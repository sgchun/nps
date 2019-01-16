VERSION <- "1.0.1"

cat("Non-Parametric Shrinkage", VERSION, "\n")

# Cut-off for small lambda
LAMBDA.CO <- 0.5
# P-value threshold for the GWAS-significant tail partition
TAIL.THR <- 5E-8

ASSERT <- function(test) {
    if (length(test) == 0) {
        stop(paste("ASSERT fail for empty conditional:",
                   deparse(substitute(test))))
    }

    if (is.na(test)) {
        stop(paste("ASSERT fail for missing value:",
                   deparse(substitute(test))))
    }
    
    if (!test) {
        stop(paste("ASSERT fail:", deparse(substitute(test))))
    }
}

#########################################################################

cargs <- commandArgs(trailingOnly=TRUE)

if (length(cargs) != 3) {
    stop("Usage: Rscript nps_split_gwassig.R <work dir> <chrom> <WINSHIFT>")
}

tempprefix <- paste(cargs[1], "/", sep='')

# Read in saved settings
args <- readRDS(paste(tempprefix, "args.RDS", sep=''))

Nt <- args[["Nt"]]
summstatfile <- args[["summstatfile"]] 
traindir <- args[["traindir"]] 
traintag <- args[["traintag"]]
WINSZ <- args[["WINSZ"]]

# Rest of command args
CHR <- as.numeric(cargs[2])

if (!(CHR %in% 1:22)) {
    stop("invalid chrom", CHR)
}

WINSHIFT <- as.numeric(cargs[3])

if (is.nan(WINSHIFT) || WINSHIFT < 0 || WINSHIFT >= WINSZ) {
    stop("Invalid shift:", cargs[3])
}

#########################################################################

# Read summary stats (discovery)
summstat <- read.delim(summstatfile, header=TRUE, stringsAsFactors=FALSE,
                       sep="\t")
# dim(summstat)

se <- sqrt(2 * summstat$reffreq * (1 - summstat$reffreq))

std.effalt <- summstat$effalt * se

chrom <- 1
snpIdx0 <- 0

for (chrom in 1:22) {

    M.chr <- sum(summstat$chr == paste('chr', chrom, sep=''))

    if (chrom > CHR) {
        break 
    }

    if (chrom < CHR) {
        snpIdx0 <- snpIdx0 + M.chr
        next
    }

#    print(chrom)
    cat("Starting from snpIdx offset", snpIdx0, "\n")

    # Scan the tail-p-value regions 
    summstat.chr <- summstat[summstat$chr == paste('chr', chrom, sep=''), ]
    std.effalt.chr <- std.effalt[(snpIdx0 + 1):(snpIdx0 + M.chr)]
    pval.chr <- summstat.chr$pval
   
    start.idx <- c()
    end.idx <- c()
    betahat.tail.chr <- rep(0, M.chr)
    
    while (min(pval.chr) < TAIL.THR) {

        pick1 <- which.min(pval.chr)[1]
        pick1.chrpos <- summstat.chr$pos[pick1]
            
        betahat.tail.chr[pick1] <- std.effalt.chr[pick1]
        
        begin <- max(which.min(pval.chr)[1] - WINSZ, 1)
        end <- min(which.min(pval.chr)[1] + WINSZ, nrow(summstat.chr))

        start.idx <- c(start.idx, begin)
        end.idx <- c(end.idx, end)
            
#        pval.chr[begin:end] <- 1
        # 500 kb window
        pval.chr[ abs(summstat.chr$pos - pick1.chrpos) < 500000 ] <- 1
    }

    print(data.frame(start=start.idx, end=end.idx,
              bp=(summstat.chr$pos[end.idx] - summstat.chr$pos[start.idx])))

    if (WINSHIFT == 0) {
        write.table(data.frame(betahat.tail.chr),
                    file=paste(tempprefix, "tail_betahat.", chrom, ".table",
                               sep=''),
                    row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")

    } 

    if (length(start.idx) == 0) {
        cat("No tail\nDone\n")

        q(save="no")
    }

    # Save tail-p-value SNPs
    X.tail <- NULL
    betahat.tail <- c()

    std.effalt.chr.adjusted <- std.effalt.chr

    for (X in 1:length(start.idx)) {
        begin <- start.idx[X]
        end <- end.idx[X]

        stdgt.file <-
            gzfile(paste(traindir, "/chrom", chrom, ".", traintag,
                         ".stdgt.gz", sep=''), open="rb")

        # seek
        X2 <- 1

        while ( (X2 + 1000) < begin ) {
            readBin(stdgt.file, "double", n=(Nt*1000))
            X2 <- X2 + 1000
        }

        while (X2 < begin) {
            readBin(stdgt.file, "double", n=Nt)
            X2 <- X2 + 1
        }
        
        ASSERT(X2 == begin)

        # read genotypes
        span <- (end - begin + 1)
        X0v <- readBin(stdgt.file, "double", n=(span * Nt))
        X0 <- matrix(X0v, nrow=Nt, ncol=span)
        rm(X0v)

        pos.cur <- summstat.chr$pos[begin:end]

        ld0 <- (t(X0) %*% X0) / (Nt - 1)

        for (J in 1:span) {
            
            pos.J <- pos.cur[J]
            
            for (K in 1:span) {

                pos.dist <- abs(pos.J - pos.cur[K])
                
                if (pos.dist > 500000) {
                    ld0[J, K] <- 0
                }
            }
        }

        # SE ~ 1/sqrt(Nt), 5 SD
        ld0[abs(ld0) < 5 / sqrt(Nt)] <- 0

        # Treat the tail effect as a fixed effect and correct it
        pval.cur <- summstat.chr$pval[begin:end]
        betahat.cur <- std.effalt.chr.adjusted[begin:end]
        pick1 <- which.min(pval.cur)[1]

        X.tail <- cbind(X.tail, X0[, pick1])
        betahat.tail <- c(betahat.tail, betahat.cur[pick1])
            
        # calculate residual effects
        tailbeta <- rep(0, span)
        tailbeta[pick1] <- betahat.cur[pick1]

        # fix for numerical error in ld0[pick1, pick1]
        ld0 <- ld0 / ld0[pick1, pick1]
            
        betahat.cur.tailfix <-
            betahat.cur - ld0 %*% as.matrix(tailbeta)

#        cat("tail fix: residual=", betahat.cur.tailfix[pick1], "\n")
#        print(ld0 %*% as.matrix(tailbeta))

        std.effalt.chr.adjusted[begin:end] <- betahat.cur.tailfix
        
        close(stdgt.file)
        gc()
    }
    
    
    # stdgt dosage 
    stdgt.file <-
        gzfile(paste(traindir, "/chrom", chrom, ".", traintag,
                     ".stdgt.gz", sep=''), open="rb")
    
    I <- 1 
    snpIdx <- 1
    X0 <- NULL

    while ((snpIdx + WINSZ) <= M.chr) {

        print(I)

        if (I == 1 & WINSHIFT > 0) {
            span <- WINSZ - WINSHIFT

            cat("First WINSZ:", span, "\n")
            
        } else {
            span <- WINSZ
        }

        X0v <- readBin(stdgt.file, "double", n=(span * Nt))
        X0 <- cbind(X0, matrix(X0v, nrow=Nt, ncol=span))
        rm(X0v)

        betahat.cur <-
            std.effalt[(snpIdx0 + snpIdx):(snpIdx0 + snpIdx + span - 1)]
        betahat.cur.adjusted <-
            std.effalt.chr.adjusted[snpIdx:(snpIdx + span - 1)]
        pos.cur <-
            summstat$pos[(snpIdx0 + snpIdx):(snpIdx0 + snpIdx + span - 1)]

        if (any(betahat.cur != betahat.cur.adjusted)) {

            cat("Update window ", I, "\n")
            
            # Load projection data
            if (WINSHIFT == 0) {
                winfilepre <-
                    paste(tempprefix, "win.", chrom, ".", I, sep='')
            } else {
                winfilepre <-
                    paste(tempprefix, "win_", WINSHIFT, ".", chrom, ".", I,
                          sep='')
            }
        
            windata <- readRDS(paste(winfilepre, ".RDS", sep=''))
        
            s0 <- windata[["eigen"]]

            ASSERT(!is.null(s0))

            Q0 <- s0$vectors        
            Q0 <- Q0[, s0$values > LAMBDA.CO]
            lambda0 <- s0$values[s0$values > LAMBDA.CO]

            etahat0 <-
                t(Q0) %*% as.matrix(betahat.cur.adjusted) / sqrt(lambda0)

            rm(windata)
            
            # pruned
            if (WINSHIFT == 0) {
                winfilepre <-
                    paste(tempprefix, "win.", chrom, ".", I, ".pruned", sep='')
            } else {
                winfilepre <-
                    paste(tempprefix, "win_", WINSHIFT, ".", chrom, ".", I,
                          ".pruned", sep='')
            }

            df.tailfix.pruned <- 
                read.delim(paste(winfilepre, ".table", sep=''),
                           sep="\t", header=TRUE, stringsAsFactors=FALSE)

            ASSERT(nrow(df.tailfix.pruned) == length(lambda0))
            ASSERT(ncol(df.tailfix.pruned) == 2)
            ASSERT(all(abs(df.tailfix.pruned$lambda - lambda0) < 0.00001 | 
                       df.tailfix.pruned$lambda == 0))
            
            df.tailfix.pruned$etahat <- etahat0
            df.tailfix.pruned$etahat[df.tailfix.pruned$lambda == 0] <- 0

            write.table(df.tailfix.pruned, 
                        file=paste(winfilepre, ".tailfix.table", sep=''),
                        quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t")

        }

        
        # move on to next iteration
        I <- I + 1
        snpIdx <- snpIdx + span
        X0 <- NULL
        gc()
    }

    # last chunk
    min.span <- max(WINSZ / 40, 10)
    
    if ((snpIdx + min.span) <= M.chr) {
        
        span <- M.chr - snpIdx + 1 
            
        X0v <- readBin(stdgt.file, "double", n=(span * Nt))
        X0 <- cbind(X0, matrix(X0v, nrow=Nt, ncol=span))
        rm(X0v)

        betahat.cur <-
            std.effalt[(snpIdx0 + snpIdx):(snpIdx0 + snpIdx + span - 1)]
        betahat.cur.adjusted <-
            std.effalt.chr.adjusted[snpIdx:(snpIdx + span - 1)]
        pos.cur <-
            summstat$pos[(snpIdx0 + snpIdx):(snpIdx0 + snpIdx + span - 1)]

        if (any(betahat.cur != betahat.cur.adjusted)) {

            cat("Update window ", I, "\n")
            
            # Load projection data
            if (WINSHIFT == 0) {
                winfilepre <-
                    paste(tempprefix, "win.", chrom, ".", I, sep='')
            } else {
                winfilepre <-
                    paste(tempprefix, "win_", WINSHIFT, ".", chrom, ".", I,
                          sep='')
            }
        
            windata <- readRDS(paste(winfilepre, ".RDS", sep=''))
        
            s0 <- windata[["eigen"]]

            ASSERT(!is.null(s0))

            Q0 <- s0$vectors        
            Q0 <- Q0[, s0$values > LAMBDA.CO]
            lambda0 <- s0$values[s0$values > LAMBDA.CO]

            etahat0 <-
                t(Q0) %*% as.matrix(betahat.cur.adjusted) / sqrt(lambda0)

            rm(windata)

            # pruned
            if (WINSHIFT == 0) {
                winfilepre <-
                    paste(tempprefix, "win.", chrom, ".", I, ".pruned", sep='')
            } else {
                winfilepre <-
                    paste(tempprefix, "win_", WINSHIFT, ".", chrom, ".", I,
                          ".pruned", sep='')
            }

            df.tailfix.pruned <- 
                read.delim(paste(winfilepre, ".table", sep=''),
                           sep="\t", header=TRUE, stringsAsFactors=FALSE)

            ASSERT(nrow(df.tailfix.pruned) == length(lambda0))
            ASSERT(ncol(df.tailfix.pruned) == 2)
            ASSERT(all(abs(df.tailfix.pruned$lambda - lambda0) < 0.00001 | 
                       df.tailfix.pruned$lambda == 0))
            
            df.tailfix.pruned$etahat <- etahat0
            df.tailfix.pruned$etahat[df.tailfix.pruned$lambda == 0] <- 0

            write.table(df.tailfix.pruned, 
                        file=paste(winfilepre, ".tailfix.table", sep=''),
                        quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t")

        }        
    }

    snpIdx0 <- snpIdx0 + M.chr

    close(stdgt.file)

    # save tail 
    if (length(betahat.tail) > 0) {
        cat("Saving trPT for tail\n")
        
        ASSERT(nrow(X.tail) == Nt)
        ASSERT(ncol(X.tail) == length(betahat.tail))

        prs.tail <- X.tail %*% as.matrix(betahat.tail)

        ASSERT(length(prs.tail) == Nt)

        if (WINSHIFT == 0) {
            saveRDS(prs.tail, file=paste(tempprefix, "trPT.", CHR, ".tail.RDS",
                                         sep=''))
        }
    }
}

cat("Done\n")

