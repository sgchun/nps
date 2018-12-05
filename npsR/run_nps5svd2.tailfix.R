library(MASS)

source("nps.1.2b/nps.R")

#########################################################################

# WINSZ <- 4000
LAMBDA.CO <- 0.5
TAIL.THR <- 5E-8

#########################################################################

cargs <- commandArgs(trailingOnly=TRUE)

if (length(cargs) != 3) {
    stop("Rscript run_nps5svd2.R <outdir> <chrom> <WINSHIFT>")
}

tempprefix <- paste(cargs[1], "/", sep='')

# Read in saved settings
args <- readRDS(paste(tempprefix, "args.RDS", sep=''))

Nt <- args[["Nt"]]
summstatfile <- args[["summstatfile"]] 
traindir <- args[["traindir"]] 
traintag <- args[["traintag"]]
betaprefix <- args[["betaprefix"]]      # known truth (only for diagnostics)
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

#tempprefix <- 
#    "/net/home/schun/gwassim2/mvn2/train/nps5.23301.4000/23301."

# WINSHIFT <- 0
# WINSHIFT <- 1000
# WINSHIFT <- 2000
# WINSHIFT <- 3000


#########################################################################

# Read summary stats (discovery)
summstat <- read.delim(summstatfile, header=TRUE, stringsAsFactors=FALSE,
                       sep="\t")
dim(summstat)

se <- sqrt(2 * summstat$reffreq * (1 - summstat$reffreq))

std.effalt <- summstat$effalt * se

if (betaprefix != ".") {

    betafile <- paste(betaprefix, "beta.chrom1.txt", sep='')

    ASSERT(file.exists(betafile))

    beta <-
        read.table(betafile, 
                   sep="\t", header=FALSE)[, 1]
    
    for (chrom in 2:22) {
        beta.file <-
            paste(betaprefix, "beta.chrom", chrom, ".txt", sep='')

        ASSERT(file.exists(betafile))

        beta0 <- 
            read.table(beta.file, sep="\t", header=FALSE)[, 1]
                   
        beta <- c(beta, beta0)
    }

    ASSERT(length(beta) == nrow(summstat))

    # from per-allele effects to standardized effects
    beta <- beta * se

    dhs <- rep(0, length(beta))

    # DHS
    dhsfile <- paste(betaprefix, "DHS.txt", sep='')

    if (file.exists(dhsfile)) {
        dhs <- 
            read.table(dhsfile,
                       sep="\t", header=FALSE)[, 1]
        ASSERT(length(beta) == length(dhs))
        mean(dhs == 1)
    }

} else {
    
    beta <- rep(0, nrow(summstat))
    dhs <- rep(0, nrow(summstat))
    
}

########################

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

    print(chrom)
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
            gzfile(paste(traindir, "/chrom", chrom, ".maf_0.05.", traintag,
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
#    stdgt.file <-
#       file(paste(traindir, "/chrom", chrom, ".maf_0.05.", traintag, ".stdgt",
#                  sep=''), open="rb")

    stdgt.file <-
        gzfile(paste(traindir, "/chrom", chrom, ".maf_0.05.", traintag,
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
            ASSERT(ncol(df.tailfix.pruned) == 4)
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

    # Min 100 SNPs
    # last chunk
    if ((snpIdx + 100) <= M.chr) {
        
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
            ASSERT(ncol(df.tailfix.pruned) == 4)
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

        #else {
        #    saveRDS(prs.tail, 
        #            file=paste(tempprefix, "win_", WINSHIFT, ".trPT.", CHR,
        #                       ".tail.RDS", sep=''))
        #                   
        #}
    }
}

cat("Done\n")

q(save="no")

