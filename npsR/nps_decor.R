VERSION <- "1.0.0"

cat("Non-Parametric Shrinkage", VERSION, "\n")

# Cut-off for small lambda
LAMBDA.CO <- 0.5

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
    stop("Usage: Rscript nps_decor.R <work dir> <chrom> <WINSHIFT>")
}

tempprefix <- paste(cargs[1], "/", sep='')

# Read in args.RDS setting
args <- readRDS(paste(tempprefix, "args.RDS", sep=''))

Nt <- args[["Nt"]]
summstatfile <- args[["summstatfile"]] 
traindir <- args[["traindir"]] 
traintag <- args[["traintag"]]
WINSZ <- args[["WINSZ"]]

# Rest of command args
CHR <- as.numeric(cargs[2])

if (!(CHR %in% 1:22)) {
    stop("invalid chrom:", CHR)
}

WINSHIFT <- as.numeric(cargs[3])

if (is.nan(WINSHIFT) || WINSHIFT < 0 || WINSHIFT >= WINSZ) {
    stop("Invalid window shift:", cargs[3])
}

#########################################################################

# Read summary stats (discovery cohort)
summstat <- read.delim(summstatfile, header=TRUE, stringsAsFactors=FALSE,
                       sep="\t")
# dim(summstat)

cat("Converting effect sizes relative to standardized genotypes...")

se <- sqrt(2 * summstat$reffreq * (1 - summstat$reffreq))

std.effalt <- summstat$effalt * se

cat(" OK\n")

########################

etahat.all <- c()
eval.all <- c()                         # lambda

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

    cat("Processing chr", chrom, "...\n")
#    cat("Starting from snpIdx offset", snpIdx0, "\n")
    
    
    # stdgt dosage (uncompressed)    
#    stdgt.file <-
#       file(paste(traindir, "/chrom", chrom, ".", traintag, ".stdgt",
#                  sep=''), open="rb")

    stdgt.file <-
        gzfile(paste(traindir, "/chrom", chrom, ".", traintag,
                     ".stdgt.gz", sep=''), open="rb")
    
    I <- 1 
    snpIdx <- 1
    X0 <- NULL

    while ((snpIdx + WINSZ) <= M.chr) {

        if (I == 1 & WINSHIFT > 0) {
            span <- WINSZ - WINSHIFT

            cat("size of first window: m=", span, "\n")
            
        } else {
            span <- WINSZ
        }

        cat(I, "\n")

        X0v <- readBin(stdgt.file, "double", n=(span * Nt))
        X0 <- cbind(X0, matrix(X0v, nrow=Nt, ncol=span))
        rm(X0v)

        pos.cur <-
            summstat$pos[(snpIdx0 + snpIdx):(snpIdx0 + snpIdx + span - 1)]
        
        betahat.cur <-
            std.effalt[(snpIdx0 + snpIdx):(snpIdx0 + snpIdx + span - 1)]

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

        s0 <- eigen(ld0, symmetric=TRUE)

        ASSERT(!is.null(s0))
        rm(ld0)

        Q0 <- s0$vectors

        lambda0 <- s0$values[s0$values > LAMBDA.CO]
        
        Q0 <- Q0[, s0$values > LAMBDA.CO]

        # in contrast to manuscript, qX0 was not divided by sqrt(lambda0) here
        # thus, backconversion code to beta scale does not divide by
        # sqrt(lambda0) either.
        qX0 <- X0 %*% Q0

        etahat0 <- t(Q0) %*% as.matrix(betahat.cur) / sqrt(lambda0)

        windata <- list()
    
        windata[["eigen"]] <- s0
        windata[["Q0.X"]] <- qX0
        windata[["etahat0"]] <- etahat0

        if (WINSHIFT == 0) {
            winfilepre <-
                paste(tempprefix, "win.", chrom, ".", I, sep='')
        } else {
            winfilepre <-
                paste(tempprefix, "win_", WINSHIFT, ".", chrom, ".", I,
                      sep='')
        }
        
        saveRDS(windata, file=paste(winfilepre, ".RDS", sep=''))
        rm(windata)

        write.table(
            data.frame(lambda=lambda0, etahat=etahat0),
            file=paste(winfilepre, ".table", sep=''),
            quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t")

        saveRDS(Q0, paste(winfilepre, ".Q.RDS", sep=''))        

        etahat.all <- c(etahat.all, etahat0)
        eval.all <- c(eval.all, lambda0)
        

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

        cat("size of last window: m=", span, "\n")
            
        X0v <- readBin(stdgt.file, "double", n=(span * Nt))
        X0 <- cbind(X0, matrix(X0v, nrow=Nt, ncol=span))
        rm(X0v)

        pos.cur <-
            summstat$pos[(snpIdx0 + snpIdx):(snpIdx0 + snpIdx + span - 1)]
        
        betahat.cur <-
            std.effalt[(snpIdx0 + snpIdx):(snpIdx0 + snpIdx + span - 1)]

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

        s0 <- eigen(ld0, symmetric=TRUE)

        ASSERT(!is.null(s0))
        rm(ld0)

        Q0 <- s0$vectors

        lambda0 <- s0$values[s0$values > LAMBDA.CO]
        
        Q0 <- Q0[, s0$values > LAMBDA.CO]

        qX0 <- X0 %*% Q0

        etahat0 <- t(Q0) %*% as.matrix(betahat.cur) / sqrt(lambda0)

        windata <- list()
    
        windata[["eigen"]] <- s0
        windata[["Q0.X"]] <- qX0
        windata[["etahat0"]] <- etahat0

        if (WINSHIFT == 0) {
            winfilepre <-
                paste(tempprefix, "win.", chrom, ".", I, sep='')
        } else {
            winfilepre <-
                paste(tempprefix, "win_", WINSHIFT, ".", chrom, ".", I,
                      sep='')
        }
        
        saveRDS(windata, file=paste(winfilepre, ".RDS", sep=''))
        rm(windata)

        write.table(
            data.frame(lambda=lambda0, etahat=etahat0),
            file=paste(winfilepre, ".table", sep=''),
            quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t")

        saveRDS(Q0, paste(winfilepre, ".Q.RDS", sep=''))        

        etahat.all <- c(etahat.all, etahat0)
        eval.all <- c(eval.all, lambda0)
    } else {
        span <- M.chr - snpIdx + 1
        cat("Ignore SNPs at the end of chromosome: m=", span, "\n")
    }

    snpIdx0 <- snpIdx0 + M.chr

    close(stdgt.file)    
}

cat("Done\n")
