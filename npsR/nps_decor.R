VERSION <- "1.1"

cat("Non-Parametric Shrinkage", VERSION, "\n")

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
LAMBDA.CO <- args[["LAMBDA.CO"]]

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
summstat.chr <- read.delim(paste(summstatfile, ".", CHR, sep=''),
                           header=TRUE, stringsAsFactors=FALSE,
                           sep="\t")
# dim(summstat.chr)


########################

etahat.all <- c()
eval.all <- c()                         # lambda


M.chr <- nrow(summstat.chr)


cat("Processing chr", CHR, "...\n")
    
## stdgt dosage (uncompressed)    
#    stdgt.file <-
#       file(paste(traindir, "/chrom", CHR, ".", traintag, ".stdgt",
#                  sep=''), open="rb")

stdgt.file <-
    gzfile(paste(traindir, "/chrom", CHR, ".", traintag,
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
        summstat.chr$pos[snpIdx:(snpIdx + span - 1)]
        
    betahat.cur <-
        summstat.chr$std.effalt[snpIdx:(snpIdx + span - 1)]

    ld0 <- (t(X0) %*% X0) / (Nt - 1)

#        for (J in 1:span) {
#
#            pos.J <- pos.cur[J]
#            
#            for (K in 1:span) {
#
#                pos.dist <- abs(pos.J - pos.cur[K])
#                
#                if (pos.dist > 500000) {
#                    ld0[J, K] <- 0
#                }
#            }
#        }

    ## SE ~ 1/sqrt(Nt), 5 SD
    ld0[abs(ld0) < 5 / sqrt(Nt)] <- 0

    s0 <- eigen(ld0, symmetric=TRUE)

    ASSERT(!is.null(s0))
    rm(ld0)
    
    Q0 <- s0$vectors

    lambda0 <- s0$values[s0$values > LAMBDA.CO]
    
    Q0 <- Q0[, s0$values > LAMBDA.CO, drop=FALSE]

    # in contrast to manuscript, qX0 was not divided by sqrt(lambda0) here
    # thus, backconversion code to beta scale does not divide by
    # sqrt(lambda0) either.
    qX0 <- X0 %*% Q0

    if (length(lambda0) == 0) {
        etahat0 <- c()
    } else {
        etahat0 <- t(Q0) %*% as.matrix(betahat.cur) / sqrt(lambda0)
    }

    windata <- list()
    
    windata[["eigen"]] <- s0
    windata[["Q0.X"]] <- qX0
    windata[["etahat0"]] <- etahat0

    winfilepre <-
        paste(tempprefix, "win_", WINSHIFT, ".", CHR, ".", I,
              sep='')
        
    saveRDS(windata, file=paste(winfilepre, ".RDS", sep=''))
    rm(windata)

    write.table(
        data.frame(lambda=lambda0, etahat=etahat0),
        file=paste(winfilepre, ".table", sep=''),
        quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t")

    saveRDS(Q0, paste(winfilepre, ".Q.RDS", sep=''))        

    etahat.all <- c(etahat.all, etahat0)
    eval.all <- c(eval.all, lambda0)
    

    ## move on to next iteration
    I <- I + 1
    snpIdx <- snpIdx + span
    X0 <- NULL
    gc()
}

## last chunk
min.span <- max(WINSZ / 40, 10)
    
if ((snpIdx + min.span) <= M.chr) {
        
    span <- M.chr - snpIdx + 1

    cat("size of last window: m=", span, "\n")
            
    X0v <- readBin(stdgt.file, "double", n=(span * Nt))
    X0 <- cbind(X0, matrix(X0v, nrow=Nt, ncol=span))
    rm(X0v)

    pos.cur <-
        summstat.chr$pos[snpIdx:(snpIdx + span - 1)]
        
    betahat.cur <-
        summstat.chr$std.effalt[snpIdx:(snpIdx + span - 1)]
    
    ld0 <- (t(X0) %*% X0) / (Nt - 1)

#        for (J in 1:span) {
#            
#            pos.J <- pos.cur[J]
#            
#            for (K in 1:span) {
#
#                pos.dist <- abs(pos.J - pos.cur[K])
#                
#                if (pos.dist > 500000) {
#                    ld0[J, K] <- 0
#                }
#            }
#        }

    ## SE ~ 1/sqrt(Nt), 5 SD
    ld0[abs(ld0) < 5 / sqrt(Nt)] <- 0
    
    s0 <- eigen(ld0, symmetric=TRUE)
    
    ASSERT(!is.null(s0))
    rm(ld0)
    
    Q0 <- s0$vectors
    
    lambda0 <- s0$values[s0$values > LAMBDA.CO]
    
    Q0 <- Q0[, s0$values > LAMBDA.CO, drop=FALSE]

    qX0 <- X0 %*% Q0

    if (length(lambda0) == 0) {
        etahat0 <- c()
    } else {
        etahat0 <- t(Q0) %*% as.matrix(betahat.cur) / sqrt(lambda0)
    }

    windata <- list()
    
    windata[["eigen"]] <- s0
    windata[["Q0.X"]] <- qX0
    windata[["etahat0"]] <- etahat0
    
    winfilepre <-
        paste(tempprefix, "win_", WINSHIFT, ".", CHR, ".", I,
              sep='')
        
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

close(stdgt.file)    

cat("Done\n")
