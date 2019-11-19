VERSION <- "1.1"

cat("Non-Parametric Shrinkage", VERSION, "\n")

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

if (length(cargs) != 2) {
    stop("Usage: Rscript nps_split_gwassig.R <work dir> <chrom>")
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

#########################################################################

# Read summary stats (discovery)
summstat <- read.delim(summstatfile, header=TRUE, stringsAsFactors=FALSE,
                       sep="\t")
# dim(summstat)

## se <- sqrt(2 * summstat$reffreq * (1 - summstat$reffreq))
## std.effalt <- summstat$effalt * se

std.effalt <-
    abs(qnorm(summstat$pval/2, lower.tail=TRUE))*sign(summstat$effalt)

summstat <- cbind(summstat, std.effalt=std.effalt)

chr.str <- paste('chr', CHR, sep='')

summstat.chr <- summstat[summstat$chr == chr.str, ]

M.chr <- nrow(summstat.chr)


#####
## Scan starting from the most associated SNPs

pval.chr <- summstat.chr$pval

start.idx <- c()
end.idx <- c()
betahat.tail.chr <- rep(0, M.chr)

# for trPT of tails 
X.trPT <- NULL
betahat.trPT <- c()
    
while (min(pval.chr) < TAIL.THR) {

    ## Focal SNP
    pick1 <- which.min(pval.chr)[1]
    pick1.chrpos <- summstat.chr$pos[pick1]

    cat(pick1.chrpos, "-", min(pval.chr), "\n")
            
    betahat.tail.chr[pick1] <- summstat.chr$std.effalt[pick1]
        
    begin <- max(which.min(pval.chr)[1] - WINSZ, 1)
    end <- min(which.min(pval.chr)[1] + WINSZ, M.chr)

    start.idx <- c(start.idx, begin)
    end.idx <- c(end.idx, end)

    ## LD-adjust the std.effalt around the focal SNP

    stdgt.file <-
        gzfile(paste(traindir, "/chrom", CHR, ".", traintag,
                     ".stdgt.gz", sep=''), open="rb")

    ## Seek to the "begin" pos
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

    ## Read training genotypes
    span <- (end - begin + 1)
    X0v <- readBin(stdgt.file, "double", n=(span * Nt))
    X0 <- matrix(X0v, nrow=Nt, ncol=span)
    rm(X0v)

    ## Get the reference LD matrix
    pos.win <- summstat.chr$pos[begin:end]

    ld0 <- (t(X0) %*% X0) / (Nt - 1)

#    for (J in 1:span) {
#            
#        pos.J <- pos.win[J]
#            
#        for (K in 1:span) {
#
#            pos.dist <- abs(pos.J - pos.win[K])
#                
#            if (pos.dist > 500000) {
#                ld0[J, K] <- 0
#            }
#        }
#    }

    # SE ~ 1/sqrt(Nt), 5 SD
    ld0[abs(ld0) < 5 / sqrt(Nt)] <- 0

    ## Treat the tail effect as a fixed effect and correct it
    pick1 <- which.min(pval.chr[begin:end])[1]
    betahat.win <- summstat.chr$std.effalt[begin:end]
    
    X.trPT <- cbind(X.trPT, X0[, pick1])
    betahat.trPT <- c(betahat.trPT, betahat.win[pick1])
            
    # calculate residual effects
    tailbeta <- rep(0, span)
    tailbeta[pick1] <- betahat.win[pick1]

    # fix for numerical error in ld0[pick1, pick1]
    ld0 <- ld0 / ld0[pick1, pick1]
            
    betahat.win.tailfix <-
        betahat.win - ld0 %*% as.matrix(tailbeta)

#        cat("tail fix: residual=", betahat.win.tailfix[pick1], "\n")
#        print(ld0 %*% as.matrix(tailbeta))

    summstat.chr$std.effalt[begin:end] <- betahat.win.tailfix

    rm(X0)
    close(stdgt.file)
    gc()

    ## Mask out LD neighbors of focal SNPs from being selected as
    ## next focal SNPs
    
    pos.mask <- pos.win[abs(ld0[pick1, ]) > 0.3]

#    cat("Mask out -",
#        paste(pval.chr[summstat.chr$pos %in% pos.mask], collapse=", "),
#        "\n")

    pval.chr[summstat.chr$pos %in% pos.mask] <- 1

}

## Save betahat info for the tail partition

print(data.frame(start=start.idx, end=end.idx,
                 bp=(summstat.chr$pos[end.idx] - summstat.chr$pos[start.idx])))

write.table(data.frame(betahat.tail.chr),
            file=paste(tempprefix, "tail_betahat.", CHR, ".table",
                       sep=''),
            row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")

## Save updated std.effalt

write.table(summstat.chr,
            file=paste(summstatfile, ".", CHR, sep=''),
            row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")

## Save tail 
if (length(betahat.trPT) > 0) {
    cat("Saving trPT for tail\n")
        
    ASSERT(nrow(X.trPT) == Nt)
    ASSERT(ncol(X.trPT) == length(betahat.trPT))

    prs.tail <- X.trPT %*% as.matrix(betahat.trPT)

    ASSERT(length(prs.tail) == Nt)

    saveRDS(prs.tail, file=paste(tempprefix, "trPT.", CHR, ".tail.RDS",
                                 sep=''))
}

cat("Done\n")

