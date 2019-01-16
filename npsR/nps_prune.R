VERSION <- "1.0.0"

cat("Non-Parametric Shrinkage", VERSION, "\n")

# Cut-off for corss-window pruning
CXWCOR.CO <- 0.3

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
    stop("Usage: Rscript nps_prune.R <work dir> <chrom> <WINSHIFT>")
}

tempprefix <- paste(cargs[1], "/", sep='')

# Read in saved settings
args <- readRDS(paste(tempprefix, "args.RDS", sep=''))
WINSZ <- args[["WINSZ"]]

CHR <- as.numeric(cargs[2])

if (!(CHR %in% 1:22)) {
    stop("invalid chrom", CHR)
}

WINSHIFT <- as.numeric(cargs[3])

if (is.nan(WINSHIFT) || WINSHIFT < 0 || WINSHIFT >= WINSZ) {
    stop("Invalid shift:", cargs[3])
}


#############################################################################

chrom <- CHR

# cat("chr", chrom, "\n")

I <- 1
snpIdx <- 1

# central block
if (WINSHIFT == 0) {
    winfilepre <-
        paste(tempprefix, "win.", chrom, ".", I, sep='')
} else {
    winfilepre <-
        paste(tempprefix, "win_", WINSHIFT, ".", chrom, ".", I,
              sep='')
}

ASSERT(file.exists(paste(winfilepre, ".RDS", sep='')))

windata <- readRDS(file=paste(winfilepre, ".RDS", sep=''))
QX0 <- windata[["Q0.X"]]
Nq <- ncol(QX0)

wintab0 <- read.delim(paste(winfilepre, ".table", sep=''),
                      header=TRUE, sep="\t")

# read block on the right
if (WINSHIFT == 0) {
    winfilepre <-
        paste(tempprefix, "win.", chrom, ".", (I + 1), sep='')
} else {
    winfilepre <-
        paste(tempprefix, "win_", WINSHIFT, ".", chrom, ".", (I + 1),
              sep='')
}

ASSERT(file.exists(paste(winfilepre, ".RDS", sep='')))
    
windata <- readRDS(file=paste(winfilepre, ".RDS", sep=''))
QX0 <- cbind(QX0, windata[["Q0.X"]])

wintabR <- read.delim(paste(winfilepre, ".table", sep=''),
                      header=TRUE, sep="\t")

# calc cross-window corr
ld.q <- cor(QX0)

ld.qx <- ld.q[1:Nq, (Nq + 1):ncol(ld.q)]
    
corcx0 <- apply(ld.qx, 1, function(x) max(abs(x)))
corcxR.idx <- apply(ld.qx, 1, function(x) which.max(abs(x))) 

# prune
for (I0 in which(corcx0 > CXWCOR.CO)) {
        
    IR <- corcxR.idx[I0]

    if (abs(wintab0$etahat[I0]) < abs(wintabR$etahat[IR])) {
            
        cat(wintab0$etahat[I0], "<", wintabR$etahat[IR], ":",
            ld.qx[I0, IR], "\n")
        
        wintab0[I0, ] <- 0 
    }
}

# re-write    
if (WINSHIFT == 0) {
    winfilepre <-
        paste(tempprefix, "win.", chrom, ".", I, ".pruned", sep='')
} else {
    winfilepre <-
        paste(tempprefix, "win_", WINSHIFT, ".", chrom, ".", I, ".pruned",
              sep='')
}
    
write.table(wintab0, 
            file=paste(winfilepre, ".table", sep=''),
            quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t")

# shift
NqL <- Nq
wintabL <- wintab0
wintab0 <- wintabR

I <- I + 2
    
# next block
if (WINSHIFT == 0) {
    winfilepre <-
        paste(tempprefix, "win.", chrom, ".", I, sep='')
} else {
    winfilepre <-
        paste(tempprefix, "win_", WINSHIFT, ".", chrom, ".", I,
              sep='')
}
    
while (file.exists(paste(winfilepre, ".RDS", sep=''))) {

    cat(I, "\n")

    # read block on the right        
    windata <- readRDS(file=paste(winfilepre, ".RDS", sep=''))
    QX0 <- cbind(QX0, windata[["Q0.X"]])
    NqR <- ncol(windata[["Q0.X"]])
    
    wintabR <- read.delim(paste(winfilepre, ".table", sep=''),
                          header=TRUE, sep="\t")

    Nq <- ncol(QX0) - NqL - NqR
    
    ld.q <- cor(QX0)
    
    ld.qx <- ld.q[(NqL + 1):(NqL + Nq),
                  c(1:NqL, (NqL + Nq + 1):ncol(ld.q))]
    
    corcx0 <- apply(ld.qx, 1, function(x) max(abs(x)))
    corcxLR.idx <- apply(ld.qx, 1, function(x) which.max(abs(x)))

    # prune
    for (I0 in which(corcx0 > CXWCOR.CO)) {
            
        ILR <- corcxLR.idx[I0]

        if (ILR <= NqL) {
            IL <- ILR
            
            if (abs(wintab0$etahat[I0]) < abs(wintabL$etahat[IL])) {
                    
                cat(wintab0$etahat[I0], "<", wintabL$etahat[IL], ":",
                    ld.qx[I0, IL], "\n")
                    
                wintab0[I0, ] <- 0 
            }
        } else {
            IR <- ILR - NqL
            
            if (abs(wintab0$etahat[I0]) < abs(wintabR$etahat[IR])) {
                
                cat(wintab0$etahat[I0], "<", wintabR$etahat[IR], ":",
                    ld.qx[I0, ILR], "\n")
                
                wintab0[I0, ] <- 0 
            }
        }
    }

    # re-write    
    if (WINSHIFT == 0) {
        winfilepre <-
            paste(tempprefix, "win.", chrom, ".", (I - 1), ".pruned",
                  sep='')
    } else {
        winfilepre <-
            paste(tempprefix, "win_", WINSHIFT, ".", chrom, ".", (I - 1),
                  ".pruned", sep='')
    }
    
    write.table(wintab0, 
                file=paste(winfilepre, ".table", sep=''),
                quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t")
        
        
    # move on to next iteration
    QX0 <- QX0[, (NqL + 1):ncol(QX0)]
        
    NqL <- Nq
    wintabL <- wintab0
    wintab0 <- wintabR
    
    I <- I + 1
    
    # next block
    if (WINSHIFT == 0) {
        winfilepre <-
            paste(tempprefix, "win.", chrom, ".", I, sep='')
    } else {
        winfilepre <-
            paste(tempprefix, "win_", WINSHIFT, ".", chrom, ".", I,
                  sep='')
    }        
}

# last block
Nq <- ncol(QX0) - NqL
    
ld.q <- cor(QX0)

ld.qx <- ld.q[(NqL + 1):(NqL + Nq), 1:NqL]
                  
corcx0 <- apply(ld.qx, 1, function(x) max(abs(x)))
corcxL.idx <- apply(ld.qx, 1, function(x) which.max(abs(x)))

# prune
for (I0 in which(corcx0 > CXWCOR.CO)) {
            
    IL <- corcxL.idx[I0]
    
    if (abs(wintab0$etahat[I0]) < abs(wintabL$etahat[IL])) {
        
        cat(wintab0$etahat[I0], "<", wintabL$etahat[IL], ":",
            ld.qx[I0, IL], "\n")
        
        wintab0[I0, ] <- 0 
    }
}

# re-write    
if (WINSHIFT == 0) {
    winfilepre <-
        paste(tempprefix, "win.", chrom, ".", (I - 1), ".pruned",
              sep='')
} else {
    winfilepre <-
        paste(tempprefix, "win_", WINSHIFT, ".", chrom, ".", (I - 1),
              ".pruned", sep='')
}
    
write.table(wintab0, 
            file=paste(winfilepre, ".table", sep=''),
            quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t")

cat("Done\n")
