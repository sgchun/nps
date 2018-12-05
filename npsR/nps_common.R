library(MASS)
library(methods)

do.nothing <- function(x) { }

VERBOSE <- print
#VERBOSE <- do.nothing

ASSERT <- function(test) {
    if (length(test) == 0) {
        stop(paste("ASSERT fail for empty conditional:",
                   deparse(substitute(test))))
    }

    if (is.na(test)) {
        stop(paste("ASSERT fail for missing value:", deparse(substitute(test))))
    }
    
    if (!test) {
        stop(paste("ASSERT fail:", deparse(substitute(test))))
    }
}

w.quant.o1 <- function(x, nbin) {

    ASSERT(all(x >= 0))
    ASSERT(all(!is.na(x)))

    x <- sort(x, decreasing=FALSE)

    w <- x / sum(x)                     # normalize
    
    sum.wx <- cumsum(w)

    cuts <- sum(w) * (1:(nbin - 1) / nbin)

    q <- c()

    for (cx in cuts) {
        xL <- max(x[sum.wx < cx])
        xR <- min(x[sum.wx > cx])
        q <- c(q, (xL + xR)/2)
    }
    
    c(min(x), q, max(x))
}

w.quant.o2 <- function(x, nbin) {

    ASSERT(all(x >= 0))
    ASSERT(all(!is.na(x)))

    x <- sort(x, decreasing=FALSE)

    w <- x / sum(x)                     # normalize
    
    sum.wx <- cumsum(w * x)

    cuts <- sum(w * x) * (1:(nbin - 1) / nbin)

    q <- c()

    for (cx in cuts) {
        xL <- max(x[sum.wx < cx])
        xR <- min(x[sum.wx > cx])
        q <- c(q, (xL + xR)/2)
    }
    
    c(min(x), q, max(x))
}

w.quant.o3 <- function(x, nbin) {

    ASSERT(all(x >= 0))
    ASSERT(all(!is.na(x)))

    x <- sort(x, decreasing=FALSE)

    w <- x / sum(x)                     # normalize
    
    sum.wx <- cumsum(w * w * x)

    cuts <- sum(w * w * x) * (1:(nbin - 1) / nbin)

    q <- c()

    for (cx in cuts) {
        xL <- max(x[sum.wx < cx])
        xR <- min(x[sum.wx > cx])
        q <- c(q, (xL + xR)/2)
    }
    
    c(min(x), q, max(x))
}

w.quant.o4 <- function(x, nbin) {

    ASSERT(all(x >= 0))
    ASSERT(all(!is.na(x)))

    x <- sort(x, decreasing=FALSE)

    w <- x / sum(x)                     # normalize
    
    sum.wx <- cumsum((w * x)** 2)

    cuts <- sum((w * x)**2) * (1:(nbin - 1) / nbin)

    q <- c()

    for (cx in cuts) {
        xL <- max(x[sum.wx < cx])
        xR <- min(x[sum.wx > cx])
        q <- c(q, (xL + xR)/2)
    }
    
    c(min(x), q, max(x))
}

regularize.ldmat <- function(ld0, pos, bandbp) {
    
    M.cur <- nrow(ld0)

    ASSERT(length(pos) == M.cur)
    
        
    for (J in 1:M.cur) {
        
        pos.J <- pos[J]
            
        for (K in 1:M.cur) {

            pos.dist <- abs(pos.J - pos[K])
            
            if (pos.dist > bandbp) {
                ld0[J, K] <- 0
            }
        }
    }

    # project to PSD matrix
    s <- eigen(ld0, symmetric=TRUE)

    Ds <- matrix(0, nrow=M.cur, ncol=M.cur)
    diag(Ds) <- pmax(s$values, 0)
    Q <- s$vectors

    ld0.psd <- Q %*% Ds %*% t(Q)

    if (!isSymmetric(ld0.psd)) {
        ld0.psd <- (ld0.psd + t(ld0.psd)) / 2
    }
    
    return(ld0.psd)
}


############################################################################
# Set up windows and calculate reference LD for each window
# Run it for each chromosome
# gtchr: genotype matrix (individual x marker, standardized to mean 0, var 1)

reference.ld <-
    function(gtchr, window.size, window.offset=0, calc.cross.window=FALSE) {

        I <- 1
        snpIdx <- 1

        M <- ncol(gtchr)                   # total SNPs
        
        M.prev <- 0                     # mark first window

        ld0.list <- list()
        ld0prev.list <- list()

        while (snpIdx <= M) {

            #VERBOSE(I)
    
            if (M.prev == 0) {
                endIdx <- min(M, snpIdx + window.size + window.offset - 1)
            } else {
                endIdx <- min(M, snpIdx + window.size - 1)
            }

            M.local <- endIdx - snpIdx + 1

            gt <- gtchr[, snpIdx:endIdx, drop=FALSE]

            ASSERT( ncol(gt) > 0 )
            ASSERT( ncol(gt) == M.local)

            ld0 <- cor(gt)

            af <- apply(gt, 1, mean)

            if (sum(af == 0 | af == 1) > 0) {
                VERBOSE(paste("filling MAF=0: ", sum(af == 0 | af == 1)))

                for (I in which(af == 0 | af == 1)) {
                    ld0[I, ] <- 0
                    ld0[, I] <- 0
                    ld0[I, I] <- 1
                }
            }

            ASSERT( sum(is.na(ld0)) == 0 )
    
            ld0.list[[I]] <- ld0

            # ld0 prev
            if (M.prev == 0 || !calc.cross.window) {
        
                ld0prev <- matrix(0, nrow=0, ncol=0)
        
            } else {
    
                gt <- gtchr[, (snpIdx - M.prev):endIdx, drop=FALSE]

                ld0prev <- cor(gt)

                ASSERT( sum(is.na(ld0prev)) == 0 )
                ASSERT( nrow(ld0prev) > 1 )
                ASSERT( nrow(ld0prev) == (M.prev + M.local))        

                ld0prev <- ld0prev[1:M.prev, (M.prev + 1):nrow(ld0prev)]

            }

            ld0prev.list[[I]] <- ld0prev

            I <- I + 1
            snpIdx <- snpIdx + M.local

            M.prev <- M.local
        }                    
        
        ASSERT(snpIdx == (M + 1))
        ASSERT(endIdx == M)

        if (!calc.cross.window) {
            ld0prev.list <- NULL
        }

        return (list(ld0=ld0.list, ld0prev=ld0prev.list))
    }



############################################################################
# Average over sliding windows

perSNP.effect.avg <- function(nps1, nps2, nps3, nps4) {

    return (apply(cbind(perSNP.effect(nps1),
                        perSNP.effect(nps2),
                        perSNP.effect(nps3),
                        perSNP.effect(nps4)), 1, mean))

}

############################################################################
# Predict risk
# gt: genotype matrix 
# (marker x individual, standardized to mean 0, var 1)

predict.phenotype <- function(gt, beta) {

    ASSERT(nrow(gt) == length(beta))

    return (apply(gt, 2, function (x) sum(x * beta)))
}
