library(MASS)

source("nps.1.2b/nps.R")

#########################################################################

# WINSZ <- 4000
LAMBDA.CO <- 0.5

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

etahat.all <- c()
eta.all <- c()
eval.all <- c()                         # lambda
dhs.all <- c()

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

        pos.cur <-
            summstat$pos[(snpIdx0 + snpIdx):(snpIdx0 + snpIdx + span - 1)]
        
        beta.cur <- beta[(snpIdx0 + snpIdx):(snpIdx0 + snpIdx + span - 1)]
    
        betahat.cur <-
            std.effalt[(snpIdx0 + snpIdx):(snpIdx0 + snpIdx + span - 1)]

        dhs.cur <- dhs[(snpIdx0 + snpIdx):(snpIdx0 + snpIdx + span - 1)]
        dhs.cur <- as.numeric(dhs.cur)

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

        ## for (J in 1:ncol(s0$vectors)) {
        ##     plot(1:nrow(Q0), Q0[, J], ty='h')
        ##     readline()
        ## }
        
        ## D1m <- matrix(0, nrow=ncol(s0$vectors), ncol=ncol(s0$vectors))
        ## diag(D1m) <- pmax(s0$values, 0)
        ## m <- Q0 %*% D1m %*% t(Q0)
        ## image(m)

        lambda0 <- s0$values[s0$values > LAMBDA.CO]
        
        Q0 <- Q0[, s0$values > LAMBDA.CO]

        # in contrast to manuscript, qX0 was not divided by sqrt(lambda0) here
        # thus, backconversion code to beta scale does not 
        qX0 <- X0 %*% Q0

        etahat0 <- t(Q0) %*% as.matrix(betahat.cur) / sqrt(lambda0)

        eta0 <- t(Q0) %*% as.matrix(beta.cur) * sqrt(lambda0)

        # dhs
        dhs0 <- apply(Q0, 2, function (x) sum(x**2 * dhs.cur))
        
        windata <- list()
    
        windata[["eigen"]] <- s0
        windata[["Q0.X"]] <- qX0
        windata[["etahat0"]] <- etahat0
        windata[["eta0"]] <- eta0
        windata[["dhs0"]] <- dhs0

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
            data.frame(lambda=lambda0, etahat=etahat0, eta=eta0, dhs=dhs0),
            file=paste(winfilepre, ".table", sep=''),
            quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t")

        saveRDS(Q0, paste(winfilepre, ".Q.RDS", sep=''))        

        etahat.all <- c(etahat.all, etahat0)
        eta.all <- c(eta.all, eta0)
        eval.all <- c(eval.all, lambda0)
        dhs.all <- c(dhs.all, dhs0)
        

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

        pos.cur <-
            summstat$pos[(snpIdx0 + snpIdx):(snpIdx0 + snpIdx + span - 1)]
        
        beta.cur <- beta[(snpIdx0 + snpIdx):(snpIdx0 + snpIdx + span - 1)]
    
        betahat.cur <-
            std.effalt[(snpIdx0 + snpIdx):(snpIdx0 + snpIdx + span - 1)]

        dhs.cur <- dhs[(snpIdx0 + snpIdx):(snpIdx0 + snpIdx + span - 1)]
        dhs.cur <- as.numeric(dhs.cur)

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

        eta0 <- t(Q0) %*% as.matrix(beta.cur) * sqrt(lambda0)

        # dhs
        dhs0 <- apply(Q0, 2, function (x) sum(x**2 * dhs.cur))
        
        windata <- list()
    
        windata[["eigen"]] <- s0
        windata[["Q0.X"]] <- qX0
        windata[["etahat0"]] <- etahat0
        windata[["eta0"]] <- eta0
        windata[["dhs0"]] <- dhs0

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
            data.frame(lambda=lambda0, etahat=etahat0, eta=eta0, dhs=dhs0),
            file=paste(winfilepre, ".table", sep=''),
            quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t")

        saveRDS(Q0, paste(winfilepre, ".Q.RDS", sep=''))        

        etahat.all <- c(etahat.all, etahat0)
        eta.all <- c(eta.all, eta0)
        eval.all <- c(eval.all, lambda0)
        dhs.all <- c(dhs.all, dhs0)
    }

    snpIdx0 <- snpIdx0 + M.chr

    close(stdgt.file)    
}

cat("Done\n")

q(save="no")

