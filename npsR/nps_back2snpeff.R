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
    stop("Usage: Rscript nps_back2snpeff.R <work dir> <chrom> <WINSHIFT>")
}

tempprefix <- paste(cargs[1], "/", sep='')

args <- readRDS(paste(tempprefix, "args.RDS", sep=''))

summstatfile <- args[["summstatfile"]] 
traindir <- args[["traindir"]]
trainfreqfile <- args[["trainfreqfile"]]
traintag <- args[["traintag"]]
WINSZ <- args[["WINSZ"]]

CHR <- as.numeric(cargs[2])

if (!(CHR %in% 1:22)) {
    stop("invalid chrom", CHR)
}

WINSHIFT <- as.numeric(cargs[3])

if (is.nan(WINSHIFT) || WINSHIFT < 0 || WINSHIFT >= WINSZ) {
    stop("Invalid shift:", cargs[3])
}


#########################################################################

# Load partition data 
part <- readRDS(paste(tempprefix, "win_", WINSHIFT, ".part.RDS", sep=''))

Nt <- part[["Nt"]]
nLambdaPT <- part[["nLambdaPT"]]
nEtaPT <- part[["nEtaPT"]]

lambda.q <- part[["lambda.q"]]
betahatH.q <- part[["betahatH.q"]]


# Read summary stats (discovery)
summstat.chr <- read.delim(paste(summstatfile, ".", CHR, sep=''),
                           header=TRUE, stringsAsFactors=FALSE,
                           sep="\t")
#dim(summstat)

# Use traing AF instead of discovery AF
trfrq.chr <- read.table(paste(trainfreqfile, ".", CHR, sep=''), header=TRUE)
tr.se.chr <- sqrt(2 * trfrq.chr$AAF * (1 - trfrq.chr$AAF))
#plot(tr.se, se, cex=0.25)
#abline(0, 1, col="red")

M.chr <- length(tr.se.chr)

ASSERT(M.chr == nrow(summstat.chr))

cat("M", "CHR", CHR, "=", M.chr, "\n")

########################

PTwt <- readRDS(paste(tempprefix, "win_", WINSHIFT, ".PTwt.RDS", sep=''))

PTwt.tail <- 0    

if (file.exists(paste(tempprefix, "PTwt.tail.RDS", sep=''))) {
    # Need handling for tail 

    PTwt.tail <- readRDS(paste(tempprefix, "PTwt.tail.RDS", sep=''))

    cat("Load PTwt.tail.RDS: ", PTwt.tail, "\n")
    
}

# get weighted betahat back

wt.betahat <- c()

I <- 1

winfilepre <-
    paste(tempprefix, "win_", WINSHIFT, ".", CHR, ".", I, sep='')

while (file.exists(paste(winfilepre, ".pruned", ".table", sep=''))) {

    print(I)
    
    tailfixfile <- paste(winfilepre, ".pruned", ".table", sep='')
                             
    wintab <- read.delim(tailfixfile, header=TRUE, sep="\t")
    
    lambda0 <- wintab$lambda
    etahat0 <- wintab$etahat
        
    Q0 <- readRDS(paste(winfilepre, ".Q.RDS", sep=''))

    etahat0 <- etahat0[lambda0 > 0]
    Q0 <- Q0[, lambda0 > 0, drop=FALSE]
    lambda0 <- lambda0[lambda0 > 0]

    ## FIXME
    etahat0 <- etahat0[lambda0 > 10]
    Q0 <- Q0[, lambda0 > 10, drop=FALSE]
    lambda0 <- lambda0[lambda0 > 10]
    # 

    Nq <- length(etahat0)

    if (Nq == 0) {
        ## No projection left
        ## move on to next iteration

        I <- I + 1

        winfilepre <-
            paste(tempprefix, "win_", WINSHIFT, ".", CHR, ".", I, sep='')

        next
    }

    wt0 <- rep(NA, Nq)
    
    for (Il in 1:nLambdaPT) {
        
        lambda.lo <- lambda.q[Il]
        lambda.hi <- lambda.q[Il+1]
        in.lambda.bin <- lambda0 > lambda.lo & lambda0 <= lambda.hi

        for (Je in 1:nEtaPT) {

            betahatH.lo <- betahatH.q[Je, Il]
            betahatH.hi <- betahatH.q[Je+1, Il]
            in.betahatH.bin <-
                (in.lambda.bin & 
                 abs(etahat0) > betahatH.lo & abs(etahat0) <= betahatH.hi)

            if (any(in.betahatH.bin)) {
                wt0[in.betahatH.bin] <- PTwt[Il, Je, 1]
            }
        }
    }

    if (any(etahat0 == 0)) {
        wt0[etahat0 == 0] <- 0
    }

    ASSERT(all(!is.na(wt0)))

#   Compared to manuscript, we did not scale qX0 with lambda^(-1/2), 
#   thus no need to scale here again, wt0 includes the factor already.
#    etahat0.adj <- etahat0 * wt0 / sqrt(lambda0)
    etahat0.adj <- etahat0 * wt0 

    wt.betahat <- c(wt.betahat, Q0 %*% as.matrix(etahat0.adj))

    ASSERT(all(!is.na(wt.betahat)))

    # move on to next iteration
    I <- I + 1

    winfilepre <-
        paste(tempprefix, "win_", WINSHIFT, ".", CHR, ".", I, sep='')
}

# pad

M.written <- length(wt.betahat)

if ((M.chr - M.written) > 0) {
    cat("Pad ", (M.chr - M.written), " SNPs with 0 at the end of chrom\n")
    
    wt.betahat <- c(wt.betahat, rep(0, M.chr - M.written))
}

ASSERT(M.chr == length(wt.betahat))


# add tail betahats
tailbetahatfile <- paste(tempprefix, "tail_betahat.", CHR, ".table",
                         sep='')

if (file.exists(tailbetahatfile)) {
    
    betahat.tail.chr <-
        read.delim(tailbetahatfile, header=FALSE, sep="\t")[, 1]

    ASSERT(length(betahat.tail.chr) == M.chr)

    wt.betahat <- wt.betahat + betahat.tail.chr * PTwt.tail

}

# se: discovery af 
#    wt.betahat <- wt.betahat / se[snpIdx0 + c(1:M.chr)]

# se: training af

ASSERT(length(tr.se.chr) == M.chr)

wt.betahat <- wt.betahat / tr.se.chr

    ## write.table(data.frame(betahat=wt.betahat),
    ##             file=paste(traindir, "/", traintag, ".win_", WINSHIFT,
    ##                 ".adjbetahat.chrom", CHR, ".txt", sep=''),
    ##             quote=FALSE, row.names=FALSE, col.names=FALSE)

    filename <- paste(tempprefix, "/", traintag, ".win_", WINSHIFT,
                      ".adjbetahat.chrom", CHR, ".txt", sep='')

    cat("Saving snpeffs:", filename, "...")
    
    write.table(data.frame(betahat=wt.betahat),
                file=filename,
                quote=FALSE, row.names=FALSE, col.names=FALSE)


cat("OK\n")

cat("Done\n")
