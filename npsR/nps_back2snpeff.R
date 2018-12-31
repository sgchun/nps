VERSION <- "1.0.0"

cat("Non-Parametric Shrinkage", VERSION, "\n")

library(MASS)

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
if (WINSHIFT == 0) {
    part <- readRDS(paste(tempprefix, "part.RDS", sep=''))
} else {
    part <- readRDS(paste(tempprefix, "win_", WINSHIFT, ".part.RDS", sep=''))
}

Nt <- part[["Nt"]]
nLambdaPT <- part[["nLambdaPT"]]
nEtaPT <- part[["nEtaPT"]]

lambda.q <- part[["lambda.q"]]
betahatH.q <- part[["betahatH.q"]]


# Read summary stats (discovery)
summstat <- read.delim(summstatfile, header=TRUE, stringsAsFactors=FALSE,
                       sep="\t")
dim(summstat)

se <- sqrt(2 * summstat$reffreq * (1 - summstat$reffreq))

std.effalt <- summstat$effalt * se

summstat <- cbind(summstat, std.effalt=std.effalt)


# Use traing AF instead of discovery AF
trfrq <- read.table(trainfreqfile, header=TRUE)
tr.se <- sqrt(2 * trfrq$AAF * (1 - trfrq$AAF))
#plot(tr.se, se, cex=0.25)
#abline(0, 1, col="red")

tr.se.chr <- tr.se[summstat$chr == paste('chr', CHR, sep='')]

cat("M", "CHR", CHR, "=", length(tr.se.chr), "\n")

########################

if (WINSHIFT == 0) {
    PTwt <- readRDS(paste(tempprefix, "PTwt.RDS", sep=''))
} else {
    PTwt <- readRDS(paste(tempprefix, "win_", WINSHIFT, ".PTwt.RDS", sep=''))
}

PTwt.tail <- 0    

if (file.exists(paste(tempprefix, "PTwt.tail.RDS", sep=''))) {
    # Need handling for tail 

    if (WINSHIFT == 0) {
        PTwt.tail <- readRDS(paste(tempprefix, "PTwt.tail.RDS", sep=''))

        cat("Load PTwt.tail.RDS: ", PTwt.tail, "\n")
        
    } else {

        PTwt.file <- paste(tempprefix, "win_", WINSHIFT, ".PTwt.tail.RDS", sep='')
    
        if (file.exists(PTwt.file)) {
            
            PTwt.tail <- readRDS(PTwt.file)

            cat("Load ", PTwt.file, ":", PTwt.tail, "\n")

        } else {
            
            PTwt.tail <- readRDS(paste(tempprefix, "PTwt.tail.RDS", sep=''))

            cat("Load PTwt.tail.RDS: ", PTwt.tail, "\n")
        }
    }
}

# get weighted betahat back

wt.betahat <- c()

I <- 1

if (WINSHIFT == 0) {
    winfilepre <-
        paste(tempprefix, "win.", CHR, ".", I, sep='')
} else {
    winfilepre <-
        paste(tempprefix, "win_", WINSHIFT, ".", CHR, ".", I, sep='')
}

while (file.exists(paste(winfilepre, ".pruned", ".table", sep=''))) {

    print(I)
    
    wintab <- read.delim(paste(winfilepre, ".pruned", ".table", sep=''),
                         header=TRUE, sep="\t")

    tailfixfile <- paste(winfilepre, ".pruned", ".tailfix.table", sep='')
                             
    if (file.exists(tailfixfile)) {
        # override
        cat("Using window data residualized on GWAS-sig SNPs: ",
            tailfixfile, "\n")
            
        wintab <- read.delim(tailfixfile, header=TRUE, sep="\t")
    }
    
    lambda0 <- wintab$lambda
    etahat0 <- wintab$etahat
        
    Q0 <- readRDS(paste(winfilepre, ".Q.RDS", sep=''))

    etahat0 <- etahat0[lambda0 > 0]
    Q0 <- Q0[, lambda0 > 0]
    lambda0 <- lambda0[lambda0 > 0]

    Nq <- length(etahat0)

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

    if (WINSHIFT == 0) {
        winfilepre <-
            paste(tempprefix, "win.", CHR, ".", I, sep='')
    } else {
        winfilepre <-
            paste(tempprefix, "win_", WINSHIFT, ".", CHR, ".", I, sep='')
    }
    
}

# pad
M.chr <- sum(summstat$chr == paste('chr', CHR, sep=''))

M.written <- length(wt.betahat)

if ((M.chr - M.written) > 0) {
    cat("Pad ", (M.chr - M.written), " SNPs with 0 at the end of chrom\n")
    
    wt.betahat <- c(wt.betahat, rep(0, M.chr - M.written))
}

ASSERT(M.chr == length(wt.betahat))


# add tail betahats
betahat.tail.chr <-
    read.delim(paste(tempprefix, "tail_betahat.", CHR, ".table",
                     sep=''), header=FALSE, sep="\t")[, 1]

ASSERT(length(betahat.tail.chr) == M.chr)

wt.betahat <- wt.betahat + betahat.tail.chr * PTwt.tail 



# se: discovery af 
#    wt.betahat <- wt.betahat / se[snpIdx0 + c(1:M.chr)]

# se: training af

ASSERT(length(tr.se.chr) == M.chr)

wt.betahat <- wt.betahat / tr.se.chr


if (WINSHIFT == 0) {

    ## write.table(data.frame(betahat=wt.betahat),
    ##             file=paste(traindir, "/", traintag, ".adjbetahat.chrom",
    ##                 CHR, ".txt", sep=''),
    ##             quote=FALSE, row.names=FALSE, col.names=FALSE)

    filename <- paste(tempprefix, "/", traintag, ".adjbetahat.chrom",
                      CHR, ".txt", sep='')

    cat("Saving snpeffs:", filename, "...")

    write.table(data.frame(betahat=wt.betahat),
                file=filename,
                quote=FALSE, row.names=FALSE, col.names=FALSE)
    
} else {
    
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

}

cat("OK\n")

cat("Done\n")
