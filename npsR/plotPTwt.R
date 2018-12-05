#library(MASS)
#library(pROC)
#source("nps.1.2b/nps.R")

cargs <- commandArgs(trailingOnly=TRUE)

if (length(cargs) < 3) {
    stop("Rscript plotPTwt.R <outdir> <plot.pdf> [<WINSHIFT>]+")
}

tempprefix <- paste(cargs[1], "/", sep='')

# Read in saved settings
args <- readRDS(paste(tempprefix, "args.RDS", sep=''))

WINSZ <- args[["WINSZ"]]


plotfile <- cargs[2]


list.WINSHIFT <- c() 

for (carg.WSHIFT in cargs[3:length(cargs)]) {

    WSHIFT <- as.numeric(carg.WSHIFT)

    if (is.nan(WSHIFT) || WSHIFT < 0 || WSHIFT > WINSZ) {
        stop(paste("Invalid shift:", carg.WSHIFT)) 
    }

    list.WINSHIFT <- c(list.WINSHIFT, WSHIFT)
}

cat("Write to: ", plotfile, "\n")
cat("Shifts by: ", paste(list.WINSHIFT, collapse=", "), "\n")

#####

#tempprefix <- "~/ukbb.t2d/npsdat/"
#WINSZ <- 4000
#list.WINSHIFT <- c(0, 1000, 2000, 3000)

WINSHIFT <- list.WINSHIFT[1]

if (WINSHIFT == 0) {
    part <- readRDS(paste(tempprefix, "part.RDS", sep=''))
} else {
    part <- readRDS(paste(tempprefix, "win_", WINSHIFT, ".part.RDS", sep=''))
}

nLambdaPT <- part[["nLambdaPT"]]
nEtaPT <- part[["nEtaPT"]]


#########################################################################

pdf(file=plotfile, useKerning=FALSE, useDingbats=FALSE)

x.end.Ix <- rep(0, nLambdaPT)
x.end.Ix.cap <- rep(NA, nLambdaPT)
y.end <- 0

for (WINSHIFT in list.WINSHIFT) {

    load(paste(tempprefix, "run_nps5svd2.", "win_", WINSHIFT, ".RData",
               sep=''))

    if (WINSHIFT == 0) {
        PTwt <- readRDS(paste(tempprefix, "PTwt.RDS", sep=''))
    } else {
        PTwt <- readRDS(paste(tempprefix, "win_", WINSHIFT, ".PTwt.RDS",
                              sep=''))
    }

    x.end.Ix <- pmax(x.end.Ix, meanBetahatH[, nEtaPT, 1] * 1.2)
    x.end.Ix.cap <- pmin(x.end.Ix.cap, betahatH.q[nEtaPT + 1, ], na.rm=TRUE)
    y.end <- max(y.end, max(abs(meanBetahatH[, , 1] * PTwt[, , 1])) * 1.2)

}

x.end <- max(x.end.Ix)

x.end.Ix <- x.end.Ix.cap



scale <- 1

PTwt.tail.file <- paste(tempprefix, "PTwt.tail.RDS", sep='')

if (file.exists(PTwt.tail.file)) {
    PTwt.tail <- readRDS(PTwt.tail.file)

    if (PTwt.tail > 0) {
        scale <- 1 / PTwt.tail
    }

    cat("Rescale y-axis by ", scale, "\n")
}

Ix <- nLambdaPT

for (Ix in nLambdaPT:1) {

    x.seq <- seq(0, min(x.end.Ix[Ix], x.end), min(x.end.Ix[Ix], x.end) / 100)
    y.seq <- rep(0, length(x.seq))

    for (WINSHIFT in list.WINSHIFT) {

        load(paste(tempprefix, "run_nps5svd2.", "win_", WINSHIFT, ".RData",
                   sep=''))
        
        if (WINSHIFT == 0) {
            PTwt <- readRDS(paste(tempprefix, "PTwt.RDS", sep=''))
        } else {
            PTwt <- readRDS(paste(tempprefix, "win_", WINSHIFT, ".PTwt.RDS",
                                  sep=''))
        }

        for (Jx in 1:nEtaPT) {
            Jx.st <- betahatH.q[Jx, Ix]
            Jx.ed <- betahatH.q[Jx + 1, Ix]
            wt <- PTwt[Ix, Jx, 1]

            seg <- x.seq >= Jx.st & x.seq < Jx.ed

            y.seq[seg] <- y.seq[seg] + wt * x.seq[seg]
        }
    }
    
    y.seq <- y.seq / 4 * scale

    # ignore the end point
    x.seq <- x.seq[1:(length(x.seq) - 1)]
    y.seq <- y.seq[1:(length(y.seq) - 1)]

    plot(x.seq, y.seq, ty='l',
         xlim=c(0, x.end),
         ylim=c(0, y.end * scale),
         main=paste("Lambda: ", Ix, "-th decile", sep=''),
         xlab="Raw eta hat", ylab="Rescaled eta hat")

}

dev.off()
