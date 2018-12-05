library(MASS)
library(pROC)

source("nps.1.2b/nps.R")

#########################################################################

# WINSZ <- 4000

#########################################################################

cargs <- commandArgs(trailingOnly=TRUE)

if (length(cargs) != 2 && length(cargs) != 3) {
    stop("Rscript run_nps5svd2.PTwt.R <outdir> <WINSHIFT> [ <plot.pdf> ]")
}

tempprefix <- paste(cargs[1], "/", sep='')

# Read in saved settings
args <- readRDS(paste(tempprefix, "args.RDS", sep=''))

WINSZ <- args[["WINSZ"]]
trainfamfile <- args[["trainfamfile"]]
trainphenofile <- args[["trainphenofile"]]


WINSHIFT <- as.numeric(cargs[2])

if (is.nan(WINSHIFT) || WINSHIFT < 0 || WINSHIFT >= WINSZ) {
    stop("Invalid shift:", cargs[2])
}

plotfile <- "plot.pdf"

if (length(cargs) == 3) {

    plotfile <- cargs[3]

    suffix <- substr(plotfile, nchar(plotfile) - 3, nchar(plotfile))

    if (suffix != ".pdf") {
        stop("Invalid pdf output file name:", cargs[3])
    }
}

#tempprefix <- 
#    "/net/home/schun/gwassim2/mvn2/train/nps5.23301.4000/"

# WINSHIFT <- 0
# WINSHIFT <- 1000
# WINSHIFT <- 2000
# WINSHIFT <- 3000

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
nDHSPT <- part[["nDHSPT"]]

lambda.q <- part[["lambda.q"]]
betahatH.q <- part[["betahatH.q"]]

#########################################################################

cat("trainfamfile:", trainfamfile, "\n")
cat("trainphenofile:", trainphenofile, "\n")

# phenotypes
trfam <- read.delim(trainfamfile, sep=" ", header=FALSE,
                    stringsAsFactors=FALSE)
trphen <- read.delim(trainphenofile, sep="\t", header=TRUE,
                     stringsAsFactors=FALSE)

rownames(trphen) <- paste(trphen$FID, trphen$IID, sep=":")
trphen <- trphen[paste(trfam[, 1], trfam[, 2], sep=":"), ]


print(length(intersect(paste(trfam[, 1], trfam[, 2], sep=":"),
                       paste(trphen$FID, trphen$IID, sep=":")
                       )))
print(sum(is.na(trphen$Outcome)))

ASSERT(all(!is.na(trphen$Outcome)))
ASSERT(all(trphen$FID == trfam[, 1]))
ASSERT(all(trphen$IID == trfam[, 2]))

trY <- trphen$Outcome

ASSERT(Nt == length(trY))

#########################################################################
# Read partitions

trPT <- array(0, dim=c(Nt, nLambdaPT, nEtaPT, nDHSPT))    

for (chrom in 1:22) {

    if (WINSHIFT == 0) {

        trPT.chr <- readRDS(paste(tempprefix, "trPT.", chrom, ".RDS", sep=''))
        
    } else {
        trPT.chr <-
            readRDS(paste(tempprefix, "win_", WINSHIFT, ".trPT.", chrom,
                          ".RDS", sep=''))
    }

    trPT <- trPT + trPT.chr
}


#########################################################################

PTwt <- array(0, dim=c(nLambdaPT, nEtaPT, nDHSPT))

for (I in 1:nLambdaPT) {
    for (J in 1:nEtaPT) {
        for (K in 1:nDHSPT) {
            
            trcaVAR <- var(trPT[trY == 1, I, J, K])
            trctVAR <- var(trPT[trY == 0, I, J, K])
            trptVAR <- (trcaVAR + trctVAR) / 2
        
            trcaMU <- mean(trPT[trY == 1, I, J, K])
            trctMU <- mean(trPT[trY == 0, I, J, K])
        
            PTwt[I, J, K] <- (trcaMU - trctMU) / trptVAR
        }
    }
}

print(sum(is.nan(PTwt)))

PTwt[is.nan(PTwt)] <- 0

print(PTwt[ , , 1])

if (WINSHIFT == 0) {
    saveRDS(PTwt, paste(tempprefix, "PTwt.RDS", sep=''))
} else {
    saveRDS(PTwt, paste(tempprefix, "win_", WINSHIFT, ".PTwt.RDS", sep=''))
}


#########################################################################

# tail partition
trPT.tail <- rep(0, Nt)

for (chrom in 1:22) {
    
    trPT.tail.file <-
            paste(tempprefix, "trPT.", chrom, ".tail.RDS", sep='')
            
    if (file.exists(trPT.tail.file)) {
            
        cat("chrom - tail - ", chrom, "\n")

        trPT.tail.chr <- readRDS(trPT.tail.file)
        
        trPT.tail <- trPT.tail + trPT.tail.chr
    }

}

PTwt.tail <- 0

if (any(trPT.tail != 0)) {
    
    trcaVAR <- var(trPT.tail[trY == 1])
    trctVAR <- var(trPT.tail[trY == 0])
    trptVAR <- (trcaVAR + trctVAR) / 2
        
    trcaMU <- mean(trPT.tail[trY == 1])
    trctMU <- mean(trPT.tail[trY == 0])
    
    PTwt.tail <- (trcaMU - trctMU) / trptVAR
}

print(PTwt.tail)

if (WINSHIFT == 0) {
    
    saveRDS(PTwt.tail, paste(tempprefix, "PTwt.tail.RDS", sep=''))
    
}


######################################################################
# Training R2

predY0 <- rep(0, Nt)

for (I in 1:nLambdaPT) {
    for (J in 1:nEtaPT) {
        for(K in 1:nDHSPT) {

            predY0 <- predY0 + PTwt[I, J, K] * trPT[, I, J, K]

        }
    }
}

predY0 <- predY0 + PTwt.tail * trPT.tail 

cat("R2obs =", cor(trY, predY0)**2, "\n")

cat("AUC :\n")
print(roc(cases=predY0[trY == 1], controls=predY0[trY == 0], ci=FALSE))


########################################################################
# For diagnostics

pdf(file=plotfile)

load(paste(tempprefix, "run_nps5svd2.", "win_", WINSHIFT, ".RData", sep=''))

for (Ix in nLambdaPT:1) {
    bh.range <- range(meanBetahatH[, , 1] * PTwt[, ,1])
    bh.sign <- sign(bh.range)[which.max(abs(bh.range))]
    
    plot(c(), c(), 
         xlim=c(0, max(meanBetahatH[, , 1]) * 1.2),
         ylim=c(0, bh.sign * max(abs(meanBetahatH[, , 1] * PTwt[, ,1])) * 1.2),
         main=paste("Lambda: ", Ix, "-th decile"),
         xlab="Raw eta hat", ylab="Re-scaled eta hat")
    
    for (Jx in 1:nEtaPT) {

        points(betahatH.q[Jx:(Jx + 1), Ix],
               betahatH.q[Jx:(Jx + 1), Ix] * PTwt[Ix, Jx, 1], ty='l')
        points(meanBetahatH[Ix, Jx, 1], meanBetahatH[Ix, Jx, 1] * PTwt[Ix, Jx, 1])
        
    }
        
#    plot(meanBetahatH[Ix, , 1], meanBetahatH[Ix, , 1] * PTwt[Ix, ,1], ty='o',
#         xlim=c(0, max(meanBetahatH[Ix, , 1])),
#         ylim=c(0, max(meanBetahatH[Ix, , 1] * PTwt[Ix, ,1])))

##    readLines(n=1) # for interactive mode
} 


## for (K in nLambdaPT:1) {

##     trPT.cor <- matrix(0, nrow=nEtaPT, ncol=nEtaPT)

##     for (I in 1:nEtaPT) {
##         for (J in 1:nEtaPT) {
            
##             trPT.cor[I, J] <- cor(trPT[, K, I, 1], trPT[, K, J, 1])
##         }
##     }

##     image(trPT.cor)

## ##    readLines(n=1) # for interactive mode
## }

dev.off()
