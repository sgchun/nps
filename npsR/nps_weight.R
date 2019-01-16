VERSION <- "1.0.1"

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

if (length(cargs) != 2) {    
    stop("Usage: Rscript nps_weight.R <work dir> <WINSHIFT>")
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

#########################################################################

cat("train fam file:", trainfamfile, "\n")
cat("train pheno file:", trainphenofile, "\n")
cat("Size of training cohort:", Nt, "\n")

# phenotypes
trfam <- read.delim(trainfamfile, sep=" ", header=FALSE,
                    stringsAsFactors=FALSE)
trphen <- read.delim(trainphenofile, sep="\t", header=TRUE,
                     stringsAsFactors=FALSE)

rownames(trphen) <- paste(trphen$FID, trphen$IID, sep=":")
trphen <- trphen[paste(trfam[, 1], trfam[, 2], sep=":"), ]


# print(length(intersect(paste(trfam[, 1], trfam[, 2], sep=":"),
#                        paste(trphen$FID, trphen$IID, sep=":")
#                        )))
# print(sum(is.na(trphen$Outcome)))

ASSERT(all(!is.na(trphen$Outcome)))
ASSERT(all(trphen$FID == trfam[, 1]))
ASSERT(all(trphen$IID == trfam[, 2]))

trY <- trphen$Outcome

ASSERT(Nt == length(trY))

if (any(trY == -9)) {
    stop("Missing outcome (-9) is not allowed")
}

#########################################################################
# Read partitions

trPT <- array(0, dim=c(Nt, nLambdaPT, nEtaPT, 1))    

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

PTwt <- array(0, dim=c(nLambdaPT, nEtaPT, 1))

for (I in 1:nLambdaPT) {
    for (J in 1:nEtaPT) {
        K <- 1
            
        trcaVAR <- var(trPT[trY == 1, I, J, K])
        trctVAR <- var(trPT[trY == 0, I, J, K])
        trptVAR <- (trcaVAR + trctVAR) / 2
        
        trcaMU <- mean(trPT[trY == 1, I, J, K])
        trctMU <- mean(trPT[trY == 0, I, J, K])
        
        PTwt[I, J, K] <- (trcaMU - trctMU) / trptVAR
    }
}

if (any(is.nan(PTwt))) {
    cat("WARNING: ", sum(is.nan(PTwt)), "partitions produced NaN\n")
}

# cat(PTwt[ , , 1])

PTwt[is.nan(PTwt)] <- 0

cat("Saving ", nLambdaPT, "x", nEtaPT, "partition weights...")

if (WINSHIFT == 0) {
    saveRDS(PTwt, paste(tempprefix, "PTwt.RDS", sep=''))
} else {
    saveRDS(PTwt, paste(tempprefix, "win_", WINSHIFT, ".PTwt.RDS", sep=''))
}

cat("OK\n")

#########################################################################

# tail partition
trPT.tail <- rep(0, Nt)

for (chrom in 1:22) {
    
    trPT.tail.file <-
        paste(tempprefix, "trPT.", chrom, ".tail.RDS", sep='')
    
    if (file.exists(trPT.tail.file)) {
        
        cat("Loading S0 partition for chrom", chrom, "...\n")

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

#    cat("Weight for S0 =", PTwt.tail, "\n")

if (WINSHIFT == 0) {

    cat("Saving S0 weight...")
    
    saveRDS(PTwt.tail, paste(tempprefix, "PTwt.tail.RDS", sep=''))

    cat("OK\n")
} 

######################################################################
# Training R2

predY0 <- rep(0, Nt)

for (I in 1:nLambdaPT) {
    for (J in 1:nEtaPT) {
        K <- 1

        predY0 <- predY0 + PTwt[I, J, K] * trPT[, I, J, K]

    }
}

predY0 <- predY0 + PTwt.tail * trPT.tail 

cat("Observed scale R2 in training =", cor(trY, predY0)**2, "\n")
cat("Done\n")

