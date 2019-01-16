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

if (length(cargs) < 2) {
    stop("Usage: Rscript nps_train_AUC.R <work dir> [ <WINSHIFT> ]+")
}

tempprefix <- paste(cargs[1], "/", sep='')

# Read in saved settings
args <- readRDS(paste(tempprefix, "args.RDS", sep=''))

WINSZ <- args[["WINSZ"]]
trainfamfile <- args[["trainfamfile"]]
trainphenofile <- args[["trainphenofile"]]


WINSHIFT.list <- as.numeric(cargs[2:length(cargs)])

if (any(is.nan(WINSHIFT.list)) || any(WINSHIFT.list < 0) ||
    any(WINSHIFT.list >= WINSZ)) {
    stop("Invalid shift:", cargs[2:length(cargs)])
}

#########################################################################

cat("trainfamfile:", trainfamfile, "\n")
cat("trainphenofile:", trainphenofile, "\n")

# phenotypes
trfam <- read.delim(trainfamfile, sep=" ", header=FALSE,
                    stringsAsFactors=FALSE)
trphen <- read.delim(trainphenofile, sep="\t", header=TRUE,
                     stringsAsFactors=FALSE)

if (ncol(trfam) != 6) {
    # re-try with tab delimination

    trfam <- read.delim(trainfamfile, sep="\t", header=FALSE,
                        stringsAsFactors=FALSE)
}

ASSERT(ncol(trfam) == 6)

rownames(trphen) <- paste(trphen$FID, trphen$IID, sep=":")
trphen <- trphen[paste(trfam[, 1], trfam[, 2], sep=":"), ]


#print(length(intersect(paste(trfam[, 1], trfam[, 2], sep=":"),
#                       paste(trphen$FID, trphen$IID, sep=":")
#                       )))
#print(sum(is.na(trphen$Outcome)))

ASSERT(all(!is.na(trphen$Outcome)))
ASSERT(all(trphen$FID == trfam[, 1]))
ASSERT(all(trphen$IID == trfam[, 2]))

trY <- trphen$Outcome

#########################################################################
# Read partitions

predY0 <- rep(0, length(trY))

# tail partition
trPT.tail <- rep(0, length(trY))

for (chrom in 1:22) {
    
    trPT.tail.file <-
        paste(tempprefix, "trPT.", chrom, ".tail.RDS", sep='')
    
    if (file.exists(trPT.tail.file)) {

        cat("Loading S0 partition for chrom", chrom, "...\n")

        trPT.tail.chr <- readRDS(trPT.tail.file)
        
        trPT.tail <- trPT.tail + trPT.tail.chr
    }

}

PTwt.tail <- readRDS(paste(tempprefix, "PTwt.tail.RDS", sep=''))

for (WINSHIFT in WINSHIFT.list) {

    cat("winshift =", WINSHIFT, "...\n")

    if (WINSHIFT == 0) {
        part <- readRDS(paste(tempprefix, "part.RDS", sep=''))
    } else {
        part <-
            readRDS(paste(tempprefix, "win_", WINSHIFT, ".part.RDS", sep=''))
    }

    Nt <- part[["Nt"]]
    nLambdaPT <- part[["nLambdaPT"]]
    nEtaPT <- part[["nEtaPT"]]

    ASSERT(Nt == length(trY))

    trPT <- array(0, dim=c(Nt, nLambdaPT, nEtaPT, 1))    

    for (chrom in 1:22) {

        if (WINSHIFT == 0) {

            trPT.chr <-
                readRDS(paste(tempprefix, "trPT.", chrom, ".RDS", sep=''))
        
        } else {
            trPT.chr <-
                readRDS(paste(tempprefix, "win_", WINSHIFT, ".trPT.", chrom,
                              ".RDS", sep=''))
        }

        trPT <- trPT + trPT.chr
    }


    if (WINSHIFT == 0) {
        PTwt <- readRDS(paste(tempprefix, "PTwt.RDS", sep=''))
    } else {
        PTwt <-
            readRDS(paste(tempprefix, "win_", WINSHIFT, ".PTwt.RDS", sep=''))
    }


    for (I in 1:nLambdaPT) {
        for (J in 1:nEtaPT) {
            K <- 1
            
            predY0 <- predY0 + PTwt[I, J, K] * trPT[, I, J, K]
        }
    }

    predY0 <- predY0 + PTwt.tail * trPT.tail 

}

library(pROC)

cat("AUC :\n")
print(roc(cases=predY0[trY == 1], controls=predY0[trY == 0], ci=TRUE))
