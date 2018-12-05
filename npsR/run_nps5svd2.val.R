library(MASS)
library(pROC)
library(DescTools)
# library(boot)

source("nps.1.2b/nps.R")


#########################################################################
# MVN 2

cargs <- commandArgs(trailingOnly=TRUE)

if (length(cargs) < 4) {
    stop("Rscript run_nps5svd2.R <outdir> <valphenofile> <K> [<WINSHIFT>]+...")
}

tempprefix <- paste(cargs[1], "/", sep='')

args <- readRDS(paste(tempprefix, "args.RDS", sep=''))

traintag <- args[["traintag"]]
valdir <- args[["valdir"]] 
valfamfile <- args[["valfamfile"]]
WINSZ <- args[["WINSZ"]]

valphenofile <- cargs[2]

K <- as.numeric(cargs[3])

if (is.na(K) || is.nan(K) || K >= 1 || K < 0) {
    stop("Invalid prevalence:", cargs[3])
}

list.WINSHIFT <- c() 

for (carg.WSHIFT in cargs[4:length(cargs)]) {

    WSHIFT <- as.numeric(carg.WSHIFT)

    if (is.nan(WSHIFT) || WSHIFT < 0 || WSHIFT > WINSZ) {
        stop(paste("Invalid shift:", carg.WSHIFT)) 
    }

    list.WINSHIFT <- c(list.WINSHIFT, WSHIFT)
}

#########################################################################
### validation

# phenotypes
vlfam <- read.delim(valfamfile, sep=" ", header=FALSE,
                    stringsAsFactors=FALSE)
vlphen <- read.delim(valphenofile, sep="\t", header=TRUE,
                     stringsAsFactors=FALSE)

rownames(vlphen) <- paste(vlphen$FID, vlphen$IID, sep=":")
vlphen <- vlphen[paste(vlfam[, 1], vlfam[, 2], sep=":"), ]

ASSERT(all(vlphen$FID == vlfam[, 1]))
ASSERT(all(vlphen$IID == vlfam[, 2]))
ASSERT(all(!is.na(vlphen$Outcome)))

# genetic risks
for (WINSHIFT in list.WINSHIFT) {

    vlY <- vlphen$Outcome

    prisk <- rep(0, length(vlY))

    for (chr in 1:22) {

        cat("chr", chr, "\n")

        # read per-chrom genetic risk file
        if (WINSHIFT == 0) {
            prisk.file <-
                paste(valdir, "/", traintag, ".predY.chrom", chr, ".txt",
                      sep='')
        } else {
            prisk.file <-
                paste(valdir, "/", traintag, ".win_", WINSHIFT,
                      ".predY.chrom", chr, ".txt", sep='')
        }
    
        prisk.chr <- read.delim(prisk.file, header=FALSE, sep="\t")[, 1]
    
        ASSERT(length(prisk.chr) == length(vlY))
    
        prisk <- prisk + prisk.chr
    
    }

    prisk <- prisk[vlY >= 0] 
    vlL <- vlphen$TotalLiability[vlY >= 0]	
    vlY <- vlY[vlY >= 0]
    
    cat("WINSHIFT =", WINSHIFT, "\n")
    
    # R2 observed scale
    cat("R2obs =", cor(vlY, prisk)**2, "\n")

    # R2 Nagelkerke
    mod <- glm(vlY ~ prisk, family=binomial(link="logit"))
    cat("R2nag =", PseudoR2(mod, "Nagelkerke"), "\n")

    # R2 liability scale
    if ("TotalLiability" %in% colnames(vlphen)) {
        cat("R2liability =", cor(vlL, prisk)**2, "\n")
    }

    # AUC
    cat("AUC:\n")
    print(roc(cases=prisk[vlY == 1], controls=prisk[vlY == 0], ci=FALSE))

    cat("------------------------------------------------\n")
}


# Average
vlY <- vlphen$Outcome

prisk <- rep(0, length(vlY))    

for (WINSHIFTx in list.WINSHIFT) {

    vlY <- vlphen$Outcome

    for (chr in 1:22) {

        # read per-chrom genetic risk file
        if (WINSHIFTx == 0) {
            prisk.file <-
                paste(valdir, "/", traintag, ".predY.chrom", chr, ".txt",
                      sep='')
        } else {
            prisk.file <-
                paste(valdir, "/", traintag, ".win_", WINSHIFTx,
                      ".predY.chrom", chr, ".txt", sep='')
        }

        prisk.chr <- read.delim(prisk.file, header=FALSE, sep="\t")[, 1]

        ASSERT(length(prisk.chr) == length(vlY))
        
        prisk <- prisk + prisk.chr
            
    }
}

cat("FINAL: averaged over WINSHIFTs\n")

prisk <- prisk[vlY >= 0] 
vlL <- vlphen$TotalLiability[vlY >= 0]	
vlY <- vlY[vlY >= 0]

cat("R2obs =", cor(vlY, prisk)**2, "\n")

mod <- glm(vlY ~ prisk, family=binomial(link="logit"))
cat("R2nag =", PseudoR2(mod, "Nagelkerke"), "\n")

if ("TotalLiability" %in% colnames(vlphen)) {
    cat("R2liability =", cor(vlL, prisk)**2, "\n")
}

cat("AUC:\n")
#print(roc(cases=prisk[vlY == 1], controls=prisk[vlY == 0], ci=FALSE))
print(roc(cases=prisk[vlY == 1], controls=prisk[vlY == 0], ci=TRUE))

df.out <- cbind(vlphen, Score=prisk)
write.table(df.out, file=paste(valphenofile, ".nps_score", sep=''),
            row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)

if (cor(vlY, prisk) < 0) {
    cat("Sign flipped: ", cor(vlY, prisk), "\n")

    prisk <- - prisk
}

cat("Tail OR:\n")

# estimate cutoff score 
cutoff.tab <- list()

if (K > 0) {
    cat("Prevalence parameter = ", K, "\n")
    cat("Sampling.")

    expected.NControls <- round( (1 - K) / K * sum(vlY == 1) )
    prisk.cases <- prisk[vlY == 1]
    prisk.controls <- prisk[vlY == 0]

    dist.0.005 <- c()
    dist.0.01 <- c()
    dist.0.05 <- c()

    idx.0.005 <- round( (sum(vlY == 1) + expected.NControls) * 0.005 )
    idx.0.01 <- round( (sum(vlY == 1) + expected.NControls) * 0.01 )
    idx.0.05 <- round( (sum(vlY == 1) + expected.NControls) * 0.05 )

    ASSERT(idx.0.005 >= 1)
    
    for (I in 1:10000) {
        if ((I %% 100) == 0) {
            cat(".")
        }
    
        # Randomly up-sample 
#        sampled.prisk <-
#            c(prisk.cases,
#              sample(prisk.controls, expected.NControls, replace=TRUE))
        sampled.prisk <-
            c(sample(prisk.cases, length(prisk.cases), replace=TRUE), 
              sample(prisk.controls, expected.NControls, replace=TRUE))

        sampled.prisk <-  sort(sampled.prisk, decreasing=TRUE)

        dist.0.005 <-
            c(dist.0.005, sampled.prisk[idx.0.005])

        dist.0.01 <-
            c(dist.0.01, sampled.prisk[idx.0.01])

        dist.0.05 <-
            c(dist.0.05, sampled.prisk[idx.0.05])
    }

    cat("\n")

    cat("0.005:", sd(dist.0.005)/sqrt(length(dist.0.005)), ", ",
        median(dist.0.005), "\n")
    cat("0.01:", sd(dist.0.01)/sqrt(length(dist.0.005)), ", ",
        median(dist.0.01), "\n")
    cat("0.05:", sd(dist.0.05)/sqrt(length(dist.0.005)), ", ",
        median(dist.0.05), "\n")

    cutoff.tab[["0.005"]] <- median(dist.0.005)
    cutoff.tab[["0.01"]] <- median(dist.0.01)
    cutoff.tab[["0.05"]] <- median(dist.0.05)    
}

btOR <- function(data, ind, cutoff) {
    d <- data[ind, ]
    
    odds1 <- sum(d$prisk >= cutoff & d$vlY == 1) /
        sum(d$prisk >= cutoff & d$vlY == 0)
    odds0 <- sum(d$prisk < cutoff & d$vlY == 1) /
        sum(d$prisk < cutoff & d$vlY == 0)

    odds1 / odds0
}

for (cutoff in c("0.01", "0.05")) {

    cutoff.prisk <- NA
    
    if (K == 0) {
        # epidemiological cohort
        cutoff.at <- round(length(prisk) * as.numeric(cutoff))
        cutoff.prisk <- sort(prisk, decreasing=TRUE)[cutoff.at]
    } else {
        
        # estimate cutoff score 

#    cutoff.at <- round(sum(vlY == 0) * cutoff)
#    cutoff.prisk <- sort(prisk[vlY == 0], decreasing=TRUE)[cutoff.at]

        cutoff.prisk <- cutoff.tab[[cutoff]]
        
    }
    
    odds1 <- sum(prisk >= cutoff.prisk & vlY == 1) /
        sum(prisk >= cutoff.prisk & vlY == 0)
    odds0 <- sum(prisk < cutoff.prisk & vlY == 1) /
        sum(prisk < cutoff.prisk & vlY == 0)

    cat(cutoff, ", Ntail=", sum(prisk >= cutoff.prisk), ", ",
        (odds1 / odds0), "\n")


#    if (TRUE) {
#
#        bt.iter <- ceiling(length(vlY) / 10000) * 10000
#
#        btdist <-
#            boot(data=data.frame(prisk=prisk, vlY=vlY), statistic=btOR, R=bt.iter,
#                 cutoff=cutoff.prisk)
#
#        print(btdist$t0)
#        print(boot.ci(btdist, type="bca"))
#    }
}
