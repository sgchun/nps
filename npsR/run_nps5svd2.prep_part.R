library(MASS)

source("nps.1.2b/nps.R")

#########################################################################

# WINSZ <- 4000

#########################################################################

cargs <- commandArgs(trailingOnly=TRUE)

if (length(cargs) != 4 && length(cargs) != 5) {
    stop("Rscript run_nps5svd2.prep_part.R <outdir> <WINSHIFT> <nLambdaPT> <nEtaPT> [ <plot.pdf> ]")
}

tempprefix <- paste(cargs[1], "/", sep='')

args <- readRDS(paste(tempprefix, "args.RDS", sep=''))

Nt <- args[["Nt"]]
WINSZ <- args[["WINSZ"]]


WINSHIFT <- as.numeric(cargs[2])

if (is.nan(WINSHIFT) || WINSHIFT < 0 || WINSHIFT >= WINSZ) {
    stop("Invalid shift:", cargs[2])
}

nLambdaPT <- as.numeric(cargs[3])
nEtaPT <- as.numeric(cargs[4])
nDHSPT <- 1

if (is.nan(nLambdaPT) || nLambdaPT < 1) {
    stop("Invalid nLambdaPT:", cargs[3])
}

if (is.nan(nEtaPT) || nEtaPT < 1) {
    stop("Invalid nEtaPT:", cargs[4])
}

plotfile <- "plot.pdf"
if (length(cargs) == 5) {
    plotfile <- cargs[5]

    suffix <- substr(plotfile, nchar(plotfile) - 3, nchar(plotfile))

    if (suffix != ".pdf") {
        stop("Invalid pdf output file name:", cargs[5])
    }
}

#tempprefix <- 
#    "/net/home/schun/gwassim2/mvn2/train/nps5.23301.4000/"


# WINSHIFT <- 0
# WINSHIFT <- 1000
# WINSHIFT <- 2000
# WINSHIFT <- 3000

# nLambdaPT <- 10
# nEtaPT <- 10
# nDHSPT <- 1

# Nt <- 5000

#############################################################################
# Read back

etahat.all <- c()
eta.all <- c()
eval.all <- c()                         # lambda
dhs.all <- c()

chrom <- 1

for (chrom in 1:22) {
    print(chrom)

    I <- 1

    if (WINSHIFT == 0) {
        winfilepre <-
            paste(tempprefix, "win.", chrom, ".", I, sep='')
    } else {
        winfilepre <-
            paste(tempprefix, "win_", WINSHIFT, ".", chrom, ".", I,
                  sep='')
    }
    
    while (file.exists(paste(winfilepre, ".pruned", ".table", sep=''))) {

        wintab <- read.delim(paste(winfilepre, ".pruned", ".table", sep=''),
                             header=TRUE, sep="\t")

        tailfixfile <- paste(winfilepre, ".pruned", ".tailfix.table", sep='')
                             
        if (file.exists(tailfixfile)) {
            # override
            cat("Overriding: ", tailfixfile, "\n")
            
            wintab <- read.delim(tailfixfile, header=TRUE, sep="\t")
        }


        lambda0 <- wintab$lambda
        etahat0 <- wintab$etahat
        eta0 <- wintab$eta
        dhs0 <- wintab$dhs

        etahat0 <- etahat0[lambda0 > 0]
        eta0 <- eta0[lambda0 > 0]
        dhs0 <- dhs0[lambda0 > 0]
        lambda0 <- lambda0[lambda0 > 0]

#        print(length(lambda0))
        
        etahat.all <- c(etahat.all, etahat0)
        eta.all <- c(eta.all, eta0)
        eval.all <- c(eval.all, lambda0)
        dhs.all <- c(dhs.all, dhs0)
        
        # move on to next iteration
        I <- I + 1

        if (WINSHIFT == 0) {
            winfilepre <-
                paste(tempprefix, "win.", chrom, ".", I, sep='')
        } else {
            winfilepre <-
                paste(tempprefix, "win_", WINSHIFT, ".", chrom, ".", I, sep='')
        }
    }
}

print(length(etahat.all))

###########################################################################

## #plot(eta.all, etahat.all, cex=0.25)
## #cor(eta.all, etahat.all)

## #plot(eval.all, (etahat.all - eta.all)**2, cex=0.25)


## lambda.q <- w.quant.o2(sqrt(eval.all), nLambdaPT)
## lambda.q <- lambda.q ** 2
## lambda.q[1] <- 0
## lambda.q[nLambdaPT + 1] <- lambda.q[nLambdaPT + 1] * 1.1

## for (I in 10:1) {
##     in.lambda.bin <- eval.all > lambda.q[I] & eval.all <= lambda.q[I + 1]

##     plot(eta.all[in.lambda.bin], etahat.all[in.lambda.bin], cex=0.25)
## #         xlim=c(-1, 1)*max(abs(eta.all)),
## #         ylim=c(-1, 1)*max(abs(etahat.all)))
##     abline(0, 1)
    
##     print(cor(eta.all[in.lambda.bin], etahat.all[in.lambda.bin]))

##     readLines(n=1)
## }

## dhsEF <- 0

## y <- etahat.all**2
## x1 <- eval.all * dhs.all
## x2 <- eval.all

## mod <- lm(y ~ x1 + x2)
## summary(mod)

## dhsEF <- 
##     mod$coefficients["x1"] / mod$coefficients["x2"]
## dhsEF

######
# Partition by lambda 
# variance scale to sd scale

lambda.all <- eval.all
## lambda.all <- eval.all * (dhs.all * dhsEF + 1)

lambda.q <- w.quant.o2(sqrt(lambda.all), nLambdaPT)
lambda.q <- lambda.q ** 2
lambda.q[1] <- 0
lambda.q[nLambdaPT + 1] <- lambda.q[nLambdaPT + 1] * 1.1

print(lambda.q)

#pdf(file=plotfile)

for (I in nLambdaPT:1) {
    in.lambda.bin <- lambda.all > lambda.q[I] & lambda.all <= lambda.q[I + 1]

#    plot(eta.all[in.lambda.bin], etahat.all[in.lambda.bin], cex=0.25)
#    abline(0, 1)

    print(sum(in.lambda.bin))

# for debugging
#    print(sd(eta.all[in.lambda.bin]))
#    print(sd(etahat.all[in.lambda.bin]))
#    print(cor(eta.all[in.lambda.bin], etahat.all[in.lambda.bin]))

    ##    readLines(n=1) # for interactive mode
}

#dev.off()

######
# Partition by eta

betahatH.q <- matrix(NA, nrow=(nEtaPT + 1), ncol=nLambdaPT)

count <- 0
nBetahatH <- array(0, dim=c(nLambdaPT, nEtaPT, nDHSPT))
meanBetahatH <- array(0, dim=c(nLambdaPT, nEtaPT, nDHSPT))

for (I in 2:length(lambda.q)) {
    etahat.all.sub <- etahat.all[lambda.all > lambda.q[I - 1] &
                                 lambda.all <= lambda.q[I]]
    etahat.all.sub <- abs(etahat.all.sub)

#    dhs.all.sub <- dhs.all[lambda.all > lambda.q[I - 1] &
#                           lambda.all <= lambda.q[I]]

    betahatH.q[, I - 1] <- w.quant.o2(etahat.all.sub, nEtaPT)
    betahatH.q[1, I - 1] <- 0
    betahatH.q[nEtaPT + 1, I - 1] <- betahatH.q[nEtaPT + 1, I - 1] * 1.1 
    
    for (J in 1:nEtaPT) {

        betahatH.lo <- betahatH.q[J, I - 1]
        betahatH.hi <- betahatH.q[J+1, I - 1]
        print(sum(etahat.all.sub > betahatH.lo &
                  etahat.all.sub <= betahatH.hi))

        nBetahatH[I - 1, J, 1] <- nBetahatH[I - 1, J, 1] +
            sum(etahat.all.sub > betahatH.lo &
                etahat.all.sub <= betahatH.hi)
        meanBetahatH[I - 1, J, 1] <- meanBetahatH[I - 1, J, 1] +
            sum(etahat.all.sub[(etahat.all.sub > betahatH.lo &
                                etahat.all.sub <= betahatH.hi)])

        count <- count + sum(nBetahatH[I - 1, J, ])
    }
}

ASSERT(count == length(lambda.all))

print(betahatH.q)

meanBetahatH <- meanBetahatH / nBetahatH

meanBetahatH[is.nan(meanBetahatH)] <- 0

# print(meanBetahatH)


# Save partition boundaries 
partdata <- list()

partdata[["Nt"]] <- Nt
partdata[["nLambdaPT"]] <- nLambdaPT
partdata[["nEtaPT"]] <- nEtaPT
partdata[["nDHSPT"]] <- nDHSPT

partdata[["lambda.q"]] <- lambda.q
partdata[["betahatH.q"]] <- betahatH.q

if (WINSHIFT == 0) {
    saveRDS(partdata, paste(tempprefix, "part.RDS", sep=''))
} else {
    saveRDS(partdata,
            paste(tempprefix, "win_", WINSHIFT, ".part.RDS", sep=''))
}

save.image(file=paste(tempprefix, "run_nps5svd2.", "win_", WINSHIFT, ".RData",
               sep=''))

cat("Done\n")
