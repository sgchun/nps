args <- commandArgs(trailingOnly=TRUE)

if (length(args) != 12) {
    stop("Rscript run_nps5svd2.params.R <summstatfile> <betaprefix> <Nt> <traindir> <trainfamfile> <trainfreqfile> <trainphenofile> <traintag> <valdir> <valfamfile> <win size> <outdir>")
}

summstatfile <- args[1]
betaprefix <- args[2]
Nt <- as.numeric(args[3])

if (is.nan(Nt) || Nt <= 1) {
    stop(paste("Invalid training cohort size:", args[3]))
}

traindir <- args[4]
trainfamfile <- args[5]
trainfreqfile <- args[6]
trainphenofile <- args[7]
traintag <- args[8]
valdir <- args[9]
valfamfile <- args[10]

WINSZ <- as.numeric(args[11])

if (is.nan(WINSZ) || WINSZ <= 100) {
    stop(paste("Invalid window size:", args[11]))
}

tempprefix <- paste(args[12], "/", sep='')


#################################################################
# SANITY CHECKS

if (!file.exists(summstatfile)) {
   stop("File does not exists:", summstatfile)
}

if (!dir.exists(traindir)) {
   stop("Dir does not exists:", traindir)
}

if (!file.exists(trainfamfile)) {
   stop("File does not exists:", trainfamfile)
}

if (!file.exists(trainfreqfile)) {
   stop("File does not exists:", trainfreqfile)
}

if (!file.exists(trainphenofile)) {
   stop("File does not exists:", trainphenofile)
}

if (!dir.exists(valdir)) {
   stop("Dir does not exists:", valdir)
}

if (!file.exists(valfamfile)) {
   stop("File does not exists:", valfamfile)
}

summstat <- read.delim(summstatfile, header=TRUE, stringsAsFactors=FALSE,
                       sep="\t")
dim(summstat)

if (length(intersect(colnames(summstat), 
                     c("chr", "pos", "ref", "alt", "reffreq", "effalt"))) 
    != 6) {
    
    stop("Missing essential columns in", summstatfile, ":", 
         paste(c("chr", "pos", "ref", "alt", "reffreq", "effalt"), 
               collapse=", "))
}

if (nchar(summstat$chr[1]) < 4 || substr(summstat$chr[1], 1, 3) != "chr") {
    stop("chr column is expected to be chr1, chr2, ...:", summstatfile)
}

summstat.chr <- substr(summstat$chr, 4, nchar(summstat$chr))

if (length(setdiff(summstat.chr, as.character(1:22))) != 0) {
   stop("expect only chr1 through chr22 in ", summstatfile, ":", 
        paste(setdiff(summstat.chr, as.character(1:22)), collapse=", "))
}

summstat.chr <- as.numeric(summstat.chr)

if (is.unsorted(summstat.chr)) {
    stop("unsorted chr:", summstatfile)
}

for (CHR in 1:22) {
    if (is.unsorted(summstat$pos[summstat.chr == CHR])) {
        stop("unsorted pos:", summstatfile)
    }

    if (any(duplicated(summstat$pos[summstat.chr == CHR]))) {
        stop("duplicated pos in chr", CHR, " :", summstatfile)
    }
}

if (any(is.na(summstat$reffreq)) || any(is.na(summstat$effalt))) {
    stop("NA is not allowed in reffreq or effalt:", summstatfile)
}

summstat.SNPID <- paste(summstat.chr, ":", summstat$pos, "_", summstat$ref, 
    "_", summstat$alt, sep='')

# training freq
trfrq <- read.table(trainfreqfile, header=TRUE, stringsAsFactors=FALSE)

if (length(intersect(colnames(trfrq), 
                     c("SNPID", "MAF"))) != 2) {

    if (length(intersect(colnames(trfrq), 
                         c("SNP", "MAF"))) == 2) {

        trfrq.SNPID <- trfrq$SNP

    } else {
        stop("Missing essential columns in ", trainfreqfile, ":", 
                   paste(c("SNPID", "MAF"), 
                         collapse=", "))
    }

} else {

    trfrq.SNPID <- trfrq$SNPID

}

if (nrow(trfrq) != nrow(summstat)) {
    stop("The number of markers does not match:", summstatfile, ", ", 
         trainfreqfile)
}

# Temporarily disabled 
#
#if (any(trfrq.SNPID != summstat.SNPID)) {
##   cat(summstat.SNPID[which(trfrq.SNPID != summstat.SNPID)[1]], "...\n")
#
#   stop("The order of SNPIDs does not match:", summstatfile, ",", 
#         trainfreqfile)
#}

trfam <- read.delim(trainfamfile, sep=" ", header=FALSE,
                    stringsAsFactors=FALSE)

trphen <- read.delim(trainphenofile, sep="\t", header=TRUE,
                     stringsAsFactors=FALSE)

if (ncol(trfam) != 6) {    
    stop(trainfamfile, " does not have standard 6 columns")
}

if (length(intersect(colnames(trphen), 
                     c("FID", "IID", "Outcome"))) 
    != 3) {
    
    stop("Missing essential columns in ", trainphenfile, ":", 
         paste(c("FID", "IID", "Outcome"), 
               collapse=", "))
}

if (Nt != nrow(trfam)) {
    stop("The number of samples in ", trainfamfile, " does not match Nt =", Nt)
}

rownames(trphen) <- paste(trphen$FID, trphen$IID, sep=":")

if (length(intersect(paste(trfam[, 1], trfam[, 2], sep=":"),
                     paste(trphen$FID, trphen$IID, sep=":")
                     )) != Nt) {
    stop("Missing phenotypes for FID/IID pairs in ", trainfamfile, ":", 
         trainphenofile)
}

if (length(setdiff(trphen$Outcome, c(0, 1))) != 0) {
    stop("Expect only 0 or 1 in Outcome:", trainphenofile, ":", 
         paste(setdiff(trphen$Outcome, c("0", "1")), collapse=", "))
}

if (sum(trphen$Outcome == 0) == 0) {
    stop("Must have controls (Outcome = 0):", trainphenofile)
}

if (sum(trphen$Outcome == 1) == 0) {
    stop("Must have controls (Outcome = 1):", trainphenofile)
}

if (valfamfile != ".") {
    vlfam <- read.delim(valfamfile, sep=" ", header=FALSE,
                        stringsAsFactors=FALSE)

    if (ncol(vlfam) != 6) {    
        stop(valfamfile, " does not have standard 6 columns")
    }
}

#summstatfile <-
#    "/net/home/schun/gwassim2/mvn2/discovery/ldpred.23301_P0.01_hsq0.5_N_DHS5xPi.std_summstats"

#betaprefix <- "/net/home/schun/gwassim2/mvn2/discovery/23301_P0.01_hsq0.5_N_DHS5xPi."

#Nt <- 5000

#traindir <-
#    "/net/home/schun/gwassim2/mvn2/train/"

#trainfamfile <-
#    "/net/home/schun/gwassim2/mvn2/train/all3.maf_0.05.2.5K_2.5K.23301_P0.01_hsq0.5_N_DHS5xPi.fam"

#trainfreqfile <-
#    "/net/home/schun/gwassim2/mvn2/train/all3.maf_0.05.2.5K_2.5K.23301_P0.01_hsq0.5_N_DHS5xPi.frq"

#trainphenofile <-
#    "/net/home/schun/gwassim2/mvn2/train/total_liability.2.5K_2.5K.23301_P0.01_hsq0.5_N_DHS5xPi.phen"

#traintag <- "2.5K_2.5K.23301_P0.01_hsq0.5_N_DHS5xPi"

#valdir <-
#    "/net/home/schun/gwassim2/mvn2/val/"

#valfamfile <-
#    "/net/home/schun/gwassim2/mvn2/val/faminfo.maf_0.05.txt"

#tempprefix <- 
#    "/net/home/schun/gwassim2/mvn2/train/nps5.23301.4000/"


##################################################################
# SAVE
args <- list()

args[["summstatfile"]] <- summstatfile
args[["betaprefix"]] <- betaprefix

args[["Nt"]] <- Nt
args[["traindir"]] <- traindir
args[["trainfamfile"]] <- trainfamfile
args[["trainfreqfile"]] <- trainfreqfile
args[["trainphenofile"]] <- trainphenofile
args[["traintag"]] <- traintag

args[["valdir"]] <- valdir
args[["valfamfile"]] <- valfamfile

args[["WINSZ"]] <- WINSZ

print(args)

dir.create(tempprefix)

saveRDS(args, file=paste(tempprefix, "args.RDS", sep=''))

cat("Done\n")
