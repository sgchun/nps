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

if (length(cargs) < 5) {
    stop("Usage: Rscript nps_val.R <work dir> <val dir> <val fam file> <val pheno file> [<WINSHIFT>]+...")
}

tempprefix <- paste(cargs[1], "/", sep='')

args <- readRDS(paste(tempprefix, "args.RDS", sep=''))

traintag <- args[["traintag"]]
WINSZ <- args[["WINSZ"]]

valdir <- cargs[2]
valfamfile <- cargs[3]
valphenofile <- cargs[4]

list.WINSHIFT <- c() 

for (carg.WSHIFT in cargs[5:length(cargs)]) {

    WSHIFT <- as.numeric(carg.WSHIFT)

    if (is.nan(WSHIFT) || WSHIFT < 0 || WSHIFT > WINSZ) {
        stop(paste("Invalid shift:", carg.WSHIFT)) 
    }

    list.WINSHIFT <- c(list.WINSHIFT, WSHIFT)
}

if (!dir.exists(valdir)) {
   stop("Dir does not exists:", valdir)
}

if (!file.exists(valfamfile)) {
   stop("File does not exists:", valfamfile)
}

if (!file.exists(valphenofile)) {
   stop("File does not exists:", valphenofile)
}

#########################################################################
### validation

# phenotypes
vlfam <- read.delim(valfamfile, sep=" ", header=FALSE,
                    stringsAsFactors=FALSE)

if (ncol(vlfam) != 6) {
    # re-try with tab delimination

    vlfam <- read.delim(valfamfile, sep="\t", header=FALSE,
                        stringsAsFactors=FALSE)
}

if (ncol(vlfam) != 6) {    
    stop(valfamfile, " does not have standard 6 columns (space or tab-delimited)")
}

if (any(duplicated(paste(vlfam[, 1], vlfam[, 2], sep=":")))) {
    stop("Duplicated FID IID combinations:", valfamfile)
}

vlphen <- read.delim(valphenofile, sep="\t", header=TRUE,
                     stringsAsFactors=FALSE)

if (length(intersect(colnames(vlphen), c("FID", "IID", "Outcome"))) != 3) {
    stop(valphenofile, " does not include standard columns: FID IID Outcome (tab-delimited")
}

if (any(duplicated(paste(vlphen$FID, vlphen$IID, sep=":")))) {
    stop("Duplicated FID IID combinations:", valphenofile)
}

rownames(vlphen) <- paste(vlphen$FID, vlphen$IID, sep=":")

# No sample IDs match 
if (length(intersect(rownames(vlphen), paste(vlfam[, 1], vlfam[, 2], sep=":")))
    == 0) {
    stop("IID/FID does not match between .fam and phenotype files")
}

# samples in .fam but not phenotype file
missing.entry <- setdiff(paste(vlfam[, 1], vlfam[, 2], sep=" "),
                           paste(vlphen$FID, vlphen$IID, sep=" "))

if (length(missing.entry) > 0) {
    cat("FID IID\n")
    cat(paste(missing.entry, collapse="\n"))
    cat("\n")
    stop("The above samples declared in ", valfamfile,
         " are missing in the phenotype file: ",
         valphenofile)
}

vlphen <- vlphen[paste(vlfam[, 1], vlfam[, 2], sep=":"), ]

vlphen$Outcome[is.na(vlphen$Outcome)] <- -9

ASSERT(all(vlphen$FID == vlfam[, 1]))
ASSERT(all(vlphen$IID == vlfam[, 2]))

cat("Validation cohort:\n")
cat("Total ", nrow(vlfam), "samples\n")
cat(sum(vlphen$Outcome == -9), " samples with missing phenotype (-9)\n")

binary.phen <- TRUE

if (length(setdiff(vlphen$Outcome, c(0, 1, -9))) != 0) {
    cat(sum(vlphen$Outcome != -9),
        " samples with non-missing QUANTITATIVE phenotype values\n")

    if (!is.numeric(vlphen$Outcome)) {
        stop("phenotype values are not numeric: Outcome :", valphenofile)
    }

    binary.phen <- FALSE

} else {
    cat(sum(vlphen$Outcome != -9),
        " samples with non-missing BINARY phenotype values\n")
    
    cat("    ", sum(vlphen$Outcome == 1), " case samples\n")
    cat("    ", sum(vlphen$Outcome == 0), " control samples\n")

    binary.phen <- TRUE
}

if ("TotalLiability" %in% colnames(vlphen)) {
    cat("Includes TotalLiability\n")
}

if ((nrow(vlphen) <= 1) || (sum(vlphen$Outcome != -9) <= 1)) {
    stop("Invalid validation cohort size: N=", sum(vlphen$Outcome != -9))
}

if (binary.phen) {
    if (sum(vlphen$Outcome == 1) <= 1) {
        stop("Too few cases: N_case=", sum(vlphen$Outcome == 1))
    }

    if (sum(vlphen$Outcome == 0) <= 1) {
        stop("Too few controls: N_control=", sum(vlphen$Outcome == 0))
    }
}

# genetic risks
for (WINSHIFT in list.WINSHIFT) {

    cat("Checking a prediction model (winshift =", WINSHIFT, ")...\n")

    vlY <- vlphen$Outcome

    prisk <- rep(0, length(vlY))

    for (chr in 1:22) {

        # cat("chr", chr, "\n")

        # read per-chrom genetic risk file
        prisk.file <-
            paste(valdir, "/", traintag, ".win_", WINSHIFT,
                  ".predY.chrom", chr, ".txt", sep='')
    
        prisk.chr <- read.delim(prisk.file, header=FALSE, sep="\t")[, 1]
    
        ASSERT(length(prisk.chr) == length(vlY))
    
        prisk <- prisk + prisk.chr
    
    }

    prisk <- prisk[vlY != -9] 
    vlY <- vlY[vlY != -9]
    
    # R2 observed scale
    cat("Observed-scale R2 =", cor(vlY, prisk)**2, "\n")

    # R2 liability scale
    if ("TotalLiability" %in% colnames(vlphen)) {
        vlL <- vlphen$TotalLiability[vlY != -9]	
        cat("Liability-scale R2 =", cor(vlL, prisk)**2, "\n")
    }
}


cat("\n\n\n")
cat("Producing a combined prediction model...")

# Combined average 
vlY <- vlphen$Outcome

prisk <- rep(0, length(vlY))    

for (WINSHIFTx in list.WINSHIFT) {

    prisk0 <- rep(0, length(vlY))

    for (chr in 1:22) {

        # read per-chrom genetic risk file
        prisk.file <-
            paste(valdir, "/", traintag, ".win_", WINSHIFTx,
                  ".predY.chrom", chr, ".txt", sep='')

        prisk.chr <- read.delim(prisk.file, header=FALSE, sep="\t")[, 1]

        ASSERT(length(prisk.chr) == length(vlY))
        
        prisk0 <- prisk0 + prisk.chr

    }

    prisk <- prisk + prisk0
}

prisk <- prisk[vlY != -9] 
vlY <- vlY[vlY != -9]

# DEPRECATED 
# if ((cor(vlY, prisk) < 0) {
#     cat("WARNING: auto-correct sign flip: ", cor(vlY, prisk), "\n")
# 
#     prisk <- - prisk
# }

filename <- paste(valphenofile, ".nps_score", sep='')
df.out <- cbind(vlphen[vlphen$Outcome != -9, ], Score=prisk)
write.table(df.out, file=filename,
            row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)
cat("OK (saved in", filename, ")\n")

cat("Observed-scale R2 =", cor(vlY, prisk)**2, "\n")

if ("TotalLiability" %in% colnames(vlphen)) {
    vlL <- vlphen$TotalLiability[vlY != -9]	
    cat("Liability-scale R2 =", cor(vlL, prisk)**2, "\n")
}

if (binary.phen) {
    if (require(pROC)) {
    
        library(pROC)

        cat("AUC:\n")
        print(roc(cases=prisk[vlY == 1], controls=prisk[vlY == 0], ci=TRUE))

    } else {
        cat("Skip AUC calculation\n")
        cat("Please install pROC package to enable this\n")
    }

    if (require(DescTools)) {
    
        library(DescTools)
    
        mod <- glm(vlY ~ prisk, family=binomial(link="logit"))
    
        cat("Nagelkerke's R2 =", PseudoR2(mod, "Nagelkerke"), "\n")
        
        print(mod)

    } else {
        cat("Skip Nagelkerke's R2 calculation\n")
        cat("Please install DescTools package to enable this\n")
    }
}

cat("Done\n")
