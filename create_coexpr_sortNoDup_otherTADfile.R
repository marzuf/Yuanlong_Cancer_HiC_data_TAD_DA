startTime <- Sys.time()
cat(paste0("> Rscript create_coexpr_sortNoDup_otherTADfile.R\n"))

suppressPackageStartupMessages(library(reshape2, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))

options(scipen=100)

# Rscript create_coexpr_sortNoDup_otherTADfile.R ENCSR079VIJ_G401_40kb TCGAkich_norm_kich

### UPDATE sortNoDup 30.06.2018
# -> sort rows of the expression table to ensure alphabetical order of the genes
#    so that after melt the 1st gene will always be the 1st in alphabetical order
# -> replace diag + upper.tri with NA, use melt with na.rm=TRUE
#    so that data saved are smaller !!!

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2)
curr_TADlist <- args[1]
curr_dataset <- args[2]
corMethod <- "pearson"


SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")
# registerDoMC(ifelse(SSHFS, 2, 40))

outFold <- file.path("CREATE_COEXPR_SORTNODUP", curr_TADlist, paste0(curr_dataset), paste0(corMethod))
dir.create(outFold, recursive=TRUE)

# PIPELINE/OUTPUT_FOLDER/GSE105566_ENCFF358MNA_Panc1/TCGApaad_wt_mutKRAS/0_prepGeneData/pipeline_geneList.Rdata
dataset_pipDir <- file.path("PIPELINE", "OUTPUT_FOLDER", curr_TADlist, curr_dataset)

script0_name <- "0_prepGeneData"

pipeline_geneList <- eval(parse(text = load(file.path(dataset_pipDir, script0_name, "pipeline_geneList.Rdata"))))

qqnormDT <- eval(parse(text = load(file.path(dataset_pipDir, script0_name, "rna_qqnorm_rnaseqDT.Rdata"))))
stopifnot(names(pipeline_geneList) %in% rownames(qqnormDT))
qqnormDT <- qqnormDT[rownames(qqnormDT) %in% names(pipeline_geneList),]
stopifnot(nrow(qqnormDT) == length(pipeline_geneList))
rownames(qqnormDT) <- pipeline_geneList[rownames(qqnormDT)]
stopifnot(setequal(pipeline_geneList, rownames(qqnormDT)))

################################################################ UPDATE 30.06.2018 -> sort the rownames to have gene1 < gene2
qqnormDT_sorted <- qqnormDT[sort(rownames(qqnormDT)),]
stopifnot(dim(qqnormDT_sorted) == dim(qqnormDT))
stopifnot(setequal(rownames(qqnormDT), rownames(qqnormDT_sorted)))
qqnormDT <- qqnormDT_sorted
qqnormDT_tmp <- qqnormDT

cor_qqnormMat <- cor(t(qqnormDT), method = corMethod)
stopifnot(nrow(cor_qqnormMat) == nrow(qqnormDT))
stopifnot(ncol(cor_qqnormMat) == nrow(qqnormDT))

###### UPDATE 30.06.2018
# coexprDT <- melt(cor_qqnormMat)
# colnames(coexprDT) <- c("gene1", "gene2", "coexpr")
# coexprDT <- coexprDT[coexprDT$gene1 != coexprDT$gene2,]

cor_qqnormMat_NA <- cor_qqnormMat
cor_qqnormMat_NA[lower.tri(cor_qqnormMat_NA, diag=T)] <- NA
coexprDT <- melt(cor_qqnormMat_NA, na.rm = T)
colnames(coexprDT) <- c("gene1", "gene2", "coexpr")
coexprDT$gene1 <- as.character(coexprDT$gene1)
coexprDT$gene2 <- as.character(coexprDT$gene2)
stopifnot(coexprDT$gene1 < coexprDT$gene2)

outFile <- file.path(outFold, "coexprDT.Rdata")
save(coexprDT, file = outFile)
cat(paste0("... written: ", outFile, "\n"))


######################################################################################
######################################################################################
######################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))
