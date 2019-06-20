startTime <- Sys.time()
cat(paste0("> Rscript create_coexpr_byCond_sortNoDup_otherTADfile.R\n"))

suppressPackageStartupMessages(library(reshape2, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))

options(scipen=100)

# Rscript create_coexpr_byCond_sortNoDup_otherTADfile.R ENCSR079VIJ_G401_40kb TCGAkich_norm_kich

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

outFold <- file.path("CREATE_COEXPR_BYCOND_SORTNODUP", curr_TADlist, paste0(curr_dataset), paste0(corMethod))
dir.create(outFold, recursive=TRUE)

# PIPELINE/OUTPUT_FOLDER/GSE105566_ENCFF358MNA_Panc1/TCGApaad_wt_mutKRAS/0_prepGeneData/pipeline_geneList.Rdata
dataset_pipDir <- file.path("PIPELINE", "OUTPUT_FOLDER", curr_TADlist, curr_dataset)

### ADDED 20.06.2019
# to retrieve the conditions
# load variable "sample1_file" and "sample2_file"
settingFile <- file.path("PIPELINE", "INPUT_FILES", curr_TADlist, paste0("run_settings_", curr_dataset, ".R"))
stopifnot(file.exists(settingFile))
source(settingFile)
stopifnot(file.exists(sample1_file))
stopifnot(file.exists(sample2_file))

cond1_ID <- eval(parse(text = load(sample1_file)))
cond2_ID <- eval(parse(text = load(sample2_file)))

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

### ADDED 20.06.2019 to split coexpression by condition
stopifnot(ncol(qqnormDT) == length(cond1_ID) + length(cond2_ID))

stopifnot( cond1_ID %in% colnames(qqnormDT) )
stopifnot( cond2_ID %in% colnames(qqnormDT) )

qqnormDT_cond1 <- qqnormDT[,cond1_ID]
qqnormDT_cond2 <- qqnormDT[,cond2_ID]
qqnormDT_tmp_cond1 <- qqnormDT_cond1
qqnormDT_tmp_cond2 <- qqnormDT_cond2

stopifnot( ncol(qqnormDT_tmp_cond2) + ncol(qqnormDT_tmp_cond1) == ncol(qqnormDT_tmp))
stopifnot( ncol(qqnormDT_cond2) + ncol(qqnormDT_cond1) == ncol(qqnormDT))

rm(qqnormDT)
rm(qqnormDT_tmp)

############################################################################## 
############################################################################## COEXPR FOR COND1
############################################################################## 

cor_qqnormMat_cond1 <- cor(t(qqnormDT_cond1), method = corMethod)
stopifnot(nrow(cor_qqnormMat_cond1) == nrow(qqnormDT_cond1))
stopifnot(ncol(cor_qqnormMat_cond1) == nrow(qqnormDT_cond1))

###### UPDATE 30.06.2018
# coexprDT_cond1 <- melt(cor_qqnormMat_cond1)
# colnames(coexprDT_cond1) <- c("gene1", "gene2", "coexpr")
# coexprDT_cond1 <- coexprDT_cond1[coexprDT_cond1$gene1 != coexprDT_cond1$gene2,]

cor_qqnormMat_cond1_NA <- cor_qqnormMat_cond1
cor_qqnormMat_cond1_NA[lower.tri(cor_qqnormMat_cond1_NA, diag=T)] <- NA
coexprDT_cond1 <- melt(cor_qqnormMat_cond1_NA, na.rm = T)
colnames(coexprDT_cond1) <- c("gene1", "gene2", "coexpr")
coexprDT_cond1$gene1 <- as.character(coexprDT_cond1$gene1)
coexprDT_cond1$gene2 <- as.character(coexprDT_cond1$gene2)
stopifnot(coexprDT_cond1$gene1 < coexprDT_cond1$gene2)

outFile <- file.path(outFold, "coexprDT_cond1.Rdata")
save(coexprDT_cond1, file = outFile)
cat(paste0("... written: ", outFile, "\n"))


############################################################################## 
############################################################################## COEXPR FOR COND2
############################################################################## 

cor_qqnormMat_cond2 <- cor(t(qqnormDT_cond2), method = corMethod)
stopifnot(nrow(cor_qqnormMat_cond2) == nrow(qqnormDT_cond2))
stopifnot(ncol(cor_qqnormMat_cond2) == nrow(qqnormDT_cond2))

###### UPDATE 30.06.2018
# coexprDT_cond2 <- melt(cor_qqnormMat_cond2)
# colnames(coexprDT_cond2) <- c("gene1", "gene2", "coexpr")
# coexprDT_cond2 <- coexprDT_cond2[coexprDT_cond2$gene1 != coexprDT_cond2$gene2,]

cor_qqnormMat_cond2_NA <- cor_qqnormMat_cond2
cor_qqnormMat_cond2_NA[lower.tri(cor_qqnormMat_cond2_NA, diag=T)] <- NA
coexprDT_cond2 <- melt(cor_qqnormMat_cond2_NA, na.rm = T)
colnames(coexprDT_cond2) <- c("gene1", "gene2", "coexpr")
coexprDT_cond2$gene1 <- as.character(coexprDT_cond2$gene1)
coexprDT_cond2$gene2 <- as.character(coexprDT_cond2$gene2)
stopifnot(coexprDT_cond2$gene1 < coexprDT_cond2$gene2)

outFile <- file.path(outFold, "coexprDT_cond2.Rdata")
save(coexprDT_cond2, file = outFile)
cat(paste0("... written: ", outFile, "\n"))


######################################################################################
######################################################################################
######################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))
