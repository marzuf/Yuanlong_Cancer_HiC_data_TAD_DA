startTime <- Sys.time()
cat(paste0("> Rscript create_angDist_byCond_sortNoDup_otherTADfile.R\n"))

# 

suppressPackageStartupMessages(library(reshape2, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))

library(lsa)
get_ang_dist <- function(vect1, vect2) {
  acos(cosine(vect1,vect2))/pi
}

corMethod <- "pearson"


options(scipen=100)

require(foreach)
require(doMC)

registerDoMC(80)

# Rscript create_coexpr_byCond_sortNoDup_otherTADfile.R ENCSR079VIJ_G401_40kb TCGAkich_norm_kich

### UPDATE sortNoDup 30.06.2018
# -> sort rows of the expression table to ensure alphabetical order of the genes
#    so that after melt the 1st gene will always be the 1st in alphabetical order
# -> replace diag + upper.tri with NA, use melt with na.rm=TRUE
#    so that data saved are smaller !!!
# curr_TADlist = "ENCSR079VIJ_G401_40kb"
# curr_dataset = "TCGAkich_norm_kich"
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2)
curr_TADlist <- args[1]
curr_dataset <- args[2]


SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")
# registerDoMC(ifelse(SSHFS, 2, 40))

outFold <- file.path("CREATE_ANGDIST_BYCOND_SORTNODUP", curr_TADlist, paste0(curr_dataset), corMethod)
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

# rm(qqnormDT)
rm(qqnormDT_tmp)

############################################################################## 
############################################################################## angDIST FOR COND1
############################################################################## 

stopifnot(rownames(qqnormDT) == rownames(qqnormDT_cond1))
stopifnot(rownames(qqnormDT) == rownames(qqnormDT_cond2))


cat("... build corMat \n")
cor_qqnormMat <- cor(t(qqnormDT), method = corMethod)
stopifnot(nrow(cor_qqnormMat) == nrow(qqnormDT))
stopifnot(ncol(cor_qqnormMat) == nrow(qqnormDT))
stopifnot(rownames(cor_qqnormMat) == colnames(cor_qqnormMat))
stopifnot(rownames(cor_qqnormMat) == rownames(qqnormDT))
rm(qqnormDT)


cat("... build corMat cond1\n")
cor_qqnormMat_cond1 <- cor(t(qqnormDT_cond1), method = corMethod)
stopifnot(nrow(cor_qqnormMat_cond1) == nrow(qqnormDT_cond1))
stopifnot(ncol(cor_qqnormMat_cond1) == nrow(qqnormDT_cond1))
stopifnot(rownames(cor_qqnormMat_cond1) == colnames(cor_qqnormMat_cond1))
stopifnot(rownames(cor_qqnormMat_cond1) == rownames(qqnormDT_cond1))
rm(qqnormDT_cond1)


cat("... build corMat cond2\n")
cor_qqnormMat_cond2 <- cor(t(qqnormDT_cond2), method = corMethod)
stopifnot(nrow(cor_qqnormMat_cond2) == nrow(qqnormDT_cond2))
stopifnot(ncol(cor_qqnormMat_cond2) == nrow(qqnormDT_cond2))
stopifnot(rownames(cor_qqnormMat_cond2) == colnames(cor_qqnormMat_cond2))
stopifnot(rownames(cor_qqnormMat_cond2) == rownames(qqnormDT_cond2))
rm(qqnormDT_cond2)


stopifnot(rownames(cor_qqnormMat) == rownames(cor_qqnormMat_cond1))
stopifnot(rownames(cor_qqnormMat) == rownames(cor_qqnormMat_cond2))


cat(" create gene comb.\n")

all_genes_cmbs <- combn(rownames(cor_qqnormMat), 2)

stopifnot(all_genes_cmbs[1,] %in% rownames(cor_qqnormMat))
stopifnot(all_genes_cmbs[2,] %in% rownames(cor_qqnormMat))


# all_genes_cmbs = all_genes_cmbs[,1:10]

cat(" create angDist - all cond.\n")

all_angDist <- apply(all_genes_cmbs, 2, function(gene_pair) {
  # cat(".")
  gene1 <- gene_pair[1]
  gene2 <- gene_pair[2]
  as.numeric(get_ang_dist(cor_qqnormMat[paste0(gene1),],
                          cor_qqnormMat[paste0(gene2),]))
})
# cat("\n")

all_angDist_DT <- data.frame(
  gene1 = all_genes_cmbs[1,],
  gene2 = all_genes_cmbs[2,],
  angDist = all_angDist,
  stringsAsFactors = FALSE
)

# all_angDist_DT <- foreach(i = seq_len(ncombs), .combine='rbind') %dopar% {
# cat("... i", "/", ncombs, "\n")
#   gene1 <- all_genes_cmbs[1,i]
#   gene2 <- all_genes_cmbs[2,i]
#   # stopifnot(gene1 %in% rownames(cor_qqnormMat))
#   # stopifnot(gene2 %in% rownames(cor_qqnormMat))
#   data.frame(
#     gene1 = gene1,
#     gene2 = gene2,
#     angDist = as.numeric(get_ang_dist(cor_qqnormMat[paste0(gene1),],
#                                       cor_qqnormMat[paste0(gene2),])),
#     stringsAsFactors = FALSE
#   )
# }

stopifnot(nrow(all_angDist_DT) == ncol(all_genes_cmbs))

outFile <- file.path(outFold, "all_angDist_DT.Rdata" )
save(all_angDist_DT, file = outFile)
cat(paste0("... written: ", outFile, "\n"))


############################################################################## 
############################################################################## angDIST FOR COND1
############################################################################## 

# all_genes_cmbs <- combn(rownames(qqnormDT_cond1), 2)


cat(" create angDist - cond. 1\n")

stopifnot(all_genes_cmbs[1,] %in% rownames(cor_qqnormMat_cond1))
stopifnot(all_genes_cmbs[2,] %in% rownames(cor_qqnormMat_cond1))


all_angDist_cond1 <- apply(all_genes_cmbs, 2, function(gene_pair) {
  # cat(".")
  gene1 <- gene_pair[1]
  gene2 <- gene_pair[2]
  as.numeric(get_ang_dist(cor_qqnormMat_cond1[paste0(gene1),],
                          cor_qqnormMat_cond1[paste0(gene2),]))
})
# cat("\n")

all_angDist_DT_cond1 <- data.frame(
  gene1 = all_genes_cmbs[1,],
  gene2 = all_genes_cmbs[2,],
  angDist = all_angDist_cond1,
  stringsAsFactors = FALSE
)

# all_angDist_DT_cond1 <- foreach(i = seq_len(ncol(all_genes_cmbs)), .combine='rbind') %dopar% {
#   gene1 <- all_genes_cmbs[1,i]
#   gene2 <- all_genes_cmbs[2,i]
#   # stopifnot(gene1 %in% rownames(cor_qqnormMat_cond1))
#   # stopifnot(gene2 %in% rownames(cor_qqnormMat_cond1))
#   data.frame(
#     gene1 = gene1,
#     gene2 = gene2,
#     angDist = as.numeric(get_ang_dist(cor_qqnormMat_cond1[paste0(gene1),],
#                                       cor_qqnormMat_cond1[paste0(gene2),])),
#     stringsAsFactors = FALSE
#     )
# }


stopifnot(nrow(all_angDist_DT_cond1) == ncol(all_genes_cmbs))

outFile <- file.path(outFold, "all_angDist_DT_cond1.Rdata" )
save(all_angDist_DT_cond1, file = outFile)
cat(paste0("... written: ", outFile, "\n"))

############################################################################## 
############################################################################## angDIST FOR COND2
############################################################################## 

# all_genes_cmbs <- combn(rownames(qqnormDT_cond2), 2)


cat(" create angDist - cond. 2\n")

stopifnot(all_genes_cmbs[1,] %in% rownames(cor_qqnormMat_cond2))
stopifnot(all_genes_cmbs[2,] %in% rownames(cor_qqnormMat_cond2))


all_angDist_cond2 <- apply(all_genes_cmbs, 2, function(gene_pair) {
  gene1 <- gene_pair[1]
  gene2 <- gene_pair[2]
  as.numeric(get_ang_dist(cor_qqnormMat_cond2[paste0(gene1),],
                          cor_qqnormMat_cond2[paste0(gene2),]))
})
all_angDist_DT_cond2 <- data.frame(
  gene1 = all_genes_cmbs[1,],
  gene2 = all_genes_cmbs[2,],
  angDist = all_angDist_cond2,
  stringsAsFactors = FALSE
)
# cat("\n")
# all_angDist_DT_cond2 <- foreach(i = seq_len(ncol(all_genes_cmbs)), .combine='rbind') %dopar% {
#   gene1 <- all_genes_cmbs[1,i]
#   gene2 <- all_genes_cmbs[2,i]
#   # stopifnot(gene1 %in% rownames(cor_qqnormMat_cond2))
#   # stopifnot(gene2 %in% rownames(cor_qqnormMat_cond2))
#   data.frame(
#     gene1 = gene1,
#     gene2 = gene2,
#     angDist = as.numeric(get_ang_dist(cor_qqnormMat_cond2[paste0(gene1),],
#                                       cor_qqnormMat_cond2[paste0(gene2),])),
#     stringsAsFactors = FALSE
#   )
# }



stopifnot(nrow(all_angDist_DT_cond2) == ncol(all_genes_cmbs))



outFile <- file.path(outFold, "all_angDist_DT_cond2.Rdata" )
save(all_angDist_DT_cond2, file = outFile)
cat(paste0("... written: ", outFile, "\n"))




######################################################################################
######################################################################################
######################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))
