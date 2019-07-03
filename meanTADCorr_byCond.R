
# Rscript meanTADCorr_byCond.R

withDiago <- FALSE
corrMethod <- "pearson"
  
script_name <- "meanTADCorr_byCond.R"

cat(paste0("> START ", script_name, "\n"))

startTime <- Sys.time()

setDir <- ""

withDiag <- FALSE

outFold <- "MEANTADCORR_BYCOND"
dir.create(outFold)

require(foreach)
require(doMC)
registerDoMC(40)

script0_name <- "0_prepGeneData"

pipOutFolder <- file.path("PIPELINE", "OUTPUT_FOLDER")

all_qqnorm_files <- list.files(pipOutFolder, recursive = TRUE, pattern="rna_qqnorm_rnaseqDT.Rdata", full.names = TRUE)
stopifnot(length(all_qqnorm_files) > 0)

qqnorm_file = all_qqnorm_files[1]

for(qqnorm_file in all_qqnorm_files) {
  
  hicds <- basename(dirname(dirname(dirname(qqnorm_file))))
  exprds <- basename(dirname(dirname(qqnorm_file)))
  
  settingFile <- file.path("PIPELINE", "INPUT_FILES", hicds, paste0("run_settings_", exprds, ".R"))
  stopifnot(file.exists(settingFile))
  source(settingFile)
  stopifnot(file.exists(sample1_file))
  stopifnot(file.exists(sample2_file))
  
  cond1_ID <- eval(parse(text = load(sample1_file)))
  cond2_ID <- eval(parse(text = load(sample2_file)))
  
  script0_name <- "0_prepGeneData"
  qqnormDTfile <- file.path(pipOutFolder, hicds, exprds,script0_name, "rna_qqnorm_rnaseqDT.Rdata")
  stopifnot(file.exists(qqnormDTfile))
  qqnormDT <- eval(parse(text = load(qqnormDTfile)))
  
  geneListFile <- file.path(pipOutFolder, hicds, exprds,script0_name, "pipeline_geneList.Rdata")
  stopifnot(file.exists(geneListFile))
  geneList <- eval(parse(text = load(geneListFile)))
  
  stopifnot(names(geneList) %in% rownames(qqnormDT))
  
  stopifnot(setequal(colnames(qqnormDT), c(cond1_ID, cond2_ID)))
  
  # Reorder the DT so that after I can subset the rows using geneList
  norm_rnaseqDT <- qqnormDT[names(geneList),]    
  
  stopifnot(all(rownames(norm_rnaseqDT) == names(geneList)))
  stopifnot(!any(duplicated(names(geneList))))
  
  g2tFile <- file.path(hicds, "genes2tad", "all_genes_positions.txt")
  stopifnot(file.exists(g2tFile))
  gene2tadDT <- read.delim(g2tFile, header=F, col.names = c("entrezID", "chromo", "start", "end", "region"), stringsAsFactors = F)
  gene2tadDT$entrezID <- as.character(gene2tadDT$entrezID)
  stopifnot(geneList %in% gene2tadDT$entrezID)
  gene2tadDT <- gene2tadDT[gene2tadDT$entrezID %in% geneList,]
  
  regionFile <- file.path(pipOutFolder, hicds, exprds,script0_name, "pipeline_regionList.Rdata")
  stopifnot(file.exists(regionFile))
  pipeline_regionList <- eval(parse(text = load(regionFile)))
  
  stopifnot(setequal(pipeline_regionList, gene2tadDT$region))
  
  all_regions <- unique(gene2tadDT$region)
  
  
  norm_rnaseqDT_cond1 <- norm_rnaseqDT[, cond1_ID]
  norm_rnaseqDT_cond2 <- norm_rnaseqDT[, cond2_ID]
  
  stopifnot(ncol(norm_rnaseqDT_cond1) == length(cond1_ID))
  stopifnot(ncol(norm_rnaseqDT_cond2) == length(cond2_ID))
  stopifnot(ncol(norm_rnaseqDT) == length(cond2_ID) + length(cond1_ID))
  
  rm(norm_rnaseqDT)
  
  
  
  ################################****************************************************************************************
  ################################********************************************* iterate over regions to compute intraTAD correlation - COND1
  ################################****************************************************************************************
  
  all_meanCorr_TAD_cond1 <- foreach(reg=all_regions, .combine='c') %dopar% {
    # cat(paste0("... doing region: ", reg, "/", length(all_regions), "\n"))
    reg_genes <- gene2tadDT$entrezID[gene2tadDT$region == reg]
    rowsToKeep <- which(geneList %in% reg_genes)
    subData <- as.data.frame(t(norm_rnaseqDT_cond1[rowsToKeep,,drop=F]))
    stopifnot(colnames(subData) == names(geneList[rowsToKeep]))
    # columns => the genes
    # rows => the samples
    ######################################################## BECAUSE I COULD HAVE DUPLICATED GENE ENTREZ ID FOR DIFFERENT NAMES !!!
    ### UPDATE: this should not happen in the latest version !!!
    # stopifnot(ncol(subData) == length(reg_genes))
    # stopifnot(ncol(subData) == length(geneList[which(geneList %in% reg_genes)]))
    stopifnot(ncol(subData) == length(reg_genes))
    #### CORRELATION
    corrMatrix_all <- cor(subData, method = corrMethod)
    # should be correlation of the genes
    ######################################################## BECAUSE I COULD HAVE DUPLICATED GENE ENTREZ ID FOR DIFFERENT NAMES !!!
    ### UPDATE: this should not happen in the latest version !!!
    # stopifnot(nrow(corrMatrix_all) == length(geneList[which(geneList %in% reg_genes)]))
    # stopifnot(ncol(corrMatrix_all) == length(geneList[which(geneList %in% reg_genes)]))
    stopifnot(nrow(corrMatrix_all) == length(reg_genes))
    stopifnot(ncol(corrMatrix_all) == length(reg_genes))
    mean(corrMatrix_all[lower.tri(corrMatrix_all, diag = withDiag)], na.rm=T)
  }
  
  cat(paste0("... end intra TAD correlation\n"))
  
  names(all_meanCorr_TAD_cond1) <- all_regions
  stopifnot(length(all_meanCorr_TAD_cond1) == length(all_regions))
  
  outFile <- file.path(outFold, hicds, exprds, "all_meanCorr_TAD_cond1.Rdata")
  dir.create(dirname(outFile), recursive = TRUE)
  save(all_meanCorr_TAD_cond1, file = outFile)
  cat(paste0("... written: ", outFile,"\n"))
  
  
  
  
  ################################****************************************************************************************
  ################################********************************************* iterate over regions to compute intraTAD correlation - COND2
  ################################****************************************************************************************
  
  all_meanCorr_TAD_cond2 <- foreach(reg=all_regions, .combine='c') %dopar% {
    # cat(paste0("... doing region: ", reg, "/", length(all_regions), "\n"))
    reg_genes <- gene2tadDT$entrezID[gene2tadDT$region == reg]
    rowsToKeep <- which(geneList %in% reg_genes)
    subData <- as.data.frame(t(norm_rnaseqDT_cond2[rowsToKeep,,drop=F]))
    stopifnot(colnames(subData) == names(geneList[rowsToKeep]))
    # columns => the genes
    # rows => the samples
    ######################################################## BECAUSE I COULD HAVE DUPLICATED GENE ENTREZ ID FOR DIFFERENT NAMES !!!
    ### UPDATE: this should not happen in the latest version !!!
    # stopifnot(ncol(subData) == length(reg_genes))
    # stopifnot(ncol(subData) == length(geneList[which(geneList %in% reg_genes)]))
    stopifnot(ncol(subData) == length(reg_genes))
    #### CORRELATION
    corrMatrix_all <- cor(subData, method = corrMethod)
    # should be correlation of the genes
    ######################################################## BECAUSE I COULD HAVE DUPLICATED GENE ENTREZ ID FOR DIFFERENT NAMES !!!
    ### UPDATE: this should not happen in the latest version !!!
    # stopifnot(nrow(corrMatrix_all) == length(geneList[which(geneList %in% reg_genes)]))
    # stopifnot(ncol(corrMatrix_all) == length(geneList[which(geneList %in% reg_genes)]))
    stopifnot(nrow(corrMatrix_all) == length(reg_genes))
    stopifnot(ncol(corrMatrix_all) == length(reg_genes))
    mean(corrMatrix_all[lower.tri(corrMatrix_all, diag = withDiag)], na.rm=T)
  }
  
  cat(paste0("... end intra TAD correlation\n"))
  
  names(all_meanCorr_TAD_cond2) <- all_regions
  stopifnot(length(all_meanCorr_TAD_cond2) == length(all_regions))
  
  outFile <- file.path(outFold, hicds, exprds, "all_meanCorr_TAD_cond2.Rdata")
  dir.create(dirname(outFile), recursive = TRUE)
  save(all_meanCorr_TAD_cond2, file = outFile)
  cat(paste0("... written: ", outFile,"\n"))
  
}

################################****************************************************************************************
################################****************************************************************************************

cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))




