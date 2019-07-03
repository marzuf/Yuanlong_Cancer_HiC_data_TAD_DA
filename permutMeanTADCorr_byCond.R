
# Rscript permutMeanTADCorr_byCond.R


# Rscript permutMeanTADCorr_byCond.R

withDiago <- FALSE

script_name <- "permutMeanTADCorr_byCond.R"

cat(paste0("> START ", script_name, "\n"))

startTime <- Sys.time()

setDir <- ""

outFold <- "PERMUTMEANTADCORR_BYCOND"
dir.create(outFold)

require(foreach)
require(doMC)
registerDoMC(40)

###################


pipOutFolder <- file.path("PIPELINE", "OUTPUT_FOLDER")

all_permut_files <- list.files(pipOutFolder, recursive = TRUE, pattern="permutationsDT.Rdata", full.names = TRUE)
stopifnot(length(all_permut_files) > 0)


permut_file = all_permut_files[1]

args <- commandArgs(trailingOnly = TRUE)

if(length(args) == 2) {
  all_permut_files <- all_permut_files[grepl(args[1], all_permut_files) & grepl(args[2], all_permut_files)]
}

cat(all_permut_files)

stopifnot(length(all_permut_files) == 1)

for(permut_file in all_permut_files) {
  
  
  
  hicds <- basename(dirname(dirname(dirname(permut_file))))
  exprds <- basename(dirname(dirname(permut_file)))
  
  
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
  
  
  cat("... ", hicds, " - ", exprds, " - load permutation data\n")
  permutationsDT <- eval(parse(text = load(permut_file)))
  
  all_regions <- sort(unique(as.character(permutationsDT[,2])))
  
  stopifnot(setequal(all_regions, pipeline_regionList))
  
  
  norm_rnaseqDT_cond1 <- norm_rnaseqDT[, cond1_ID]
  norm_rnaseqDT_cond2 <- norm_rnaseqDT[, cond2_ID]
  
  stopifnot(ncol(norm_rnaseqDT_cond1) == length(cond1_ID))
  stopifnot(ncol(norm_rnaseqDT_cond2) == length(cond2_ID))
  stopifnot(ncol(norm_rnaseqDT) == length(cond2_ID) + length(cond1_ID))
  
  rm(norm_rnaseqDT)
  
  ################################****************************************************************************************
  ####################################################### COMPUTE MEAN INTRA CORR BY TAD FOR THE PERMUTATIONS  ----> COND1
  ################################****************************************************************************************
  
  cat("... start intraCorr permutDT for cond1 \n")
  
  intraTADcorr_permDT_allReg_cond1 <- foreach(i_col = 1:ncol(permutationsDT), .combine='cbind') %dopar% {
    #intraTADcorr_permDT_allReg <- foreach(i_col = 1:2, .combine='cbind') %dopar% {
    
    cat(paste0("... ", exprds, " - ", hicds, ", COND1 - intraTAD correlation for permutation: ", i_col, "/", ncol(permutationsDT), "\n"))
    
    g2t_permDT <- data.frame(entrezID = rownames(permutationsDT), 
                             region = permutationsDT[,i_col], stringsAsFactors = F)
    g2t_permDT$entrezID <- as.character(g2t_permDT$entrezID)
    g2t_permDT$region <-  as.character(g2t_permDT$region)
    
    permutCorr <- sapply(unique(all_regions), function(reg) {
      reg_genes <- g2t_permDT$entrezID[g2t_permDT$region == reg]
      stopifnot( reg_genes %in% geneList)
      subData <- as.data.frame(t(norm_rnaseqDT_cond1[which(geneList %in% reg_genes),,drop=F]))
      
      # THE REGIONS HERE ARE UNFILTERED, SO IT IS POSSIBLE THAT THERE ARE SOME REGIONS WITH ONLY 1 GENE
      # (these regions have not been used for the logReg_meanTAD etc.)
      if(ncol(subData) < 2)
        return(NA)
      # columns => the genes
      # rows => the samples
      stopifnot(nrow(subData) == ncol(norm_rnaseqDT_cond1))
      ################################################# BECAUSE OF POSSIBLE DUPLICATED GENE ENTREZ ID
      # UPDATE:  now should not have duplicates !!!
      stopifnot(ncol(subData) == length(reg_genes))
      # stopifnot(ncol(subData) == length(geneList[which(geneList %in% reg_genes)]))
      
      #### ALL CORRELATION
      corrMatrix_all <- cor(subData)
      # should be correlation of the genes
      ################################################# BECAUSE OF POSSIBLE DUPLICATED GENE ENTREZ ID
      # UPDATE:  now should not have duplicates !!!
      stopifnot(nrow(corrMatrix_all) == length(reg_genes))
      stopifnot(ncol(corrMatrix_all) == length(reg_genes))
      # stopifnot(ncol(corrMatrix_all) == length(geneList[which(geneList %in% reg_genes)]))
      # stopifnot(nrow(corrMatrix_all) == length(geneList[which(geneList %in% reg_genes)]))
      
      meanCorr_all <- mean(corrMatrix_all[lower.tri(corrMatrix_all, diag = withDiago)], na.rm=T)
      
    })
    # permutCorr <- rbindlist(permutCorr)
    curr_permutDT <- data.frame(permutCorr)
    stopifnot(ncol(curr_permutDT) == 1)
    colnames(curr_permutDT) <- paste0("result", i_col-1)
    stopifnot(all(rownames(curr_permutDT) == all_regions))
    curr_permutDT
  }
  cat("... end intraCorr permutDT \n")
  
  cat(paste0("*** DONE: ", script_name, "\n"))
  #stop("-- ok\n")
  
  meanCorr_permDT_cond1 <- as.data.frame(intraTADcorr_permDT_allReg_cond1)
  stopifnot(ncol(meanCorr_permDT_cond1) == ncol(permutationsDT))  
  colnames(meanCorr_permDT_cond1) <- paste0("permutation",  c(1:ncol(permutationsDT)))
  stopifnot(nrow(meanCorr_permDT_cond1) == length(all_regions))
  rownames(meanCorr_permDT_cond1) <- all_regions
  
  outFile <- file.path(outFold, hicds, exprds, "meanCorr_permDT_cond1.Rdata")
  dir.create(dirname(outFile), recursive = TRUE)
  save(meanCorr_permDT_cond1, file = outFile)
  cat(paste0("... written: ", outFile, "\n"))
  
  
  ################################****************************************************************************************
  ####################################################### COMPUTE MEAN INTRA CORR BY TAD FOR THE PERMUTATIONS  ----> COND2
  ################################****************************************************************************************
  
  cat("... start intraCorr permutDT for cond2  \n")
  
  intraTADcorr_permDT_allReg_cond2 <- foreach(i_col = 1:ncol(permutationsDT), .combine='cbind') %dopar% {
    #intraTADcorr_permDT_allReg <- foreach(i_col = 1:2, .combine='cbind') %dopar% {
    
    cat(paste0("... ", exprds, " - ", hicds, ", COND2 - intraTAD correlation for permutation: ", i_col, "/", ncol(permutationsDT), "\n"))
    
    g2t_permDT <- data.frame(entrezID = rownames(permutationsDT), 
                             region = permutationsDT[,i_col], stringsAsFactors = F)
    g2t_permDT$entrezID <- as.character(g2t_permDT$entrezID)
    g2t_permDT$region <-  as.character(g2t_permDT$region)
    
    permutCorr <- sapply(unique(all_regions), function(reg) {
      reg_genes <- g2t_permDT$entrezID[g2t_permDT$region == reg]
      stopifnot( reg_genes %in% geneList)
      subData <- as.data.frame(t(norm_rnaseqDT_cond2[which(geneList %in% reg_genes),,drop=F]))
      
      # THE REGIONS HERE ARE UNFILTERED, SO IT IS POSSIBLE THAT THERE ARE SOME REGIONS WITH ONLY 1 GENE
      # (these regions have not been used for the logReg_meanTAD etc.)
      if(ncol(subData) < 2)
        return(NA)
      # columns => the genes
      # rows => the samples
      stopifnot(nrow(subData) == ncol(norm_rnaseqDT_cond2))
      ################################################# BECAUSE OF POSSIBLE DUPLICATED GENE ENTREZ ID
      # UPDATE:  now should not have duplicates !!!
      stopifnot(ncol(subData) == length(reg_genes))
      # stopifnot(ncol(subData) == length(geneList[which(geneList %in% reg_genes)]))
      
      #### ALL CORRELATION
      corrMatrix_all <- cor(subData)
      # should be correlation of the genes
      ################################################# BECAUSE OF POSSIBLE DUPLICATED GENE ENTREZ ID
      # UPDATE:  now should not have duplicates !!!
      stopifnot(nrow(corrMatrix_all) == length(reg_genes))
      stopifnot(ncol(corrMatrix_all) == length(reg_genes))
      # stopifnot(ncol(corrMatrix_all) == length(geneList[which(geneList %in% reg_genes)]))
      # stopifnot(nrow(corrMatrix_all) == length(geneList[which(geneList %in% reg_genes)]))
      
      meanCorr_all <- mean(corrMatrix_all[lower.tri(corrMatrix_all, diag = withDiago)], na.rm=T)
      
    })
    # permutCorr <- rbindlist(permutCorr)
    curr_permutDT <- data.frame(permutCorr)
    stopifnot(ncol(curr_permutDT) == 1)
    colnames(curr_permutDT) <- paste0("result", i_col-1)
    stopifnot(all(rownames(curr_permutDT) == all_regions))
    curr_permutDT
  }
  cat("... end intraCorr permutDT \n")
  
  cat(paste0("*** DONE: ", script_name, "\n"))
  #stop("-- ok\n")
  
  meanCorr_permDT_cond2 <- as.data.frame(intraTADcorr_permDT_allReg_cond2)
  stopifnot(ncol(meanCorr_permDT_cond2) == ncol(permutationsDT))  
  colnames(meanCorr_permDT_cond2) <- paste0("permutation",  c(1:ncol(permutationsDT)))
  stopifnot(nrow(meanCorr_permDT_cond2) == length(all_regions))
  rownames(meanCorr_permDT_cond2) <- all_regions
  
  outFile <- file.path(outFold, hicds, exprds, "meanCorr_permDT_cond2.Rdata")
  dir.create(dirname(outFile), recursive = TRUE)
  save(meanCorr_permDT_cond2, file = outFile)
  cat(paste0("... written: ", outFile, "\n"))
  
}

################################****************************************************************************************
################################****************************************************************************************

cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))








