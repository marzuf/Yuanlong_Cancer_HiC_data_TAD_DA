source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

script_name <- "cmp_meanCorr_topdomTADs.R"
cat("... start ", script_name, "\n")
startTime <- Sys.time()

script0_name <- "0_prepGeneData"
script3_name <- "3_runMeanTADLogFC"
script4_name <- "4_runMeanTADCorr"
script8_name <- "8c_runAllDown"
script11_name <- "11_runEmpPvalCombined"
script9_name <- "9_runEmpPvalMeanTADLogFC"
script10_name <- "10_runEmpPvalMeanTADCorr"

pipFolder <- file.path("PIPELINE","OUTPUT_FOLDER")
td_folder <- file.path("..", "Cancer_HiC_data_TAD_DA")

td_pipFolder <- file.path(td_folder, pipFolder)

# PIPELINE/OUTPUT_FOLDER/ENCSR444WCZ_A549_40kb/TCGAluad_norm_luad/4_runMeanTADCorr/all_meanCorr_TAD.Rdata

all_files <- list.files(pipFolder, recursive = TRUE, pattern="all_meanCorr_TAD.Rdata", full.names = FALSE)

curr_file = "Panc1_rep12_40kb/TCGApaad_wt_mutKRAS/4_runMeanTADCorr/all_meanCorr_TAD.Rdata"
for(curr_file in all_files) {
  
  curr_ds <- dirname(dirname(curr_file))
  
  hicds <- dirname(curr_ds)
  td_hicds <- file.path(td_folder, hicds)
  
  ### YUANLONG DATA
  ### GENES and TADs INFO
  yl_tad_DT_file <- file.path(hicds, "genes2tad", "all_assigned_regions.txt")
  stopifnot(file.exists(yl_tad_DT_file))
  yl_tad_DT <- read.delim(yl_tad_DT_file, header=F, 
                          col.names=c("chromo", "region", "start", "end"), stringsAsFactors = FALSE)
  stopifnot(is.numeric(yl_tad_DT$start))
  stopifnot(is.numeric(yl_tad_DT$end))
  yl_tad_DT <- yl_tad_DT[grepl("_TAD", yl_tad_DT$region),,drop=FALSE] 
  
  yl_g2tFile <- file.path(hicds, "genes2tad", "all_genes_positions.txt")
  stopifnot(file.exists(yl_g2tFile))
  yl_g2t_DT <- read.delim(yl_g2tFile, header=F, 
                          col.names = c("entrezID",  "chromo", "start", "end", "region"), stringsAsFactors = FALSE)
  yl_g2t_DT$entrezID <- as.character(yl_g2t_DT$entrezID)
  stopifnot(yl_pipeline_geneList %in% yl_g2t_DT$entrezID)
  
  yl_g2t_DT <- yl_g2t_DT[yl_g2t_DT$entrezID %in% yl_pipeline_geneList,]
  stopifnot(yl_g2t_DT$region %in% yl_pipeline_tadList)
  stopifnot(yl_pipeline_tadList %in% yl_g2t_DT$region)
  # retrieve # genes per TAD
  yl_tadNgenes <- setNames(as.numeric(table(yl_g2t_DT$region)),names(table(yl_g2t_DT$region)))
  
  # retrieve size of the TADs
  yl_tad_DT <- yl_tad_DT[yl_tad_DT$region %in% yl_pipeline_tadList,]
  stopifnot(!duplicated(yl_tad_DT$region))
  yl_tad_DT$tad_size <- yl_tad_DT$end-yl_tad_DT$start+1
  yl_tadSize <- setNames(yl_tad_DT$tad_size, yl_tad_DT$region)
  
  # KEEP ONLY THE TADs USED IN THE PIPELINE
  stopifnot(dir.exists(file.path(pipFolder,curr_ds, script0_name)))
  yl_tadListFile <- file.path(pipFolder,curr_ds, script0_name, "pipeline_regionList.Rdata")
  stopifnot(file.exists(yl_tadListFile))
  yl_pipeline_tadList <- eval(parse(text = load(yl_tadListFile))) # not adjusted
  # RETRIEVE THE GENES USED IN THE PIPELINE - script0
  yl_geneListFile <- file.path(pipFolder,curr_ds, script0_name, "pipeline_geneList.Rdata")
  stopifnot(file.exists(yl_geneListFile))
  yl_pipeline_geneList <- eval(parse(text = load(yl_geneListFile))) # not adjusted
  # LOAD TAD MEAN CORRELATION DATA
  yl_corrData <- eval(parse(text = load(file.path(pipFolder, curr_file))))
  
  # LOAD TAD CONCORDANCE
  yl_concordData <- eval(parse(text  = load(file.path(pipFolder, curr_ds, script8_name, "all_obs_prodSignedRatio.Rdata"))))
  
  # LOAD TAD ratioFC
  yl_ratioDownData <- eval(parse(text  = load(file.path(pipFolder, curr_ds, script8_name, "all_obs_ratioDown.Rdata"))))
  
  # LOAD TAD meanFC
  yl_meanFCdata <- eval(parse(text  = load(file.path(pipFolder, curr_ds, script3_name, "all_meanLogFC_TAD.Rdata"))))
  
  # LOAD TAD emp pval
  yl_pvalData <- eval(parse(text  = load(file.path(pipFolder, curr_ds, script11_name, "emp_pval_combined.Rdata"))))
  yl_pvalData <- p.adjust(yl_pvalData, method="BH")
  
  # LOAD TAD corr pval
  yl_corrPvalData <- eval(parse(text  = load(file.path(pipFolder, curr_ds, script10_name, "emp_pval_meanCorr.Rdata"))))
  yl_corrPvalData <- p.adjust(yl_corrPvalData, method="BH")
  
  # LOAD TAD fc pval
  yl_fcPvalData <- eval(parse(text  = load(file.path(pipFolder, curr_ds, script9_name, "emp_pval_meanLogFC.Rdata"))))
  yl_fcPvalData <- p.adjust(yl_fcPvalData, method="BH")
  
  
  ### TOPDOM DATA
  ### GENES and TADs INFO
  td_tad_DT_file <- file.path(td_folder, hicds, "genes2tad", "all_assigned_regions.txt")
  stopifnot(file.exists(td_tad_DT_file))
  td_tad_DT <- read.delim(td_tad_DT_file, header=F, 
                          col.names=c("chromo", "region", "start", "end"), stringsAsFactors = FALSE)
  stopifnot(is.numeric(td_tad_DT$start))
  stopifnot(is.numeric(td_tad_DT$end))
  td_tad_DT <- td_tad_DT[grepl("_TAD", td_tad_DT$region),,drop=FALSE] 
  
  td_g2tFile <- file.path(td_folder, hicds, "genes2tad", "all_genes_positions.txt")
  stopifnot(file.exists(td_g2tFile))
  td_g2t_DT <- read.delim(td_g2tFile, header=F, col.names = c("entrezID",  "chromo", "start", "end", "region"), stringsAsFactors = FALSE)
  td_g2t_DT$entrezID <- as.character(td_g2t_DT$entrezID)
  stopifnot(td_pipeline_geneList %in% td_g2t_DT$entrezID)
  
  td_g2t_DT <- td_g2t_DT[td_g2t_DT$entrezID %in% td_pipeline_geneList,]
  stopifnot(td_g2t_DT$region %in% td_pipeline_tadList)
  stopifnot(td_pipeline_tadList %in% td_g2t_DT$region)
  # retrieve # genes per TAD
  td_tadNgenes <- setNames(as.numeric(table(td_g2t_DT$region)),names(table(td_g2t_DT$region)))
  
  # retrieve size of the TADs
  td_tad_DT <- td_tad_DT[td_tad_DT$region %in% td_pipeline_tadList,]
  stopifnot(!duplicated(td_tad_DT$region))
  td_tad_DT$tad_size <- td_tad_DT$end-td_tad_DT$start+1
  td_tadSize <- setNames(td_tad_DT$tad_size, td_tad_DT$region)
  
  # KEEP ONLY THE TADs USED IN THE PIPELINE
  stopifnot(dir.exists(file.path(td_pipFolder,curr_ds, script0_name)))
  td_tadListFile <- file.path(td_pipFolder,curr_ds, script0_name, "pipeline_regionList.Rdata")
  stopifnot(file.exists(td_tadListFile))
  td_pipeline_tadList <- eval(parse(text = load(td_tadListFile))) # not adjusted
  # RETRIEVE THE GENES USED IN THE PIPELINE - script0
  td_geneListFile <- file.path(td_pipFolder,curr_ds, script0_name, "pipeline_geneList.Rdata")
  stopifnot(file.exists(td_geneListFile))
  td_pipeline_geneList <- eval(parse(text = load(td_geneListFile))) # not adjusted
  # LOAD TAD MEAN CORRELATION DATA
  td_corrData <- eval(parse(text  = load(file.path(topdomFolder, curr_file))))

  # LOAD TAD CONCORDANCE
  td_concordData <- eval(parse(text  = load(file.path(topdomFolder, curr_ds, script8_name, "all_obs_prodSignedRatio.Rdata"))))
  
  # LOAD TAD ratioFC
  td_ratioDownData <- eval(parse(text  = load(file.path(topdomFolder, curr_ds, script8_name, "all_obs_ratioDown.Rdata"))))
  
  
  # LOAD TAD meanFC
  td_meanFCdata <- eval(parse(text  = load(file.path(topdomFolder, curr_ds, script3_name, "all_meanLogFC_TAD.Rdata"))))
  
  
  # LOAD TAD emp pval combined
  td_pvalData <- eval(parse(text  = load(file.path(topdomFolder, curr_ds, script11_name, "emp_pval_combined.Rdata"))))
  td_pvalData <- p.adjust(td_pvalData, method="BH")
  
  # LOAD TAD corr pval
  td_corrPvalData <- eval(parse(text  = load(file.path(topdomFolder, curr_ds, script10_name, "emp_pval_meanCorr.Rdata"))))
  td_corrPvalData <- p.adjust(td_corrPvalData, method="BH")
  
  # LOAD TAD fc pval
  td_fcPvalData <- eval(parse(text  = load(file.path(topdomFolder, curr_ds, script9_name, "emp_pval_meanLogFC.Rdata"))))
  td_fcPvalData <- p.adjust(td_fcPvalData, method="BH")
  

  
  
  
  ### 1) YL VS. TOPDOM MEAN TAD CORRELATION  
  plot_multiDens_argList(list(
    yl_meanCorr = yl_corrData,
    td_meanCorr = td_corrData
  ),
  my_xlab = "mean intra-TAD correlation")
  
  ### 2) YL VS. TOPDOM # genes in TADs
  plot_multiDens_argList(list(
    yl_nGenesTADs = yl_tadNgenes,
    td_nGenesTADs = td_tadNgenes
  ), 
  my_xlab = "# genes/TAD")
  
  ### 3) YL VS. TOPDOM size of TADs
  plot_multiDens_argList(list(
    yl_nGenesTADs = yl_tadSize/1000,
    td_nGenesTADs = td_tadSize/1000
  ), 
  my_xlab = "TAD size (kb)")
  
  ### 4) YL VS. TOPDOM intraTAD concordance
  plot_multiDens_argList(list(
    yl_FCC = yl_concordData,
    td_FCC = td_concordData
  ), 
  my_xlab = "FCC score")
  
  ### 5) YL VS. TOPDOM ratioDown
  
  plot_multiDens_argList(list(
    yl_ratioDown = yl_ratioDownData,
    td_ratioDown = td_ratioDownData
  ), 
  my_xlab = "TAD ratioDown")
  
  
  ### 6) YL VS. TOPDOM meanFC
  
  plot_multiDens_argList(list(
    yl_meanFC = yl_meanFCdata,
    td_meanFC = td_meanFCdata
  ), 
  my_xlab = "mean TAD FC")
  
  
  ### 7) YL VS. TOPDOM emp pval
  plot_multiDens_argList(list(
    yl_adjEmpPval = -log10(yl_pvalData),
    td_adjEmpPval = -log10(td_pvalData)
  ), 
  my_xlab = "TAD emp. pval (-log10)")
  
  
  
  densplot(x = -log10(yl_pvalData), y = yl_meanFCdata,
           xlab = "-log10 emp. pval combined",
           ylab = "mean FC")
  mtext(side=3, text = "YL data")
  
  
  densplot(x = -log10(yl_pvalData), y = -log10(yl_fcPvalData),
           xlab = "-log10 emp. pval combined",
           ylab = "-log10 emp. pval meanFC")
  mtext(side=3, text = "YL data")
  
  
  densplot(x = -log10(yl_pvalData), y = yl_corrData,
           xlab = "-log10 emp. pval combined",
           ylab = "intra-TAD corr.")
  mtext(side=3, text = "YL data")
  
  densplot(x = -log10(yl_pvalData), y = -log10(yl_corrPvalData),
           xlab = "-log10 emp. pval combined",
           ylab = "-log10 emp. pval intraCorr")
  mtext(side=3, text = "YL data")
  
  
  densplot(x = -log10(td_pvalData), y = td_meanFCdata,
           xlab = "-log10 emp. pval combined",
           ylab = "mean FC")
  mtext(side=3, text = "TD data")
  
  densplot(x = -log10(td_pvalData), y = -log10(td_fcPvalData),
           xlab = "-log10 emp. pval combined",
           ylab = "-log10 emp. pval meanFC")
  mtext(side=3, text = "TD data")
  
  
  densplot(x = -log10(td_pvalData), y = td_corrData,
           xlab = "-log10 emp. pval combined",
           ylab = "intra-TAD corr.")
  mtext(side=3, text = "TD data")
  

  densplot(x = -log10(td_pvalData), y = -log10(td_corrPvalData),
           xlab = "-log10 emp. pval combined",
           ylab = "-log10 emp. pval intraCorr")
  mtext(side=3, text = "TD data")
  
  
}


txt <- paste0(startTime, "\n", Sys.time(), "\n")
cat(paste0(txt))
cat(paste0("*** DONE: ", script_name, "\n"))




