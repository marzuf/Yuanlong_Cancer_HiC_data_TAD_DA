# Rscript cmp_meanCorr_topdomTADs.R

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

script_name <- "cmp_meanCorr_topdomTADs.R"
cat("... start ", script_name, "\n")
startTime <- Sys.time()

require(foreach)
require(doMC)

registerDoMC(40)

outFolder <- "CMP_MEANCORR_TOPDOMTADS"
dir.create(outFolder)

script0_name <- "0_prepGeneData"
script3_name <- "3_runMeanTADLogFC"
script4_name <- "4_runMeanTADCorr"
script8_name <- "8c_runAllDown"
script11_name <- "11_runEmpPvalCombined"
script9_name <- "9_runEmpPvalMeanTADLogFC"
script10_name <- "10_runEmpPvalMeanTADCorr"

pipFolder <- file.path("PIPELINE","OUTPUT_FOLDER")
td_folder <- file.path("..", "Cancer_HiC_data_TAD_DA")
topdomFolder <- file.path(td_folder, pipFolder)

td_pipFolder <- file.path(td_folder, pipFolder)

plotType <- "png"
myHeight <- 400
myWidth <- 600

buildAllData <- FALSE

# PIPELINE/OUTPUT_FOLDER/ENCSR444WCZ_A549_40kb/TCGAluad_norm_luad/4_runMeanTADCorr/all_meanCorr_TAD.Rdata

all_files <- list.files(pipFolder, recursive = TRUE, pattern="all_meanCorr_TAD.Rdata", full.names = FALSE)

curr_file = "Panc1_rep12_40kb/TCGApaad_wt_mutKRAS/4_runMeanTADCorr/all_meanCorr_TAD.Rdata"

if(buildAllData){

  all_data <- foreach(curr_file = all_files) %dopar% {
    
    curr_ds <- dirname(dirname(curr_file))
    
    hicds <- dirname(curr_ds)
    exprds <- basename(curr_ds)
    
    
    
    if(hicds == "ENCSR489OCU_NCI-H460_40kb") {
      hicds_td <- "NCI-H460_40kb"
    } else {
      hicds_td <- hicds
    }
    
    curr_ds_td <- file.path(hicds_td, exprds)
    
    td_hicds <- file.path(td_folder, hicds_td)
    
    ### YUANLONG DATA
    # KEEP ONLY THE TADs USED IN THE PIPELINE
    stopifnot(dir.exists(file.path(pipFolder,curr_ds, script0_name)))
    yl_tadListFile <- file.path(pipFolder,curr_ds, script0_name, "pipeline_regionList.Rdata")
    stopifnot(file.exists(yl_tadListFile))
    yl_pipeline_tadList <- eval(parse(text = load(yl_tadListFile))) # not adjusted
    # RETRIEVE THE GENES USED IN THE PIPELINE - script0
    yl_geneListFile <- file.path(pipFolder,curr_ds, script0_name, "pipeline_geneList.Rdata")
    stopifnot(file.exists(yl_geneListFile))
    yl_pipeline_geneList <- eval(parse(text = load(yl_geneListFile))) # not adjusted
    
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
    # KEEP ONLY THE TADs USED IN THE PIPELINE
    stopifnot(dir.exists(file.path(td_pipFolder,curr_ds_td, script0_name)))
    td_tadListFile <- file.path(td_pipFolder,curr_ds_td, script0_name, "pipeline_regionList.Rdata")
    stopifnot(file.exists(td_tadListFile))
    td_pipeline_tadList <- eval(parse(text = load(td_tadListFile))) # not adjusted
    # RETRIEVE THE GENES USED IN THE PIPELINE - script0
    td_geneListFile <- file.path(td_pipFolder,curr_ds_td, script0_name, "pipeline_geneList.Rdata")
    stopifnot(file.exists(td_geneListFile))
    td_pipeline_geneList <- eval(parse(text = load(td_geneListFile))) # not adjusted
    
    ### GENES and TADs INFO
    td_tad_DT_file <- file.path(td_folder, hicds_td, "genes2tad", "all_assigned_regions.txt")
    stopifnot(file.exists(td_tad_DT_file))
    td_tad_DT <- read.delim(td_tad_DT_file, header=F, 
                            col.names=c("chromo", "region", "start", "end"), stringsAsFactors = FALSE)
    stopifnot(is.numeric(td_tad_DT$start))
    stopifnot(is.numeric(td_tad_DT$end))
    td_tad_DT <- td_tad_DT[grepl("_TAD", td_tad_DT$region),,drop=FALSE] 
    
    td_g2tFile <- file.path(td_folder, hicds_td, "genes2tad", "all_genes_positions.txt")
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
    
    # LOAD TAD MEAN CORRELATION DATA
    td_corrData <- eval(parse(text  = load(file.path(topdomFolder, curr_ds_td, script4_name, "all_meanCorr_TAD.Rdata"))))
    
    # LOAD TAD CONCORDANCE
    td_concordData <- eval(parse(text  = load(file.path(topdomFolder, curr_ds_td, script8_name, "all_obs_prodSignedRatio.Rdata"))))
    
    # LOAD TAD ratioFC
    td_ratioDownData <- eval(parse(text  = load(file.path(topdomFolder, curr_ds_td, script8_name, "all_obs_ratioDown.Rdata"))))
    
    
    # LOAD TAD meanFC
    td_meanFCdata <- eval(parse(text  = load(file.path(topdomFolder, curr_ds_td, script3_name, "all_meanLogFC_TAD.Rdata"))))
    
    
    # LOAD TAD emp pval combined
    td_pvalData <- eval(parse(text  = load(file.path(topdomFolder, curr_ds_td, script11_name, "emp_pval_combined.Rdata"))))
    td_pvalData <- p.adjust(td_pvalData, method="BH")
    
    # LOAD TAD corr pval
    td_corrPvalData <- eval(parse(text  = load(file.path(topdomFolder, curr_ds_td, script10_name, "emp_pval_meanCorr.Rdata"))))
    td_corrPvalData <- p.adjust(td_corrPvalData, method="BH")
    
    # LOAD TAD fc pval
    td_fcPvalData <- eval(parse(text  = load(file.path(topdomFolder, curr_ds_td, script9_name, "emp_pval_meanLogFC.Rdata"))))
    td_fcPvalData <- p.adjust(td_fcPvalData, method="BH")
    
  
    ### 1) YL VS. TOPDOM MEAN TAD CORRELATION  
    outFile <- file.path(outFolder, paste0(dirname(curr_ds), "_", basename(curr_ds), "_", "TD_vs_YL_meanIntraCorr_density.", plotType))
    do.call(plotType, list(outFile, height=myHeight, width=myWidth))
    plot_multiDens_argList(list(
      yl_meanCorr = yl_corrData,
      td_meanCorr = td_corrData
    ),
    my_xlab = "mean intra-TAD correlation")
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
    
    ### 2) YL VS. TOPDOM # genes in TADs
    outFile <- file.path(outFolder, paste0(dirname(curr_ds), "_", basename(curr_ds), "_", "TD_vs_YL_nbrGenesByTAD_density.", plotType))
    do.call(plotType, list(outFile, height=myHeight, width=myWidth))
    plot_multiDens_argList(list(
      yl_nGenesTADs = yl_tadNgenes,
      td_nGenesTADs = td_tadNgenes
    ), 
    my_xlab = "# genes/TAD")
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
    
    ### 3) YL VS. TOPDOM size of TADs
    outFile <- file.path(outFolder, paste0(dirname(curr_ds), "_", basename(curr_ds), "_", "TD_vs_YL_TADsize_density.", plotType))
    do.call(plotType, list(outFile, height=myHeight, width=myWidth))
    plot_multiDens_argList(list(
      yl_TADsize = yl_tadSize/1000,
      td_TADsize = td_tadSize/1000
    ), 
    my_xlab = "TAD size (kb)")
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
    
    ### 4) YL VS. TOPDOM intraTAD concordance
    outFile <- file.path(outFolder, paste0(dirname(curr_ds), "_", basename(curr_ds), "_", "TD_vs_YL_intraConcord_density.", plotType))
    do.call(plotType, list(outFile, height=myHeight, width=myWidth))
    plot_multiDens_argList(list(
      yl_FCC = yl_concordData,
      td_FCC = td_concordData
    ), 
    my_xlab = "FCC score")
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
    
    ### 5) YL VS. TOPDOM ratioDown
    outFile <- file.path(outFolder, paste0(dirname(curr_ds), "_", basename(curr_ds), "_", "TD_vs_YL_ratioDown_density.", plotType))
    do.call(plotType, list(outFile, height=myHeight, width=myWidth))
    plot_multiDens_argList(list(
      yl_ratioDown = yl_ratioDownData,
      td_ratioDown = td_ratioDownData
    ), 
    my_xlab = "TAD ratioDown")
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
    
    
    ### 6) YL VS. TOPDOM meanFC
    outFile <- file.path(outFolder, paste0(dirname(curr_ds), "_", basename(curr_ds), "_", "TD_vs_YL_meanFC_density.", plotType))
    do.call(plotType, list(outFile, height=myHeight, width=myWidth))
    plot_multiDens_argList(list(
      yl_meanFC = yl_meanFCdata,
      td_meanFC = td_meanFCdata
    ), 
    my_xlab = "mean TAD FC")
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
    
    
    ### 7) YL VS. TOPDOM emp pval
    outFile <- file.path(outFolder, paste0(dirname(curr_ds), "_", basename(curr_ds), "_", "TD_vs_YL_empPvalComb_density.", plotType))
    do.call(plotType, list(outFile, height=myHeight, width=myWidth))
    plot_multiDens_argList(list(
      yl_adjEmpPval = -log10(yl_pvalData),
      td_adjEmpPval = -log10(td_pvalData)
    ), 
    my_xlab = "TAD emp. pval (-log10)")
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
    
    # densplot empPvalCombined vs. meanFC - YL
    outFile <- file.path(outFolder, paste0(dirname(curr_ds), "_", basename(curr_ds), "_", "YL_empPvalComb_vs_meanFC.", plotType))
    do.call(plotType, list(outFile, height=myHeight, width=myHeight))
    densplot(x = -log10(yl_pvalData), y = yl_meanFCdata,
             xlab = "-log10 emp. pval combined",
             ylab = "mean FC")
    mtext(side=3, text = "YL data")
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
    
    # densplot empPvalCombined vs. empPval meanFC - YL
    outFile <- file.path(outFolder, paste0(dirname(curr_ds), "_", basename(curr_ds), "_", "YL_empPvalComb_vs_empPvalFC.", plotType))
    do.call(plotType, list(outFile, height=myHeight, width=myHeight))
    densplot(x = -log10(yl_pvalData), y = -log10(yl_fcPvalData),
             xlab = "-log10 emp. pval combined",
             ylab = "-log10 emp. pval meanFC")
    mtext(side=3, text = "YL data")
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
    
    # densplot empPvalCombined vs. intraCorr - YL
    outFile <- file.path(outFolder, paste0(dirname(curr_ds), "_", basename(curr_ds), "_", "YL_empPvalComb_vs_meanIntraCorr.", plotType))
    do.call(plotType, list(outFile, height=myHeight, width=myHeight))
    densplot(x = -log10(yl_pvalData), y = yl_corrData,
             xlab = "-log10 emp. pval combined",
             ylab = "intra-TAD corr.")
    mtext(side=3, text = "YL data")
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
    
    # densplot empPvalCombined vs. empPvalIntraCorr - YL
    outFile <- file.path(outFolder, paste0(dirname(curr_ds), "_", basename(curr_ds), "_", "YL_empPvalComb_vs_empPvalCorr.", plotType))
    do.call(plotType, list(outFile, height=myHeight, width=myHeight))
    densplot(x = -log10(yl_pvalData), y = -log10(yl_corrPvalData),
             xlab = "-log10 emp. pval combined",
             ylab = "-log10 emp. pval intraCorr")
    mtext(side=3, text = "YL data")
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
    
    # densplot empPvalCombined vs. meanFC - TD
    outFile <- file.path(outFolder, paste0(dirname(curr_ds), "_", basename(curr_ds), "_", "TD_empPvalComb_vs_meanFC.", plotType))
    do.call(plotType, list(outFile, height=myHeight, width=myHeight))
    densplot(x = -log10(td_pvalData), y = td_meanFCdata,
             xlab = "-log10 emp. pval combined",
             ylab = "mean FC")
    mtext(side=3, text = "TD data")
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
    
    # densplot empPvalCombined vs. empPval meanFC - TD
    outFile <- file.path(outFolder, paste0(dirname(curr_ds), "_", basename(curr_ds), "_", "TD_empPvalComb_vs_empPvalFC.", plotType))
    do.call(plotType, list(outFile, height=myHeight, width=myHeight))
    densplot(x = -log10(td_pvalData), y = -log10(td_fcPvalData),
             xlab = "-log10 emp. pval combined",
             ylab = "-log10 emp. pval meanFC")
    mtext(side=3, text = "TD data")
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
    
    # densplot empPvalCombined vs. intraCorr - TD
    outFile <- file.path(outFolder, paste0(dirname(curr_ds), "_", basename(curr_ds), "_", "TD_empPvalComb_vs_meanIntraCorr.", plotType))
    do.call(plotType, list(outFile, height=myHeight, width=myHeight))
    densplot(x = -log10(td_pvalData), y = td_corrData,
             xlab = "-log10 emp. pval combined",
             ylab = "intra-TAD corr.")
    mtext(side=3, text = "TD data")
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
    
    # densplot empPvalCombined vs. empPvalIntraCorr - TD
    outFile <- file.path(outFolder, paste0(dirname(curr_ds), "_", basename(curr_ds), "_", "TD_empPvalComb_vs_empPvalCorr.", plotType))
    do.call(plotType, list(outFile, height=myHeight, width=myHeight))
    densplot(x = -log10(td_pvalData), y = -log10(td_corrPvalData),
             xlab = "-log10 emp. pval combined",
             ylab = "-log10 emp. pval intraCorr")
    mtext(side=3, text = "TD data")
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
    
    
    list(
    yl_meanCorr = yl_corrData,
    td_meanCorr = td_corrData,
    yl_nGenesTADs = yl_tadNgenes,
    td_nGenesTADs = td_tadNgenes,
    yl_TADsizeKb = yl_tadSize/1000,
    td_TADsizeKb = td_tadSize/1000,
    yl_FCC = yl_concordData,
    td_FCC = td_concordData,
    yl_ratioDown = yl_ratioDownData,
    td_ratioDown = td_ratioDownData,
    yl_meanFC = yl_meanFCdata,
    td_meanFC = td_meanFCdata,
    yl_adjEmpPvalComb = yl_pvalData,
    td_adjEmpPvalComb = td_pvalData,
    yl_adjEmpPvalFC = yl_fcPvalData,
    td_adjEmpPvalFC = td_fcPvalData,
    yl_adjEmpPvalCorr = yl_corrPvalData,
    td_adjEmpPvalCorr = td_corrPvalData
    )
    
    
  } # end-foreach iterating over datasets
  names(all_data) <- dirname(dirname(all_files))
  outFile <- file.path(outFolder, "all_data.Rdata")
  save(all_data, file = outFile)
  cat(paste0("... written: ", outFile, "\n"))
} else { # end-if build table
  outFile <- file.path(outFolder, "all_data.Rdata")
  all_data <- eval(parse(text = load(outFile)))
}

nDS <- length(all_data)

all_vars <- c("meanCorr", "nGenesTADs", "TADsizeKb", "FCC", "ratioDown", "meanFC", "adjEmpPvalComb", "adjEmpPvalFC", "adjEmpPvalCorr")
curr_var = all_vars[1]

for(curr_var in all_vars) {
  
  td_var <- unlist(lapply(all_data, function(x) x[[paste0("td_", curr_var)]]))
  yl_var <- unlist(lapply(all_data, function(x) x[[paste0("yl_", curr_var)]]))
  
  outFile <- file.path(outFolder, paste0("all_datasets_TD_vs_YL_", curr_var, ".", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  plot_multiDens(
    list(YL_TADs = yl_var,
    TD_TADs = td_var),
    my_xlab = paste0(curr_var) ,
    plotTit = paste0("YL vs. TD data: ", curr_var)
  )
  mtext(side=3, text = paste0("(nDS = ", nDS, ")"))
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  if(grepl("adjEmpPval", curr_var)) {
    outFile <- file.path(outFolder, paste0("all_datasets_TD_vs_YL_", curr_var, "_log10.", plotType))
    do.call(plotType, list(outFile, height=myHeight, width=myWidth))
    plot_multiDens(
      list(YL_TADs = -log10(yl_var),
           TD_TADs = -log10(td_var)),
      my_xlab = paste0(curr_var, "[-log10]") ,
      plotTit = paste0("YL vs. TD data: ", curr_var)
    )
    mtext(side=3, text = paste0("(nDS = ", nDS, ")"))
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
  }
}


x_var <- "adjEmpPvalComb"
y_var <- "adjEmpPvalFC"

dataSource="td"
for(dataSource in c("td", "yl")) {
  myx <- unlist(lapply(all_data, function(x) x[[paste0(dataSource, "_", x_var)]]))
  myy <- unlist(lapply(all_data, function(x) x[[paste0(dataSource, "_", y_var)]]))
  
  outFile <- file.path(outFolder, paste0("all_datasets_", toupper(dataSource), "_", y_var, "_vs_", x_var, "_log10.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myHeight))
  densplot(
    x = -log10(myx),
    y = -log10(myy),
    xlab = paste0(x_var, " [-log10]") ,
    ylab = paste0(y_var, " [-log10]") ,
    main = paste0(y_var, " vs. ", x_var, " - ", toupper(dataSource))
  )
  mtext(side=3, text = paste0("(nDS = ", nDS, ")"))
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
}

y_var <- "adjEmpPvalCorr"

for(dataSource in c("td", "yl")) {
  myx <- unlist(lapply(all_data, function(x) x[[paste0(dataSource, "_", x_var)]]))
  myy <- unlist(lapply(all_data, function(x) x[[paste0(dataSource, "_", y_var)]]))
  
  outFile <- file.path(outFolder, paste0("all_datasets_", toupper(dataSource), "_", y_var, "_vs_", x_var, "_log10.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myHeight))
  densplot(
    x = -log10(myx),
    y = -log10(myy),
    xlab = paste0(x_var, " [-log10]") ,
    ylab = paste0(y_var, " [-log10]") ,
    main = paste0(y_var, " vs. ", x_var, " - ", toupper(dataSource))
  )
  mtext(side=3, text = paste0("(nDS = ", nDS, ")"))
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
}




td_var <- unlist(lapply(all_data, function(x) x[[paste0("td_", curr_var)]]))
yl_var <- unlist(lapply(all_data, function(x) x[[paste0("yl_", curr_var)]]))

outFile <- file.path(outFolder, paste0("all_datasets_TD_vs_YL_", curr_var, ".", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot_multiDens(
  list(YL_TADs = yl_var,
       TD_TADs = td_var),
  my_xlab = paste0(curr_var) ,
  plotTit = paste0("YL vs. TD data: ", curr_var)
)
mtext(side=3, text = paste0("(nDS = ", nDS, ")"))
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))




x_var <- "adjEmpPvalCorr"
y_var <- "adjEmpPvalFC"

dataSource="td"
for(dataSource in c("td", "yl")) {
  myx <- unlist(lapply(all_data, function(x) x[[paste0(dataSource, "_", x_var)]]))
  myy <- unlist(lapply(all_data, function(x) x[[paste0(dataSource, "_", y_var)]]))
  
  outFile <- file.path(outFolder, paste0("all_datasets_", toupper(dataSource), "_", y_var, "_vs_", x_var, "_log10.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myHeight))
  densplot(
    x = -log10(myx),
    y = -log10(myy),
    xlab = paste0(x_var, " [-log10]") ,
    ylab = paste0(y_var, " [-log10]") ,
    main = paste0(y_var, " vs. ", x_var, " - ", toupper(dataSource))
  )
  mtext(side=3, text = paste0("(nDS = ", nDS, ")"))
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
}


txt <- paste0(startTime, "\n", Sys.time(), "\n")
cat(paste0(txt))
cat(paste0("*** DONE: ", script_name, "\n"))




