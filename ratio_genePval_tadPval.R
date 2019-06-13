# Rscript ratio_genePval_tadPval.R

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

script_name <- "ratio_genePval_tadPval.R"
cat("... start ", script_name, "\n")
startTime <- Sys.time()

require(foreach)
require(doMC)

registerDoMC(40)

outFolder <- "RATIO_GENEPVAL_TADPVAL"
dir.create(outFolder)

script0_name <- "0_prepGeneData"
script1_name <- "1_runGeneDE"
# script3_name <- "3_runMeanTADLogFC"
# script4_name <- "4_runMeanTADCorr"
# script8_name <- "8c_runAllDown"
script11_name <- "11_runEmpPvalCombined"
# script9_name <- "9_runEmpPvalMeanTADLogFC"
# script10_name <- "10_runEmpPvalMeanTADCorr"

yl_pipFolder <- file.path("PIPELINE","OUTPUT_FOLDER")
td_folder <- file.path("..", "Cancer_HiC_data_TAD_DA")
td_pipFolder <- file.path(td_folder, yl_pipFolder)

plotType <- "png"
myHeight <- 400
myWidth <- 600

axisCex <- 1.2

buildTable <- TRUE

# PIPELINE/OUTPUT_FOLDER/ENCSR444WCZ_A549_40kb/TCGAluad_norm_luad/4_runMeanTADCorr/all_meanCorr_TAD.Rdata

all_files <- list.files(yl_pipFolder, recursive = TRUE, pattern="emp_pval_combined.Rdata", full.names = FALSE)

stopifnot(length(all_files) > 0)

curr_file = "Panc1_rep12_40kb/TCGApaad_wt_mutKRAS/4_runMeanTADCorr/all_meanCorr_TAD.Rdata"

curr_file <- all_files[1]

if(buildTable) {
  
  
  all_ranks_DT <- foreach(curr_file = all_files, .combine='rbind') %dopar% {
    
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
    
    ### YUANLONG DATA - list of genes
    yl_geneListFile <- file.path(yl_pipFolder,curr_ds, script0_name, "pipeline_geneList.Rdata")
    stopifnot(file.exists(yl_geneListFile))
    # stopifnot(yl_pipeline_geneList %in% yl_g2t_DT$entrezID) # => TRUE
    yl_pipeline_geneList <- eval(parse(text = load(yl_geneListFile))) # not adjusted
    
    yl_g2tFile <- file.path(hicds, "genes2tad", "all_genes_positions.txt")
    stopifnot(file.exists(yl_g2tFile))
    yl_g2t_DT <- read.delim(yl_g2tFile, header=F, 
                            col.names = c("entrezID",  "chromo", "start", "end", "region"), stringsAsFactors = FALSE)
    yl_g2t_DT$entrezID <- as.character(yl_g2t_DT$entrezID)
    stopifnot(yl_pipeline_geneList %in% yl_g2t_DT$entrezID)
    yl_g2t_DT <- yl_g2t_DT[yl_g2t_DT$entrezID %in% yl_pipeline_geneList,]
    yl_g2t <- setNames( as.character(yl_g2t_DT$region),  as.character(yl_g2t_DT$entrezID))
    
    
    ### YL DATA - COMBINED PVAL
    yl_pvalData <- eval(parse(text  = load(file.path(yl_pipFolder, curr_ds, script11_name, "emp_pval_combined.Rdata"))))
    yl_pvalData <- p.adjust(yl_pvalData, method="BH")
    yl_tadRank <- rank(yl_pvalData, ties="min")
    
    ### YL GENE DE DATA
    yl_geneDE <- eval(parse(text  = load(file.path(yl_pipFolder, curr_ds, script1_name, "DE_topTable.Rdata")))) 
    yl_geneDE$genes <- as.character(yl_geneDE$genes)
    stopifnot(names(yl_pipeline_geneList) %in% yl_geneDE$genes)
    yl_geneDE <- yl_geneDE[yl_geneDE$genes %in% names(yl_pipeline_geneList),,drop=FALSE]
    # stopifnot((yl_pipeline_geneList) %in% yl_geneDE$genes) # FALSE
    stopifnot(!duplicated(yl_geneDE$genes))
    yl_geneDE$geneEntrez <- sapply(yl_geneDE$genes, function(x) {
      tmp <- pipeline_geneList[names(pipeline_geneList) == x]
      stopifnot(length(tmp) == 1)
      tmp
      })
    stopifnot(!duplicated(yl_geneDE$geneEntrez))
    stopifnot((yl_pipeline_geneList) %in% yl_geneDE$geneEntrez) # 
    stopifnot(yl_geneDE$geneEntrez %in% yl_g2t_DT$entrezID) # 
    yl_genePval <- setNames(yl_geneDE$adj.P.Val,  as.character(yl_geneDE$geneEntrez))
    
    
    ### TD DATA - list of genes
    td_geneListFile <- file.path(td_pipFolder,curr_ds_td, script0_name, "pipeline_geneList.Rdata")
    stopifnot(file.exists(td_geneListFile))
    td_pipeline_geneList <- eval(parse(text = load(td_geneListFile))) # not adjusted
    
    td_g2tFile <- file.path(td_folder, hicds_td, "genes2tad", "all_genes_positions.txt")
    stopifnot(file.exists(td_g2tFile))
    td_g2t_DT <- read.delim(td_g2tFile, header=F, col.names = c("entrezID",  "chromo", "start", "end", "region"), stringsAsFactors = FALSE)
    td_g2t_DT$entrezID <- as.character(td_g2t_DT$entrezID)
    stopifnot(td_pipeline_geneList %in% td_g2t_DT$entrezID)
    td_g2t_DT <- td_g2t_DT[td_g2t_DT$entrezID %in% td_pipeline_geneList,]
    stopifnot(td_pipeline_geneList %in% td_g2t_DT$entrezID)
    td_g2t <- setNames( as.character(td_g2t_DT$region),  as.character(td_g2t_DT$entrezID))
    
    # td_LOAD TAD emp pval
    td_pvalData <- eval(parse(text  = load(file.path(td_pipFolder, curr_ds_td, script11_name, "emp_pval_combined.Rdata"))))
    td_pvalData <- p.adjust(td_pvalData, method="BH")
    td_tadRank <- rank(td_pvalData, ties="min")
    
    ### TD GENE DE DATA
    td_geneDE <- eval(parse(text  = load(file.path(td_pipFolder, curr_ds_td, script1_name, "DE_topTable.Rdata")))) 
    td_geneDE$genes <- as.character(td_geneDE$genes)
    stopifnot(names(td_pipeline_geneList) %in% td_geneDE$genes)
    td_geneDE <- td_geneDE[td_geneDE$genes %in% names(td_pipeline_geneList),,drop=FALSE]
    # stopifnot((td_pipeline_geneList) %in% td_geneDE$genes) # FALSE
    stopifnot(!duplicated(td_geneDE$genes))
    td_geneDE$geneEntrez <- sapply(td_geneDE$genes, function(x) {
      tmp <- pipeline_geneList[names(pipeline_geneList) == x]
      stopifnot(length(tmp) == 1)
      tmp
    })
    stopifnot(!duplicated(td_geneDE$geneEntrez))
    stopifnot((td_pipeline_geneList) %in% td_geneDE$geneEntrez) # 
    stopifnot(td_geneDE$geneEntrez %in% td_g2t_DT$entrezID) # 
    td_genePval <- setNames(td_geneDE$adj.P.Val,  as.character(td_geneDE$geneEntrez))
    

        
    ### COMMON GENES
    common_geneList <- intersect(yl_pipeline_geneList, td_pipeline_geneList)
    
    stopifnot(common_geneList %in% names(yl_g2t))
    stopifnot(yl_g2t %in% names(yl_tadRank))
    
    stopifnot(common_geneList %in% names(td_g2t))
    stopifnot(td_g2t %in% names(td_tadRank))
    
    stopifnot(common_geneList %in% names(td_genePval))
    stopifnot(common_geneList %in% names(yl_genePval))
    
    data.frame(
      hicds = hicds,
      exprds = exprds,
      geneID = common_geneList,
      yl_rank = as.numeric(yl_tadRank[yl_g2t[common_geneList]]),
      yl_adjPvalComb = as.numeric(yl_pvalData[yl_g2t[common_geneList]]),
      yl_geneAdjPval = as.numeric(yl_genePval[common_geneList]),
      td_rank = as.numeric(td_tadRank[td_g2t[common_geneList]]),
      td_adjPvalComb = as.numeric(td_pvalData[td_g2t[common_geneList]]),
      td_geneAdjPval = as.numeric(td_genePval[common_geneList]),
      stringsAsFactors = FALSE
    )
  } # end-foreach iterating over datasets
  outFile <- file.path(outFolder, "all_ranks_DT.Rdata")
  save(all_ranks_DT, file=outFile)
  cat(paste0("... written: ", outFile, "\n"))
} else {
  outFile <- file.path(outFolder, "all_ranks_DT.Rdata")
  all_ranks_DT <- eval(parse(text = load(outFile)))
}

nDS <- length(unique(paste0(all_ranks_DT$hicds, "_", all_ranks_DT$exprds)))


all_ranks_DT$yl_pvalRatio <- -log10(all_ranks_DT$yl_adjPvalComb)/-log10(all_ranks_DT$yl_geneAdjPval)
all_ranks_DT$td_pvalRatio <-  -log10(all_ranks_DT$td_adjPvalComb)/-log10(all_ranks_DT$td_geneAdjPval)

### YL VS TD - GENE PVAL 
# myx <- all_ranks_DT$td_pvalRatio
# myx <- myx[myx >=  quantile(myx, 0.05) & myx <  quantile(myx, 0.95)]
# myy <- all_ranks_DT$yl_pvalRatio
# myy <- myy[myy >=  quantile(myy, 0.05) & myy <  quantile(myy, 0.95)]
# outFile <- file.path(outFolder, paste0("yl_vs_td_pvalRatio.", plotType))
# do.call(plotType, list(outFile, height=myHeight, width=myWidth))
# plot_multiDens(
#   list(myx,myy)
# )
myx <- all_ranks_DT$td_pvalRatio
myy <- all_ranks_DT$yl_pvalRatio
toKeep <- which(myx >=  quantile(myx, 0.05) & myx <  quantile(myx, 0.95) & 
  myy >=  quantile(myy, 0.05) & myy <  quantile(myy, 0.95))
myx <- myx[toKeep]
myy <- myy[toKeep]
outFile <- file.path(outFolder, paste0("yl_vs_td_pvalRatio.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myHeight))
densplot(
  x = myx,
  xlab = paste0("TD - -log10(pval.) ratio"),
  y = myy,
  ylab = paste0("YL - -log10(pval.) ratio"),
  main = "YL vs. TD pval. ratio (no outlier)",
  cex = 0.7,
  cex.lab = axisCex,
  cex.axis = axisCex
)  
mtext(side=3, text = paste0("(nDS = ", nDS, ")"))
addCorr(x=myx, legPos="topleft",
        y=myy, bty='n')
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))



### YL VS TD - GENE PVAL 
myx <- all_ranks_DT$td_geneAdjPval
myy <- all_ranks_DT$yl_geneAdjPval
outFile <- file.path(outFolder, paste0("yl_vs_td_genePval.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myHeight))
densplot(
  x = myx,
  xlab = paste0("TD - gene adj. pval."),
  y = myy,
  ylab = paste0("YL - gene adj. pval."),
  main = "YL vs. TD gene adj. pval.",
  cex = 0.7,
  cex.lab = axisCex,
  cex.axis = axisCex
)  
mtext(side=3, text = paste0("(nDS = ", nDS, ")"))
addCorr(x=myx, legPos="topleft",
        y=myy, bty='n')
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))



### YL VS TD - LOG10 GENE PVAL 
myx <- log10(all_ranks_DT$td_geneAdjPval)
myy <- log10(all_ranks_DT$yl_geneAdjPval)
outFile <- file.path(outFolder, paste0("yl_vs_td_log10genePval.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myHeight))
densplot(
  x = -myx,
  xlab = paste0("TD - gene adj. pval. [-log10]"),
  y = -myy,
  ylab = paste0("YL - gene adj. pval. [-log10]"),
  main = "YL vs. TD gene adj. pval.",
  cex = 0.7,
  cex.lab = axisCex,
  cex.axis = axisCex
)  
mtext(side=3, text = paste0("(nDS = ", nDS, ")"))
addCorr(x=myx, legPos="topleft",
        y=myy, bty='n')
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


### YL - LOG10 GENE PVAL  vs TAD pval
myx <- log10(all_ranks_DT$yl_adjPvalComb)
myy <- log10(all_ranks_DT$yl_geneAdjPval)
outFile <- file.path(outFolder, paste0("yl_log10_tadPval_vs_genePval.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myHeight))
densplot(
  x = -myx,
  xlab = paste0("YL - TAD adj. pval. [-log10]"),
  y = -myy,
  ylab = paste0("YL - gene adj. pval. [-log10]"),
  main = "YL gene vs. TAD adj. pval.",
  cex = 0.7,
  cex.lab = axisCex,
  cex.axis = axisCex
)  
mtext(side=3, text = paste0("(nDS = ", nDS, ")"))
addCorr(x=myx, legPos="topleft",
        y=myy, bty='n')
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


### YL - LOG10 GENE PVAL  vs TAD pval
myx <- log10(all_ranks_DT$td_adjPvalComb)
myy <- log10(all_ranks_DT$td_geneAdjPval)
outFile <- file.path(outFolder, paste0("td_log10_tadPval_vs_genePval.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myHeight))
densplot(
  x = -myx,
  xlab = paste0("TD - TAD adj. pval. [-log10]"),
  y = -myy,
  ylab = paste0("TD - gene adj. pval. [-log10]"),
  main = "TD gene vs. TAD adj. pval.",
  cex = 0.7,
  cex.lab = axisCex,
  cex.axis = axisCex
)  
mtext(side=3, text = paste0("(nDS = ", nDS, ")"))
addCorr(x=myx, legPos="topleft",
        y=myy, bty='n')
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


# # densplot empPvalCombined vs. empPvalIntraCorr - TD
# outFile <- file.path(outFolder, paste0(dirname(curr_ds), "_", basename(curr_ds), "_", "TD_empPvalComb_vs_empPvalCorr.", plotType))
# do.call(plotType, list(outFile, height=myHeight, width=myWidth))
# densplot(x = -log10(td_pvalData), y = -log10(td_corrPvalData),
#          xlab = "-log10 emp. pval combined",
#          ylab = "-log10 emp. pval intraCorr")
# mtext(side=3, text = "TD data")
# foo <- dev.off()
# cat(paste0("... written: ", outFile, "\n"))

txt <- paste0(startTime, "\n", Sys.time(), "\n")
cat(paste0(txt))
cat(paste0("*** DONE: ", script_name, "\n"))




