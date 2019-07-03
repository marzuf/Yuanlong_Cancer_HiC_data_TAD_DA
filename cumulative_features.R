# Rscript cumulative_features.R

script_name <- "cumulative_features.R"

startTime <- Sys.time()

cat("> START cumulative_features.R \n")

require(foreach)
require(doMC)

registerDoMC(40)

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

pipFolder <- file.path("PIPELINE", "OUTPUT_FOLDER")

cexAxis <- cexLab <- 1.2
plotType <- "png"
myHeight <- ifelse(plotType == "png", 400, 7)
myWidth <- ifelse(plotType == "png", 600, 10)

outFold <- "CUMULATIVE_FEATURES"
dir.create(outFold)


coexprFolder <- "CREATE_COEXPR_SORTNODUP"
corMet <- "pearson"

script0_name <- "0_prepGeneData"
script1_name <- "1_runGeneDE"
script3_name <- "3_runMeanTADLogFC"
script4_name <- "4_runMeanTADCorr"

##################### FOLD CHANGE
all_tadFC_files <- list.files(pipFolder, pattern = "all_meanLogFC_TAD.Rdata", full.names = TRUE, recursive = TRUE)
stopifnot(length(all_tadFC_files) > 0)

fc_file = all_tadFC_files[1]

foo <- foreach(fc_file = all_tadFC_files) %dopar% {
  
  hicds <- basename(dirname(dirname(dirname((fc_file)))))
  exprds <- basename(dirname(dirname(fc_file)))
  
  tad_fc <- eval(parse(text = load(fc_file)))
  
  DE_DT_file <- file.path(pipFolder, hicds, exprds, script1_name, "DE_topTable.Rdata")  
  stopifnot(file.exists(DE_DT_file))
  DE_DT <- eval(parse(text = load(DE_DT_file)))
  DE_DT$genes <- as.character(DE_DT$genes)
  
  geneList_file <- file.path(pipFolder, hicds, exprds, script0_name, "pipeline_geneList.Rdata")
  stopifnot(file.exists(geneList_file))
  geneList <- eval(parse(text = load(geneList_file)))
  
  stopifnot(rownames(DE_DT) == DE_DT$genes)
  stopifnot(names(geneList) %in% DE_DT$genes)  
  
  DE_DT <- DE_DT[DE_DT$genes %in% names(geneList),]
  
  gene_fc <- setNames(DE_DT$logFC, DE_DT$genes)
  
  tad_fc_sorted <- sort(tad_fc, decreasing=TRUE)
  gene_fc_sorted <- sort(gene_fc, decreasing=TRUE)
  
  tad_fc_sorted_norm <- tad_fc_sorted/max(tad_fc_sorted)
  gene_fc_sorted_norm <- gene_fc_sorted/max(gene_fc_sorted)
  
  tad_fc_cumsum <- cumsum(tad_fc_sorted)
  gene_fc_cumsum <- cumsum(gene_fc_sorted)
  
  tad_fc_cumsum_norm <- tad_fc_cumsum/max(cumsum(tad_fc_sorted))
  gene_fc_cumsum_norm <-  gene_fc_cumsum/max(cumsum(gene_fc_sorted))
  
  tad_fc_ranks <- c(1:length(tad_fc_cumsum))
  gene_fc_ranks <- c(1:length(gene_fc_cumsum))
  
  tad_fc_ranks_norm <- tad_fc_ranks/length(tad_fc_ranks)
  gene_fc_ranks_norm <- gene_fc_ranks/length(gene_fc_ranks)
  
  mySub <- paste0("nTADs = ", length(tad_fc_ranks_norm), " - nGenes = ", length(gene_fc_ranks_norm))
  
  outFile <- file.path(outFold, paste0(hicds, "_", exprds, "_", "cumulative_fc_tads.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  plot(
    x = tad_fc_ranks,
    y = tad_fc_cumsum,
    type="l",
    cex.axis = cexAxis,
    cex.lab = cexLab,
    xlab = "TAD ranks",
    ylab = "cumulative FC",
    main = paste0(hicds, " - ", exprds)
  )
  mtext(side=3, text=mySub, font=3)
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  outFile <- file.path(outFold, paste0(hicds, "_", exprds, "_", "cumulative_fc_genes.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  plot(
    x = gene_fc_ranks,
    y = gene_fc_cumsum,
    type="l",
    cex.axis = cexAxis,
    cex.lab = cexLab,
    xlab = "gene ranks",
    ylab = "cumulative FC",
    main = paste0(hicds, " - ", exprds)
  )
  mtext(side=3, text=mySub, font=3)
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  outFile <- file.path(outFold, paste0(hicds, "_", exprds, "_", "cumulative_fc_bothGenesTADs.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  plot(
    x = tad_fc_ranks_norm,
    y = tad_fc_cumsum_norm,
    type="l",
    cex.axis = cexAxis,
    cex.lab = cexLab,
    xlab = "relative rank",
    ylab = "relative FC cumsum",
    main = paste0(hicds, " - ", exprds)
  )
  mtext(side=3, text=mySub, font=3)
  points(
    x = gene_fc_ranks_norm,
    y = gene_fc_cumsum_norm,
    type = "l",
    col="red"
  )
  legend("topleft", lty=1, col = c("black", "red"), legend = c("TAD", "gene"), bty="n")
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  ############################################################################################# CORRELATION DATA
  
  meanCorr_file <- file.path(pipFolder, hicds, exprds, script4_name, "all_meanCorr_TAD.Rdata")
  stopifnot(file.exists(meanCorr_file))
  tad_corr <- eval(parse(text = load(meanCorr_file)))
  
  cat("... load coexprDT\n")
  
  coexprFile <- file.path(coexprFolder, hicds, exprds, corMet, "coexprDT.Rdata")
  stopifnot(file.exists(coexprFile))  
  coexprDT <- eval(parse(text = load(coexprFile)))
  gene_corr <- coexprDT$coexpr
  
  tad_corr_sorted <- sort(tad_corr, decreasing = TRUE)
  gene_corr_sorted <- sort(gene_corr, decreasing = TRUE)
  
  tad_corr_cumsum <- cumsum(tad_corr_sorted)
  gene_corr_cumsum <- cumsum(gene_corr_sorted)
  
  tad_corr_cumsum_norm <- tad_corr_cumsum/max(tad_corr_cumsum)
  gene_corr_cumsum_norm <- gene_corr_cumsum/max(gene_corr_cumsum)
  
  tad_corr_ranks <- 1:length(tad_corr_sorted)
  gene_corr_ranks <- 1:length(gene_corr_sorted)
  
  tad_corr_ranks_norm <- tad_corr_ranks/max(tad_corr_ranks)  
  gene_corr_ranks_norm <- gene_corr_ranks/max(gene_corr_ranks)
  
  mySub <- paste0("nTADs = ", length(tad_corr_ranks_norm), " - nGenes = ", length(gene_corr_ranks_norm))
  
  
  
  outFile <- file.path(outFold, paste0(hicds, "_", exprds, "_", "cumulative_corr_tads.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  plot(
    x = tad_corr_ranks,
    y = tad_corr_cumsum,
    type="l",
    cex.axis = cexAxis,
    cex.lab = cexLab,
    xlab = "TAD ranks",
    ylab = "cumulative corr.",
    main = paste0(hicds, " - ", exprds)
  )
  mtext(side=3, text=mySub, font=3)
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  outFile <- file.path(outFold, paste0(hicds, "_", exprds, "_", "cumulative_corr_genes.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  plot(
    x = gene_corr_ranks,
    y = gene_corr_cumsum,
    type="l",
    cex.axis = cexAxis,
    cex.lab = cexLab,
    xlab = "gene ranks",
    ylab = "cumulative corr.",
    main = paste0(hicds, " - ", exprds)
  )
  mtext(side=3, text=mySub, font=3)
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  outFile <- file.path(outFold, paste0(hicds, "_", exprds, "_", "cumulative_corr_bothGenesTADs.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  plot(
    x = tad_corr_ranks_norm,
    y = tad_corr_cumsum_norm,
    type="l",
    cex.axis = cexAxis,
    cex.lab = cexLab,
    xlab = "relative rank",
    ylab = "relative corr. cumsum",
    main = paste0(hicds, " - ", exprds)
  )
  mtext(side=3, text=mySub, font=3)
  points(
    x = gene_corr_ranks_norm,
    y = gene_corr_cumsum_norm,
    type = "l",
    col="red"
  )
  legend("topleft", lty=1, col = c("black", "red"), legend = c("TAD", "gene"), bty="n")
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
}


# ######################################################################################
# ######################################################################################
# ######################################################################################
cat("*** DONE ", script_name, "\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))



