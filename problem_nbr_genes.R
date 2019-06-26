# Rscript problem_nbr_genes.R

script_name <- "problem_nbr_genes.R"

startTime <- Sys.time()

cat("> START coexpr_DE_queryTAD.R \n")

require(foreach)
script0_name <- "0_prepGeneData"


outFold <- "PROBLEM_NBR_GENES"
dir.create(outFold)

plotType <- "png"
myWidth <- 400
myHeight <- 400
cexPlot <- 1.2

pipDir <- file.path("PIPELINE", "OUTPUT_FOLDER")
mainDir_0 <- file.path("..", "Cancer_HiC_data_TAD_DA") 
pipDir_0 <- file.path(mainDir_0, "PIPELINE", "OUTPUT_FOLDER")

g2tFolder <- "genes2tad"
g2tSuffix <-  "all_genes_positions.txt"
posSuffix <- "all_assigned_regions.txt"

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

all_geneList_files <- list.files(paste0(pipDir), recursive = TRUE, pattern="pipeline_geneList.Rdata", full.names = TRUE)
stopifnot(length(all_geneList_files) > 0)

g_file = all_geneList_files[1]
# all_geneList_files = all_geneList_files[1]
genes_nbr_DT <- foreach(g_file = all_geneList_files,.combine='rbind') %do% {
  curr_file <- g_file
  
  # stopifnot(file.exists(curr_file))
  curr_hicds <- basename(dirname(dirname(dirname(g_file))))
  curr_exprds <- basename(dirname(dirname(g_file)))
  
  cat(curr_file, "\n")
  
  gList <- eval(parse(text = load(curr_file)))
  
  curr_file <- file.path(pipDir, curr_hicds, curr_exprds, script0_name, "pipeline_regionList.Rdata")
  # stopifnot(file.exists(curr_file))
  cat(curr_file, "\n")
  
  rList <- eval(parse(text = load(curr_file)))
  
  if(curr_hicds == "ENCSR489OCU_NCI-H460_40kb") {
    curr_hicds_0 <- "NCI-H460_40kb"
  } else if(curr_hicds == "GSE75070_MCF-7_shNS_40kb"){
    curr_hicds_0 <- "MCF-7_40kb"
  } else if(curr_hicds == "GSE118514_RWPE1_40kb" | curr_hicds == "GSE58752_liver_40kb") {
    return(NULL)
  } else {
    curr_hicds_0 <- curr_hicds
  }
  
  curr_file_0 <- file.path(pipDir_0, curr_hicds_0, curr_exprds, script0_name, "pipeline_geneList.Rdata")
  # stopifnot(file.exists(curr_file_0))
  gList_0 <- eval(parse(text = load(curr_file_0)))
  
  curr_file_0 <- file.path(pipDir_0, curr_hicds_0, curr_exprds, script0_name, "pipeline_regionList.Rdata")
  # stopifnot(file.exists(curr_file_0))
  rList_0 <- eval(parse(text = load(curr_file_0)))
  
  
  
  g2tFile <- file.path(curr_hicds, g2tFolder, g2tSuffix)
  stopifnot(file.exists(g2tFile))
  g2t_DT <- read.delim(g2tFile, header=F, col.names = c("entrezID",  "chromo", "start", "end", "region"), stringsAsFactors = FALSE)
  g2t_DT$entrezID <- as.character(g2t_DT$entrezID)
  g2t_nbr <- sum(grepl("_TAD", g2t_DT$region))
  
  
  g2tFile_0 <- file.path(mainDir_0, curr_hicds_0, g2tFolder, g2tSuffix)
  stopifnot(file.exists(g2tFile_0))
  g2t_DT_0 <- read.delim(g2tFile_0, header=F, col.names = c("entrezID",  "chromo", "start", "end", "region"), stringsAsFactors = FALSE)
  g2t_DT_0$entrezID <- as.character(g2t_DT_0$entrezID)
  g2t_nbr_0 <- sum(grepl("_TAD", g2t_DT_0$region))
  
  xGenes <- sum(g2t_DT$chromo == "chrX" & grepl("_TAD", g2t_DT$region))
  xGenes_0 <- sum(g2t_DT_0$chromo == "chrX" & grepl("_TAD", g2t_DT_0$region))
  
  
  posFile <- file.path(curr_hicds, g2tFolder, posSuffix)
  stopifnot(file.exists(posFile))
  pos_DT <- read.delim(posFile, header=F, 
                          col.names=c("chromo", "region", "start", "end"), stringsAsFactors = FALSE)
  stopifnot(is.numeric(pos_DT$start))
  stopifnot(is.numeric(pos_DT$end))
  tadpos_DT <- pos_DT[grepl("_TAD", pos_DT$region),]
  tadpos_nbr <- sum(grepl("_TAD", tadpos_DT$region))
  stopifnot(tadpos_nbr == nrow(tadpos_DT))
  tadpos_DT$tadSize <- tadpos_DT$end - tadpos_DT$start + 1
  tadpos_cvr <- sum(tadpos_DT$tadSize)
  
  posFile_0 <- file.path(mainDir_0, curr_hicds_0, g2tFolder, posSuffix)
  stopifnot(file.exists(posFile_0))
  pos_DT_0 <- read.delim(posFile_0, header=F, 
                       col.names=c("chromo", "region", "start", "end"), stringsAsFactors = FALSE)
  stopifnot(is.numeric(pos_DT_0$start))
  stopifnot(is.numeric(pos_DT_0$end))
  tadpos_DT_0 <- pos_DT_0[grepl("_TAD", pos_DT_0$region),]
  tadpos_nbr_0 <- sum(grepl("_TAD", tadpos_DT_0$region))
  stopifnot(tadpos_nbr_0 == nrow(tadpos_DT_0))
  tadpos_DT_0$tadSize <- tadpos_DT_0$end - tadpos_DT_0$start + 1
  tadpos_cvr_0 <- sum(tadpos_DT_0$tadSize)
  
  chr <- unique(g2t_DT$chromo[grep("_TAD", g2t_DT$region)])
  chr_0 <- unique(g2t_DT_0$chromo[grep("_TAD", g2t_DT_0$region)])
  
  missingChr <- paste0(chr_0[!chr_0 %in% chr], collapse = ",")
  
  data.frame(
    hicds = curr_hicds,
    exprds = curr_exprds,
    nGenes = length(gList),
    nGenes_0 = length(gList_0),
    nRegions = length(rList),
    nRegions_0 = length(rList_0),
    g2t_nbr = g2t_nbr,
    g2t_nbr_0 = g2t_nbr_0,
    tadpos_nbr = tadpos_nbr,
    tadpos_nbr_0 = tadpos_nbr_0,
    tadpos_cvr = tadpos_cvr,
    tadpos_cvr_0 = tadpos_cvr_0,
    missingChr = missingChr,
    xGenes = xGenes,
    xGenes_0 = xGenes_0,
    stringsAsFactors = FALSE
  )
}

genes_nbr_DT$nGenes_Diff <- genes_nbr_DT$nGenes_0 - genes_nbr_DT$nGenes
genes_nbr_DT$coverDiff <- genes_nbr_DT$tadpos_cvr_0 - genes_nbr_DT$tadpos_cvr

######################################################################### # genes diff vs. # genes chrX
xvar <- "nGenes_Diff"
yvar <- "xGenes_0"
outFile <- file.path(outFold, paste0(yvar, "_", xvar, ".", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot(x = genes_nbr_DT[,xvar],
         y = genes_nbr_DT[,yvar],
         xlab = paste0(xvar),
         ylab = paste0(yvar),
         cex.axis = cexPlot,
         cex.lab = cexPlot,
     pch=16, cex = 0.7
         )
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


######################################################################### # genes vs. # genes v0
xvar <- "nGenes"
yvar <- "nGenes_0"
outFile <- file.path(outFold, paste0(yvar, "_", xvar, ".", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot(x = genes_nbr_DT[,xvar],
     y = genes_nbr_DT[,yvar],
     xlab = paste0(xvar),
     ylab = paste0(yvar),
     cex.axis = cexPlot,
     cex.lab = cexPlot,
     main  = "# pipeline genes",
     pch=16, cex = 0.7
)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

outFile <- file.path(outFold, paste0(yvar, "_", xvar, "_densplot.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot_multiDens(list(
     nGenes = genes_nbr_DT[,xvar],
     nGenes_0 = genes_nbr_DT[,yvar]),
     my_xlab = paste0(xvar),
     # my_ylab = paste0(yvar),
     plotTit =  "# pipeline genes"
)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

######################################################################### # genes vs. # assigned genes v0
xvar <- "g2t_nbr"
yvar <- "g2t_nbr_0"
outFile <- file.path(outFold, paste0(yvar, "_", xvar, ".", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot(x = genes_nbr_DT[,xvar],
     y = genes_nbr_DT[,yvar],
     xlab = paste0(xvar),
     ylab = paste0(yvar),
     cex.axis = cexPlot,
     cex.lab = cexPlot,
     main = "# assigned genes",
     pch=16, cex = 0.7
)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

outFile <- file.path(outFold, paste0(yvar, "_", xvar, "_densplot.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot_multiDens(list(
     g2t_nbr = genes_nbr_DT[,xvar],
     g2t_nbr_0 = genes_nbr_DT[,yvar]),
     my_xlab = paste0(xvar),
     # my_ylab = paste0(yvar),
     plotTit =  "# assigned genes"
)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


######################################################################### TAD cvr vs. TAD cvr v0
xvar <- "tadpos_cvr"
yvar <- "tadpos_cvr_0"
outFile <- file.path(outFold, paste0(yvar, "_", xvar, ".", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot(x = log10(genes_nbr_DT[,xvar]),
     y = log10(genes_nbr_DT[,yvar]),
     xlab = paste0(xvar, " (log10)"),
     ylab = paste0(yvar, "(log10)"),
     cex.axis = cexPlot,
     cex.lab = cexPlot,
     pch=16, cex = 0.7
)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

outFile <- file.path(outFold, paste0(yvar, "_", xvar, "_densplot.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot_multiDens(list(
  tadpos_cvr =  log10(genes_nbr_DT[,xvar]),
  tadpos_cvr_0 = log10(genes_nbr_DT[,yvar])
  ),
    my_xlab = paste0(xvar, " (log10)")
     # my_ylab = paste0(yvar, "(log10)")
)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

######################################################################### # assigned TAD vs. # assigned TAD v0
xvar <- "tadpos_nbr"
yvar <- "tadpos_nbr_0"
outFile <- file.path(outFold, paste0(yvar, "_", xvar, ".", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot(x = log10(genes_nbr_DT[,xvar]),
     y = log10(genes_nbr_DT[,yvar]),
     xlab = paste0(xvar, " (log10)"),
     ylab = paste0(yvar, "(log10)"),
     cex.axis = cexPlot,
     cex.lab = cexPlot,
     pch=16, cex = 0.7
)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

outFile <- file.path(outFold, paste0(yvar, "_", xvar, ".", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot_multiDens(list(
     tadpos_nbr = log10(genes_nbr_DT[,xvar]),
     tadpos_nbr_0 = log10(genes_nbr_DT[,yvar])
     ),
     my_xlab = paste0(xvar, " (log10)")
     # my_ylab = paste0(yvar, "(log10)")
)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

######################################################################### # pipeline TAD vs. # pipeline TAD v0
xvar <- "nRegions"
yvar <- "nRegions_0"
outFile <- file.path(outFold, paste0(yvar, "_", xvar, ".", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot(x = log10(genes_nbr_DT[,xvar]),
     y = log10(genes_nbr_DT[,yvar]),
     xlab = paste0(xvar, " (log10)"),
     ylab = paste0(yvar, "(log10)"),
     cex.axis = cexPlot,
     cex.lab = cexPlot,
     pch=16, cex = 0.7
)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


outFile <- file.path(outFold, paste0(yvar, "_", xvar, "_densplot.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot_multiDens(list(
     nRegions = log10(genes_nbr_DT[,xvar]),
     nRegions_0 = log10(genes_nbr_DT[,yvar])),
     my_xlab = paste0(xvar, " (log10)")
     # my_ylab = paste0(yvar, "(log10)")
)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))




# ######################################################################################
# ######################################################################################
# ######################################################################################
cat(paste0("*** DONE: ", script_name, "\n"))
cat(paste0(startTime, "\n", Sys.time(), "\n"))





