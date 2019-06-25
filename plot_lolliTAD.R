#!/usr/bin/Rscript

startTime <- Sys.time()

script_name <- "plot_lolliTAD.R"
cat(paste0("> START ", script_name,  "\n"))

# Rscript plot_lolliTAD.R TCGAskcm_wt_mutCTNNB1 ENCSR312KHQ_SK-MEL-5_40kb chr1_TAD200


SSHFS <- FALSE
setDir <- ifelse(SSHFS, "/media/electron", "")
setDir <- ifelse(SSHFS, "/media/electron", "")

exprds <- "TCGAluad_mutKRAS_mutEGFR"
hicds <- "NCI-H460_40kb"
args <- commandArgs(trailingOnly = T)
stopifnot(length(args) > 2)
exprds <- args[1] 
hicds <- args[2] 
all_TADs <- args[3:length(args)]

mainDir <- "."
pipelineDir <- file.path(mainDir, "PIPELINE")
settingDir <- file.path(pipelineDir, "INPUT_FILES")
pipOutDir <-  file.path(pipelineDir, "OUTPUT_FOLDER")

pipOutDir <- file.path(pipOutDir, hicds, exprds)

stopifnot(dir.exists(hicds))
stopifnot(dir.exists(settingDir))
stopifnot(dir.exists(pipOutDir))

pipScriptDir <- paste0(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2")
caller <- "TopDom"
mainSettingsFile <- file.path(paste0(pipScriptDir, "_",  caller), "main_settings.R") # needed for entrezDT_file
cat(mainSettingsFile)
stopifnot(file.exists(mainSettingsFile))

outFold <- "PLOT_LOLLITAD"
dir.create(outFold, recursive=TRUE)

logfile <- file.path(outFold, paste0(exprds, "_", hicds, "_", "plot_lollitad_logfile.txt"))
if(!SSHFS) file.remove(logfile)
if(SSHFS) logfile <- "" 

# PIPELINE/INPUT_FILES/NCI-H460_40kb/run_settings_TCGAluad_mutKRAS_mutEGFR.R
settingF <- file.path(settingDir, hicds, paste0("run_settings_", exprds, ".R"))
stopifnot(file.exists(settingF))
  
script0_name <- "0_prepGeneData"
script1_name <- "1_runGeneDE"
script3_name <- "3_runMeanTADLogFC"


source(mainSettingsFile)
#source("run_settings.R")
source(paste0(pipScriptDir, "/", "TAD_DE_utils.R"))
source(settingF)
source(file.path("..", "Cancer_HiC_data_TAD_DA", "utils_fct.R"))
suppressPackageStartupMessages(library(grid, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 
suppressPackageStartupMessages(library(tools, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 

plotType <- "svg"
myHeight <- ifelse(plotType == "png", 480, 7)
myWidth <- ifelse(plotType == "png", 600, 10)

################################****************************************************************************************
####################################################### PREPARE INPUT
################################****************************************************************************************

# INPUT DATA
txt <- paste0("... gene2tadDT_file = ", gene2tadDT_file, "\n")
printAndLog(txt, logFile = logfile)
gene2tadDT <- read.delim(gene2tadDT_file, header=F, col.names = c("entrezID", "chromo", "start", "end", "region"), stringsAsFactors = F)
gene2tadDT$entrezID <- as.character(gene2tadDT$entrezID)

txt <- paste0("... entrezDT_file = ", entrezDT_file, "\n")
printAndLog(txt, logFile = logfile)
entrez2symbDT <- read.delim(entrezDT_file, header=T, stringsAsFactors=F)
entrez2symbDT <- entrez2symbDT[,c("entrezID", "symbol")]
colnames(entrez2symbDT) <- c("entrezID", "geneName")
entrez2symbDT$entrezID <- as.character(entrez2symbDT$entrezID)

meanTADlogFC <- eval(parse(text = load(paste0(pipOutDir, "/", script3_name, "/", "all_meanLogFC_TAD.Rdata"))))

########################################################################################
########### PREPARE THE GENES TO PLOT THE logFC 
########################################################################################
samp1 <- eval(parse(text = load(paste0(setDir, "/", sample1_file))))
samp2 <- eval(parse(text = load(paste0(setDir, "/", sample2_file))))

# UPDATE SELECT THE GENES ACCORDING TO THE SETTINGS PREPARED IN 0_PREPGENEDATA
rnaseqDT <- eval(parse(text = load(paste0(pipOutDir, "/", script0_name, "/rna_rnaseqDT.Rdata"))))
initList <- eval(parse(text = load(paste0(pipOutDir, "/", script0_name, "/rna_geneList.Rdata"))))
geneList <- eval(parse(text = load(paste0(pipOutDir, "/", script0_name, "/pipeline_geneList.Rdata"))))

txt <- paste0(toupper(script_name), "> Start with # genes: ", length(geneList), "/", length(initList), "\n")
printAndLog(txt, logfile)

rnaseqDT <- rnaseqDT[names(geneList),]    
stopifnot(all(rownames(rnaseqDT) == names(geneList)))
stopifnot(!any(duplicated(names(geneList))))
stopifnot(!any(duplicated(geneList)))

DE_topTable <- eval(parse(text = load(paste0(pipOutDir, "/", script1_name, "/DE_topTable.Rdata"))))
stopifnot(!any(duplicated(names(geneList))))
DE_topTable <- DE_topTable[DE_topTable$genes %in% names(geneList),]
stopifnot(nrow(DE_topTable) > 0)

gene2tadDT <- gene2tadDT[gene2tadDT$entrezID %in% geneList,]

DE_topTable$genes <- unlist(sapply(DE_topTable$genes, function(x) geneList[x]))
rownames(DE_topTable) <- NULL

# ! duplicated row.names are not allowed !
dupEntrez <- unique(geneList[duplicated(geneList)])
geneList <- geneList[! geneList %in% dupEntrez]
stopifnot(!any(duplicated(geneList)))

DE_topTable <- DE_topTable[!DE_topTable$genes %in% dupEntrez,]
stopifnot(!any(duplicated(DE_topTable$genes)))

initNrow <- nrow(rnaseqDT)
rnaseqDT <- rnaseqDT[which(rownames(rnaseqDT) %in% names(geneList)),]
txt <- paste0(toupper(script_name), "> Discard duplicated symbol, retain: ", nrow(rnaseqDT), "/", initNrow , " genes\n")
printAndLog(txt, logfile)


stopifnot(all(rownames(rnaseqDT) == names(geneList)))
stopifnot(is.numeric(rnaseqDT[1,1]))
#if(applyVoomAndCPM) {
if(inputDataType == "raw" | inputDataType == "RSEM") {
  log2_rnaseqDT <- log2(rnaseqDT + 0.0001) 
} else{
  log2_rnaseqDT <- rnaseqDT
}
meanExpr <- rowMeans(log2_rnaseqDT, na.rm=T)
names(meanExpr) <- unlist(sapply(names(meanExpr), function(x) geneList[x]))

################################****************************************************************************************
####################################################### DO THE PLOTS 
################################****************************************************************************************
###################################################################### PLOT LOLLI TAD
 
# retrieve which condition is cond1, i.e. the one that is more expressed when logFC is positive
geneHighestLogFC <- names(geneList[geneList == DE_topTable$genes[which.max(DE_topTable$logFC)] ])

samp1_vect <- rnaseqDT[geneHighestLogFC,samp1, drop=F]
stopifnot(dim(samp1_vect) == c(1,length(samp1)))
samp2_vect <- rnaseqDT[geneHighestLogFC,samp2, drop=F]
stopifnot(dim(samp2_vect) == c(1,length(samp2)))

plot_cond1 <- ifelse(as.numeric(rowMeans(samp1_vect, na.rm = T)) > as.numeric(rowMeans(samp2_vect, na.rm = T)), cond1, cond2)
plot_cond2 <- ifelse(as.numeric(rowMeans(samp1_vect, na.rm = T)) > as.numeric(rowMeans(samp2_vect, na.rm = T)), cond2, cond1)

n_topTAD_toplot <- length(all_TADs)

outWidth <- 20
outHeight <- min(c(7 * n_topTAD_toplot/2, 49))

if(outHeight < 7) outHeight <- 7

vect_plot <- list()
for(i_plot in 1:n_topTAD_toplot)  {
  tad_to_plot <- all_TADs[i_plot]
  vect_plot[[i_plot]] <- plot_lolliTAD(TAD_to_plot = tad_to_plot,
                                       meanExpr_vect = meanExpr, 
                                       DE_table = DE_topTable,
                                       g2t_table = gene2tadDT, 
                                       id2name_table=entrez2symbDT, 
                                       geneList = geneList,
                                       textLeft =  meanTADlogFC[tad_to_plot] > 0,
                                       orderBy = "logFC", cond1=plot_cond1, cond2=plot_cond2)
}
all_plots <- do.call(grid.arrange, c(vect_plot,  list(ncol=2, top=textGrob(paste(exprds),
                                                                           gp=gpar(fontsize=20,font=2)))))
# cat("myHeight =", outHeight, "\n")

outFile <- file.path(outFold, paste0(exprds, "_", hicds, "_", paste0(all_TADs, collapse = "_"), ".", plotType))
ggsave(filename = outFile, all_plots, width=outWidth, height = outHeight)
cat("... written: ", outFile, "\n")


###################################################################### 
txt <- paste0(startTime, "\n", Sys.time(), "\n")
printAndLog(txt, logfile)
cat(paste0("*** DONE: ", script_name, "\n"))


