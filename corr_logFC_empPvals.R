
# Rscript corr_logFC_empPvals.R 

startTime <- Sys.time()

cat("> START corr_logFC_empPvals.R \n")

SSHFS <- FALSE

require(foreach)
require(doMC)

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

plotType <- "png"
myHeight <- ifelse(plotType=="png", 500, 7)
myWidth <- myHeight
myHeightDensity <- myHeight
myWidthDensity <- ifelse(plotType=="png", 600, 10)

myWidthBarplot <- ifelse(plotType=="png", 600, 10)
myHeightBarplot <- ifelse(plotType=="png", 400, 7)

plotCex <- 1.4

registerDoMC(ifelse(SSHFS, 2, 40))

build_signifTADs_allDS_data <- TRUE

setDir <- ifelse(SSHFS, "/media/electron", "")

outFolder <- file.path("CORR_LOGFC_EMPPVALS")
dir.create(outFolder, recursive = TRUE)

checkFile <- file.path(outFolder, paste0("check_file_traceback", ".txt"))
file.remove(checkFile)

logFile=""

nPermut <- 10000
minEmpPval <- 1/(nPermut+1)
minEmpPval

# stopifnot(!is.na(signifThresh))
# stopifnot(signifThresh > 0)

#PIPELINE/OUTPUT_FOLDER/GSE105318_DLD1_40kb/TCGAcoad_msi_mss//emp_pval_combined.Rdata

script0_name <- "0_prepGeneData"
script3_name <- "3_runMeanTADLogFC"
script4_name <- "4_runMeanTADCorr"
script8c_name <- "8c_runAllDown"
script9_name <- "9_runEmpPvalMeanTADLogFC"
script9v2_name <- "9v2_runEmpPvalWilcoxStat"
script10_name <- "10_runEmpPvalMeanTADCorr"
script10v2_name <- "10v2_runEmpPvalMeanTADCorr"
script10b_name <- "10b_runEmpPvalProdSignedRatio"
script11_name <- "11_runEmpPvalCombined"

pipOutFolder <- file.path("PIPELINE", "OUTPUT_FOLDER")

all_hicexpr_ds <- unname(unlist(sapply(list.files(pipOutFolder, full.names = TRUE), 
                                       function(x) file.path(basename(x),  list.files(x)))))


all_hicexpr_ds <- list.files(pipOutFolder, recursive = TRUE, pattern="emp_pval_combined.Rdata", full.names = FALSE)
all_hicexpr_ds <- dirname(dirname(all_hicexpr_ds))

stopifnot(dir.exists(file.path(pipOutFolder, all_hicexpr_ds)))

ds=all_hicexpr_ds[1]

# all_hicexpr_ds=all_hicexpr_ds[1]

if(build_signifTADs_allDS_data){
  cat("... start building allPvals_allDS_DT data \n")
  
  allPvals_allDS_DT <- foreach(ds = all_hicexpr_ds, .combine='rbind') %dopar% {
    cat("... start: ", ds, "\n")
    
    hicds <- dirname(ds)
    exprds <- basename(ds)
    stopifnot(dir.exists(hicds))
    dsPipOutDir <- file.path(pipOutFolder, ds)
    stopifnot(dir.exists(dsPipOutDir))
    
    g2tFile <- file.path(hicds, "genes2tad", "all_genes_positions.txt")
    stopifnot(file.exists(g2tFile))
    g2t_DT <- read.delim(g2tFile, header=F, 
                         col.names = c("entrezID",  "chromo", "start", "end", "region"), stringsAsFactors = FALSE)
    g2t_DT$entrezID <- as.character(g2t_DT$entrezID)
    
    # RETRIEVE gene list - script0
    stopifnot(dir.exists(file.path(dsPipOutDir, script0_name)))
    geneListFile <- file.path(dsPipOutDir, script0_name, "pipeline_geneList.Rdata")
    stopifnot(file.exists(geneListFile))
    pipeline_geneList <- eval(parse(text = load(geneListFile))) # not adjusted
    
    regionListFile <- file.path(dsPipOutDir, script0_name, "pipeline_regionList.Rdata")
    stopifnot(file.exists(regionListFile))
    pipeline_regionList <- eval(parse(text = load(regionListFile))) # not adjusted
    
    stopifnot(pipeline_geneList %in% g2t_DT$entrezID)
    g2t_DT <- g2t_DT[g2t_DT$entrezID %in% pipeline_geneList,]
    stopifnot(length(pipeline_geneList) == nrow(g2t_DT))
    tad_nGenes <- setNames(as.numeric(table(g2t_DT$region)), as.character(names(table(g2t_DT$region))))
    
    # RETRIEVE logFC pval - script9
    stopifnot(dir.exists(file.path(dsPipOutDir, script9_name)))
    tad_pvalFCFile <- file.path(dsPipOutDir, script9_name, "emp_pval_meanLogFC.Rdata")
    stopifnot(file.exists(tad_pvalFCFile))
    tad_pvalFC <- eval(parse(text = load(tad_pvalFCFile))) # not adjusted
    adj_tad_pvalFC <- sort(p.adjust(tad_pvalFC, method="BH"))
    stopifnot(setequal(names(adj_tad_pvalFC), pipeline_regionList))
        
    # RETRIEVE logFC values - script3
    stopifnot(dir.exists(file.path(dsPipOutDir, script3_name)))
    tad_valuesFCFile <- file.path(dsPipOutDir, script3_name, "all_meanLogFC_TAD.Rdata")
    stopifnot(file.exists(tad_valuesFCFile))
    tad_valuesFC <- eval(parse(text = load(tad_valuesFCFile))) # not adjusted
    tad_valuesFC <- sort(tad_valuesFC)
    stopifnot(setequal(names(tad_valuesFC), pipeline_regionList))
    
    # RETRIEVE TADcorr pval - script10
    stopifnot(dir.exists(file.path(dsPipOutDir, script10_name)))
    tad_pvalCorrFile <- file.path(dsPipOutDir, script10_name, "emp_pval_meanCorr.Rdata")
    stopifnot(file.exists(tad_pvalCorrFile))
    tad_pvalCorr <- eval(parse(text = load(tad_pvalCorrFile)))
    adj_tad_pvalCorr <- sort(p.adjust(tad_pvalCorr, method="BH"))
    tad_pvalCorr <- sort(tad_pvalCorr)
    stopifnot(setequal(names(tad_pvalCorr), pipeline_regionList))
    
    # RETRIEVE TADcorr values - script4
    stopifnot(dir.exists(file.path(dsPipOutDir, script4_name)))
    tad_valuesCorrFile <- file.path(dsPipOutDir, script4_name, "all_meanCorr_TAD.Rdata")
    stopifnot(file.exists(tad_valuesCorrFile))
    tad_valuesCorr <- eval(parse(text = load(tad_valuesCorrFile)))
    tad_valuesCorr <- sort(tad_valuesCorr)
    stopifnot(setequal(names(tad_valuesCorr), pipeline_regionList))
    
    # RETRIEVE TADfcc values - script8c
    stopifnot(dir.exists(file.path(dsPipOutDir, script8c_name)))
    tad_valuesFCCfile <- file.path(dsPipOutDir, script8c_name, "all_obs_prodSignedRatio.Rdata")
    stopifnot(file.exists(tad_valuesFCCfile))
    tad_valuesFCC <- eval(parse(text = load(tad_valuesFCCfile)))
    tad_valuesFCC <- sort(tad_valuesFCC)
    stopifnot(setequal(names(tad_valuesFCC), pipeline_regionList))
    
    # # RETRIEVE SIGNIF. TADs comb pval
    stopifnot(dir.exists(file.path(dsPipOutDir, script11_name)))
    tad_pvalFile <- file.path(dsPipOutDir, script11_name, "emp_pval_combined.Rdata")
    stopifnot(file.exists(tad_pvalFile))
    tad_pvals <- eval(parse(text = load(tad_pvalFile)))
    adj_tad_pvalComb <- sort(p.adjust(tad_pvals, method="BH"))
    tad_pvalComb <- sort(tad_pvals)
    stopifnot(setequal(names(tad_pvalComb), pipeline_regionList))
    
    # RETRIEVE PVAL WILCOX
    stopifnot(dir.exists(file.path(dsPipOutDir, script9v2_name)))
    tad_pvalWilcoxFile <- file.path(dsPipOutDir, script9v2_name, "emp_pval_wilcoxStat.Rdata")
    stopifnot(file.exists(tad_pvalWilcoxFile))
    tad_pvalWilcox <- eval(parse(text = load(tad_pvalWilcoxFile))) # not adjusted
    adj_tad_pvalWilcox <- sort(p.adjust(tad_pvalWilcox, method="BH"))
    tad_pvalWilcox <- sort(tad_pvalWilcox)
    stopifnot(setequal(names(tad_pvalWilcox), pipeline_regionList))
    stopifnot(setequal(names(adj_tad_pvalWilcox), pipeline_regionList))
    
    
    all_tads <- pipeline_regionList
    
    outDT <- data.frame(
      hicds=hicds,
      exprds=exprds,
      tad=all_tads,
      tad_id=paste(hicds, exprds, all_tads, sep="_"),
      
      nGenes = tad_nGenes[all_tads],
      
      
      
      adj_pvalFC = as.numeric(adj_tad_pvalFC[all_tads]),
      adj_pvalCorr = as.numeric(adj_tad_pvalCorr[all_tads]),
      adj_pvalComb = as.numeric(adj_tad_pvalComb[all_tads]),

      adj_pvalWilcox = as.numeric(adj_tad_pvalWilcox[all_tads]),
      
      
      valuesFC = as.numeric(tad_valuesFC[all_tads]),
      valuesCorr = as.numeric(tad_valuesCorr[all_tads]),
      valuesFCC = as.numeric(tad_valuesFCC[all_tads]),
      
      stringsAsFactors = FALSE
    )
    outDT
  }
  outFile <- file.path(outFolder, "allPvals_allDS_DT.Rdata")
  save(allPvals_allDS_DT, file = outFile)
  cat(paste0("... written: ", outFile, "\n"))
  
} else { # if not to build
  outFile <- file.path(outFolder, "allPvals_allDS_DT.Rdata")
  allPvals_allDS_DT <- eval(parse(text = load(outFile)))
}

head(allPvals_allDS_DT)

all_ds <- unique(paste0(allPvals_allDS_DT$hicds, "_", allPvals_allDS_DT$exprds))


# ######################################################################################
# ###################################################################################### ajd pval Comb vs adj pval Corr
# ###################################################################################### 


myTit <- paste0("adj. empPval. meanIntraTADcorr vs. adj. empPval. Comb.")
mySub <- paste0("(nDS = ", length(all_ds), " - YL TADs)")

var1 <- "adj_pvalComb"
var2 <- "adj_pvalCorr"

myx <- allPvals_allDS_DT[,var1]
myy <- allPvals_allDS_DT[,var2]

outFile <- file.path(outFolder, paste0(var2, "_", var1, "_", "densplot.", plotType))
do.call(plotType, list(outFile,  height=myHeight, width=myWidth))
densplot(x = myx,
         y = myy,
         main = myTit,
         xlab = var1,
         ylab = var2,
         cex.axis = plotCex, cex.lab=plotCex, cex.main=plotCex)
mtext(text=mySub, side=3)
addCorr(x=myx, legPos="topleft",
        y=myy, bty='n')
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

myx <- -log10(allPvals_allDS_DT[,var1])
myy <- -log10(allPvals_allDS_DT[,var2])

outFile <- file.path(outFolder, paste0(var2, "_", var1, "_", "log10_densplot.", plotType))
do.call(plotType, list(outFile,  height=myHeight, width=myWidth))
densplot(x = myx,
         y = myy,
         main = myTit,
         xlab = paste0(var1, "[-log10]"),
         ylab = paste0(var2, "[-log10]"),
         cex.axis = plotCex, cex.lab=plotCex, cex.main=plotCex)
mtext(text=mySub, side=3)
addCorr(x=myx, legPos="topleft",
        y=myy, bty='n')
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

# ######################################################################################
# ######################################################################################  ajd pval Corr vs adj pval FC
# ###################################################################################### 

myTit <- paste0("adj. empPval. meanFC vs. adj. empPval. meanIntraTADcorr")
mySub <- paste0("(nDS = ", length(all_ds), " - YL TADs)")
var1 <- "adj_pvalCorr"
var2 <- "adj_pvalFC"

myx <- (allPvals_allDS_DT[,var1])
myy <- (allPvals_allDS_DT[,var2])

outFile <- file.path(outFolder, paste0(var2, "_", var1, "_", "densplot.", plotType))
do.call(plotType, list(outFile,  height=myHeight, width=myWidth))
densplot(x = myx,
         y = myy,
         main = myTit,
         xlab = var1,
         ylab = var2,
         cex.axis = plotCex, cex.lab=plotCex, cex.main=plotCex)
mtext(text=mySub, side=3)
addCorr(x=myx, legPos="topleft",
        y=myy, bty='n')
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))
  
myx <- -log10(allPvals_allDS_DT[,var1])
myy <- -log10(allPvals_allDS_DT[,var2])


outFile <- file.path(outFolder, paste0(var2, "_", var1, "_log10_", "densplot.", plotType))
do.call(plotType, list(outFile,  height=myHeight, width=myWidth))
densplot(x = myx,
         y = myy,
         main = myTit,
         xlab = paste0(var1, "[-log10]"),
         ylab = paste0(var2, "[-log10]"),
         cex.axis = plotCex, cex.lab=plotCex, cex.main=plotCex)
mtext(text=mySub, side=3)
addCorr(x=myx, legPos="topleft",
        y=myy, bty='n')
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


########################################################################################### 
########################################################################################### adj pvalCorr <-> adj pvalFC
########################################################################################### 


myTit <- paste0("adj. empPval. wilcox vs. adj. empPval. Corr.")
mySub <- paste0("(nDS = ", length(all_ds), ")")
var1 <- "adj_pvalCorr"
var2 <- "adj_pvalWilcox"
outFile <- file.path(outFolder, paste0(var2, "_", var1, "_", "densplot.", plotType))
do.call(plotType, list(outFile,  height=myHeight, width=myWidth))
myx <- (allPvals_allDS_DT[,var1]) 
myy <- (allPvals_allDS_DT[,var2]) 
densplot(x = myx,
         y = myy,
         main = myTit,
         xlab = var1,
         ylab = var2,
         cex.axis = plotCex, cex.lab=plotCex, cex.main=plotCex)
mtext(text=mySub, side=3)
addCorr(x=myx, legPos="topleft",
        y=myy, bty='n')
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

outFile <- file.path(outFolder, paste0(var2, "_", var1, "_log10_", "densplot.", plotType))
do.call(plotType, list(outFile,  height=myHeight, width=myWidth))
myx <- -log10(allPvals_allDS_DT[,var1]) 
myy <- -log10(allPvals_allDS_DT[,var2]) 
densplot(x = myx,
         y = myy,
         main = myTit,
         xlab = paste0(var1, "[-log10]"),
         ylab = paste0(var2, "[-log10]"),
         cex.axis = plotCex, cex.lab=plotCex, cex.main=plotCex)
mtext(text=mySub, side=3)
addCorr(x=myx, legPos="topleft",
        y=myy, bty='n')
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))



# ######################################################################################
# ###################################################################################### ajd pval Comb vs adj pval FC
# ###################################################################################### 

myTit <- paste0("adj. empPval. meanFC vs. adj. empPval. Comb.")
mySub <- paste0("(nDS = ", length(all_ds), " - YL TADs)")
var1 <- "adj_pvalComb"
var2 <- "adj_pvalFC"

myx <- (allPvals_allDS_DT[,var1])
myy <- (allPvals_allDS_DT[,var2])

outFile <- file.path(outFolder, paste0(var2, "_", var1, "_", "densplot.", plotType))
do.call(plotType, list(outFile,  height=myHeight, width=myWidth))
densplot(x = myx,
         y = myy,
         main = myTit,
         xlab = var1,
         ylab = var2,
         cex.axis = plotCex, cex.lab=plotCex, cex.main=plotCex)
mtext(text=mySub, side=3)
addCorr(x=myx, legPos="topleft",
        y=myy, bty='n')
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

myx <- -log10(allPvals_allDS_DT[,var1])
myy <- -log10(allPvals_allDS_DT[,var2])


outFile <- file.path(outFolder, paste0(var2, "_", var1, "_log10_", "densplot.", plotType))
do.call(plotType, list(outFile,  height=myHeight, width=myWidth))
densplot(x = myx,
         y = myy,
         main = myTit,
         xlab = paste0(var1, "[-log10]"),
         ylab = paste0(var2, "[-log10]"),
         cex.axis = plotCex, cex.lab=plotCex, cex.main=plotCex)
mtext(text=mySub, side=3)
addCorr(x=myx, legPos="topleft",
        y=myy, bty='n')
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))



# ######################################################################################
# ######################################################################################
# ######################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))

