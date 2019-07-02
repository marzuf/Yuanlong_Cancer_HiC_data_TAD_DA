
# Rscript cmp_empPvals_all.R 

startTime <- Sys.time()

cat("> START cmp_empPvals_all.R \n")

SSHFS <- FALSE

require(foreach)
require(doMC)

source("utils_fct.R")

plotType <- "png"
myHeight <- ifelse(plotType=="png", 500, 7)
myWidth <- myHeight
myHeightDensity <- myHeight
myWidthDensity <- ifelse(plotType=="png", 600, 10)

myWidthBarplot <- ifelse(plotType=="png", 600, 10)
myHeightBarplot <- ifelse(plotType=="png", 400, 7)

plotCex <- 1.4

nColBreaks <- 10

registerDoMC(ifelse(SSHFS, 2, 40))

build_signifTADs_allDS_data <- TRUE

setDir <- ifelse(SSHFS, "/media/electron", "")

outFolder <- file.path("CMP_EMPPVALS_ALL_2")
dir.create(outFolder, recursive = TRUE)

# checkFile <- file.path(outFolder, paste0("check_file_traceback_", toptad_id, "_", bottad_id, ".txt"))
# file.remove(checkFile)
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

tad_match_folder <- file.path("INTERSECT_topTADs_ACROSSDS", "top3")
stopifnot(dir.exists(tad_match_folder))

all_bestMatchDT_file <- file.path(tad_match_folder, "all_bestMatchDT.Rdata")
stopifnot(file.exists(all_bestMatchDT_file))

signifTADs_allDS_data_file <- file.path(tad_match_folder, "signifTADs_allDS_data.Rdata")
stopifnot(file.exists(signifTADs_allDS_data_file))

all_matchDT_file <- file.path(tad_match_folder, "all_matchDT.Rdata")
stopifnot(file.exists(all_matchDT_file))

ratio_matchingSignifTAD_DT_file <- file.path(tad_match_folder, "ratio_matchingSignifTAD_DT.Rdata")
stopifnot(file.exists(ratio_matchingSignifTAD_DT_file))

hicds_exprds_asMatch_DT_file <- file.path(tad_match_folder, "hicds_exprds_asMatch_DT.Rdata")
stopifnot(file.exists(hicds_exprds_asMatch_DT_file))


pipOutFolder <- file.path("PIPELINE", "OUTPUT_FOLDER")

all_hicexpr_ds <- unname(unlist(sapply(list.files(pipOutFolder, full.names = TRUE), function(x) file.path(basename(x),  list.files(x)))))
stopifnot(dir.exists(file.path(pipOutFolder, all_hicexpr_ds)))

ds=all_hicexpr_ds[1]

ds="ENCSR862OGI_RPMI-7951_40kb/TCGAskcm_lowInf_highInf"

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
    g2t_DT <- read.delim(g2tFile, header=F, col.names = c("entrezID",  "chromo", "start", "end", "region"), stringsAsFactors = FALSE)
    g2t_DT$entrezID <- as.character(g2t_DT$entrezID)
    
    # RETRIEVE gene list - script0
    stopifnot(dir.exists(file.path(dsPipOutDir, script0_name)))
    geneListFile <- file.path(dsPipOutDir, script0_name, "pipeline_geneList.Rdata")
    stopifnot(file.exists(geneListFile))
    pipeline_geneList <- eval(parse(text = load(geneListFile))) # not adjusted
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
    tad_pvalFC <- sort(tad_pvalFC)
    
    # RETRIEVE logFC values - script3
    stopifnot(dir.exists(file.path(dsPipOutDir, script3_name)))
    tad_valuesFCFile <- file.path(dsPipOutDir, script3_name, "all_meanLogFC_TAD.Rdata")
    stopifnot(file.exists(tad_valuesFCFile))
    tad_valuesFC <- eval(parse(text = load(tad_valuesFCFile))) # not adjusted
    tad_valuesFC <- sort(tad_valuesFC)
    
    
    # RETRIEVE Wilcox pval - script9v2
    #>>>>>>>>>>>> AT THE END; WILL NEED TO ADD
    # stopifnot(dir.exists(file.path(dsPipOutDir, script9v2_name)))
    tad_pvalWilcoxFile <- file.path(dsPipOutDir, script9v2_name, "emp_pval_wilcoxStat.Rdata")
    # TEMP! : might not all exist for the moment
    if(file.exists(tad_pvalWilcoxFile)){
      stopifnot(file.exists(tad_pvalWilcoxFile))
      tad_pvalWilcox <- eval(parse(text = load(tad_pvalWilcoxFile))) # not adjusted
      adj_tad_pvalWilcox <- sort(p.adjust(tad_pvalWilcox, method="BH"))
      tad_pvalWilcox <- sort(tad_pvalWilcox)
    } else {
      adj_tad_pvalWilcox <- NA
      tad_pvalWilcox <- NA
    }
    
    # RETRIEVE TADcorr pval - script10
    stopifnot(dir.exists(file.path(dsPipOutDir, script10_name)))
    tad_pvalCorrFile <- file.path(dsPipOutDir, script10_name, "emp_pval_meanCorr.Rdata")
    stopifnot(file.exists(tad_pvalCorrFile))
    tad_pvalCorr <- eval(parse(text = load(tad_pvalCorrFile)))
    adj_tad_pvalCorr <- sort(p.adjust(tad_pvalCorr, method="BH"))
    tad_pvalCorr <- sort(tad_pvalCorr)
    
    # RETRIEVE TADcorr pval - script10b
    #>>>>>>>>>>>> AT THE END; WILL NEED TO ADD
    # stopifnot(dir.exists(file.path(dsPipOutDir, script10b_name)))
    tad_pvalFCCfile <- file.path(dsPipOutDir, script10b_name, "emp_pval_prodSignedRatio.Rdata")
    if(file.exists(tad_pvalFCCfile)){
      stopifnot(file.exists(tad_pvalFCCfile))
      tad_pvalFCC <- eval(parse(text = load(tad_pvalFCCfile)))
      adj_tad_pvalFCC <- sort(p.adjust(tad_pvalFCC, method="BH"))
      tad_pvalFCC <- sort(tad_pvalFCC)
      
    } else{
      adj_tad_pvalFCC <- NA
      tad_pvalFCC <- NA
      
    }
    
    
    # RETRIEVE TADcorr pval - script10v2
    #>>>>>>>>>>>> AT THE END; WILL NEED TO ADD
    # stopifnot(dir.exists(file.path(dsPipOutDir, script10v2_name)))
    tad_pvalCorrV2File <- file.path(dsPipOutDir, script10v2_name, "emp_pval_meanCorr.Rdata")
    if(file.exists(tad_pvalCorrV2File)) {
      stopifnot(file.exists(tad_pvalCorrV2File))
      tad_pvalCorrV2 <- eval(parse(text = load(tad_pvalCorrV2File)))
      adj_tad_pvalCorrV2 <- sort(p.adjust(tad_pvalCorrV2, method="BH"))
      tad_pvalCorrV2 <- sort(tad_pvalCorrV2)
    } else {
      adj_tad_pvalCorrV2 <- NA
      tad_pvalCorrV2 <- NA
    }
    
    # RETRIEVE TADcorr values - script4
    stopifnot(dir.exists(file.path(dsPipOutDir, script4_name)))
    tad_valuesCorrFile <- file.path(dsPipOutDir, script4_name, "all_meanCorr_TAD.Rdata")
    stopifnot(file.exists(tad_valuesCorrFile))
    tad_valuesCorr <- eval(parse(text = load(tad_valuesCorrFile)))
    tad_valuesCorr <- sort(tad_valuesCorr)
    
    # RETRIEVE TADfcc values - script8c
    stopifnot(dir.exists(file.path(dsPipOutDir, script8c_name)))
    tad_valuesFCCfile <- file.path(dsPipOutDir, script8c_name, "all_obs_prodSignedRatio.Rdata")
    stopifnot(file.exists(tad_valuesFCCfile))
    tad_valuesFCC <- eval(parse(text = load(tad_valuesFCCfile)))
    tad_valuesFCC <- sort(tad_valuesFCC)
    
    
                        # # RETRIEVE SIGNIF. TADs 
                        stopifnot(dir.exists(file.path(dsPipOutDir, script11_name)))
                        tad_pvalFile <- file.path(dsPipOutDir, script11_name, "emp_pval_combined.Rdata")
                        stopifnot(file.exists(tad_pvalFile))
                        tad_pvals <- eval(parse(text = load(tad_pvalFile)))
                        adj_tad_pvalComb <- sort(p.adjust(tad_pvals, method="BH"))
                        tad_pvalComb <- sort(tad_pvals)
    
    
                        all_tads <- Reduce(union, list(names(tad_pvalFC), names(tad_pvalCorr), 
                                                           names(adj_tad_pvalFC), names(adj_tad_pvalCorr),
                                                           # #>>>>>>>>>>>> AT THE END; WILL NEED TO ADD
                                                           names(tad_pvalFCC), names(adj_tad_pvalFCC),
                                                           names(tad_valuesFCC),
                                                           names(tad_pvalComb),
                                                           names(adj_tad_pvalComb),
                                                           names(tad_valuesFC), names(tad_valuesCorr),
                                                           
                                                           #>>>>>>>>>>>> AT THE END; WILL NEED TO ADD
                                                           names(tad_pvalCorrV2), names(tad_pvalWilcox),
                                                           names(adj_tad_pvalCorrV2), names(adj_tad_pvalWilcox)
                        ))                        
                        
    all_tads <- Reduce(intersect, list(names(tad_pvalFC), names(tad_pvalCorr),
                                       names(adj_tad_pvalFC), names(adj_tad_pvalCorr),
                                                                     # #>>>>>>>>>>>> AT THE END; WILL NEED TO ADD
                                                                     names(tad_pvalFCC), names(adj_tad_pvalFCC),
                                       names(tad_valuesFCC),
                                                                       names(tad_pvalComb),
                                                                       names(adj_tad_pvalComb),
                                       names(tad_valuesFC), names(tad_valuesCorr),

                                       #>>>>>>>>>>>> AT THE END; WILL NEED TO ADD
                                        names(tad_pvalCorrV2), names(tad_pvalWilcox),
                                        names(adj_tad_pvalCorrV2), names(adj_tad_pvalWilcox)
    ))

    stopifnot(all_tads %in% names(tad_nGenes))

    stopifnot(length(all_tads) == length(tad_pvalFC))
    stopifnot(length(all_tads) == length(tad_pvalCorr))
    stopifnot(length(all_tads) == length(tad_pvalComb))
    stopifnot(length(all_tads) == length(adj_tad_pvalFC))
    stopifnot(length(all_tads) == length(adj_tad_pvalCorr))
    stopifnot(length(all_tads) == length(adj_tad_pvalComb))
    stopifnot(length(all_tads) == length(tad_valuesFC))
    stopifnot(length(all_tads) == length(tad_valuesCorr))

    #>>>>>>>>>>>> AT THE END; WILL NEED TO ADD
    stopifnot(length(all_tads) == length(tad_valuesFCC))
    stopifnot(length(all_tads) == length(adj_tad_pvalFCC))

    #>>>>>>>>>>>> AT THE END; WILL NEED TO ADD
    stopifnot(length(all_tads) == length(tad_pvalCorrV2))
    stopifnot(length(all_tads) == length(tad_pvalWilcox))
    
    
    outDT <- data.frame(
      hicds=hicds,
      exprds=exprds,
      tad=all_tads,
      tad_id=paste(hicds, exprds, all_tads, sep="_"),
      
      nGenes = tad_nGenes[all_tads],
      
      adj_pvalFC = as.numeric(adj_tad_pvalFC[all_tads]),
      adj_pvalWilcox = as.numeric(adj_tad_pvalWilcox[all_tads]),
      adj_pvalCorr = as.numeric(adj_tad_pvalCorr[all_tads]),
      adj_pvalCorrV2 = as.numeric(adj_tad_pvalCorrV2[all_tads]),
      adj_pvalFCC = as.numeric(adj_tad_pvalFCC[all_tads]),
      
                                                      adj_pvalComb = as.numeric(adj_tad_pvalComb[all_tads]),
      pvalFC = as.numeric(tad_pvalFC[all_tads]),
      pvalWilcox = as.numeric(tad_pvalWilcox[all_tads]),
      pvalCorr = as.numeric(tad_pvalCorr[all_tads]),
      pvalCorrV2 = as.numeric(tad_pvalCorrV2[all_tads]),
      pvalFCC = as.numeric(tad_pvalFCC[all_tads]),
      
                                                          pvalComb = as.numeric(tad_pvalComb[all_tads]),
      valuesFC = as.numeric(tad_valuesFC[all_tads]),
      valuesCorr = as.numeric(tad_valuesCorr[all_tads]),
      valuesFCC = as.numeric(tad_valuesFCC[all_tads]),
      
      stringsAsFactors = FALSE
    )
    #>>>>>>>>>>>> AT THE END; WILL NEED TO ADD
    stopifnot(!is.na(outDT))
    outDT
  }
  outFile <- file.path(outFolder, "allPvals_allDS_DT.Rdata")
  save(allPvals_allDS_DT, file = outFile)
  cat(paste0("... written: ", outFile, "\n"))
  
} else { # if not to build
  outFile <- file.path(outFolder, "allPvals_allDS_DT.Rdata")
  allPvals_allDS_DT <- eval(parse(text = load(outFile)))
}




# stop("--ok\n")




head(allPvals_allDS_DT)
# load("TAD_LOGFC_INTRACORR_DIST/allPvals_allDS_DT.Rdata")
# hicds             exprds adj_pvalFC adj_pvalCorr adj_pvalComb
# 1 ENCSR079VIJ_G401_40kb TCGAkich_norm_kich 0.02913459  0.002099790 2.915007e-05
# 2 ENCSR079VIJ_G401_40kb TCGAkich_norm_kich 0.02913459  0.002950338 3.456350e-05



load("CMP_EMPPVALS/allPvals_allDS_DT.Rdata")
head(allPvals_allDS_DT)


not_to_plot <- c("hicds", "exprds", "tad", "tad_id")

vars_to_plot <- colnames(allPvals_allDS_DT)[!colnames(allPvals_allDS_DT) %in% not_to_plot]

myTit <- "all TCGA datasets - all TADs"
mySub <- paste0("(# TADs = ", length(allPvals_allDS_DT$tad_id) , ")")

# 
# for(curr_var in vars_to_plot) {
#   
#   ############################################## DENSITY
#   outFile <- file.path(outFolder, paste0(curr_var, "_hist.", plotType))
#   do.call(plotType, list(outFile,  height=myHeightDensity, width=myWidthDensity))
#   
#   hist(allPvals_allDS_DT[, curr_var],
#        main = curr_var)
#   mtext(text = paste0(myTit, "; ", mySub), side=3)
#   foo <- dev.off()
#   cat(paste0("... written: ", outFile, "\n"))
#   
#   
# }
# 
# for(curr_var in vars_to_plot) {
#   
#   ############################################## DENSITY
#   outFile <- file.path(outFolder, paste0(curr_var, "_hist_log10.", plotType))
#   do.call(plotType, list(outFile,  height=myHeightDensity, width=myWidthDensity))
#   
#   hist(log10(allPvals_allDS_DT[, curr_var]),
#        main = curr_var)
#   mtext(text = paste0(myTit, "; ", mySub), side=3)
#   foo <- dev.off()
#   cat(paste0("... written: ", outFile, "\n"))
#   
#   
# }
# 
# 
# 
# 
# 
# 
# all_pairs <- combn(vars_to_plot, m=2)
# 
# mySub <- paste0("(tot. DS: ", nrow(unique(allPvals_allDS_DT[,c("hicds", "exprds")])) , 
#                 "; tot. TADs: ",  nrow(allPvals_allDS_DT), ")")
# 
# 
# foo <- foreach(i_col = seq_len(ncol(all_pairs))) %dopar% {
#   
#   var1 <- all_pairs[1,i_col]
#   var2 <- all_pairs[2,i_col]
#   
#   myx <- allPvals_allDS_DT[, var1]
#   myy <- allPvals_allDS_DT[, var2]
#   
#   myTit <- paste0(var2, " vs. ", var1)
#   
#   ############################################## DENSPLOT 
#   
#   outFile <- file.path(outFolder, paste0(var2, "_", var1, "_", "densplot.", plotType))
#   do.call(plotType, list(outFile,  height=myHeight, width=myWidth))
#   densplot(x = myx, y = myy, 
#            main = myTit,
#            xlab = var1,
#            ylab = var2,
#            cex.axis = plotCex, cex.lab=plotCex, cex.main=plotCex)
#   mtext(text=mySub, side=3)
#   foo <- dev.off()
#   cat(paste0("... written: ", outFile, "\n"))
#   
#   ############################################## DENSPLOT LOG10
#   outFile <- file.path(outFolder, paste0(var2, "_", var1, "_", "densplot_log10.", plotType))
#   do.call(plotType, list(outFile,  height=myHeight, width=myWidth))
#   densplot(x = log10(myx), 
#            y = log10(myy), 
#            main = myTit,
#            xlab = paste0(var1, " [log10]"),
#            ylab = paste0(var2, " [log10]"),
#            cex.axis = plotCex, cex.lab=plotCex, cex.main=plotCex)
#   mtext(text=mySub, side=3)
#   foo <- dev.off()
#   cat(paste0("... written: ", outFile, "\n"))
#   
#   ############################################## DENSPLOT 
#   ############################################## DENSPLOT LOG10
# }

barplotCol <- "skyblue"

all_vars <- colnames(allPvals_allDS_DT)[grep("^pval.+", colnames(allPvals_allDS_DT))]

xaxiscols <- c("red")
nPermut <- 10000
minEmpPval <- 1/(nPermut+1)



var <- all_vars[1]
for(var in all_vars){
  # full_var <- paste0("pval", var)
  full_var <- var
  new_var <- paste0("nbrMinPval_", full_var)
  myylab <- paste0("# min ", full_var, "(# perm.=", nPermut, ")")
  myTit <- paste0(var, ": # of TADs with min. pval")
  mySub <- paste0("(tot.: ", sum(allPvals_allDS_DT[, full_var] == minEmpPval), "/", nrow(allPvals_allDS_DT), ")")
  aggDT <- aggregate(as.formula(paste0(full_var, " ~ hicds + exprds")), data = allPvals_allDS_DT, 
                     function(x) sum(x == minEmpPval, na.rm=T))
  colnames(aggDT)[ colnames(aggDT) == full_var] <- new_var
  
  aggDT$dataset <- paste0(aggDT$hicds, "\n", aggDT$exprds)
  aggDT <- aggDT[order(aggDT[,new_var]),]
  
  
  source("colors_utils.R")
  col_labs <- aggDT$dataset
  col_labs <- gsub(".+\n(.+)", "\\1", col_labs)
  col_labs_types <- sapply(col_labs, function(x) as.character(cancer_subAnnot[x]))
  col_labs_types_colors <- sapply(col_labs_types, function(x) as.character(cancer_subColors[x]))
  
  # outFile <- file.path(outFolder, paste0(full_var, "_nbrTAD_minPval_barplot.", plotType))
  outFile <- file.path(outFolder, paste0(full_var, "_nbrTAD_minPval_barplot.", "svg"))
  
  # do.call(plotType, list(outFile,  height=myHeightBarplot, width=myWidthBarplot*1.2))
  do.call("svg", list(outFile,  height=8, width=13))
  
  
  par(oma=c(5,2,2,2))
  
  x <- barplot(aggDT[, new_var],
               # names = aggDT$dataset, 
               col = barplotCol,
               las = 2,
               main = myTit,
               ylab = myylab,
               cex.main = plotCex,
               cex.lab = plotCex,
               cex.axis = plotCex)
  
  mtext(text = mySub, side=3)
  # axis(1, at=x, las=2, labels = aggDT$dataset,cex.axis=0.6,  
  #      col.lab = col_labs_types_colors
  #      ) # for the text label: change cex.axis
  
  allcols <- unique(col_labs_types_colors)
  
  for(i in allcols){
    idx <- which(col_labs_types_colors == i)
    axis(1, at=x[idx,], las=2, labels = aggDT$dataset[idx],cex.axis=0.6,  
         col.axis = i
    )
  }
  
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  
}




stop("--ok\n")









































#============================================================================================
#============================================================================================ plot density for all columns
#============================================================================================

myTit <- "all TCGA datasets - all TADs"
mySub <- paste0("(# TADs = ", length(allPvals_allDS_DT$tad_id) , ")")


############################################## # min Pval

barplotCol <- "skyblue"

all_vars <- c("FC", "Corr")

var <- all_vars[1]
for(var in all_vars){
  full_var <- paste0("pval", var)
  new_var <- paste0("nbrMinPval_", full_var)
  myylab <- paste0("# min ", full_var, "(# perm.=", nPermut, ")")
  myTit <- paste0(var, ": # of TADs with min. pval")
  mySub <- paste0("(tot.: ", sum(allPvals_allDS_DT[, full_var] == minEmpPval), "/", nrow(allPvals_allDS_DT), ")")
  aggDT <- aggregate(as.formula(paste0(full_var, " ~ hicds + exprds")), data = allPvals_allDS_DT, 
                     function(x) sum(x == minEmpPval))
  colnames(aggDT)[ colnames(aggDT) == full_var] <- new_var
  
  aggDT$dataset <- paste0(aggDT$hicds, "\n", aggDT$exprds)
  aggDT <- aggDT[order(aggDT[,new_var]),]
  
  outFile <- file.path(outFolder, paste0(full_var, "_nbrTAD_minPval_barplot.", plotType))
  do.call(plotType, list(outFile,  height=myHeightBarplot, width=myWidthBarplot))
  
  x <- barplot(aggDT[, new_var],
               # names = aggDT$dataset, 
               col = barplotCol,
               las = 2,
               main = myTit,
               ylab = myylab,
               cex.main = plotCex,
               cex.lab = plotCex,
               cex.axis = plotCex)
  
  mtext(text = mySub, side=3)
  axis(1, at=x, las=2, labels = aggDT$dataset,cex.axis=0.7 ) # for the text label: change cex.axis
  
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  
}





all_vars <- colnames(allPvals_allDS_DT)
all_vars <- all_vars[!all_vars %in% c("hicds", "exprds", "tad", "tad_id")]

curr_var <- all_vars[1]
for(curr_var in all_vars) {
  
  ############################################## DENSITY
  outFile <- file.path(outFolder, paste0(curr_var, "_density.", plotType))
  do.call(plotType, list(outFile,  height=myHeightDensity, width=myWidthDensity))
  
  plot(density(allPvals_allDS_DT[, curr_var]),
       main = curr_var)
  mtext(text = paste0(myTit, "; ", mySub), side=3)
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  
}

comb_vars <- combn(all_vars, 2)

i=1
for(i in 1:ncol(comb_vars)){
  var1 <- comb_vars[1,i]
  var2 <- comb_vars[2,i]
  
  myx <- allPvals_allDS_DT[, var1]
  myy <- allPvals_allDS_DT[, var2]
  
  myTit <- paste0(var2, " vs. ", var1)
  
  ############################################## DENSPLOT
  outFile <- file.path(outFolder, paste0(var2, "_", var1, "_", "densplot.", plotType))
  do.call(plotType, list(outFile,  height=myHeight, width=myWidth))
  densplot(x = myx, y = myy, 
           main = myTit,
           xlab = var1,
           ylab = var2,
           cex.axis = plotCex, cex.lab=plotCex, cex.main=plotCex)
  mtext(text=mySub, side=3)
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  ############################################## IF PVAL -> DENSPLOT LOG10
  if(grepl("pval", var1)) {
    myx_log10 <- -log10(myx)
    var1_lab <- paste0(var1, " [-log10]")
  }
  if(grepl("pval", var2)) {
    myy_log10 <- -log10(myy)
    var2_lab <- paste0(var2, " [-log10]")
  }
  outFile <- file.path(outFolder, paste0(var2, "_", var1, "_", "densplot_log10.", plotType))
  do.call(plotType, list(outFile,  height=myHeight, width=myWidth))
  densplot(x = myx_log10, y = myy_log10, 
           main = myTit,
           xlab = var1_lab,
           ylab = var2_lab,
           cex.axis = plotCex, cex.lab=plotCex, cex.main=plotCex)
  mtext(text=mySub, side=3)
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  ############################################## IF PVAL + COMB -> DENSPLOT, COLOR-CODE VALUE
  if( xor(grepl("Comb", var1), grepl("Comb", var2)) & xor(grepl("FC|Corr", var1), grepl("FC|Corr", var2)) ) {
    
    valueCol <- ifelse(grepl("FC", var1) | grepl("FC", var2), "valuesFC", "valuesCorr")
    tmpCols <- colorRampPalette(c('red','blue'))(nColBreaks)[as.numeric(cut(allPvals_allDS_DT[, valueCol],breaks = nColBreaks))]
    
    outFile <- file.path(outFolder, paste0(var2, "_", var1, "_", "colValues.", plotType))
    do.call(plotType, list(outFile,  height=myHeight, width=myWidth))
    plot(x = myx, y = myy, 
         main = myTit,
         xlab = var1,
         ylab = var2,
         col = tmpCols,
         pch=16, cex  = 0.7,
         cex.axis = plotCex, cex.lab=plotCex, cex.main=plotCex)
    mtext(text=mySub, side=3)
    legend("topright", title=paste0(valueCol),
           legend=c(levels(cut(allPvals_allDS_DT[, valueCol],breaks = nColBreaks))),
           col =colorRampPalette(c('red','blue'))(nColBreaks),pch=20, bty="n")
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
    
    outFile <- file.path(outFolder, paste0(var2, "_", var1, "_", "colValues_log10.", plotType))
    do.call(plotType, list(outFile,  height=myHeight, width=myWidth))
    plot(x = -log10(myx), y = -log10(myy), 
         main = myTit,
         xlab = paste0(var1, " [-log10]"),
         ylab = paste0(var2, " [-log10]"),
         col = tmpCols,
         pch=16, cex  = 0.7,
         cex.axis = plotCex, cex.lab=plotCex, cex.main=plotCex)
    mtext(text=mySub, side=3)
    legend("bottomright", title=paste0(valueCol),
           legend=c(levels(cut(allPvals_allDS_DT[, valueCol],breaks = nColBreaks))),
           col =colorRampPalette(c('red','blue'))(nColBreaks),pch=20, bty="n")
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
    
    
  } 
  
  
}





# ######################################################################################
# ######################################################################################
# ######################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))

