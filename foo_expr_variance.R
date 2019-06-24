library(foreach)
library(doMC)

# Rscript expr_variance.R

SSHFS <- FALSE
setDir <- ifelse(SSHFS, "~/media/electron", "")

registerDoMC(ifelse(SSHFS, 2, 40))

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

varianceDT_file <- "GENE_VARIANCE/LOG2FPKM/all_ds_geneVarDT.Rdata"
stopifnot(file.exists(varianceDT_file))
load(varianceDT_file)

nTop <- 1000

plotType <- "svg"
myHeight <- 7
myWidth <- 10

plotCex <- 1.2

outFold <- "EXPR_VARIANCE"
dir.create(outFold, recursive = TRUE)

ds_pipFolder <- file.path("PIPELINE", "OUTPUT_FOLDER")
ds_settingFolder <- file.path("PIPELINE", "INPUT_FILES")

ctype <- c("luad", "lihc")

allData_allDs <- foreach(dataset = ctype) %do% {
  
  if(dataset == "luad"){
    exprds1 <- "TCGAluad_norm_luad"
    exprds2 <- "TCGAluad_mutKRAS_mutEGFR"
    exprds3 <- "TCGAluad_nonsmoker_smoker"
    exprds4 <- "TCGAluad_wt_mutKRAS"
    allDS <- c(exprds1, exprds2, exprds3, exprds4)
    hicds <- "NCI-H460_40kb"
    
  }else if(dataset == "lihc"){
    exprds1 <- "TCGAlihc_norm_lihc"  
    exprds2 <- "TCGAlihc_wt_mutCTNNB1"
    allDS <- c(exprds1, exprds2)
    hicds <- "GSE105381_HepG2_40kb"
  }

  varName <- "meanMostVar"
  
  allData <- foreach(ds = allDS) %dopar% {
    
    calc_var <- all_ds_geneVarDT[all_ds_geneVarDT$hicds == hicds & all_ds_geneVarDT$exprds == ds, varName]

    setting_file <- file.path(ds_settingFolder, hicds, paste0("run_settings_", ds, ".R"))
    stopifnot(file.exists(setting_file))
    source(setting_file)
    
    samp1_ID <- eval(parse(text=load(file.path(setDir, sample1_file))))
    samp2_ID <- eval(parse(text=load(file.path(setDir, sample2_file))))
    
    ##### FPKM
    
    fpkm_file <- file.path(ds_pipFolder, hicds, ds, "0_prepGeneData", paste0("rna_", "fpkm", "DT.Rdata"))
    stopifnot(file.exists(fpkm_file))
    
    fpkmDT <- eval(parse(text = load(fpkm_file)))
    samp1_DT <- fpkmDT[, samp1_ID]
    samp2_DT <- fpkmDT[, samp2_ID]
    
    stopifnot(nrow(samp1_DT) > 0)
    stopifnot(nrow(samp2_DT) > 0)
    
    samp1_var <- apply(samp1_DT, 1, var, na.rm=T)
    samp2_var <- apply(samp2_DT, 1, var, na.rm=T)
    
    stopifnot( length(intersect(samp1_ID, samp2_ID)) == 0)
    
    list(
      calc_var = calc_var,
      samp1_ID = samp1_ID,
      samp2_ID = samp2_ID,
      cond1 = cond1,
      cond2 = cond2,
      # nSamp1 = length(samp1_ID),
      # nSamp2 = length(samp2_ID),
      # sampIntersect = length(intersect(samp1_ID, samp2_ID)),
      # sampUnion = length(union(samp1_ID, samp2_ID)),
      samp1_var = samp1_var,
      samp2_var = samp2_var
    )
    
  }
  names(allData) <- allDS
  allData
}
names(allData_allDs) <- ctype

outFile <- file.path(outFold, "allData_allDs.Rdata")
save(allData_allDs, file = outFile)
cat(paste0("... written: ", outFile, "\n"))


# load("CMP_VARIANCE_NORM_MUT_V2/allData_allDs.Rdata")
# 
# cancerT="luad"
for(cancerT in ctype) {

  sub_data <- allData_allDs[[cancerT]]
  
  # all_var1 <- lapply(sub_data, function(x) log10(x[["samp1_var"]]))
  # all_var2 <- lapply(sub_data, function(x) log10(x[["samp2_var"]]))
  # 
  all_var1 <- lapply(sub_data, function(x) {
    tmp <- sort(na.omit(x[["samp1_var"]]), decreasing=T)
    tmp <- tmp[1:nTop]
    log10(tmp)
  })
  
  all_var2 <- lapply(sub_data, function(x) {
    tmp <- sort(na.omit(x[["samp2_var"]]), decreasing=T)
    tmp <- tmp[1:nTop]
    log10(tmp)
    })
  
  cond1 <- unlist(lapply(sub_data, function(x) x[["cond1"]]))
  cond2 <- unlist(lapply(sub_data, function(x) x[["cond2"]]))
  all_conds <- c(cond1, cond2)
  i_norm_cond <- which(all_conds == "norm")
  stopifnot(length(i_norm_cond) == 1)

  samp1_ID <-  lapply(sub_data, function(x) x[["samp1_ID"]])
  samp2_ID <-  lapply(sub_data, function(x) x[["samp2_ID"]])
  all_samps <- c(samp1_ID, samp2_ID)
  
  norm_samps <- all_samps[[i_norm_cond]]
  
  nSamp1_ID <- unlist(lapply(sub_data, function(x) length(x[["samp1_ID"]])))
  nSamp2_ID <- unlist(lapply(sub_data, function(x) length(x[["samp2_ID"]])))
  all_nSamp <- c(nSamp1_ID, nSamp2_ID)
  
  nNorm <- all_nSamp[i_norm_cond]
  
  all_nSamp_inNorm <- lapply(all_samps, function(x) {
    sum(x %in% unlist(norm_samps))
  })
  
  mytit <- paste0(toupper(cancerT))
  mysub <- paste0("top ", nTop, " most var genes ")
  myleg <- paste0(
    all_conds," (",  all_nSamp, "; ", all_nSamp_inNorm, "/", nNorm, ")"
  )
  
  myxlab <- "gene var. FPKM (all genes)"
  
  outFile <- file.path(outFold, paste0(cancerT, "_top", nTop, "_mostVar.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  plot_multiDens_argList(c(all_var1, all_var2),
                 # legTxt = c("cond (nSamp; nSampInNorm/nNormSamp)", myleg),
                 legTxt = c( myleg),
                 my_xlab = myxlab,
                 plotTit = mytit,
                 cex.axis = plotCex,
                 cex.lab = plotCex
                )
  mtext(text=mysub, side=3)
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
}




