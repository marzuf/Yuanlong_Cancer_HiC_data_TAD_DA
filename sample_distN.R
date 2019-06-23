# Rscript sample_distN.R

require(foreach)
source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

myHeight <- 400
myWidth <- 400

outFolder <- "SAMPLE_DISTN"
dir.create(outFolder)

sampleNbrFile <- "CREATE_SAMPLE_AROUND_TADS_WITHDISTN/all_ds_sample_around_TADs.Rdata"
sampleKbFile <- "CREATE_SAMPLE_AROUNDKB_TADS_WITHDISTN/2000000/all_ds_sample_aroundKb_TADs.Rdata"
stopifnot(file.exists(sampleNbrFile))
stopifnot(file.exists(sampleKbFile))

sampleNbr <- eval(parse(text = load(sampleNbrFile)))
sampleKb <- eval(parse(text = load(sampleKbFile)))

sampleNbr_DT <- foreach(i = names(sampleNbr), .combine='rbind') %dopar% {
  hicds = dirname(i)
  exprds = basename(i)
  nGenesNbr = unlist(lapply(sampleNbr[[i]], function(x) x[["nGenes"]]))
  maxDistNbr = unlist(lapply(sampleNbr[[i]], function(x) x[["maxDist"]]))
  
  data.frame(
    hicds=hicds,
    exprds=exprds,
    nGenesNbr=nGenesNbr,
    maxDistNbr=maxDistNbr,
  stringsAsFactors=FALSE
  )
  
}
outFile <- file.path(outFolder, "maxDist_vs_nGenes_betweenNbr.png")
png(outFile, height=myHeight, width=myWidth)
densplot(x = sampleNbr_DT$nGenesNbr, 
         y = sampleNbr_DT$maxDistNbr/1000,
         xlab = "nGenes",
         ylab = "maxDist (kb)",
         main = "maxDist vs. nGenes")
mtext(side=3, text = "betweenNbr", font=3)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


sampleKb_DT <- foreach(i = names(sampleNbr), .combine='rbind') %dopar% {
  hicds = dirname(i)
  exprds = basename(i)
  nGenesKb = unlist(lapply(sampleKb[[i]], function(x) x[["nGenes"]]))
  maxDistKb = unlist(lapply(sampleKb[[i]], function(x) x[["maxDist"]]))
  
  data.frame(
    hicds=hicds,
    exprds=exprds,
    nGenesKb=nGenesKb,
    maxDistKb=maxDistKb,
    stringsAsFactors=FALSE
  )
  
}

outFile <- file.path(outFolder, "maxDist_vs_nGenes_betweenKb.png")
png(outFile, height=myHeight, width=myWidth)
densplot(x = sampleKb_DT$nGenesKb,
         xlab = "nGenes",
         ylab = "maxDist (kb)",
     y = sampleKb_DT$maxDistKb/1000,
     main = "maxDist vs. nGenes")
mtext(side=3, text = "betweenKb", font=3)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))
