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
plot(sampleNbr_DT$nGenesNbr, sampleNbr_DT$maxDistNbr)

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

plot(sampleKb_DT$nGenesKb, sampleKb_DT$maxDistKb)
