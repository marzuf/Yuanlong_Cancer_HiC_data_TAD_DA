# Rscript coexpr_cond1_cond2_permut.R

hicds <- "ENCSR401TBQ_Caki2_40kb"
exprds <- "TCGAkich_norm_kich"


outFolder <- "COEXPR_COND1_COND2_PERMUT"
dir.create(outFolder)

permutMeanCorrFolder <- "PERMUTMEANTADCORR_BYCOND"

meanCorrFolder <- "MEANTADCORR_BYCOND"

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

plotType <- "png"
myHeight <- ifelse(plotType == "png", 400, 7)
myWidth <- myHeight

cond1_permDT_file <- file.path(permutMeanCorrFolder, hicds, exprds, "meanCorr_permDT_cond1.Rdata")
stopifnot(file.exists(cond1_permDT_file))

cond1_obsCorr_file <- file.path(meanCorrFolder, hicds, exprds, "all_meanCorr_TAD_cond1.Rdata")
stopifnot(file.exists(cond1_obsCorr_file))
cond1_obs <- eval(parse(text = load(cond1_obsCorr_file)))

cond2_permDT_file <- file.path(permutMeanCorrFolder, hicds, exprds, "meanCorr_permDT_cond2.Rdata")
stopifnot(file.exists(cond2_permDT_file))

cond2_obsCorr_file <- file.path(meanCorrFolder, hicds, exprds, "all_meanCorr_TAD_cond2.Rdata")
stopifnot(file.exists(cond2_obsCorr_file))
cond2_obs <- eval(parse(text = load(cond2_obsCorr_file)))

stopifnot(length(cond1_obs) == length(cond2_obs))
stopifnot(names(cond1_obs) == names(cond2_obs))

myTit <- paste0(hicds, " - ", exprds, " (obs.)")
mySub <- paste0("obs. mean TAD expr. corr. (n=", length(cond1_obs) , ")")

outFile <- file.path(outFolder, paste0(hicds, "_", exprds, "_obs_cond2_vs_cond1.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
densplot(
  x = cond1_obs,
  y = cond2_obs,
  xlab = "cond1",
  ylab = "cond2",
  main = paste0(myTit)
)
mtext(side = 3, text = mySub)
curve(1*x, add=T, lty=2, col = "grey")
addCorr(x=cond1_obs, y=cond2_obs)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

cat("... load permDT cond1\n")
cond1_permDT <- eval(parse(text = load(cond1_permDT_file)))
# head(cond1_permDT)

cat("... load permDT cond2\n")
cond2_permDT <- eval(parse(text = load(cond2_permDT_file)))
# head(cond2_permDT)

stopifnot(rownames(cond1_permDT) == rownames(cond2_permDT))
stopifnot(colnames(cond1_permDT) == colnames(cond2_permDT))


########################################################################################


nSub <- 100

sub_cond1_permDT <- cond1_permDT[,1:nSub]

save(sub_cond1_permDT, file="sub_cond1_permDT.Rdata")
sub_cond2_permDT <- cond2_permDT[,1:nSub]
save(sub_cond2_permDT, file= "sub_cond2_permDT.Rdata")

myTit <- paste0(hicds, " - ", exprds, " (permut.)")
mySub <- paste0("perm. mean TAD expr. corr. (nPerm=", nSub , " x ", nrow(sub_cond1_permDT), ")")

outFile <- file.path(outFolder, paste0(hicds, "_", exprds, "_permut", nSub, "_cond2_vs_cond1.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
densplot(
  x = unlist(c(sub_cond1_permDT)),
  y = unlist(c(sub_cond2_permDT)),
  xlab = "cond1",
  ylab = "cond2",
  main = paste0(myTit)
)
addCorr(x=unlist(c(sub_cond1_permDT)), y=unlist(c(sub_cond2_permDT)))
mtext(side = 3, text = mySub)
curve(1*x, add=T, lty=2, col = "grey")
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

########################################################################################

nSub <- 1000
sub_cond1_permDT <- cond1_permDT[,1:nSub]
sub_cond2_permDT <- cond2_permDT[,1:nSub]

myTit <- paste0(hicds, " - ", exprds, " (permut.)")
mySub <- paste0("perm. mean TAD expr. corr. (nPerm=", nSub , " x ", nrow(sub_cond1_permDT), ")")

outFile <- file.path(outFolder, paste0(hicds, "_", exprds, "_permut", nSub, "_cond2_vs_cond1.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
densplot(
  x = unlist(c(sub_cond1_permDT)),
  y = unlist(c(sub_cond2_permDT)),
  xlab = "cond1",
  ylab = "cond2",
  main = paste0(myTit)
)
addCorr(x=unlist(c(sub_cond1_permDT)), y=unlist(c(sub_cond2_permDT)))
mtext(side = 3, text = mySub)
curve(1*x, add=T, lty=2, col = "grey")
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

########################################################################################



myTit <- paste0(hicds, " - ", exprds, " (permut.)")
mySub <- paste0("perm. mean TAD expr. corr. (all perm=", ncol(sub_cond1_permDT) , " x ", nrow(sub_cond1_permDT), ")")

outFile <- file.path(outFolder, paste0(hicds, "_", exprds, "_permut", "All", "_cond2_vs_cond1.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
densplot(
  x = unlist(c(cond1_permDT)),
  y = unlist(c(cond2_permDT)),
  xlab = "cond1",
  ylab = "cond2",
  main = paste0(myTit)
)
addCorr(x=unlist(c(cond1_permDT)), y=unlist(c(cond2_permDT)))
mtext(side = 3, text = mySub)
curve(1*x, add=T, lty=2, col = "grey")
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))







