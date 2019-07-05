require(foreach)
require(doMC)

startTime <- Sys.time()

# Rscript cmp_empPval_rank.R

# source("/media/electron/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2/TAD_DE_utils.R")
source("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2/TAD_DE_utils.R")
source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

# > load("EMPPVAL_MEANTADCORR_RANK/Panc1_rep12_40kb/TCGApaad_wt_mutKRAS/emp_pval_meanLogFC_rank.Rdata")
all_fc_files <- list.files(".", pattern="emp_pval_meanLogFC_rank.Rdata", full.names = TRUE, recursive = TRUE)

corrFolder <- "EMPPVAL_MEANTADCORR_RANK"

stopifnot(length(all_fc_files) > 0)

outFold <- "CMP_EMPPVAL_RANK"
dir.create(outFold)


all_dt <- foreach(fc_file = all_fc_files, .combine='rbind') %dopar% {
  
  stopifnot(file.exists(fc_file))
  
  hicds <- basename(dirname(dirname(fc_file)))
  exprds <- basename(dirname(fc_file))
  
  fc_pval <- eval(parse(text = load(fc_file)))
  
  corr_file <- file.path(corrFolder, hicds, exprds, "emp_pval_meanCorr_rank.Rdata")
  stopifnot(file.exists(corr_file))
  
  corr_pval <- eval(parse(text = load(corr_file)))
  
  stopifnot(setequal(names(corr_pval), names(fc_pval)))
  
  all_tads <- names(corr_pval)
  
  data.frame(
    hicds = hicds,
    exprds = exprds,
    region = all_tads,
    empPval_logFC = fc_pval[all_tads],
    empPval_meanCorr = corr_pval[all_tads],
    stringsAsFactors = FALSE
  )
  
  
}

all_dt$adj_empPval_logFC <- p.adjust(all_dt$empPval_logFC, method="BH")
all_dt$adj_empPval_meanCorr <- p.adjust(all_dt$empPval_meanCorr, method="BH")


emp_pval_combined <- unlist(sapply(1:nrow(all_dt), function(x)
  stouffer(c(all_dt[x,"empPval_meanCorr"], all_dt[x, "empPval_logFC"]), two.tails = TRUE)))

all_dt$empPval_combined <- emp_pval_combined
all_dt$adj_empPval_combined <- p.adjust(all_dt$empPval_combined, method="BH")


range(all_dt$empPval_logFC)
range(all_dt$empPval_meanCorr)
range(all_dt$empPval_combined)

range(all_dt$adj_empPval_combined)
range(all_dt$adj_empPval_logFC)
range(all_dt$adj_empPval_meanCorr)

nDS <- length(unique(paste0(all_dt$hicds, "_", all_dt$exprds)))

xvar <- "empPval_meanCorr"
yvar <- "empPval_logFC"

stopifnot(xvar %in% colnames(all_dt))
stopifnot(yvar %in% colnames(all_dt))

outFile <- file.path(outFold, paste0(yvar, "_vs_", xvar, ".png"))
do.call("png", list(outFile, height=400, width=400))
densplot(
  x = -log10(all_dt[,xvar]),
  y = -log10(all_dt[,yvar]),
  xlab = paste0(xvar, " (-log10)"),
  ylab = paste0(yvar, " (-log10)"),
  cex.axis=1.2,
  cex.lab=1.2,
  main = paste0("empPval rank FC vs. Corr")
)
mtext(side = 3, text = paste0("nDS = ", nDS, ")"))
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


xvar <- "adj_empPval_meanCorr"
yvar <- "adj_empPval_logFC"

stopifnot(xvar %in% colnames(all_dt))
stopifnot(yvar %in% colnames(all_dt))

outFile <- file.path(outFold, paste0(yvar, "_vs_", xvar, ".png"))
do.call("png", list(outFile, height=400, width=400))
densplot(
  x = -log10(all_dt[,xvar]),
  y = -log10(all_dt[,yvar]),
  xlab = paste0(xvar, " (-log10)"),
  ylab = paste0(yvar, " (-log10)"),
  cex.axis=1.2,
  cex.lab=1.2,
  main = paste0("adj. empPval rank FC vs. Corr")
)
mtext(side = 3, text = paste0("nDS = ", nDS, ")"))
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


xvar <- "empPval_logFC"
yvar <- "adj_empPval_combined"

stopifnot(xvar %in% colnames(all_dt))
stopifnot(yvar %in% colnames(all_dt))

outFile <- file.path(outFold, paste0(yvar, "_vs_", xvar, ".png"))
do.call("png", list(outFile, height=400, width=400))
densplot(
  x = -log10(all_dt[,xvar]),
  y = -log10(all_dt[,yvar]),
  xlab = paste0(xvar, " (-log10)"),
  ylab = paste0(yvar, " (-log10)"),
  cex.axis=1.2,
  cex.lab=1.2,
  main = paste0("adj. empPval rank FC vs. Corr")
)
mtext(side = 3, text = paste0("nDS = ", nDS, ")"))
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


xvar <- "empPval_meanCorr"
yvar <- "adj_empPval_combined"

stopifnot(xvar %in% colnames(all_dt))
stopifnot(yvar %in% colnames(all_dt))

outFile <- file.path(outFold, paste0(yvar, "_vs_", xvar, ".png"))
do.call("png", list(outFile, height=400, width=400))
densplot(
  x = -log10(all_dt[,xvar]),
  y = -log10(all_dt[,yvar]),
  xlab = paste0(xvar, " (-log10)"),
  ylab = paste0(yvar, " (-log10)"),
  cex.axis=1.2,
  cex.lab=1.2,
  main = paste0("adj. empPval rank FC vs. Corr")
)
mtext(side = 3, text = paste0("nDS = ", nDS, ")"))
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


################################****************************************************************************************
################################****************************************************************************************

cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))









