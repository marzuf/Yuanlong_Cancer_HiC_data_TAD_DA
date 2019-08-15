# Rscript top_scoring_ratioDown.R
# source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

# 3 => done in Rscript coexpr_DE_queryTAD.R
# chekcer que withinCoexpr mm chose que meanTADcorr !!! [manque pval comb!!!]

script_name <- "top_scoring_ratioDown.R"
cat("> Start ", script_name, "\n")
startTime <- Sys.time()

require(foreach)
require(doMC)

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

registerDoMC(40)

plotType <- "png"
myHeight <- ifelse(plotType=="png", 400, 7)
myWidth <- myHeight
plotCex <- 1.4

myCexAxis <- 1.2
myCexLab <- 1.2


outFolder <- "TOP_SCORING_RATIODOWN"
dir.create(outFolder)

dataFolder <- "COEXPR_BETWEEN_WITHIN_ALL"

pipOutFolder <- file.path("PIPELINE", "OUTPUT_FOLDER")

all_pvalComb_files <- list.files(pipOutFolder, recursive = TRUE, pattern="emp_pval_combined.Rdata", full.names = FALSE)
stopifnot(length(all_pvalComb_files) > 0)

all_fc_files <- list.files(pipOutFolder, recursive = TRUE, pattern="all_meanLogFC_TAD.Rdata", full.names = FALSE)
stopifnot(length(all_fc_files) > 0)

all_meanCorr_files <- list.files(pipOutFolder, recursive = TRUE, pattern="all_meanCorr_TAD.Rdata", full.names = FALSE)
stopifnot(length(all_meanCorr_files) > 0)

dataFile <- file.path(dataFolder, "allData_within_between_coexpr.Rdata")
stopifnot(file.exists(dataFile))
allData_within_between_coexpr <- eval(parse(text = load(dataFile)))

all_ratioDown_files <- list.files(pipOutFolder, recursive = TRUE, pattern="all_obs_ratioDown.Rdata", full.names = FALSE)
stopifnot(length(all_ratioDown_files) > 0)


### BUILD THE LOGFC TABLE
fc_file = all_fc_files[1]
fc_DT <- foreach(fc_file = all_fc_files, .combine = 'rbind') %dopar% {
  curr_file <- file.path(pipOutFolder, fc_file)
  stopifnot(file.exists(curr_file))
  tad_fc <- eval(parse(text = load(curr_file)))
  dataset <- dirname(dirname(fc_file))
  data.frame(
    dataset = dataset,
    region = names(tad_fc),
    meanFC = as.numeric(tad_fc),
    stringsAsFactors = FALSE
  )
}

### BUILD THE MEANCORR TABLE
meanCorr_file = all_meanCorr_files[1]
meanCorr_DT <- foreach(meanCorr_file = all_meanCorr_files, .combine = 'rbind') %dopar% {
  curr_file <- file.path(pipOutFolder, meanCorr_file)
  stopifnot(file.exists(curr_file))
  tad_meanCorr <- eval(parse(text = load(curr_file)))
  dataset <- dirname(dirname(meanCorr_file))
  data.frame(
    dataset = dataset,
    region = names(tad_meanCorr),
    meanCorr = as.numeric(tad_meanCorr),
    stringsAsFactors = FALSE
  )
}


### BUILD THE ratio down TABLE
rd_file = all_ratioDown_files[1]
rD_DT <- foreach(rd_file = all_ratioDown_files, .combine = 'rbind') %dopar% {
  curr_file <- file.path(pipOutFolder,rd_file)
  stopifnot(file.exists(curr_file))
  tad_rd <- eval(parse(text = load(curr_file)))
  dataset <- dirname(dirname(rd_file))
  data.frame(
    dataset = dataset,
    region = names(tad_rd),
    ratioDown = as.numeric(tad_rd),
    stringsAsFactors = FALSE
  )
}

### BUILD THE PVAL COMBINED TABLE
pvalcomb_file = all_pvalComb_files[1]
pvalComb_DT <- foreach(pvalcomb_file = all_pvalComb_files, .combine = 'rbind') %dopar% {
  curr_file <- file.path(pipOutFolder,pvalcomb_file)
  stopifnot(file.exists(curr_file))
  tad_pvalComb <- eval(parse(text = load(curr_file)))
  dataset <- dirname(dirname(pvalcomb_file))
  adj_pvalComb <- p.adjust(tad_pvalComb, method="BH")
  stopifnot(setequal(names(tad_pvalComb), names(adj_pvalComb)))
  data.frame(
    dataset = dataset,
    region = names(tad_pvalComb),
    pvalComb = as.numeric(tad_pvalComb),
    adj_pvalComb = as.numeric(adj_pvalComb[names(tad_pvalComb)]),
    stringsAsFactors = FALSE
  )
}



all_DT <- merge(fc_DT, meanCorr_DT, by =c("dataset", "region"), all = TRUE)
stopifnot(nrow(fc_DT) == nrow(meanCorr_DT))
stopifnot(nrow(fc_DT) == nrow(all_DT))
stopifnot(!is.na(all_DT))

all_DT <- merge(all_DT, pvalComb_DT, by =c("dataset", "region"), all = TRUE)
stopifnot(nrow(all_DT) == nrow(pvalComb_DT))
stopifnot(!is.na(all_DT))

all_DT <- merge(all_DT, rD_DT, by =c("dataset", "region"), all = TRUE)
stopifnot(nrow(all_DT) == nrow(rD_DT))
stopifnot(!is.na(all_DT))


# sort the TADs by decreasing withinCoexpr
# plot level of coexpr within and between on the same plot

tad_coexpr_DT <- data.frame(
  dataset = as.character(unlist(lapply(1:length(allData_within_between_coexpr), function(i) {
    ds_name <- names(allData_within_between_coexpr)[i]
    ds_name <- gsub("^CREATE_COEXPR_SORTNODUP/", "", ds_name)
    ds_name <- gsub("/pearson/coexprDT.Rdata$", "", ds_name)
    rep(ds_name, length(allData_within_between_coexpr[[i]]))
  }))),
  region = as.character(unlist(lapply(1:length(allData_within_between_coexpr), function(i) {
    names(allData_within_between_coexpr[[i]])
  }))),
  
  withinCoexpr = as.numeric(unlist(lapply(allData_within_between_coexpr, 
                                          function(sublist) lapply(sublist, function(x) x[["withinCoexpr"]])))),
  betweenAllCoexpr = as.numeric(unlist(lapply(allData_within_between_coexpr, 
                                              function(sublist) lapply(sublist, function(x) x[["betweenAllCoexpr"]])))),
  betweenKbCoexpr = as.numeric(unlist(lapply(allData_within_between_coexpr, 
                                             function(sublist) lapply(sublist, function(x) x[["betweenKbCoexpr"]])))),
  betweenNbrCoexpr = as.numeric(unlist(lapply(allData_within_between_coexpr, 
                                              function(sublist) lapply(sublist, function(x) x[["betweenNbrCoexpr"]])))),
  withinCoexpr_cond1 = as.numeric(unlist(lapply(allData_within_between_coexpr, 
                                                function(sublist) lapply(sublist, function(x) x[["withinCoexpr_cond1"]])))),
  betweenAllCoexpr_cond1 = as.numeric(unlist(lapply(allData_within_between_coexpr, 
                                                    function(sublist) lapply(sublist, function(x) x[["betweenAllCoexpr_cond1"]])))),
  betweenKbCoexpr_cond1 = as.numeric(unlist(lapply(allData_within_between_coexpr, 
                                                   function(sublist) lapply(sublist, function(x) x[["betweenKbCoexpr_cond1"]])))),
  betweenNbrCoexpr_cond1 = as.numeric(unlist(lapply(allData_within_between_coexpr, 
                                                    function(sublist) lapply(sublist, function(x) x[["betweenNbrCoexpr_cond1"]])))),
  withinCoexpr_cond2 = as.numeric(unlist(lapply(allData_within_between_coexpr, 
                                                function(sublist) lapply(sublist, function(x) x[["withinCoexpr_cond2"]])))),
  betweenAllCoexpr_cond2 = as.numeric(unlist(lapply(allData_within_between_coexpr,
                                                    function(sublist) lapply(sublist, function(x) x[["betweenAllCoexpr_cond2"]])))),
  betweenKbCoexpr_cond2 = as.numeric(unlist(lapply(allData_within_between_coexpr, 
                                                   function(sublist) lapply(sublist, function(x) x[["betweenKbCoexpr_cond2"]])))),
  betweenNbrCoexpr_cond2 = as.numeric(unlist(lapply(allData_within_between_coexpr, 
                                                    function(sublist) lapply(sublist, function(x) x[["betweenNbrCoexpr_cond2"]])))),
  
  stringsAsFactors = FALSE
)

stopifnot(nrow(tad_coexpr_DT) == nrow(all_DT))

all_DT <- merge(all_DT, tad_coexpr_DT, by=c("dataset", "region"), all = TRUE)

stopifnot(nrow(tad_coexpr_DT) == nrow(all_DT))

#stopifnot(!is.na(tad_coexpr_DT))

outFile <- file.path(outFolder, "all_DT.Rdata")
save(all_DT, file = outFile)
cat("... written: ", outFile, "\n")


xvar <- "ratioDown"

all_y <- c("withinCoexpr", "meanFC", "meanCorr", "pvalComb", "adj_pvalComb")

for(yvar in all_y) {
  
  myx <- tad_coexpr_fc_DT[, paste0(xvar)]
  stopifnot(length(myx) > 0)
  
  myy <-  tad_coexpr_fc_DT[, paste0(yvar)]
  stopifnot(length(myy) > 0)
  
  outFile <- file.path(outFolder, paste0(yvar, "_vs_", xvar, "_densplot.", plotType ))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  densplot(
    y = myy,
    x = myx,
    cex.lab=myCexLab,
    cex.axis=myCexAxis,
    ylab=paste0(yvar),
    xlab = paste0(xvar),
    cex = 0.5,
    main = paste0(yvar, " vs. ", xvar)
  )
  addCorr(x = myx, myy = myy, bty="n")
  
  foo <- dev.off()
  cat("... written: ", outFile, "\n")
  
  
  
}



stop("--ok\n")



# ######################################################################################
# ######################################################################################
# ######################################################################################
cat(paste0("*** DONE: ", script_name, "\n"))
cat(paste0(startTime, "\n", Sys.time(), "\n"))





