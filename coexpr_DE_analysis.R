# Rscript coexpr_DE_analysis.R

script_name <- "coexpr_DE_analysis.R"

startTime <- Sys.time()

cat("> START coexpr_DE_analysis.R \n")

SSHFS <- FALSE

require(tools)
require(foreach)
require(doMC)
registerDoMC(ifelse(SSHFS, 2, 40))

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")
source("subtype_cols.R")

plotType <- "png"
myHeight <- ifelse(plotType=="png", 500, 7)
myWidth <- myHeight
plotCex <- 1.4

myCexAxis <- 1.2
myCexLab <- 1.2

windowSizeBp <- 500*10^3
options(scipen=100)


outFolder <- "COEXPR_DE_ANALYSIS"
dir.create(outFolder, recursive=TRUE)

dataFolder <- "COEXPR_BETWEEN_WITHIN_ALL"

corrMet <- "pearson"
pipOutFolder <- file.path("PIPELINE", "OUTPUT_FOLDER")
stopifnot(dir.exists(pipOutFolder))

all_fc_files <- list.files(pipOutFolder, recursive = TRUE, pattern="all_meanLogFC_TAD.Rdata", full.names = FALSE)
stopifnot(length(all_fc_files) > 0)

all_pvalcomb_files <- list.files(pipOutFolder, recursive = TRUE, pattern="emp_pval_combined.Rdata", full.names = FALSE)
stopifnot(length(all_pvalcomb_files) > 0)

dataFile <- file.path(dataFolder, "allData_within_between_coexpr.Rdata")
stopifnot(file.exists(dataFile))
allData_within_between_coexpr <- eval(parse(text = load(dataFile)))

all_domainScore_files <- list.files(".", recursive = TRUE, pattern="_final_domains_withScore.txt", full.names = FALSE)
stopifnot(length(all_domainScore_files) > 0)

all_ratioDown_files <- list.files(pipOutFolder, recursive = TRUE, pattern="all_obs_ratioDown.Rdata", full.names = FALSE)
stopifnot(length(all_ratioDown_files) > 0)

stopifnot(length(all_ratioDown_files) == length(all_fc_files))
stopifnot(length(all_fc_files) == length(all_pvalcomb_files) )
  
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
tad_coexpr_DT <- tad_coexpr_DT[order(tad_coexpr_DT$withinCoexpr, decreasing = TRUE),]
tad_coexpr_DT$TADrank <- 1:nrow(tad_coexpr_DT)

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

### BUILD THE CPTMT SCORE TABLE
score_file = all_domainScore_files[1]
score_DT <- foreach(score_file = all_domainScore_files, .combine = 'rbind') %dopar% {
  curr_file <- file.path(score_file)
  stopifnot(file.exists(curr_file))
  curr_DT <- read.delim(curr_file, header=F, 
                        col.names = c("chromo", "start", "end", "region", "score"))
  curr_DT$dataset <- dirname(dirname(curr_file))
  curr_DT
}

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


### BUILD THE LOGFC TABLE
pvalcomb_file = all_pvalcomb_files[1]
pvalcomb_DT <- foreach(pvalcomb_file = all_pvalcomb_files, .combine = 'rbind') %dopar% {
  curr_file <- file.path(pipOutFolder, pvalcomb_file)
  stopifnot(file.exists(curr_file))
  tad_pvalcomb <- eval(parse(text = load(curr_file)))
  # adj pval comb
  adj_tad_pvalcomb <- p.adjust(tad_pvalcomb, method="BH")
  dataset <- dirname(dirname(pvalcomb_file))
  data.frame(
    dataset = dataset,
    region = names(adj_tad_pvalcomb),
    adjPvalComb = as.numeric(adj_tad_pvalcomb),
    stringsAsFactors = FALSE
  )
}



############################################################################## MERGE THE TABLES

fc_adjpval_DT <- merge(pvalcomb_DT, fc_DT, by=c("dataset", "region"))
stopifnot(nrow(fc_adjpval_DT) == nrow(fc_DT))
stopifnot(nrow(fc_DT) == nrow(pvalcomb_DT))
stopifnot(!is.na(fc_adjpval_DT))

fc_rd_DT <- merge(rD_DT, fc_adjpval_DT, by=c("dataset", "region"))
stopifnot(nrow(fc_rd_DT) == nrow(rD_DT))
stopifnot(nrow(fc_DT) == nrow(rD_DT))
stopifnot(!is.na(fc_rd_DT))

stopifnot(fc_rd_DT$dataset %in% tad_coexpr_DT$dataset)
stopifnot(tad_coexpr_DT$dataset %in% fc_rd_DT$dataset)

tad_coexpr_fc_DT <- merge(tad_coexpr_DT, fc_rd_DT, by=c("dataset", "region"))
tad_coexpr_fc_DT <- tad_coexpr_fc_DT[order(tad_coexpr_fc_DT$withinCoexpr, decreasing = TRUE),]
tad_coexpr_fc_DT$TADrank <- 1:nrow(tad_coexpr_fc_DT)


tad_coexpr_fc_DT$withinDiffCond1Cond2 <- (tad_coexpr_fc_DT$withinCoexpr_cond1 - tad_coexpr_fc_DT$withinCoexpr_cond2)
tad_coexpr_fc_DT$withinRatioCond1Cond2 <- (tad_coexpr_fc_DT$withinCoexpr_cond1 / tad_coexpr_fc_DT$withinCoexpr_cond2)
tad_coexpr_fc_DT$withinChangeratioCond1Cond2 <- (tad_coexpr_fc_DT$withinCoexpr_cond2 - tad_coexpr_fc_DT$withinCoexpr_cond1)/tad_coexpr_fc_DT$withinCoexpr_cond1

tad_coexpr_fc_DT$withinBetweenAllDiffCond1 <- (tad_coexpr_fc_DT$withinCoexpr_cond1 - tad_coexpr_fc_DT$betweenAllCoexpr_cond1) 
tad_coexpr_fc_DT$withinBetweenAllDiffCond2 <- (tad_coexpr_fc_DT$withinCoexpr_cond2 - tad_coexpr_fc_DT$betweenAllCoexpr_cond2) 


tad_coexpr_fc_DT$withinBetweenNbrDiffCond1 <- (tad_coexpr_fc_DT$withinCoexpr_cond1 - tad_coexpr_fc_DT$betweenNbrCoexpr_cond1) 
tad_coexpr_fc_DT$withinBetweenNbrDiffCond2 <- (tad_coexpr_fc_DT$withinCoexpr_cond2 - tad_coexpr_fc_DT$betweenNbrCoexpr_cond2) 


tad_coexpr_fc_DT$withinBetweenKbDiffCond1 <- (tad_coexpr_fc_DT$withinCoexpr_cond1 - tad_coexpr_fc_DT$betweenKbCoexpr_cond1) 
tad_coexpr_fc_DT$withinBetweenKbDiffCond2 <- (tad_coexpr_fc_DT$withinCoexpr_cond2 - tad_coexpr_fc_DT$betweenKbCoexpr_cond2) 


tad_coexpr_fc_DT$withinBetweenDiffAll <- (tad_coexpr_fc_DT$withinCoexpr - tad_coexpr_fc_DT$betweenAllCoexpr) 
tad_coexpr_fc_DT$withinBetweenRatioAll <- (tad_coexpr_fc_DT$withinCoexpr / tad_coexpr_fc_DT$betweenAllCoexpr) 

tad_coexpr_fc_DT$withinBetweenDiffNbr <- (tad_coexpr_fc_DT$withinCoexpr - tad_coexpr_fc_DT$betweenNbrCoexpr) 
tad_coexpr_fc_DT$withinBetweenRatioNbr <- (tad_coexpr_fc_DT$withinCoexpr / tad_coexpr_fc_DT$betweenNbrCoexpr) 

tad_coexpr_fc_DT$withinBetweenDiffKb <- (tad_coexpr_fc_DT$withinCoexpr - tad_coexpr_fc_DT$betweenKbCoexpr) 
tad_coexpr_fc_DT$withinBetweenRatioKb <- (tad_coexpr_fc_DT$withinCoexpr / tad_coexpr_fc_DT$betweenKbCoexpr) 


# select only that with + coexpr in both condition -> I can take logFC

tad_coexpr_fc_DT$withinBetwNbrLogFC <- log10(tad_coexpr_fc_DT$withinCoexpr/tad_coexpr_fc_DT$betweenNbrCoexpr)

tad_coexpr_fc_DT$withinBetwNbrCond1LogFC <- log10(tad_coexpr_fc_DT$withinCoexpr_cond1/tad_coexpr_fc_DT$betweenNbrCoexpr_cond1)
tad_coexpr_fc_DT$withinBetwNbrCond2LogFC <- log10(tad_coexpr_fc_DT$withinCoexpr_cond2/tad_coexpr_fc_DT$betweenNbrCoexpr_cond2)

tad_coexpr_fc_DT$withinCond2WithinCond1LogFC <- log10(tad_coexpr_fc_DT$withinCoexpr_cond2/tad_coexpr_fc_DT$withinCoexpr_cond1)

tad_coexpr_fc_DT$adjPvalComb_log10 <- -log10(tad_coexpr_fc_DT$adjPvalComb)


# for each, plot i) densplot; ii) plot with color for subtypes
tad_coexpr_fc_DT$cmps <- basename(tad_coexpr_fc_DT$dataset)

colDT <- data.frame(
  cmps = names(all_cmps),
  cmpType = all_cmps,
stringsAsFactors = FALSE
)
colDT$cmpCol <- all_cols[colDT$cmpType]
stopifnot(!is.na(colDT))

stopifnot(tad_coexpr_fc_DT$cmps %in% colDT$cmps)

tad_coexpr_fc_DT <- merge(tad_coexpr_fc_DT, colDT, by = "cmps", all.x = TRUE, all.y = FALSE)

stopifnot(!is.na(tad_coexpr_fc_DT$cmpCol))


outFile <- file.path(outFolder, "tad_coexpr_fc_DT.Rdata")
save(tad_coexpr_fc_DT, file = outFile)
cat(paste0("... written: ", outFile, "\n"))


myplot_densplot <- function(xvar, yvar, addCurve=FALSE, dt = tad_coexpr_fc_DT, outPrefix="") {
  
  stopifnot(xvar %in% colnames(dt))
  myx <- dt[,xvar]
  stopifnot(yvar %in% colnames(dt))
  myy <- dt[,yvar]
  mycols <- dt$cmpCol
  
  outFile <- file.path(outFolder, paste0(outPrefix, "", yvar, "_vs_", xvar, "_densplot.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  densplot(y=myy,
           x=myx,
           pch=16,
           cex=0.7,
           cex.axis = myCexAxis,
           cex.lab = myCexLab,
           xlab=xvar,
           ylab=yvar,
           main = paste0(yvar, " vs. ", xvar)
  )
  addCorr(x=myx,y=myy,legPos="topleft", bty='n')
  if(addCurve) {
    curve(1*x, lty=2, col="grey", add = TRUE)
  }
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
}

myplot_colplot <- function(xvar, yvar, mycols, addCurve = FALSE, dt = tad_coexpr_fc_DT, outPrefix="") {
  
  stopifnot(xvar %in% colnames(dt))
  myx <- dt[,xvar]
  stopifnot(yvar %in% colnames(dt))
  myy <- dt[,yvar]
  
  outFile <- file.path(outFolder, paste0(outPrefix, "", yvar, "_vs_", xvar, "_colplot.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  plot(y=myy,
           x=myx,
           pch=16,
           cex=0.7,
           cex.axis = myCexAxis,
           cex.lab = myCexLab,
           xlab=xvar,
           ylab=yvar,
       col=mycols,
           main = paste0(yvar, " vs. ", xvar)
  )
  addCorr(x=myx,y=myy,legPos="topleft", bty='n')
  if(addCurve) {
    curve(1*x, lty=2, col="grey",add=TRUE)
  }
  addSubtypeLeg(bty="n")
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
}

mycols <- tad_coexpr_fc_DT$cmpCol

##########################
### adj. pval. combined and withinCoexpr diff cond1 cond2
##########################

yvar <- "adjPvalComb_log10"
xvar <- "withinDiffCond1Cond2"

myplot_densplot(xvar,yvar)
myplot_colplot(xvar,yvar,mycols)

##########################
### adj. pval. combined and withinBetweenCoexpr diff (kb)
##########################

yvar <- "adjPvalComb_log10"
xvar <- "withinBetweenDiffKb"

myplot_densplot(xvar,yvar)
myplot_colplot(xvar,yvar,mycols)


##########################
### adj. pval. combined and withinBetweenCoexpr diff (nbr)
##########################

yvar <- "adjPvalComb_log10"
xvar <- "withinBetweenDiffNbr"

myplot_densplot(xvar,yvar)
myplot_colplot(xvar,yvar,mycols)


##########################
### detect "disruption" and DE: those with high FC in expression and high FC in coexpression => FC coexpr vs. FC expr.
##########################

yvar <- "withinBetwNbrLogFC"
xvar <- "meanFC"

myplot_densplot(xvar,yvar)
myplot_colplot(xvar,yvar,mycols)

##########################
### detect "disruption" and DE: those with high/low ratioDown and high FC in coexpression => FC coexpr vs. ratio Down
##########################

yvar <- "withinBetwNbrLogFC"
xvar <- "ratioDown"

myplot_densplot(xvar,yvar)
myplot_colplot(xvar,yvar,mycols)



##########################
### detect "disruption": those with change in coexpr => coexpr cond2 vs coexpr cond1
##########################

yvar <- "withinCoexpr_cond2"
xvar <- "withinCoexpr_cond1"

myplot_densplot(xvar,yvar, addCurve = TRUE)
myplot_colplot(xvar,yvar,mycols, addCurve = TRUE)

mycols_sub <- tad_coexpr_fc_DT$cmpCol[tad_coexpr_fc_DT$cmpType == "subtypes"]
mycols_tumor <- tad_coexpr_fc_DT$cmpCol[tad_coexpr_fc_DT$cmpType == "norm_vs_tumor"]
mycols_mut <- tad_coexpr_fc_DT$cmpCol[tad_coexpr_fc_DT$cmpType == "wt_vs_mut"]

myplot_densplot(xvar,yvar, addCurve = TRUE, 
                dt = tad_coexpr_fc_DT[tad_coexpr_fc_DT$cmpType == "subtypes",] , outPrefix = "subtypes_")
myplot_colplot(xvar,yvar,mycols_sub, addCurve = TRUE,
               dt = tad_coexpr_fc_DT[tad_coexpr_fc_DT$cmpType == "subtypes",] , outPrefix = "subtypes_")

myplot_densplot(xvar,yvar, addCurve = TRUE, 
                dt = tad_coexpr_fc_DT[tad_coexpr_fc_DT$cmpType == "norm_vs_tumor",] , outPrefix = "norm_vs_tumor_")
myplot_colplot(xvar,yvar,mycols_tumor, addCurve = TRUE,
               dt = tad_coexpr_fc_DT[tad_coexpr_fc_DT$cmpType == "norm_vs_tumor",] , outPrefix = "norm_vs_tumor_")

myplot_densplot(xvar,yvar, addCurve = TRUE, 
                dt = tad_coexpr_fc_DT[tad_coexpr_fc_DT$cmpType == "wt_vs_mut",] , outPrefix = "wt_vs_mut_")
myplot_colplot(xvar,yvar,mycols_mut, addCurve = TRUE,
               dt = tad_coexpr_fc_DT[tad_coexpr_fc_DT$cmpType == "wt_vs_mut",] , outPrefix = "wt_vs_mut_")


##########################
### detect "disruption": those with change (or remain cohesive ?) in coexpr => betw-within cond2 vs betw-within cond1
##########################
yvar <- "withinBetweenNbrDiffCond2"
xvar <- "withinBetweenNbrDiffCond1"

myplot_densplot(xvar,yvar, addCurve = TRUE)
myplot_colplot(xvar,yvar,mycols, addCurve = TRUE)

myplot_densplot(xvar,yvar, addCurve = TRUE, 
                dt = tad_coexpr_fc_DT[tad_coexpr_fc_DT$cmpType == "subtypes",] , outPrefix = "subtypes_")
myplot_colplot(xvar,yvar,mycols_sub, addCurve = TRUE,
               dt = tad_coexpr_fc_DT[tad_coexpr_fc_DT$cmpType == "subtypes",] , outPrefix = "subtypes_")

myplot_densplot(xvar,yvar, addCurve = TRUE, 
                dt = tad_coexpr_fc_DT[tad_coexpr_fc_DT$cmpType == "norm_vs_tumor",] , outPrefix = "norm_vs_tumor_")
myplot_colplot(xvar,yvar,mycols_tumor, addCurve = TRUE,
               dt = tad_coexpr_fc_DT[tad_coexpr_fc_DT$cmpType == "norm_vs_tumor",] , outPrefix = "norm_vs_tumor_")

myplot_densplot(xvar,yvar, addCurve = TRUE, 
                dt = tad_coexpr_fc_DT[tad_coexpr_fc_DT$cmpType == "wt_vs_mut",] , outPrefix = "wt_vs_mut_")
myplot_colplot(xvar,yvar,mycols_mut, addCurve = TRUE,
               dt = tad_coexpr_fc_DT[tad_coexpr_fc_DT$cmpType == "wt_vs_mut",] , outPrefix = "wt_vs_mut_")


##########################
### detect coordinated change in expression and DE: those with high FC and high coexpr => coexpr vs. FC expr
##########################
yvar <- "withinBetwNbrLogFC"
xvar <- "meanFC"

myplot_densplot(xvar,yvar)
myplot_colplot(xvar,yvar,mycols)

##########################
### detect those that function as regulatory unit and DE: those with high FC and high diff. between-within => betw-within vs. FC expr
##########################

yvar <- "withinCoexpr"
xvar <- "meanFC"

myplot_densplot(xvar,yvar)
myplot_colplot(xvar,yvar,mycols)

yvar <- "withinCoexpr"
xvar <- "ratioDown"

myplot_densplot(xvar,yvar)
myplot_colplot(xvar,yvar,mycols)



yvar <- "withinBetweenDiffNbr"
xvar <- "meanFC"

myplot_densplot(xvar,yvar)
myplot_colplot(xvar,yvar,mycols)


##########################
# CPTMT SCORE
##########################

score_DT$hicds <- score_DT$dataset
tad_coexpr_fc_DT$hicds  <- dirname(tad_coexpr_fc_DT$dataset)

all_DT <- merge(score_DT, tad_coexpr_fc_DT, by = c("hicds", "region"))

tad_coexpr_fc_DT <- all_DT

head(all_DT)

stopifnot(!is.na(all_DT$score))

yvar <- "meanFC"
xvar <- "score"
myplot_densplot(xvar,yvar)


yvar <- "withinCoexpr"
xvar <- "score"
myplot_densplot(xvar,yvar)

yvar <- "betweenAllCoexpr"
xvar <- "score"
myplot_densplot(xvar,yvar)

yvar <- "betweenKbCoexpr"
xvar <- "score"
myplot_densplot(xvar,yvar)

yvar <- "betweenNbrCoexpr"
xvar <- "score"
myplot_densplot(xvar,yvar)



# ######################################################################################
# ######################################################################################
# ######################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))





