# Rscript coexpr_DE_queryTAD.R

script_name <- "coexpr_DE_queryTAD.R"

startTime <- Sys.time()

cat("> START coexpr_DE_queryTAD.R \n")

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


outFolder <- "COEXPR_DE_QUERYTAD"
dir.create(outFolder, recursive=TRUE)

dataFolder <- "COEXPR_BETWEEN_WITHIN_ALL"

corrMet <- "pearson"
pipOutFolder <- file.path("PIPELINE", "OUTPUT_FOLDER")
stopifnot(dir.exists(pipOutFolder))

all_fc_files <- list.files(pipOutFolder, recursive = TRUE, pattern="all_meanLogFC_TAD.Rdata", full.names = FALSE)
stopifnot(length(all_fc_files) > 0)


dataFile <- file.path(dataFolder, "allData_within_between_coexpr.Rdata")
stopifnot(file.exists(dataFile))
allData_within_between_coexpr <- eval(parse(text = load(dataFile)))

all_domainScore_files <- list.files(".", recursive = TRUE, pattern="_final_domains_withScore.txt", full.names = FALSE)
stopifnot(length(all_domainScore_files) > 0)

all_ratioDown_files <- list.files(pipOutFolder, recursive = TRUE, pattern="all_obs_ratioDown.Rdata", full.names = FALSE)
stopifnot(length(all_ratioDown_files) > 0)

  
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


fc_rd_DT <- merge(rD_DT, fc_DT, by=c("dataset", "region"))
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


#################################################################
################################################################# HIGH FC and HIGH WITHIN COEXPR
#################################################################


# by cmpTypes
all_cmps <- unique(tad_coexpr_fc_DT$cmpType)


rankingVars <- c("withinCoexpr_rank", "meanFC_rank", "withinCoexpr_meanFC_avgRank")

nTop <- 5


cmp=all_cmps[1]
for(cmp in all_cmps){
  
  cmp_DT <- tad_coexpr_fc_DT[tad_coexpr_fc_DT$cmps == cmp,]
  
  cmp_DT$withinCoexpr_rank <- rank(-cmp_DT$withinCoexpr, ties="min") # rank: highest rank = highest coexpr value
  cmp_DT$meanFC_rank <- rank(-abs(cmp_DT$meanFC), ties="min") # rank: highest rank = highest coexpr value
  
  cmp_DT$withinCoexpr_meanFC_avgRank <- (cmp_DT$withinCoexpr_rank+cmp_DT$meanFC_rank)/2
  
  # plot(
  #   x= cmp_DT$withinCoexpr_rank,
  #   y= cmp_DT$meanFC_rank
  # )
  
  for(var in rankingVars) {
    
    
    
    
    
  }
  
  
}












# ######################################################################################
# ######################################################################################
# ######################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))





