# Rscript coexpr_DE_queryTAD.R

script_name <- "coexpr_DE_queryTAD.R"

startTime <- Sys.time()

cat("> START coexpr_DE_queryTAD.R \n")

SSHFS <- FALSE

setDir <- ifelse(SSHFS, "/media/electron", "")

require(tools)
require(foreach)
require(doMC)
registerDoMC(ifelse(SSHFS, 2, 40))

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")
source("subtype_cols.R")

source("look_tads_fct.R")

source("plot_lolliTAD_funct.R")

nTop <- 5

plotType <- "png"
myHeight <- ifelse(plotType=="png", 400, 7)
myWidth <- myHeight
plotCex <- 1.4

myCexAxis <- 1.2
myCexLab <- 1.2

windowSizeBp <- 500*10^3
options(scipen=100)


outFolder <- paste0("COEXPR_DE_QUERYTAD_", nTop)
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

script0_name <- "0_prepGeneData"
  

entrez2symb_dt <- read.delim(file.path(setDir,
                                       "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt"),
                             header=T, stringsAsFactors = FALSE)
entrez2symb_dt$entrezID <- as.character(entrez2symb_dt$entrezID)



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
# load(outFile)



tad_coexpr_fc_DT$meanFC_rank <- rank(-abs(tad_coexpr_fc_DT$meanFC), ties="min") # rank: highest rank = highest coexpr value


tad_coexpr_fc_DT$withinBetweenDiffNbr_rank <- rank(-tad_coexpr_fc_DT$withinBetweenDiffNbr, ties = "min") # rank: highest rank = highest coexpr value
tad_coexpr_fc_DT$withinBetweenDiffNbr_meanFC_avgRank <- (tad_coexpr_fc_DT$withinBetweenDiffNbr_rank+tad_coexpr_fc_DT$meanFC_rank)/2

tad_coexpr_fc_DT$withinBetweenDiffKb_rank <- rank(-tad_coexpr_fc_DT$withinBetweenDiffKb, ties = "min") # rank: highest rank = highest coexpr value
tad_coexpr_fc_DT$withinBetweenDiffKb_meanFC_avgRank <- (tad_coexpr_fc_DT$withinBetweenDiffKb_rank+tad_coexpr_fc_DT$meanFC_rank)/2


tad_coexpr_fc_DT$withinCoexpr_rank <- rank(-tad_coexpr_fc_DT$withinCoexpr, ties="min") # rank: highest rank = highest coexpr value
tad_coexpr_fc_DT$withinCoexpr_meanFC_avgRank <- (tad_coexpr_fc_DT$withinCoexpr_rank+tad_coexpr_fc_DT$meanFC_rank)/2





#################################################################
################################################################# HIGH FC and HIGH WITHIN COEXPR
#################################################################

# iterate
# subset_vars (subsetVar) -> subset_levels (var) -> rankingVars (rank_var)

subset_vars <- c("cmpType", "all")
# cmpType = subtypes/norm vs tumor/wt vs mut
rankingVars <- c( "meanFC_rank",
                  "withinCoexpr_rank", "withinCoexpr_meanFC_avgRank",
                  "withinBetweenDiffNbr_rank", "withinBetweenDiffNbr_meanFC_avgRank",
                  "withinBetweenDiffKb_rank", "withinBetweenDiffKb_meanFC_avgRank",
                  "withinCoexpr_rank", "withinCoexpr_meanFC_avgRank")


#### > START ITERATING AND OUTPUTING


subsetVar = subset_vars[1]
for(subsetVar in subset_vars) {
  if(subsetVar == "all") {
    subset_levels <- ""
  } else {
    subset_levels <- unique(tad_coexpr_fc_DT[, subsetVar])  
  }
  var = subset_levels[1]
  for(var in subset_levels) {
    if(subsetVar == "all") {
      stopifnot(var == "")
      curr_DT <- tad_coexpr_fc_DT
    } else {
        curr_DT <- tad_coexpr_fc_DT[tad_coexpr_fc_DT[, paste0(subsetVar)] == var,]
    }
    stopifnot(nrow(curr_DT) > 0)
    rank_var = rankingVars[1]
    for(rank_var in rankingVars){
      stopifnot(rank_var %in% colnames(curr_DT))
      ranked_curr_DT <- curr_DT[order(curr_DT[,rank_var]),]
      
      
      outFile <- file.path(outFolder, paste0(subsetVar, "_", var, "_", rank_var, "_nTop", nTop, "_densplot.", plotType ))
      do.call(plotType, list(outFile, height=myHeight, width=myWidth))
      densplot(
        y = curr_DT$withinCoexpr,
        x = curr_DT$meanFC,
        cex.lab=myCexLab,
        cex.axis=myCexAxis,
        ylab="withinCoexpr",
        xlab = "meanFC",
        cex = 0.5,
        main = paste0(subsetVar, " - ", var, " - ", rank_var)
      )
      points(x = ranked_curr_DT$meanFC[1:nTop],
             y = ranked_curr_DT$withinCoexpr[1:nTop],
            pch = 4,
            col="red",
            cex=1.5
             )
      text(x = ranked_curr_DT$meanFC[1:nTop],
             y = ranked_curr_DT$withinCoexpr[1:nTop],
           labels = paste0(ranked_curr_DT$dataset[1:nTop],"\n", ranked_curr_DT$region[1:nTop]),
             pch = 4,
             col="red",
             cex=0.7
      )
      foo <- dev.off()
      cat(paste0("... written: ", outFile, "\n"))
      
      
      i=1
      # topTable_DT <- foreach(i = 1:nTop, .combine='rbind') %dopar% {
      plotList <- list()
        topTable_DT <- foreach(i = 1:nTop, .combine='rbind') %do% {
        i_hicds <- dirname(ranked_curr_DT$dataset[i])
        i_exprds <- basename(ranked_curr_DT$dataset[i])
        i_tad <- basename(ranked_curr_DT$region[i])
        
        
        plotList[[i]] <- plot_lolliTAD_ds(exprds = i_exprds,
                         hicds = i_hicds,
                         all_TADs = i_tad,
                         orderByLolli = "startPos")
        
        
        tmpDT <- get_tad_symbols_DT(hicds = i_hicds,
                           exprds = i_exprds,
                           region = i_tad,
                           mainPipFolder = pipOutFolder
                           )
        otherCols <- colnames(tmpDT)
        tmpDT$iTop <- i
        tmpDT[,c("iTop", otherCols)]
      } # end-foreach iterating over the "n" top TADs
      # subset_vars (subsetVar) -> subset_levels (var) -> rankingVars (rank_var)
      outFile <- file.path(outFolder, paste0(subsetVar, "_", var, "_", rank_var, "_nTop", nTop, ".txt" ))
      cat(outFile, "\n")
      write.table(topTable_DT, col.names=TRUE, row.names=FALSE, sep="\t", append=FALSE, quote=FALSE, file = outFile)
      cat(paste0("... written: ", outFile, "\n"))
      
      mytit <- paste0(subsetVar, " - ", var, " - ", rank_var, " - ", nTop)
      all_plots <- do.call(grid.arrange, c(plotList,  list(ncol=ifelse(nTop == 1, 1, 2), top=textGrob(mytit, gp=gpar(fontsize=20,font=2)))))
      

      outFile <- file.path(outFolder, paste0(subsetVar, "_", var, "_", rank_var, "_nTop", nTop, ".", plotType ))
      
      outWidthGG <- 20
      outHeightGG <- min(c(7 * nTop/2, 49))
      
      ggsave(filename = outFile, all_plots, width=outWidthGG, height = outHeightGG)
      cat("... written: ", outFile, "\n")
      
        
      
    } # end-for iterating over rankingVar # => rank_var in rankingVars (e.g. "withinCoexpr_rank")
  } # end-for iterating over levels of data subset #   for(var in subset_levels) (e.g. "subtypes")
} # end-for iterating over subset data # for(subsetVar in subset_vars) (e.g. "cmpType")






# ######################################################################################
# ######################################################################################
# ######################################################################################
cat(paste0("*** DONE: ", script_name, "\n"))
cat(paste0(startTime, "\n", Sys.time(), "\n"))





