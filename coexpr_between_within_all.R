
# Rscript coexpr_between_within_all.R

script_name <- "coexpr_between_within_all.R"

startTime <- Sys.time()

cat("> START coexpr_between_within_all.R \n")

SSHFS <- FALSE

buildData <- TRUE

require(tools)
require(foreach)
require(doMC)
registerDoMC(ifelse(SSHFS, 2, 80))

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

plotType <- "png"
myHeight <- ifelse(plotType=="png", 500, 7)
myWidth <- myHeight
plotCex <- 1.4

myCexAxis <- 1.2
myCexLab <- 1.2

windowSizeBp <- 500*10^3
options(scipen=100)


outFolder <- "COEXPR_BETWEEN_WITHIN_ALL"
dir.create(outFolder, recursive=TRUE)

corrMet <- "pearson"
script0_name <- "0_prepGeneData"
mypattern <- "coexprDT.Rdata"
pipOutFolder <- file.path("PIPELINE", "OUTPUT_FOLDER")
stopifnot(dir.exists(pipOutFolder))

all_fc_files <- list.files(pipOutFolder, recursive = TRUE, pattern="all_meanLogFC_TAD.Rdata", full.names = FALSE)
stopifnot(length(all_fc_files) > 0)

coexprFiles  <- list.files("CREATE_COEXPR_SORTNODUP", recursive = TRUE, pattern = mypattern, full.names = TRUE)
stopifnot(length(coexprFiles) > 0)

coexprByCondFolder <- "CREATE_COEXPR_BYCOND_SORTNODUP"
stopifnot(dir.exists(coexprByCondFolder))

coexprDataFile = coexprFiles[1]
hicds = "K562_40kb"
exprds = "TCGAlaml_wt_mutFLT3"
coexprDataFile = file.path("CREATE_COEXPR_SORTNODUP", 
                            hicds, 
                            paste0(exprds),
                            corrMet,
                            "coexprDT.Rdata")



#### LOAD INFORMATION ABOUT GENES AROUND TAD (SAME NUMBER)
aroundTADfile <- file.path("CREATE_SAMPLE_AROUND_TADS", "all_ds_sample_around_TADs.Rdata")
stopifnot(file.exists(aroundTADfile))
all_ds_sample_around_TADs <- eval(parse(text = load(aroundTADfile)))


#### LOAD INFORMATION ABOUT GENES AROUND TAD (WINDOW)

aroundKbTADfile <- file.path("CREATE_SAMPLE_AROUNDKB_TADS", paste0(windowSizeBp), "all_ds_sample_aroundKb_TADs.Rdata")
stopifnot(file.exists(aroundKbTADfile))
all_ds_sample_aroundKb_TADs <- eval(parse(text = load(aroundKbTADfile)))


coexprFilesOutNames <- coexprFiles

coexprDataFile = coexprFiles[1]

stopifnot(length(coexprFiles) == length(all_ds_sample_aroundKb_TADs))
stopifnot(length(coexprFiles) == length(all_ds_sample_around_TADs))

if(buildData) {

  allData_within_between_coexpr <- foreach(coexprDataFile = coexprFiles) %do% {
    
    cat("... start ", coexprDataFile, " \n")
    
    stopifnot(file.exists(coexprDataFile))
    
    hicds <- basename(dirname(dirname(dirname(coexprDataFile))))
    stopifnot(dir.exists(hicds))
    exprds <- basename(dirname(dirname(coexprDataFile)))
    
    curr_coexprByCondFolder <- file.path(coexprByCondFolder, hicds, exprds, corrMet)
    stopifnot(dir.exists(curr_coexprByCondFolder))

    coexprDT_cond1File <- file.path(curr_coexprByCondFolder , "coexprDT_cond1.Rdata")
    stopifnot(file.exists(coexprDT_cond1File))
    coexprDT_cond2File <- file.path(curr_coexprByCondFolder , "coexprDT_cond2.Rdata")
    stopifnot(file.exists(coexprDT_cond2File))
    
    
    cat("... start ", hicds, " - ", exprds, "\n")
    
    cat("... load coexprDT ...\n")
    coexprDT <- eval(parse(text = load(coexprDataFile)))
    coexprDT$gene1 <- as.character(coexprDT$gene1)
    coexprDT$gene2 <- as.character(coexprDT$gene2)
    
    
    cat("... load coexprDT_cond1 ...\n")
    coexprDT_cond1 <- eval(parse(text = load(paste0(coexprDT_cond1File))))
    coexprDT_cond1$gene1 <- as.character(coexprDT_cond1$gene1)
    coexprDT_cond1$gene2 <- as.character(coexprDT_cond1$gene2)
    
    
    cat("... load coexprDT_cond2 ...\n")
    coexprDT_cond2 <- eval(parse(text = load(paste0(coexprDT_cond2File))))
    coexprDT_cond2$gene1 <- as.character(coexprDT_cond2$gene1)
    coexprDT_cond2$gene2 <- as.character(coexprDT_cond2$gene2)
    

    g2tFile <- file.path(hicds, "genes2tad", "all_genes_positions.txt")
    stopifnot(file.exists(g2tFile))
    g2t_DT <- read.delim(g2tFile, header=F, 
                         col.names = c("entrezID",  "chromo", "start", "end", "region"), stringsAsFactors = FALSE)
    g2t_DT$entrezID <- as.character(g2t_DT$entrezID)
    
    dsPipOutDir <- file.path("PIPELINE", "OUTPUT_FOLDER", hicds, exprds)
    stopifnot(dir.exists(dsPipOutDir))
    stopifnot(dir.exists(file.path(dsPipOutDir, script0_name)))
    geneListFile <- file.path(dsPipOutDir, script0_name, "pipeline_geneList.Rdata")
    stopifnot(file.exists(geneListFile))
    pipeline_geneList <- eval(parse(text = load(geneListFile))) # 
    
    regionListFile <- file.path(dsPipOutDir, script0_name, "pipeline_regionList.Rdata")
    stopifnot(file.exists(regionListFile))
    pipeline_regionList <- eval(parse(text = load(regionListFile))) # 
    
    stopifnot(pipeline_geneList %in% g2t_DT$entrezID)
    g2t_DT <- g2t_DT[g2t_DT$entrezID %in% pipeline_geneList & g2t_DT$region %in% pipeline_regionList, ]
    stopifnot(g2t_DT$entrezID %in% coexprDT$gene1 | g2t_DT$entrezID %in% coexprDT$gene2 )
    stopifnot(pipeline_regionList %in% g2t_DT$region)
    
    reg = pipeline_regionList[1]

    within_between_coexpr_data <- foreach(reg = pipeline_regionList) %dopar% {
      
      cat(paste0("...... start ", hicds, " - ", exprds, " - ", reg, " \n"))
      tad_genes <- g2t_DT$entrezID[g2t_DT$region == reg]
      
      ####################################################################
      #### ALL CONDITIONS
      ####################################################################
      tad_coexprDT <- coexprDT[coexprDT$gene1 %in% tad_genes | 
                                 coexprDT$gene2 %in% tad_genes,]
    
      stopifnot(nrow(tad_coexprDT) > 0)  
      stopifnot(!is.na(tad_coexprDT$coexpr))
      stopifnot(is.numeric(tad_coexprDT$coexpr))
      stopifnot(tad_genes %in% tad_coexprDT$gene1 | tad_genes %in% tad_coexprDT$gene2 )
      
      ### WITHIN
      within_coexprDT <- tad_coexprDT[tad_coexprDT$gene1 %in% tad_genes & tad_coexprDT$gene2 %in% tad_genes,]
      
      ### BETWEEN: ALL
      # betweenAll: 1 of the 2 is not in TAD
      betweenAll_coexprDT <- tad_coexprDT[ ! (tad_coexprDT$gene1 %in% tad_genes & tad_coexprDT$gene2 %in% tad_genes),]
      stopifnot(nrow(within_coexprDT) + nrow(betweenAll_coexprDT) == nrow(tad_coexprDT))
      
      stopifnot(betweenAll_coexprDT$gene1 %in% tad_genes | 
                  betweenAll_coexprDT$gene2 %in% tad_genes )
      
      stopifnot( ! betweenAll_coexprDT$gene1 %in% tad_genes | 
                  ! betweenAll_coexprDT$gene2 %in% tad_genes )
      
      ### BETWEEN: GENES AROUND - NUMBER
      betweenGenes_windowNbr <- all_ds_sample_around_TADs[[file.path(hicds, exprds)]][[paste0(reg)]]
      stopifnot(length(betweenGenes_windowNbr) > 0)
      
      betweenNbr_coexprDT <- tad_coexprDT[ tad_coexprDT$gene1 %in% betweenGenes_windowNbr | 
                                             tad_coexprDT$gene2 %in% betweenGenes_windowNbr,] 
      stopifnot( betweenNbr_coexprDT$gene1 %in% betweenGenes_windowNbr |
                   betweenNbr_coexprDT$gene1 %in% tad_genes )
      
      stopifnot( betweenNbr_coexprDT$gene2 %in% betweenGenes_windowNbr |
                   betweenNbr_coexprDT$gene2 %in% tad_genes )
      
      stopifnot( ! (betweenNbr_coexprDT$gene1 %in% tad_genes & 
                      betweenNbr_coexprDT$gene2 %in% tad_genes) ) 
                      
      stopifnot(betweenGenes_windowNbr %in% betweenAll_coexprDT$gene1 | 
                  betweenGenes_windowNbr %in% betweenAll_coexprDT$gene2)
      
      stopifnot(!betweenGenes_windowNbr %in% within_coexprDT$gene1 &
                  ! betweenGenes_windowNbr %in% within_coexprDT$gene2)
      
      ### BETWEEN: GENES AROUND - KB
      betweenGenes_windowKb <- all_ds_sample_aroundKb_TADs[[file.path(hicds, exprds)]][[paste0(reg)]]
      # stopifnot(length(betweenGenes_windowKb) > 0) => not always TRUE
      
      betweenKb_coexprDT <- tad_coexprDT[ tad_coexprDT$gene1 %in% betweenGenes_windowKb | 
                                             tad_coexprDT$gene2 %in% betweenGenes_windowKb,] 
      stopifnot( betweenKb_coexprDT$gene1 %in% betweenGenes_windowKb |
                   betweenKb_coexprDT$gene1 %in% tad_genes )
      
      stopifnot( betweenKb_coexprDT$gene2 %in% betweenGenes_windowKb |
                   betweenKb_coexprDT$gene2 %in% tad_genes )
      
      stopifnot( ! (betweenKb_coexprDT$gene1 %in% tad_genes & 
                      betweenKb_coexprDT$gene2 %in% tad_genes) ) 
      
      stopifnot(betweenGenes_windowKb %in% betweenAll_coexprDT$gene1 | 
                  betweenGenes_windowKb %in% betweenAll_coexprDT$gene2)
      
      stopifnot(!betweenGenes_windowKb %in% within_coexprDT$gene1 &
                  ! betweenGenes_windowKb %in% within_coexprDT$gene2)
      
      ####################################################################
      #### COND 1
      ####################################################################
      
      tad_coexprDT_cond1 <- coexprDT_cond1[coexprDT_cond1$gene1 %in% tad_genes | 
                                 coexprDT_cond1$gene2 %in% tad_genes,]
      
      stopifnot(nrow(tad_coexprDT_cond1) > 0)  
      stopifnot(!is.na(tad_coexprDT_cond1$coexpr))
      stopifnot(is.numeric(tad_coexprDT_cond1$coexpr))
      stopifnot(tad_genes %in% tad_coexprDT_cond1$gene1 | tad_genes %in% tad_coexprDT_cond1$gene2 )
      
      ### WITHIN
      
      within_coexprDT_cond1 <- tad_coexprDT_cond1[tad_coexprDT_cond1$gene1 %in% tad_genes & tad_coexprDT_cond1$gene2 %in% tad_genes,]
      
      
      
      ### BETWEEN: ALL
      
      
      betweenAll_coexprDT_cond1 <- tad_coexprDT_cond1[ ! (tad_coexprDT_cond1$gene1 %in% tad_genes & tad_coexprDT_cond1$gene2 %in% tad_genes),]
      
      stopifnot(nrow(within_coexprDT_cond1) + nrow(betweenAll_coexprDT_cond1) == nrow(tad_coexprDT_cond1))
      
      
      stopifnot(betweenAll_coexprDT_cond1$gene1 %in% tad_genes | 
                  betweenAll_coexprDT_cond1$gene2 %in% tad_genes )
      
      stopifnot( ! betweenAll_coexprDT_cond1$gene1 %in% tad_genes | 
                   ! betweenAll_coexprDT_cond1$gene2 %in% tad_genes )
      
      
      
      
      ### BETWEEN: GENES AROUND - NBR
    
      
      
      betweenNbr_coexprDT_cond1 <- tad_coexprDT_cond1[ tad_coexprDT_cond1$gene1 %in% betweenGenes_windowNbr | 
                                             tad_coexprDT_cond1$gene2 %in% betweenGenes_windowNbr,] 
      stopifnot( betweenNbr_coexprDT_cond1$gene1 %in% betweenGenes_windowNbr |
                   betweenNbr_coexprDT_cond1$gene1 %in% tad_genes )
      
      stopifnot( betweenNbr_coexprDT_cond1$gene2 %in% betweenGenes_windowNbr |
                   betweenNbr_coexprDT_cond1$gene2 %in% tad_genes )
      
      stopifnot( ! (betweenNbr_coexprDT_cond1$gene1 %in% tad_genes & 
                      betweenNbr_coexprDT_cond1$gene2 %in% tad_genes) ) 
      
      stopifnot(betweenGenes_windowNbr %in% betweenAll_coexprDT$gene1 | 
                  betweenGenes_windowNbr %in% betweenAll_coexprDT$gene2)
      
      stopifnot(!betweenGenes_windowNbr %in% within_coexprDT$gene1 &
                  ! betweenGenes_windowNbr %in% within_coexprDT$gene2)
      
      
      
      
      
      ### BETWEEN: GENES AROUND - KB
      
      betweenKb_coexprDT_cond1 <- tad_coexprDT_cond1[ tad_coexprDT_cond1$gene1 %in% betweenGenes_windowKb | 
                                            tad_coexprDT_cond1$gene2 %in% betweenGenes_windowKb,] 
      stopifnot( betweenKb_coexprDT_cond1$gene1 %in% betweenGenes_windowKb |
                   betweenKb_coexprDT_cond1$gene1 %in% tad_genes )
      
      stopifnot( betweenKb_coexprDT_cond1$gene2 %in% betweenGenes_windowKb |
                   betweenKb_coexprDT_cond1$gene2 %in% tad_genes )
      
      stopifnot( ! (betweenKb_coexprDT_cond1$gene1 %in% tad_genes & 
                      betweenKb_coexprDT_cond1$gene2 %in% tad_genes) ) 
      
      stopifnot(betweenGenes_windowKb %in% betweenAll_coexprDT$gene1 | 
                  betweenGenes_windowKb %in% betweenAll_coexprDT$gene2)
      
      stopifnot(!betweenGenes_windowKb %in% within_coexprDT$gene1 &
                  ! betweenGenes_windowKb %in% within_coexprDT$gene2)
      
      
      
      
      ####################################################################
      #### COND 2
      ####################################################################

      tad_coexprDT_cond2 <- coexprDT_cond2[coexprDT_cond2$gene1 %in% tad_genes | 
                                             coexprDT_cond2$gene2 %in% tad_genes,]
      
      stopifnot(nrow(tad_coexprDT_cond2) > 0)  
      stopifnot(!is.na(tad_coexprDT_cond2$coexpr))
      stopifnot(is.numeric(tad_coexprDT_cond2$coexpr))
      stopifnot(tad_genes %in% tad_coexprDT_cond2$gene1 | tad_genes %in% tad_coexprDT_cond2$gene2 )
      
      ### WITHIN
      within_coexprDT_cond2 <- tad_coexprDT_cond2[tad_coexprDT_cond2$gene1 %in% tad_genes & tad_coexprDT_cond2$gene2 %in% tad_genes,]
      
      ### BETWEEN: ALL
      
      betweenAll_coexprDT_cond2 <- tad_coexprDT_cond2[ ! (tad_coexprDT_cond2$gene1 %in% tad_genes & tad_coexprDT_cond2$gene2 %in% tad_genes),]
      
      stopifnot(nrow(within_coexprDT_cond2) + nrow(betweenAll_coexprDT_cond2) == nrow(tad_coexprDT_cond2))
      
      stopifnot(betweenAll_coexprDT_cond2$gene1 %in% tad_genes | 
                  betweenAll_coexprDT_cond2$gene2 %in% tad_genes )
      
      stopifnot( ! betweenAll_coexprDT_cond2$gene1 %in% tad_genes | 
                   ! betweenAll_coexprDT_cond2$gene2 %in% tad_genes )
      
      
      ### BETWEEN: GENES AROUND - NBR
      
      
      betweenNbr_coexprDT_cond2 <- tad_coexprDT_cond2[ tad_coexprDT_cond2$gene1 %in% betweenGenes_windowNbr | 
                                             tad_coexprDT_cond2$gene2 %in% betweenGenes_windowNbr,] 
      stopifnot( betweenNbr_coexprDT_cond2$gene1 %in% betweenGenes_windowNbr |
                   betweenNbr_coexprDT_cond2$gene1 %in% tad_genes )
      
      stopifnot( betweenNbr_coexprDT_cond2$gene2 %in% betweenGenes_windowNbr |
                   betweenNbr_coexprDT_cond2$gene2 %in% tad_genes )
      
      stopifnot( ! (betweenNbr_coexprDT_cond2$gene1 %in% tad_genes & 
                      betweenNbr_coexprDT_cond2$gene2 %in% tad_genes) ) 
      
      stopifnot(betweenGenes_windowNbr %in% betweenAll_coexprDT$gene1 | 
                  betweenGenes_windowNbr %in% betweenAll_coexprDT$gene2)
      
      stopifnot(!betweenGenes_windowNbr %in% within_coexprDT$gene1 &
                  ! betweenGenes_windowNbr %in% within_coexprDT$gene2)
      
      
      
      
      ### BETWEEN: GENES AROUND - KB
      
      
      betweenKb_coexprDT_cond2 <- tad_coexprDT_cond2[ tad_coexprDT_cond2$gene1 %in% betweenGenes_windowKb | 
                                            tad_coexprDT_cond2$gene2 %in% betweenGenes_windowKb,] 
      stopifnot( betweenKb_coexprDT_cond2$gene1 %in% betweenGenes_windowKb |
                   betweenKb_coexprDT_cond2$gene1 %in% tad_genes )
      
      stopifnot( betweenKb_coexprDT_cond2$gene2 %in% betweenGenes_windowKb |
                   betweenKb_coexprDT_cond2$gene2 %in% tad_genes )
      
      stopifnot( ! (betweenKb_coexprDT_cond2$gene1 %in% tad_genes & 
                      betweenKb_coexprDT_cond2$gene2 %in% tad_genes) ) 
      
      stopifnot(betweenGenes_windowKb %in% betweenAll_coexprDT$gene1 | 
                  betweenGenes_windowKb %in% betweenAll_coexprDT$gene2)
      
      stopifnot(!betweenGenes_windowKb %in% within_coexprDT$gene1 &
                  ! betweenGenes_windowKb %in% within_coexprDT$gene2)
      
      
      
      
      
      
      #####################################################################################################
      
      list(withinCoexpr = mean(within_coexprDT$coexpr),
           withinCoexpr_cond1 = mean(within_coexprDT_cond1$coexpr),
           withinCoexpr_cond2 = mean(within_coexprDT_cond2$coexpr),
           betweenAllCoexpr = mean(betweenAll_coexprDT$coexpr),
           betweenAllCoexpr_cond1 = mean(betweenAll_coexprDT_cond1$coexpr),
           betweenAllCoexpr_cond2 = mean(betweenAll_coexprDT_cond2$coexpr),
           betweenKbCoexpr = mean(betweenKb_coexprDT$coexpr),
           betweenKbCoexpr_cond1 = mean(betweenKb_coexprDT_cond1$coexpr),
           betweenKbCoexpr_cond2 = mean(betweenKb_coexprDT_cond2$coexpr),
           betweenNbrCoexpr = mean(betweenNbr_coexprDT$coexpr),
           betweenNbrCoexpr_cond1 = mean(betweenNbr_coexprDT_cond1$coexpr),
           betweenNbrCoexpr_cond2 = mean(betweenNbr_coexprDT_cond2$coexpr)
           )
    } # end-foreach iterating over TADs
    names(within_between_coexpr_data) <- pipeline_regionList
    within_between_coexpr_data
  } # end-foreach iterating over dataset
  names(allData_within_between_coexpr) <- coexprFilesOutNames
  
  outFile <- file.path(outFolder, "allData_within_between_coexpr.Rdata")
  save(allData_within_between_coexpr, file = outFile)
  cat(paste0("... written: ", outFile, "\n"))
  
} else { # end if-buildData
  outFile <- file.path(outFolder, "allData_within_between_coexpr.Rdata")
  allData_within_between_coexpr <- eval(parse(text = load(outFile)))
}

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

tad_coexpr_fc_DT <- merge(tad_coexpr_DT, fc_DT, by=c("dataset", "region"))
tad_coexpr_fc_DT <- tad_coexpr_fc_DT[order(tad_coexpr_fc_DT$withinCoexpr, decreasing = TRUE),]
tad_coexpr_fc_DT$TADrank <- 1:nrow(tad_coexpr_fc_DT)


tad_coexpr_fc_DT$withinDiffCond1Cond2 <- (tad_coexpr_fc_DT$withinCoexpr_cond1 - tad_coexpr_fc_DT$withinCoexpr_cond2)
tad_coexpr_fc_DT$withinRatioCond1Cond2 <- (tad_coexpr_fc_DT$withinCoexpr_cond1 / tad_coexpr_fc_DT$withinCoexpr_cond2)
tad_coexpr_fc_DT$withinChangeratioCond1Cond2 <- (tad_coexpr_fc_DT$withinCoexpr_cond2 - tad_coexpr_fc_DT$withinCoexpr_cond1)/tad_coexpr_fc_DT$withinCoexpr_cond1


tad_coexpr_fc_DT$withinBetweenDiffAll <- (tad_coexpr_fc_DT$withinCoexpr - tad_coexpr_fc_DT$betweenAllCoexpr) 
tad_coexpr_fc_DT$withinBetweenRatioAll <- (tad_coexpr_fc_DT$withinCoexpr / tad_coexpr_fc_DT$betweenAllCoexpr) 

tad_coexpr_fc_DT$withinBetweenDiffNbr <- (tad_coexpr_fc_DT$withinCoexpr - tad_coexpr_fc_DT$betweenNbrCoexpr) 
tad_coexpr_fc_DT$withinBetweenRatioNbr <- (tad_coexpr_fc_DT$withinCoexpr / tad_coexpr_fc_DT$betweenNbrCoexpr) 

tad_coexpr_fc_DT$withinBetweenDiffKb <- (tad_coexpr_fc_DT$withinCoexpr - tad_coexpr_fc_DT$betweenKbCoexpr) 
tad_coexpr_fc_DT$withinBetweenRatioKb <- (tad_coexpr_fc_DT$withinCoexpr / tad_coexpr_fc_DT$betweenKbCoexpr) 


outFile <- file.path(outFolder, paste0("multidens_withinCoexpr.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myHeight*1.2))
plot_multiDens(list(
  withinCoexpr = tad_coexpr_DT[,"withinCoexpr"],
  withinCoexpr_cond1 = tad_coexpr_DT[,"withinCoexpr_cond1"],
  withinCoexpr_cond2 = tad_coexpr_DT[,"withinCoexpr_cond2"]
), my_xlab = "withinCoexpr. corr."
)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


outFile <- file.path(outFolder, paste0("multidens_betweenCoexpr_all.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myHeight*1.2))
plot_multiDens(list(
  betweenAllCoexpr = tad_coexpr_DT[,"betweenAllCoexpr"],
  betweenNbrCoexpr = tad_coexpr_DT[,"betweenNbrCoexpr"],
  betweenKbCoexpr = tad_coexpr_DT[,"betweenKbCoexpr"],
  betweenAllCoexpr_cond1 = tad_coexpr_DT[,"betweenAllCoexpr_cond1"],
  betweenNbrCoexpr_cond1 = tad_coexpr_DT[,"betweenNbrCoexpr_cond1"],
  betweenKbCoexpr_cond1 = tad_coexpr_DT[,"betweenKbCoexpr_cond1"],
  betweenAllCoexpr_cond2 = tad_coexpr_DT[,"betweenAllCoexpr_cond2"],
  betweenNbrCoexpr_cond2 = tad_coexpr_DT[,"betweenNbrCoexpr_cond2"],
  betweenKbCoexpr_cond2= tad_coexpr_DT[,"betweenKbCoexpr_cond2"]
), my_xlab = "betweenCoexpr. corr. (all cond.)",
plotTit = "all between coexpr."
)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


outFile <- file.path(outFolder, paste0("multidens_betweenAllCoexpr_all.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myHeight*1.2))
plot_multiDens(list(
  betweenAllCoexpr = tad_coexpr_DT[,"betweenAllCoexpr"],
  # betweenNbrCoexpr = tad_coexpr_DT[,"betweenNbrCoexpr"],
  # betweenKbCoexpr = tad_coexpr_DT[,"betweenKbCoexpr"],
  betweenAllCoexpr_cond1 = tad_coexpr_DT[,"betweenAllCoexpr_cond1"],
  # betweenNbrCoexpr_cond1 = tad_coexpr_DT[,"betweenNbrCoexpr_cond1"],
  # betweenKbCoexpr_cond1 = tad_coexpr_DT[,"betweenKbCoexpr_cond1"],
  betweenAllCoexpr_cond2 = tad_coexpr_DT[,"betweenAllCoexpr_cond2"]#,
  # betweenNbrCoexpr_cond2 = tad_coexpr_DT[,"betweenNbrCoexpr_cond2"],
  # betweenKbCoexpr_cond2= tad_coexpr_DT[,"betweenKbCoexpr_cond2"]
), my_xlab = "betweenCoexpr. corr. (all cond.)",
plotTit = "betweenAllCoexpr"
)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

outFile <- file.path(outFolder, paste0("multidens_betweenKbCoexpr_all.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myHeight*1.2))
plot_multiDens(list(
  # betweenAllCoexpr = tad_coexpr_DT[,"betweenAllCoexpr"],
  # betweenNbrCoexpr = tad_coexpr_DT[,"betweenNbrCoexpr"],
  betweenKbCoexpr = tad_coexpr_DT[,"betweenKbCoexpr"],
  # betweenAllCoexpr_cond1 = tad_coexpr_DT[,"betweenAllCoexpr_cond1"],
  # betweenNbrCoexpr_cond1 = tad_coexpr_DT[,"betweenNbrCoexpr_cond1"],
  betweenKbCoexpr_cond1 = tad_coexpr_DT[,"betweenKbCoexpr_cond1"],
  # betweenAllCoexpr_cond2 = tad_coexpr_DT[,"betweenAllCoexpr_cond2"]#,
  # betweenNbrCoexpr_cond2 = tad_coexpr_DT[,"betweenNbrCoexpr_cond2"],
  betweenKbCoexpr_cond2= tad_coexpr_DT[,"betweenKbCoexpr_cond2"]
), my_xlab = "betweenCoexpr. corr. (all cond.)",
plotTit = "betweenKbCoexpr"
)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


outFile <- file.path(outFolder, paste0("multidens_betweenNbrCoexpr_all.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myHeight*1.2))
plot_multiDens(list(
  # betweenAllCoexpr = tad_coexpr_DT[,"betweenAllCoexpr"],
  betweenNbrCoexpr = tad_coexpr_DT[,"betweenNbrCoexpr"],
  # betweenKbCoexpr = tad_coexpr_DT[,"betweenKbCoexpr"],
  # betweenAllCoexpr_cond1 = tad_coexpr_DT[,"betweenAllCoexpr_cond1"],
  betweenNbrCoexpr_cond1 = tad_coexpr_DT[,"betweenNbrCoexpr_cond1"],
  # betweenKbCoexpr_cond1 = tad_coexpr_DT[,"betweenKbCoexpr_cond1"],
  # betweenAllCoexpr_cond2 = tad_coexpr_DT[,"betweenAllCoexpr_cond2"]#,
  betweenNbrCoexpr_cond2 = tad_coexpr_DT[,"betweenNbrCoexpr_cond2"]#,
  # betweenKbCoexpr_cond2= tad_coexpr_DT[,"betweenKbCoexpr_cond2"]
), my_xlab = "betweenCoexpr. corr. (all cond.)",
plotTit = "betweenNbrCoexpr"
)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

outFile <- file.path(outFolder, paste0("multidens_betweenCoexpr.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myHeight*1.2))
plot_multiDens(list(
  betweenAllCoexpr = tad_coexpr_DT[,"betweenAllCoexpr"],
  betweenNbrCoexpr = tad_coexpr_DT[,"betweenNbrCoexpr"],
  betweenKbCoexpr = tad_coexpr_DT[,"betweenKbCoexpr"]#,
  # betweenAllCoexpr_cond1 = tad_coexpr_DT[,"betweenAllCoexpr_cond1"],
  # betweenNbrCoexpr_cond1 = tad_coexpr_DT[,"betweenNbrCoexpr_cond1"],
  # betweenKbCoexpr_cond1 = tad_coexpr_DT[,"betweenKbCoexpr_cond1"],
  # betweenAllCoexpr_cond2 = tad_coexpr_DT[,"betweenAllCoexpr_cond2"]#,
  # betweenNbrCoexpr_cond2 = tad_coexpr_DT[,"betweenNbrCoexpr_cond2"]#,
  # betweenKbCoexpr_cond2= tad_coexpr_DT[,"betweenKbCoexpr_cond2"]
), my_xlab = "betweenCoexpr. corr. (all cond.)",
plotTit = "betweenCoexpr - cond1+2"
)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

outFile <- file.path(outFolder, paste0("multidens_betweenCoexpr_cond1.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myHeight*1.2))
plot_multiDens(list(
  # betweenAllCoexpr = tad_coexpr_DT[,"betweenAllCoexpr"],
  # betweenNbrCoexpr = tad_coexpr_DT[,"betweenNbrCoexpr"],
  # betweenKbCoexpr = tad_coexpr_DT[,"betweenKbCoexpr"]#,
  betweenAllCoexpr_cond1 = tad_coexpr_DT[,"betweenAllCoexpr_cond1"],
  betweenNbrCoexpr_cond1 = tad_coexpr_DT[,"betweenNbrCoexpr_cond1"],
  betweenKbCoexpr_cond1 = tad_coexpr_DT[,"betweenKbCoexpr_cond1"]#,
  # betweenAllCoexpr_cond2 = tad_coexpr_DT[,"betweenAllCoexpr_cond2"]#,
  # betweenNbrCoexpr_cond2 = tad_coexpr_DT[,"betweenNbrCoexpr_cond2"]#,
  # betweenKbCoexpr_cond2= tad_coexpr_DT[,"betweenKbCoexpr_cond2"]
), my_xlab = "betweenCoexpr. corr. (cond1)",
plotTit = "betweenCoexpr - cond1"
)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

outFile <- file.path(outFolder, paste0("multidens_betweenCoexpr_cond2.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myHeight*1.2))
plot_multiDens(list(
  # betweenAllCoexpr = tad_coexpr_DT[,"betweenAllCoexpr"],
  # betweenNbrCoexpr = tad_coexpr_DT[,"betweenNbrCoexpr"],
  # betweenKbCoexpr = tad_coexpr_DT[,"betweenKbCoexpr"]#,
  # betweenAllCoexpr_cond1 = tad_coexpr_DT[,"betweenAllCoexpr_cond1"],
  # betweenNbrCoexpr_cond1 = tad_coexpr_DT[,"betweenNbrCoexpr_cond1"],
  # betweenKbCoexpr_cond1 = tad_coexpr_DT[,"betweenKbCoexpr_cond1"]#,
  betweenAllCoexpr_cond2 = tad_coexpr_DT[,"betweenAllCoexpr_cond2"],
  betweenNbrCoexpr_cond2 = tad_coexpr_DT[,"betweenNbrCoexpr_cond2"],
  betweenKbCoexpr_cond2= tad_coexpr_DT[,"betweenKbCoexpr_cond2"]
), my_xlab = "betweenCoexpr. corr. (cond2)", 
plotTit = "betweenCoexpr - cond2"
)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))
all_y <- colnames(tad_coexpr_fc_DT)
all_y <- all_y[! all_y %in% c("dataset", "region", "meanFC", "TADrank")]
xvar <- "meanFC"
myx <- tad_coexpr_fc_DT[, xvar] 

for(yvar in all_y) {
    myy <- tad_coexpr_fc_DT[, yvar] 
    outFile <- file.path(outFolder, paste0(yvar, "_vs_", xvar, ".", plotType))
    do.call(plotType, list(outFile, height=myHeight, width=myWidth))
    densplot(y=myy,
             x=myx,
             pch=16,
             cex=0.7,
             cex.axis = myCexAxis,
             cex.lab = myCexLab,
             # col="black",
             xlab=xvar,
             ylab=yvar,
             main = paste0(yvar, " vs. ", xvar)
    )
    addCorr(x=myx,y=myy,legPos="topleft", bty='n')
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
    
    
    outFile <- file.path(outFolder, paste0("density_", yvar, ".", plotType))
    do.call(plotType, list(outFile, height=myHeight, width=myWidth*1.2))
    plot(density(na.omit(myy)),
             pch=16,
             cex=0.7,
             cex.axis = myCexAxis,
             cex.lab = myCexLab,
             # col="black",
             ylab=yvar,
             main = paste0(yvar)
    )
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
    
    
    
}
  
  



suffix=""
betwType="All"
for(suffix in c("", "_cond1", "_cond2")) {
  
  var1 <- paste0("withinCoexpr", suffix)
  stopifnot(var1 %in% colnames(tad_coexpr_DT))
  
  
  for(betwType in c("All", "Nbr", "Kb")) {
    
    var2 <- paste0("between", betwType, "Coexpr", suffix)
    
    stopifnot(var2 %in% colnames(tad_coexpr_DT))
    
    myy1 <- tad_coexpr_DT[, var1] 
    myy2 <- tad_coexpr_DT[, var2] 
    
    if(betwType == "Kb") betwType <- paste0(betwType, windowSizeBp)
    
    outFile <- file.path(outFolder, paste0("within_vs_between", betwType, suffix, ".", plotType))
    do.call(plotType, list(outFile, height=myHeight, width=myWidth))
    plot(y=myy1,
         x=tad_coexpr_DT$TADrank,
         ylim=c(range(myy1, myy2, na.rm = TRUE)),
         pch=16,
         cex=0.7,
         cex.axis = myCexAxis,
         cex.lab = myCexLab,
         col="black",
         xlab="TAD",
         ylab=paste0("mean pairwise ", toTitleCase(corrMet), " corr."),
         main = paste0("Between vs. within expression correlation")
    )
    points(y=myy2,
           x=tad_coexpr_DT$TADrank,
           pch=16,
           cex=0.7,
           col="red"
    )
    
    
    if(suffix == "") {
      mtext(side=3, text = paste0("Between: ", betwType, " - cond1+cond2"), font = 3)  
    } else {
      mtext(side=3, text = paste0("Between: ", betwType, " - ", gsub("_", "", suffix)), font = 3)
    }
    
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
    
  }
  
}


plotvar="withinCoexpr"

for(plotvar in c("withinCoexpr", "betweenAllCoexpr", "betweenKbCoexpr", "betweenNbrCoexpr")) {
  
  xvar <- paste0(plotvar, "_cond1")
  stopifnot(xvar %in% colnames(tad_coexpr_DT))
  
  yvar <- paste0(plotvar, "_cond2")
  stopifnot(yvar %in% colnames(tad_coexpr_DT))
  
  myx <- tad_coexpr_DT[, xvar] 
  myy <- tad_coexpr_DT[, yvar] 
  
  outFile <- file.path(outFolder, paste0("cond2_vs_cond1_", plotvar, ".", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  densplot(y=myy,
       x=myx,
       pch=16,
       cex=0.7,
       cex.axis = myCexAxis,
       cex.lab = myCexLab,
      # col="black",
       xlab="cond1",
       ylab="cond2",
       main = paste0("Expr. corr.: ", plotvar, " - cond2 vs. cond1")
  )
  addCorr(x=myx,y=myy,legPos="topleft", bty='n')
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
}

xvar <- "withinCoexpr"
stopifnot(xvar %in% colnames(tad_coexpr_DT))

all_y <- c("betweenAllCoexpr", "betweenKbCoexpr", "betweenNbrCoexpr") 

myx <- tad_coexpr_DT[, xvar] 

for(yvar in all_y) {
  
  stopifnot(yvar %in% colnames(tad_coexpr_DT))
  

  myy <- tad_coexpr_DT[, yvar] 
  
  outFile <- file.path(outFolder, paste0(yvar, "_vs_", xvar, ".", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  densplot(y=myy,
           x=myx,
           pch=16,
           cex=0.7,
           cex.axis = myCexAxis,
           cex.lab = myCexLab,
           # col="black",
           xlab=paste0(xvar),
           ylab=paste0(yvar),
           main = paste0("Expr. corr.: ", yvar, " vs. ", xvar)
  )
  addCorr(x=myx,y=myy,legPos="topleft", bty='n')
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
}



# ######################################################################################
# ######################################################################################
# ######################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))





