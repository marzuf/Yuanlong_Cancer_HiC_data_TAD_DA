
# Rscript cosine_between_within_all_diag1.R

script_name <- "cosine_between_within_all_diag1.R"

startTime <- Sys.time()

cat("> START cosine_between_within_all_diag1.R \n")

SSHFS <- FALSE

buildData <- TRUE

require(tools)
require(foreach)
require(doMC)
registerDoMC(ifelse(SSHFS, 2, 80))

library(lsa)
get_ang_dist <- function(vect1, vect2) {
  acos(cosine(vect1,vect2))/pi
}

library(igraph)


source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

plotType <- "png"
myHeight <- ifelse(plotType=="png", 500, 7)
myWidth <- myHeight
plotCex <- 1.4

myCexAxis <- 1.2
myCexLab <- 1.2

windowSizeBp <- 500*10^3
options(scipen=100)


outFolder <- "COSINE_BETWEEN_WITHIN_ALL_DIAG1"
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

# stopifnot(length(coexprFiles) == length(all_ds_sample_aroundKb_TADs))
# stopifnot(length(coexprFiles) == length(all_ds_sample_around_TADs))

if(buildData) {

  allData_ang_dist <- foreach(coexprDataFile = coexprFiles) %do% {
    
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
    
    # convert data frame to symmetric matrix (functions from igraph)

    ang_dist_data <- foreach(reg = pipeline_regionList) %dopar% {
      
      cat(paste0("...... start ", hicds, " - ", exprds, " - ", reg, " \n"))
      tad_genes <- g2t_DT$entrezID[g2t_DT$region == reg]
      
      ####################################################################
      #### ALL CONDITIONS
      ####################################################################
      
      ### WITHIN TAD ONLY
      tad_coexprDT <- coexprDT[coexprDT$gene1 %in% tad_genes & 
                                 coexprDT$gene2 %in% tad_genes,]
    
      stopifnot(nrow(tad_coexprDT) > 0)  
      stopifnot(!is.na(tad_coexprDT$coexpr))
      stopifnot(is.numeric(tad_coexprDT$coexpr))
      stopifnot(tad_genes %in% tad_coexprDT$gene1 | tad_genes %in% tad_coexprDT$gene2 )
      
      tad_coexprMatrix <- as.matrix(get.adjacency(graph.data.frame(tad_coexprDT, directed=FALSE), 
                                                  attr=colnames(tad_coexprDT)[3], sparse=FALSE))
      stopifnot(tad_genes %in% rownames(tad_coexprMatrix) & tad_genes %in% colnames(tad_coexprMatrix))
      stopifnot(rownames(tad_coexprMatrix) %in% tad_genes &  colnames(tad_coexprMatrix) %in% tad_genes)
      stopifnot( rownames(tad_coexprMatrix) == colnames(tad_coexprMatrix))
      ### CORRECTED HERE -> coexpression with itself should be 1 not zero !
      diag(tad_coexprMatrix) <- 1
      

      tad_genes_cmbs <- combn(x = tad_genes, m = 2)
      
      tad_ang_dist <- apply(tad_genes_cmbs, 2, function(x) {
        stopifnot(x[1] %in% rownames(tad_coexprMatrix))
        stopifnot(x[2] %in% rownames(tad_coexprMatrix))
        get_ang_dist(vect1 = tad_coexprMatrix[x[1],], vect2 = tad_coexprMatrix[x[2],])
      })
      
      tad_meanAngDist <- mean(tad_ang_dist)
      
      ### WITH ALL GENES
      all_coexprDT <- coexprDT[coexprDT$gene1 %in% tad_genes |
                                 coexprDT$gene2 %in% tad_genes,]
      
      stopifnot(nrow(all_coexprDT) > 0)  
      stopifnot(!is.na(all_coexprDT$coexpr))
      stopifnot(is.numeric(all_coexprDT$coexpr))
      stopifnot(tad_genes %in% all_coexprDT$gene1 | tad_genes %in% all_coexprDT$gene2 )
      
      all_coexprMatrix <- as.matrix(get.adjacency(graph.data.frame(all_coexprDT, directed=FALSE), 
                                                  attr=colnames(all_coexprDT)[3], sparse=FALSE))
      stopifnot(tad_genes %in% rownames(all_coexprMatrix) & tad_genes %in% colnames(all_coexprMatrix))
      stopifnot( rownames(all_coexprMatrix) == colnames(all_coexprMatrix))
      ### CORRECTED HERE -> coexpression with itself should be 1 not zero !
      diag(all_coexprMatrix) <- 1
      
      tad_genes_cmbs <- combn(x = tad_genes, m = 2)
      
      all_ang_dist <- apply(tad_genes_cmbs, 2, function(x) {
        stopifnot(x[1] %in% rownames(all_coexprMatrix))
        stopifnot(x[2] %in% rownames(all_coexprMatrix))
        get_ang_dist(vect1 = all_coexprMatrix[x[1],], vect2 = all_coexprMatrix[x[2],])
      })
      
      all_meanAngDist <- mean(all_ang_dist)
      
      ####################################################################
      #### COND 1
      ####################################################################
      
      ### WITHIN TAD ONLY
      tad_coexprDT_cond1 <- coexprDT_cond1[coexprDT_cond1$gene1 %in% tad_genes & 
                                 coexprDT_cond1$gene2 %in% tad_genes,]
      
      stopifnot(nrow(tad_coexprDT_cond1) > 0)  
      stopifnot(!is.na(tad_coexprDT_cond1$coexpr))
      stopifnot(is.numeric(tad_coexprDT_cond1$coexpr))
      stopifnot(tad_genes %in% tad_coexprDT_cond1$gene1 | tad_genes %in% tad_coexprDT_cond1$gene2 )
      
      
      tad_coexprMatrix_cond1 <- as.matrix(get.adjacency(graph.data.frame(tad_coexprDT_cond1, directed=FALSE), 
                                                  attr=colnames(tad_coexprDT_cond1)[3], sparse=FALSE))
      stopifnot(tad_genes %in% rownames(tad_coexprMatrix_cond1) & tad_genes %in% colnames(tad_coexprMatrix_cond1))
      stopifnot(rownames(tad_coexprMatrix_cond1) %in% tad_genes &  colnames(tad_coexprMatrix_cond1) %in% tad_genes)
      stopifnot( rownames(tad_coexprMatrix_cond1) == colnames(tad_coexprMatrix_cond1))
      ### CORRECTED HERE -> coexpression with itself should be 1 not zero !
      diag(tad_coexprMatrix_cond1) <- 1
      
      tad_genes_cmbs_cond1 <- combn(x = tad_genes, m = 2)
      
      all_ang_dist_cond1 <- apply(tad_genes_cmbs_cond1, 2, function(x) {
        stopifnot(x[1] %in% rownames(tad_coexprMatrix_cond1))
        stopifnot(x[2] %in% rownames(tad_coexprMatrix_cond1))
        get_ang_dist(vect1 = tad_coexprMatrix_cond1[x[1],], vect2 = tad_coexprMatrix_cond1[x[2],])
      })
      
      tad_meanAngDist_cond1 <- mean(all_ang_dist_cond1)
      
      ### WITH ALL GENES
      all_coexprDT_cond1 <- coexprDT_cond1[coexprDT_cond1$gene1 %in% tad_genes |
                                             coexprDT_cond1$gene2 %in% tad_genes,]
      
      stopifnot(nrow(all_coexprDT_cond1) > 0)  
      stopifnot(!is.na(all_coexprDT_cond1$coexpr))
      stopifnot(is.numeric(all_coexprDT_cond1$coexpr))
      stopifnot(tad_genes %in% all_coexprDT_cond1$gene1 | tad_genes %in% all_coexprDT_cond1$gene2 )
      
      all_coexprMatrix_cond1 <- as.matrix(get.adjacency(graph.data.frame(all_coexprDT_cond1, directed=FALSE), 
                                                  attr=colnames(all_coexprDT_cond1)[3], sparse=FALSE))
      stopifnot(tad_genes %in% rownames(all_coexprMatrix_cond1) & tad_genes %in% colnames(all_coexprMatrix_cond1))
      stopifnot( rownames(all_coexprMatrix_cond1) == colnames(all_coexprMatrix_cond1))
      ### CORRECTED HERE -> coexpression with itself should be 1 not zero !
      diag(all_coexprMatrix_cond1) <- 1
      
      tad_genes_cmbs <- combn(x = tad_genes, m = 2)
      
      all_ang_dist_cond1 <- apply(tad_genes_cmbs, 2, function(x) {
        stopifnot(x[1] %in% rownames(all_coexprMatrix_cond1))
        stopifnot(x[2] %in% rownames(all_coexprMatrix_cond1))
        get_ang_dist(vect1 = all_coexprMatrix_cond1[x[1],], vect2 = all_coexprMatrix_cond1[x[2],])
      })
      
      all_meanAngDist_cond1 <- mean(all_ang_dist_cond1)
      
      ####################################################################
      #### COND 2
      ####################################################################
      
      tad_coexprDT_cond2 <- coexprDT_cond2[coexprDT_cond2$gene1 %in% tad_genes & 
                                             coexprDT_cond2$gene2 %in% tad_genes,]
      
      stopifnot(nrow(tad_coexprDT_cond2) > 0)  
      stopifnot(!is.na(tad_coexprDT_cond2$coexpr))
      stopifnot(is.numeric(tad_coexprDT_cond2$coexpr))
      stopifnot(tad_genes %in% tad_coexprDT_cond2$gene1 | tad_genes %in% tad_coexprDT_cond2$gene2 )
      
      
      tad_coexprMatrix_cond2 <- as.matrix(get.adjacency(graph.data.frame(tad_coexprDT_cond2, directed=FALSE), 
                                                        attr=colnames(tad_coexprDT_cond2)[3], sparse=FALSE))
      stopifnot(tad_genes %in% rownames(tad_coexprMatrix_cond2) & tad_genes %in% colnames(tad_coexprMatrix_cond2))
      stopifnot(rownames(tad_coexprMatrix_cond2) %in% tad_genes &  colnames(tad_coexprMatrix_cond2) %in% tad_genes)
      stopifnot( rownames(tad_coexprMatrix_cond2) == colnames(tad_coexprMatrix_cond2))
      ### CORRECTED HERE -> coexpression with itself should be 1 not zero !
      diag(tad_coexprMatrix_cond2) <- 1
      
      tad_genes_cmbs_cond2 <- combn(x = tad_genes, m = 2)
      
      all_ang_dist_cond2 <- apply(tad_genes_cmbs_cond2, 2, function(x) {
        stopifnot(x[1] %in% rownames(tad_coexprMatrix_cond2))
        stopifnot(x[2] %in% rownames(tad_coexprMatrix_cond2))
        get_ang_dist(vect1 = tad_coexprMatrix_cond2[x[1],], vect2 = tad_coexprMatrix_cond2[x[2],])
      })
      
      tad_meanAngDist_cond2 <- mean(all_ang_dist_cond2)
      
      
      ### WITH ALL GENES
      all_coexprDT_cond2 <- coexprDT_cond2[coexprDT_cond2$gene1 %in% tad_genes |
                                             coexprDT_cond2$gene2 %in% tad_genes,]
      
      stopifnot(nrow(all_coexprDT_cond2) > 0)  
      stopifnot(!is.na(all_coexprDT_cond2$coexpr))
      stopifnot(is.numeric(all_coexprDT_cond2$coexpr))
      stopifnot(tad_genes %in% all_coexprDT_cond2$gene1 | tad_genes %in% all_coexprDT_cond2$gene2 )
      
      all_coexprMatrix_cond2 <- as.matrix(get.adjacency(graph.data.frame(all_coexprDT_cond2, directed=FALSE), 
                                                        attr=colnames(all_coexprDT_cond2)[3], sparse=FALSE))
      stopifnot(tad_genes %in% rownames(all_coexprMatrix_cond2) & tad_genes %in% colnames(all_coexprMatrix_cond2))
      stopifnot( rownames(all_coexprMatrix_cond2) == colnames(all_coexprMatrix_cond2))
      ### CORRECTED HERE -> coexpression with itself should be 1 not zero !
      diag(all_coexprMatrix_cond2) <- 1
      
      
      tad_genes_cmbs <- combn(x = tad_genes, m = 2)
      
      all_ang_dist_cond2 <- apply(tad_genes_cmbs, 2, function(x) {
        stopifnot(x[1] %in% rownames(all_coexprMatrix_cond2))
        stopifnot(x[2] %in% rownames(all_coexprMatrix_cond2))
        get_ang_dist(vect1 = all_coexprMatrix_cond2[x[1],], vect2 = all_coexprMatrix_cond2[x[2],])
      })
      
      all_meanAngDist_cond2 <- mean(all_ang_dist_cond2)
      
      ################################################################################
      
      
      list(
        tad_meanAngDist = tad_meanAngDist,
        tad_meanAngDist_cond1 = tad_meanAngDist_cond1,
        tad_meanAngDist_cond2 = tad_meanAngDist_cond2,
        all_meanAngDist = all_meanAngDist,
        all_meanAngDist_cond1 = all_meanAngDist_cond1,
        all_meanAngDist_cond2 = all_meanAngDist_cond2
      )
    } # end-foreach iterating over TADs
    names(ang_dist_data) <- pipeline_regionList
    ang_dist_data
  } # end-foreach iterating over dataset
  names(allData_ang_dist) <- coexprFilesOutNames
  
  outFile <- file.path(outFolder, "allData_ang_dist.Rdata")
  save(allData_ang_dist, file = outFile)
  cat(paste0("... written: ", outFile, "\n"))
  
} else { # end if-buildData
  outFile <- file.path(outFolder, "allData_ang_dist.Rdata")
  allData_ang_dist <- eval(parse(text = load(outFile)))
}




# sort the TADs by decreasing withinCoexpr
# plot level of coexpr within and between on the same plot

tad_angDist_DT <- data.frame(
  dataset = as.character(unlist(lapply(1:length(allData_ang_dist), function(i) {
    ds_name <- names(allData_ang_dist)[i]
    ds_name <- gsub("^CREATE_COEXPR_SORTNODUP/", "", ds_name)
    ds_name <- gsub("/pearson/coexprDT.Rdata$", "", ds_name)
    rep(ds_name, length(allData_ang_dist[[i]]))
  }))),
  region = as.character(unlist(lapply(1:length(allData_ang_dist), function(i) {
    names(allData_ang_dist[[i]])
  }))),
  
  tad_meanAngDist = as.numeric(unlist(lapply(allData_ang_dist, 
                                          function(sublist) lapply(sublist, function(x) x[["tad_meanAngDist"]])))),
  
  
  tad_meanAngDist_cond1 = as.numeric(unlist(lapply(allData_ang_dist, 
                                             function(sublist) lapply(sublist, function(x) x[["tad_meanAngDist_cond1"]])))),
  
  tad_meanAngDist_cond2 = as.numeric(unlist(lapply(allData_ang_dist, 
                                             function(sublist) lapply(sublist, function(x) x[["tad_meanAngDist_cond2"]])))),
  
  all_meanAngDist = as.numeric(unlist(lapply(allData_ang_dist, 
                                             function(sublist) lapply(sublist, function(x) x[["all_meanAngDist"]])))),
  
  all_meanAngDist_cond1 = as.numeric(unlist(lapply(allData_ang_dist, 
                                             function(sublist) lapply(sublist, function(x) x[["all_meanAngDist_cond1"]])))),
  
  all_meanAngDist_cond2 = as.numeric(unlist(lapply(allData_ang_dist, 
                                             function(sublist) lapply(sublist, function(x) x[["all_meanAngDist_cond2"]])))),
  

  
  stringsAsFactors = FALSE
)

tad_angDist_DT$tadAngDistDiffCond1Cond2 <- (tad_angDist_DT$tad_meanAngDist_cond1 - tad_angDist_DT$tad_meanAngDist_cond2)

tad_angDist_DT$tadAngDistRatioCond1Cond2 <- (tad_angDist_DT$tad_meanAngDist_cond1 / tad_angDist_DT$tad_meanAngDist_cond2)

tad_angDist_DT$allAngDistDiffCond1Cond2 <- (tad_angDist_DT$all_meanAngDist_cond1 - tad_angDist_DT$all_meanAngDist_cond2)

tad_angDist_DT$allAngDistRatioCond1Cond2 <- (tad_angDist_DT$all_meanAngDist_cond1 / tad_angDist_DT$all_meanAngDist_cond2)

tad_angDist_DT$tadAngDistChangeratioCond1Cond2 <- (tad_angDist_DT$tad_meanAngDist_cond2 - tad_angDist_DT$tad_meanAngDist_cond1)/tad_angDist_DT$tad_meanAngDist_cond1
tad_angDist_DT$allAngDistChangeratioCond1Cond2 <- (tad_angDist_DT$all_meanAngDist_cond2 - tad_angDist_DT$all_meanAngDist_cond1)/tad_angDist_DT$all_meanAngDist_cond1
cat("helloC\n")

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

tad_angDist_fc_DT <- merge(tad_angDist_DT, fc_DT, by=c("dataset", "region"))

cat("helloD\n")
##### build the coexpr table

coexprFolder <- "COEXPR_BETWEEN_WITHIN_ALL"
coexprFile <- file.path(coexprFolder, "allData_within_between_coexpr.Rdata")
stopifnot(file.exists(coexprFile))
allData_within_between_coexpr <- eval(parse(text = load(coexprFile)))

cat("helloE\n")

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

tad_angDist_fc_coexpr_DT <- merge(tad_angDist_fc_DT, tad_coexpr_DT, by=c("dataset", "region"))

cat("helloF\n")

tad_angDist_fc_coexpr_DT$withinDiffCond1Cond2 <- (tad_angDist_fc_coexpr_DT$withinCoexpr_cond1 - tad_angDist_fc_coexpr_DT$withinCoexpr_cond2)
tad_angDist_fc_coexpr_DT$withinRatioCond1Cond2 <- (tad_angDist_fc_coexpr_DT$withinCoexpr_cond1 / tad_angDist_fc_coexpr_DT$withinCoexpr_cond2)
tad_angDist_fc_coexpr_DT$withinChangeratioCond1Cond2 <- (tad_angDist_fc_coexpr_DT$withinCoexpr_cond2 - tad_angDist_fc_coexpr_DT$withinCoexpr_cond1)/tad_angDist_fc_coexpr_DT$withinCoexpr_cond1


tad_angDist_fc_coexpr_DT$withinBetweenDiffAll <- (tad_angDist_fc_coexpr_DT$withinCoexpr - tad_angDist_fc_coexpr_DT$betweenAllCoexpr) 
tad_angDist_fc_coexpr_DT$withinBetweenRatioAll <- (tad_angDist_fc_coexpr_DT$withinCoexpr / tad_angDist_fc_coexpr_DT$betweenAllCoexpr) 

tad_angDist_fc_coexpr_DT$withinBetweenDiffNbr <- (tad_angDist_fc_coexpr_DT$withinCoexpr - tad_angDist_fc_coexpr_DT$betweenNbrCoexpr) 
tad_angDist_fc_coexpr_DT$withinBetweenRatioNbr <- (tad_angDist_fc_coexpr_DT$withinCoexpr / tad_angDist_fc_coexpr_DT$betweenNbrCoexpr) 

tad_angDist_fc_coexpr_DT$withinBetweenDiffKb <- (tad_angDist_fc_coexpr_DT$withinCoexpr - tad_angDist_fc_coexpr_DT$betweenKbCoexpr) 
tad_angDist_fc_coexpr_DT$withinBetweenRatioKb <- (tad_angDist_fc_coexpr_DT$withinCoexpr / tad_angDist_fc_coexpr_DT$betweenKbCoexpr) 



outFile <- file.path(outFolder, paste0("multidens_tad_meanAngDist.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myHeight*1.2))
plot_multiDens(list(
  tad_meanAngDist = tad_angDist_fc_coexpr_DT[,"tad_meanAngDist"],
  tad_meanAngDist_cond1 = tad_angDist_fc_coexpr_DT[,"tad_meanAngDist_cond1"],
  tad_meanAngDist_cond2 = tad_angDist_fc_coexpr_DT[,"tad_meanAngDist_cond2"]
), my_xlab = "TAD mean ang. dist.",
legPos="topleft"
)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


outFile <- file.path(outFolder, paste0("multidens_all_meanAngDist.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myHeight*1.2))
plot_multiDens(list(
  all_meanAngDist = tad_angDist_fc_coexpr_DT[,"all_meanAngDist"],
  all_meanAngDist_cond1 = tad_angDist_fc_coexpr_DT[,"all_meanAngDist_cond1"],
  all_meanAngDist_cond2 = tad_angDist_fc_coexpr_DT[,"all_meanAngDist_cond2"]
), my_xlab = "all mean ang. dist.",
legPos="topleft"
)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))




outFile <- file.path(outFolder, paste0("multidens_all_tad_meanAngDist.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myHeight*1.2))
plot_multiDens(list(
  all_meanAngDist = tad_angDist_fc_coexpr_DT[,"all_meanAngDist"],
  tad_meanAngDist = tad_angDist_fc_coexpr_DT[,"tad_meanAngDist"]
), my_xlab = "all+tad mean ang. dist.",
legPos="topleft"
)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

outFile <- file.path(outFolder, paste0("multidens_all_tad_meanAngDist_cond1.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myHeight*1.2))
plot_multiDens(list(
  all_meanAngDist_cond1 = tad_angDist_fc_coexpr_DT[,"all_meanAngDist_cond1"],
  tad_meanAngDist_cond1 = tad_angDist_fc_coexpr_DT[,"tad_meanAngDist_cond1"]
), my_xlab = "all+tad mean ang. dist. (cond1)",
legPos="topleft"
)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


outFile <- file.path(outFolder, paste0("multidens_all_tad_meanAngDist_cond2.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myHeight*1.2))
plot_multiDens(list(
  all_meanAngDist_cond2 = tad_angDist_fc_coexpr_DT[,"all_meanAngDist_cond2"],
  tad_meanAngDist_cond2 = tad_angDist_fc_coexpr_DT[,"tad_meanAngDist_cond2"]
), my_xlab = "all+tad mean ang. dist. (cond2)",
legPos="topleft"
)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))




outFile <- file.path(outFolder, paste0("multidens_all_meanAngDist.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myHeight*1.2))
plot_multiDens(list(
  all_meanAngDist = tad_angDist_fc_coexpr_DT[,"all_meanAngDist"],
  tad_meanAngDist = tad_angDist_fc_coexpr_DT[,"tad_meanAngDist"],
  all_meanAngDist_cond1 = tad_angDist_fc_coexpr_DT[,"all_meanAngDist_cond1"],
  tad_meanAngDist_cond1 = tad_angDist_fc_coexpr_DT[,"tad_meanAngDist_cond1"],
  all_meanAngDist_cond2 = tad_angDist_fc_coexpr_DT[,"all_meanAngDist_cond2"],
  tad_meanAngDist_cond2 = tad_angDist_fc_coexpr_DT[,"tad_meanAngDist_cond2"]
), my_xlab = "all mean ang. dist.",
legPos="topleft"
)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))




cat("\nhelloA\n")
colnames(tad_angDist_fc_coexpr_DT)
cat("\nhelloA\n")

all_y <- c("tad_meanAngDist", "tad_meanAngDist_cond1", "tad_meanAngDist_cond2", 
           "all_meanAngDist", "all_meanAngDist_cond1", "all_meanAngDist_cond2")
xvar <- "withinCoexpr"
myx <- tad_angDist_fc_coexpr_DT[, xvar] 

for(yvar in all_y) {
    myy <- tad_angDist_fc_coexpr_DT[, yvar] 
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
    addCorr(x=myx,y=myy,legPos="topright", bty='n')
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
  
# ######################################################################################
# ######################################################################################
# ######################################################################################
cat(paste0("*** DONE: ", script_name, "\n"))
cat(paste0(startTime, "\n", Sys.time(), "\n"))





