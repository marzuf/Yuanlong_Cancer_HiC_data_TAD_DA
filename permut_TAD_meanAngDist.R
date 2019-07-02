# Rscript permut_TAD_meanAngDist.R 

script_name <- "permut_TAD_meanAngDist.R"

startTime <- Sys.time()

cat("> START permut_TAD_meanAngDist.R \n")

SSHFS <- FALSE

require(foreach)
require(doMC)

registerDoMC(40)

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

outFold <- "PERMUT_TAD_MEANANGDIST"
dir.create(outFold)

angDistFolder <- "CREATE_ANGDIST_BYCOND_SORTNODUP"
stopifnot(dir.exists(angDistFolder))

all_angDist_files <- list.files(angDistFolder, pattern = "all_angDist_DT.Rdata", recursive = TRUE, full.names = TRUE)
stopifnot(length(all_angDist_files) > 0)

pipOutfolder <- file.path("PIPELINE", "OUTPUT_FOLDER")
stopifnot(dir.exists(pipOutfolder))

script0_name <- "0_prepGeneData"
script5_name <- "5_runPermutationsMedian"

angDist_file = all_angDist_files[1]


all_angDist_files[1] = all_angDist_files[1]

for(angDist_file in all_angDist_files){
  
  hicds <- basename(dirname(dirname(dirname(angDist_file))))
  exprds <- basename(dirname(dirname(angDist_file)))
  
  
  regionFile <- file.path(pipOutfolder, hicds, exprds, script0_name, "pipeline_regionList.Rdata")
  stopifnot(file.exists(regionFile))
  regionList <- eval(parse(text = load(regionFile)))
  
  
  geneFile <- file.path(pipOutfolder, hicds, exprds, script0_name, "pipeline_geneList.Rdata")
  stopifnot(file.exists(geneFile))
  geneList <- eval(parse(text = load(geneFile)))
  
  g2tFile <- file.path(hicds, "genes2tad", "all_genes_positions.txt")
  stopifnot(file.exists(g2tFile))
  g2t_DT <- read.delim(g2tFile, header=F, col.names = c("entrezID",  "chromo", "start", "end", "region"), stringsAsFactors = FALSE)
  g2t_DT$entrezID <- as.character(g2t_DT$entrezID)
  
  stopifnot(regionList %in% g2t_DT$region)
  stopifnot(geneList %in% g2t_DT$entrezID)
  
  g2t_DT <- g2t_DT[g2t_DT$entrezID %in% geneList,]
  
  angDist_DT <- eval(parse(text = load(angDist_file)))
  
  permDT_file <- file.path(pipOutfolder, hicds, exprds, script5_name, "permutationsDT.Rdata")
  permutationsDT <- eval(parse(text = load(permDT_file)))
  
  all_regions <- sort(unique(as.character(permutationsDT[,2])))
  
  
  stopifnot(is.character(angDist_DT$gene1))
  stopifnot(is.character(angDist_DT$gene2))
  
  stopifnot(geneList %in% angDist_DT$gene1 | geneList %in% angDist_DT$gene2)
  
  stopifnot(setequal(unique(g2t_DT$region), regionList))
  stopifnot(setequal(all_regions, regionList))
  

  ####################################################### COMPUTE MEAN INTRA CORR BY TAD FOR THE PERMUTATIONS
  
  cat("... start intraCorr permutDT \n")
  
  permutationsDT <- permutationsDT[,1:11]
  
  tadMeanAngDist_permDT_allReg <- foreach(i_col = 1:ncol(permutationsDT), .combine='cbind') %dopar% {
    
    cat(paste0("... angDist for permutation: ", i_col, "/", ncol(permutationsDT), "\n"))
    
    g2t_permDT <- data.frame(entrezID = rownames(permutationsDT), 
                             region = permutationsDT[,i_col], stringsAsFactors = F)
    g2t_permDT$entrezID <- as.character(g2t_permDT$entrezID)
    g2t_permDT$region <-  as.character(g2t_permDT$region)
    
    permutAngDist <- sapply(regionList, function(reg) {
      reg_genes <- g2t_permDT$entrezID[g2t_permDT$region == reg]
      
      
      tad_angDist_DT <- angDist_DT[angDist_DT$gene1 %in% reg_genes & angDist_DT$gene2 %in% reg_genes,]
      mean(tad_angDist_DT$angDist)
      
      
    })
    # permutAngDist <- rbindlist(permutAngDist)
    curr_permutDT <- data.frame(permutAngDist)
    stopifnot(ncol(curr_permutDT) == 1)
    colnames(curr_permutDT) <- paste0("result", i_col-1)
    stopifnot(all(rownames(curr_permutDT) == regionList))
    curr_permutDT
  }
  cat("... end angDist permutDT \n")
  
  meanAngDist_permDT <- as.data.frame(tadMeanAngDist_permDT_allReg)
  stopifnot(ncol(meanAngDist_permDT) == ncol(permutationsDT))  
  colnames(meanAngDist_permDT) <- paste0("permutation",  c(1:ncol(permutationsDT)))
  stopifnot(nrow(meanAngDist_permDT) == length(regionList))
  rownames(meanAngDist_permDT) <- regionList
  

  outFile <-  file.path(outFold, hicds, exprds, "meanAngDist_permDT.Rdata")
  dir.create(dirname(outFile), recursive = TRUE)
  save(meanAngDist_permDT, file= outFile)
  cat(paste0("... written: ", outFile, "\n"))
  

  
  
}


######################################################################################
######################################################################################
######################################################################################
cat("*** DONE: ", script_name, "\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))

































