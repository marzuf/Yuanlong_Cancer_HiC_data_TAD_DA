# Rscript TAD_meanAngDist.R 

script_name <- "TAD_meanAngDist.R"

startTime <- Sys.time()

cat("> START TAD_meanAngDist.R \n")

SSHFS <- FALSE

require(foreach)
require(doMC)

registerDoMC(40)

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

outFold <- "TAD_MEANANGDIST"
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
  
  stopifnot(is.character(angDist_DT$gene1))
  stopifnot(is.character(angDist_DT$gene2))
  
  stopifnot(geneList %in% angDist_DT$gene1 | geneList %in% angDist_DT$gene2)
  
  stopifnot(setequal(unique(g2t_DT$region), regionList))
  
  all_tad_meanAngDist <- foreach(tad = regionList, .combine='c') %dopar% {
    tad_genes <- as.character(g2t_DT$entrezID[g2t_DT$region == tad])
    tad_angDist_DT <- angDist_DT[angDist_DT$gene1 %in% tad_genes & angDist_DT$gene2 %in% tad_genes,]
    mean(tad_angDist_DT$angDist)
  }
  names(all_tad_meanAngDist) <- regionList
  
  
  outFile <- file.path(outFold, hicds, exprds, "all_tad_meanAngDist.Rdata")
  dir.create(dirname(outFile), recursive = TRUE)
  save(all_tad_meanAngDist, file = outFile)
  cat(paste0("... written: ", outFile, "\n"))
}


######################################################################################
######################################################################################
######################################################################################
cat("*** DONE: ", script_name, "\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))


