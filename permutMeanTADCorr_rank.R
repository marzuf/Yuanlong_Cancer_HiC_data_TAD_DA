

stop("--- use directly empPval_meanTADCorr_rank.R \n")

# Rscript permutMeanTADCorr_rank.R

# Rscript permutMeanTADCorr_rank.R Panc1_rep12_40kb TCGApaad_wt_mutKRAS

withDiago <- FALSE

script_name <- "permutMeanTADCorr_rank.R"

cat(paste0("> START ", script_name, "\n"))

startTime <- Sys.time()

setDir <- ""

outFold <- "PERMUTMEANTADCORR_RANK"
dir.create(outFold)

require(foreach)
require(doMC)
registerDoMC(40)

script0_name <- "0_prepGeneData"
script7_name <- "7_runPermutationsMeanTADCorr"

# rank 1 -> high correlation

tiesMeth <- "min"


###################
pipOutFolder <- file.path("PIPELINE", "OUTPUT_FOLDER")

all_permut_files <- list.files(pipOutFolder, recursive = TRUE, pattern="meanCorr_permDT.Rdata", full.names = TRUE)
stopifnot(length(all_permut_files) > 0)

permut_file = all_permut_files[1]

args <- commandArgs(trailingOnly = TRUE)

if(length(args) == 2) {
  all_permut_files <- all_permut_files[grepl(args[1], all_permut_files) & grepl(args[2], all_permut_files)]
  stopifnot(length(all_permut_files) == 1)
}

for(permut_file in all_permut_files) {
  
  hicds <- basename(dirname(dirname(dirname(permut_file))))
  exprds <- basename(dirname(dirname(permut_file)))

  cat("... ", hicds, " - ", exprds, "\n")
  cat("...... load permDT", "\n")  
  meanCorr_permDT <- eval(parse(text = load(permut_file)))
  
  
  # take -x to have decreasing ranking
  rank_meanCorr_permDT <- apply(meanCorr_permDT, 2 , function(x) rank(-x, ties=tiesMeth))
  
  
  max1_idx <- which(rownames(meanCorr_permDT) == names(which(rank_meanCorr_permDT[,1] == 1)))
  stopifnot(meanCorr_permDT[-max1_idx,1] <  rank_meanCorr_permDT[max1_idx,1])
  
  outFile <- file.path(outFold, hicds, exprds, "rank_meanCorr_permDT.Rdata")
  dir.create(dirname(outFile), recursive = TRUE)
  save(rank_meanCorr_permDT, file = outFile)
  cat(paste0("... written: ", outFile, "\n"))
  
}

################################****************************************************************************************
################################****************************************************************************************

cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))








