# Rscript empPval_meanTADCorr_rank.R

# Rscript empPval_meanTADCorr_rank.R Panc1_rep12_40kb TCGApaad_wt_mutKRAS

withDiago <- FALSE

script_name <- "empPval_meanTADCorr_rank.R"

cat(paste0("> START ", script_name, "\n"))

startTime <- Sys.time()

setDir <- ""

outFold <- "EMPPVAL_MEANTADCORR_RANK"
dir.create(outFold)

require(foreach)
require(doMC)
registerDoMC(40)

script0_name <- "0_prepGeneData"
script4_name <- "4_runMeanTADCorr"
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

  # observed data
  obs_meanCorr_TAD_file <- file.path(pipOutFolder, hicds, exprds, script4_name, "all_meanCorr_TAD.Rdata")
  stopifnot(file.exists(obs_meanCorr_TAD_file))
  meanCorr_obs <- eval(parse(text = load(obs_meanCorr_TAD_file)))

  # take -x to have decreasing ranking
  meanCorr_obs_rank <- rank(-meanCorr_obs, ties = tiesMeth)
  all_regions <- names(meanCorr_obs_rank)
  
  cat("...... load permDT", "\n")  
  meanCorr_permDT <- eval(parse(text = load(permut_file)))
  stopifnot( setequal(all_regions, rownames(meanCorr_permDT) ))
  
  # take -x to have decreasing ranking
  rank_meanCorr_permDT <- apply(meanCorr_permDT, 2 , function(x) rank(-x, ties=tiesMeth))
  max1_idx <- which(rownames(meanCorr_permDT) == names(which(rank_meanCorr_permDT[,1] == 1)))
  stopifnot(meanCorr_permDT[-max1_idx,1] <  rank_meanCorr_permDT[max1_idx,1])
  
  stopifnot( setequal(all_regions, rownames(rank_meanCorr_permDT) ))

  
  
  # compute the empirical p-value  
  
  emp_pval_meanCorr_rank <- unlist(foreach(reg = all_regions, .combine='c') %dopar% {
    # select the random values for this region
    all_shuff_meanCorr_rank <- rank_meanCorr_permDT[paste0(reg),]
    stopifnot(length(all_shuff_meanCorr_rank) == ncol(rank_meanCorr_permDT))
    # select the observed value for this region
    obs_rank <- meanCorr_obs_rank[reg]
    stopifnot(length(obs_rank) == 1 )
    # the number of times the permut rank is smaller than the observed one
    # => the pvalue will be high if the permutation often a smaller rank
    emp_pval <- sum(all_shuff_meanCorr_rank <= obs_rank) 
    # emp_pval/length(shuff_intraTADcorr_region)
    ### ADD THE +1 => THIS IS FOR THE OBSERVED VALUE -> AVOID THE 0 VALUES !!!
    (emp_pval+1)/(length(all_shuff_meanCorr_rank)+1)
  })
  names(emp_pval_meanCorr_rank) <- all_regions
  
  stopifnot(all(emp_pval_meanCorr_rank > 0 & emp_pval_meanCorr_rank <= 1 ))
  
  
  outFile <- file.path(outFold, hicds, exprds, "emp_pval_meanCorr_rank.Rdata")
  dir.create(dirname(outFile), recursive = TRUE)
  save(emp_pval_meanCorr_rank, file = outFile)
  cat(paste0("... written: ", outFile, "\n"))
  
}

################################****************************************************************************************
################################****************************************************************************************

cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))








