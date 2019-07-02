# Rscript empPval_TAD_meanAngDist.R

script_name <- "empPval_TAD_meanAngDist.R"

cat("> Start : ", script_name, "\n")

startTime <- Sys.time()


# maximal value 1 
# if smaller than 1 -> more similar
# frequency that by chance the distance is smaller than the observed distance

suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
registerDoMC(40) 

outFold <- "EMPPVAL_TAD_MEANANGDIST"
dir.create(outFold)

obsAngDistFolder <- "TAD_MEANANGDIST"
stopifnot(dir.exists(obsAngDistFolder))


permAngDistFolder <- "PERMUT_TAD_MEANANGDIST"
stopifnot(dir.exists(permAngDistFolder))

all_obsDist_files <- list.files(obsAngDistFolder, pattern = "all_tad_meanAngDist.Rdata", recursive = TRUE, full.names = TRUE)
stopifnot(length(all_angDist_files) > 0)

obs_file = all_obsDist_files[1]

for(obs_file in all_obsDist_files) {
  
  hicds <- basename(dirname(dirname(obs_file)))
  exprds <- basename(dirname(obs_file))
  
  perm_file <- file.path(permAngDistFolder, hicds, exprds, "meanAngDist_permDT.Rdata")  
  stopifnot(file.exists(perm_file))
  
  all_obs_meanAngDist <- eval(parse(text = load(obs_file)))
  
  
  permDT <- eval(parse(text = load(perm_file)))
  all_regions <- rownames(permDT)
  
  regionListFile <- file.path(pipFolder, hicds, exprds, script0_name, "pipeline_geneList.Rdata")
  stopifnot(file.exists(regionListFile))
  regionList <- eval(parse(text = load(regionListFile)))
  
  stopifnot(setequal(all_regions, regionList))
  stopifnot(setequal(all_regions, names(all_obs_meanAngDist)))
  
  stopifnot(regionList %in% names(all_obs_meanAngDist))
  
  emp_pval_meanAngDist <- foreach(reg = all_regions, .combine='c') %dopar% {
    
    perm_angDist_reg <- permDT[paste0(reg),]
    
    obs_meanAngDist <- all_obs_meanAngDist[reg]
    
    stopifnot(length(obs_meanAngDist) == 1)
    
    # if there is a lot of permutations where the distance is smaller
    # than the observed -> pval will be high (not signif)
    emp_pval <- sum(perm_angDist_reg <= obs_meanAngDist)
    
    (emp_pval+1)/(length(perm_angDist_reg)+1)
    
    
  }
  names(emp_pval_meanAngDist) <- all_regions
  
  
  stopifnot(all(emp_pval_meanAngDist > 0 & emp_pval_meanAngDist <= 1 ))
  
  outFile <- file.path(outFold, hicds, exprds, "emp_pval_meanAngDist.Rdata")
  dir.create(dirname(outFile), recursive = TRUE)
  
  save(emp_pval_meanAngDist, file = outFile)
  cat(paste0("... written: ", outFile, "\n"))
  
  
}
  
  
##########################################################################
cat(paste0("*** DONE: ", script_name, "\n"))
cat(paste0(startTime, "\n", Sys.time(), "\n"))

