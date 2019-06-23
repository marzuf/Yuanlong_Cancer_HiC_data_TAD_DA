
# Rscript extract_cptmtScore.R

script_name <- "extract_cptmtScore.R"
startTime <- Sys.time()

options(scipen=100)

SSHFS <- FALSE
prefixDir <- ifelse(SSHFS, "/media/electron", "")
allData_file <- file.path(prefixDir, "/mnt/ndata/Yuanlong/1.Projects/2.PROFILE/0.Scripts/Generate_TADs_for_Marie/bin_comp_CELL_LINEs_list.Rdata")
stopifnot(file.exists(allData_file))

nCpu <- ifelse(SSHFS, 2, 40)
suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
registerDoMC(nCpu)

binSize <- 40000
binSizeKb <- 40

allData <- eval(parse(text = load(allData_file)))

cancerDS <- names(allData)[grep("^mega_", names(allData))]

ds = cancerDS[1]
# cancerDS=cancerDS[1]
foo <- foreach(ds = cancerDS) %dopar% {
  
  cat(paste0("... start DS: ", ds, "\n"))
  
  folderName <- gsub("^mega_", "", ds) 
  outFolder <- file.path(paste0(folderName, "_",binSizeKb , "kb"), "CPTMT_SCORE")
  dir.create(outFolder, recursive = TRUE)
  
  allChr_tad_DT <- allData[[ds]]
  stopifnot( colnames(allChr_tad_DT) == c("chr", "bin_index", "CD_rank"))
  allChr_tad_DT$chr <- paste0("chr", as.character(allChr_tad_DT$chr))
  
  all_chrs <- unique(allChr_tad_DT$chr)
  
  chromo="chr1"
  all_chromoDT <- foreach(chromo = all_chrs, .combine='rbind') %do% {
    
    cat(paste0("...... start chromo: ", chromo, "\n"))
    
    chr_DT <- allChr_tad_DT[allChr_tad_DT$chr == chromo, ,drop=FALSE]
    chr_DT$CD_score <- as.numeric(chr_DT$CD_rank)
    stopifnot(!is.na(chr_DT$CD_score))
    chr_DT$CD_rank <- factor(chr_DT$CD_rank, levels = unique(chr_DT$CD_rank))
    chr_DT$domainNbr <- as.numeric(chr_DT$CD_rank)
    stopifnot( diff(chr_DT$domainNbr) == 0 | diff(chr_DT$domainNbr) == 1 )
    
    all_bins_tadPos_DT <- chr_DT
    all_bins_tadPos_DT$CD_rank <- NULL
    head( all_bins_tadPos_DT)
    
    stopifnot( diff(as.numeric(all_bins_tadPos_DT$bin_index)) > 0 )
    
    tadPos_DT <- do.call(rbind, by(all_bins_tadPos_DT, all_bins_tadPos_DT$domainNbr, function(x) {
      regionIdx <- unique(x$domainNbr)
      cdScore <- unique(x$CD_score)
      stopifnot(length(cdScore) == 1)
      stopifnot(!is.na(cdScore))
      stopifnot(length(regionIdx) == 1)
      data.frame(
        chr = chromo,
        region = paste0(chromo, "_TAD", regionIdx),
        startIdx = min(as.numeric(as.character(x$bin_index))),
        endIdx = max(as.numeric(as.character(x$bin_index))),
        score = cdScore,
        stringsAsFactors = FALSE
      )
    }))
    stopifnot(nrow(tadPos_DT) == length(unique(all_bins_tadPos_DT$domainNbr)))
    tadPos_DT$startIdx <- as.numeric(as.character(tadPos_DT$startIdx))
    stopifnot(!is.na(tadPos_DT$startIdx))
    tadPos_DT$endIdx <- as.numeric(as.character(tadPos_DT$endIdx))
    stopifnot(!is.na(tadPos_DT$endIdx))
    stopifnot(tadPos_DT$endIdx >= tadPos_DT$startIdx)
    # pos_start = (bin_index - 1)*40E3 + 1; pos_end = bin_index*40E3
    # e.g. bin1 -> 1 - 40'000
    tadPos_DT$start <- (tadPos_DT$startIdx - 1)*binSize + 1 
    tadPos_DT$end <- (tadPos_DT$endIdx)*binSize 
    outDT <- tadPos_DT[,c("chr", "start", "end", "region", "score"), drop=FALSE]
    stopifnot(is.numeric(outDT$start))
    stopifnot(is.numeric(outDT$end))
    stopifnot(outDT$chr == chromo)
    outDT <- outDT[order(outDT$start, outDT$end),,drop=FALSE]
    stopifnot(outDT$end > outDT$start)
    
    stopifnot(outDT$end[-nrow(outDT)] < outDT$start[-1])
    
    # output file name should like GSE105194_cerebellum_40kb/FINAL_DOMAINS/GSE105194_cerebellum_chr3_KR_40kb_final_domains.txt
    # outFile <- file.path(outFolder, paste0(folderName, "_", chromo, "_", "YL", "_", binSizeKb, "kb_final_domains_withScore.txt"))
    # write.table(outDT, file = outFile, col.names = FALSE, row.names = FALSE, append = FALSE, sep="\t", quote = FALSE)
    # cat(paste0("... written: ", outFile, "\n"))
    outDT
  } # end foreach-iterating over chromo    
  
  outFile <- file.path(outFolder, paste0(folderName, "_", "all_chromo", "_", "YL", "_", binSizeKb, "kb_final_domains_withScore.txt"))
  write.table(all_chromoDT, file = outFile, col.names = FALSE, row.names = FALSE, append = FALSE, sep="\t", quote = FALSE)
  cat(paste0("... written: ", outFile, "\n"))
  
} # end for-iterating over ds



# check -> cmp with the assigned region

all_domainScore_files <- list.files(".", recursive = TRUE, pattern="_final_domains_withScore.txt", full.names = FALSE)
stopifnot(length(all_domainScore_files) > 0)

score_file = all_domainScore_files[1]
score_DT <- foreach(score_file = all_domainScore_files, .combine = 'rbind') %dopar% {
  curr_file <- file.path(score_file)
  stopifnot(file.exists(curr_file))
  curr_DT <- read.delim(curr_file, header=F, 
                        col.names = c("chromo", "start", "end", "region", "score"))
  curr_DT$dataset <- dirname(dirname(curr_file))
  curr_DT
}

all_assigned_files <- list.files(".", recursive = TRUE, pattern="all_assigned_regions.txt", full.names = FALSE)
stopifnot(length(all_assigned_files) > 0)

assigned_DT <- foreach(assigned_file = all_assigned_files, .combine = 'rbind') %dopar% {
  curr_file <- file.path(assigned_file)
  stopifnot(file.exists(curr_file))
  curr_DT <- read.delim(curr_file, header=F, 
                        col.names = c("chromo", "region", "start", "end"))
  curr_DT$dataset <- dirname(dirname(curr_file))
  curr_DT
}
assigned_DT <- assigned_DT[grep("_TAD", assigned_DT$region),]

colnames(score_DT)[colnames(score_DT) == "start"] <- "start_scores"
colnames(score_DT)[colnames(score_DT) == "end"] <- "end_scores"

colnames(assigned_DT)[colnames(assigned_DT) == "start"] <- "start_assigned"
colnames(assigned_DT)[colnames(assigned_DT) == "end"] <- "end_assigned"

stopifnot(nrow(score_DT) == nrow(assigned_DT))

all_DT <- merge(score_DT, assigned_DT, by = c("chromo", "region", "dataset"))

stopifnot(!is.na(all_DT))

stopifnot(nrow(score_DT) == nrow(all_DT))

stopifnot(all_DT$start_scores == all_DT$start_assigned)
stopifnot(all_DT$end_scores == all_DT$end_assigned)

########################################################################################
########################################################################################
########################################################################################
txt <- paste0(startTime, "\n", Sys.time(), "\n")
cat(paste0(txt))
cat(paste0("*** DONE: ", script_name, "\n"))


