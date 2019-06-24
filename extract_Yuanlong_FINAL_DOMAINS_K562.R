
# Rscript extract_Yuanlong_FINAL_DOMAINS_K562.R

script_name <- "extract_Yuanlong_FINAL_DOMAINS_K562.R"
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

# cancerDS <- names(allData)[grep("^mega_", names(allData))]
cancerDS <- "AWS_K562"

ds = cancerDS[1]
# cancerDS=cancerDS[1]


outFolderPrefix <- "K562"

foo <- foreach(ds = cancerDS) %dopar% {
  
  cat(paste0("... start DS: ", ds, "\n"))
  
  # folderName <- gsub("^mega_", "", ds) 
  # outFolder <- file.path(paste0(folderName, "_",binSizeKb , "kb"), "FINAL_DOMAINS")
  # dir.create(outFolder, recursive = TRUE)
  folderName <- outFolderPrefix
  outFolder <- file.path(paste0(outFolderPrefix, "_",binSizeKb , "kb"), "FINAL_DOMAINS")
  dir.create(outFolder, recursive = TRUE)
  
  allChr_tad_DT <- allData[[ds]]
  stopifnot( colnames(allChr_tad_DT) == c("chr", "bin_index", "CD_rank"))
  allChr_tad_DT$chr <- paste0("chr", as.character(allChr_tad_DT$chr))
  
  all_chrs <- unique(allChr_tad_DT$chr)
  
  chromo="chr1"
  for(chromo in all_chrs) {
    
    cat(paste0("...... start chromo: ", chromo, "\n"))
    
    chr_DT <- allChr_tad_DT[allChr_tad_DT$chr == chromo, ,drop=FALSE]
    chr_DT$CD_rank <- factor(chr_DT$CD_rank, levels = unique(chr_DT$CD_rank))
    chr_DT$domainNbr <- as.numeric(chr_DT$CD_rank)
    stopifnot( diff(chr_DT$domainNbr) == 0 | diff(chr_DT$domainNbr) == 1 )
    
    all_bins_tadPos_DT <- chr_DT
    all_bins_tadPos_DT$CD_rank <- NULL
    head( all_bins_tadPos_DT)
    
    stopifnot( diff(as.numeric(all_bins_tadPos_DT$bin_index)) > 0 )
    
    tadPos_DT <- do.call(rbind, by(all_bins_tadPos_DT, all_bins_tadPos_DT$domainNbr, function(x) {
      regionIdx <- unique(x$domainNbr)
      stopifnot(length(regionIdx) == 1)
      data.frame(
        chr = chromo,
        region = paste0(chromo, "_TAD", regionIdx),
        startIdx = min(as.numeric(as.character(x$bin_index))),
        endIdx = max(as.numeric(as.character(x$bin_index))),
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
    outDT <- tadPos_DT[,c("chr", "start", "end"), drop=FALSE]
    stopifnot(is.numeric(outDT$start))
    stopifnot(is.numeric(outDT$end))
    stopifnot(outDT$chr == chromo)
    outDT <- outDT[order(outDT$start, outDT$end),,drop=FALSE]
    stopifnot(outDT$end > outDT$start)
    
    stopifnot(outDT$end[-nrow(outDT)] < outDT$start[-1])
    
    # output file name should like GSE105194_cerebellum_40kb/FINAL_DOMAINS/GSE105194_cerebellum_chr3_KR_40kb_final_domains.txt
    outFile <- file.path(outFolder, paste0(folderName, "_", chromo, "_", "YL", "_", binSizeKb, "kb_final_domains.txt"))
    write.table(outDT, file = outFile, col.names = FALSE, row.names = FALSE, append = FALSE, sep="\t", quote = FALSE)
    cat(paste0("... written: ", outFile, "\n"))

  } # end for-iterating over chromo    
} # end for-iterating over ds

txt <- paste0(startTime, "\n", Sys.time(), "\n")
cat(paste0(txt))
cat(paste0("*** DONE: ", script_name, "\n"))


