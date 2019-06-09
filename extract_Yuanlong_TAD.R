SSHFS <- FALSE
prefixDir <- ifelse(SSHFS, "/media/electron", "")
allData_file <- file.path(prefixDir, "/mnt/ndata/Yuanlong/1.Projects/2.PROFILE/0.Scripts/Generate_TADs_for_Marie/bin_comp_CELL_LINEs_list.Rdata")
stopifnot(file.exists(allData_file))

binSize <- 40000

allData <- eval(parse(text = load(allData_file)))

cancerDS <- names(allData)[grep("^mega_", names(allData))]

ds = cancerDS[1]
  
foreach(ds = cancerDS) %dopar% {
  
  
  folderName <- gsub("^mega_", "", ds) 
  
  dir.create(folderName)
  
  allChr_tad_DT <- allData[[ds]]
  stopifnot( colnames(allChr_tad_DT) == c("chr", "bin_index", "CD_rank"))
  allChr_tad_DT$chr <- paste0("chr", as.character(allChr_tad_DT$chr))
  
  all_chrs <- unique(allChr_tad_DT$chr)
  
  chromo="chr1"
  for(chromo in all_chrs) {
  
    chr_DT <- allChr_tad_DT[allChr_tad_DT$chr == chromo, ,drop=FALSE]
    chr_DT$CD_rank <- factor(chr_DT$CD_rank, levels = unique(chr_DT$CD_rank))
    chr_DT$domainNbr <- as.numeric(chr_DT$CD_rank)
    stopifnot( diff(chr_DT$domainNbr) == 0 | diff(chr_DT$domainNbr) == 1 )
    
    all_bins_tadPos_DT <- chr_DT
    all_bins_tadPos_DT$CD_rank <- NULL
    head( all_bins_tadPos_DT)
    
    tadPos_DT <- do.call(rbind, by(all_bins_tadPos_DT, all_bins_tadPos_DT$domainNbr, function(x) {
      
      regionIdx <- unique(x$domainNbr)
      stopifnot(length(regionIdx) == 1)
      
      data.frame(
        chr = chromo,
        region = paste0(chromo, "_TAD", regionIdx),
        startIdx = min(x$bin_index), 
        endIdx = max(x$bin_index),
        stringsAsFactors = FALSE
        )
    
    }))
    
    
    stopifnot(nrow(tadPos_DT) == length(unique(all_bins_tadPos_DT$domainNbr)))
    tadPos_DT$startIdx <- as.numeric(as.character(tadPos_DT$startIdx))
    stopifnot(!is.na(tadPos_DT$startIdx))
    tadPos_DT$endIdx <- as.numeric(as.character(tadPos_DT$endIdx))
    stopifnot(!is.na(tadPos_DT$endIdx))
    # pos_start = (bin_index - 1)*40E3 + 1; pos_end = bin_index*40E3
    # e.g. bin1 -> 1 - 40'000
    tadPos_DT$start <- (tadPos_DT$startIdx - 1)*binSize + 1 
    tadPos_DT$end <- (tadPos_DT$endIdx)*binSize 
    
    # allStarts <- tadPos_DT$end + 1
    # allStarts <- allStarts[1:(length(allStarts) -1)]
    # allEnds <- tadPos_DT$start- 1
    # allEnds <- allEnds[-1]
    
    allStarts <- c(1,tadPos_DT$end + 1)
    allStarts <- allStarts[1:(length(allStarts) -1)]
    allEnds <- tadPos_DT$start- 1
    
    stopifnot( length(allStarts) == length(allEnds) )
    
    
    bdStarts <- allStarts[allStarts < allEnds]
    bdEnds <- allEnds[allStarts < allEnds]
    
    nBD <- length(bdStarts)
    stopifnot( nBD == length(bdEnds))
    
    if(nBD > 0) {
      
      
      
      stopifnot(bdEnds > bdStarts)
      
      bd_DT <- data.frame(
        chr = chromo,
        region = paste0(chromo, "_BOUND", 1:nBD),
        start = bdStarts,
        end = bdEnds,
        stringsAsFactors = FALSE
      )
      
      all_regions_DT <- rbind(
        tadPos_DT[,c("chr", "region", "start", "end")],
        bd_DT
      )
        
    
      
    } else {
      all_regions_DT <- tadPos_DT[,c("chr", "region", "start", "end")]
    }

    stopifnot(nrow(all_regions_DT) > 0)


    stopifnot(is.numeric(all_regions_DT$start))
    stopifnot(is.numeric(all_regions_DT$end))
    
    all_regions_DT <- all_regions_DT[order(all_regions_DT$start, all_regions_DT$end),]
    
    tmpStarts <- all_regions_DT$end + 1
    stopifnot(tmpStarts[-length(tmpStarts)] == all_regions_DT$start[-1])    
    stopifnot(all_regions_DT$start[1] == 1)
  }
  
  
}