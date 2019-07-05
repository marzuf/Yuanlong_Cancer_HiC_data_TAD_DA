
# ISSUE WITH THE WRONG DATA -> for a given dataset, all chromosomes were the same
# setDir <- "/media/electron"
# load(file.path(setDir, "/mnt/ndata/Yuanlong/1.Projects/2.PROFILE/0.Scripts/Generate_TADs_for_Marie/bin_comp_CELL_LINEs_list.Rdata"))
# 
# all_ds <- names(bin_comp_CELL_LINEs_list)
# all_ds <- all_ds[grepl("mega_", all_ds)]
# 
# for(ds in all_ds) {
#  cat("... ", ds, "\n") 
#   ref_chr <- 1
#   
#   currdt <- bin_comp_CELL_LINEs_list[[paste0(ds)]]
# 
#   ref_ranks <- currdt$CD_rank[currdt$chr == ref_chr]
#   ref_idxs <- currdt$bin_index[currdt$chr == ref_chr]
#     
#   other_chr <- 2:22
#   
#   for(curr_chr in other_chr) {
#     cat("....... ", ds, " - ", curr_chr, "\n")   
#     curr_ranks <- currdt$CD_rank[currdt$chr == curr_chr]
#     curr_idxs <- currdt$bin_index[currdt$chr == curr_chr]
#     
#     stopifnot(all(ref_ranks == curr_ranks))
#     stopifnot(all(ref_idxs == curr_idxs))
#     # break
#   }
#   # break
# }
# 
# currdt <- bin_comp_CELL_LINEs_list[["mega_ENCSR079VIJ_G401"]]
# 
# currdt <- bin_comp_CELL_LINEs_list[["mega_ENCSR549MGQ_T47D"]]
# 
# currdt_1 <- currdt[currdt$chr == 1,]
# stopifnot(nrow(currdt_1) > 0)
# 
# currdt_21 <- currdt[currdt$chr == 21,]
# stopifnot(nrow(currdt_21) > 0)
# 
# all(currdt_1$CD_rank == currdt_21$CD_rank)
# all(currdt_1$index == currdt_21$bin_index)




# Rscript check_YL_TADs.R

require(foreach)
require(doMC)
require(ggplot2)
require(hrbrthemes)
require(ggthemes)
startTime <- Sys.time()

script_name <- "check_YL_TADs.R"

cat("> START ", script_name, "\n")

outFold <- "CHECK_YL_TADS"
dir.create(outFold)

buildData <- FALSE

SSHFS <- FALSE
setDir <- ifelse(SSHFS, "/media/electron", "")

registerDoMC(ifelse(SSHFS, 2, 40))

binKb <- 40

tdFolder <- file.path("..", "Cancer_HiC_data_TAD_DA")

plotType <- "png"
myHeight <- myWidth <- 7

### PREP YL data
ylData_file <- file.path(setDir, "/mnt/ndata/Yuanlong/1.Projects/2.PROFILE/0.Scripts/Generate_TADs_for_Marie/bin_comp_CELL_LINEs_list.Rdata")
stopifnot(file.exists(ylData_file))
ylData <- eval(parse(text = load(ylData_file)))

### PREP all genome gene data
genome_geneFile <- file.path(setDir, "mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt")
geneDT <- read.delim(genome_geneFile, header=T, stringsAsFactors = F)
colnames(geneDT) <- c("entrezID", "chromo", "start", "end", "assembly", "strand")
geneDT <- geneDT[,c("entrezID", "chromo", "start", "end", "strand")]
stopifnot(is.numeric(geneDT$start) & is.numeric(geneDT$end))
geneDT <- geneDT[order(geneDT$chromo, geneDT$start),] 

# for each chromo -> find the biggest end
chromo_geneEnd <- c(by(geneDT, geneDT$chromo, function(x) max(x$end)))
stopifnot(!is.na(chromo_geneEnd))

chromo_binEnd <- ceiling(chromo_geneEnd/(binKb*1000))
stopifnot(chromo_geneEnd <= chromo_binEnd*(binKb*1000))
stopifnot(chromo_binEnd > 0)


all_hicds <- list.files(".", pattern = "all_assigned_regions.txt", full.names = TRUE, recursive = TRUE)
all_hicds <- basename(dirname(dirname(all_hicds)))
stopifnot(length(all_hicds) > 0)

##### > ITERATE OVERDATASETS

if(buildData) {
  hicds = "K562_40kb"
  
  all_hicds = all_hicds[1]
  
  all_ds_TAD_DT <- foreach(hicds = all_hicds, .combine='rbind') %dopar% {
    
    cat("... start ", hicds, "\n")
    
    if(hicds == "GSE75070_MCF-7_shNS_40kb") {
      hicds_0 <- "MCF-7_40kb"
    } else if(hicds == "ENCSR489OCU_NCI-H460_40kb") {
      hicds_0 <- "NCI-H460_40kb"
    } else {
      hicds_0 <- hicds
    }
    
    ## -----> YL Data
    yl_idx <- grep(gsub(paste0("_", binKb, "kb"), "", hicds), names(ylData))
    stopifnot(length(yl_idx) == 1)
    hic_ylDT <- ylData[[yl_idx]]
    stopifnot("chr" %in% colnames(hic_ylDT))
    stopifnot(nrow(hic_ylDT) > 0)
    hic_ylDT$chr <- paste0("chr", as.character(hic_ylDT$chr))
    colnames(hic_ylDT)[colnames(hic_ylDT) == "chr"] <- "chromo"
    hic_ylDT$bin_index <- as.numeric(as.character(hic_ylDT$bin_index))
    
    assignedFile <- file.path(hicds, "genes2tad", "all_assigned_regions.txt")
    assignedDT <- read.delim(assignedFile, header=F, stringsAsFactors = F, col.names=c("chromo", "region", "start", "end"))
    stopifnot(is.numeric(assignedDT$start))
    stopifnot(is.numeric(assignedDT$end))
    tad_assignedDT <- assignedDT[grepl("_TAD", assignedDT$region),]
    bd_assignedDT <- assignedDT[grepl("_BOUND", assignedDT$region),]
    
    g2tFile <- file.path(hicds, "genes2tad", "all_genes_positions.txt")
    stopifnot(file.exists(g2tFile))
    g2tDT <- read.delim(g2tFile, header=F, stringsAsFactors = F, col.names=c("entrezID", "chromo", "start", "end","region"))
    stopifnot(is.numeric(g2tDT$start))
    stopifnot(is.numeric(g2tDT$end))
    tad_geneDT <- g2tDT[grepl("_TAD", g2tDT$region),]
    bd_geneDT <- g2tDT[grepl("_BOUND", g2tDT$region),]
    stopifnot( ( nrow(bd_geneDT) + nrow(tad_geneDT) ) == nrow(g2tDT))
    
    ## -----> TopDom Data
    assignedFile_0 <- file.path(tdFolder, hicds_0, "genes2tad", "all_assigned_regions.txt")
    if(file.exists(assignedFile_0)) {
      avTopDomData <- TRUE
      assignedDT_0 <- read.delim(assignedFile_0, header=F, stringsAsFactors = F, col.names=c("chromo", "region", "start", "end"))
      stopifnot(is.numeric(assignedDT_0$start))
      stopifnot(is.numeric(assignedDT_0$end))
      tad_assignedDT_0 <- assignedDT_0[grepl("_TAD", assignedDT_0$region),]
      bd_assignedDT_0 <- assignedDT_0[grepl("_BOUND", assignedDT_0$region),]
      
      g2tFile_0 <- file.path(tdFolder, hicds_0, "genes2tad", "all_genes_positions.txt")
      stopifnot(file.exists(g2tFile_0))
      g2tDT_0 <- read.delim(g2tFile_0, header=F, stringsAsFactors = F, col.names=c("entrezID", "chromo", "start", "end","region"))
      stopifnot(is.numeric(g2tDT_0$start))
      stopifnot(is.numeric(g2tDT_0$end))
      tad_geneDT_0 <- g2tDT_0[grepl("_TAD", g2tDT_0$region),]
      bd_geneDT_0 <- g2tDT_0[grepl("_BOUND", g2tDT_0$region),]
      stopifnot( ( nrow(bd_geneDT_0) + nrow(tad_geneDT_0) ) == nrow(g2tDT_0))
    } else {
      avTopDomData <- FALSE
    }
    
    all_chromo <- unique(as.character(hic_ylDT$chromo))
    chromo = "chr1"
    
    # all_chromo = all_chromo[1]
    
    chr_ds_TAD_DT <- foreach(chromo = all_chromo, .combine='rbind') %do% {
      
      cat("... start ", hicds, " - ", chromo, "\n")
      
      stopifnot(chromo %in% names(chromo_binEnd))
      chr_all_ref_bins <- c(1:chromo_binEnd[paste0(chromo)])
      
      tadFile <- file.path(hicds, "FINAL_DOMAINS", paste0(gsub(paste0("_", binKb, "kb"), "", hicds), "_", chromo, "_", "YL_", binKb, "kb_final_domains.txt"))
      stopifnot(file.exists(tadFile))
      tadDT <- read.delim(tadFile, header=F, stringsAsFactors = F, col.names = c("chromo", "start", "end"))
      stopifnot(is.numeric(tadDT$start))
      stopifnot(is.numeric(tadDT$end))
      
      if(avTopDomData) {
        
        # if(hicds == "K562_40kb")  {
        #   
        # } else {
        #   
        # }
        tadFile_0 <- file.path(tdFolder, hicds_0, "FINAL_DOMAINS", paste0(gsub(paste0("_", binKb, "kb"), "", hicds_0), "_", chromo, "_", "KR_", binKb, "kb_final_domains.txt"))
        stopifnot(file.exists(tadFile_0))
        tadDT_0 <- read.delim(tadFile_0, header=F, stringsAsFactors = F, col.names = c("chromo", "start", "end"))
        stopifnot(is.numeric(tadDT_0$start))
        stopifnot(is.numeric(tadDT_0$end))
        tadDT_0$end <- tadDT_0$end
        
        chr_tad_assignedDT_0 <- tad_assignedDT_0[tad_assignedDT_0$chromo == chromo,]
        tmp_chr_tad_assignedDT_0 <- chr_tad_assignedDT_0
        chr_tad_assignedDT_0$region <- NULL
        #-> check that the same number of TADs assigned than retrieved from YL data (same TADs should be kept in the assigned_regions file)
        stopifnot(nrow(chr_tad_assignedDT_0) == nrow(tadDT_0))
        stopifnot(all.equal(chr_tad_assignedDT_0, tadDT_0, check.attributes=FALSE)) # rownames differ if check.attributes = TRUE
        chr_tad_assignedDT_0 <- tmp_chr_tad_assignedDT_0
        
        # to have similar to YL -> all binIdx covered with tad
        bin_tadDT_0 <- tadDT_0
        # bin_tadDT_0$start_bin <- floor((bin_tadDT_0$start - 1)/(binKb*1000)) + 1 # same as:
        bin_tadDT_0$start_bin <- (bin_tadDT_0$start - 1) %/% (binKb*1000) + 1
        bin_tadDT_0$end_bin <- bin_tadDT_0$end / (binKb * 1000)
        stopifnot( (((bin_tadDT_0$start_bin-1) * binKb * 1000) + 1) == bin_tadDT_0$start)
        stopifnot( (bin_tadDT_0$end_bin * binKb * 1000) == bin_tadDT_0$end)
        
        all_bins_0 <- unlist(apply(bin_tadDT_0, 1, function(x) x["start_bin"]:x["end_bin"]))
        stopifnot(!duplicated(all_bins_0))
        
        # foo_tadDT_0 <- data.frame(start = c(1, 81, 161, 241), 
        #                           end = c(80,  160, 200, 280))
        # foo_tadDT_0$start_bin <- (foo_tadDT_0$start - 1) %/% (40) + 1
        # foo_tadDT_0$end_bin <- foo_tadDT_0$end / (40)
        
        ####################################################################################################################
        # LOOK AT THE NUMBER OF BINS LOST (TopDom data)
        ####################################################################################################################
        missing_bins_0 <- chr_all_ref_bins[ ! chr_all_ref_bins %in% all_bins_0]  
        
        ####################################################################################################################
        # # of genes in TADs and number of genes in boundaries (YL data)
        ####################################################################################################################
        chr_tad_geneDT_0 <- tad_geneDT_0[tad_geneDT_0$chromo == chromo,]
        stopifnot(nrow(tad_geneDT_0) > 0)
        
        chr_bd_geneDT_0 <- bd_geneDT_0[bd_geneDT_0$chromo == chromo,]
        stopifnot(nrow(bd_geneDT_0) > 0)
        
        chr_bd_assignedDT_0 <- bd_assignedDT_0[bd_assignedDT_0$chromo == chromo,]
        stopifnot(nrow(chr_bd_assignedDT_0) > 0)
        
        chr_tad_assignedDT_0 <- tad_assignedDT_0[tad_assignedDT_0$chromo == chromo,]
        stopifnot(nrow(chr_tad_assignedDT_0) > 0)
        
        nGenesTAD_0 <- nrow(chr_tad_geneDT_0)
        nGenesBound_0 <- nrow(chr_bd_geneDT_0)
        
        meanTADsize_0 <- mean(chr_tad_assignedDT_0$end - chr_tad_assignedDT_0$start + 1)
        medianTADsize_0 <- median(chr_tad_assignedDT_0$end - chr_tad_assignedDT_0$start + 1)
        
        nAssignedBound_TD <- nrow(chr_bd_assignedDT_0)
        nAssignedTAD_TD <- nrow(chr_tad_assignedDT_0)
        nBinsCovered_TD <- length(all_bins_0)
        missingBins_TD <- length(missing_bins_0)
        nGenesTAD_TD <- nGenesTAD_0
        nGenesBound_TD <- nGenesBound_0
        meanSizeTAD_TD <- meanTADsize_0
        medianSizeTAD_TD <- medianTADsize_0
        
      } else {   # TopDom data not available
        nAssignedBound_TD <- NA
        nAssignedTAD_TD <- NA
        nBinsCovered_TD <- NA
        missingBins_TD <- NA
        nGenesTAD_TD <- NA
        nGenesBound_TD <- NA
        meanSizeTAD_TD <- NA
        medianSizeTAD_TD <- NA
      }
      
      chr_tad_geneDT <- tad_geneDT[tad_geneDT$chromo == chromo,]
      stopifnot(nrow(chr_tad_geneDT) > 0)
      
      chr_bd_geneDT <- bd_geneDT[bd_geneDT$chromo == chromo,]
      stopifnot(nrow(chr_bd_geneDT) > 0)
      
      chr_assignedDT <- assignedDT[assignedDT$chromo == chromo,]
      stopifnot(nrow(chr_assignedDT) > 0)
      
      chr_tad_assignedDT <- tad_assignedDT[tad_assignedDT$chromo == chromo,]
      stopifnot(nrow(chr_tad_assignedDT) > 0)
      
      chr_bd_assignedDT <- bd_assignedDT[bd_assignedDT$chromo == chromo,]
      stopifnot(nrow(chr_bd_assignedDT) > 0)
      
      tmp_chr_tad_assignedDT <- chr_tad_assignedDT
      chr_tad_assignedDT$region <- NULL
      #-> check that the same number of TADs assigned than retrieved from YL data (same TADs should be kept in the assigned_regions file)
      stopifnot(nrow(chr_tad_assignedDT) == nrow(tadDT))
      stopifnot(all.equal(chr_tad_assignedDT, tadDT, check.attributes=FALSE)) # rownames differ if check.attributes = TRUE
      
      # keep back the region column...
      chr_tad_assignedDT <- tmp_chr_tad_assignedDT
      
      ####################################################################################################################
      # LOOK AT THE NUMBER OF BINS LOST (YL data)
      ####################################################################################################################
      chr_hic_ylDT <- hic_ylDT[hic_ylDT$chromo == chromo, ]
      stopifnot(nrow(chr_hic_ylDT) > 0)
      stopifnot(!duplicated(chr_hic_ylDT$bin_index))
      all_bins <- chr_hic_ylDT$bin_index
      missing_bins <- chr_all_ref_bins[ ! chr_all_ref_bins %in% all_bins]
      ####################################################################################################################
      # # of genes in TADs and number of genes in boundaries (YL data)
      ####################################################################################################################
      nGenesTAD <- nrow(chr_tad_geneDT)
      nGenesBound <- nrow(chr_bd_geneDT)
      
      
      meanTADsize <- mean(chr_tad_assignedDT$end - chr_tad_assignedDT$start + 1)
      
      save(chr_tad_assignedDT, file = paste0(chromo, "_", "chr_tad_assignedDT.Rdata"))
      save(chr_hic_ylDT, file = paste0(chromo, "_", "chr_hic_ylDT.Rdata"))
      
      head(chr_tad_assignedDT)
      
      medianTADsize <- median(chr_tad_assignedDT$end - chr_tad_assignedDT$start + 1)
      
      data.frame(
        hicds = hicds,
        chromo = chromo,
        
        nAssignedBound_YL = nrow(chr_bd_assignedDT),
        nAssignedBound_TD = nAssignedBound_TD,
        
        nAssignedTAD_YL = nrow(chr_tad_assignedDT),
        nAssignedTAD_TD = nAssignedTAD_TD,
        
        missingBins_YL = length(missing_bins),
        missingBins_TD = missingBins_TD,
        
        nBinsCovered_YL = length(all_bins),
        nBinsCovered_TD = nBinsCovered_TD,
        
        nGenesTAD_YL = nGenesTAD,
        nGenesTAD_TD = nGenesTAD_TD,
        
        nGenesBound_YL= nGenesBound,
        nGenesBound_TD = nGenesBound_TD,
        
        meanSizeTAD_YL = meanTADsize,
        meanSizeTAD_TD = meanSizeTAD_TD,
        
        medianSizeTAD_YL = medianTADsize,
        medianSizeTAD_TD = medianSizeTAD_TD,
        
        stringsAsFactors = FALSE
      )
    } # end-foreach iterating over chromo
    chr_ds_TAD_DT
  } # end-foreach iterating over datasets
  outFile <- file.path(outFold, "all_ds_TAD_DT.Rdata")
  save(all_ds_TAD_DT, file = outFile)
  cat(paste0("... written: ", outFile, "\n"))
} else {  # not buildData
  outFile <- file.path(outFold, "all_ds_TAD_DT.Rdata")
  all_ds_TAD_DT <- eval(parse(text = load(outFile)))
}



#################################################################################################################### START PLOTING

# plot the different variables to compare YL and TD

# rehsape dataframe
# xx <- all_ds_TAD_DT[1:10,c("hicds", "chromo", "nAssignedBound_YL", "nAssignedBound_TD")]
# reshape(xx, ids = c("hicds", "chromo"), times = c("nAssignedBound_YL", "nAssignedBound_TD"), 
#           timevar = "feature", direction="long", varying=list(c("nAssignedBound_YL", "nAssignedBound_TD")),
#         idvar=c("hicds", "chromo"), v.names=c("nAssigned"))
# 
# 
# xx <- all_ds_TAD_DT[1:10,c("hicds", "chromo", "nAssignedBound_YL", "nAssignedBound_TD", "missingBins_YL", "missingBins_TD")]
# reshape(xx, ids = c("hicds", "chromo"), times = c("nAssignedBound_YL", "nAssignedBound_TD", "missingBins_YL", "missingBins_TD" ), 
#         timevar = "feature", direction="long", 
#         varying=list(c("nAssignedBound_YL", "nAssignedBound_TD",  "missingBins_YL", "missingBins_TD")),
#         idvar=c("hicds", "chromo"), v.names=c("nbr"))
# 
# 
# all_vars <- colnames(all_ds_TAD_DT)[! colnames(all_ds_TAD_DT) %in% c("hicds", "chromo")]
# 
# xx <- all_ds_TAD_DT[1:10,]
# 
# yy <- reshape(xx, 
#         ids = c("hicds", "chromo"), idvar = c("hicds", "chromo"),
#         times = all_vars,
#         timevar = "feature",
#         direction = "long",
#         v.names="nbr",
#         varying = list(all_vars)
#         )
# unique(yy$feature)



idVars <- c("hicds", "chromo")
varVars <- colnames(all_ds_TAD_DT)[! colnames(all_ds_TAD_DT) %in% c("hicds", "chromo")]

reshaped_DT <- reshape(all_ds_TAD_DT, 
                       ids = idVars, 
                       idvar = idVars,
                       times = varVars, 
                       timevar = "feature", # name of the new column that contains the category
                       direction = "long",
                       v.names="nbr", # name of the new column that contains the values
                       varying = list(varVars))

unique(reshaped_DT$feature)

reshaped_DT$featureName <- gsub("(.+)_.+", "\\1", reshaped_DT$feature)
reshaped_DT$featureDataset <- gsub(".+_(.+)", "\\1", reshaped_DT$feature)

plotList <- list(
  c("nAssignedBound", "nAssignedTAD"),
  c("missingBins", "nBinsCovered"),
  c("nGenesBound", "nGenesTAD"),
  c("meanSizeTAD", "medianSizeTAD")
                 )

for(toplot_vars in plotList) {
  
  if( all(toplot_vars == c("nAssignedBound", "nAssignedTAD")) | all(toplot_vars ==  c("nGenesBound", "nGenesTAD"))) {
    myylab <- "# genes"
  } else if(all(toplot_vars ==  c("meanSizeTAD", "medianSizeTAD"))) {
    myylab <- "TAD size (bp)"
  } else if(all(toplot_vars == c("missingBins", "nBinsCovered"))) {
    myylab <- "# bins"
  } else {
    myylab <- NA
  }
  stopifnot(!is.na(myylab))
  
  sub_DT <- reshaped_DT[reshaped_DT$featureName %in% toplot_vars,]
  stopifnot(nrow(sub_DT) > 0)
  
  nDS <- length(unique(sub_DT$hicds))

  
  outFile <- file.path(outFold, paste0(paste0(toplot_vars, collapse = "_"), "_YL_TD_boxplot.", plotType))
  p <- ggplot(sub_DT, aes(x = featureName, y= nbr, fill = featureDataset, col=featureDataset))  +
    geom_boxplot()+
    scale_fill_ipsum() +
    scale_color_ipsum() +
    theme_ipsum_rc(grid="XY", axis_title_just="c", axis_title_size=14) + 
    labs(x="", 
         y=paste0(myylab),
         title="Comparison YL vs. TopDom TADs",
         # caption="Brought to you by the letter 'g'",
         subtitle=paste0("nDS = ", nDS),
         colour = "",
         fill = ""
         ) 
  ggsave(p, filename = outFile, height = myHeight, width = myWidth)
  cat(paste0("... written: ", outFile, "\n"))
    
  # THEMES FROM hrbrthemes package
  # theme_ft_rc(grid="XY", axis_title_just="c") # dark
  # theme_ipsum(grid="XY", axis_title_just="c")
  #theme_ipsum_ps(grid="XY", axis_title_just="c") !
  # theme_ipsum_pub(grid="XY", axis_title_just="c")
  # theme_ipsum_tw(grid="XY", axis_title_just="c")
  # theme_modern_rc(grid="XY", axis_title_just="c") # dark
  
  # THEMES FROM ggthemes package
  # theme_base()
  # theme_calc()
  # theme_economist()  # colored
  # theme_economist_white() # colored bg
  # theme_excel() # colored bg
  # theme_few()# colored bg
  # theme_fivethirtyeight()# colored bg
  # theme_foundation()
  # theme_gdocs()   
  # theme_hc()
  # theme_igray()# colored bg
  # theme_map    ()        
  # theme_pander()
  # theme_par()
  # theme_solarized()      # colored bg
  # theme_solarized_2() # colored bg
  # theme_solid()
  # theme_stata()       
  # theme_tufte() # nice font !
  
  
}

# ######################################################################################
# ######################################################################################
# ######################################################################################
cat(paste0("*** DONE: ", script_name, "\n"))
cat(paste0(startTime, "\n", Sys.time(), "\n"))








