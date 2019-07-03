# Rscript comp_coexprCorr_coexprAngDist.R

script_name <- "comp_coexprCorr_coexprAngDist.R"

buildDT <- FALSE

startTime <- Sys.time()
cat(paste0("> Rscript comp_coexprCorr_coexprAngDist.R\n"))

# angDistFiles <- list.files("CREATE_ANGDIST_BYCOND_SORTNODUP", full.names = TRUE, pattern = "all_angDist_DT.Rdata", recursive = TRUE)
# stopifnot(length(angDistFiles) > 0)
# angDistFiles_cond1 <- list.files("CREATE_ANGDIST_BYCOND_SORTNODUP", full.names = TRUE, pattern = "all_angDist_DT_cond1.Rdata", recursive = TRUE)
# stopifnot(length(angDistFiles_cond1) > 0)
# angDistFiles_cond2 <- list.files("CREATE_ANGDIST_BYCOND_SORTNODUP", full.names = TRUE, pattern = "all_angDist_DT_cond2.Rdata", recursive = TRUE)
# stopifnot(length(angDistFiles_cond2) > 0)
# 
# 
# coexprFiles <- list.files("CREATE_COEXPR_SORTNODUP", full.names = TRUE, pattern = "coexprDT.Rdata", recursive = TRUE)
# stopifnot(length(coexprFiles) > 0)
# coexprFiles_cond1 <- list.files("CREATE_COEXPR_BYCOND_SORTNODUP", full.names = TRUE, pattern = "coexprDT_cond1.Rdata", recursive = TRUE)
# stopifnot(length(coexprFiles_cond1) > 0)
# coexprFiles_cond2 <- list.files("CREATE_COEXPR_BYCOND_SORTNODUP", full.names = TRUE, pattern = "coexprDT_cond2.Rdata", recursive = TRUE)
# stopifnot(length(coexprFiles_cond2) > 0)
# 
# 
# sameTADfiles <- list.files("CREATE_SAME_TAD_SORTNODUP", full.names = TRUE, pattern = "all_TAD_pairs.Rdata", recursive = TRUE)

################################################################################################################################
################################################################################################################################
################################################################################################################################

options(scipen=100)

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")
require(foreach)
require(doMC)

registerDoMC(40)

outFold <- file.path("COMP_COEXPRCORR_COEXPRANGDIST")
dir.create(outFold, recursive = TRUE)
plotType <- "png"
myHeight <- ifelse(plotType == "png", 400, 7)
myWidth <- ifelse(plotType == "png", 600, 10)


sameTADfolder <- "CREATE_SAME_TAD_SORTNODUP"
stopifnot(dir.exists(sameTADfolder))

angDistFolder <- "CREATE_ANGDIST_BYCOND_SORTNODUP"
stopifnot(dir.exists(angDistFolder))

coexprFolder <- "CREATE_COEXPR_SORTNODUP"
stopifnot(dir.exists(coexprFolder))

coexprByCondFolder <-"CREATE_COEXPR_BYCOND_SORTNODUP"
stopifnot(dir.exists(coexprByCondFolder))

corMet <- "pearson"

all_datasets <- list.dirs(file.path("PIPELINE", "OUTPUT_FOLDER"), recursive = TRUE)

hicds <- basename(dirname(dirname(all_datasets)))
exprds <- basename(dirname(all_datasets))


pip_datasets <- file.path(hicds, exprds)
pip_datasets <- pip_datasets[dir.exists(hicds) & 
                               pip_datasets != file.path("PIPELINE", "OUTPUT_FOLDER") &
                               pip_datasets != file.path(".", "PIPELINE")
                               ]
pip_datasets <- unique(pip_datasets)

dataset = pip_datasets[1]

# pip_datasets = pip_datasets[1]

if(buildDT) {
  all_corr_angDist_DT <- foreach(dataset = pip_datasets) %dopar% {
    
    hicds <- dirname(dataset)
    stopifnot(dir.exists(hicds))
    exprds <- basename(dataset)
    
    cat(paste0("... START ", hicds, " - ", exprds, "\n"))
    
    # CREATE_ANGDIST_BYCOND_SORTNODUP/ENCSR346DCU_LNCaP_40kb/TCGAprad_norm_prad/pearson/all_angDist_DT_cond2.Rdata
    # CREATE_COEXPR_BYCOND_SORTNODUP/K562_40kb/TCGAlaml_wt_mutFLT3/pearson/coexprDT_cond1.Rdata
    # CREATE_SAME_TAD_SORTNODUP/ENCSR079VIJ_G401_40kb/all_TAD_pairs.Rdata
    
    angDistFile <- file.path(angDistFolder, hicds, exprds, corMet, "all_angDist_DT.Rdata")
    stopifnot(file.exists(angDistFile))
    angDistFile_cond1 <- file.path(angDistFolder, hicds, exprds, corMet, "all_angDist_DT_cond1.Rdata")
    stopifnot(file.exists(angDistFile_cond1))
    angDistFile_cond2 <- file.path(angDistFolder, hicds, exprds, corMet, "all_angDist_DT_cond2.Rdata")
    stopifnot(file.exists(angDistFile_cond2))
    
    cat(paste0("...... load ", angDistFile, "\n"))
    angDist_DT <- eval(parse(text = load(angDistFile)))
    # tmp_angDist_DT <- angDist_DT[1:10,]
    # save(tmp_angDist_DT, file ="tmp_angDist_DT.Rdata")
    cat(paste0("...... load ", angDistFile_cond1, "\n"))
    angDist_DT_cond1 <- eval(parse(text = load(angDistFile_cond1)))
    # tmp_angDist_DT_cond1 <- angDist_DT_cond1[1:10,]
    # save(tmp_angDist_DT_cond1, file ="tmp_angDist_DT_cond1.Rdata")
    cat(paste0("...... load ", angDistFile_cond2, "\n"))
    angDist_DT_cond2 <- eval(parse(text = load(angDistFile_cond2)))
    # tmp_angDist_DT_cond2 <- angDist_DT_cond2[1:10,]
    # save(tmp_angDist_DT_cond2, file ="tmp_angDist_DT_cond2.Rdata")
    
    coexprCorrFile <- file.path(coexprFolder, hicds, exprds, corMet, "coexprDT.Rdata")
    stopifnot(file.exists(coexprCorrFile))
    coexprCorrFile_cond1 <- file.path(coexprByCondFolder, hicds, exprds, corMet, "coexprDT_cond1.Rdata")
    stopifnot(file.exists(coexprCorrFile_cond1))
    coexprCorrFile_cond2 <- file.path(coexprByCondFolder, hicds, exprds, corMet, "coexprDT_cond2.Rdata")
    stopifnot(file.exists(coexprCorrFile_cond2))
    
    cat(paste0("...... load ", coexprCorrFile, "\n"))
    coexpr_DT <- eval(parse(text = load(coexprCorrFile)))
    # tmp_coexpr_DT <- coexpr_DT[1:10,] 
    # save(tmp_coexpr_DT, file ="tmp_coexpr_DT.Rdata")
    cat(paste0("...... load ", coexprCorrFile_cond1, "\n"))
    coexpr_DT_cond1 <- eval(parse(text = load(coexprCorrFile_cond1)))
    # tmp_coexpr_DT_cond1 <- coexpr_DT_cond1[1:10,]
    # save(tmp_coexpr_DT_cond1, file ="tmp_coexpr_DT_cond1.Rdata")
    cat(paste0("...... load ", coexprCorrFile_cond2, "\n"))
    coexpr_DT_cond2 <- eval(parse(text = load(coexprCorrFile_cond2)))
    # tmp_coexpr_DT_cond2 <- coexpr_DT_cond2[1:10,]
    # save(tmp_coexpr_DT_cond2, file ="tmp_coexpr_DT_cond2.Rdata")
    
    tadFile <- file.path(sameTADfolder, hicds, "all_TAD_pairs.Rdata")
    stopifnot(file.exists(tadFile))
    cat(paste0("...... load ", tadFile, "\n"))
    tad_DT <- eval(parse(text = load(tadFile)))
    stopifnot(grepl("_TAD", tad_DT$region))
    # tmp_tad_DT <- tad_DT[1:10,]
    # save(tmp_tad_DT, file ="tmp_tad_DT.Rdata")
    
    # load("tmp_angDist_DT.Rdata")
    # load("tmp_angDist_DT_cond1.Rdata")
    # load("tmp_angDist_DT_cond2.Rdata")
    # load("tmp_coexpr_DT.Rdata")
    # load("tmp_coexpr_DT_cond1.Rdata")
    # load("tmp_coexpr_DT_cond2.Rdata")
    # load("tmp_tad_DT.Rdata")
    
    
    cat("...... merge angDistDT and coexprDT", "\n" )
    stopifnot(angDist_DT$gene2 > angDist_DT$gene1)
    stopifnot(coexpr_DT$gene2 > coexpr_DT$gene1)
    angDist_coexpr_DT <- merge(angDist_DT, coexpr_DT, by = c("gene1", "gene2"), all = TRUE)
    stopifnot(!is.na(angDist_coexpr_DT))
    
    cat("...... merge angDist_coexpr_DT and tadDT", "\n" )
    angDist_coexpr_tad_DT <- merge(angDist_coexpr_DT, tad_DT, by = c("gene1", "gene2"), all = TRUE)
    # angDist_coexpr_tad_DT$sameTAD <- ! is.na(angDist_coexpr_tad_DT$region)
    angDist_coexpr_tad_DT$sameTADTxt <- ifelse(! is.na(angDist_coexpr_tad_DT$region), "sameTAD", "diffTAD")
    
    outFile <- file.path(outFold, paste0(hicds, "_", exprds, "_", "coexprAngDist_sameDiffTAD_boxplot.", plotType))
    do.call(plotType, list(outFile, height=myHeight, width=myWidth))
    boxplot(angDist_coexpr_tad_DT$angDist ~ angDist_coexpr_tad_DT$sameTADTxt,
            ylab = "angDist",
            main = paste0(hicds, " - ", exprds))
    mtext(side = 3, text = paste0("pairwise coexpr. ang. dist."), font = 3)
    
    t_test <- t.test(x=angDist_coexpr_tad_DT$angDist[angDist_coexpr_tad_DT$sameTADTxt == "sameTAD"],
                     y = angDist_coexpr_tad_DT$angDist[angDist_coexpr_tad_DT$sameTADTxt == "diffTAD"])
    legText <- paste0("t-test pval. =\n", sprintf("%.2f", t_test$p.value))
    legend("bottomright",bty="n", legend=legText )
    
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
    
    outFile <- file.path(outFold, paste0(hicds, "_", exprds, "_", "coexprCorr_sameDiffTAD_boxplot.", plotType))
    do.call(plotType, list(outFile, height=myHeight, width=myWidth))
    boxplot(angDist_coexpr_tad_DT$coexpr ~ angDist_coexpr_tad_DT$sameTADTxt,
            ylab = paste0(corMet, " corr."),
            main = paste0(hicds, " - ", exprds))
    mtext(side = 3, text = paste0("pairwise coexpr. corr."), font = 3)
    
    t_test <- t.test(x=angDist_coexpr_tad_DT$coexpr[angDist_coexpr_tad_DT$sameTADTxt == "sameTAD"],
                     y = angDist_coexpr_tad_DT$coexpr[angDist_coexpr_tad_DT$sameTADTxt == "diffTAD"])
    legText <- paste0("t-test pval. =\n", sprintf("%.2f", t_test$p.value))
    legend("bottomright",bty="n", legend=legText )
    
    
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
    
    outFile <- file.path(outFold, paste0(hicds, "_", exprds, "_", "coexprCorr_coexprAngDist_sameDiffTAD_densplot.", plotType))
    myx <- angDist_coexpr_tad_DT$angDist
    myy <- angDist_coexpr_tad_DT$coexpr
    do.call(plotType, list(outFile, height=myHeight, width=myHeight))
    densplot( y = myy,
              x = myx,
              ylab = paste0(corMet, " corr."),
              xlab = paste0("ang. dist."),
              main = paste0(hicds, " - ", exprds))
    addCorr(x=myx,y=myy,legPos="topright", bty='n')
    mtext(side = 3, text = ("pairwise coexpr. corr. vs. ang. dist."), font = 3)
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
    
    outFile <- file.path(outFold, paste0(hicds, "_", exprds, "_", "coexprCorr_coexprAngDist_sameDiffTAD_density.", plotType))
    do.call(plotType, list(outFile, height=myHeight, width=myWidth))
    plot_multiDens(
      list(
        sameTAD_coexprCorr = angDist_coexpr_tad_DT$coexpr[angDist_coexpr_tad_DT$sameTADTxt == "sameTAD"],
        diffTAD_coexprCorr = angDist_coexpr_tad_DT$coexpr[angDist_coexpr_tad_DT$sameTADTxt == "diffTAD"],
        sameTAD_coexprAngDist = angDist_coexpr_tad_DT$angDist[angDist_coexpr_tad_DT$sameTADTxt == "sameTAD"],
        diffTAD_coexprAngDist = angDist_coexpr_tad_DT$angDist[angDist_coexpr_tad_DT$sameTADTxt == "diffTAD"]
      ),
      legPos="topleft",
      plotTit = paste0(hicds, " - ", exprds))
    mtext(side = 3, text = ("pairwise coexpr. corr. and ang. dist."), font = 3)
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
    
    data.frame(
      cbind(data.frame(hicds = hicds, exprds = exprds),
            angDist_coexpr_tad_DT),
      stringsAsFactors = FALSE
    )
    
    # break
    
  } # end-foreach iterating over datasets
  # all_corr_angDist_DT
  
  outFile <- file.path(outFold, "all_corr_angDist_DT.Rdata")
  save(all_corr_angDist_DT, file = outFile)
  cat(paste0("... written: ", outFile, "\n"))
} else {
  outFile <- file.path(outFold, "all_corr_angDist_DT.Rdata")
  stopifnot(file.exists(outFile))
  cat("... load all_corr_angDist_DT \n")
  all_corr_angDist_DT <- eval(parse(text = load(outFile)))
}

nDS <- length(unique(paste0(all_corr_angDist_DT$hicds, all_corr_angDist_DT$exprds)))
myTit <- paste0("all datasets (n=", nDS, ")")

outFile <- file.path(outFold, paste0("all_", "coexprAngDist_sameDiffTAD_boxplot.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
boxplot(all_corr_angDist_DT$angDist ~ all_corr_angDist_DT$sameTADTxt,
        ylab = "angDist",
        main = myTit)
mtext(side = 3, text = paste0("pairwise coexpr. ang. dist."), font = 3)

t_test <- t.test(x=all_corr_angDist_DT$angDist[all_corr_angDist_DT$sameTADTxt == "sameTAD"],
                 y = all_corr_angDist_DT$angDist[all_corr_angDist_DT$sameTADTxt == "diffTAD"])
legText <- paste0("t-test pval. =\n", sprintf("%.2f", t_test$p.value))
legend("bottomright",bty="n", legend=legText )

foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

outFile <- file.path(outFold, paste0("all_", "coexprCorr_sameDiffTAD_boxplot.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
boxplot(all_corr_angDist_DT$coexpr ~ all_corr_angDist_DT$sameTADTxt,
        ylab = paste0(corMet, " corr."),
        main = myTit)
mtext(side = 3, text = paste0("pairwise coexpr. corr."), font = 3)

t_test <- t.test(x=all_corr_angDist_DT$coexpr[all_corr_angDist_DT$sameTADTxt == "sameTAD"],
                 y = all_corr_angDist_DT$coexpr[all_corr_angDist_DT$sameTADTxt == "diffTAD"])
legText <- paste0("t-test pval. =\n", sprintf("%.2f", t_test$p.value))
legend("bottomright",bty="n", legend=legText )


foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


######################################################################################
######################################################################################
######################################################################################
cat("*** DONE: ", script_name, "\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))
