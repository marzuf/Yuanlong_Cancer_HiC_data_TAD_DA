#!/usr/bin/bash


######## #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 


# LIVER
# ./run_pipeline.sh GSE105381_HepG2_40kb TCGAlihc_norm_lihc   # 10.06.2019
# ./run_pipeline.sh GSE105381_HepG2_40kb TCGAlihc_wt_mutCTNNB1   # 12.06.2019

# ./run_pipeline.sh GSE58752_liver_40kb TCGAlihc_norm_lihc # ADDED 20.06.2019
# ./run_pipeline.sh GSE58752_liver_40kb TCGAlihc_wt_mutCTNNB1 # ADDED 20.06.2019


# BREAST
# ./run_pipeline.sh MCF-7_40kb TCGAbrca_lum_bas   # does not exist 
# ./run_pipeline.sh ENCSR549MGQ_T47D_40kb TCGAbrca_lum_bas   # 10.06.2019

# ./run_pipeline.sh GSE75070_MCF-7_shNS_40kb TCGAbrca_lum_bas # ADDED 20.06.2019


# KIDNEY
# ./run_pipeline.sh ENCSR079VIJ_G401_40kb TCGAkich_norm_kich  # 10.06.2019
# ./run_pipeline.sh ENCSR401TBQ_Caki2_40kb TCGAkich_norm_kich # 11.06.2019


# SKIN
# ./run_pipeline.sh ENCSR312KHQ_SK-MEL-5_40kb TCGAskcm_lowInf_highInf # 10.06.2019
# ./run_pipeline.sh ENCSR862OGI_RPMI-7951_40kb TCGAskcm_lowInf_highInf  # 11.06.2019

# ./run_pipeline.sh ENCSR312KHQ_SK-MEL-5_40kb TCGAskcm_wt_mutBRAF # 10.06.2019
# ./run_pipeline.sh ENCSR862OGI_RPMI-7951_40kb TCGAskcm_wt_mutBRAF # 11.06.2019  

# ./run_pipeline.sh ENCSR312KHQ_SK-MEL-5_40kb TCGAskcm_wt_mutCTNNB1		# 10.06.2019
# ./run_pipeline.sh ENCSR862OGI_RPMI-7951_40kb TCGAskcm_wt_mutCTNNB1		# 11.06.2019 


# LUNG
# ./run_pipeline.sh ENCSR444WCZ_A549_40kb TCGAluad_norm_luad   # 10.06.2019
# ./run_pipeline.sh NCI-H460_40kb TCGAluad_norm_luad   # does not exist
# ./run_pipeline.sh ENCSR489OCU_NCI-H460_40kb TCGAluad_norm_luad   # 12.06.2019

# ./run_pipeline.sh ENCSR444WCZ_A549_40kb TCGAluad_mutKRAS_mutEGFR # 10.06.2019
# ./run_pipeline.sh NCI-H460_40kb TCGAluad_mutKRAS_mutEGFR  # does not exist 
# ./run_pipeline.sh ENCSR489OCU_NCI-H460_40kb TCGAluad_mutKRAS_mutEGFR  # 11.06.2019 

# ./run_pipeline.sh ENCSR444WCZ_A549_40kb TCGAluad_nonsmoker_smoker # 11.06.2019
# ./run_pipeline.sh NCI-H460_40kb TCGAluad_nonsmoker_smoker  # does not exist
# ./run_pipeline.sh ENCSR489OCU_NCI-H460_40kb TCGAluad_nonsmoker_smoker  # 11.06.2019 

# ./run_pipeline.sh ENCSR444WCZ_A549_40kb TCGAluad_wt_mutKRAS # 11.06.2019
# ./run_pipeline.sh NCI-H460_40kb TCGAluad_wt_mutKRAS # does not exist
# ./run_pipeline.sh ENCSR489OCU_NCI-H460_40kb TCGAluad_wt_mutKRAS # 11.06.2019 

# ./run_pipeline.sh ENCSR444WCZ_A549_40kb TCGAlusc_norm_lusc  # 10.06.2019
# ./run_pipeline.sh NCI-H460_40kb TCGAlusc_norm_lusc # does not exist
# ./run_pipeline.sh ENCSR489OCU_NCI-H460_40kb TCGAlusc_norm_lusc # 11.06.2019 

	
# PANCREAS
# ./run_pipeline.sh Panc1_rep12_40kb TCGApaad_wt_mutKRAS # 10.06.2019
	

# PROSTATE	
# ./run_pipeline.sh ENCSR346DCU_LNCaP_40kb TCGAprad_norm_prad # 10.06.2019
# ./run_pipeline.sh GSE73782_PC3_40kb TCGAprad_norm_prad # 12.06.2019

# ./run_pipeline.sh GSE118514_RWPE1_40kb TCGAprad_norm_prad # ADDED 20.06.2019

# GBM
# ./run_pipeline.sh GSE105194_cerebellum_40kb TCGAgbm_classical_mesenchymal # 12.06.2019
# ./run_pipeline.sh GSE105194_cerebellum_40kb TCGAgbm_classical_neural # 12.06.2019
# ./run_pipeline.sh GSE105194_cerebellum_40kb TCGAgbm_classical_proneural # 12.06.2019
# ./run_pipeline.sh GSE105194_cerebellum_40kb TCGAlgg_IDHwt_IDHmutnc # 12.06.2019

# ./run_pipeline.sh GSE105194_spinal_cord_40kb TCGAgbm_classical_mesenchymal # 12.06.2019
# ./run_pipeline.sh GSE105194_spinal_cord_40kb TCGAgbm_classical_neural # 12.06.2019
# ./run_pipeline.sh GSE105194_spinal_cord_40kb TCGAgbm_classical_proneural # 12.06.2019
# ./run_pipeline.sh GSE105194_spinal_cord_40kb TCGAlgg_IDHwt_IDHmutnc # 12.06.2019


# COLORECTAL
# ./run_pipeline.sh GSE105318_DLD1_40kb TCGAcoad_msi_mss # 10.06.2019

# LYMPHOBLAST
# ./run_pipeline.sh K562_40kb TCGAlaml_wt_mutFLT3 # 10.06.2019



start_time=$(date -R)    
set -e

if [[ $# != 2 ]]; then
    echo "invalid # of arguments"
    exit 1
fi

hic_dataset="$1"
expr_dataset="$2"

echo "*** START ***"
echo "... > Hi-C dataset: $hic_dataset"
echo "... > Gene expression dataset: $expr_dataset"

#********************** HARD-CODED SETTINGS FOR THE PIPELINE ********************************************
step1=1     # prepare setting file
step2=1    # run the pipeline

TAD_DE_pipSteps=( "0cleanInputTCGA" "1cleanInput" "2" "2v2" "3" "4" "5" "6" "7" "8c" "9" "10" "11" "13cleanInput" "14f2" "170revision2EZH2" )
#TAD_DE_pipSteps=( "7" "8c" "9" "10" "11" "13cleanInput" "14f2" "170revision2EZH2" )
#TAD_DE_pipSteps=( "0cleanInputTCGA" )
#TAD_DE_pipSteps=( "0cleanInputTCGA" "1cleanInput" "2" "3" )
#TAD_DE_pipSteps=( "4" "5" "6" "7" "8c")

# new steps 03.03.2019 à refaire: 2v2, 5v2, 6v2, 7v2, 9v2, 10v2, 11v2
#TAD_DE_pipSteps=( "6v2" "9v2" )
#TAD_DE_pipSteps=( "5v2" )
#TAD_DE_pipSteps=( "6v2onlyW" )
#TAD_DE_pipSteps=( "7v2" )
#TAD_DE_pipSteps=( "9v2" )
#TAD_DE_pipSteps=( "10v2" "10b" )

# ./run_pipeline.sh GSE105381_HepG2_40kb TCGAlihc_norm_lihc   # 


# ./run_pipeline.sh ENCSR079VIJ_G401_40kb TCGAkich_norm_kich   # 

# ./run_pipeline.sh ENCSR444WCZ_A549_40kb TCGAluad_norm_luad

runDir="/mnt/etemp/marie/Yuanlong_Cancer_HiC_data_TAD_DA"

TAD_DE_pipDir="/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom"
TAD_DE_script="./zzz_run_given_step_given_data_v2.sh"
old_inputFolder="/mnt/ed4/marie/scripts/TAD_DE_pipeline/SETTING_FILES_cleanInput"
nPermut="10000"
#nPermut="10"
ncpu="40"

Rexec=`which Rscript`

new_inputFolder="$runDir/PIPELINE/INPUT_FILES/$hic_dataset"
mkdir -p $new_inputFolder
outputFolder="$runDir/PIPELINE/OUTPUT_FOLDER/$hic_dataset/$expr_dataset"
mkdir -p $outputFolder

echo "!!! IMPORTANT HARD-CODED SETTINGS !!!"
echo "... ! step1 = $step1"
echo "... ! step2 = $step2"

echo "... ! old_inputFolder = $old_inputFolder"
echo "... ! nPermut = $nPermut"
echo "... ! ncpu = $ncpu"

if [[ $step2 -eq 1 ]]; then
	echo "... ! TAD_DE_pipDir = $TAD_DE_pipDir"
	echo "... ! TAD_DE_script = $TAD_DE_script"
	echo "... ! TAD_DE_pipStep(s): ${TAD_DE_pipSteps[*]}"
fi

###################################### FUNCTION DEFINITIONS

function mvBack {
  echo "... go back to my folder"
  cd $runDir  
}
trap mvBack EXIT

runCMD() {
  echo "> $1"
  eval $1
}

checkFile() {
if [[ ! -f  $2 ]]; then
    echo "... $1 ($2) does not exist !"
    exit 1
fi
}
############################################################################

old_setting_file="$old_inputFolder/run_settings_${expr_dataset}.R"

checkFile old_setting_file $old_setting_file

new_setting_file="$new_inputFolder/run_settings_${expr_dataset}.R"

runCMD "cp $old_setting_file $new_inputFolder"

# !!! TO CHECK FORMAT ZZZZZ
# $hic_dataset/genes2tad/all_assigned_regions.txt
# MCF-7/genes2tad/all_assigned_regions.txt
# TADposDT.txt
#chr1    chr1_TAD1       750001  1300000
new_TADpos_file="$hic_dataset/genes2tad/all_assigned_regions.txt"

# !!! TO CHECK FORMAT ZZZZZ
# $hic_dataset/genes2tad/all_genes_positions.txt
# MCF-7/genes2tad/all_genes_positions.txt
#gene2tadDT.txt
#12893       chr1    761586  762902  chr1_TAD1
new_gene2tad_file="$hic_dataset/genes2tad/all_genes_positions.txt"

checkFile new_TADpos_file $new_TADpos_file

checkFile new_gene2tad_file $new_gene2tad_file

if [[ "$step1" -eq 1 ]] ; then

	cat >> ${new_setting_file} <<- EOM

			# > file edited: `date -R` 

			# path to output folder:
			pipOutFold <- "${outputFolder}"

			# OVERWRITE THE DEFAULT SETTINGS FOR INPUT FILES - use TADs from the current Hi-C dataset 
			TADpos_file <- paste0(setDir, "`realpath $new_TADpos_file`")
							#chr1    chr1_TAD1       750001  1300000
							#chr1    chr1_TAD2       2750001 3650000
							#chr1    chr1_TAD3       3650001 4150000

			gene2tadDT_file <- paste0(setDir, "`realpath $new_gene2tad_file`")
							#LINC00115       chr1    761586  762902  chr1_TAD1
							#FAM41C  chr1    803451  812283  chr1_TAD1
							#SAMD11  chr1    860260  879955  chr1_TAD1
							#NOC2L   chr1    879584  894689  chr1_TAD1

			# overwrite main_settings.R: nCpu <- 25
			nCpu <- ${ncpu}

			# *************************************************************************************************************************
			# ************************************ SETTINGS FOR PERMUTATIONS (5#_, 8c_)
			# *************************************************************************************************************************

			# number of permutations
			nRandomPermut <- $nPermut
			gene2tadAssignMethod <- "maxOverlap"
			nRandomPermutShuffle <- $nPermut
			step8_for_permutGenes <- TRUE
			step8_for_randomTADsFix <- FALSE
			step8_for_randomTADsGaussian <- FALSE
			step8_for_randomTADsShuffle <- FALSE
			step14_for_randomTADsShuffle <- FALSE

	EOM
	echo "WRITTEN: ${new_setting_file}"


echo "> END STEP1:" $(date -R)

checkFile new_setting_file $new_setting_file

fi # end if STEP1

###################################################################################################################################################
### STEP2: MOVE TO THE TAD DE PIPELINE DIRECTORY TO LAUNCH THE TAD DE PIPELINE
###################################################################################################################################################
if [[ "$step2" -eq 1 ]] ; then
	echo "> START STEP2:" $(date -R)
	runCMD "cd $TAD_DE_pipDir"

	echo $TAD_DE_script ${new_setting_file} `echo ${TAD_DE_pipSteps[*]}`
	$TAD_DE_script ${new_setting_file} `echo ${TAD_DE_pipSteps[*]}`

	cd $runDir
	echo "> END STEP2:" $(date -R)
fi

###################################################################################################################################################
########## END ####################################################################################################################################
###################################################################################################################################################
echo "*** DONE"
end_time=$(date -R)    
echo $start_time
echo $end_time
exit 0

