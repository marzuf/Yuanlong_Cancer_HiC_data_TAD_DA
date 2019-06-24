#!/usr/bin/bash

# ./all_prep_data_for_AUC_coexprDist_coexpr_byCond.sh

start_time=$(date -R)   

scriptCoexpr="create_coexpr_byCond_sortNoDup_otherTADfile.R"

# should have hicds + expr ds !!!

all_TAD_files_ds=(
## DS launched 20.01.2018:
#"ENCSR079VIJ_G401_40kb TCGAkich_norm_kich"
## 2
#"ENCSR312KHQ_SK-MEL-5_40kb TCGAskcm_lowInf_highInf"
#"ENCSR312KHQ_SK-MEL-5_40kb TCGAskcm_wt_mutBRAF"
#"ENCSR312KHQ_SK-MEL-5_40kb TCGAskcm_wt_mutCTNNB1"
## 3
#"ENCSR346DCU_LNCaP_40kb TCGAprad_norm_prad"
## 4
#"ENCSR401TBQ_Caki2_40kb TCGAkich_norm_kich"
## 5
#"ENCSR444WCZ_A549_40kb TCGAluad_mutKRAS_mutEGFR"
#"ENCSR444WCZ_A549_40kb TCGAluad_nonsmoker_smoker"
#"ENCSR444WCZ_A549_40kb TCGAluad_norm_luad"
#"ENCSR444WCZ_A549_40kb TCGAluad_wt_mutKRAS"
#"ENCSR444WCZ_A549_40kb TCGAlusc_norm_lusc"
## 6
#"ENCSR489OCU_NCI-H460_40kb TCGAluad_mutKRAS_mutEGFR"
#"ENCSR489OCU_NCI-H460_40kb TCGAluad_nonsmoker_smoker"
#"ENCSR489OCU_NCI-H460_40kb TCGAluad_norm_luad"
#"ENCSR489OCU_NCI-H460_40kb TCGAluad_wt_mutKRAS"
#"ENCSR489OCU_NCI-H460_40kb TCGAlusc_norm_lusc"
## 7
#"ENCSR549MGQ_T47D_40kb TCGAbrca_lum_bas"
## 8
#"ENCSR862OGI_RPMI-7951_40kb TCGAskcm_lowInf_highInf"
#"ENCSR862OGI_RPMI-7951_40kb TCGAskcm_wt_mutBRAF"
#"ENCSR862OGI_RPMI-7951_40kb TCGAskcm_wt_mutCTNNB1"
## 9
#"GSE105194_cerebellum_40kb TCGAgbm_classical_mesenchymal"
#"GSE105194_cerebellum_40kb TCGAgbm_classical_neural"
#"GSE105194_cerebellum_40kb TCGAgbm_classical_proneural"
#"GSE105194_cerebellum_40kb TCGAlgg_IDHwt_IDHmutnc"
## 10
#"GSE105194_spinal_cord_40kb TCGAgbm_classical_mesenchymal"
#"GSE105194_spinal_cord_40kb TCGAgbm_classical_neural"
#"GSE105194_spinal_cord_40kb TCGAgbm_classical_proneural"
#"GSE105194_spinal_cord_40kb TCGAlgg_IDHwt_IDHmutnc"
## 11
#"GSE105318_DLD1_40kb TCGAcoad_msi_mss"
## 12
#"GSE105381_HepG2_40kb TCGAlihc_norm_lihc"
#"GSE105381_HepG2_40kb TCGAlihc_wt_mutCTNNB1"
## 13
#"GSE73782_PC3_40kb TCGAprad_norm_prad"
## 14
#"Panc1_rep12_40kb TCGApaad_wt_mutKRAS"

### TODO - added 20.06.2019
#"GSE118514_RWPE1_40kb TCGAprad_norm_prad"
#"GSE58752_liver_40kb TCGAlihc_norm_lihc"
#"GSE58752_liver_40kb TCGAlihc_wt_mutCTNNB1"
#"GSE75070_MCF-7_shNS_40kb TCGAbrca_lum_bas"


### TODO - added 24.06.2019
"K562_40kb TCGAlaml_wt_mutFLT3"

)


# Rscript create_coexpr_sortNoDup_otherTADfile.R "ENCSR346DCU_LNCaP"GSE73782_PC3_ICE_40kb TCGAprad_norm_prad"

for ds in "${all_TAD_files_ds[@]}"; do
    
    echo Rscript $scriptCoexpr $ds
    Rscript $scriptCoexpr $ds
done


###################################################################################################################################################
########## END ####################################################################################################################################
###################################################################################################################################################

echo "*** DONE"
end_time=$(date -R)    
echo $start_time
echo $end_time
exit 0

