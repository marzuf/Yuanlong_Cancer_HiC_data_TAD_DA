#!/usr/bin/bash

# ./all_prep_data_for_AUC_coexprDist_sameFam.sh

start_time=$(date -R)   

scriptPrepFam="prep_gene_families_TAD_data_otherTADfile.R"

scriptSameFam="create_sameFamily_sortNoDup_otherFamFile.R"


# DS launched 24.06.2019:
all_TAD_files_ds=(
# 1
"ENCSR079VIJ_G401_40kb"
# 2
"ENCSR312KHQ_SK-MEL-5_40kb"
# 3
"ENCSR346DCU_LNCaP_40kb"
# 4
"ENCSR401TBQ_Caki2_40kb"
# 5
"ENCSR444WCZ_A549_40kb"
# 6
"ENCSR489OCU_NCI-H460_40kb"
# 7
"ENCSR549MGQ_T47D_40kb"
# 8
"ENCSR862OGI_RPMI-7951_40kb"
# 9
"GSE105194_cerebellum_40kb"
# 10
"GSE105194_spinal_cord_40kb"
# 11
"GSE105318_DLD1_40kb"
# 12
"GSE105381_HepG2_40kb"
# 13
"GSE73782_PC3_40kb"
# 14
"Panc1_rep12_40kb"
# 15
"GSE118514_RWPE1_40kb"
# 16
"GSE58752_liver_40kb"
# 17
"GSE75070_MCF-7_shNS_40kb" 
# 18 # TODO - added 24.06.2019
"K562_40kb"
)



for ds in "${all_TAD_files_ds[@]}"; do

    echo Rscript $scriptPrepFam $ds
    Rscript $scriptPrepFam $ds
    
    echo Rscript $scriptSameFam $ds
    Rscript $scriptSameFam $ds
    
done


###################################################################################################################################################
########## END ####################################################################################################################################
###################################################################################################################################################

echo "*** DONE"
end_time=$(date -R)    
echo $start_time
echo $end_time
exit 0

