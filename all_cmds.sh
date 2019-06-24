#!/bin/bash

#all_datasets=(
#"ENCSR079VIJ_G401"
#"ENCSR105KFX_SK-N-DZ"
#"ENCSR312KHQ_SK-MEL-5"
#"ENCSR346DCU_LNCaP"
#"ENCSR401TBQ_Caki2"
#"ENCSR444WCZ_A549"
#"ENCSR489OCU_NCI-H460"
#"ENCSR549MGQ_T47D"
#"ENCSR834DXR_SK-N-MC"
#"ENCSR862OGI_RPMI-7951"
#"GSE105194_cerebellum"
#"GSE105194_spinal_cord"
#"GSE105318_DLD1"
#"GSE105381_HepG2"
#"GSE118514_RWPE1"
#"GSE58752_liver"
#"GSE73782_PC3"
#"GSE75070_MCF-7_shNS"
#"GSM1826481_SK-N-SH"
#"GSM2334832_RPMI-8226_HindIII"
#"GSM2334834_U266_MobI"
#"Panc1_rep12"
#)

all_datasets=(
"K562"
)

for ds in ${all_datasets[@]}; do

    echo $ds
    ./3_assign_genes.sh $ds

done 


#./run_pipeline.sh GSE105381_HepG2_40kb TCGAlihc_wt_mutCTNNB1   						# =>  done 22.01.2019




