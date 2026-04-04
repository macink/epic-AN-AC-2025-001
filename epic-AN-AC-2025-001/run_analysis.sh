#!/bin/bash

submit=$1

###############################################################
### Batch production for real data using file list          ###
###############################################################

if [ "$submit" -eq 2 ]; then

    # Uncomment for the corresponding file list needed
    
    configs=(sartre1.39-1.0_coherent_phi_eAu_bsat_10x100_ab)
    #configs=(BeAGLE1.03.02-1.1_phi_eAu_10x100_q2_1to10000_hiAcc_run)
    #configs=(sartre1.39-1.1_coherent_rho_eAu_bsat_10x100_q2_1to20_hiAcc)
    #configs=(BeAGLE1.03.02-1.0_DIS_eAu_10x100_q2_1to10_ab_run)

    pwd=$PWD

    for config in "${configs[@]}"; do
        echo "[i] Running config = $config"
        
        mode="unknown"

        if [[ "$config" == *"DIS"* ]]; then
            mode="DIS"
        elif [[ "$config" == *"rho"* ]]; then
            mode="rho"
        elif [[ "$config" == *"coherent_phi"* ]]; then
            mode="phi_coh"
        elif [[ "$config" == *"BeAGLE"* ]] && [[ "$config" == *"phi"* ]]; then
            mode="phi_incoh"
    fi
    echo "[i] Auto‑detected mode = $mode"

    cutMode="allCuts"
    #cutMode="noCuts"
    vetoes="allVetoes"
    #vetoes="noVetoes"
    pid="withPID"
    #pid="noPID"
    MCcuts="no"
    #MCcuts="yes"

    # modify to your own directory 
    odir=/eic/u/macink/EICreconOutputReader/analysis/$config
    logdir=$odir/log

    mkdir -pv "$odir"
    rm -rf "$odir"/*
    mkdir -pv "$logdir"

    # modify to your own directory
    echo "Checking files: " $(ls /eic/u/macink/EICreconOutputReader/$config/subList*)

    executable=job_run.sh
    cp -v ${executable} $odir/.
    cp -v analysis $odir/.
    cp -v pleaseIncludeMe.h $odir/.
    cp -v diffractive_vm_full_analysis.cxx $odir/.
    cp -v run_DiffractiveVM.sh $odir/.

    condor_file=CondorFile_$config
    echo "" > ${condor_file}
    echo "Universe    = vanilla" >> ${condor_file}
    echo "Executable  = ${odir}/${executable}" >> ${condor_file}
    echo "GetEnv      = True" >> ${condor_file}

    # modify to your own directory
    files=$(ls /eic/u/macink/EICreconOutputReader/$config/subList*)

    for file in $files; do
        listNum=$(basename ${file} | sed "s/.list//g" | cut -f 2 -d _)
        LogFile=${logdir}/log_${listNum}.out
        ErrFile=${logdir}/log_${listNum}.err
        OutFile=${odir}/output_${listNum}.root

        echo "" >> ${condor_file}
        echo "Output    = ${LogFile}" >> ${condor_file}
        echo "Error     = ${ErrFile}" >> ${condor_file}
        echo "Arguments = ${odir} ${file} ${mode} ${cutMode} ${vetoes} ${pid} ${MCcuts} ${OutFile}" >> ${condor_file}
        echo "Queue" >> ${condor_file}
    done

    mv ${condor_file} $odir/.
    cd $odir
    condor_submit ${condor_file}
    cd $pwd

done
fi
