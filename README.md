# epic-AN-AC-2025-001
Diffractive phi ePIC analysis note.

This repository contains the code to reproduce the diffractive phi analysis and plots from the analysis note found at: 

## Getting input file lists

The analysis code is set up to take in filelists as input. To geneate filelists, run
```./prep_file_list.sh``` from eic-shell. You must select the data directory and the corresponding basename and files you wish to generate.

## Run analysis locally:

From eic-shell 

```./eic-shell --version 25.12.0-stable```

run `diffractive_vm_full_analysis.cxx(input_file.list,mode,cutMode,vetoes,pid,output)` where the options of the input parameters are defined in `run_analysis.sh`

For example:

```root -b 'diffractive_vm_full_analysis.cxx("subList_000.list","phi_coh","allCuts","allVetoes","withPID","test")'```

This will create a root file, test_output.root that can be ran through the plot macros to reproduce the analysis note plots.

## Run analysis on condor:

In `run_analysis.sh`, modify the config filename for the analysis you are performing (i.e. coherent phi production, incoherent phi production, coherent rho production, or DIS). Adjust the input parameters as necessary. Modify directory locations to correspond to your analysis directory. `job_run.sh` and `run_DiffractiveVM.sh` also need to have the directories updated to match your own.

To submit the jobs to condor, you must not be in eic-shell:

```./run_analysis 2```

This will create a list of output root files in the specified output folder.

## Reproduce analysis note plots

`AnalysisNote_plots_diffractive_phi.C` contains the macros to reproduce all of the plots in the analysis note. The input root file used to generate the plots must be updated accordingly.


