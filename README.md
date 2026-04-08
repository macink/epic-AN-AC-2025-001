# epic-AN-AC-2025-001
Diffractive phi ePIC analysis note.

This repository contains the code to reproduce the diffractive phi analysis and plots from the analysis note found here: [https://zenodo.org/records/17944278](https://zenodo.org/records/19439688)

The folder analysis_root_files contains all of root files resulting from the analysis of coherent diffractive phi and rho (with and without PID), incoherent phi, and DIS. These files can be ran in the `AnalysisNote_plots_diffractive_phi.C` macro to generate plots from the analysis note.

## Getting input file lists

The analysis code is set up to take in filelists as input. To geneate filelists, run
```./prep_file_list.sh``` from eic-shell. You must select the data directory and the corresponding basename and files you wish to generate.

## Run analysis locally:

From eic-shell (**Note:** if you do not have eic-shell installed, you can find instructions for installation here: https://eic.github.io/tutorial-setting-up-environment/)

```./eic-shell --version 25.12.0-stable```

run `diffractive_vm_full_analysis.cxx("input_file.list","mode","cutMode","vetoes","pid","MCcuts","output")` where the options of the input parameters are defined in `run_analysis.sh`

For example:

```root -b 'diffractive_vm_full_analysis.cxx("subList_000.list","phi_coh","allCuts","allVetoes","withPID","no","test")'```

This will create a root file, test_output.root from subList_000.list containing diffractive phi root files with all cuts, vetoes, and PID implemented in the analysis. This test_output.root file can be ran through the plot macro 'AnalysisNote_plots_diffractive_phi.C' to reproduce the analysis note plots.

## Run analysis on condor:

In `run_analysis.sh`, modify the config filename for the analysis you are performing (i.e. coherent phi production, incoherent phi production, coherent rho production, or DIS). Adjust the input parameters as necessary. Modify directory locations to correspond to your analysis directory. `job_run.sh` and `run_DiffractiveVM.sh` also need to have the directories updated to match your own.

To submit the jobs to condor, you must not be in eic-shell:

```./run_analysis 2```

This will create a list of output root files in the specified output folder.

## Reproduce analysis note plots

`AnalysisNote_plots_diffractive_phi.C` contains the macros to reproduce all of the plots in the analysis note. The input root file used to generate the plots must be updated accordingly.


