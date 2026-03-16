# epic-AN-AC-2025-001
Diffractive phi ePIC analysis note.

This repository contains the code to reproduce the diffractive phi analysis and plots from the analysis note found at: 


## Run analysis local:

From eic-shell 

```./eic-shell --version 25.12.0-stable```

run `diffractive_vm_full_analysis.cxx(input_file.list,mode,cutMode,vetoes,pid,output.root)` where the options of the input parameters are defined in run_analysis.sh

For example:

''''diffractive_vm_full_analysis.cxx("subList_000.list","phi_coh","allCuts","allVetoes","withPID","testOutput.root")''''

