#!/bin/bash

input=$1
mode=$2
cutMode=$3
vetoes=$4
pid=$5
MCcuts=$6

# modify to your own directory
output_dir="/eic/u/macink/EICreconOutputReader/output"

mkdir -p "${output_dir}"

output="${output_dir}/$(basename ${input})"
echo "Running analysis for file: ${input}"
echo "Output at [ ${output}_output.root ]"
    
# Run the diffractive VM analysis
root -b -q diffractive_vm_full_analysis.cxx\(\"${input}\",\"${mode}\",\"${cutMode}\",\"${vetoes}\",\"${pid}\",\"${MCcuts}\",\"${output}\"\)
echo "Finished processing file: ${input}"
echo
    
