#!/bin/bash

cd ${1}

# modify to your own directory
/eic/u/macink/eic-shell --version 25.12.0-stable << EOF
./run_DiffractiveVM.sh ${2} ${3} ${4} ${5} ${6} ${7}
exit
EOF