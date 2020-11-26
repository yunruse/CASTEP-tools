#!/bin/bash

rm -r $1
cp -r aimd_structures $1
#python3.8 supercell.py -fcvs4 $1/*.cell
python3.8 supercell.py -kfcvs6 -N $2 -p "$1/all.param" -r "$1/all.sh" $1/*.cell
