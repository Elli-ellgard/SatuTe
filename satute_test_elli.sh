#!/usr/bin/env bash

set -e

path=$(pwd)
iqtree_path=/home/elgert/IQ-TREE/iqtree-2.2.2.4-Linux/bin/iqtree2

echo "-----------------------------------"
echo "Different tests for input handling"
echo "-----------------------------------"
echo ""

echo "-----------------------------------"
echo "TEST 1: empty directory"
cd $path/test/case_empty
rm -r *
python3 satute_cli.py  -iqtree $iqtree_path -dir $path/test/case_empty/

echo "-----------------------------------"
echo ""


echo "-----------------------------------"
echo "TEST 2a: only fasta alignment file"

echo "-----------------------------------"
echo ""


echo "-----------------------------------"
echo "TEST 2b: only phylip alignment file"

echo "-----------------------------------"
echo ""


echo "-----------------------------------"
echo  "TEST 3: directory with complete iqtree output"
echo "-----------------------------------"
echo ""




