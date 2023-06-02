!/usr/bin/env bash

#set -e

path=$(pwd)
iqtree_path=/home/elgert/IQ-TREE/iqtree-2.2.2.4-Linux/bin/iqtree2

echo "-----------------------------------"
echo "Different tests for input handling"
echo "-----------------------------------"
echo ""

echo "-----------------------------------"
echo "TEST 1: empty directory"
echo""
    DIR=$path/test/case_empty
    if [ -d "$DIR" ]
    then
        if [ "$(ls -A $DIR)" ]; then
             echo "Take action $DIR is not Empty"
             rm -r $DIR/*
        else
            echo "$DIR is Empty"
	fi
    else
        echo "Directory $DIR not found."
    fi
    
    python3 satute_cli.py  -iqtree $iqtree_path -dir $DIR

echo "-----------------------------------"
echo ""


echo "-----------------------------------"
echo "TEST 2a: only fasta alignment file"
echo ""
    DIR=$path/test/case_fasta
    python3 satute_cli.py  -iqtree $iqtree_path -dir $DIR
    mv $DIR/example.fasta $path
    rm -r $DIR
    mkdir $DIR
    mv $path/example.fasta  $DIR

echo "-----------------------------------"
echo ""


echo "-----------------------------------"
echo "TEST 2b: only phylip alignment file"
echo ""

echo "-----------------------------------"
echo ""


echo "-----------------------------------"
echo  "TEST 3: directory with complete iqtree output"
echo ""

echo "-----------------------------------"
echo ""




