#/usr/bin/python3
#-----------------------------------------------------
import sys
import argparse
import numpy as np
import regex as re
import pandas as pd
import ete3
from ete3 import Tree
import csv
import os
import scipy
import scipy.linalg
import subprocess
#-----------------------------------------------------


def parse_input():
    """Parsing input arguments passed to the script"""
    #################################################
    print("="*100)
    print("")
    print("Satute - asymptotic test for branch saturation")
    print("")
    print("="*100)
    print("Author: ")
    print("Citation: ")
    print("="*100)
    print("")
    ##################################################

    parser = argparse.ArgumentParser(description='Satute')

    # parameters ---------------------------------------
    parser.add_argument("-o",help="output files prefix",default="satute_output_file",metavar="<file_name>")
    parser.add_argument("-dir",help="path to input directory",default="no directory",metavar="<file_name>")
    parser.add_argument("-iqtree",help="path to IQ-TREE",default="iqtree2",metavar="<file_name>")
    parser.add_argument("-nr", help="number of rate categories", type=int, default=1,metavar="<num>")

    #parsing arguments
    args = parser.parse_args(sys.argv[1:])
    #------------------------------------
    
    d=vars(args)
    file=open(args.o+".log","w")
    file.write("-"*100+'\n')
    file.write("Satute"+'\n')
    file.write("-"*100+'\n')
    file.write("Command line: ")
    file.write(str(sys.argv).replace("'","").replace(",","").replace("[","").replace("]","")+'\n')
    file.write("-"*100+'\n')
    file.write("Printing parameter values:"+'\n')
    file.write("-"*100+'\n')
    for i in d.keys():
        file.write(str(i)+":"+str(d[i])+'\n')
    file.write("-"*100+'\n')
    file.close()

    #-------------
    check_input(d["iqtree"], d["dir"], d["nr"], d["o"])


def check_input(iqtree,idir,nr,o):
    """ 
    Check if all parameters are correct and which files the input directory includes
    """
    print("Check available input")
    if os.path.isdir(idir):
        if os.path.isfile(idir + ".+.fasta"):
            print("fasta file exists")
        elif os.path.isfile(idir + ".+.phy"):
            print("phy file exists")
        
    else:
        sys.exit("INPUT_ERROR: input directory does not exist")


if __name__ == '__main__':
    parse_input()
    #check_input()
    #satute()
