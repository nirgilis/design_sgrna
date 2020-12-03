#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  3 10:53:04 2019

@author: pjvonk
"""
# this program requires the following things:
# biopython installed
# bowtie and bowtie_build in PATH (NOTE: this is bowtie, not bowtie2)
import os
import shutil

# Define variables
input_file = "./data/raw/input_sequences.fasta"
genome_assembly = "./data/raw/Schco3.assembly.fasta"
database = "Schco3"


# Confirm that bowtie is in the PATH. If it is, then extract the directory
# of the executable. This is needed by CCTOP. 
# NOTE: this is bowtie, not bowtie2!
bowtie_path = shutil.which("bowtie")
if not bowtie_path:
    raise Exception("the executable 'bowtie' is not in your PATH variable")
bowtie_location = os.path.dirname(bowtie_path)

# run bowtie
os.system("bowtie-build %s ./data/temp/%s" % (genome_assembly, database))

# run cctop
# TODO: change to format {}
cctop_command = "python ./bin/cctop_standalone/CCTop.py "\
                "--input {input_file} "\
                "--index ./data/temp/{database} "\
                "--bowtie {bowtie_location} "\
                "--pam NGG "\
                "--sgRNA5 NN "\
                "--totalMM 4 "\
                "--output ./data/temp/".format(
                input_file=input_file,
                database=database,
                bowtie_location=bowtie_location)
                
os.system(cctop_command)

# parse results and output the required files
# TODO: Write code to parse output of CCTop
def reverse_complement(protospacer):
    reverse_protospacer = ""
    
    for nt in protospacer.strip()[::-1]:
        if nt == "A":
            reverse_protospacer.append("T")
        elif nt == "T":
            reverse_protospacer.append("A")
        elif nt == "C":
            reverse_protospacer.append("G")
        elif nt == "G":
            reverse_protospacer.append("C")
        else:
            raise Exception("invalid character in protospacer. Please check your input sequences")
    return(reverse_protospacer)

def design_primers_sgRNA(protospacer):
    fw_addition = "TAATACGACTCACTATAG"
    rv_addition = "TTCTAGCTCTAAAAC"
    if protospacer[0] == "G":
        fw_primer = fw_addition + protospacer
    else:
        fw_primer = fw_addition + "G" + "protospacer"
    rc_protospacer = reverse_complement(protospacer)
    rv_primer = rv_addition + rc_protospacer
    return(fw_primer, rv_primer)