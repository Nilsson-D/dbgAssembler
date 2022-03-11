#!/usr/bin/env python3

"""
Title: test_run_viruses.py
Created on Fri 11-03-2022
Author: Daniel Nilsson

Description:
    Tests the dbgAssembler on smaller genomes


List of functions:
    None
    
List of "non standard" modules:
    - networkx
    - matplotlib
    - pygraphviz
   
Procedure:
    1. reads in three small virus gemomes with varying size and one larger genome.
    2. calls the dbgAssembler
    3. outputs result to a test file <result_test_virus>

Usage:
    python test_run_viruses.py 
"""

from difflib import SequenceMatcher #to compare output and input sequence
import glob #to get the data
import sys # import required module
from pathlib import Path #use pathlib to create a result directory 
import time
# append the path of the parent directory so it can be run in the root dir
sys.path.append(".")

from Assembler import DbgSolver #get the assembler
from Assembler.Utilities.readfastafiles import readFasta_returnDict #to read in fasta file

genomes = glob.glob("data/viruses/*.fna")
dir_path = "test/test_viruses/"
output_file = "test/test_viruses/result_test_virus.fna"
output_test_log = "test/test_viruses/log_test_virus.txt"

#create directory from root if not existing
if not Path(dir_path).exists():
    Path(dir_path).mkdir()

k = 200

with open(output_file, "w") as output_fasta, open(output_test_log, "w") as output_log:

    for genome in genomes:
        
        fasta = readFasta_returnDict(genome)
        
        #get the sequences:
        sequences = "" 
        
        for seq in fasta.values():
            sequences += seq
        
        genome_size = len(sequences)
        
        start = time.time()
        solver = DbgSolver(seq, k) #call the class
        new_seq = solver.assemble() #get the reconstructed sequence 
        poportion_correct = 100*SequenceMatcher(None, new_seq, sequences).ratio()
        exe_time = round(time.time()-start,4)
        
        
        output_fasta.write(">%s\n" % genome)
        output_fasta.write("%s\n\n" % new_seq)
        
        output_log.write("Input genome: %s \n" % list(fasta.keys())[0])
        output_log.write("Genome size: %s \n" % genome_size)
        output_log.write("kmer size used: %s \n" % k)
        output_log.write("Poportion of sequence assembled correctly: %s%%\n" % poportion_correct)
        output_log.write("Execution time: %s \n\n" % exe_time)