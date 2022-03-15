#!/usr/bin/env python3
"""
Title: test_run_viruses.py
Created on Fri 11-03-2022
Author: Daniel Nilsson

Description:
    Tests the dbgAssembler on smaller virus genomes

List of functions:
    None
    
List of "non standard" modules:
    - difflib

   
Procedure:
    1. reads in three small virus gemomes with varying size and one larger genome.
    2. calls the dbgAssembler
    3. calculates the hamming distance between the two sequences
    3. outputs result to a test file <result_test_virus>

Usage:
    python test_run_viruses.py 
"""

import glob #to get the data
from pathlib import Path #use pathlib to create a result directory 
import time #to time the runs

from Assembler import DbgSolver #get the assembler
from Assembler.Utilities.input_manager import InputManager #to read in fasta file


script_dir = Path( __file__ ).parent.absolute() #get the path to the dir where the test script is placed
genomes = glob.glob(f"{script_dir}/../test_data/*.fna") #fetch the genome files

#set the paths for the output
dir_path = f"{script_dir}/test_viruses"

#create directory from root if not existing
if not Path(dir_path).exists():
    Path(dir_path).mkdir()

ks = [x for x in range(3,34,3)] #create kmer size

def ham_dist(seq1, seq2):
    """
    calculates the hamming distance
    """
    if len(seq1) == len(seq2):
        hamming_distance = sum(nuc1 != nuc2 for nuc1,nuc2 in zip(seq1, seq2))
        
        portion_corr = 100*(len(seq1) - hamming_distance)/len(seq1)
    else:
        portion_corr = 100*(1-(len(seq2) - sum(nuc1 == nuc2 for nuc1,nuc2 in zip(seq1, seq2[:len(seq1)])))/len(seq2))
        
    return portion_corr


for genome in genomes: #for each genome (fasta file)
    name = Path(genome).name
    for k in ks: #run all sizes of k
            
        output_file = f"{dir_path}/k_{k}_{name}"       
        output_test_log = f"{dir_path}/log_k_{name}.txt"
        
        
        with open(output_file, "w") as output_fasta, open(output_test_log, "a") as output_log: #open the output files
        
            sequence = "".join(InputManager(genome).read_fasta())
            
                       
            
            start = time.time() #get the current time
            solver = DbgSolver(genome, k) #call the solver
            new_seq = solver.assemble() #get the reconstructed sequence 
            exe_time = round(time.time()-start,4) #get the execution time and poportion correctly assembled sequence
            
            name_id = genome.split("/")[-1] #format the header
            header = f"{name_id}_kmer_size:{k}"
            
            genome_size = len(sequence) #get the genome size for log info
                       
            #calc statistics            
            poportion_correct = ham_dist(new_seq, sequence) 
            
            #write the aseemblies to a file
            output_fasta.write(">Scaffold_1") 
            output_fasta.write("%s\n" % new_seq)
            
            
            #Write information to the log file
            output_log.write("Input genome: %s \n" % Path(genome).name) #first the input file
            output_log.write("Genome size: %s \n" % genome_size) #the size of the genome
            output_log.write("kmer size used: %s \n" % k) #kmer size used
            output_log.write("Poportion of sequence assembled correctly: %s%%\n" % poportion_correct) #poportion correct
            output_log.write("Execution time: %s \n\n" % exe_time) #execution time
            
            