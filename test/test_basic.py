#!/usr/bin/env python3

"""
Title: test_basic.py
Created on Wen 09-03-2022
Author: Daniel Nilsson

Description:
    A test script for dbgAssembler.py
    Generates 10 DNA sequences with length 10000 for de- and reconstruction.
    The DNA strings are used as input to the dbgAssembler class. 
    It will test the assembler on both unbiased sequences (all nucleotides are random)
    and on biased sequences (e.i higher GC-content).
    Prints to standard output whatever the reconstruction was successful or not.

List of functions:
    rnd_sequences
    
    
List of "non standard" modules:
    numpy - for making input probabilities for each nucleotide poosible
    SequenceMatcher
   
Procedure:
    1. Create two lists of dna sequences, one with gc-content 50% and onw with 90%
    2. Call the DbgGraph class on all the sequences and check whatever the sequences are reconstructed
    correctly.
    3. Print to standard output the 

Usage:
    python test_basic.py 
    
"""
from pathlib import Path #use pathlib to create a result directory
from difflib import SequenceMatcher #to compare output and input sequence
from numpy.random import choice #import choice for creating random sequences


from Assembler import DbgSolver


def rnd_sequences(times, length, bias=0.5):
    """
    Create >times> number of sequences with length >length> and with a gc-content of <bias>
    First create the probabilities and then loop <times> times and create the sequences
    by sampling a sequence of length <length> 
    """
    #create the probabilities. Divide by 2 as the prob must sum to 1
    a, t = (1 - bias)/2, (1 - bias)/2
    c, g = (bias)/2, (bias)/2
    
    probs = [a, t, c, g] #store the probabilities in a vector
    nucleotides = ["A", "T", "C", "G"] #create a list for the nucleotides
    sequences = list() #create a list to store the sequences
    
    
    #loop <times> times
    for _ in range(times):
        
        #create a sequence, sampling from nucleotides with a size <length> and
        #sample each nucleotide with the corresponding probability in p
        sequence = choice(nucleotides, size = length, p = probs)
        
        #join the sequence to create a string and then add to the sequences list
        sequences.append("".join(sequence))
    
    #return the list
    return sequences
    



if __name__ == "__main__":
    """
    If the script is run as main script
    """
length = 5000
times = 10
ks = [x for x in range(11,53,2)] #create kmer sizes
#generate unbiased sequences     
unbiased_seqs = rnd_sequences(times, length)

poportion_unbiased = dict() #to store the stats of the assembly  
result_unbiased_seqs = list() #to store the resulting sequences

#Print what sequences that will be run
print("Running sequences with GC-content: 50% \n\n")

for k in ks:
    poportion_unbiased[f"{k}"] = list()
    for run, seq in enumerate(unbiased_seqs, 1): #iterate over each sequence and record the run
        g = DbgSolver(seq, k) #call the class
        new_seq = g.assemble() #get the reconstructed sequence 
        
        result_unbiased_seqs.append(list(new_seq.values())[0]) # as the input is one entire sequence, we expect only one contig back
        
        poportion_correct = 100*SequenceMatcher(None, list(new_seq.values())[0], seq).ratio() #calculate the poportion correctly positioned nucs
        poportion_unbiased[f"{k}"].append(poportion_correct) #add to dict
    
          


#generate biased sequences   
gc=0.90

#generate biased sequences 
biased_seqs = rnd_sequences(times, length, gc)

poportion_biased = dict() #to store the stats of the assembly 
result_biased_seqs = list()#to store the resulting sequences

print(f"Running sequences with GC-content: {100*gc}%\n\n")

for k in ks:
    poportion_biased[f"{k}"] = list()
    for run, seq in enumerate(biased_seqs, 1): #iterate over each sequence and record the run
        g = DbgSolver(seq, k) #call the class
        new_seq = g.assemble() #get the reconstructed sequence 
        
        result_biased_seqs.append(list(new_seq.values())[0]) # as the input is one entire sequence, we expect only one contig back
        
        poportion_correct = 100*SequenceMatcher(None, list(new_seq.values())[0], seq).ratio() #calculate the poportion correctly positioned nucs
        poportion_biased[f"{k}"].append(poportion_correct) #add to dict

        
    
    
    
script_dir = Path( __file__ ).parent.absolute() #get the path to the dir where the test script is placed


#set the output paths
dir_path = f"{script_dir}/test_basic"
output_file = f"{dir_path}/result_test_basic.fna"
input_seqs = f"{dir_path}/inputSequences_test_basic.fna"
output_test_log = f"{dir_path}/log_test_basic.txt"

#create directory from root if not existing
if not Path(dir_path).exists():
    Path(dir_path).mkdir()
    

with open(input_seqs, "w") as fileObj: #write the input sequences to a file
    i = 1
    for unb_seq in unbiased_seqs: #First write the unbiased
        fileObj.write(f"Run nr {i}_unbiased\n")
        fileObj.write(f"{unb_seq}\n")
        i += 1
        
    i = 1    
    for b_seq in biased_seqs: #write the biased
        fileObj.write(f"Run nr {i}_biased\n")
        fileObj.write(f"{b_seq}\n")
        i += 1  
        

with open(output_file, "w") as fileObj: #write the resulting sequences to a file      
    i = 1
    for unb_seq in result_unbiased_seqs: #First write the unbiased
        fileObj.write(f"Run nr {i}_unbiased\n")
        fileObj.write(f"{unb_seq}\n")
        i += 1
        
    i = 1    
    for b_seq in result_biased_seqs: #write the biased
        fileObj.write(f"Run nr {i}_biased\n")
        fileObj.write(f"{b_seq}\n")
        i += 1     
    

with open(output_test_log, "w") as fileObj: #write a log file     

    i = 1
    #write the the unbiased
    fileObj.write(f"GC-content: 50%, length: {length}bp\n")
    for k, pop in poportion_unbiased.items():
        average = round(sum(pop)/len(pop),2) #calc average correct assembly
        fileObj.write(f"For kmer size: {k}, average poportion of correct assembly: {average}\n")
        i += 1

    i = 1
    #write the the biased
    fileObj.write(f"\n\nGC-content: 90%, length: {length}bp\n")
    for k, pop in poportion_biased.items():
        average = round(sum(pop)/len(pop),2) #calc average correct assembly
        fileObj.write(f"For kmer size: {k}, average poportion of correct assembly: {average}\n")
        i += 1
