#!/usr/bin/env python3

"""
Title: test.py
Created on Wen 09-03-2022
Author: Daniel Nilsson

Description:
    A test script for dbgAssembler.py
    Generates 10 DNA sequences with length 1000 for de- and reconstruction.
    The DNA strings are used as input to the dbgAssembler class. 
    It will test the assembler on both unbiased sequences (all nucleotides are random)
    and on biased sequences (e.i higher GC-content).
    Prints to standard output whatever the reconstruction was successful or not.

List of functions:
    rnd_sequences
    biased_sequences
    
List of "non standard" modules:
    numpy - for making input probabilities for each nucleotide poosible
   
Procedure:
    1. Create two lists of dna sequences, one with gc-content 50% and onw with 90%
    2. Call the DbgGraph class on all the sequences and check whatever the sequences are reconstructed
    correctly.
    3. Print to standard output the 

Usage:
    python test.py 
"""


# import required module
import sys
  
# append the path of the
# parent directory so it can be run in the root dir
sys.path.append(".")


from numpy.random import choice #import choice for creating random sequences
from Assembler import DbgSolver


def rnd_sequences(times = 10, length = 1000, bias=0.5):
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


#generate unbiased sequences     
seqs = rnd_sequences()

#create variables for counting the failed and successful assemblies
correct = 0
failure = 0
#Print what sequences that will be run
print("Running sequences with GC-content: 50%")
for run, seq in enumerate(seqs, 1): #iterate over each sequence and record the run
    g = DbgSolver(seq) #call the class
    new_seq = g.assemble() #get the reconstructed sequence 
    
    #print some information for each run
    print(f"Run nr: {run}")
    print(f"Input sequence: {seq} \n\nOutput sequence: {new_seq}")
    print(f"Correct assembly: {new_seq == seq}\n-----------------------------------------------------------\n")
    
    #count the failed and successful assemblies 
    if new_seq == seq:
        correct += 1 
    else:
        failure += 1
        
      

print("\n\n\n\n\n\n")


#generate biased sequences   
gc=0.90

#generate biased sequences 
biased_seqs = rnd_sequences(bias = gc)

#create variables for counting the failed and successful assemblies
biased_correct = 0
biased_failure = 0
#Print what sequences that will be run
print(f"Running sequences with GC-content: {100*gc}%")
for run, seq in enumerate(biased_seqs, 1): #iterate over each sequence and record the run
    g = DbgSolver(seq) #call the class
    new_seq = g.assemble() #get the reconstructed sequence 
    
    #print some information for each run
    print(f"Run nr: {run}")
    print(f"Input sequence: {seq} \n\nOutput sequence: {new_seq}")
    print(f"Correct assembly: {new_seq == seq}\n-----------------------------------------------------------\n")
    
    #count the failed and successful assemblies 
    if new_seq == seq:
        biased_correct += 1 
    else:
        biased_failure += 1
        
 
#Print overall statistics for biased and unbiased sequences        
print(f"Assembly statistics for gc-content: 50% \nCorrect assemblies: {correct}\nFailed assemblies: {failure}\n\n")
print(f"Assembly statistics for gc-content: {100*gc}% \nCorrect assemblies: {biased_correct}\nFailed assemblies: {biased_failure}")
