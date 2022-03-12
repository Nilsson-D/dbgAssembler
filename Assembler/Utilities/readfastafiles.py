#!/usr/bin/env python3

"""
Title: readfastafiles.py
Created on Thur 09-03-2022
Author: Daniel Nilsson

Description:
    Utilities for reading in a fasta file
    * Create a dictionary where the key corresponds to the sequence id and the value corresponds to the sequence
    * convert fastq to fasta


List of functions:
    - readFastQ_returnDict_Fasta   
    - readFasta_returnDict
    
"""

import sys
from Assembler.Utilities.check_valid_sequence import check_valid_sequence
#script for reading fasta files


def readFasta_returnDict(input_fasta, seq_type = "DNA", allowN = True, isAligned = False):    
    """
    This function takes a fasta file (DNA) and write the ids and sequences
    to a dictionary, where the keys are the ids and the values are the
    corresponding values.
    
    """
    
    print(f"\nReading input fasta file: '{input_fasta}'") #message to inform the user
    
    with open(input_fasta, "r") as fasta_file: #open fasta in read mode
        
        fasta_dict = {} #A dictionary to hold the ids and their sequence
        sequence = "" #the variable to hold our complete sequence
        
        line_number = 0 #counts the number of lines
        
        #Checking if first line starts with > to ensure correct FASTA format as the header might only contain the characters A C G T
        if not fasta_file.readline().startswith('>'): #just an initial check
             print(f"""\nError: Specified file: '{input_fasta}' is not in FASTA format (DNA)
                   \nThe file should start with a fasta id (>id)""") #give the user an error message
             sys.exit() #exit 
        
        #Reset line pointer to the first line again
        fasta_file.seek(0) 
        
        
        for line in fasta_file: #iterate over each line in the fasta file
            line_number += 1 #for each line, increase the line count
            
            if line.startswith(">"): #check if the line starts with a '>' as this would be our id/header
               
              if line_number > 1: #checks if the line nr is larger than 1. If so we can add the header to the dictionary
                                  #as the first line will be a header and the line number is 1, this wont be executed.
                                  #thus, we are escaping the fact that we are referring to the 'header' before its assignment
                    fasta_dict[header] = sequence #add the header as a key to the dictionary with the sequence as the value 
               
              header = line.strip().replace(">", "") #take the id and strip from whitespaces and > 
              sequence = "" #reset the sequence after added
              
              if header in fasta_dict.keys(): #Check for duplicates in the fasta file
                  print(f"Error: fasta file: {input_fasta} contains duplicates") #print error message
                  sys.exit() #exit
                               
                
            else: #If the line does not start with a >. 
                  #This means that the line is a sequence
                line = line.strip()
                #Check if the  DNA sequence is valid an only has the characters -ACTGactg
                if not check_valid_sequence(line, seq_type, allowN, isAligned): #if any other character than -ACTGactg is found
                    print(f"""\nError: Specified file: '{input_fasta}' contains invalid characters in sequence with the id '{header}'
                          """) #print invalid DNA format
                    sys.exit() #exit 
               
                "".join([sequence, line.strip().upper()]) #append to the sequence, strip the new line and make the sequence to uppercase letters
                
        fasta_dict[header] = sequence #this is needed as the last header and sequence wont be added otherwise
                                      #because we are only adding a key value pair to the dict when a new header is found
                                      #and no new header is found at the end of the file
    
    
    for fasta_id, fasta_seq in fasta_dict.items(): #check so that every id has a sequence
        if not fasta_seq: #if a sequence is empty then warn the user and exit
            print(f"Error: missing sequence for fasta id '{fasta_id}'") #print a message
            sys.exit() #exit
            
    print(f"Done reading input fasta file: '{input_fasta}'")  #message to inform the user                           
    return fasta_dict #return the dictionary of sequence ids and sequences


def readFastQ_returnDict_Fasta(fastq_1, seq_type = "DNA", allowN = True):
    """
    Convert a fastq file to a fasta file, paired reads can be specified. The result is stored in a dictionary without the > as the first
    character in the header line
    """
    #create the headers for each line
    info = ['header', 'sequence', 'optional', 'quality']
    n = 4 #each id consists of 4 lines
    fq_1 = dict() #dictionary for forward reads

    with open(fastq_1, 'r') as f1: #start with the fastq with forward reads
        lines = []  #create a list to store the lines
        for line in f1:
            lines.append(line.rstrip()) #append each line to the list
            if len(lines) == n: #when the length of the list is 4
                record = {k: v for k, v in zip(info, lines)} #assign each line to the given headers in the info list
                header = "_".join(record.pop("header").split(" "))[1:] #remove the header and assign to a separate variable
                record.pop("optional") #remove lines with unecessary info
                record.pop("quality") #remove lines with unecessary info

                if not check_valid_sequence(record["sequence"], seq_type, allowN): #if any other character than -ACTGactg is found
                    print(f"""\nError: Specified file: '{fastq_1}' contains invalid characters in sequence with the id '{header}'
                          """) #print invalid DNA format
                    sys.exit() #exit

                fq_1[header] = record["sequence"]  #assign the sequence as the value to the header
                lines = [] #reset the lines list


        return fq_1
        
if __name__ == "__main__":

    pass