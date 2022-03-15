#!/usr/bin/env python3
"""
Title: input_manager.py
Created on Mon 14-03-2022
Author: Daniel Nilsson

Description:
    Utilities for reading in a fasta or fastq file and checking for errors
    * Create a dictionary where the key corresponds to the sequence id and the value corresponds to the sequence
    * convert fastq to fasta


List of functions:
    -

List of classes: 
    -InputManager
        -functions:
            -correct_fasta
            -read_fasta

            
   
Procedure:
    1. The class takes a path to a file as input
    2. dependning on the uses, either fastq or fasta can be parsed and checked for errors.
    3. If reading the files. A generator object is returned.

"""

import re
from pathlib import Path

class InputManager():
    """
    Class for managing fasta and fastq files.
    It handles reading of the files and format checking
    """
    def __init__(self, file):
        
        if Path(file).exists():
            self.file = file
        else:
            raise Exception(f"Error: No such file: {file}")
        
        
    def correct_fasta(self, sequence, allowN=False):
        """
        Quick check whatever the fasta sequence is in the correct format. sequences containg Ns can be used
        if the user specied it.
        """
        
        if allowN:
            valid_nucs = r"[^NACGT]+" #create a pattern to (not) search for
            
        else:
            valid_nucs = r"[^ACGT]+" #create a pattern to (not) search for
        
        valid = True 

        if re.match(valid_nucs, sequence, re.IGNORECASE): #check if there are invalid characters
            valid = False
            
        return valid #will return true if the input is correct
        
    
    
    def read_fasta(self, allowN=False):
        """
        Read in a fasta file and output only the sequence
        
        """
        
        with open(self.file, "r") as fasta:
            header = None   #store the header for checking if it is valid
            sequence = list() #the variable to hold our complete sequence
            
            for line in fasta:
                if line.startswith(">"): #check if the line starts with a '>' as this would be our id/header      
                    if header: #this is fine for the first line as im only interested in the sequence
                        yield "".join(sequence) #only use the sequence as output
                            
                            
                    header = line.strip() #get the new header 
                    sequence = [] #reset the sequence
                    
                else:
                    line = line.strip().upper() #if not a header
                    
                    if not self.correct_fasta(line, allowN): #check if the line is a DNA sequence
                        raise Exception(f"Error: Found invalid characters in {self.file}")
                    
                    sequence.append(line) #add to the total sequence
                    
            if header: #last sequence wont be included in the loop
                yield "".join(sequence) #only use the sequence as output
                
                 

if __name__ == "__main__":

    pass  