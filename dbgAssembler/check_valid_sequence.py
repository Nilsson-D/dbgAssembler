#!/usr/bin/env python3

"""
Title: check_valid_sequence.py
Created on Mon 07-03-2022
Author: Daniel Nilsson

Description:
    A script to store helper functions for the dbgAssembler regarding 
    checking if the input sequence is valid. 

List of functions:
    check_valid_sequence()
             
   
Procedure:
    1. get the sequence, the type, and if the sequence is from an alignment
    2. check whatever the sequence contains invalid characters
    3. return a boolean whatever the sequence is valid or not
    
Usage:
    - 
    used as module
"""

import re #allows searching for invalid characters 



def check_valid_sequence(fasta_line, seq_type = "DNA", allowN = True, isAligned = False):
    """
    Check whatever a dna, rna or protein sequence has the correct characters.
    It should be specified if the sequence is from an alignment 
    """
    seq_type = seq_type.upper() #seq_type to uppercase so this wont create a problem
    valid_fasta = True #create a boolean to return whatever the string/sequence has valid characters
    
    
    if seq_type == "DNA": #if DNA is specified
        #check whatever the sequence is from an alignment
        if isAligned:      
            if allowN:
                valid_fasta_chrs = r"[^NACGT-]+"
            else:
                valid_fasta_chrs = r"[^ACGT-]+"
        else:
            if allowN:
                valid_fasta_chrs = r"[^NACGT]+"
            else:
                valid_fasta_chrs = r"[^ACGT]+"
        
        if fasta_line == None:
            fasta_line = ""
        #search for all characters except valid_fasta_chrs    
        result = re.findall(valid_fasta_chrs, fasta_line.strip(), re.IGNORECASE) #ignore case when searching
        if result:
            valid_fasta = False #return true if any other chr is found
            
        
    elif seq_type == "PROTEIN":  #if protein is specified
    #check whatever the sequence is from an alignment
        if isAligned: 
            valid_fasta_chrs = r"[^ABCDEFGHIJKLMNOPQRSTUVWYZX*-]+"
        else:
            valid_fasta_chrs = r"[^ABCDEFGHIJKLMNOPQRSTUVWYZX*]+"
            
        #search for all characters except valid_fasta_chrs    
        result = re.findall(valid_fasta_chrs, fasta_line.strip(), re.IGNORECASE) #ignore case when searching
        if result:
            valid_fasta = False #return true if any other chr is found

    
    elif seq_type == "RNA": #if RNA is specified
        #check whatever the sequence is from an alignment
        if isAligned:         
            valid_fasta_chrs = r"[^ACGU-]+"
        else:
            valid_fasta_chrs = r"[^ACGU]+"
            
        #search for all characters except valid_fasta_chrs
        result = re.findall(valid_fasta_chrs, fasta_line.strip(), re.IGNORECASE) #ignore case when searching
        if result:
            valid_fasta = False #return true if any other chr is found
    
    return valid_fasta #return boolean
