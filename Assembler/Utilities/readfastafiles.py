#!/usr/bin/env python3

"""
Title: readfastafiles.py
Created on Thur 09-03-2022
Author: Daniel Nilsson

Description:
    Utilities for reading in a fasta file
    * Either get a tuple containing two lists. One list containing the ids and the other containg the sequences
    * Or create a dictionary where the key corresponds to the sequence id and the value corresponds to the sequence


List of functions:
    - readFasta_returnTuple   
    - readFasta_returnDict
    
"""

import sys
from Assembler.Utilities.check_valid_sequence import check_valid_sequence
#script for reading fasta files

def readFasta_returnTuple(input_fasta, seq_type = "DNA", allowN = True, isAligned = False):
    """
    
    This function takes a FASTA file as input and stores the ids in one list and
    the sequence in one list, where id with index x will have the corresponding
    sequence with index x in the sequence list
    
    """
    print("\n\nReading the input fasta")  #print start of execution
    
    with open(input_fasta, "r") as in_fasta: #open the fasta file 
       
        #create variables to store the fasta ids and sequences in
        ids = [] #holds the fasta ids
        sequences = [] #holds the complete sequence
        seqFragments = [] #this is to store each sequence stretching over more than one line
        
        
        #Checking if first line starts with > to ensure correct FASTA format as the header might only contain the characters A C G T
        if not in_fasta.readline().startswith('>'): #just an initial check
             print(f"\nError: Specified file: '{input_fasta}' is not in DNA FASTA format") #give the user an error message
             sys.exit() #exit 
        
        #Reset line pointer to the first line again
        in_fasta.seek(0) 
        
        
        #parse the FASTA file     
        for fasta_line in in_fasta:
                              
            if fasta_line.startswith('>'): #check if the line starts with a '>'
                fasta_line = fasta_line[1:] #remove the > sign
                
                if fasta_line in ids: #Check for duplicates in the fasta file
                    print(f"Error: fasta file: {input_fasta} contains duplicates") #print error message
                    sys.exit() #exit
                    
                
                ids.append(fasta_line.rstrip("\n")) #Strip the new line from the end of the sequence id and
                                                    #append the id to the list of ids
                
                # found start of next sequence
                if seqFragments: #if seqFragment is not empty
                    
                    
                    sequence = ''.join(seqFragments).upper() #join the seqFragment to the sequence
                    sequences.append(sequence) #append the complete sequence to the sequence list
                    seqFragments = [] #reset the seqFragment for the next sequence
           
            else: #if the line does not start with a '>'
                
                fasta_line.strip()
                #Checking if the DNA sequence is valid
                if not check_valid_sequence(fasta_line, seq_type, allowN, isAligned): #if any other character than ACGT is found
                    print(f"\nError: Specified file: '{input_fasta}' contains invalid characters") #print not valid DNA format
                    sys.exit() #exit 
                    
                    
                # found more of existing sequence if the line does not start with a '>'
                seq = fasta_line.rstrip() # remove the new line character
                seqFragments.append(seq) #then append more of the sequence to seqFragment
            
            
        if seqFragments: #at the last line no new id that starts with '>' is found 
            # therefore we need to append the last sequence fragment to the sequence list
            sequence = ''.join(seqFragments).upper() #first join the seqFragments and the
            sequences.append(sequence) #append the complete sequence to the list of sequences
    
   
    for seq in sequences: #this is just a check to see if there is a empty string in the list of sequences
        if not seq: #if a empty string is found then an id does not have a sequence which creates a jump 
                    #in the index, not linking the index positions to each other in the id list and sequence list
                    #Therefore, this must be checked.
            print(f"""\nError: Number of sequence ids does not match the number of sequences 
              in specified file: '{input_fasta}'. """) #give the user an error message
            sys.exit() #exit
                    
        
    print("\n\nDone reading the input fasta") #print the end of execution
    return ids, sequences #return the id and sequence list




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
               
                sequence += line.strip().upper() #append to the sequence, strip the new line and make the sequence to uppercase letters
                
        fasta_dict[header] = sequence #this is needed as the last header and sequence wont be added otherwise
                                      #because we are only adding a key value pair to the dict when a new header is found
                                      #and no new header is found at the end of the file
    
    
    for fasta_id, fasta_seq in fasta_dict.items(): #check so that every id has a sequence
        if not fasta_seq: #if a sequence is empty then warn the user and exit
            print(f"Error: missing sequence for fasta id '{fasta_id}'") #print a message
            sys.exit() #exit
            
    print(f"Done reading input fasta file: '{input_fasta}'")  #message to inform the user                           
    return fasta_dict #return the dictionary of sequence ids and sequences




    
if __name__ == "__main__":
    """
    
    If the script is run as main.
    Here is the number of arguments provided checked and if they are valid.
    If everything looks fine then the functions are called  
    
    """ 
    pass