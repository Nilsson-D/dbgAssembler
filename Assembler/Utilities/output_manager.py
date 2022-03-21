#!/usr/bin/env python3
"""
Title: output_manager.py
Created on Fri 11-03-2022
Author: Daniel Nilsson

Description:
    This program handles the writing of the output from dbgAssembler to a file
    
List of functions:
    output_results

Procedure:
    1. gets the resulting output from dbgAssembler.py
    2. creates a fasta file fore the sequence
    3. creates a log file for the input file used and k value

"""


from pathlib import Path #used to check the path
from datetime import date #to get the date

from Assembler.Utilities import pathchecker #used for creating the directories

#fetch the date for the logfile and defualt directory
today = date.today()
current_date = today.strftime("%d_%m_%Y")

def output_results(input_arg, k, output_seq, output_file, directory):
    """
    Get the arguments from dbgAssembler and outputs to a fasta file and a log file
    input_arg is the input fasta from the user, k is the kmer size used, output_seq is the assembled sequence,
    output_file is the name of the output file and lastly, directory is the name of the output directory
    """
    
    logfile = f"log_run_{current_date}.txt"  #create a log file
   
    #create directories if necessary
    if not Path(directory).exists():
        pathchecker.createDirs(directory)
    

    #first create the file for the assembly
    output_file = directory + "/"+ output_file
    with open(output_file, "w") as output_result:
        for i,seq in enumerate(output_seq,1):
            output_result.write(f">Scaffold_{i}\n") #create a header
            output_result.write(f"{seq}\n") #Write out the sequence
            
 
    #then create a log file 
    logfile = directory + "/"+ logfile
    with open(logfile, "w") as output_log: 
        output_log.write("Input file: %s\n" % input_arg) #write the input file
        output_log.write("Size of kmer used: %s\n" % k) #write the kmer size

    print(f"Written to output directory: {directory}")       
            

if __name__ == "__main__":

    pass


