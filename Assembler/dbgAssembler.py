#!/usr/bin/env python3

"""
Title: dbgAssembler.py
Created on Mon 07-03-2022
Author: Daniel Nilsson

Description:
    This program is the main script for calling the assembler. 
    It parses the command line for inputs and calls the DbgSolver to create the assembly.

List of functions:
    option_handler()

List of "non standard" modules:
    -
   
Procedure:
    1. calls the dbg.py which creates the de bruijn graph and ouputs the sequence
    2. outputs the assembled sequence to the function in output_manager to create
       output files.
Usage:
    python dbgAssembler.py -i <input_file> -k <kmer_size> [optional] -o <output_file> [optional]

"""


import sys 
from pathlib import Path #for checking the existence path
import argparse #enables user input
from datetime import date #to get the date for creating the deafult names for the output files

from Assembler.Utilities import pathchecker #used for creating the directories
from Assembler.Utilities.output_manager import output_results #to write result to a file
from Assembler.de_bruijn_graph import DbgSolver #import the script that deconstructs and reconstructs the sequence


def option_handler():
    """
    When called as main, get the input from the user
    """
    today = date.today()
    current_date = today.strftime("%d_%m_%Y")
    
    #create default parameters
    output = f"dbgAssembler_run_{current_date}.fna"
    directory = f"dbgAssembler_{current_date}"
    k = 31
    
    
    
    #create a parser for the command line
    parser = argparse.ArgumentParser(usage="""%(prog)s -i <input_file> -k <kmer_size> [optional] -o <output_file> [optional] \nType -h/--help for the help message""",
                                description="This program takes an one-line fasta file (DNA) as input and breaks the sequence into kmers of size k. Then reassembles the string using a de Bruijn graph based approach")
   

    #Here are all the parameters that the user can specify           
    parser.add_argument("-i", metavar='<input file>',  help="path to fasta file", required=True)
    
    parser.add_argument("-k",  metavar='<kmer size>', type=int, help="kmer size (default: 31, max: 251)", default=k)
    
    parser.add_argument("-o", metavar='<output file>', type=str,  help="name of output file, default: dbgAssembler_run{current_date}.fna", default=output)
    
    parser.add_argument("-d", metavar='<directory>', type=str,  help="name of output directory to create, default: dbgAssembler_{current_date}", default=directory)
    
    parser.add_argument("-n", metavar='<y/n>',  help="if y, allow Ns in sequence", default=False)
    
    #assign the parsed input to a variable
    args = parser.parse_args()
    
    #assign each input to a new variable
    input_file = args.i
    k = args.k 
    allowN = args.n
    
    #variables for output files
    output = args.o 
    directory = args.d
        
    
    if allowN.lower() == "y":
        allowN = True 
    elif allowN.lower() == "n":
        allowN = False
    else:
        raise Exception("Error: Invalid input for oprtion -n")
    
    if not pathchecker.check_path_overwrite(directory):
        print("Exiting...")
        sys.exit()  
        
    
    if Path(input_file).is_file():  #check whatever the input file is an actual file..   
        
        #call the DgbSolver for t
        solver = DbgSolver(input_file, k, allowN)
 
                                      
    else: #raise an error in the 
        raise Exception(f"Error: {input_file} does not exist")
                      
    
    #get the reconstructed sequence 
    new_seq = solver.assemble()
       
    #write to output file     
    output_results(input_file, solver.k, new_seq, output, directory) #call the function to write to output files


    

if __name__ == "__main__":
    option_handler()