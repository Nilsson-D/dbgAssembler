#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Title: output_manager.py
Created on Fri 11-03-2022
Author: Daniel Nilsson

Description:
    This program handles the writing of the output from dbgAssembler to a file
    
List of functions:
    output_results

"""

import sys
from pathlib import Path #used to check the path
from datetime import date #to get the date

from Assembler.Utilities import pathchecker #used for creating the directories

#featch the date for the logfile and defualt directory
today = date.today()
current_date = today.strftime("%d_%m_%Y")

def output_results(input_arg, k, output_seq, poportion, 
                   output_file, directory, onlySeq = False):
    """
    Get the arguments from dbgAssembler
    """
    
    logfile = f"log_run_{current_date}.txt"  #create a log file
   
    #create directories if necessary
    if not Path(directory).exists():
        pathchecker.createDirs(directory)
    
    else: #check if the directory name already exists and if the user wants to overwrite it
       if not pathchecker.check_path_overwrite(directory):
         print("Exiting...")
         sys.exit()  

    #first create the file for the assembly
    output_file = directory + "/"+ output_file
    with open(output_file, "w") as output_result:
        output_result.write(">Scaffold_1\n") #create a header
        output_result.write("%s\n" % output_seq) #Write out the sequence
            
    #check whatever it is a sequence or a fasta file
    if not onlySeq:   
        #then create a log file 
        logfile = directory + "/"+ logfile
        with open(logfile, "w") as output_log: 
            output_log.write("Input file: %s\n" % input_arg) #write the input file
            output_log.write("Size of kmer used: %s\n" % k) #write the kmer size
            output_log.write("Poportion of sequence assembled correctly: %s%%\n" % poportion) #and the porportion that is correct
    
    else:       
        #then create a log file 
        logfile = directory + "/"+ logfile
        with open(logfile, "w") as output_log: 
            output_log.write("Input sequence: %s\n" % input_arg) #write the input file
            output_log.write("Size of kmer used: %s\n" % k) #write the kmer size
            output_log.write("Poportion of sequence assembled correctly: %s%%\n" % poportion) #and the porportion that is correct       
            
            




