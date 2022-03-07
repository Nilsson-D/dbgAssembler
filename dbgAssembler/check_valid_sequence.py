#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 28 08:03:55 2022

@author: daniel
"""

import re

#script for reading fasta files


class Error(Exception):
    """Base class for other exceptions"""
    pass

class WrongFormat(Error):
    """Raised when the input is in wrong format"""
    pass




def check_valid_sequence(fasta_line, seq_type = "DNA", isAligned = False):
    seq_type = seq_type.upper()
    valid_fasta = True
    
    
    if seq_type == "DNA":
        if isAligned:         
            valid_fasta_chrs = r"[^ACGT-]+"
        else:
            valid_fasta_chrs = r"[^ACGT]+"
            
            
        result = re.findall(valid_fasta_chrs, fasta_line.strip(), re.IGNORECASE)
        if result:
            valid_fasta = False
            
        
    elif seq_type == "PROTEIN": 
        if isAligned: 
            valid_fasta_chrs = r"[^ABCDEFGHIJKLMNOPQRSTUVWYZX*-]+"
        else:
            valid_fasta_chrs = r"[^ABCDEFGHIJKLMNOPQRSTUVWYZX*]+"
            
        result = re.findall(valid_fasta_chrs, fasta_line.strip(), re.IGNORECASE)
        if result:
            valid_fasta = False

        
    elif seq_type == "RNA":
        if isAligned:         
            valid_fasta_chrs = r"[^ACGU-]+"
        else:
            valid_fasta_chrs = r"[^ACGU]+"
    
        result = re.findall(valid_fasta_chrs, fasta_line.strip(), re.IGNORECASE)
        if result:
            valid_fasta = False
    
    return valid_fasta
