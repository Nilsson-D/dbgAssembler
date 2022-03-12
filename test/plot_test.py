#!/usr/bin/env python3

"""
Title: plot_test.py
Created on Thur 10-03-2022
Author: Daniel Nilsson

Description:
    Visulize the de Bruijn graph created by the DbgSolver using a predefined DNA sequence.
    Networkx and pygraphviz are used later on to plot and the figure 


List of functions:
    None
    
List of "non standard" modules:
    - networkx
    - matplotlib
    - pygraphviz
   
Procedure:
    1. Create a Dna sequence
    2. Call DbgSolver (the assembler) using the sequence as input with a defined kmer size
    3. Create the de Bruijn graph and convert to a networkx object
    4. Set layout and convert to pygraphviz graph and save figure


Usage:
    python plot_test.py 
"""


# import requi9red module

from pathlib import Path #use pathlib to create a result directory

#Following pages are for visualizing the de Bruijn graph
import networkx as nx
import matplotlib.pyplot as plt
import pygraphviz
 

#Import the script
from Assembler import DbgSolver


k = 3 #set a value for the kmer size 
sequence = "CCAGTCAGCGGGATGTACGGGATT" #create a sequence to visualize

assembler = DbgSolver(sequence, k) #Create a DbgSolver object
dbg = assembler.create_graph() #get the de Bruijn graph

G = nx.MultiDiGraph(dbg) #convert to a networkx graph

#set layout 
G.graph['edge'] = {'splines': 'curved'}
G.graph = {"rankdir":"LR",
           "size":"10"}

#convert to a pygraphviz graph for saving the figure
A = nx.nx_agraph.to_agraph(G) 
A.graph_attr["label"] = f"Sequence: {sequence}, k = {k}"
A.layout('dot')                                                                


script_dir = Path( __file__ ).parent.absolute() #get the path to the dir where the test script is placed

#Create directory and output file
path = f"{script_dir}/test_plot"
filename_fig = "test_dbg_plot.png"

#create directory from root if not existing
if not Path(path).exists():
    Path(path).mkdir()

#draw and save the figure
A.draw(f'{path}/{filename_fig}')  

