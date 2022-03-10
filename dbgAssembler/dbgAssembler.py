#!/usr/bin/env python3

"""
Title: dbgAssembler.py
Created on Mon 07-03-2022
Author: Daniel Nilsson

Description:

List of functions:

List of "non standard" modules:
   
Procedure:
    1. Create kmers of specified size
    2. Create de bruijn graph
        * edges being the kmer
        * nodes k-1mers of the kmers
        * only uniq kmers will be used
        * repeated kmers are counted
    4. Check Eulerian path
    5. Use Depth-first search algorithm to traverse the edges and nodes
        * traverse the same edge only once
    6. When a dead end is found, put the node at the last position for the construction of
       the new sequence, then backtrack until a node with an unused edge is found and start over the path search

Usage:
    python dbgAssembler.py -i <sequence> -k <kmer_size>

"""


from check_valid_sequence import check_valid_sequence
from readfastafiles import readFasta_returnDict


from pathlib import Path
import random
import numpy as np
import argparse

class DbgSolver:
    """
    class to create kmers of length k using a string sequence.    
    """
    def __init__(self, sequence, k):
        """
        intialize variables
        """
        
        if not check_valid_sequence(sequence):
            raise Exception("Error: DNA sequence contains ambiguous characters")
        
        
        self.sequence = sequence
        self.k = k
        
    def create_kmers(self):
        """
        Split the sequence in kmers of size k and add to the dictionary self.edges. 
        Associate each kmer with its repeat count
        
        This creates |Text| k + 1 edges
        """
        edges = list()

        #get length of the sequence
        sequence_length = len(self.sequence)

        #iterate over the length to get the start index and the length of the sequence 
        #starting at position kto get the end index for slicing
        for start, end in zip(range(sequence_length), range(self.k, sequence_length+1)):
            #slice the sequence in kmers using the indices:
            kmer = self.sequence[start:end]

            edges.append(kmer)

        return edges

    def create_graph(self):
        """
        Create the nodes by splitting up the edges in pairs of k-1mers
        This is an initial graph which do not have uniq k-1mers but will be used when creating the 
        final de Bruijn graph as an adjacency matrix
        """
        edges = self.create_kmers()
        graph = dict()
        
        
        for node in edges:
            left, right =  node[:-1], node[1:]
            
            if left not in graph:
                graph[left] = list()
            
            if right not in graph:
                graph[right] = list()
                
            #As left and right are from the same kmer, this will always be true in terms of 
            #overlap
            graph[left].append(right)
            
        return graph
    
    def create_adj_matrix(self):
        """
        Represent the de Bruijn graph as an adjacency matrix 
        Here we merge repeated nodes to an unique node 
        """
        graph = self.create_graph() #get the intial graph

        adjMatrix = np.zeros([len(graph.keys()),len(graph.keys())], dtype=int)
        node_index = list()
        for i in list(graph.keys()):
            node_index.append(i)
            
        
        for i in list(graph.keys()):
            row = node_index.index(i)
        
            for j in graph[i]:
                col = node_index.index(j)
                adjMatrix[row,col] +=1
                
        return adjMatrix, node_index

    def count_indegree_outdegree(self):
        """
        Count the number of outgoing and ingoing edges for each node
        """

        adj_matrix, nodes = self.create_adj_matrix()
        
        in_out = dict(zip(nodes, [None]*len(nodes)))

        adj_matrix = np.array(adj_matrix)
        for key, i in zip(nodes, range(adj_matrix.shape[0])):    
            in_degrees = sum(adj_matrix.T[i])
            out_degrees = sum(adj_matrix[i])
            in_out[key] = [in_degrees, out_degrees]

        return in_out
    
    def get_nodes(self):
        """
        The problem is to reconstruct a given string. This is not a circular path.
        Thus, we need to find the start and end node which should be the unbalanced nodes
        """
        indegree_outdegree = self.count_indegree_outdegree()
        
        for node in indegree_outdegree.items():
            in_deg, out_deg = node[1][0], node[1][1]
            
            if in_deg < out_deg:
                start = node[0]                   
                    
        return start
    
    def dfs(self):
        """
        Use a depth-first search algorithm to traverse each edge to a new node
        When a stop is reached, put it in the list and go back to a node that has edges 
        not yet travered
        """
        
        node = self.get_nodes()  
        edge_count = self.count_indegree_outdegree()
        graph, nodes = self.create_adj_matrix()
        
        path = list()
        
        def traversal(node):
            #while there still are unused bridges 
            while edge_count[node][1] != 0: 
                
                
                row = nodes.index(node) #get the row index from the kmer dictionary                            
                connections = np.where(graph[row] > 0)   # find the indices where the values is larger than 0
                                                         # These are the unused edges
                                                         
                connected_node = random.choice(connections[0]) #pick a random node to traverse to
                
                #remove the the edge connecting the two nodes in the graph
                graph[row, connected_node] -= 1
                
                #decrease the out degree of the node by one
                edge_count[node][1] -= 1
                
    
                new_node = nodes[connected_node]
                            
                traversal(new_node)
                   
            path.append(node)   

        traversal(node)    
        
    

        
        return path[::-1]      
    
    
    def get_sequence(self):
        sequence = ""
        first = True
        
        path = self.dfs()
        
        for node in path:
            if first:
                sequence += node
                first = False
            else:    
                sequence += node[-1]
                
        return(sequence)
    

def main():
    """
    When called as main, get the input of 
    """
    
    #create a parser for the command line
    parser = argparse.ArgumentParser(usage="""%(prog)s python dbgAssembler.py <sequence> -k <kmer_size> [optional] \nType -h/--help for the help message""",
                                 description="This program takes a dna string as input and breaks it to kmers of size k. Then reassembles the string using a De Bruijn graph styled manner")
    
    #create arguments for inputs: sequence and k
    parser.add_argument("-i", metavar='<input_sequence>', type=str,  help="a DNA string", default=None)
        
    #create arguments for inputs: sequence and k
    parser.add_argument("-f", metavar='<fasta_file>',  help="path to fasta file", default=None)
    
    parser.add_argument("-k",  metavar='<kmer_size>', type=int, help="kmer size (default 1/3 of the sequence length, max: 256)")
    
    #assign the parsed input to a variable
    args = parser.parse_args()
    
    #assign each input to a new variable
    sequence = args.i
    fasta_file = args.f
    k = args.k 
    
    if not sequence and not fasta_file:
        raise Exception("Only one of -i and -f can be specified at a time")
    
    
    if fasta_file != None :
        if  Path(fasta_file).is_file():    
              adv_test = readFasta_returnDict(fasta_file)
              
              #get the sequences:
              sequences = "" 
              
              for seq in adv_test.values():
                  sequences += seq
              #call the DgbSolver
              if k: 
                  graph = DbgSolver(sequences, k)
              else:
                 graph = DbgSolver(sequences)  
           
        else:
            raise Exception(f"Error: {fasta_file} is not a file or does not exist")
        
    #check if k is specified and create a DbgGraph object with the k or without if not specified 
    elif k: 
        graph = DbgSolver(sequence, k)
    else:
       graph = DbgSolver(sequence) 
    
    
    #get the reconstructed sequence and print the sequence to standard output
    new_seq = graph.get_sequence()
   # print(f"\nReassembly: {new_seq}")
    
    #Tell the user if the reassembly was succesful or not by comparing 
    #the input sequence with the reconstructed
    print(f"Successful assembly: {sequence == new_seq}")
    
    





if __name__ == "__main__":
    """
    If the script is run as main
    """
    main()
    
  #  test_seq = "ACTGACGTACGTACGTGTG"
   # my_graph = DbgSolver(test_seq, 4) 
    #my_graph.create_graph()
   # my_graph.get_nodes()
    
   # adj_m, ind_m = my_graph.create_adj_matrix()
   # adj_m
   # ind_m
    
    
   # new_seq = my_graph.get_sequence()
    
   # print(new_seq)    
   # print(test_seq)
   # print(new_seq == test_seq)
