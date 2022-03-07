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
from check_valid_sequence import WrongFormat
import random
import numpy as np


class DbgSolver:
    """
    class to create kmers of length k using a string sequence.    
    """
    def __init__(self, sequence, k):
        """
        intialize variables
        """
        
        if not check_valid_sequence(sequence):
            raise WrongFormat("Error: DNA sequence contains ambiguous characters")
        
        
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
        node_index = dict()
        for i, node in enumerate(list(graph.keys())):
            node_index[node] = i
            
        
        for i in list(graph.keys()):
            row = node_index[i]
        
            for j in graph[i]:
                col = node_index[j]
                adjMatrix[row,col] +=1
                
        return adjMatrix, node_index

    def count_indegree_outdegree(self):
        """
        Count the number of outgoing and ingoing edges for each node
        """

        adj_matrix, nodes = self.create_adj_matrix()
        
        in_out = dict(zip(nodes.keys(),[None]*len(nodes.keys())))

        adj_matrix = np.array(adj_matrix)
        for key, i in zip(nodes.keys(), range(adj_matrix.shape[0])):    
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
        
        #function to traverse the graph
        def traversal(node):
            #while there still are unused bridges 
            while edge_count[node][1] != 0: 
                row = nodes[node] #get the row index from the kmer dictionary
                                
                connections = np.where(graph[row] > 0)   # find the indices where the values is over 0
                                                         # These are the unused edges
                                                         
                connected_node = random.choice(connections[0]) #pick a random node to traverse to
                
                #remove the the edge connecting the two nodes in the graph
                graph[row, connected_node] -= 1
                
                #decrease the out degree of the node by one
                edge_count[node][1] -= 1
                
                #To get the node in 
                for key, value in nodes.items():
                    if value == connected_node:
                        new_node = key
                        
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
    

if __name__ == "__main__":
    """
    If the script is run as main
    """
    test_seq = "ACTGACGTACGTACGTGTG"
    my_graph = DbgSolver(test_seq, 4) 
    my_graph.create_graph()
    
    new_seq = my_graph.get_sequence()
    
    print(new_seq)    
    print(test_seq)
    print(new_seq == test_seq)
