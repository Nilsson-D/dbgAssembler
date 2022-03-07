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

import random
from check_valid_sequence import check_valid_sequence
from check_valid_sequence import WrongFormat

class DbgGraph:
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
        self.graph = dict()

        
    def create_uniq_kmers(self):
        """
        Split the sequence in kmers of size k and add to the dictionary self.edges. 
        Associate each kmer with its repeat count
        """
        edges = dict()
        
        #get length of the sequence
        sequence_length = len(self.sequence)
        
        #iterate over the length to get the start index and the length of the sequence 
        #starting at position kto get the end index for slicing
        for start, end in zip(range(sequence_length), range(self.k, sequence_length+1)):
            #slice the sequence in kmers using the indices:
            kmer = self.sequence[start:end]
            
            if kmer not in edges.keys():
                edges[kmer] = 1
            else:
                edges[kmer] += 1
        
        print(f"Creating edges: {edges}")
        return edges
        
    
    def create_graph(self, edges):
        """
        Create the node by splitting up the edges in pairs of k-1mers
        """
        
        for node in edges.keys():
            left, right =  node[:-1], node[1:]
            
            if left not in self.graph:
                self.graph[left] = set()
            
            if right not in self.graph:
                self.graph[right] = set()
                
            #As left and right are from the same kmer, this will always be true in terms of 
            #overlap
            self.graph[left].add(right)
            
        print(f"Creating graph with nodes: {self.graph}")
        return self.graph       
    
    
    def count_edges(self, edges):
        """
        Count the number of outgoing and ingoing edges for each node
        """
        
        nodes = dict()
        
        for node in self.graph.keys():
            
            nodes[node] = {"in": 0, "out": 0} 
            
            nodes[node]["out"] = len(self.graph[node])
            
            #get the in edges by checking the keys in which 
            #the node is the value for
            for connected_nodes in self.graph.values():
                if node in connected_nodes:
                    nodes[node]["in"] += 1  


        for edge, count in edges.items():
            left, right =  edge[:-1], edge[1:]
            if count > 1:
                nodes[left]["out"] += count -1
                nodes[left]["out"] += count -1
                
                nodes[right]["in"] += count -1
                nodes[right]["in"] += count -1
                
        return nodes
    
    
    def is_eulerian(self, edges):
        """
        Check if an eulerian path exists
        The conditions are that we must have one node with:
            * outdegree - indegree = 1 and 
            * indegree - outdegree = 1
        
        The node with one more outdegree than indgree will be our start node when searching
        """
        exists = False
        start = set()
        end = set()
        nodes = self.count_edges(edges)
        
        for node, degrees in nodes.items(): 
            diff = degrees["out"] - degrees["in"] 
            if diff == 1:
                start.add(node)
            elif diff == -1:
                end.add(node)
                
        if len(start) == 1 and len(end) == 1:
            exists = True
           
        if exists:
            return "".join(start)
        else:
            return exists
    
    
    def dfs(self, edges):
        """
        Use a depth-first search algorithm to traverse each edge to a new node
        """
        
        node = self.is_eulerian(edges)  
        edge_count = self.count_edges(edges)
        path = list()
        
        
        def traversal(node):
            
            #choose new node at random
            while edge_count[node]["out"] != 0: 
               
               connected_node = random.choice(list(self.graph[node]))
               self.graph[node].remove(connected_node)
               edge_count[node]["out"]  -= 1
               traversal(connected_node)
               
            path.append(node)                 
    
        traversal(node)
        
        return path
            
    
    def print_sequence(self):
        sequence = ""
        first = True
        
        edges = self.create_uniq_kmers()
        
        self.create_graph(edges)
        self.is_eulerian(edges)
        path = self.dfs(edges)
        
        for node in path[::-1]:
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
    test_seq = "ACTGACGTACGTACGTGTGAAGCTAGTCGCGCATTACGGGTTGAAACGTTGGTT"
    my_graph = DbgGraph(test_seq, 3) 
    
    edges = my_graph.create_uniq_kmers()
    my_graph.create_graph(edges)
    my_graph.count_edges(edges)
    
    new_seq = my_graph.print_sequence()
    
    print(new_seq)    
    print(test_seq)
    print(new_seq == test_seq)
