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
    def __init__(self, sequence, k=None):
        """
        intialize variables
        """
        
        if not check_valid_sequence(sequence):
            raise WrongFormat("Error: DNA sequence contains ambiguous characters")
        
        self.sequence = sequence
        
        if not k:
            k = int(len(self.sequence)/3)
        self.k = k
        
    def __create_kmers(self):
        """
        Split the sequence in kmers of size k and add to the dictionary edges. 
        Associate each kmer with its repeat count
        
        Nr of kmers = L-k+1, L is the length of the input sequence, k is the kmer size
        """
        edges = dict()
        
        #get length of the sequence
        sequence_length = len(self.sequence)
        
        #iterate over the length to get the start index and the length of the sequence 
        #starting at position k to get the end index for slicing
        for start, end in zip(range(sequence_length), range(self.k, sequence_length+1)):
            #slice the sequence in kmers using the indices:
            kmer = self.sequence[start:end]
            
            if kmer not in edges.keys():
                edges[kmer] = 1
            else:
                edges[kmer] += 1
        
        return edges
        
    
    def create_graph(self):
        """
        Create the node by splitting up the edges in pairs of k-1mers
        """
        edges = self.__create_kmers()
        graph = dict()
        
        for node, repeat in edges.items():

            left, right =  node[:-1], node[1:]
            
            if left not in graph:
                graph[left] = dict()
            
            if right not in graph:
                graph[right] = dict()
                
            #As left and right are from the same kmer, this will always be true in terms of 
            #overlap
            if right not in graph[left].keys():              
                graph[left].update({right:repeat})
                
        self.graph = graph
        return graph       
    
    
    def count_edges(self):
        """
        Count the number of outgoing and ingoing edges for each node
        """
        edges = self.__create_kmers()
        nodes = dict() #create a dictionary to store each node with its in and out degreees
        
        for node in self.graph.keys(): #loop over each node
            
            nodes[node] = {"in": 0, "out": 0} #create nested dictionaries to store in and out degrees
            
            #count the out degree
            for connected_node in self.graph[node]:
                edges = self.graph[node][connected_node] #The number of repeats will be the number of 
                                                         #edges for one connected node
                nodes[node]["out"] += edges #add the edges to one connected node to get the out degree of the current node
            
            #get the in edges by checking the keys in which 
            #the node is the value for
            for connected_nodes in self.graph.values():
                if node in connected_nodes: #if the node is present as a value for another node, then
                                            #the number of out edges for the other node will 
                                            #be the number of in edges for the current node
                    nodes[node]["in"] += connected_nodes[node] #Sum all in edges to get the in degree of the current node
              
        return nodes
    
    
    def get_nodes(self):
        """
            If an eulerian trail is present, then we have one node with:
                
            * outdegree - indegree = 1 (start) and one with
            * indegree - outdegree = 1 (end)
        
            If instead we have an eulerian cycle, all nodes are balanced:  outdegree = indegree
            Then start at a random node
        """
        start = set()
        end = set()
        nodes = self.count_edges()
        
        for node, degrees in nodes.items(): 
            diff = degrees["out"] - degrees["in"] 
            if diff == 1:
                start.add(node)
            elif diff == -1:
                end.add(node)
                
        if len(start) == 1 and len(end) == 1:
            return "".join(start)
        
        elif len(start) == 0 and len(end) == 0:
            random.choice(nodes.keys())
           
    
    def dfs(self):
        """
        Use a depth-first search algorithm to traverse each edge to a new node
        """
        
        graph = self.create_graph()
        node = self.get_nodes()
        edge_count = self.count_edges()
        path = list()

        def traversal(node):

            #choose new node at random
            while edge_count[node]["out"] != 0: 
                
               #Get the connected nodes which have an untraversed connected edge from the current node
               connected_nodes = list()
               for next_node, edges in graph[node].items():
                   if edges > 0:
                     connected_nodes.append(next_node)  
                     
               connected_node = random.choice(connected_nodes)
               
               graph[node][connected_node] -= 1
               
               edge_count[node]["out"]  -= 1
               
               traversal(connected_node)
               
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
                sequence += node[self.k-2:]  
                
        return(sequence)
            

    

if __name__ == "__main__":
    """
    If the script is run as main script
    """
     
    pass
        
