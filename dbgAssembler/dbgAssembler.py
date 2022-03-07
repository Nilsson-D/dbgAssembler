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
        
        return edges
        
    
    def create_graph(self):
        """
        Create the node by splitting up the edges in pairs of k-1mers
        """
        edges = self.create_uniq_kmers()
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
            
        return graph       
    
    
    def count_edges(self, graph):
        """
        Count the number of outgoing and ingoing edges for each node
        """
        edges = self.create_uniq_kmers()
        nodes = dict() #create a dictionary to store each node with its in and out degreees
        
        for node in graph.keys(): #loop over each node
            
            nodes[node] = {"in": 0, "out": 0} #create nested dictionaries to store in and out degrees
            
            #count the out degree
            for connected_node in graph[node]:
                edges = graph[node][connected_node] #The number of repeats will be the number of 
                                                         #edges for one connected node
                nodes[node]["out"] += edges #add the edges to one connected node to get the out degree of the current node
            
            #get the in edges by checking the keys in which 
            #the node is the value for
            for connected_nodes in graph.values():
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
           
    def simplification(self, g):
        """
        Merge nodes with where one node only has one out degree
        """
        graph = self.create_graph()
        in_out_degree = self.count_edges(g)
        
        to_remove = set()
        to_update = dict()
        for node, degrees in in_out_degree.items():
            if degrees["in"] == degrees["out"] == 1:
                
                next_node = "".join(graph[node].keys())
                merged_node = node + next_node[-1]
                             
                to_remove.add(node)
                to_update[next_node] = merged_node
                for prev_node, connected_nodes in graph.items():
                    if node in connected_nodes:                      
                        graph[prev_node][merged_node] = graph[prev_node].pop(node)
        

        #To do: merge nodes with one out degree and update graph                              
        print(graph)
        for next_node, merged_node in to_update.items():
            
            if next_node in graph.keys():
                graph[merged_node] = graph.pop(next_node)  
            
        return graph
    
    def dfs(self):
        """
        Use a depth-first search algorithm to traverse each edge to a new node
        """

        graph = self.simplification()
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
                sequence += node[-1]
                
        return(sequence)
            

    

if __name__ == "__main__":
    """
    If the script is run as main
    """
    test_seq = "ACTGACGTACGTACGTGTGAAGCTAGTCGCGCATTACGGGTTGAAACGTTGGTT"
    my_graph = DbgGraph(test_seq, 3) 
    
    edges = my_graph.create_uniq_kmers()
    
    
    g = my_graph.create_graph()  
    g
    
    g = my_graph.simplification(g)
    g
    my_graph.count_edges(g)
    my_graph.get_nodes()   
    
    
    new_seq = my_graph.get_sequence()
    new_seq
    
    count = 0
    while new_seq != test_seq:
       new_seq = my_graph.get_sequence()
       count += 1
        
    print(count)