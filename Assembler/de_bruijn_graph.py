#!/usr/bin/env python3

"""
Title: de_bruijn_graph.py
Created on Mon 07-03-2022
Author: Daniel Nilsson

Description:
    This program takes a DNA sequenc as input and an optional integer k.
    It then splits the string in overlapping kmers of size k, produces a De Bruijn graph
    of k-1 mers. Then finds an Eulerian path and outputs the reconstructed sequence.
    
    OBS, this is not a memory efficient script. It cannot handle too large sequences

List of functions:
    main()

List of classes: 
    DbgGraph - functions:
        __create_kmers()
        __create_graph()
        __count_edges()
        __get_nodes()
        traversal()
        assemble()
             

List of "non standard" modules:
    -
   
Procedure:
    1. Create kmers of specified size from the input sequence
        * repeated kmers are counted
    2. Create the de bruijn graph
        * edges will represent the kmer
        * nodes are the prefix and suffix of the kmer
    4. Check Eulerian path/cycle
    5. Use Hierholzer's algorithm to traverse all the edges exactly once
        * traverse the same edge only once
    6. When a dead end is found, put the node at the last position for the construction of
       the new sequence, then backtrack until a node with an unused edge is found and start over the path search

"""


import random #enables random choice of edge traversal
import sys

from Assembler.Utilities.input_manager import InputManager #enables reading the fasta file

class DbgSolver:
    """
    class to create kmers of length k from a DNA sequence, create De Bruijn graph using 
    k-1 mers, and reconstruction of the DNA sequence, by finding a eulerian path
        
    """
    def __init__(self, input_file, k=31, allowN=False):
        """
        intialize variables, k is optional
        """
        max_k = 251 #set a max kmer size
        min_k = 3
               
        self.sequence = "".join(InputManager(input_file).read_fasta(allowN)) #make the sequence to contain only uppercase letters
        
        if len(self.sequence) < 10:
            print("Error: Sequence length to short. Minimum length: 10 bp")
            sys.exit()
                    
        #handle value of k
        if not k: #if k is not specified, give k value 1/3 of the sequence length
            k = int(len(self.sequence)/3)
            
        if k > len(self.sequence): #if the specified k is longer than the length of the sequence
            k = int(len(self.sequence)/3) #then set it to be 1/3 of the sequence length
            
            if k > max_k: #set an upper limit of k to max_k 
                print("Max kmer size reached. Using size {max_k}")
                k=max_k
                
            else: #if k is not bigger than max_k
                print(f"kmer size cannot be longer than the sequence\n1/3 of the sequence length: {k} will instead be used")
            
        if k > max_k: #if user specified k > max_k or 1/3 of the length is > 256, set k to max_k
            print(f"Max kmer size reached. Using size {max_k}")
            k=max_k
        
        if k < min_k:
            print(f"k-mer size cannot be smaller than 3. Using size: {min_k}")
            k=min_k
            
        self.k = k #intialize k   
        self.min_contig_size = 150 #min contig length
        
    def __create_kmers(self):
        """
        Split the sequence in kmers of size k and add to the dictionary edges. 
        Associate each kmer with its repeat count
     
        Nr of kmers = L-k+1, L is the length of the input sequence, k is the kmer size
        """
        print("Creating kmers")
        
        edges = list()
        
        #get length of the sequence
        sequence_length = len(self.sequence)
        
        #iterate over the length to get the start index and the length of the sequence 
        #starting at position kto get the end index for slicing
        for start, end in zip(range(sequence_length), range(self.k, sequence_length+1)):
            #slice the sequence in kmers using the indices:
            kmer = self.sequence[start:end]
            
            edges.append(kmer) 
        
        print("Done creating kmers")
        
        return edges

        
    
    def create_graph(self):
        """
        Create the node by splitting up the edges in pairs of k-1mers
        These will be used to create the De Bruijn graph
        """               
        edges = self.__create_kmers()
        print("Constructing graph")
        graph = dict() #create a empty dictionary to be  with nodes 
        
        for node in edges: #get the edge and the repeat count
            left, right =  node[:-1], node[1:] #create the nodes by creating two sequences with k-1 in overlap  

            if left not in graph: #if the left node is not in the graph, assign it as a key
                graph[left] = list()
            
            if right not in graph: #if the right node is not in the graph, assign it as a key
                graph[right] = list()
                
            #As left and right are from the same kmer, this will always be true in terms of overlap
            
            graph[left].append(right) #assign the suffix to the prefix as this will always be true. In this graph, a key can hold multiple copies of one node corresponding to the number of edges between the two nodes
            

                     
        
        print("Done with constructing graph")
        self.graph = graph #create a graph for the class
        return graph  #also return a graph if the case one would only want to use the graph
        
    
    def __count_edges(self):
        """
        Represent the de Bruijn graph as an adjacency matrix 
        Here we merge repeated nodes to an unique node 
        """

        nodes = dict() #create a dictionary to store each node with its in and out degreees
        for node in self.graph: #loop over each node
            nodes[node] = {"in": 0, "out": 0} #create nested dictionaries to store in and out degrees
            
            #count the out degree
            nodes[node]["out"] = len(self.graph[node]) #add the edges to one connected node to get the out degree of the current node
            
            #get the in edges by checking the keys in which 
            #the node is the value for
            for connected_nodes in self.graph.values():
                if node in connected_nodes: #if the node is present as a value for another node, then
                                            #the number of out edges for the other node will 
                                            #be the number of in edges for the current node
                    edges = connected_nodes.count(node)
                    nodes[node]["in"] += edges #Sum all in edges to get the in degree of the current node
        
        return nodes 
  
    
    def __get_nodes(self):
        """
            If an eulerian trail is present, then we have one node with:
                
            * outdegree - indegree = 1 (start) and one with
            * indegree - outdegree = 1 (end)
        
            If instead we have an eulerian cycle, all nodes are balanced:  outdegree = indegree
            Then start at a random node
        """
        #create two sets to store the potential start and end node
        start = list()
        end = list()
        
        nodes = self.__count_edges()
        
        for node, degrees in nodes.items(): #start by looping over each node and its degree
            diff = degrees["out"] - degrees["in"] 
            if diff == 1: #if there is one more out degree, append to the start list
                start.append(node) #add the node to keep track of the unbalanced node
                
            elif diff == -1: #if there is one more in degree, append to the end list
                end.append(node) #add the node to keep track of the unbalanced node
                
        if len(start) == 1 and len(end) == 1: #if there are two unblanced nodes, we have an eulerian path. return the start value if that is the case
            return start[0]
        
        
        elif len(start) == 0 and len(end) == 0: #Check if the sequence is circular.
            node = random.choice(list(self.graph.keys())) #choose a random node    
            return node
        
        else: #This is for when there is no eulerian path or cicle. Make an attempt to solve it anyways
            node = random.choice(list(self.graph.keys())) #choose a random node      
            return node    
            
    def __traversal(self):
         """
         Traverse the de Bruijn graph by chosing random edges from the current node.
         When an end is reached (a node with no outgoing edges) append the node to the and backtrack the path
         """

         unvisited = list() #create a list to track the nodes which have untraversed edges
         path = list() #create a list to store the path taken through the de bruijn graph
         #assign the start node
         node = self.__get_nodes()            
         
         unvisited.append(node)        
         
         #While there are unvisited edges in graph:
         while unvisited: #while there still is untraversed edges
              node = unvisited[-1] #assign the last element in the unvisited list to be the next node
              neighbors = self.graph[node] #get the neighbors
              if len(neighbors) == 0: #if the out degree is 0, we have reached an end

                  path.append(node) #append the node to the path (starting from the last node)
                  unvisited.pop(-1) #remove the node from the unvisited list as there is no edges left

                
              #if there still is remaining edges
              else:                
                  next_node = random.choice(neighbors) #pick one random node from its neighbors, 
                                                       #the point of picking a random node instead of the always the first one in the list
                                                       #is to simulate an unkown path taking. When using real sequenced data, the order wont be
                                                       #as straightforward as using the first added edge (neighbors[0]) to the path...  
                  self.graph[node].remove(next_node) #remove the edge coming that we have taken
        
                  unvisited.append(next_node) #add the next node to the  unvisited list as this will be the next starting point
                   
            
                    
         return path[::-1] #reverse the order of the path as the order in which the element have been added are in the reverse order (backtracking)
       
 
    def assemble(self):
        """
        Takes the eulerian path/cycle and assembles it
        It will continue to call the 
        """
        
        self.create_graph() #create the graph
 
        print("Creating assembly") 
        
        contigs = list()
        
        
        while any(self.graph.values()): #create all paths to get as many contigs as possible
            
            path = self.__traversal()

            sequence = ""
            if len(path) > self.min_contig_size: #define a minimum contig length 
                sequence = "".join([sequence, path[0]])

                sequence += "".join([node[-1] for node in path[1:]])
            
                                                                                     
                contigs.append(sequence)                                             
                                                                                                                                                                             
        print("Done :)")
        return contigs
    


if __name__ == "__main__":
    """
    If the script is run as main
    """
    pass