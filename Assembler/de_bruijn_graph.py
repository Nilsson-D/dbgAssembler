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
    The assembler cannot handle circular genomes well (though, it is possible to input)
    The implementations assumes a linear genome

List of functions:
    main()

List of classes: 
    DbgGraph - functions:
        __create_kmers()
        create_graph()
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
    5. Use Depth-first search algorithm to traverse the edges and nodes
        * traverse the same edge only once
    6. When a dead end is found, put the node at the last position for the construction of
       the new sequence, then backtrack until a node with an unused edge is found and start over the path search

"""


import random #enables random choice of edge traversal


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
               
        self.sequence = "".join(InputManager(input_file).read_fasta(allowN)) #make the sequence to contain only uppercase letters
        
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
            
        self.k = k #intialize k
        
        self.graph = dict() 
        
    def __create_kmers(self):
        """
        Split the sequence in kmers of size k and add to the dictionary edges. 
        Associate each kmer with its repeat count
     
        Nr of kmers = L-k+1, L is the length of the input sequence, k is the kmer size
        """
        print("Creating kmers")
        
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
        
        print("Done creating kmers")
        self.sequence = "" #reset the sequence
        
        self.edges = edges
        return edges

        
    
    def create_graph(self):
        """
        Create the node by splitting up the edges in pairs of k-1mers
        These will be used to create the De Bruijn graph
        """               
        self.__create_kmers()
        print("Constructing graph")
        
        graph = dict() #create a empty dictionary to be  with nodes 
        
        for node, repeat in self.edges.items(): #get the edge and the repeat count

            left, right =  node[:-1], node[1:] #create the nodes by creating two sequences with k-1 in overlap 
            
            if left not in graph: #if the left node is not in the graph, assign it as a key
                graph[left] = dict()
            
            if right not in graph: #if the right node is not in the graph, assign it as a key
                graph[right] = dict()
                
            #As left and right are from the same kmer, this will always be true in terms of overlap
            if right not in graph[left].keys():              
                graph[left].update({right:repeat}) #assign the repeat count as this represent the number of edges to th node
        
        
        print("Done with constructing graph")
        self.graph = graph #create a graph for the class
        return graph  #also return a graph if the case one would only want to use the graph
        
    
    def __count_edges(self):
        """
        Represent the de Bruijn graph as an adjacency matrix 
        Here we merge repeated nodes to an unique node 
        """

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
              
        self.degrees = nodes  
  
    
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
        
        for node, degrees in self.degrees.items(): #start by looping over each node and its degree
            diff = degrees["out"] - degrees["in"] 
            if diff == 1: #if there is one more out degree, append to the start list
                start.append(node) #add the node to keep track of the unbalanced node
                
            elif diff == -1: #if there is one more in degree, append to the end list
                end.append(node) #add the node to keep track of the unbalanced node
                
        if len(start) == 1 and len(end) == 1: #if there are two unblanced nodes, we have an eulerian path. return the start value if that is the case
            return start[0]
        
        else: #else, just return a random node
            node = random.choice(list(self.graph.keys())) #choose a random node      
            return node
    
            
    def __traversal(self):
         """
         Traverse the de Bruijn graph by chosing random edges from the current node.
         When an end is reached (a node with no outgoing edges) append the node to the and backtrack the path
         """

         self.create_graph() #create the graph
         unvisited = list() #create a list to track the nodes which have untraversed edges
         path = list() #create a list to store the path taken through the de bruijn graph
         
         #assign the start node
         self.__count_edges()
         node = self.__get_nodes()    
         
         print("Finding path...")
         
         unvisited.append(node)
         #While there are unvisited edges in graph:
         while unvisited: #while there still is untraversed edges
              node = unvisited[-1] #assign the last element in the unvisited list to be the next node
              neighbors = list(self.graph[node].keys())#get the neighbors
              if self.degrees[node]["out"] == 0: #if the out degree is 0, we have reached an end
                  path.append(node) #append the node to the path (starting from the last node)
                  unvisited.pop(-1) #remove the node from the unvisited list as there is no edges left
                  
  
              #if there still is remaining edges
              else:
                  self.degrees[node]["out"] -= 1 #remove one from the total count of the outgoing edge for the node
                  next_node = random.choice(neighbors) #pick one random node from its neighbors                                          
                  unvisited.append(next_node) #add the next node to the  unvisited list as this will be the next starting point
                   
            
                    
         return path[::-1] #revrese the order of the path as the order in which the element have been added are in the reverse order (backtracking)
       
 
    def assemble(self):
        """
        Takes the eulerian path/cycle and assembles it
        It will continue to call the 
        """
        
        c = self.__traversal() #start by calling traversal function 
        print("Creating assembly") 
       
        sequence = "" # and add to the sequence and set the first to false
        first = True
        #Missing step:  convert the cycle c into the sequence:
        for node in c: #iterate through the path list
            if first: #check for the first element 
                sequence = "".join([sequence, node])  # and add to the sequence and set the "first" to false
                first = False
            else: #Add only the last character of all the other nodes that are not the start
                sequence = "".join([sequence, node[-1]]) #as the path represents the overlap by one, we cannot add the enitre node
                               
        
        print("Done :)")
        return sequence
      



if __name__ == "__main__":
    """
    If the script is run as main
    """
    pass
    
