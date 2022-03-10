#!/usr/bin/env python3

"""
Title: dbgAssembler.py
Created on Mon 07-03-2022
Author: Daniel Nilsson

Description:
    This program takes a DNA sequenc as input and an optional integer k.
    It then splits the string in overlapping kmers of size k, produces a De Bruijn graph
    of k-1 mers. Then finds an Eulerian path and outputs the reconstructed sequence

List of functions:
    main()

List of classes: 
    DbgGraph - functions:
        __create_kmers()
        create_graph()
        count_edges()
        get_nodes()
        dfs(): 
            nested func - traversal()
        get_sequence()
             

List of "non standard" modules:
    check_valid_sequence - to handle incorrect input format of the sequence (should be DNA)
   
Procedure:
    1. Create kmers of specified size
    2. Create de bruijn graph
        * edges being the kmer
        * nodes k-1mers of the kmers
        * only uniq kmers will be used
        * repeated kmers are counted
    4. Check Eulerian path/cycle
    5. Use Depth-first search algorithm to traverse the edges and nodes
        * traverse the same edge only once
    6. When a dead end is found, put the node at the last position for the construction of
       the new sequence, then backtrack until a node with an unused edge is found and start over the path search

Usage:
    python dbgAssembler.py <sequence> -k <kmer_size> [optional]

"""
  
import random #enables random choice of edge traversal
import argparse #enables user input
from check_valid_sequence import check_valid_sequence #to check for correct input

class DbgGraph:
    """
    class to create kmers of length k from a DNA sequence, create De Bruijn graph using 
    k-1 mers, and reconstruction of the DNA sequence, by finding a eulerian path
        
    """
    def __init__(self, sequence, k=None):
        """
        intialize variables, k is optional
        """
        
        if not check_valid_sequence(sequence): #check for correct sequence input
            #raise error otherwise
            raise Exception("Error: DNA sequence contains ambiguous characters") 
        
        self.sequence = sequence.upper() #make the sequence to contain only uppercase letters
        
        #handle value of k
        if not k: #if k is not specified, give k value 1/3 of the sequence length
            k = int(len(self.sequence)/3)
            
        if k > len(self.sequence): #if the specified k is longer than the length of the sequence
            k = int(len(self.sequence)/3) #then set it to be 1/3 of the sequence length
            
            if k > 256: #set an upper limit of k to 256 
                print("Max kmer size set reached. Using size 256")
                k=256
                
            else: #if k is not bigger than 256
                print(f"kmer size cannot be longer than the sequence\n1/3 of the sequence length: {k} will instead be used")
            
        if k > 256: #if user specified k > 256 or 1/3 of the length is > 256, set k to 256
            print("Max kmer size set reached. Using size 256")
            k=256
            
        self.k = k #intialize k
        
    def __create_kmers(self):
        """
        Split the sequence in kmers of size k and add to the dictionary edges. 
        Associate each kmer with its repeat count
        
        Nr of kmers = L-k+1, L is the length of the input sequence, k is the kmer size
        """
        edges = dict() #represent all edges by a dictionary
                       #This is to keep track of the repeat of the kmer
        
        #get length of the sequence
        sequence_length = len(self.sequence)
        
        #iterate over the length to get the start index and the length of the sequence 
        #starting at position k to get the end index for slicing
        for start, end in zip(range(sequence_length), range(self.k, sequence_length+1)):
            #slice the sequence in kmers using the indices:
            kmer = self.sequence[start:end]
            
            if kmer not in edges.keys(): #if the kmer is not present, add to dict and give the repeat count 1
                edges[kmer] = 1
            else: #if already present, increase the repeat count by 1
                edges[kmer] += 1
        
        return edges #return the edges
        
    
    def create_graph(self):
        """
        Create the node by splitting up the edges in pairs of k-1mers
        These will be used to create the De Bruijn graph
        """
        
        edges = self.__create_kmers() #firstly, get the edges (kmers)
        graph = dict() #store the nodes in a dictionary
                       #the De Bruijn graph will be represented as an adjacency list with slight modification taking repeats into account
        
        for node, repeat in edges.items(): #get the edge and the repeat count

            left, right =  node[:-1], node[1:] #create the nodes which will create two nodes per kmer, 
                                               #the nodes will have a k-1 overlap to each other
            
            #populate the adjacency list, all nodes needs to be present as a key
            #in the dictionary
            if left not in graph: 
                graph[left] = dict()
            
            if right not in graph:
                graph[right] = dict()
                
            #As left and right are from the same kmer, this will always be true in terms of 
            #overlap, also set the next node (value) to have the value of the repeat.
            #This represents the number of edges pointing towards the connected node
            if right not in graph[left].keys():              
                graph[left].update({right:repeat})
                
        self.graph = graph
        return graph       
    
    
    def count_edges(self, graph):
        """
        Count the number of outgoing and ingoing edges for each node
        """
        edges = self.__create_kmers()
        nodes = dict() #create a dictionary to store each node with its in and out degrees
        
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
    
    
    def get_nodes(self, edge_count):
        """
            If an eulerian trail is present, then we have one node with:
                
            * outdegree - indegree = 1 (start) and one with
            * indegree - outdegree = 1 (end)
        
            If instead we have an eulerian cycle, all nodes are balanced:  outdegree = indegree
            Then start at a random node
        """
        #create two sets to store the potential start and end node
        start = set() 
        end = set()
        
        
        for node, degrees in edge_count.items():  #start by looping through the in and out degree
            diff = degrees["out"] - degrees["in"] #check the difference
            if diff == 1: #if 1, we have a start node as the out degree is larger
                start.add(node)
            elif diff == -1: #if -1, we have a higher in degree than out, meaning it is a end node
                end.add(node)
        
        #if the sets start and end have length 1, we have a eulerian trial, then get the start node
        if len(start) == 1 and len(end) == 1:
            return "".join(start) #get the string value
        
        #if all nodes are balanced, then choose one at random to return as we do not know the start
        #for an eulerian cycle
        elif len(start) == 0 and len(end) == 0:
            return random.choice(edge_count.keys())
           
    
    def dfs(self):
        """
        Use a depth-first search algorithm with backtracking to traverse each edge to a new node.
        When reaching a node with no edges left, add it to a list representing the path. Then 
        go back to find a node with unused edges and start over until next stop and so on.
        The order of the list needs to be reversed as the first element to be put in is the end node
        """
        
        graph = self.create_graph() #create the Bruijn graph
        edge_count = self.count_edges(graph) #fecth the degree for each node
        node = self.get_nodes(edge_count) #get the start node
        
        path = list() #create a list to hold the path

        #A recursive func for the traversal
        def traversal(node): #input is the node to traverse from

            #while the node has unused edges, jump to one of the connected nodes
            while edge_count[node]["out"] != 0: 
                
               #Get the connected nodes which have an untraversed connected edge from the current node
               connected_nodes = list() #create a list to store possible nodes to jump to
               for next_node, edges in graph[node].items(): #get the connected nodes and their edges
                   if edges > 0: #if the edge count is not 0, we still have unused edges
                     connected_nodes.append(next_node)  #add these nodes to the connected_nodes list
                     
               connected_node = random.choice(connected_nodes) #choose a connected node at random
               
               #now remove the edge that connects the node and its connected_node
               graph[node][connected_node] -= 1
               
               #decrease the nodes out degree by one as we have removed one edge
               edge_count[node]["out"]  -= 1
               
               #call the function on the next node to repeat the process
               traversal(connected_node)
               
            #if the out degree of the node is 0, append the node to the path
            path.append(node)       
            
        #call the traversal func
        traversal(node)
        
               
        return path[::-1] #return the reversed order of the path
            

    def get_sequence(self):
        """
        This function reconstructs the sequence from the eulerian path by
        concatenating each elements last character in the list to a string. All 
        the characters of the first element int the list will be used as this is the start        
        """
        
        sequence = "" #create an empty string to use for the concatenation
        first = True #create a boolean to track the first element in the list
        
        path = self.dfs() #get the eulerian path 
        
        for node in path: #iterate through the path list
            if first: #check for the first element 
                sequence += node # and add to the sequence and set the first to false
                first = False
            else: #Add only the last character of all the other nodes that are not the first
                sequence += node[-1]  
                
        return(sequence) #return the reconstructed sequence
            


def main():
    """
    When called as main, get the input of 
    """
    
    #create a parser for the command line
    parser = argparse.ArgumentParser(usage="""%(prog)s python dbgAssembler.py <sequence> -k <kmer_size> [optional] \nType -h/--help for the help message""",
                                 description="This program takes a dna string as input and breaks it to kmers of size k. Then reassembles the string using a De Bruijn graph styled manner")
    
    #create arguments for inputs: sequence and k
    parser.add_argument("i", metavar='<input_sequence>', type=str,  help="a DNA string")
    
    parser.add_argument("-k",  metavar='<kmer_size>', type=int, help="kmer size (default 1/3 of the sequence length, max: 256)")
    
    #assign the parsed input to a variable
    args = parser.parse_args()
    
    #assign each input to a new variable
    sequence = args.i
    k = args.k 
    
    #check if k is specified and create a DbgGraph object with the k or without if not specified 
    if k: 
        graph = DbgGraph(sequence, k)
    else:
       graph = DbgGraph(sequence) 
    
    
    #get the reconstructed sequence and print the sequence to standard output
    new_seq = graph.get_sequence()
    print(f"\nReassembly: {new_seq}")
    
    #Tell the user if the reassembly was succesful or not by comparing 
    #the input sequence with the reconstructed
    print(f"Successful assembly: {sequence == new_seq}")
    

if __name__ == "__main__":
    """
    If the script is run as main script
    """    
    main() #call main
        
