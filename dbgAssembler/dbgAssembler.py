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
        __create_graph()
        __count_edges()
        __get_nodes()
        traversal()
        assemble()
             

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

from pathlib import Path
import random #enables random choice of edge traversal
import argparse #enables user input
from check_valid_sequence import check_valid_sequence #to check for correct input
from readfastafiles import readFasta_returnDict



class DbgSolver:
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
                print("Max kmer size reached. Using size 256")
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
        edges = list() #represent all edges by a dictionary
                       #This is to keep track of the repeat of the kmer

        #get length of the sequence
        sequence_length = len(self.sequence)

        #iterate over the length to get the start index and the length of the sequence 
        #starting at position k to get the end index for slicing
        for start, end in zip(range(sequence_length), range(self.k, sequence_length+1)):
            #slice the sequence in kmers using the indices:
            kmer = self.sequence[start:end]

            edges.append(kmer)

        return edges

        
    
    def __create_graph(self):
        """
        Create the node by splitting up the edges in pairs of k-1mers
        These will be used to create the De Bruijn graph
        """
                
        edges = self.__create_kmers() #firstly, get the edges (kmers)
        graph = dict() #store the nodes in a dictionary
                       #the De Bruijn graph will be represented as an adjacency list with slight modification taking repeats into account
        
        for node in edges: #get the edge and the repeat count

            left, right =  node[:-1], node[1:] #create the nodes which will create two nodes per kmer, 
                                               #the nodes will have a k-1 overlap to each other
            
            #populate the adjacency list, all nodes needs to be present as a key
            #in the dictionary
            if left not in graph: 

                graph[left] = list()
            
            if right not in graph:
                graph[right] = list()
                
            #As left and right are from the same kmer, this will always be true in terms of 
            #overlap
            graph[left].append(right)
            
        self.graph = graph
    
    
    def __count_edges(self):
        """
        Represent the de Bruijn graph as an adjacency matrix 
        Here we merge repeated nodes to an unique node 
        """

        nodes = dict() #create a dictionary to store each node with its in and out degrees
        
        for node in self.graph.keys(): #loop over each node
            
            nodes[node] = {"in": 0, "out": 0} #create nested dictionaries to store in and out degrees
            
            edges = len(self.graph[node]) #The number of repeats will be the number of 
                                                         #edges for one connected node
            nodes[node]["out"] = edges #add the edges to one connected node to get the out degree of the current node
            
            #get the in edges by checking the keys in which 
            #the node is the value for
            for connected_nodes in self.graph.values():
                if node in connected_nodes: #if the node is present as a value for another node, then
                                            #the number of out edges for the other node will 
                                            #be the number of in edges for the current node
                    nodes[node]["in"] += connected_nodes.count(node) #Sum all in edges to get the in degree of the current node
                    
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
        start = dict() 
        end = dict()
        
        self.__count_edges()
        
        for node, degrees in self.degrees.items():  #start by looping through the in and out degree
            diff = degrees["out"] - degrees["in"] #check the difference
            if diff == 1: #if 1, we have a start node as the out degree is larger
                return node, self.graph[node]
            elif diff == -1: #if -1, we have a higher in degree than out, meaning it is a end node
                end[node] = self.graph[node]
        
        #if the sets start and end have length 1, we have a eulerian trial, then get the start node
        if len(start) == 1 and len(end) == 1:
            return start #get the string value
        
        #if all nodes are balanced, then choose one at random to return as we do not know the start
        #for an eulerian cycle
        elif len(start) == 0 and len(end) == 0:
            node = random.choice(list(self.degrees.edge_count.keys()) )
            return node, self.graph[node]         
    
            
    def __traversal(self):
         """
         Form a cycle Cycle by randomly walking in Graph
         """
         #Put all edges into the unvisited edges:
         self.__create_graph()
         unvisited = self.graph.copy()
              
         #assign the start node
         start_node, value = self.__get_nodes()
         
         
         #The get a random neighbor to the start node
         node = random.choice(value)
         
         #assign the first node to the path/cycle depending on the of all nodes are balanced
         cycle = [start_node, value[unvisited[start_node].index(node)]]
         
         
         unvisited[start_node].pop(unvisited[start_node].index(node)) #pop chosen connected node

         if len(value) == 1: #if the node only has one edge, we will set it as explored. Thus, we remove it from our unvisited list
             unvisited.pop(start_node)
             
         #While there are unvisited edges in graph:
         while unvisited:
              #Check if the last node in the cycle still is present in the unvisited graph:
              if cycle[-1] in unvisited:
                   neighbors = unvisited.pop(cycle[-1]) #if so remove the last node
                   if len(neighbors) > 0: #check if the last node has any connected nodes
                        if len(neighbors) > 1:  #if more than one edge
                             #Put back the remaining unvisited edges:
                             unvisited[cycle[-1]] = neighbors[1:]       
                             
                        #Add the connected node whatever there only is one node or multiple connected nodes
                        cycle.append(neighbors[0])
    
              #If the last node in the cycle does not have any remaining edges.
              #This is only applicable if we have an euileran cycle
              else:                   
                   for i in range(len(cycle)): #iterate over the nodes in the path
                        if cycle[i] in unvisited: #when we find one node with unvisited nodes, get the connected nodes
                             neighbors = unvisited.pop(cycle[i]) #get the connected nodes from the graph
                             if len(neighbors) > 0: #check how many edges there are
                             
                                  cycle = cycle[:i] + cycle[i:] #this is needed if we have a cycle instead of a path
                                  if len(neighbors) > 1: #if more than one edge
                                       #Put back the remaining unvisited edges:
                                       unvisited[cycle[-1]] = neighbors[1:]
                                  #Add the connected node whatever there only is one node or multiple connected nodes
                                  cycle.append(neighbors[0])
                             break #break the for loop and start over at the last node
    
         return cycle

        
    def assemble(self):
        """
        Takes kthe eulerian path/cycle and assembles it
        """
           
        c = self.__traversal()
        
        sequence = "" # and add to the sequence and set the first to false
        first = True
        #Missing step:  convert the cycle c into the sequence:
        for node in c: #iterate through the path list
               
            if first: #check for the first element 
                sequence += node # and add to the sequence and set the first to false
                first = False
            else: #Add only the last character of all the other nodes that are not the first
                sequence += node[-1]  
                

        return(sequence)
    


def main():
    """
    When called as main, get the input of 
    """
   
    #create a parser for the command line
    parser = argparse.ArgumentParser(usage="""%(prog)s python dbgAssembler.py <sequence> -k <kmer_size> [optional] \nType -h/--help for the help message""",
                                description="This program takes a dna string as input and breaks it to kmers of size k. Then reassembles the string using a De Bruijn graph styled manner")
   
    
    #create arguments for inputs: sequence, file and k
    parser.add_argument("-i", metavar='<input_sequence>', type=str,  help="a DNA string", default=None)
        
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
    

