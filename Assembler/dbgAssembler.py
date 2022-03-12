#!/usr/bin/env python3

"""
Title: dbgAssembler.py
Created on Mon 07-03-2022
Author: Daniel Nilsson

Description:
    This program takes a DNA sequenc as input and an optional integer k.
    It then splits the string in overlapping kmers of size k, produces a De Bruijn graph
    of k-1 mers. Then finds an Eulerian path and outputs the reconstructed sequence.
    
    OBS, this is not a memory efficient script. It cannot handle to large sequences

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

from pathlib import Path #check existing path
import random #enables random choice of edge traversal
import argparse #enables user input
from datetime import date #to get the date

from Assembler.Utilities.check_valid_sequence import check_valid_sequence #to check for correct input
from Assembler.Utilities.readfastafiles import readFastQ_returnDict_Fasta #to read in fasta file
from Assembler.Utilities.readfastafiles import readFasta_returnDict #to read in fasta file
from Assembler.Utilities.output_manager import output_results #to write result to a file

class DbgSolver:
    """
    class to create kmers of length k from a DNA sequence, create De Bruijn graph using 
    k-1 mers, and reconstruction of the DNA sequence, by finding a eulerian path
        
    """
    def __init__(self, sequence, k=None):
        """
        intialize variables, k is optional
        """
        max_k = 256
        if not check_valid_sequence(sequence): #check for correct sequence input
            #raise error otherwise
            raise Exception("Error: DNA sequence contains ambiguous characters") 
        
        self.sequence = sequence.upper() #make the sequence to contain only uppercase letters
        
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
        
        self.sequence = "" #reset the sequence
        return edges

        
    
    def create_graph(self):
        """
        Create the node by splitting up the edges in pairs of k-1mers
        These will be used to create the De Bruijn graph
        """
                
        edges = self.__create_kmers() #get the edges (kmers)
        graph = dict() #create a empty dictionary to be  with nodes 
        
        for node, repeat in edges.items(): #get the edge and the repeat count

            left, right =  node[:-1], node[1:] #create the nodes by creating two sequences with k-1 in overlap 
            
            if left not in graph: #if the left node is not in the graph, assign it as a key
                graph[left] = dict()
            
            if right not in graph: #if the right node is not in the graph, assign it as a key
                graph[right] = dict()
                
            #As left and right are from the same kmer, this will always be true in terms of overlap
            if right not in graph[left].keys():              
                graph[left].update({right:repeat}) #assign the repeat count as this represent the number of edges to th node
        
        self.graph = graph #create a graph for the class
        return graph  #also return a graph if the case one would only want to use the graph
        
    
    def count_edges(self):
        """
        Represent the de Bruijn graph as an adjacency matrix 
        Here we merge repeated nodes to an unique node 
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
              
        self.degrees = nodes  

    
    
    def get_nodes(self):
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
                start.append((node, list(self.graph[node].keys()))) #add the node and the neighbors as a tuple
                
            elif diff == -1: #if there is one more in degree, append to the end list
                end.append((node, list(self.graph[node].keys()))) #add the node and the neighbors as a tuple
                
        if len(start) == 1 and len(end) == 1: #if there are two unblanced nodes, we have an eulerian path. return the start value if that is the case
            return start[0]
        
        else: #else, just return a random node
            node = random.choice(list(self.graph.keys())) #choose a random node
            neighbors = random.choice(list(self.graph[node].keys()))  #add the node and the neighbors as a tuple       
            return node, neighbors
    
            
    def __traversal(self):
         """
         Traverse the de Bruijn graph by chosing random edges from the current node.
         When an end is reached (a node with no outgoing edges) append the node to the and backtrack the path
         """

         self.create_graph() #create the graph
         
         unvisited = list() #create a list to track the nodes which have untraversed edges
         path = list() #create a list to store the path taken through the de bruijn graph
         
         #assign the start node
         self.count_edges()
         node, neighbors = self.get_nodes()        
         unvisited.append(node)
         #While there are unvisited edges in graph:
         while unvisited: #while there still is untraversed nodes
              node = unvisited[-1] #assign the last element in the unvisited list to be the next node
              neighbors = list(self.graph[node].keys())#get the neighbors
              if self.degrees[node]["out"] == 0: #if the out degree is 0, we have reached an end
                  path.append(node) #append the node to the path (starting from the last node)
                  unvisited.pop(-1) #remove the node from the unvisited list
                  
  
              #if there still is remaning edges
              else:
                  self.degrees[node]["out"] -= 1 #remove one from the total count of outgoing edge
                  next_node = random.choice(neighbors) #get one random node from the neighbors   
                  
                  self.graph[node][next_node] -= 1  #remove the chosen neighbors from the current node
                  if self.graph[node][next_node] == 0:
                      self.graph[node].pop(next_node)   
                      
                  unvisited.append(next_node) #add the next node to the  unvisited list as this will be the next starting point
                   
            
                    
         return path[::-1] #revrese the order of the path as the order in which the element have been added are in the reverse order (backtracking)
       
 
    def assemble(self):
        """
        Takes the eulerian path/cycle and assembles it
        It will continue to call the 
        """
          
        scaffold = dict() #a dictionary to hol the assembled sequences
        i = 1 #index for the header
        
        c = self.__traversal() #start by calling traversal function 
       
        sequence = "" # and add to the sequence and set the first to false
        first = True
        #Missing step:  convert the cycle c into the sequence:
        for node in c: #iterate through the path list
            if first: #check for the first element 
                sequence = "".join([sequence, node])  # and add to the sequence and set the first to false
                first = False
            else: #Add only the last character of all the other nodes that are not the first
                sequence = "".join([sequence, node[-1]]) 
                
                
        scaffold[f"Scaffold_{i}"] = sequence #add the sequence to the dictionary for output
        i += 1
        
        
        while any(self.graph.values()): #while there still are untraversed edges    
            c = self.__traversal() #call traversal function again
                
            sequence = "" # and add to the sequence and set the first to false
            first = True
            #Missing step:  convert the cycle c into the sequence:
            for node in c: #iterate through the path list
                   
                if first: #check for the first element 
                   sequence = "".join([sequence, node]) # and add to the sequence and set the first to false
                   first = False
                   
                else: #Add only the last character of all the other nodes that are not the first
                   sequence = "".join([sequence,node[-1]]) 
                    
            scaffold[f"Scaffold_{i}"] = sequence #add the sequence to the dictionary for output
            i += 1

        return scaffold
      
"".join(["ads", "dsad"[-1]]) 

def main():
    """
    When called as main, get the input from the user
    """
    today = date.today()
    current_date = today.strftime("%d_%m_%Y")
    output = f"dbgAssembler_run_{current_date}.fna"
    directory = f"dbgAssembler_{current_date}"
    #create a parser for the command line
    parser = argparse.ArgumentParser(usage="""%(prog)s python dbgAssembler.py -i <sequence> -f <input_file> -k <kmer_size> [optional] -o <output_file> [optional] \nType -h/--help for the help message""",
                                description="This program takes a dna string or a fasta file as input and breaks the sequence into kmers of size k. Then reassembles the string using a De Bruijn graph approach")
   
    
    #create arguments for inputs: sequence, file and k
    parser.add_argument("-i", metavar='<input_sequence>', type=str,  help="a DNA string", default=None)
        
    parser.add_argument("-f", metavar='<input_file>',  help="path to fasta file", default=None)
    
    parser.add_argument("-t", metavar='<fasta/fastq>',  help="Specify whatever the input file is in fasta or fastq format", default=None)
    
    parser.add_argument("-k",  metavar='<kmer_size>', type=int, help="kmer size (default 1/3 of the sequence length, max: 256)")
    
    parser.add_argument("-o", metavar='<output_file>', type=str,  help="path to outputfile, default: dbgAssembler_run{current_date}.fna", default=output)
    
    parser.add_argument("-d", metavar='<directory>', type=str,  help="name of output directory to create, default: dbgAssembler_{current_date}", default=directory)
    
    #assign the parsed input to a variable
    args = parser.parse_args()
    
    #assign each input to a new variable
    sequence = args.i
    input_file = args.f
    type_f = args.t  
    k = args.k 
    output = args.o 
    directory = args.d
    
    #if a sequence is specified, then the format cannot be specified
    if sequence and type_f:
        raise Exception("Error: format type cannot be specified with a string (sequence) as input")    
 
    #both sequence and a file cant be used together
    if sequence and input_file:
        raise Exception("Only one of -i and -f can be specified at a time")
    
    
    elif input_file:
        if not type_f: #check if the format is specified
            raise Exception("Error: The format of the input file must be specified") #give a message that the type format must be specified
        
        #check whatever the format is given correctly
        elif type_f.lower() != "fastq" and type_f.lower() != "fasta":
            raise Exception("Error: Invalid format") #give a message that the format is invalid
            
        if Path(input_file).is_file():  #check whatever the input file is an actual file..   
            if type_f.lower() == "fasta": #if it is in fasta format
              
              #get the sequences:
              sequences = "" 
              
              #create one large sequence
              for seq in readFasta_returnDict(input_file).values():
                 sequences = "".join([sequences, seq])
                  
              #call the DgbSolver with k or without
              if k: 
                  graph = DbgSolver(sequences, k)
              else:
                 graph = DbgSolver(sequences)  
                 
         
            elif type_f.lower() == "fastq": #if it is in fastq format
            
                #create an empty string to hold the one-line sequence
                sequences = ""

                #convert to fasta and create one large sequence
                for seq in readFastQ_returnDict_Fasta(input_file).values(): #
                    sequences = "".join([sequences, seq])

                #call the DgbSolver with k or without
                if k:
                    graph = DbgSolver(sequences, k)
                else:
                   graph = DbgSolver(sequences)     
                   
        else: #raise an error in the 
            raise Exception(f"Error: {input_file} is not a file or does not exist")
               
        
    elif not sequence: #if neither a sequence or a file is present
        raise Exception("Error: No input specified")
        
    #check if k is specified and create a DbgGraph object with the k or without if not specified 
    elif k: 
        graph = DbgSolver(sequence, k)
    else:
       graph = DbgSolver(sequence) 
    
    
    #get the reconstructed sequence and print the sequence to standard output
    new_seq = graph.assemble()
    
    
    
    #Tell the user if the reassembly was succesful or not by comparing 
    #the input sequence with the reconstructed
    #write to output file
    if not sequence and input_file:    
        output_results(input_file, graph.k, new_seq, output, directory) #call the function to write to output files
        
    elif sequence and not input_file: #specify if the input was a file or a single string     
            output_results(sequence, graph.k, new_seq, output, directory, onlySeq = True)  #call the function to write to output files    

    



if __name__ == "__main__":
    """
    If the script is run as main
    """
    main()
    
 
    
