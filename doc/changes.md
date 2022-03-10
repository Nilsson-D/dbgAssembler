## Documentation for dbgAssembler.py


#### Representation of de Bruijn Graph:
Three different structures were tried:
  * adjacency matrix
  * adjacency list, using a dictionary
    - Each repeat was stored in an additional dictionary as an integer with the cooresponding node
      - Ex {ATG:{TGC:2}}, The node ATG has two out going edges connected to the node TGC
  * adjacency list, using a dictionary
    - Each repeat was displayed as multiple elements in a list for outgoing edges of a node
      - Ex {ATG:[TGC, TCG]} The node ATG has two out going edges connected to the node TGC

The two first alternatives were implemented with a recursive function. This did not work for longer sequences. Third alternative implemented without a recursive function and had a eaiser implementation
to follow. Thus, the third alternative was chosen.

#### Eulerian path/cycle
To do:
  * Change choice of node to be random in the loop
    - Do this by sample a random node
    - get the index of the node and slice the list to exclude the chosen node.
      - obs! Do not remove from graph, this changes the graph for later iterations

  * Create a script to compare larger output
    - map func?
    - regular expression iter with set to check index?

  * Create threading if possible

  * Fix test script import
  * Create advanced test script for larger sequences
  * Create test script for visualization
  * Correct the format checker
