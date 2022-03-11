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

The two first alternatives were implemented with a recursive function. This did not work for longer sequences as maximum recursion depth was exceeded. Possible solution: increase recursion limit.

Third alternative implemented without a recursive function and had a eaiser implementation
to follow. Thus, the third alternative was chosen.

Last change, changed path walking to be completly random. Compressed the code and
no recursive func is used.
 * Test on genome size 5kb gives a 100% correct assembly
 * Test on 170 kb gives a ~96% correct assembly

#### Eulerian path/cycle
To do:
  * Create a script to compare larger output
  * Create a Utility file for writing output to a separate folder
  * Create a test script for genome files <= 20 kb
