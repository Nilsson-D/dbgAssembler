# Simple Genome Assembler

### Description
dbgAssembler is a simple genome/sequence assembler for showing the utility of
a de bruijn graph based approach. dbgAssembler uses a string as an input or a fasta file to
generate kmers of the sequence and reassembles it by finding the eulierian path for the de bruijn graph.

### Dependencies
* pygraphviz 1.9
* pandas     1.4.1
* numpy      1.22.2
* networkx   2.5
* matplotlib 3.5.1

### Test data
The following test data was used for evaluating the assembly
 * human corona virus - RefSeq Assembly Accession: GCF_009858895.2
 * human c.fna - RefSeq Assembly Accession: GCF_001274345.1
 * Human adenovirus - NCBI Reference Sequence: NC_012959.1
 * Human gammaherpesvirus 4 (Epstein-Barr virus) - NCBI Reference Sequence: NC_007605.1

 fastq reads
 * novel aviadenovirus Oriolus adenovirus - SRR17301511
