# dbgAssembler

### Description
dbgAssembler is a simple genome/sequence assembler using a de Bruijn graph based approach.
The dbgAssembler uses a one-line sequence fasta file as input to generate kmers of the sequence and reassembles it by finding the Eulierian path in the de Bruijn graph.

### Dependencies
* python v3.9.10

### Test data
Test data for evaluating the assembler
 * human corona virus - RefSeq Assembly Accession: GCF_009858895.2
 * human papilloma virus - RefSeq Assembly Accession: GCF_001274345.1
 * Human adenovirus - NCBI Reference Sequence: NC_012959.1
 * Human gammaherpesvirus 4 - NCBI Reference Sequence: NC_007605.1


### Installation
The script can be installed by following instructions:
```
git clone https://github.com/Nilsson-D/dbgAssembler.git
cd dbgAssembler
python setup.py install
```

After installation, the dbgAssembler is located in the Assembler folder

### Example runs
```
Assembler/dbgAssembler.py -h
usage: dbgAssembler.py -i <input_file> -k <kmer_size> [optional] -o <output_file> [optional]
Type -h/--help for the help message

This program takes an one-line fasta file (DNA) as input and breaks the
sequence into kmers of size k. Then reassembles the string using a de Bruijn
graph based approach

optional arguments:
  -h, --help        show this help message and exit
  -i <input file>   path to fasta file
  -k <kmer size>    kmer size (default: 31, max: 251)
  -o <output file>  name of output file, default:
                    dbgAssembler_run{current_date}.fna
  -d <directory>    name of output directory to create, default:
                    dbgAssembler_{current_date}
  -n <y/n>          if y, allow Ns in sequence

```

Running with defualt parameters

```
dbgAssembler.py -i <input_file>
```


Running with k-mer size 15
```
dbgAssembler.py -i <input_file> -k 15
```

A test run can be done by running the script test_run_viruses in the test folder.
If deciding to run the test script, keep in mind that the assembler runs for
different k-mer sizes. Each run for gammaherpesvirus 4 will take about 12 minutes (see example evaluation below).
```
test/test_run_viruses
```

If only testing the script for one size of k, the normal command can be run for one of the
test data found in the folder test_data
```
dbgAssembler.py -i test_data/<input.fna> -k 15
```


The output is a directory containing two files:
  * a fasta file for the assembly
  * a log file with the information about the input file and k-mer size


### Example result
Here are some examples how correct the assemblies are.
The assemblies are aligned against the reference genomes using mummers (v3.23) dnadiff.
The assembler has problems with lower kmers and even though a complete alignment is obtained, the assembly had to be broken up into contigs, which is not optimal (Fig 1). A higher kmer size is needed for larger genomes. This is even clear when evaluating the variants between the assembly and the reference (Fig 2). Though the dbgAssembler did not manage to complete assemble the gamma herpesvirus. Lastly, the execution time increases quickly with the sequence length as well (Fig 3). Thus, larger genomes will be time-consuming.

![](/example_figures/Aligned_seq_and_contigs.png)

Fig 1. The percentage of aligned bases to the reference genome and the number of contigs that were needed to align the assembly to reference genome evaluated by mummers dnadiff for different k-mer sizes.


![](/example_figures/variants.png)

Fig 2. The number of variants detected when aligning the assembly to the reference genome evaluated by mummers dnadiff for differ-ent k-mer sizes.

![](/example_figures/execution_time.png)

Fig 3. The execution time for each genome run using different k-mer sizes

### Author
Daniel Nilsson
