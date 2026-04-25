# Hypergraph-motif-counting

This is the official implementation accompanying

Mornie B., Colle D., Audenaert P., Pickavet M. (2026) **Efficient Counting of Higher-Order Network Motifs in Large Hypergraphs**,  
submitted to The 18th International Conference on Advances in Social Networks Analysis and Mining (ASONAM 2026).


### Installation
In order to build our application, clone the repository and run make:
```bash
git clone https://github.com/bmornie/Hypergraph-motif-counting.git
cd Hypergraph-motif-counting
make
```
This will create 2 executables: `motif_counting` and `compute_abundance`.

### Motif counting
Basic usage:
```bash
./motif_counting <input_path>
```
The input should be in hyperedge list format. Each line is a hyperedge, written as integers (nodes) separated by spaces.  
An example of a valid input file is:
```bash
0 1 5
5 7
3 7 0
```
Duplicate hyperedges and hyperedges of size >4 are ignored.  

There are 2 optional arguments:
- `-s <int>`: Motif size. Can be either 3 or 4 (default: 4)
- `-o <path>`: Output file. Default output file is `[name inputfile]_result.txt`  
  
The output is a list of frequencies of all higher-order motifs of size 3 or 4 in the input hypergraph, e.g.:
```bash
Type 1:	730438
Type 2:	2290461
...
Type 171: 70
```
The file `motif_definitions.txt` contains a representative of each type.

