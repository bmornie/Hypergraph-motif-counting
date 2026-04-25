# Hypergraph-motif-counting

This is the official implementation accompanying

Mornie B., Colle D., Audenaert P., Pickavet M. (2026) **Efficient Counting of Higher-Order Network Motifs in Large Hypergraphs**,  
submitted to The 18th International Conference on Advances in Social Networks Analysis and Mining (ASONAM 2026).


### Requirements
- GCC ≥ 13 (with C++23 support)
- GNU Make

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
- `-s <int>`: Motif size. Can be either 3 or 4 (default: 4).
- `-o <path>`: Output file. Default output file is `[name inputfile]_result.txt`.
  
The output contains a list of frequencies of all higher-order motifs of size 3 or 4 in the input hypergraph, e.g.:
```bash
Type 1:  76607
Type 2:  93241
...
Type 171:  28
```
The file `motif_definitions.txt` contains a representative of each type in a list-of-lists format.

### Compute abundance
Basic usage:
```bash
./compute_abundance <input_path>
```
Optional arguments:
- `-s <int>`: Motif size. Can be either 3 or 4 (default: 4).
- `-m <str>`: Null model. Can be either `conf` (configuration model) or `sp` (size-preserving model). Default is `conf`.
- `-n <int>`: Number of random hypergraph samples (default: 100).
- `-e <float>`: Value of $\epsilon$ used for calculating the abundance, see Eq. (2) in the paper (default: 4).
- `-t <int>`: Number of random hyperedge swaps used to construct samples under the configuration model. If not provided, this is set to 10x the number of hyperedges in the input.
- `-o <path>`: Output file. Default output file is `[name inputfile]_abundance_[null model].txt`.

The output contains a list of frequencies and corresponding abundance values for each motif type, e.g.:
```bash
Type 1:  76607  -0.0558137
Type 2:  93241  0.0702636
...
Type 171:  28  0.867414
```

### Datasets
The `data/` folder contains all datasets used in our experiments in the correct format.  

For each dataset, precomputed abundance profiles under both null models are stored in `results_correlation/`. These can be reproduced with `compute_abundance`. For example:
```bash
./compute_abundance data/hospital.txt -m sp
```
will compute the abundances (order 4 motifs) for the Hospital dataset under the size-preserving model, and store the results in `hospital_abundance_sp.txt`.  
Figure 4 from the paper can be reproduced by running the python script `data/correlation_matrix.py`.

### Dependencies
This project uses `ankerl::unordered_dense` (included in `utilities/`), for storage of hypergraphs.  
See its repository for details and license. No additional installation is required.

### License
This project is licensed under the GNU Affero General Public License v3.0 (AGPL-3.0).  
See the `LICENSE` file for details.
