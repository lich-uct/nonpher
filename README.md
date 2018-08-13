# Nonpher
This Python package contains functions for generating hard-to-synthesize (HS) structures as well as several functions for
calculating molecular complexity. Nonpher utilizes morphing algorithm as is implemented by Molpher-lib, which is based on
the cheminformatic toolkit RDKit.

A starting structure is randomly, step-by-step transformed to more complex one (usually but not necessarily as it is a
random process). HS structures can be useful when someone wants to classify them, because they are not easily accessible
from any database.
 
## Instalation
### Prerequisities
#### Supported platforms:
* Linux 64-bit (as molpher-lib is compiled only for it, but otherwise the package is platform independent)

#### Dependencies
* RDKit
* Molpher-lib >=0.0.0b2 (now it works with RDKit 2018.3.1 and libboost 1.65.1)

### Installation with Anaconda
Nonpher is distributed as a conda package. At the moment, this is the preferred way to install and use the library.
All you need to do is get the full Anaconda distribution or its lightweight variant, Miniconda. It is essentially a
Python distribution, package manager and virtual environment in one and makes setting up a development environment
for any project very easy. After installing Anaconda/Miniconda (and environment preparing) you can run the following in
the Linux terminal:
```bash
conda install -c rdkit -c lich nonpher
```

### Installation with setup.py
Anyway you have installed RDKit and Molpher-lib, you can download/clone Nonpher and install it 
from its directory with:
```bash
python setup.py install
```

## Quick start
To generate HS structures, Nonpher requires starting structures, which should be relatively easy-to-synthesize (ES)
and should somehow be standardized, at least removed ions and charges. To quick start, you can use Nonpher with prepared
simple script which only needs an input (stdin or file) with structures in csv format. First column is for structure IDs,
the second for structures themselves encoded in SMILES, all other columns are skipped. If you use input file, you can also
specify output file which will be in the same format: ID, SOURCE_STRUCTURE and generated HS structure.
To run the script just write: 

```bash
python -m nonpher.nonpher [-H] [INPUT_FILE [OUTPUT_FILE]]
```

where parameter -H is for skipping the first line (header)

