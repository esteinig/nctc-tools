# nctc-tools

Pipeline facilitating access to complete bacterial reference assemblies from public repositories:
* maintain complete and curated reference genomes (PacBio) from [NCTC3000](http://www.sanger.ac.uk/resources/downloads/bacteria/nctc/)
* type local assemblies or species with [mlst](https://github.com/tseemann/mlst) and [Abricate](https://github.com/tseemann/abricate)
* collect analysis summaries for local assemblies or species from NCTC3000

## Data Usage

Please refer to data usage guidelines from [NCTC3000](http://www.sanger.ac.uk/resources/downloads/bacteria/nctc/):

>This sequencing centre plans on publishing the completed and annotated sequences in a peer-reviewed journal as soon as possible. Permission of the principal investigator should be obtained before publishing analyses of the sequence/open reading frames/genes on a chromosome or genome scale. See our data sharing policy.

>**Please note**: these are pre-submission assemblies that should not be treated as final versions. Assemblies contain both chromosomal and plasmid contigs.

## Setup

Conda:

```
pip install conda
```

Clone this repository recursively to include latest version of [Abricate](https://github.com/tseemann/abricate):

```
git clone --recursive https://github.com/esteinig/nctc-tools
```

Create environment with dependencies:

```
conda env create -f nctc-tools/env/nctc_tools.yaml
```

Make `./nctc-tools/nctc_tools.py` executable and export to PATH.

## Run

Activate environment before executing tasks:

```
source activate nctc-tools
```

Print help message for global options and tasks:

```
nctc_tools.py --help

nctc_tools.py make --help
nctc_tools.py type --help
nctc_tools.py collect --help
```

## Quick Start

```
source activate nctc

nctc_tools.py --path ./ref_db --species "Escherichia coli" make --chromosomes

cd ./ref_db

nctc_tools.py --species "Escherichia coli" type -v vfdb -r resfinder
nctc_tools.py --species "Escherichia coli" collect --cnv --csv -o ../reports

# snakemake at nctc_tools/pipes
nctc_tools.py --species "Escherichia coli" type --cluster

cd ..

# fasta files at path, not project
nctc_tools.py --path ./mrsa type --minid 90 

source deactivate
```

## Tasks

[`make`](): create species repository on which to operate

[`type`](): perform general typing methods for resistance, virulence, mlst

[`collect`](): summarise typing outputs

[`update`](): update species repository and summarise updates

## Contact

eikejoachim.steinig@my.jcu.edu.au
