# nctc-tools

Pipeline facilitating access to complete bacterial reference assemblies from public repositories:
* maintain complete and curated reference genomes (PacBio) from [NCTC3000](http://www.sanger.ac.uk/resources/downloads/bacteria/nctc/)
* type local assemblies or species with [mlst](https://github.com/tseemann/mlst) and [Abricate](https://github.com/tseemann/abricate)
* collect analysis summaries for local assemblies or species from NCTC3000

## Data Usage

Please refer to data usage guidelines from [NCTC3000](http://www.sanger.ac.uk/resources/downloads/bacteria/nctc/):

>This sequencing centre plans on publishing the completed and annotated sequences in a peer-reviewed journal as soon as possible. Permission of the principal investigator should be obtained before publishing analyses of the sequence/open reading frames/genes on a chromosome or genome scale. See our data sharing policy.

Assemblies will be able to be updated with task [`update`]()

>**Please note**: these are pre-submission assemblies that should not be treated as final versions. Assemblies contain both chromosomal and plasmid contigs.

## Setup

Dependencies:

* [Anaconda](https://www.continuum.io/DOWNLOADS) for Python 3

Clone this repository recursively to include latest version of [Abricate](https://github.com/tseemann/abricate):

```
git clone --recursive https://github.com/esteinig/nctc-tools
```


Next update will include installation via `conda`, for now you can use `install.sh`:

**Script**:

```
chmod +x ./nctc-tools/install.sh
bash ./nctc-tools/install.sh
```

**Manual**:

Create environment with dependencies:

```
conda env create -f nctc-tools/env/nctc_tools.yaml
```


Make script executable and export to PATH, then:

**Run**

Activate environment before executing tasks:

```
source activate nctc-tools
```

Print help message for global options and task options:

```
nctc_tools.py --help

nctc_tools.py make --help
nctc_tools.py type --help
nctc_tools.py collect --help
```

## Quick Start

```
source activate nctc

# make project at --path
nctc_tools.py --path ./ref_db --species "Escherichia coli" make --chromosomes

# specify full path via --path or cd into project
cd ./ref_db

# run tasks within --path
nctc_tools.py --species "Escherichia coli" type -v vfdb -r resfinder
nctc_tools.py --species "Escherichia coli" collect --cnv --csv -o ../reports

# snakemake at nctc_tools/pipes
nctc_tools.py --species "Escherichia coli" type --cluster

# return to parent dir
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
