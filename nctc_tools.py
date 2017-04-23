#!/usr/bin/env python

"""

nctc_tools.py

Establishes NCTC3000 library by parsing the NCTC3000 website, updating a local data structure, downloading
assemblies, performing genotyping, resistance / virulence factor detection and pangenome analysis on the
latest sequences.

@NCTC3000

http://www.sanger.ac.uk/resources/downloads/bacteria/nctc/

Install:

% git clone --recursive https://github.com/esteinig/nctc_tools

% conda env create -f ./nctc_tools/envs/nctc.yaml

Run:

% source activate nctc

% nctc_tools.py --project ./nctc_db --species "Staphylococcus aureus" make --chromosomes
% nctc_tools.py --project ./nctc_db --species "Staphylococcus aureus" type --minid 90
% nctc_tools.py --project ./nctc_db --species "Staphylococcus aureus" collect -o reports --cnv

% nctc_tools.py --project ./nctc_db --species "Staphylococcus aureus" update

% # Cluster with Snakemake (might take longer than local runs depending on queue)

% nctc_tools.py --project ./nctc_db --species "Staphylococcus aureus" type --cluster

% # local type and collect where files at user_path/*.fasta:

% nctc_tools.py --user_path ./mrsa type --minid 90 --resistance_db resfinder

"""

from nctc.core import NCTC3000
from nctc.utils import CommandLine


def main():

    cmd_line = CommandLine()

    args = vars(cmd_line.args)
    nctc = NCTC3000(project_path=args["project"], species=args["species"], force=args["overwrite"],
                    user=args["user_path"])

    if args["subparser"] == "make":

        nctc.parse_website()
        nctc.parse_species(strict=args["chromosomes"], force=args["force"])

    if args["subparser"] == "type":

        if not args["user_path"]:
            nctc.load_project()

        nctc.analyse(user_path=args["user_path"], cluster=args["cluster"], mlst=args["mlst"],
                     resistance=args["resistance"], resistance_db=args["resistance_db"], virulence=args["virulence"],
                     virulence_db=args["virulence_db"], minid=args["minid"], mincov=args["mincov"],
                     force=args["force"])

    if args["subparser"] == "collect":

        if not args["user_path"]:
            nctc.load_project()

        nctc.summarise(out_path=args["output"], user_path=args["user_path"], merge=args["merge"], cnv=args["cnv"],
                       cnv_mincov=args["cnv_mincov"], force=args["force"], resistance_db=args["resistance_db"],
                       virulence_db=args["virulence_db"], csv=args["csv"])

    if args["subparser"] == "update":

        nctc.load_project()

main()



