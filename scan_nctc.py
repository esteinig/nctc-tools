#!/usr/bin/env python

"""

scanNCTC

Establishes NCTC3000 library by parsing the NCTC3000 website, updating a local data structure, downloading
assemblies, performing genotyping, resistance / virulence factor detection and pangenome analysis on the
latest sequences.

@NCTC3000

http://www.sanger.ac.uk/resources/downloads/bacteria/nctc/

Install:

% git clone https://gitbub.com/esteinig/scanNCTC
% git clone https://github.com/tseemann/abricate ./scanNCTC/bin/

% conda env create -f ./scanNCTC/envs/nctc.yaml

Run:

% source activate nctc

% nctc make --project ./nctc_db --species staphylococcus_aureus --db nctc
% nctc analysis --species staphylococcus_aureus --min_id 90 --db resfinder,vfdb

% nctc summary --species staphylococcus_aureus --report ./saureus_report.tab
% nctc update --project ./nctc_db --species staphylococcus_aureus

Cluster with Snakemake (might take longer depending on queue on HPC)

% nctc analyse --species staphylococcus_aureus --cluster --background


"""

from nctc.core import NCTC3000
from nctc.utils import CommandLine


def main():

    cmd_line = CommandLine()

    args = vars(cmd_line.args)

    nctc = NCTC3000(project_path=args["project"], species=args["species"], force=args["overwrite"])

    if args["subparser"] == "make":

        nctc.parse_website()
        nctc.parse_species(strict=args["chromosomes"], force=args["force"])

    if args["subparser"] == "type":

        nctc.load_project()

        nctc.analyse(user_path=args["user_path"], cluster=args["cluster"], mlst=args["mlst"],
                     resistance=args["resistance"], resistance_db=args["resistance_db"], virulence=args["virulence"],
                     virulence_db=args["virulence_db"], minid=args["minid"], mincov=args["mincov"],
                     force=args["force"])

    if args["subparser"] == "collect":

        nctc.load_project()

        nctc.summarise(out_path=args["output"], user_path=args["user_path"], merge=args["merge"], cnv=args["cnv"],
                       cnv_mincov=args["cnv_mincov"], force=args["force"], resistance_db=args["resistance_db"],
                       virulence_db=args["virulence_db"], csv=args["csv"])

    if args["subparser"] == "update":

        nctc.load_project()

main()



