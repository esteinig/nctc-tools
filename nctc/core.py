import os
import bs4
import json
import time
import pandas
import shutil
import requests
import urllib.request

from numpy import nan

from BCBio import GFF
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord

from nctc.utils import *
from nctc.snek import Snek

# Somewhere in report summary a SettingWithCopyWarning is thrown.
# Disabled warning for now, but need to fix.
pandas.options.mode.chained_assignment = None

FASTA_EXTENSIONS = (".fa", ".fasta")


class NCTC3000:

    """

    Class for parsing HTML table of closed bacterial genomes from NCTC3000. Creates data table for local copy
    of all (or a subset by species) genome assemblies from NCTC3000. Generates .fasta and .gbk files, then
    analyses with mlst (https://github.com/tseemann/mlst) and abricate (https://github.com/tseemann/abricate) for
    virulence and antibiotic resistance prediction. Summarize output and update projects and species.

    """

    def __init__(self, path=os.getcwd(), species="Staphylococcus aureus", force=False, verbose=True):

        self.date = time.strftime("%d-%m-%Y-%H-%M")

        self.path = os.path.abspath(path)
        self.project = True

        self.species = species
        self.species_name = self.species.replace(" ", "_").lower()
        self.species_path = os.path.join(self.path, self.species_name)

        self.force = force
        self.verbose = verbose

        # Infrastructure

        self.config_path = os.path.join(self.path, ".nctc")

        # Analysis with NCTC Programs and Snakemake

        self.nctc = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))

        self.abricate = os.path.join(self.nctc, "bin", "abricate", "bin", "abricate")

        self.env_file = os.path.join(self.nctc, "env", "nctc.json")
        self.snake_file = os.path.join(self.nctc, "pipe", "nctc.snk")
        self.config_file = os.path.join(self.nctc, "pipe", "nctc.json")

        self.config_run_file = None

        # HTML NCTC3000

        self.url = "http://www.sanger.ac.uk/resources/downloads/bacteria/nctc/"

        self.tables = dict()

        self.files = dict()
        self.report_files = dict()

    def _load_project(self):

        files = [file for file in os.listdir(self.config_path) if file.endswith(".tab")]

        table_history = self._load_history(files)

        try:
            last = sorted(table_history.keys(), key=lambda x: time.strptime(x, "%d-%m-%Y-%H-%M"))[-1]
        except IndexError:
            raise ValueError("There is no history to load in", self.config_path)

        last_table = table_history[last]

        stamp("Loading last check point for table", last_table["name"],
              "from project created at", last_table["created"])

        self.tables[last_table["name"]] = NCTCTable(name=last_table["name"], dataframe=last_table["dataframe"],
                                                    links=last_table["links"])

    def _load_history(self, files):

        table_history = dict()

        for file in files:
            components = file.replace(".tab", "").split("_")

            date = components[-2]
            table_name = components[0]
            table_type = components[-1]

            if date not in table_history.keys():
                table_history[date] = {"name": table_name, "created": date}

            if table_type == "links":
                table_history[date]["links"] = os.path.join(self.config_path, file)
            else:
                table_history[date]["dataframe"] = os.path.join(self.config_path, file)

        return table_history

    def _setup_project(self):

        """

        Called upon initialisation of class NCTC3000

        1. If self.path exists and contains path/.nctc, create a new project directory.

        2. If self.path exists and does not contain path/.nctc, it is a user directory.

        3. If self.path does not exist, create new project directory.


        """

        if os.path.exists(self.path) and os.path.exists(self.config_path):
            stamp("Operating on project at", self.path)
            self._load_project()
        elif os.path.exists(self.path) and not os.path.exists(self.config_path):
            stamp("Could not find project at", self.path)
            exit(1)
        elif not os.path.exists(self.path):
            stamp("Creating project at", self.path)
            self._create_project()
        else:
            stamp("Hmmmmm...")

    def make(self, strict=True, force=False):

        """

        Main function for creating a new project and maintaining it.

        """

        if self.verbose:
            self._print_parameters(locals())

        self._setup_project()

        self.parse_website()
        self.parse_species(strict=strict, force=force)

    def _create_project(self):

        gff_path = os.path.join(self.species_path, "gff")
        fasta_path = os.path.join(self.species_path, "fasta")
        genbank_path = os.path.join(self.species_path, "genbank")

        for path in (self.path, self.species_path, self.config_path,
                     gff_path, fasta_path, genbank_path):
            os.makedirs(path)

    def parse_website(self):

        stamp("Parsing NCTC300 website at Sanger Institute")

        site = requests.get(self.url)
        html = bs4.BeautifulSoup(site.text, "lxml")

        for table in html.find_all('table'):
            nctc_table = NCTCTable(name="nctc_assemblies", table=table)

            self.tables[nctc_table.name] = nctc_table

            nctc_table.summarise()

            stamp("Storing table", nctc_table.name, "from date", self.date, "at", self.config_path)

            nctc_table.store(self.config_path)

    def parse_species(self, name="genomes", assembly="manual", strict=True, force=False):

        stamp("Fetching assemblies for", self.species, "from NCTC3000")

        bacteria = self.tables[name]

        if strict:
            species_indices = bacteria.dataframe[(bacteria.dataframe["Species"] == self.species) &
                                                 (bacteria.dataframe["Chromosomes"] > 0)].index
        else:
            species_indices = bacteria.dataframe[(bacteria.dataframe["Species"] == self.species)].index

        gff_files = download_assemblies(gff_path=os.path.join(self.species_path, "gff"),
                                        links=get_links(bacteria, species_indices, assembly), force=force)

        gff_to_fasta(gff_files, fasta_path=os.path.join(self.species_path, "fasta"),
                     genbank_path=os.path.join(self.species_path, "genbank"), verbose=self.verbose)

    @staticmethod
    def _print_parameters(parameters):

        """

        http://stackoverflow.com/questions/10724495/getting-all-arguments-and-values-passed-to-a-python-function

        :param parameters: within function call to local()
        :return:

        """

        parameters.pop("self")

        for key, value in sorted(parameters.items()):
            stamp("{key} = {value}".format(key=key, value=value))

    def collect(self, out_path, user_path=None, merge="outer", resistance_db="resfinder", virulence_db="vfdb",
                  cnv=True, cnv_mincov=95, single=False, csv=True, force=False):

        out_path = os.path.abspath(out_path)
        
        if self.verbose:
            self._print_parameters(locals())

        if os.path.exists(out_path) and force:
            shutil.rmtree(out_path)

        os.makedirs(out_path)

        _, _, mlst_path, res_path, vir_path = self._get_analysis_paths(resistance_db, virulence_db, mode="collect")

        frames = parse_reports(mlst_path, res_path, vir_path, resistance_db=resistance_db, virulence_db=virulence_db)

        cleaned_frames = clean_reports(frames)
        merged_frames = merge_frames(cleaned_frames, merge=merge)

        if cnv:
            cnv_frames = get_copy_numbers(cleaned_frames, cnv_mincov)
            cnv_merged = merge_frames(cnv_frames, left=cleaned_frames["mlst"], merge=merge)

            cleaned_frames.update(cnv_frames)
            merged_frames.update(cnv_merged)

        if single:
            write_frames = cleaned_frames
        else:
            write_frames = merged_frames

        for output, frame in write_frames.items():
            if csv:
                sep = ","
                file_path = os.path.join(out_path, output + ".csv")
            else:
                sep = "\t"
                file_path = os.path.join(out_path, output + ".tab")

            frame.to_csv(file_path, sep=sep, na_rep=".")
            stamp("Written summary to file:", file_path)

    def type(self, resistance_db="resfinder", virulence_db="vfdb", minid=90, mincov=95,
             cluster=True, mlst=True, resistance=True, virulence=True, force=False):

        """

        Analyse a collection of species reference genomes using

            mlst: https://github.com/tseemann/mlst
            abricate: https://github.com/tseemann/abricate

        If executed in a cluster environment the function uses Snakemake and the wrapper class Snek.

        Recommend local runs as they are quick enough for single species, cluster if all are required.

        """

        if self.verbose:
            self._print_parameters(locals())

        analysis_path, fasta_path, mlst_path, res_path, vir_path = self._get_analysis_paths(resistance_db,
                                                                                            virulence_db, mode="type")

        if force:
            stamp("May the force be with you.")

        stamp("Running analyses on files (.fasta) in", fasta_path)

        if cluster:

            self._modify_config(virulence_db, resistance_db, minid, mincov, force)

            snek = Snek(path=analysis_path, snake_file=self.snake_file,
                        config_file=self.config_run_file, executable="snakemake", force=force)

            snek.run()

        else:

            if mlst:
                mlst = MLST(target_path=fasta_path, out_dir=mlst_path, exec_path="mlst", force=force)

                mlst.run(min_id=minid, min_cov=mincov)

            if resistance:
                res = Abricate(target_path=fasta_path, out_dir=res_path, exec_path=self.abricate, force=force)
                res.run(db=resistance_db, min_id=minid)

            if virulence:
                vir = Abricate(target_path=fasta_path, out_dir=vir_path, exec_path=self.abricate, force=force)
                vir.run(db=virulence_db, min_id=minid)

    def _get_analysis_paths(self, res_db, vir_db, mode="type"):

        """

        Called by analyse and summarise to define paths available for analysis and paths containing analysis reports.
        If user_path is not defined, return path to species and fasta directories of a project.

        :param user_path: path to analysis directory, requires files at user_path/fasta
        :return: user_path, fasta_path

        """

        analysis_path = None
        fasta_path = None

        if os.path.exists(self.path) and os.path.exists(self.config_path):
            stamp("Operating on species", self.species, "in project", self.path)
            analysis_path = os.path.join(self.species_path, "analysis")
            fasta_path = os.path.join(self.species_path, "fasta")

        elif os.path.exists(self.path) and not os.path.exists(self.config_path):
            stamp("Could not find project at", self.path, "- assuming user path.")

            analysis_path = os.path.join(self.path, "analysis")
            fasta_path = self.path

            if mode == "type":
                fasta_files = [file for file in os.listdir(self.path) if file.endswith(FASTA_EXTENSIONS)]

                if len(fasta_files) == 0:
                    stamp("Could not detect any files (.fasta, .fa) in path", self.path)
                    exit(1)
                else:
                    stamp('Found files (.fasta, .fa) in path', self.path)
            elif mode == "collect":
                if not os.path.exists(analysis_path):
                    stamp("Could not find analysis directory at", analysis_path)
                    stamp("Perhaps you need to run typing first?")
                    exit(1)

        elif not os.path.exists(self.path):
            stamp("Could not find path", self.path)
            exit(1)
        else:
            stamp("Hmmmmm...")
            exit(1)

        # Generate paths to analysis results

        mlst_path = os.path.join(analysis_path, "mlst")
        res_path = os.path.join(analysis_path, "resistance", res_db)
        vir_path = os.path.join(analysis_path, "virulence", vir_db)

        # Return user or species analysis and fasta paths

        return analysis_path, fasta_path, mlst_path, res_path, vir_path

    def _modify_config(self, vf_db, res_db, min_id, min_cov, force):

        with open(self.config_file, "r") as infile:
            config = json.load(infile)

        config["snek"]["force"] = force
        config["abricate_path"] = self.abricate
        config["env"] = self.env_file
        config["abricate_vir_db"] = vf_db
        config["abricate_res_db"] = res_db
        config["minid"] = min_id
        config["mlst_mincov"] = min_cov

        self.config_run_file = os.path.join(self.project_config_path, "run_" + os.path.basename(self.config_file))

        with open(self.config_run_file, "w") as outfile:
            json.dump(config, outfile)