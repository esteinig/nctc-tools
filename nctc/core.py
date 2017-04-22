import os
import bs4
import json
import time
import pandas
import shutil
import requests
import urllib.request

from numpy import nan, int64

from BCBio import GFF
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord

from nctc.snek import Snek
from nctc.utils import stamp, MLST, Abricate, NCTCTable

# Somewhere in report summary a SettingWithCopyWarning is thrown.
# Disabled warning for now but need to fix.
pandas.options.mode.chained_assignment = None


class NCTC3000:

    """

    Class for parsing HTML table of closed bacterial genomes from NCTC3000. Creates data table for local copy
    of all (or a subset by species) genome assemblies from NCTC3000. Generates .fasta and .gbk files, then
    analyses with mlst (https://github.com/tseemann/mlst) and abricate (https://github.com/tseemann/abricate) for
    virulence and antibiotic resistance prediction. Summarize output and update projects and species.

    """

    def __init__(self, project_path="NCTC3000", species="Staphylococcus aureus", force=False):

        self.project_path = os.path.abspath(project_path)
        self.project = os.path.basename(self.project_path)

        self.force = force

        self.date = time.strftime("%d-%m-%Y-%H-%M")
        self.species = species

        # Infrastructure

        self.species_path = os.path.join(self.project_path, self.species.replace(" ", "_"))

        self.species_gff_path = os.path.join(self.species_path, "gff")
        self.species_fasta_path = os.path.join(self.species_path, "fasta")
        self.species_genbank_path = os.path.join(self.species_path, "genbank")

        # Analysis with NCTC Programs and Snakemake

        self.nctc = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))

        self.abricate = os.path.join(self.nctc, "bin", "abricate", "bin", "abricate")

        self.env_file = os.path.join(self.nctc, "env", "nctc.json")
        self.snake_file = os.path.join(self.nctc, "pipe", "nctc.snk")
        self.config_file = os.path.join(self.nctc, "pipe", "nctc.json")

        self.mlst_path = None
        self.res_path = None
        self.vir_path = None

        self.config_run_file = None

        # NCTC3000

        self.url = "http://www.sanger.ac.uk/resources/downloads/bacteria/nctc/"

        self.tables = dict()
        self.files = dict()

        self.project_config_path = os.path.join(self.project_path, ".config")

        self.report_files = dict()

        # Initialisation Functions

        if self.force:
            print("Warning. Overwriting project", self.project_path, "in three seconds...")
            for i in range(3):
                time.sleep(1)

            shutil.rmtree(self.project_path)

        self._setup_project()

    def load_project(self):

        table_history = dict()

        files = [file for file in os.listdir(self.project_config_path)
                 if file.startswith(self.project) and file.endswith(".tab")]

        for file in files:

            components = file.replace(".tab", "").split("_")

            date = components[-2]
            project = components[0]
            table_name = components[1]
            table_type = components[-1]

            if date not in table_history.keys():
                table_history[date] = {"project": project, "name": table_name, "created": date}

            if table_type == "links":
                table_history[date]["links"] = os.path.join(self.project_config_path, file)
            else:
                table_history[date]["dataframe"] = os.path.join(self.project_config_path, file)

        try:
            last = sorted(table_history.keys(), key=lambda x: time.strptime(date, "%d-%m-%Y-%H-%M"))[-1]
        except IndexError:
            raise ValueError("There is no history to load in", self.project_config_path)

        last_table = table_history[last]

        stamp("Loading last check point from project", last_table["project"], "created at", last_table["created"])

        self.tables[last_table["name"]] = NCTCTable(project=last_table["project"], name=last_table["name"],
                                                    dataframe=last_table["dataframe"], links=last_table["links"])

    def _setup_project(self):

        for path in (self.project_path, self.species_path, self.species_gff_path,
                     self.species_fasta_path, self.species_genbank_path,
                     self.project_config_path):
            os.makedirs(path, exist_ok=True)

    def parse_website(self):

        stamp("Parsing NCTC300 website at Sanger Institute")

        site = requests.get(self.url)
        html = bs4.BeautifulSoup(site.text, "lxml")

        for table in html.find_all('table'):
            nctc_table = NCTCTable(project=self.project_path, table=table)

            self.tables[nctc_table.name] = nctc_table

            nctc_table.summarise()

            stamp("Storing table:", nctc_table.name)

            nctc_table.store(os.path.join(self.project_config_path, nctc_table))

    def parse_species(self, name="genomes", assembly="manual", strict=True, force=False):

        stamp("Fetching assemblies:", self.species)

        bacteria = self.tables[name]

        if strict:

            species_indices = bacteria.dataframe[(bacteria.dataframe["Species"] == self.species) &
                                                 (bacteria.dataframe["Chromosomes"] > 0)].index

        else:
            species_indices = bacteria.dataframe[(bacteria.dataframe["Species"] == self.species)].index

        gff_files = self._download_assemblies(self._get_links(bacteria, species_indices, assembly), force=force)

        self._gff_to_fasta(gff_files)

    @staticmethod
    def _print_parameters(parameters):

        """

        http://stackoverflow.com/questions/10724495/getting-all-arguments-and-values-passed-to-a-python-function

        :param args:
        :return:
        """

        parameters.pop("self")

        for key, value in sorted(parameters.items()):
            stamp("{key} = {value}".format(key=key, value=value))

    def summarise(self, out_path, user_path=None, merge="outer", resistance_db="resfinder", virulence_db="vfdb",
                  cnv=True, cnv_mincov=95, single=False, csv=True, force=False, verbose=True):

        out_path = os.path.abspath(out_path)
        
        if verbose:
            self._print_parameters(locals())

        if os.path.exists(out_path) and force:
            shutil.rmtree(out_path)

        os.makedirs(out_path)

        self._get_analysis_paths(user_path)

        frames = self._parse_reports(resistance_db=resistance_db, virulence_db=virulence_db)

        cleaned_frames = self._clean_reports(frames)
        merged_frames = self._merge_frames(cleaned_frames, merge=merge)

        if cnv:
            cnv_frames = self._get_copy_numbers(cleaned_frames, cnv_mincov)
            cnv_merged = self._merge_frames(cnv_frames, left=cleaned_frames["mlst"], merge=merge)

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

    @staticmethod
    def _merge_frames(frames, left=None, merge="outer"):

        merged = {}
        for output in ("virulence", "resistance"):
            if left is None:
                merged_df = pandas.merge(left=frames["mlst"], right=frames[output], left_on="ID", right_on="ID",
                                         how=merge)
            else:
                merged_df = pandas.merge(left=left, right=frames[output], left_on="ID", right_on="ID",
                                         how=merge)

            merged[output] = merged_df

        return merged

    def _get_copy_numbers(self, frames, cnv_mincov):

        cnv_frames = {}

        for output, dataframe in frames.items():
            if output in ("virulence", "resistance"):
                cnv_frame = {}
                for column in dataframe.columns:
                    if column not in ("ID", "Count"):
                        cnv_frame[column] = self._get_cnv_column(output, dataframe, column, cnv_mincov)
                    else:
                        cnv_frame[column] = dataframe[column]

                cnv_frames[output + "_cnv"] = pandas.DataFrame(cnv_frame)

        # Return dictionary without key "mlst", just copy numebr variation for Abricate
        return cnv_frames

    @staticmethod
    def _get_cnv_column(output, dataframe, column, cnv_mincov):

        cnv_column = []
        for entry in dataframe[column].str.split(";"):
            if len(entry) == 1 and entry[0] == ".":
                cnv_column.append(nan)
            else:
                cov_filtered = []
                for e in entry:
                    try:
                        e = float(e)
                        if e >= cnv_mincov:
                            cov_filtered.append(e)
                    except ValueError:
                        print("Could not convert coverage in dataframe for", output, "to float.")

                cnv_column.append(len(cov_filtered))

        return cnv_column

    def _clean_reports(self, frames):

        clean = {}

        for output, frame in frames.items():
            ids = frame.ix[:, 0]
            cleaned = frame[ids.str.endswith(".fasta")]
            cleaned_ids = cleaned.iloc[:, 0].map(lambda x: self._extract_name(x))
            cleaned = self._clean_output_reports(output, cleaned)
            cleaned.insert(0, "ID", cleaned_ids)
            cleaned.reindex()

            clean[output] = cleaned

        return clean

    def _clean_output_reports(self, output, dataframe):
        if output in ("virulence", "resistance"):
            keep = ~dataframe.columns.str.contains("|".join(['Unnamed:', '#FILE']))
            cleaned = dataframe.iloc[:, keep]
            cleaned.rename(columns={'NUM_FOUND': 'Count'}, inplace=True)

        elif output == "mlst":
            cleaned = dataframe.iloc[:, 1:]
            cleaned.rename(columns={1: 'Species', 2: 'MLST'}, inplace=True)
            cleaned.rename(columns=self._get_mlst_header(cleaned), inplace=True)
        else:
            cleaned = None

        return cleaned

    @staticmethod
    def _get_mlst_header(df):

        allele_cols = {}
        i = 1
        for col in df.columns:
            if type(col) is not int64:
                allele_cols[col] = col
            else:
                allele_cols[col] = "MLST_Allele_" + str(i)
                i += 1

        return allele_cols

    @staticmethod
    def _extract_name(entry):

        file_name = os.path.basename(entry)
        name, extension = os.path.splitext(file_name)

        return name

    def analyse(self, user_path=None, resistance_db="resfinder", virulence_db="vfdb", minid=90, mincov=95,
                cluster=True, mlst=True, resistance=True, virulence=True, force=False, verbose=True):

        """

        Analyse a collection of species reference genomes using

            mlst: https://github.com/tseemann/mlst
            abricate: https://github.com/tseemann/abricate

        If executed in a cluster environment the function uses Snakemake and the wrapper class Snek.

        Recommend local runs as they are quick enough for single species, cluster if all are required.

        """

        if verbose:
            self._print_parameters(locals())

        analysis_path, fasta_path = self._get_analysis_paths(user_path)

        stamp("Analysing data for", self.species)

        if cluster:
            self._modify_config(virulence_db, resistance_db, minid, mincov, force)

            snek = Snek(path=analysis_path, snake_file=self.snake_file,
                        config_file=self.config_run_file, executable="snakemake", force=force)

            snek.run()

        else:

            if mlst:
                mlst = MLST(target_path=fasta_path, out_dir=self.mlst_path, exec_path="mlst",
                            force=force)
                mlst.run(min_id=minid, min_cov=mincov)
            else:
                stamp("MLST results exist at", self.mlst_path)

            if resistance:
                res = Abricate(target_path=fasta_path, out_dir=self.res_path, exec_path=self.abricate,
                               force=force)
                res.run(db=resistance_db, min_id=minid, ext=".fasta")
            else:
                stamp("Resfinder results exist at", self.res_path)

            if virulence:
                vir = Abricate(target_path=fasta_path, out_dir=self.vir_path, exec_path=self.abricate,
                               force=force)
                vir.run(db=virulence_db, min_id=minid, ext=".fasta")
            else:
                stamp("Virulence results exist at", self.vir_path)

    def _get_analysis_paths(self, user_path):

        """

        Called by analyse and summarise to define paths available for analysis and paths containing analysis reports.
        If user_path is not defined, return path to species and fasta directories of a project.

        :param user_path: path to analysis directory, requires files at user_path/fasta
        :return: user_path, fasta_path

        """

        if user_path is not None:
            analysis_path = user_path
            fasta_path = os.path.join(user_path, "fasta")
            stamp("Inititating local analysis paths at user path:", analysis_path)
        else:
            analysis_path = self.species_path
            fasta_path = self.species_fasta_path
            stamp("Inititating project analysis paths at species path:", self.species_path)

        # Generate paths to analysis results, make available for parsing report functions:

        self.mlst_path = os.path.join(analysis_path, "mlst")
        self.res_path = os.path.join(analysis_path, "res")
        self.vir_path = os.path.join(analysis_path, "vir")

        # Return user or species analysis and fasta paths

        return analysis_path, fasta_path

    def _parse_reports(self, resistance_db="resfinder", virulence_db="vfdb"):

        files = {"mlst": os.path.join(self.mlst_path, "mlst.report"),
                 "resistance": os.path.join(self.res_path, resistance_db + ".report"),
                 "virulence": os.path.join(self.vir_path, virulence_db + ".report")}

        frames = {"mlst": nan,
                  "resistance": nan,
                  "virulence": nan}

        for output, file in files.items():
            if os.path.exists(file):
                stamp("Found output file", file)
                if output == "mlst":
                    frames[output] = pandas.read_csv(file, sep="\t", header=None)
                else:
                    frames[output] = pandas.read_csv(file, sep="\t")
            else:
                stamp("Could not find output file, skipping", file)

        return frames

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

    def _gff_to_fasta(self, gff_files):

        for file in gff_files:

            file_name, ext = os.path.splitext(file)
            base_name = os.path.basename(file_name)

            fasta = os.path.join(self.species_fasta_path, base_name + ".fasta")
            genbank = os.path.join(self.species_genbank_path, base_name + ".gbk")

            if not os.path.exists(fasta) or not os.path.exists(genbank):

                records = self._parse_gff(file, file_name)

                if not os.path.exists(fasta):
                    SeqIO.write(records, fasta, format="fasta")
                if not not os.path.exists(genbank):
                    SeqIO.write(records, genbank, format="genbank")

            self.files[base_name] = {"gff": os.path.realpath(file),
                                     "fasta": os.path.realpath(fasta),
                                     "genbank": os.path.realpath(genbank)}

    @staticmethod
    def _parse_gff(file, file_name):

        records = []

        with open(file, "r") as infile:

            for i, rec in enumerate(GFF.parse(infile)):

                # Enumerates the contigs (can be chromosome, plasmid and unidentified)
                # based on total number of contigs (not type)
                rec_id = rec.id + "_" + str(i + 1)

                if len(rec_id) > 15:
                    rec_id = "contig_" + "_" + str(i + 1)

                seq_record = SeqRecord(Seq(str(rec.seq), IUPAC.unambiguous_dna), id=rec_id,
                                       description=os.path.basename(file_name),
                                       features=rec.features)

                records.append(seq_record)

        return records

    def _get_links(self, bacteria, species_indices, assembly):

        links = bacteria.links.ix[species_indices, assembly + "_gff"].dropna()

        file_names = self._get_file_names(bacteria.dataframe.ix[links.index])

        stamp("There are", len(links), "valid entries for", self.species, "on the website for NCTC3000")

        files = pandas.concat([file_names, links], axis=1)
        files.columns = ["File", "Link"]

        return files

    def _download_assemblies(self, links, force=False):

        exist = 0
        downloaded = 0
        out_paths = []
        for entry in links.itertuples():
            file, url = entry.File, entry.Link

            out_path = os.path.join(self.species_gff_path, file)
            out_paths.append(out_path)

            if os.path.exists(out_path) and not force:
                exist += 1
            else:
                stamp("Downloading", file, "from", url)
                urllib.request.urlretrieve(url, out_path)
                downloaded += 1

        stamp("Found", exist, "files and downloaded", downloaded, "new assemblies for",
              self.species, "in project", self.project)

        return out_paths

    @staticmethod
    def _get_file_names(link_entries):

        species_names = link_entries["Species"].str.replace(" ", "_")
        strain_names = species_names.str.cat([link_entries["Strain"]], sep="_")

        return strain_names.str.cat([".gff" for _ in range(len(link_entries))])
