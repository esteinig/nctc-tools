"""

NCTC Utils:
    -- Command Line Parser
    -- NCTC Table for parsing HTML
    -- NCTC Programs for local Abricate / MLST / Roary

"""

import bs4
import json
import urllib

from numpy import nan

from BCBio import GFF
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord

import os
import time
import pandas
import shutil
import textwrap
import argparse

from numpy import nan, int64

from subprocess import call
from subprocess import CalledProcessError

FASTA_EXTENSIONS = (".fa", ".fasta")


class CommandLine:

    def __init__(self):

        parser = argparse.ArgumentParser()

        parser.add_argument("--version", "-v", default=False, action="store_true",
                            help="print version and exit")

        parser.add_argument("--path", "-p", type=str, required=False,
                            default=os.getcwd(), dest="path", help="path to project or fasta directory")

        parser.add_argument("--species", "-s", type=str, required=False,
                            default="Staphylococcus aureus", dest="species", help="species name")

        parser.add_argument("--verbose", "-d", default=False, dest="verbose", action="store_true",
                            help="print details")

        subparsers = parser.add_subparsers(help='Command-line interface for scanNCTC')

        make_parser = subparsers.add_parser("make", help="create local repository of NCTC3000")

        make_parser.add_argument("--force", "-f", default=False, dest="force", action="store_true",
                                 help="force download of assemblies")
        make_parser.add_argument("--chromosomes", "-c",  default=False, dest="chromosomes",
                                 action="store_true", help="only manual assemblies with identified chromosomes")

        make_parser.set_defaults(subparser='make')

        analysis_parser = subparsers.add_parser("type", help="type assemblies")

        analysis_parser.add_argument("--cluster", default=False, dest="cluster", action="store_true",
                                     help="enable cluster execution")
        analysis_parser.add_argument("--resistance_db", "-r", type=str, default="resfinder", dest="resistance_db",
                                     help="resistance database")
        analysis_parser.add_argument("--virulence_db", "-v", type=str, default="vfdb", dest="virulence_db",
                                     help="virulence database")
        analysis_parser.add_argument("--mincov", type=float, default=90, dest="mincov",
                                     help="minimum coverage for mlst")
        analysis_parser.add_argument("--minid", type=float, default=90, dest="minid",
                                     help="minimum coverage for abricate and mlst")
        analysis_parser.add_argument("--no-mlst", default=True, dest="mlst", action="store_false",
                                     help="disable mlst")
        analysis_parser.add_argument("--no-resistance", default=True, dest="resistance", action="store_false",
                                     help="disable resistance")
        analysis_parser.add_argument("--no-virulence", default=True, dest="virulence", action="store_false",
                                     help="disable virulence")
        analysis_parser.add_argument("--force", "-f", default=False, dest="force", action="store_true",
                                     help="force analysis of assemblies")

        analysis_parser.set_defaults(subparser='type')

        collect_parser = subparsers.add_parser("collect", help="collect results from analyses")

        collect_parser.add_argument("--output", "-o", type=str, default="report.tab", required=True,
                                    dest="output", help="summary output file")
        collect_parser.add_argument("--merge", "-m", type=str, default="outer", required=False,
                                    dest="merge", help="merge on mlst and abricate reports")
        collect_parser.add_argument("--cnv", default=False, dest="cnv", action="store_true",
                                    help="summarise copy number variation in abricate reports")
        collect_parser.add_argument("--cnv_mincov", type=float, default=95, dest="cnv_mincov", required=False,
                                    help="minimum coverage to consider copy in copy number variation")
        collect_parser.add_argument("--single", default=False, dest="single", action="store_true",
                                    help="write unmerged reports")
        collect_parser.add_argument("--csv", default=False, dest="csv", action="store_true",
                                    help="write reports in csv")
        collect_parser.add_argument("--force", "-f", default=False, dest="force", action="store_true",
                                    help="force collection of reports")
        collect_parser.add_argument("--resistance_db", "-r", type=str, default="resfinder", dest="resistance_db",
                                    help="resistance database used in analysis")
        collect_parser.add_argument("--virulence_db", "-v", type=str, default="vfdb", dest="virulence_db",
                                    help="virulence database used in analysis")

        collect_parser.set_defaults(subparser='collect')

        pangenome_parser = subparsers.add_parser("pangenome")

        pangenome_parser = subparsers.add_parser("pangenome")

        self.args = parser.parse_args()


class NCTCTable:

    """ Class parsing a HTML Soup table into a Pandas dataframe """

    def __init__(self, name, table=None, dataframe=pandas.DataFrame(), links=pandas.DataFrame()):

        """

        Initialize class NCTCTable.

        Calls:

        _transform()

        :param table: HTML soup string containing table section of NCTC3000.

        """

        self.name = name
        self.table = table

        self.dataframe = dataframe
        self.links = links

        if self.table is not None:
            self._transform()

    def store(self, path):

        filename = os.path.join(path, self.name + "_" + time.strftime("%d-%m-%Y-%H-%M"))

        filename_table = filename + "_data.tab"
        filename_links = filename + "_links.tab"

        self.dataframe.to_csv(filename_table, sep="\t")
        self.links.to_csv(filename_links, sep="\t")

    def _transform(self):

        """

        Private method called upon initiation of NCTCTable.

        Assigns the HTML table descriptor ('summary') as name to the class (self.name)
        and extracts data from the HTML table (self.dataframe), including FTP links to
        assembled genomes (self.links).

        Calls:

        _transform_links()
        _clean_tables()

        """

        self.name = self.table["summary"]

        if self.name == "NCTC genomes":
            self.name = "genomes"

        data = [[td.text for td in row.findAll("td")]
                for row in self.table.findAll("tr")]

        ftp_links = [[a["href"] for a in row.find_all("a", href=True)
                      if a["href"].startswith("ftp")]
                     for row in self.table.findAll("tr")]

        links = [self._transform_links(links) for links in ftp_links]

        self.dataframe = pandas.DataFrame(data[1:], columns=[tag.text for tag in self.table.findAll("th")])

        self.links = pandas.DataFrame(links[1:], columns=links[0])

        self._clean_tables()

        assert len(self.dataframe) == len(self.links)

    def summarise(self, top=10):

        """

        Prints a short summary of the current NCTC3000.

        :param top: Number of species and their genome counts to show in summary.
        :return: message: Summary message as formatted string.

        """

        stamp("Here is a summary of the project:")

        message = textwrap.dedent("""

        NCTC3000, {date}
        -------------------------------------------

        Database: {db_name}

        Species:                       {n_species}
        Entries:                       {n_entries}
        Entries with Assembly:         {n_assemblies}
        Entries with Chromosomes:      {n_chromosomes}

        Top {top_number}:

        {top_species_counts}

        -------------------------------------------

        """).format(date=time.strftime("%d-%m-%Y %H:%M"),
                    db_name=self.name,
                    n_species=len(self.dataframe["Species"].unique()),
                    n_entries=len(self.dataframe),
                    n_assemblies=self.dataframe["Manual Assembly"].count(),
                    n_chromosomes=len(self.dataframe[self.dataframe["Chromosomes"]
                                      .astype(float) > 0]),
                    top_number=top,
                    top_species_counts=self.dataframe["Species"].value_counts()[:top].to_string())

        print(message)

        return message

    def _clean_tables(self):

        """

        Need better cleaning of Tables...

        """

        self.dataframe = self.dataframe.replace("Pending", nan)

        self.dataframe.columns = ["Species", "Strain", "Sample", "Runs", "Automated Assembly", "Manual Assembly",
                                  "Chromosomes", "Plasmids", "Unidentified"]

        self.dataframe[["Chromosomes", "Plasmids", "Unidentified"]] = \
            self.dataframe[["Chromosomes", "Plasmids", "Unidentified"]].astype(float)

    @staticmethod
    def _transform_links(links):

        """

        Private static method to transform link data from HTML table into a Pandas dataframe.

        :param links: List of links extracted from HTML table on NCTC3000.
        :return: Dictionary of links from HTML table, missing links are None.

        """

        link_dict = {"manual_gff": nan,
                     "manual_embl": nan,
                     "automatic_gff": nan,
                     "automatic_embl": nan}

        if not links:
            return link_dict
        else:
            for link in links:
                if "manual" and "gff" in link:
                    link_dict["manual_gff"] = link
                if "manual" and "embl" in link:
                    link_dict["manual_embl"] = link
                if "automatic" and "gff" in link:
                    link_dict["automatic_gff"] = link
                if "automatic" and "embl" in link:
                    link_dict["automatic_embl"] = link

        return link_dict


class NCTCProgram:

    """

    Superclass NCTCProgram - minimal wrappers for Abricate and MLST

            htpps://github.com/tseemann/mlst
            htpps://github.com/tseemann/abricate

    """

    def __init__(self, target_path=os.getcwd(), out_dir="mlst", exec_path="mlst", force=False):

        self.out_dir = out_dir
        self.exec_path = exec_path
        self.target_path = target_path

        self.output = None

        if not self.target_path.endswith("/"):
            self.target_path += "/"

        if not self.out_dir.endswith("/"):
            self.out_dir += "/"

        self.go = True

        if os.path.exists(out_dir):
            stamp("Analysis exists at", self.out_dir)
            if force:
                stamp("Force activated, removing directory", self.out_dir)
                shutil.rmtree(self.out_dir)
            else:
                self.go = False

        if self.go:
            stamp("Creating analysis directory", self.out_dir)
            os.makedirs(out_dir)


class MLST(NCTCProgram):

    def run(self, min_cov=80, min_id=90):

        if self.go:
            output = os.path.join(self.out_dir, "mlst.report")

            cmd = " ".join([self.exec_path, "--mincov", str(min_cov), "--minid", str(min_id),
                            self.target_path + "*", ">", output])

            try:
                stamp("Running mlst with database", "on", self.target_path)
                call(cmd, shell=True, stdout=open(os.devnull, 'w'), stderr=open(os.devnull, 'w'))
            except CalledProcessError:
                shutil.rmtree(self.out_dir)
                stamp("Could not run MLST. Removing directory", self.out_dir)
            except KeyboardInterrupt:
                stamp("Interrupted by user, removing analysis directory", self.out_dir)
                shutil.rmtree(self.out_dir, ignore_errors=True)

            return output


class Abricate(NCTCProgram):

    def _get_file_paths(self, file, db):

        file_name, extension = os.path.splitext(file)
        file_path = os.path.join(self.target_path, file)
        output = os.path.join(self.out_dir, file_name + "_" + db + ".tab")

        return file_path, output

    def run(self, db="resfinder", min_id=80):

        if self.go:
            files = [file for file in os.listdir(self.target_path) if file.endswith(FASTA_EXTENSIONS)]

            fail = []
            for file in files:

                file_path, output = self._get_file_paths(file, db)

                cmd = " ".join([self.exec_path, "--db", db, "--minid", str(min_id), "--nopath", file_path, ">", output])

                try:
                    stamp("Running Abricate with database", db, "on", file)
                    call(cmd, shell=True, stdout=open(os.devnull, 'w'), stderr=open(os.devnull, 'w'))
                except CalledProcessError:
                    stamp("Could not run Abricate for file", file)
                    fail.append(file)
                    continue
                except KeyboardInterrupt:
                    stamp("Interrupted by user, removing analysis directory", self.out_dir)
                    shutil.rmtree(self.out_dir, ignore_errors=True)
                finally:
                    if len(fail) == len(files):
                        stamp("Failed to run all files in Abricate, removing directory", self.out_dir)
                        shutil.rmtree(self.out_dir, ignore_errors=True)

            summary_file = self.summarize(db=db)

            return summary_file

    def summarize(self, db):

        if self.go:
            summary = os.path.join(self.out_dir, db + ".report")

            cmd = " ".join([self.exec_path, "--summary", self.out_dir + "*.tab", ">", summary])

            try:
                stamp("Writing results from database", db, "to", summary)
                call(cmd, shell=True, stdout=open(os.devnull, 'w'), stderr=open(os.devnull, 'w'))
            except CalledProcessError:
                raise

            return summary


def stamp(*args):

    print(str(time.strftime("[%H:%M:%S]")) + " " + " ".join([str(arg) for arg in args]))


def get_copy_numbers(frames, cnv_mincov):
    cnv_frames = {}

    for output, dataframe in frames.items():
        if "virulence" in output or "resistance" in output:
            cnv_frame = {}
            for column in dataframe.columns:
                if column not in ("ID", "Count"):
                    cnv_frame[column] = _get_cnv_column(output, dataframe, column, cnv_mincov)
                else:
                    cnv_frame[column] = dataframe[column]

            cnv_frames[output + "_cnv"] = pandas.DataFrame(cnv_frame)

    # Return dictionary without key "mlst", just copy number variation for Abricate
    return cnv_frames


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


def clean_reports(frames):

    clean = {}

    for output, frame in frames.items():
        ids = frame.ix[:, 0]
        cleaned = frame[ids.str.endswith(FASTA_EXTENSIONS)]
        cleaned_ids = cleaned.iloc[:, 0].map(lambda x: _extract_file_name(x))
        cleaned = _clean_output_reports(output, cleaned)
        cleaned.insert(0, "ID", cleaned_ids)
        cleaned.reindex()

        clean[output] = cleaned

    return clean


def _clean_output_reports(output, dataframe):

    if "virulence" in output or "resistance" in output:
        keep = ~dataframe.columns.str.contains("|".join(['Unnamed:', '#FILE']))
        cleaned = dataframe.iloc[:, keep]
        cleaned.rename(columns={'NUM_FOUND': 'Count'}, inplace=True)

    elif output == "mlst":
        cleaned = dataframe.iloc[:, 1:]
        cleaned.rename(columns={1: 'Species', 2: 'MLST'}, inplace=True)
        cleaned.rename(columns=_get_mlst_header(cleaned), inplace=True)
    else:
        cleaned = None

    return cleaned


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


def _extract_file_name(entry):
    file_name = os.path.basename(entry)
    name, extension = os.path.splitext(file_name)

    return name


def merge_frames(frames, left=None, merge="outer"):

    merged = {}

    for output in frames.keys():
        if "virulence" in output or "resistance" in output:
            if left is None:
                merged_df = pandas.merge(left=frames["mlst"], right=frames[output], left_on="ID", right_on="ID",
                                         how=merge)
            else:
                merged_df = pandas.merge(left=left, right=frames[output], left_on="ID", right_on="ID",
                                         how=merge)

            merged[output] = merged_df

    return merged


def parse_reports(mlst_path, res_path, vir_path, resistance_db="resfinder", virulence_db="vfdb"):

    files = {"mlst": os.path.join(mlst_path, "mlst.report"),
             "resistance_" + resistance_db: os.path.join(res_path, resistance_db + ".report"),
             "virulence_" + virulence_db: os.path.join(vir_path, virulence_db + ".report")}

    frames = {}

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


def download_assemblies(gff_path, links, force=False):
    exist = 0
    downloaded = 0
    out_paths = []

    for entry in links.itertuples():
        file, url = entry.File, entry.Link

        out_path = os.path.join(gff_path, file)

        if os.path.exists(out_path) and not force:
            exist += 1
        else:
            stamp("Downloading", file, "from", url)
            try:
                urllib.request.urlretrieve(url, out_path)
            except KeyboardInterrupt:
                stamp("Interrupted by user, deleting last file at", out_path)
                if ospath.exist(out_path):
                    os.remove(out_path)
            downloaded += 1
            
        out_paths.append(out_path)

    stamp("Found", exist, "files and downloaded", downloaded, "new assemblies to", gff_path)

    return out_paths


def gff_to_fasta(gff_files, fasta_path, genbank_path, verbose=True):
    files = {}

    for file in gff_files:

        if verbose:
            stamp("Converting file", file)

        file_name, ext = os.path.splitext(file)
        base_name = os.path.basename(file_name)

        fasta = os.path.join(fasta_path, base_name + ".fasta")
        genbank = os.path.join(genbank_path, base_name + ".gbk")

        if not os.path.exists(fasta) or not os.path.exists(genbank):

            records = _parse_gff(file, file_name)

            if not os.path.exists(fasta):
                SeqIO.write(records, fasta, format="fasta")
            if not not os.path.exists(genbank):
                SeqIO.write(records, genbank, format="genbank")

        files[base_name] = {"gff": os.path.realpath(file),
                            "fasta": os.path.realpath(fasta),
                            "genbank": os.path.realpath(genbank)}

    return files


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


def get_links(bacteria, species_indices, assembly):
    links = bacteria.links.ix[species_indices, assembly + "_gff"].dropna()

    file_names = get_file_names(bacteria.dataframe.ix[links.index])

    files = pandas.concat([file_names, links], axis=1)
    files.columns = ["File", "Link"]

    return files


def get_file_names(link_entries):
    species_names = link_entries["Species"].str.replace(" ", "_")
    strain_names = species_names.str.cat([link_entries["Strain"]], sep="_")

    return strain_names.str.cat([".gff" for _ in range(len(link_entries))])
