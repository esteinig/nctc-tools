"""

NCTC Utils:
    -- Command Line Parser
    -- NCTC Table for parsing HTML
    -- NCTC Programs for local Abricate / MLST / Roary

"""

import os
import time
import pandas
import shutil
import textwrap
import argparse

from numpy import nan

from subprocess import check_output
from subprocess import CalledProcessError


class CommandLine:

    def __init__(self):

        parser = argparse.ArgumentParser()

        parser.add_argument("--version", "-v", help="print version and exit")

        parser.add_argument("--project", "-p", type=str, required=False,
                            default="NCTC3000", dest="project", help="project path")

        parser.add_argument("--species", "-s", type=str, default="Staphylococcus aureus", required=False,
                            dest="species", help="species name")

        parser.add_argument("--force", "-f", default=False, dest="overwrite", action="store_true",
                            help="force overwrite of project")

        parser.add_argument("--verbose", "-d", default=False, dest="verbose", action="store_true",
                            help="print details")

        subparsers = parser.add_subparsers(help='Command-line interface for scanNCTC')

        make_parser = subparsers.add_parser("make")

        make_parser.add_argument("--force", "-f", default=False, dest="force", action="store_true",
                                 help="force download of assemblies")

        make_parser.add_argument("--chromosomes", "-c",  default=False, dest="chromosomes",
                                 action="store_true", help="only assemblies with identified chromosomes")

        make_parser.set_defaults(subparser='make')

        analysis_parser = subparsers.add_parser("analyse")

        analysis_parser.add_argument("--user_path", "-u", type=str, default=None, required=False,
                                    dest="output", help="user path with user_path/fasta")
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

        analysis_parser.set_defaults(subparser='analyse')

        collect_parser = subparsers.add_parser("collect")

        collect_parser.add_argument("--output", "-o", type=str, default="report.tab", required=True,
                                    dest="output", help="summary output file")
        collect_parser.add_argument("--user_path", "-u", type=str, default=None, required=False,
                                    dest="user_path", help="user path with user_path/fasta")
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

        update_parser = subparsers.add_parser("update")

        self.args = parser.parse_args()


class NCTCTable:

    """ Class parsing a HTML Soup table into a Pandas dataframe """

    def __init__(self, project, table=None, dataframe=pandas.DataFrame(), links=pandas.DataFrame(), name=None):

        """

        Initialize class NCTCTable.

        Calls:

        _transform()

        :param table: HTML soup string containing table section of NCTC3000.

        """

        self.project = project
        self.table = table
        self.name = name

        self.dataframe = dataframe
        self.links = links

        if self.table is not None:
            self._transform()

    def store(self, path):

        filename = os.path.join(path, self.project + "_" + self.name + "_" + time.strftime("%d-%m-%Y-%H-%M"))

        filename_table = os.path.join(filename + "_data.tab")
        filename_links = os.path.join(filename + "_links.tab")

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


def stamp(*args):

    print(str(time.strftime("[%H:%M:%S]")) + "\t" + " ".join([str(arg) for arg in args]))


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

        if force:
            shutil.rmtree(self.out_dir)

        os.makedirs(out_dir, exist_ok=True)


class MLST(NCTCProgram):

    def run(self, min_cov=80, min_id=90):

        output = os.path.join(self.out_dir, "mlst.report")

        cmd = " ".join([self.exec_path, "--mincov", str(min_cov), "--minid", str(min_id),
                        self.target_path + "*", ">", output])

        try:
            check_output(cmd, shell=True)
        except CalledProcessError:
            shutil.rmtree(self.out_dir)
            stamp("Could not run MLST. Removing directory", self.out_dir)
        except KeyboardInterrupt:
            stamp("Interrupted by user, removing analysis directory", self.out_dir)
            shutil.rmtree(self.out_dir)

        return output


class Abricate(NCTCProgram):

    def _get_file_paths(self, file):

        file_name, extension = os.path.splitext(file)
        file_path = os.path.join(self.target_path, file)
        output = os.path.join(self.out_dir, file_name + ".tab")

        return file_path, output

    def run(self, db="resfinder", min_id=80, ext=".fasta"):

        files = [file for file in os.listdir(self.target_path) if file.endswith(ext)]

        fail = []
        for file in files:

            file_path, output = self._get_file_paths(file)

            cmd = " ".join([self.exec_path, "--db", db, "--minid", str(min_id), "--nopath", file_path, ">", output])

            try:
                check_output(cmd, shell=True)
            except CalledProcessError:
                stamp("Could not run Abricate for file", file)
                fail.append(file)
                continue
            except KeyboardInterrupt:
                stamp("Interrupted by user, removing analysis directory", self.out_dir)
                shutil.rmtree(self.out_dir)
            finally:
                if len(fail) == len(files):
                    stamp("Failed to run all files in Abricate, removing directory", self.out_dir)
                    shutil.rmtree(self.out_dir)

        summary_file = self.summarize(db=db)

        return summary_file

    def summarize(self, db):

        summary = os.path.join(self.out_dir, db + ".report")

        cmd = " ".join([self.exec_path, "--summary", self.out_dir + "*.tab", ">", summary])

        try:
            check_output(cmd, shell=True)
        except CalledProcessError:
            raise

        return summary



