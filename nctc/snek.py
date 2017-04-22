"""

Convenience wrapper for Snakemake on cluster environments.

Notes:

    Configuration file for Snakemake needsto have one sub-dictionary accessible by

        "snek"      dictionary with command-line configuration parameters for Snakemake

    Snek configuration is extracted / copied to the working directory.

    Files used as initial input for Snakemake (in Snakefile) must be present in the working directory, e.g.
    for NCTC there needs to be:

        fasta/{sample}.fasta

Tested on Cheetah, CDU.

"""

import os
import json
import shutil

from subprocess import check_output, CalledProcessError


class Snek:

    def __init__(self, path=os.getcwd(), snake_file="nctc.snk", config_file="nctc.json",
                 executable="snakemake", major=3, minor=11, force=False):

        self.current = os.getcwd()

        self.path = path
        self.force = force

        self.snake_file = snake_file
        self.config_file = config_file

        self.executable = executable

        self.major = major
        self.minor = minor

        self.config = {
            "force": self.force,
            "jobs": 300,
            "latency": 300,
            "conda": True,
            "cluster": "qsub -l pmem=1gb -l walltime=2:00:00 -l nodes=1:ppn=1 -V -S /bin/bash"
        }

        self._check_version()
        self._configure_snek()

    def _configure_snek(self):

        self.snake_path = os.path.join(self.path, os.path.basename(self.snake_file))
        self.config_path = os.path.join(self.path, os.path.basename(self.config_file))

        if os.path.exists(self.snake_path):
            os.remove(self.snake_path)
        if os.path.exists(self.config_path):
            os.remove(self.config_path)

        shutil.copy(self.snake_file, self.snake_path)

        with open(self.config_file) as infile:
            config = json.load(infile)

        # This is particular format of a config file for Snakemake (.json) containing
        # config params used in Snakefile (dictionary), with one entry with key
        # "snek", which contains as values the higher parameters for executing snakemake through command-line.

        # The entry is removed from dictionary and configures the class when initialised, while the dictionary for
        # Snakemake (config file given basename) is written to the config path, which is the targe tpath this class
        # is initialised with. This allows execution of the Snakefile in the location it references in relative paths
        # in the rules of the Snakefile (e.g. glob_wildcards("fasta/{sample}.fasta").

        self.config = config.pop("snek")

        with open(self.config_path, "w") as outfile:
            json.dump(config, outfile)

    def _check_version(self):

        try:
            version = check_output([self.executable, "--version"])

            components = [int(component) for component in str(version, "utf-8").split(".")]

            if components[0] < self.major:
                raise ValueError("You need Snakemake major version >= 3")
            elif components[1] < self.minor:
                raise ValueError("You need Snakemake minor version >= 11")
            else:
                print("Version of Snakemake is", str(version, "utf-8"))

        except CalledProcessError:
            print("Could not run Snakemake.")

    def _make_command(self):

        cmd = [self.executable, "--cluster", self.config["cluster"], "--jobs", str(self.config["jobs"]),
               "--latency-wait", str(self.config["latency"]), "-s", self.snake_path]

        if self.config["conda"]:
            cmd += ["--use-conda"]

        if self.config["force"]:
            cmd += ["--forceall"]

        return cmd

    def _make_run(self):

        log_dir = os.path.join(self.path, "logs")

        if not os.path.exists(log_dir):
            os.makedirs(log_dir)
        else:
            print("Replacing log directory...")

            shutil.rmtree(log_dir)
            os.makedirs(log_dir)

        return log_dir

    def _clean_successful_run(self, run_log, log_dir):

        print("Snek was run successfully.")

        submission_logs = [file for file in os.listdir(self.path)
                           if file.startswith("snakejob")]

        for log in submission_logs:
            shutil.move(os.path.join(self.path, log), os.path.join(log_dir, log))

        with open(os.path.join(log_dir, "snakemake.log"), "w") as out_file:
            out_file.write(run_log)

    def _clean_failed_run(self):

        """ Needs error handling on failed runs, identify samples, re-run."""

        print("Error in running Snakemake, please check submission logs.")

    def run(self):

        os.chdir(self.path)

        log_dir = self._make_run()

        cmd = self._make_command()

        try:
            run_log = check_output(cmd)
            run_log = str(run_log, "utf-8")
            self._clean_successful_run(run_log, log_dir)
        except CalledProcessError:
            self._clean_failed_run()
        finally:
            os.chdir(self.current)
