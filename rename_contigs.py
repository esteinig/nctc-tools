import os
import shutil

files = [file for file in os.listdir(os.getcwd()) if file.endswith(".fa")]

target_path = "/home/esteinig/analysis/ST772/renamed_assemblies"

translate = {}
with open("rename_contigs.csv") as infile:
    for line in infile:
        content = line.strip().split(",")
        translate[content[0]] = content[1]

for file in files:
    try:
        new_file = os.path.join(target_path, translate[file] + ".fasta")
        shutil.copyfile(os.path.join(os.getcwd(), file), new_file)
    except KeyError:
        print("Could not find file", file)

