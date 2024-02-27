import os
import subprocess
import shutil
######## CHECK DIRECTORIES
desktopdirectory = os.path.expanduser("~/Desktop")
pathforworkingfolder = os.path.join(desktopdirectory, "PROBE")
if not os.path.exists(pathforworkingfolder):
    os.mkdir(pathforworkingfolder)
subfolders = ["USERFILES", "RESULTS", "BLAST-Databases",
              "BLAST-Search-Results", "App-Files", "QUAST", "Overlap-Genomes"]
pathdictionary = {"PROBE": pathforworkingfolder}
for folder in subfolders:
    subfolderpath = os.path.join(pathforworkingfolder, folder)
    if not os.path.exists(subfolderpath):
        os.makedirs(subfolderpath)
    pathdictionary[folder] = subfolderpath
####### CHECK IF GENOMES HAVE BEEN MOVED
are_genomes_already_there = any(file.endswith(".fasta") for file in os.listdir(pathdictionary["App-Files"]))
if not are_genomes_already_there:
    for filename in os.listdir(pathdictionary['USERFILES']):
        if filename.endswith(".fasta"):
            src_file = os.path.join(pathdictionary['USERFILES'], filename)
            dest_file = os.path.join(pathdictionary["App-Files"], filename)
            if os.path.isfile(src_file):
                shutil.copy(src_file, dest_file)
####### GET AND READ TSV OF BARCODES AND NAMES, RENAME GENOMES
tsv_pathlist = [file for file in os.listdir(pathdictionary["USERFILES"]) if file.endswith(".tsv") and "Name" in file]
NameIDDict = {}
with open(os.path.join(pathdictionary["USERFILES"] ,tsv_pathlist[0]), 'r') as input:
    for liner in input:
        barcode, name = liner.strip().split("\t")
        NameIDDict[barcode] = name
for barcode, name in NameIDDict.items():
    currentfilename = f'{barcode}.fasta'
    newfilename = f'{name}.fasta'
    current_filepath = os.path.join(pathdictionary['App-Files'], currentfilename)
    new_filepath = os.path.join(pathdictionary['App-Files'], newfilename)
    if os.path.exists(current_filepath):
        subprocess.run(['mv', current_filepath, new_filepath])
######### RUN COMPLETE QUAST ANALYSIS
listofGAdirsinautoblast = [x for x in os.listdir(pathdictionary["App-Files"]) if not x.endswith(".json")]
fullpathlist = [os.path.join(pathdictionary["App-Files"], GA) for GA in listofGAdirsinautoblast if GA !='.DS_Store']
subprocess.run(['quast.py', '-t 1', *fullpathlist, '-o', pathdictionary['QUAST']])