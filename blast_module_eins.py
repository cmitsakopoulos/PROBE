import os
import json
import subprocess
import datetime
import sys
######### CHECKING DIRECTORIES
desktopdirectory = os.path.expanduser("~/Desktop")
pathforworkingfolder = os.path.join(desktopdirectory, "PROBE")
if not os.path.exists(pathforworkingfolder):
    os.mkdir(pathforworkingfolder)
subfolders = ["USERFILES", "RESULTS", "BLAST-Databases", "BLAST-Search-Results", "App-Files", "QUAST", "Overlap-Genomes"]
pathdictionary = {"PROBE": pathforworkingfolder}
for folder in subfolders:
    subfolderpath = os.path.join(pathforworkingfolder, folder)
    if not os.path.exists(subfolderpath):
        os.makedirs(subfolderpath)
    pathdictionary[folder] = subfolderpath
######### SOURCING TSV
print("> Sourcing genome barcode and name .tsv file...")
tsvpathlist = [file for file in os.listdir(pathdictionary["USERFILES"]) if file.endswith(".tsv") and "Name" in file]
NameIDDict = {}
with open(os.path.join(pathdictionary["USERFILES"] ,tsvpathlist[0]), 'r') as input:
    for liner in input:
        barcode, name = liner.strip().split("\t")
        NameIDDict[barcode] = name
######## RUNNING QUAST AND PRINTING RECEIPT
querydatapath = os.path.join(pathdictionary["App-Files"], "querygenes1.json")
resultsoutput = pathdictionary["BLAST-Search-Results"]
databasedirectory = pathdictionary["BLAST-Databases"]
print("> Module one in operation...")
x=0
with open(querydatapath) as info:
    whattododict = json.load(info)
for queryname, querypath in whattododict.items():
    for barcode, name in NameIDDict.items():
        outputfile = os.path.join(resultsoutput, f"{name}_result_{queryname}.txt")
        if not os.path.exists(outputfile):
            x+=1
            print(f"> Generating {outputfile}")
            cmd = ['blastn', '-query', querypath, '-db', barcode, '-outfmt', '6', '-evalue', '1e-50', '-out', outputfile]
            subprocess.run(cmd, check=True, cwd= databasedirectory) 
        else:
            print(f"{outputfile} already exists...moving on")
if not x == 0:
    currentime = datetime.datetime.now() 
    receiptdict = {names: f"BLAST ran with this query at {currentime}" for names in whattododict.keys()}
    printpath = os.path.join(pathdictionary["RESULTS"], f"1st_BLASTnRECEIPT_{currentime}.json")
    with open(printpath, "w") as receipt:
        json.dump(receiptdict,receipt)  
    print("> Module one finished operations, receipt printed in RESULTS.")
else:
    print("> Module one discovered that you have already ran BLAST with the query genes and genomes, \n that you have deposited in USERFILES. In case you add more query genes, run this function again and \n the program will automatically run BLAST only with the additional genes to conserve time.")
sys.exit()