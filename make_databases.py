import subprocess
import os
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
############# MAKE DATABASES
print("> Initiating BLAST Database generation from user provided genomes...")
genomefiles = [x for x in os.listdir(pathdictionary["USERFILES"]) if x.endswith('.fasta')]
for file in genomefiles: 
    inputpath = os.path.join(pathdictionary["USERFILES"], file)
    filename = os.path.splitext(file)[0]
    outputpath = os.path.join(pathdictionary["BLAST-Databases"], filename)
    subprocess.run(['makeblastdb', '-in', inputpath, '-dbtype', 'nucl', '-out', outputpath], check=True)