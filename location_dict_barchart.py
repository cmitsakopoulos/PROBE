import os
import json
import pandas as pd
import seaborn as sns 
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
################## CHECK DIRECTORIES ################
desktopdirectory = os.path.expanduser("~/Desktop")
pathforworkingfolder = os.path.join(desktopdirectory, "MINAGIASGENE")
if not os.path.exists(pathforworkingfolder):
    os.mkdir(pathforworkingfolder)
subfolders = ["USERFILES", "RESULTS", "BLAST-Databases", "BLAST-Search-Results", "App-Files", "QUAST", "Overlap-Genomes"]
pathdictionary = {"MINAGIAS": pathforworkingfolder}
for folder in subfolders:
    subfolderpath = os.path.join(pathforworkingfolder, folder)
    if not os.path.exists(subfolderpath):
        os.makedirs(subfolderpath)
    pathdictionary[folder] = subfolderpath
######### START LOCATION BASED INTERPRETATION
querydatapath = os.path.join(pathdictionary["App-Files"], f"query_sequences.json")
tsvpathlist = [file for file in os.listdir(pathdictionary["USERFILES"]) if file.endswith(".tsv") and "Path" in file]
with open(querydatapath, 'r') as info:
    whattododict = json.load(info)
frame = pd.read_csv(os.path.join(pathdictionary["USERFILES"], tsvpathlist[0]), sep='\t')
bacloc = {}
for row in frame.values:
    bacloc[row[0]] = row[1]  
######### PREPARE TEMP DICT 
pathovar_empty_dict = {}
for pat in bacloc.values():
    pathovar_empty_dict[pat] = ""
######## SORT BACS WITH LOCATION KEY
pathovar_count = {}
for pathovar in pathovar_empty_dict.keys():
    bacs = []
    for bac, local in bacloc.items():
        if local == pathovar:
            bacs.append(bac)
    pathovar_empty_dict[pathovar] = bacs
    pathovar_count[pathovar] = len(bacs)
######### COUNT PRESENCE WITH LOCATION
all_pathovars_dict = {}
for querygene in whattododict.keys():
    geolocdict = {}
    for file in os.listdir(pathdictionary["App-Files"]):
        if file.endswith(".json"):
            filepath = os.path.join(pathdictionary["App-Files"], file)
            if f"{querygene}_Presence" in filepath:
                with open(filepath, 'r') as presence:
                    temppresencedict = json.load(presence)
                    for path, bacterias in pathovar_empty_dict.items():
                        newlist = []
                        sumlist = []
                        for bacteria in bacterias:
                            newlist.append(temppresencedict[bacteria])
                        for item in newlist:
                            if "inconclusive" in item or "present" in item:
                                sumlist.append(1)
                            else:
                                sumlist.append(0)
                        summed = sum(sumlist)
                        geolocdict[path] = summed 
                    if sum(geolocdict.values()) > 0:
                        numberoflocs = len(geolocdict.keys())
                        sns.set(style="whitegrid")
                        plt.figure(figsize=(16, 14))
                        secondsavelocation = os.path.join(pathdictionary["RESULTS"], f"{querygene}distribution_pathovars_BLASTn.pdf")
                        plt.bar(geolocdict.keys(), geolocdict.values(), color='black', width=0.5)
                        plt.title(f"{querygene} presence distribution by X. euvesicatoria pathovar", fontstyle="italic", fontweight="bold", fontsize=22)
                        plt.xticks(rotation=70, fontsize=12)
                        plt.setp(plt.gca().get_xticklabels(), ha='right')
                        plt.ylabel('Number of samples with gene presence', fontstyle="italic", fontweight="bold", fontsize=20, labelpad=10)
                        plt.savefig(secondsavelocation, format='pdf')
                        plt.close()
                        all_pathovars_dict[querygene] = [x for x in geolocdict.values()]   
                    else:
                        newsave = os.path.join(pathdictionary["RESULTS"], "IMPORTANT_BLASTn.txt")
                        with open(newsave, 'a') as do:
                            do.write(f"{querygene} is absent from every genome provided. \n")
#### PRINT CONCATENATED TSV
save_help = os.path.join(pathdictionary["RESULTS"], f"ALL_pathovars_BLASTn.tsv")
all_pathovars_dict["Pathovars"] = [y for y in geolocdict.keys()]
all_dataframe = pd.DataFrame(list(all_pathovars_dict.items()), columns=['Key', 'Value'])
all_dataframe.to_csv(save_help, sep='\t', index=False)
print(all_dataframe)


                        



