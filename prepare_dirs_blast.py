import os
import json
######## CHECK DIRS
print("> Setting up for BLAST")
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
print("> Obtaining query gene info...")
######## READ QUERIES AND RECORD IT
queryfolder = pathdictionary["USERFILES"]
querysequencepaths = {}
for file in os.listdir(queryfolder):
    if file.endswith('.fna'):  
        path = os.path.join(queryfolder, file)
        queryname = os.path.splitext(file)[0]
        querysequencepaths[queryname] = path
pathforappdata = pathdictionary["App-Files"]
querydatafile = "query_sequences.json"
pathforquerydata = os.path.join(pathforappdata, querydatafile)
with open(pathforquerydata, 'w') as file:
    json.dump(querysequencepaths, file)
print("> User query genes obtained")
##### SPLIT QUERY SEQUENCE PATHS INTO DIFFERENT JSONS
queryloadONE = {}
queryloadTWO = {}
tempcountlist = []
halfqueries = []
print("> Prioritising and sharing workload")
with open(pathforquerydata) as querydictfull:
    allqueries = json.load(querydictfull)
    for query in allqueries.keys():
        tempcountlist.append(query)
halfrun = int(len(tempcountlist) / 2) 
for x in range(0, halfrun):
    halfqueries.append(tempcountlist[x])
for query in halfqueries:
    for queryname, path in allqueries.items():
        if query == queryname:
            path = allqueries[queryname]
            queryloadONE[queryname] = path
halfqueriesset = set(halfqueries)
allqueriesset = set(tempcountlist)
otherhalfqueries = allqueriesset - halfqueriesset
for query in otherhalfqueries:
    for queryname, path in allqueries.items():
        if query == queryname:
            addpath = allqueries[query]
            queryloadTWO[query] = addpath
print("> Finishing up prioritisation...")
jsonname1, jsonname2 = "querygenes1.json", "querygenes2.json"
printpath1 = os.path.join(pathdictionary["App-Files"], jsonname1)
printpath2 = os.path.join(pathdictionary["App-Files"], jsonname2)
with open(printpath1, "w") as ooga:
    json.dump(queryloadONE, ooga)
with open(printpath2, "w") as ooga:
    json.dump(queryloadTWO, ooga)
####### SOURCE TSV OF ASSEMBLY BARCODES AND BACTERIAL NAMES
tsvpathlist = [file for file in os.listdir(pathdictionary["USERFILES"]) if file.endswith(".tsv") and "Name" in file]
NameIDDict = {}
with open(os.path.join(pathdictionary["USERFILES"] ,tsvpathlist[0]), 'r') as input:
    for liner in input:
        barcode, name = liner.strip().split("\t")
        NameIDDict[barcode] = name
####### CHECK THAT THE TSV IS CORRECT
meow = {}
for barcoden in NameIDDict.keys():
    meow[barcoden] = "y"
    for item in os.listdir(pathdictionary["USERFILES"]):
        if item.endswith(".fasta") and barcoden in item:
            meow[barcoden] = "x"
####### REPORT FINDINGS
outlier_list = []
for key, value in meow.items():
    if value == "y":
        outlier_list.append(f"{NameIDDict[key]} sample, with assembly barcode: {key}")
if outlier_list:
    print(f"A list has been generated for name and assembly barcode errors: \n {outlier_list}")
    just_names = [bacname for bacname in NameIDDict.values()]
    duplicates_count_dict = {}
    report = []
    for info in just_names:
        if info not in duplicates_count_dict.keys():
          duplicates_count_dict[info] = 1
        else:
          duplicates_count_dict[info] += 1
    for bacterium, count in duplicates_count_dict.items():
        if count > 1:
            report.append(bacterium)
    if report:
        print(f"The following bacteria have duplicate names, please ammend: {report}")
        print(144)
    else:
        print("No tsv duplicates!")
else:
    print("No tsv errors!") ####### CHECK TSV FOR DUPLICATES
    just_names = [bacname for bacname in NameIDDict.values()]
    duplicates_count_dict = {}
    report = []
    for info in just_names:
        if info not in duplicates_count_dict.keys():
          duplicates_count_dict[info] = 1
        else:
          duplicates_count_dict[info] += 1
    for bacterium, count in duplicates_count_dict.items():
        if count > 1:
            report.append(bacterium)
    if report:
        print(f"The following bacteria have duplicate names, please ammend: {report}")
        print(144)
    else:
        print("No tsv duplicates!")
    