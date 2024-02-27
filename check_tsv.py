import os
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
###### GATHER AND SAVE NAME BARCODE TSV
tsv_pathlist = [file for file in os.listdir(pathdictionary["USERFILES"]) if file.endswith(".tsv") and "Name" in file]
NameIDDict = {}
with open(os.path.join(pathdictionary["USERFILES"] ,tsv_pathlist[0]), 'r') as input:
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
else:
    print("No tsv errors!") ####### CHECK TSV FOR DUPLICATES
    just_names = [bacname for bacname in NameIDDict.values()]
    duplicates_count_dict = {}
    report = []
    for info in just_names:
        duplicates_count_dict[info] += 1
    for bacterium, count in duplicates_count_dict.items():
        if count > 1:
            report.append(bacterium)
    if report:
        print(f"The following bacteria have duplicate names, please ammend: {report}")
    else:
        print("No tsv duplicates!")
