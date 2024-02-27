import os
import json
import pandas as pd
import seaborn as sns 
import matplotlib.pyplot as plt
from statistics import stdev
from matplotlib.backends.backend_pdf import PdfPages
################## CHECK DIRECTORIES ################
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
#################GET ID NAMES#######################
tsvpathlist = [file for file in os.listdir(pathdictionary["USERFILES"]) if file.endswith(".tsv") and "Name" in file]
NameIDDict = {}
with open(os.path.join(pathdictionary["USERFILES"] ,tsvpathlist[0]), 'r') as input:
    for liner in input:
        barcode, name = liner.strip().split("\t")
        NameIDDict[barcode] = name
####### READ QUERIES #### INTERPRET RESULTS ####### PRINT
querydatapath = os.path.join(pathdictionary["App-Files"], f"query_sequences.json")
resultsoutput = pathdictionary["BLAST-Search-Results"]
databasedirectory = pathdictionary["BLAST-Databases"]
dataforheatmap = {}
with open(querydatapath, 'r') as info:
    whattododict = json.load(info)
for query, querypath in whattododict.items():
    Presencedict = {}
    Percentagedict = {}
    listofpresenceforspecficgene = []
    multiplehitsdict = {names: f"{query} less than three hits" for names in NameIDDict.values()}
    sequence_length = 0
    with open(querypath, 'r') as file:
        data = file.readlines()
        sequence = "".join(line.strip() for line in data if not line.startswith(">")).replace("\n", "")
        sequence_length = len(sequence)
    queryresults = pathdictionary["BLAST-Search-Results"]
    for bacname in NameIDDict.values():
        if bacname not in Presencedict.keys():
            for file in os.listdir(queryresults):
                filepath = os.path.join(queryresults, file)
                if bacname in filepath:
                    if query in filepath:
                        filesize = os.path.getsize(filepath)
                        if filesize < 1:
                            Presencedict[bacname] = f"{query} absent"
                            Percentagedict[bacname] = 0
                            continue
                        else:
                            with open(filepath, 'r') as filedata:
                                lineslist = filedata.readlines()
                                topresult = lineslist[0]
                                linesplit = topresult.split('\t')
                                percidentity, percoverlap = float(linesplit[2]), (float(linesplit[3]) / sequence_length) * 100
                                if percidentity > 97 and percoverlap > 97:
                                    Presencedict[bacname] = f"{query} present"
                                    Percentagedict[bacname] = round(percidentity, 1)
                                    continue
                                else:                
                                    Presencedict[bacname] = f"{query} inconclusive"
                                    Percentagedict[bacname] = round(percidentity, 1)
                                    continue
                    else:
                        continue
                else:
                    continue
        else:
            continue
        listofpresenceforspecficgene.append(Percentagedict[bacname])
    dataforheatmap[query] = listofpresenceforspecficgene
    presencefile = f"{query}_Presence.json"
    outputpath1 = os.path.join(pathdictionary["App-Files"], presencefile)
    with open(outputpath1, 'w') as file:
        json.dump(Presencedict, file)  
###### RE ARRANGE ANNOTATION BASED ON HIGHEST PRESENCE
re_arranged_data_for_heatmap = {}
temp_dict_heatmap = {}
gene_stats_dict = {}
for key, listing in dataforheatmap.items():
    temp_sum = sum(listing)
    temp_dict_heatmap[key] = temp_sum
just_sums = [x for x in temp_dict_heatmap.values()]
just_sums.sort(reverse=True)
for index in just_sums:
    for gene in dataforheatmap.keys():
        if temp_dict_heatmap[gene] == index:
            re_arranged_data_for_heatmap[gene] = dataforheatmap[gene]
            gene_stats_dict[gene] = dataforheatmap[gene]
stylista = []
for genename, summed_value in temp_dict_heatmap.items():
    if summed_value == 0:
        stylista.append(genename)
        del re_arranged_data_for_heatmap[genename]
        del gene_stats_dict[genename]
notice_path = os.path.join(pathdictionary["RESULTS"], "BLASTn_HEATMAP_WARNING.txt")
application_component_path = os.path.join(pathdictionary["App-Files"], "ABSENTGENES.json")
with open (notice_path, 'w') as write:
    write.write(f"The following genes are not present (0 hits, 0% overlap) in any of the genomes provided. \nThey have therefore been omitted from the figure: \n {stylista}")
application_component = {}
for item in stylista:
    application_component[item] = 0
with open(application_component_path, 'w') as record:
    json.dump(application_component, record)
###### ANOTTATE AND SAVE HEATMAP
genomes = [x for x in NameIDDict.values()]
height = len(genomes) - 2
length = len([x for x in whattododict.keys()]) + 2
savelocation = os.path.join(pathdictionary["RESULTS"], "GENE_PRESENCE_HEATMAP_BLASTn.pdf")
re_arranged_data_for_heatmap["Genome"] = genomes
df = pd.DataFrame(re_arranged_data_for_heatmap)
df = df.set_index('Genome')
plt.figure(figsize=(length, height))
sns.heatmap(df, annot=True, fmt=".2f", annot_kws={"size": 8}, cmap='inferno')
plt.xlabel("Genes", fontstyle='italic', fontweight='bold', fontsize=18, labelpad=20)
plt.ylabel("Genomes", fontstyle='italic', fontweight='bold', fontsize=18, labelpad=20)
plt.savefig(savelocation, format='pdf', bbox_inches='tight')
plt.close()
##### MAKE NEW DICT WITH STATS
data_frame_dict = {}
for x, values in gene_stats_dict.items():
    summed = 0
    floats = []
    for value in values:
        summed += float(value)
        floats.append(float(value))
    temp_list = []
    average = summed / len(values)
    average_appropriate = round(average, 1)
    temp_list.append(average_appropriate)
    st_dev = stdev(floats)
    st_dev_fixed = round(st_dev, 1)
    temp_list.append(st_dev_fixed)
    data_frame_dict[x] = temp_list
###### GENERATE TABLE WITH TRANSPOSED STATS FOR GENES
df_for_stats = pd.DataFrame.from_dict(data_frame_dict, orient='index', columns=['Mean PSS (AVG_%)', 'PSS Standard Deviation'])
target = os.path.join(pathdictionary["RESULTS"], "GENE_PRESENCE_STATS_BLASTn.pdf")
with PdfPages(target) as pdf:
    fig, ax = plt.subplots(figsize=(5, 9))
    ax.axis('tight')
    ax.axis('off')
    table_object = ax.table(cellText=df_for_stats.values, colLabels=df_for_stats.columns, rowLabels=df_for_stats.index, cellLoc='center', loc='center')
    for (column, row), cell in table_object.get_celld().items():
        if column == 0:  
            cell.set_text_props(fontweight='bold', fontsize="13")
        if row == -1:
            cell.set_text_props(fontweight='bold', fontsize="13")
    plt.title('Percentage Sequence Similarity (PSS) Between Xeu Genomes and T3SS Genes using BLASTn', fontsize=13, pad=5)
    pdf.savefig(fig, bbox_inches='tight')
    plt.close()
save_help = os.path.join(pathdictionary["RESULTS"], "BLASTn_Presence_Tabulated.tsv")
with open(save_help, 'w') as new:
    for key_name, contents in data_frame_dict.items():
        new.write(f"{key_name}\t{contents[0]}\t{contents[1]}\n")
