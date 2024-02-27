import os 
import json
import subprocess
############### CHECK DIRECTORIES
desktopdirectory = os.path.expanduser("~/Desktop")
pathforworkingfolder = os.path.join(desktopdirectory, "PROBE")
if not os.path.exists(pathforworkingfolder):
    os.mkdir(pathforworkingfolder)
subfolders = ["USERFILES", "RESULTS", "BLAST-Databases", "BLAST-Search-Results", "App-Files", "QUAST", "Overlap-Genomes", "ITERATIVE"]
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
###### SOURCE FULL QUERY SEQUENCE PATHS
queryfolder = pathdictionary["USERFILES"]
with open(os.path.join(pathdictionary["App-Files"], "querygenes1.json")) as meow:
   querysequencepaths = json.load(meow)
###### SOURCE GENES WHICH ARE PRESENT
pre_determined_absence = os.path.join(pathdictionary["App-Files"], "ABSENTGENES.json")
with open(pre_determined_absence) as just_genes:
   absence_dict = json.load(just_genes)
for x in absence_dict.keys():
    if x in querysequencepaths.keys():
        del querysequencepaths[x]
###### IDENTIFY OVERLAP STARTS AND ENDS, CONTIG NAMES -> SAVE LOCALLY FOR NEXT FEATURE, IN CASE OF COMPUTER FAILURE
blastresults = pathdictionary["BLAST-Search-Results"]
tempcalcdict = {}
begin_new_metric = {}
for query in querysequencepaths.keys():
    list_overlaps = []
    for genomename in NameIDDict.values():
        for file in os.listdir(blastresults):
            filepath = os.path.join(blastresults, file)
            if file.endswith(".txt"):
                if genomename in filepath:
                    if query in filepath:
                        if os.path.getsize(filepath) > 1:
                            with open(filepath) as filedata:
                                lineslist = filedata.readlines()
                                topresult = lineslist[0]
                                linesplit = topresult.split('\t')
                                contig, overlapstart, overlapend = linesplit[1], int(linesplit[8]), int(linesplit[9])
                                templist = []
                                templist.append(contig)
                                templist.append(overlapstart)
                                templist.append(overlapend)
                                list_overlaps.append(overlapstart)
                                list_overlaps.append(overlapend)
                                tempcalcdict[genomename] = templist
                        else:
                            tempcalcdict[genomename] = ""
    newfilename = f"{query}_overlaps.json"
    newfilepath = os.path.join(pathdictionary["App-Files"], newfilename)
    with open(newfilepath, "w") as newfile:
        json.dump(tempcalcdict, newfile)  
########## 
########## MODULATED OVERLAP RIP COMPONENT
for querygene in querysequencepaths.keys():
    notice_text = os.path.join(pathdictionary["ITERATIVE"], f"REMOVED_GENOMES_{querygene}_TREE.txt")
    concat_file = os.path.join(pathdictionary["ITERATIVE"], f"{querygene}_CONCAT.fasta")
    if not os.path.exists(concat_file):
        for appfile in os.listdir(pathdictionary["App-Files"]):
            furtherfilepath = os.path.join(pathdictionary["App-Files"], appfile)
            if f"{querygene}_overlaps" in furtherfilepath:
                overlap_sequence = ""
                with open(furtherfilepath) as data:
                    overlapdict = json.load(data)  
        for genomefile in os.listdir(pathdictionary["USERFILES"]):
            if genomefile.endswith(".fasta"):
                genetic_code = ""
                fullinfo = ""
                overlap_sequence = ""
                reversedoverlapsequence = ""
                translatefirst = ""
                orthodoxtranslation = ""
                inb_tween = genomefile.strip(".fasta")
                genome = NameIDDict[inb_tween]
                fixedname = genome.replace(" ", "_")
                print(f"> Working with {querygene}...searching for overlap sequence on {genome}")
                filepath_genome = os.path.join(pathdictionary["USERFILES"], genomefile)                     
                if not overlapdict[genome]:
                    print(f"> {querygene} did not have a BLAST hit with: {genome}. Moving on...")
                    with open(notice_text, 'a') as save:
                        save.write(f"{genome} was not added to the phylogenetic tree as it did not have a {querygene} homologue!")
                else:
                    indexes = overlapdict[genome]
                    contig_name, start, end = indexes[0], indexes[1], indexes[2]
                    list_of_overlap_points = []
                    list_of_overlap_points.append(indexes[1])
                    list_of_overlap_points.append(indexes[2])
                    with open(filepath_genome) as genetic_code_source:
                        are_we_good = False
                        contig_position = f">{contig_name}"
                        for line in genetic_code_source:
                            if contig_position in line: #Initiate genetic code saving
                                are_we_good = True
                                continue
                            elif are_we_good == True: #Save code from contig, stop on the next one
                                if line.startswith(">"):
                                        are_we_good = False
                                else:
                                    genetic_code += line.strip()
                    criteria = end - start
                    if criteria > 0: #If the start and end position of the overlapping sequence on the genome has a positive difference
                        begin, finish = list_of_overlap_points[0], list_of_overlap_points[1]
                        overlap_sequence += genetic_code[begin-1:finish-1]
                        fullinfo = f">{fixedname}\n{overlap_sequence}\n"
                        with open(concat_file, 'a') as dumpy:
                            dumpy.write(fullinfo)
                    else:
                        if criteria < 0: #If the BLAST hit is on the complementary strand of the genome assembly contig -> translate
                            list_of_overlap_points.sort()
                            begin, finish = list_of_overlap_points[0], list_of_overlap_points[1]
                            overlap_sequence += genetic_code[begin-1:finish-1]
                            reversedoverlapsequence = overlap_sequence[::-1]
                            translatefirst = reversedoverlapsequence.replace("T", "Y").replace("G", "M")
                            orthodoxtranslation = translatefirst.replace("A", "T").replace("C", "G").replace("Y", "A").replace("M", "C")
                            fullinfo = f">{fixedname}\n{orthodoxtranslation}\n"
                            with open(concat_file, 'a') as mega_dumpy:
                                mega_dumpy.write(fullinfo)
    # Begin MSA production, following the for loop for creating the concat file                        
    mouse = pathdictionary["ITERATIVE"]
    targetoutput = os.path.join(mouse, f"{querygene}_HOMOLOGUES_ALIGNED.fasta")
    print("> Preparing MSA file...")
    if not os.path.exists(targetoutput):
        print("> Generating output directory...")
        with open(targetoutput, "w") as write:
            pass
    subprocess.run(f"mafft --maxiterate 1000 --thread -1 {concat_file} > {targetoutput}", shell=True)
    print("> MSA complete, moving onto NWK generation.")
    output_for_nwk = os.path.join(pathdictionary["RESULTS"], f"{querygene}")
    subprocess.run(f"raxml-ng --msa {targetoutput} --model GTR+G+I --prefix {output_for_nwk}", shell=True)
    for file in os.listdir(pathdictionary["RESULTS"]):
       if "bestTree" in file and f"{querygene}" in file:
          filepa = os.path.join(pathdictionary["RESULTS"], file)
          rename = os.path.join(pathdictionary["RESULTS"], f"{querygene}.nwk")
          os.rename(filepa, rename)
    print("> NWK complete!")
