import os 
import json
import subprocess
############### CHECK DIRECTORIES
desktopdirectory = "/home/student"
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
#################GET ID NAMES#######################
tsvpathlist = [file for file in os.listdir(pathdictionary["USERFILES"]) if file.endswith(".tsv") and "Name" in file]
NameIDDict = {}
with open(os.path.join(pathdictionary["USERFILES"] ,tsvpathlist[0]), 'r') as input:
    for liner in input:
        barcode, name = liner.strip().split("\t")
        NameIDDict[barcode] = name
###### SOURCE FULL QUERY SEQUENCE PATHS
queryfolder = pathdictionary["USERFILES"]
querysequencepaths = {}
for file in os.listdir(queryfolder):
    if file.endswith('.fna'):  
        path = os.path.join(queryfolder, file)
        queryname = os.path.splitext(file)[0]
        querysequencepaths[queryname] = path
###############################################################
blastresults = pathdictionary["BLAST-Search-Results"]
tempcalcdict = {}
for query in querysequencepaths.keys():
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
                                templist.append(overlapend)
                                templist.append(overlapstart)
                                tempcalcdict[genomename] = templist
                        else:
                            tempcalcdict[genomename] = ""
    newfilename = f"{query}_overlaps.json"
    newfilepath = os.path.join(pathdictionary["App-Files"], newfilename)
    with open(newfilepath, "w") as newfile:
        json.dump(tempcalcdict, newfile)  
########## MODULATED OVERLAP RIP COMPONENT
if len(os.listdir(pathdictionary["Overlap-Genomes"])) <= 1:
    buffer = "n" * 5
    for genomefile in os.listdir(pathdictionary["USERFILES"]):
        if genomefile.endswith(".fasta"):
            inb_tween = genomefile.strip(".fasta")
            genome = NameIDDict[inb_tween]
            newname = f"{genome}_overlap_sequence.fasta"
            newgenomeassembly = os.path.join(pathdictionary["Overlap-Genomes"], newname)
            file_path = os.path.join(pathdictionary["USERFILES"], genomefile)
            for querygene in querysequencepaths.keys():
                genetic_code = ""
                print(f"> Working with {querygene}...searching for overlap sequence on {genome}")
                for appfile in os.listdir(pathdictionary["App-Files"]):
                    furtherfilepath = os.path.join(pathdictionary["App-Files"], appfile)
                    if querygene in furtherfilepath:
                        if "overlap" in furtherfilepath:
                            overlap_sequence = ""
                            with open(furtherfilepath) as data:
                                overlapdict = json.load(data)                                
                            if not overlapdict[genome]:
                                with open(newgenomeassembly, "a") as overlaprip:
                                    overlaprip.write("")
                                    print(f"> {querygene} not found on genome.")
                            else:
                                indexes = overlapdict[genome]
                                contig_name, start, end = indexes[0], indexes[1], indexes[2]
                                list_of_overlap_points = []
                                list_of_overlap_points.append(indexes[1])
                                list_of_overlap_points.append(indexes[2])
                                with open(file_path) as genetic_code_source:
                                    are_we_good = False
                                    contig_position = f">{contig_name}"
                                    for line in genetic_code_source:
                                        if contig_position in line:
                                            are_we_good = True
                                            continue
                                        elif are_we_good == True:
                                            if line.startswith(">"):
                                                 are_we_good = False
                                            else:
                                             genetic_code += line.strip()
                                criteria = end - start
                                if criteria > 0:
                                    begin, finish = list_of_overlap_points[0], list_of_overlap_points[1]
                                    overlap_sequence += genetic_code[begin-1:finish-1]
                                    with open(newgenomeassembly, "a") as overlaprip:
                                        overlaprip.write(buffer + overlap_sequence + buffer)
                                        print(f"> Added {querygene} overlap sequence!")
                                else:
                                    if criteria < 0:
                                        list_of_overlap_points.sort()
                                        begin, finish = list_of_overlap_points[0], list_of_overlap_points[1]
                                        overlap_sequence += genetic_code[begin-1:finish-1]
                                        reversedoverlapsequence = overlap_sequence[::-1]
                                        translatefirst = reversedoverlapsequence.replace("T", "Y").replace("G", "M")
                                        orthodoxtranslation = translatefirst.replace("A", "T").replace("C", "G").replace("Y", "A").replace("M", "C")
                                        with open(newgenomeassembly, "a") as overlaprip:
                                            overlaprip.write(buffer + orthodoxtranslation + buffer) 
                                            print(f"> Added {querygene} translated overlap sequence!") 
else:
  print("> 'Overlap genomes' already exist, moving on...")
###### GENERATE COMBINED FILE OF SEQUENCES
concatfile = os.path.join(pathdictionary["RESULTS"], "COMBINED.fasta")
if not os.path.exists(concatfile):
    print("> Combining the overlap genomes...please wait.")
    for file in os.listdir(pathdictionary["Overlap-Genomes"]):
        for bacterium in NameIDDict.values():
            thepath = os.path.join(pathdictionary["Overlap-Genomes"], file)
            if bacterium in thepath:
                with open(thepath) as source:
                    largesequence = source.readline()
                    sequencename = f">{bacterium}"
                    fixedname = sequencename.replace(" ", "_")
                    fullinfo = f"{fixedname}\n{largesequence}\n"
                    with open(concatfile, 'a') as final:
                        final.write(fullinfo)
    print("> Done!")
else:
    print("> No problems with the requirements for MSA generation, moving on...")
######## GENERATE MSA FILE
results = pathdictionary["RESULTS"]
targetoutput = os.path.join(results, "ALIGNED.fasta")
print("> Preparing MSA file...")
if not os.path.exists(targetoutput):
    print("> Generating output directory...")
    with open(targetoutput, "w") as write:
        pass
if os.path.getsize(targetoutput) < 10: 
 subprocess.run(f"mafft --maxiterate 1000 --thread -1 {concatfile} > {targetoutput}", shell=True)
else:
    pass
########## GENERATE NWK FILE
teliko = results + "/NEWTREE"
print("> MSA complete, moving onto NWK generation.")
if not teliko in os.listdir(results):
   subprocess.run(f"raxml-ng --search1 --msa {targetoutput} --model GTR+G+I --prefix {teliko}", shell=True)
   for file in os.listdir(results):
       if "bestTree" in file:
          filepa = os.path.join(results, file)
          rename = os.path.join(results, "NEWTREE.nwk")
          os.rename(filepa, rename)
else:
    print("> The program has identified you have already ran this module with prior results. \nTo avoid problems with data generation, either move youre previous results to another directory, \nor choose to run by deleting previous results (Y/y).")
 
#raxml-ng --msa /Users/chrismitsacopoulos/Desktop/XopAM_HOMOLOGUES_ALIGNED.fasta --model GTR+G+I --prefix /Users/chrismitsacopoulos/Desktop/XopAM_TREE
#mafft --maxiterate 1000 --thread -1 
