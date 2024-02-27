from ete3 import Tree, NodeStyle, TreeStyle, RectFace, TextFace, faces
import os
import json
import pandas as pd
# CHECK DIRS
desktopdirectory = os.path.expanduser("~/Desktop")
pathforworkingfolder = os.path.join(desktopdirectory, "MINAGIASGENE")
if not os.path.exists(pathforworkingfolder):
    os.mkdir(pathforworkingfolder)
subfolders = ["USERFILES", "RESULTS", "BLAST-Databases",
              "BLAST-Search-Results", "App-Files", "QUAST", "PHYLOGENETICS"]
pathdictionary = {"MINAGIAS": pathforworkingfolder}
for folder in subfolders:
    subfolderpath = os.path.join(pathforworkingfolder, folder)
    if not os.path.exists(subfolderpath):
        os.makedirs(subfolderpath)
    pathdictionary[folder] = subfolderpath
# SORT HIGHEST PRESENCE
blastresults = pathdictionary["App-Files"]
tempcalcdict = {}
for file in os.listdir(blastresults):
    if file.endswith(".json"):
        path = os.path.join(blastresults, file)
        if "Presence" in path:
            with open(path) as data:
                presencedict = json.load(data)
                gonlist = [y for y in presencedict.values()]
                genename = str(gonlist[0]).strip(" present").strip(
                    " inconclusive").strip(" absent")
                presencelist = [x for x in presencedict.values()]
                numberofpresence = presencelist.count(f"{genename} present")
                tempcalcdict[genename] = numberofpresence
dictrearranged = {}
descendinglist = []
descendinglist = [x for x in tempcalcdict.values()]
descendinglist.sort(reverse=True)
for x in descendinglist:
    for key, value in tempcalcdict.items():
        if value == x:
            dictrearranged[key] = value
# GET PATHOVAR DICT
tsvpathlist = [file for file in os.listdir(
    pathdictionary["USERFILES"]) if file.endswith(".tsv") and "Path" in file]
frame = pd.read_csv(os.path.join(
    pathdictionary["USERFILES"], tsvpathlist[0]), sep='\t')
bacloc = {}
for row in frame.values:
    bacloc[row[0]] = row[1]
pathovars_of_interest = {}
for key_, value_ in bacloc.items():
    if value_ == "euvesicatoria" or value_ == "unknown":
        pathovars_of_interest[key_] = value_
# READ THE PRESENCE JSONS, USED IN FOR LOOP ITERATION
def readresultsjson(queryname, results):
    presencedict = {}
    for file in os.listdir(results):
        if file.endswith(".json"):
            if f"{queryname}_Presence" in file:
                jsonpath = os.path.join(results, file)
    with open(jsonpath, newline='') as jsonfile:
        presencedict = json.load(jsonfile)
    return presencedict
# MSA TREE ANNOTATION, VISUALISATION, PRINTING
x = 0
newickpath = "/Users/chrismitsacopoulos/Desktop/REVISED.nwk"
newtree = Tree(newickpath)
tree = TreeStyle()
nodes = newtree.get_leaves()
topnode = nodes[0]
printname = "SNP_presence_tree_BLASTn.png"
print_dir = os.path.join(pathdictionary["PHYLOGENETICS"], printname)
# Automatic appending should follow, colours for legend and annotation should match
for name in dictrearranged.keys():
    x += 1
    presencedict = readresultsjson(name, pathdictionary["App-Files"])
    tree.show_leaf_name = True
    tree.force_topology = True
    tree.draw_guiding_lines = True
    tree.guiding_lines_color = '#000000'
    for node in newtree.traverse():
        treestyle = NodeStyle()
        treestyle.vt_line_width = 200
        treestyle.hz_line_width = 200
        new_style = NodeStyle()
        new_style["bgcolor"] = "light blue"
        second_style = NodeStyle()
        second_style["bgcolor"] = "light yellow"
        node.set_style(treestyle)
        if node.is_leaf():
            if node != topnode:
                if node.name in presencedict:
                    node.add_features(genepresence=presencedict.get(node.name))
                    nodepresence = presencedict.get(node.name)
                    if "present" in nodepresence:
                        heatmap = "black"
                    elif "inconclusive" in nodepresence:
                        heatmap = "grey"
                    else:
                        heatmap = "white"
                    heatmapface = RectFace(
                        width=9, height=9, fgcolor="black", bgcolor=heatmap)
                    node.add_face(heatmapface, column=x, position="aligned")
                    if node.name in pathovars_of_interest.keys():
                        if pathovars_of_interest[node.name] == "euvesicatoria":
                            node.set_style(new_style)
                        elif pathovars_of_interest[node.name] == "unknown":
                            node.set_style(second_style)
            else:
                if node.name in pathovars_of_interest.keys():
                    if pathovars_of_interest[node.name] == "euvesicatoria":
                            node.set_style(new_style)
                    elif pathovars_of_interest[node.name] == "unknown":
                            node.set_style(second_style)
                node.add_features(genepresence=presencedict.get(node.name))
                nodepresence = presencedict.get(node.name)
                if "present" in nodepresence:
                    heatmap = "black"
                elif "inconclusive" in nodepresence:
                    heatmap = "grey"
                else:
                    heatmap = "white"
                topnodetext = TextFace(f"{name}", fstyle="bold", fsize=12)
                topnodetext.rotation = 90
                topnode.add_face(topnodetext, column=x, position="aligned")
                heatmapface = RectFace(
                    width=9, height=9, fgcolor="black", bgcolor=heatmap)
                node.add_face(heatmapface, column=x, position="aligned")
                spacer = RectFace(width=9, height=9,
                                  fgcolor="white", bgcolor="white")
                topnode.add_face(spacer, column=x, position="aligned")
                spacer = RectFace(width=9, height=9,
                                  fgcolor="white", bgcolor="white")
                topnode.add_face(spacer, column=x, position="aligned")
                spacer = RectFace(width=9, height=9,
                                  fgcolor="white", bgcolor="white")
                topnode.add_face(spacer, column=x, position="aligned")
tree.legend.add_face(
    RectFace(width=6, height=6, fgcolor="black", bgcolor="black"), column=x)
tree.legend.add_face(TextFace("Present: Near Identical Hit (>=97%)", fsize=10,
                     ftype="Times New Roman"), column=x+1)
tree.legend.add_face(
    RectFace(width=6, height=6, fgcolor="black", bgcolor="grey"), column=x)
tree.legend.add_face(TextFace("Present: Query/Genome overlap significant (<97%)",
                     fsize=10, ftype="Times New Roman"), column=x+1)
tree.legend.add_face(
    RectFace(width=6, height=6, fgcolor="black", bgcolor="white"), column=x)
tree.legend.add_face(TextFace("Not present: No overlap between Query/Genome (N/A)", fsize=10,
                     ftype="Times New Roman"), column=x+1)
tree.legend_position = 3
tree.show_scale = True
newtree.render(print_dir, dpi=3000, tree_style=tree)
