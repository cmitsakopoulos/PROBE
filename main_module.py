import os
import click
import subprocess
import time
import sys
import multiprocessing

APPNAME = "PROBE"
VERSION = "1.0 beta"
desktopdirectory = os.path.expanduser("~/Desktop")
pathforworkingfolder = os.path.join(desktopdirectory, "PROBE")
if not os.path.exists(pathforworkingfolder):
    os.mkdir(pathforworkingfolder)
subfolders = ["USERFILES", "RESULTS", "BLAST-Databases", "BLAST-Search-Results",
              "App-Files", "QUAST", "ITERATIVE", "PHYLOGENETICS"]
pathdictionary = {"PROBE": pathforworkingfolder}
for folder in subfolders:
    subfolderpath = os.path.join(pathforworkingfolder, folder)
    if not os.path.exists(subfolderpath):
        os.makedirs(subfolderpath)
    pathdictionary[folder] = subfolderpath

def maincomponent():
    while True:
        click.echo(click.style(f"""Welcome to {APPNAME} v{VERSION} \n
        Choose your desired function accordingly by typing the function's number and pressing enter. \n
        -> List of automated functions;
            1) Make Databases of '.fasta' genomes
            2) Run BLAST with '.fna' queries on genomes
            3) Run QUAST on genomes
            4) Generate '.fasta' MSA through maximum likelyhood(MAFFT) -> Generate GTR+G based tree in NWK(RAxML)
            5) Generate pathovar based graphical representation for data (need extra .tsv with sample names and location) \n """, fg="green", bold=True))
        answer = click.prompt("Respond with desired function number, 'e' to shut down: ", type=click.Choice(
            ["1", "2", "3", "4", "5", "e"], case_sensitive=False))
        if answer == "1":
            DATABASES()
            break
        elif answer == "2":
            ONLYBLAST()
            break
        elif answer == "3":
            QUAST()
            break
        elif answer == "4":
            ONLYTREE()
            break
        elif answer == "5":
            location_analysis()
            break
        elif answer == "e":
            click.echo(click.style(
                f"Exiting application, thank you for choosing {APPNAME}!", fg="blue", bold=True))
            break

def DATABASES():
    programmescwd = os.path.dirname(os.path.realpath(__file__))
    contents = os.listdir(pathdictionary["BLAST-Databases"])
    if len(contents) > 1:
        click.echo(click.style(
            "BLAST Databases already exist, do you want to remove them and start anew? (MUST)", fg="green", bold=True))
        answer = click.prompt("Enter (Y/y) to delete previous results and re-run,\n type 'b' to go back to options, 'e' to shut down: ",
                              type=click.Choice(["y", "n", "e", "b"], case_sensitive=False))
        if answer == "Y" or answer == "y":
            click.echo(click.style(
                "> Deleting previous BLAST databases...", fg="red", bold=True))
            for file in os.listdir(pathdictionary["BLAST-Databases"]):
                filepath = os.path.join(
                    pathdictionary["BLAST-Databases"], file)
                if os.path.isfile(filepath):
                    print(f"> Deleting {filepath}")
                    os.remove(filepath)
            subprocess.run(['python3', 'make_databases.py'],
                           cwd=programmescwd, check=True)
            click.echo(click.style("> Done!", fg="blue", bold=True))
            maincomponent()
        elif answer == "e":
            click.echo(click.style(
                f"> Exiting application, thank you for choosing {APPNAME}!", fg="blue", bold=True))
            sys.exit()
        elif answer == "b":
            click.echo(click.style("> Going back...", fg="blue", bold=True))
            time.sleep(2)
            maincomponent()
    else:
        click.echo(click.style(
            "> Starting database creation...please wait", fg="blue", bold=True))
        subprocess.run(['python3', 'make_databases.py'],
                       cwd=programmescwd, check=True)
        click.echo(click.style("> Done!", fg="green", bold=True))
        maincomponent()

def QUAST():
    programmescwd = os.path.dirname(os.path.realpath(__file__))
    contents = os.listdir(pathdictionary["QUAST"])
    if len(contents) > 1:
        click.echo(click.style(
            "QUAST has been run previously, do you want to delete prior results and start anew?", fg="green", bold=True))
        answer = click.prompt("Enter (Y/y) to delete previous results and re-run,\n type 'b' to go back to options, 'e' to shut down: ",
                              type=click.Choice(["y", "n", "e", "b"], case_sensitive=False))
        if answer == "Y" or answer == "y":
            click.echo(click.style(
                "> Deleting previous QUAST results...", fg="red", bold=True))
            for file in os.listdir(pathdictionary["QUAST"]):
                filepath = os.path.join(pathdictionary["QUAST"], file)
                if os.path.isfile(filepath):
                    print(f"> Deleting {filepath}")
                    os.remove(filepath)
            click.echo(click.style(
                "> Starting QUAST...", fg="blue", bold=True))
            subprocess.run(['python3', 'run_quast.py'],
                           cwd=programmescwd, check=True)
            click.echo(click.style("> Done!", fg="green", bold=True))
            maincomponent()
        elif answer == "e":
            click.echo(click.style(
                f"> Exiting application, thank you for choosing {APPNAME}!", fg="blue", bold=True))
            sys.exit()
        elif answer == "b":
            click.echo(click.style("> Going back...", fg="blue", bold=True))
            time.sleep(2)
            maincomponent()
    else:
        click.echo(click.style("> Starting QUAST...", fg="blue", bold=True))
        subprocess.run(['python3', 'run_quast.py'],
                       cwd=programmescwd, check=True)
        click.echo(click.style("> Done!", fg="green", bold=True))
        maincomponent()

def ONLYBLAST():
    programmescwd = os.path.dirname(os.path.realpath(__file__))
    contents = os.listdir(pathdictionary["BLAST-Search-Results"])
    if len(contents) > 1:
        click.echo(click.style(
            "BLAST has been run previously, you can run and keep prior results, or start anew:", fg="green", bold=True))
        answer = click.prompt("Enter (D/d) to delete previous results and re-run,\n enter (K/k) to run without deleting previous results, \n type 'b' to go back to options or 'e' to shut down: ",
                              type=click.Choice(["k", "d", "e", "b"], case_sensitive=False))
        if answer == "D" or answer == "d":
            click.echo(click.style(
                "> Deleting previous QUAST results...", fg="red", bold=True))
            for file in os.listdir(pathdictionary["BLAST-Search-Results"]):
                filepath = os.path.join(
                    pathdictionary["BLAST-Search-Results"], file)
                if os.path.isfile(filepath):
                    print(f"> Deleting {filepath}")
                    os.remove(filepath)
            click.echo(click.style(
                "> Preparing for BLAST...", fg="blue", bold=True))
            try:
                checker = subprocess.run(['python3', 'prepare_dirs_blast.py'], cwd=programmescwd, capture_output=True, check=True)
                if 144 in checker.stdout:
                    print("> Returning to options...")
                    maincomponent()
                else:
                    print("> No TSV errors, moving on...")
            except subprocess.CalledProcessError as e:
                print(f"OOPS! > An error occurred: {e}")
            click.echo(click.style(
                "> Running BLAST...please wait", fg="blue", bold=True))
            #### RUN SIMULTANEOUS BLAST
            run_simult()
            ##### CONTINUE
            click.echo(click.style(
                "> Interpreting BLAST results...please wait", fg="blue", bold=True))
            subprocess.run(['python3', 'calculate_results.py'], cwd=programmescwd)
            click.echo(click.style(
                "> Done! Going back to options...", fg="green", bold=True))
            time.sleep(2)
            maincomponent()
        elif answer == "K" or answer == 'k':
            click.echo(click.style(
                "> Preparing for BLAST...", fg="blue", bold=True))
            try:
                checker = subprocess.run(['python3', 'prepare_dirs_blast.py'], cwd=programmescwd, capture_output=True, check=True)
                if 144 in checker.stdout:
                    print("> Returning to options...")
                    maincomponent()
                else:
                    print("> No TSV errors, moving on...")
                    
            except subprocess.CalledProcessError as e:
                print(f"OOPS! > An error occurred: {e}")
            click.echo(click.style(
                "> Running BLAST...please wait", fg="blue", bold=True))
            #### RUN SIMULTANEOUS BLAST
            run_simult()
            ##### CONTINUE
            subprocess.run(
                ['python3', 'calculate_results.py'], cwd=programmescwd, check=True)
            click.echo(click.style(
                "> Done! Going back to options...", fg="green", bold=True))
            time.sleep(2)
            maincomponent()
        elif answer == "e":
            click.echo(click.style(
                f"Exiting application, thank you for choosing {APPNAME}!", fg="blue", bold=True))
            sys.exit()
        elif answer == "b":
            click.echo(click.style("Going back...", fg="blue", bold=True))
            time.sleep(2)
            maincomponent()
    else:
        click.echo(click.style("Preparing for BLAST...", fg="blue", bold=True))
        subprocess.run(['python3', 'prepare_dirs_blast.py'],
                       cwd=programmescwd, check=True)
        click.echo(click.style(
            "> Running BLAST, please wait...", fg="blue", bold=True))
        try:
            checker = subprocess.run(['python3', 'prepare_dirs_blast.py'], cwd=programmescwd, capture_output=True, check=True)
            if 144 in checker.stdout:
                print("> Returning to options...")
                maincomponent()
            else:
                print("> No TSV errors, moving on...")
        except subprocess.CalledProcessError as e:
                print(f"OOPS! > An error occurred: {e}")
        #### RUN SIMULTANEOUS BLAST
        run_simult()
        ##### CONTINUE
        click.echo(click.style(
            "> Moving onto result interpretation...", fg="blue", bold=True))
        subprocess.run(['python3', 'calculate_results.py'], cwd=programmescwd)
        click.echo(click.style(
            "> Done! Going back to options...", fg="green", bold=True))
        maincomponent()
            
def ONLYTREE():
    programmescwd = os.path.dirname(os.path.realpath(__file__))
    click.echo(click.style(
        "The phylogenetics function has been run before, \n do you want to delete prior results and start again?", fg="green", bold=True))
    answer = click.prompt("By entering (I/i), you can create separate trees for each of your genome-gene homologues.\n Else, type 'b' to go back to options, 'e' to shut down: ",
                            type=click.Choice(["I", "i", "n", "e", "b"], case_sensitive=False))
    if answer == "I" or answer == "i":
        click.echo(click.style(
            "> Deleting previous phylo results...", fg="red", bold=True))
        for file in os.listdir(pathdictionary["ITERATIVE"]):
            filepath = os.path.join(pathdictionary["ITERATIVE"], file)
            if os.path.isfile(filepath):
                print(f"> Deleting {filepath}")
                os.remove(filepath)
        for file in os.listdir(pathdictionary["PHYLOGENETICS"]):
            filepath = os.path.join(pathdictionary["PHYLOGENETICS"], file)
            if os.path.isfile(filepath):
                print(f"> Deleting {filepath}")
                os.remove(filepath)
        click.echo(click.style(
            "> Working on MSA and NWK generation, please wait...", fg="blue", bold=True))
        accelerate_phylogenomics()
        click.echo(click.style(
            "> Done! Going back to options...", fg="green", bold=True))
        maincomponent()
    elif answer == "K" or answer == "k":
        click.echo(click.style(
            "> Looking through previous results, please wait...", fg="blue", bold=True))
        subprocess.run(['python3', 'produce_MSA_NWK_complete.py'],
                        cwd=programmescwd, check=True)
        click.echo(click.style(
            "> Done! Going back to options...", fg="green", bold=True))
        maincomponent()
    elif answer == "e":
        click.echo(click.style(
            f"> Exiting application, thank you for choosing {APPNAME}!", fg="blue", bold=True))
        sys.exit()
    elif answer == "b":
        click.echo(click.style("> Going back...", fg="blue", bold=True))
        time.sleep(1)
        maincomponent()

def location_analysis():
    programmescwd = os.path.dirname(os.path.realpath(__file__))
    contents = [x for x in os.listdir(
        pathdictionary["RESULTS"]) if x.endswith(".pdf")]
    if ".pdf" in contents:
        click.echo(click.style(
            "The location graph function has been run before,\
                  \n do you want to delete prior results and start again?",
            fg="green",
            bold=True))
        answer = click.prompt("Enter (Y/y) to overwrite previous results and add to them \n type 'b' to go back to options, 'e' to shut down: ",
                              type=click.Choice(["y", "n", "e", "b"], case_sensitive=False))
        if answer == "Y" or answer == "y":
            click.echo(click.style("> Making graphs...", fg="blue", bold=True))
            subprocess.run(['python3', 'location_dict_barchart.py'],
                           cwd=programmescwd, check=True)
            click.echo(click.style(
                "> Done! Going back to options...", fg="green", bold=True))
            maincomponent()
        elif answer == "e":
            click.echo(click.style(
                f"> Exiting application, thank you for choosing {APPNAME}!", fg="blue", bold=True))
            sys.exit()
        elif answer == "b":
            click.echo(click.style("> Going back...", fg="blue", bold=True))
            time.sleep(2)
            maincomponent()
    else:
        click.echo(click.style(
            "> Working on drawing graphs...please wait", fg="blue", bold=True))
        subprocess.run(['python3', 'location_dict_barchart.py'],
                       cwd=programmescwd, check=True)
        click.echo(click.style(
            "> Done! Going back to options...", fg="green", bold=True))
        maincomponent()

######## MULTIPROCESSING COMPONENTS
def first_blast_module():
    programmescwd = os.path.dirname(os.path.realpath(__file__))
    insanity = subprocess.Popen(['python3', 'blast_module_eins.py'], 
                            cwd=programmescwd)
    insanity.wait()

def second_blast_module():
    programmescwd = os.path.dirname(os.path.realpath(__file__))
    insanity2 = subprocess.Popen(['python3', 'blast_module_zwei.py'],
                        cwd=programmescwd)
    insanity2.wait()

def run_simult():
    blast1 = multiprocessing.Process(target=first_blast_module)
    blast2 = multiprocessing.Process(target=second_blast_module)
    blast1.start()
    blast2.start()
    blast1.join()
    blast2.join()

def phylo_first():
    programmescwd = os.path.dirname(os.path.realpath(__file__))
    insanity_extreme = subprocess.Popen(['python3', 'produce_MSA_NWK_iterative.py'], 
                            cwd=programmescwd)
    insanity_extreme.wait()

def phylo_second():
    programmescwd = os.path.dirname(os.path.realpath(__file__))
    insanity_crazy = subprocess.Popen(['python3', 'produce_MSA_NWK_iterative_2.py'], 
                            cwd=programmescwd)
    insanity_crazy.wait()

def accelerate_phylogenomics():
    phylo1 = multiprocessing.Process(target=phylo_first)
    phylo2 = multiprocessing.Process(target=phylo_second)
    phylo1.start()
    phylo2.start()
    phylo1.join()
    phylo2.join()

if __name__ == "__main__":
    maincomponent()