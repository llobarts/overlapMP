"""
This module contains the very first functions to run PSI-BLAST.
"""
import os
from argparse import ArgumentParser
import subprocess


def parsearguments() -> dict:
    """
    It reads the arguments from de command line.

    Returns:
        dict: Contains the values for the parameters
    """
    # Args processing
    parser = ArgumentParser()
    parser.add_argument(
    '-f', '--fasta',
    type=str,
    dest='fastafile',
    action='store',
    required=False,
    default= "Example/initialseq.fa",
    help='File that contains the families in format fasta,\
        where headers are their names')
    parser.add_argument(
    '-d', '--database',
    type=str,
    dest='db',
    action='store',
    required=True,
    help='Absolute path of the database')
    parser.add_argument(
    '-i', '--maxiterations',
    type=int,
    dest='max_iter',
    action='store',
    required=False,
    default=25,
    help='Maximum number of iterations for runnining PSI-BLAST (default 25)')
    parser.add_argument(
    '-e', '--evalue',
    dest='evalue',
    action='store',
    required=False,
    type=float,
    default=0.001,
    help='Evalue to run PSIBLAST (default 0.001)')
    parser.add_argument(
    '-t', '--threads',
    type=int,
    dest='threads',
    action='store',
    required=False,
    default=8,
    help='Number of threads (dafault 8)')
    parser.add_argument(
    '-p', '--psiblastpath',
    type=str,
    dest='psiblastpath',
    action='store',
    required=True,
    help='Path of the binary executable of psi-blast')
    parser.add_argument(
    '-n', '--nameproject',
    type=str,
    dest='nameproject',
    action='store',
    required=False,
    default= 'newproject',
    help='Name of the project')

    args = parser.parse_args()
    
    return vars(args)


def writecommand(psiblast, evalue, max_iter, db, threads) -> list:
    """
    Joins the parameters introduced by the user and creates the complete
    command.

    Args:
        evalue (float): Evalues to run PSI-BLAST
        max_iter (int): Maximum number of iterations 
        db (str): Path of the database to run PSI-BLAST
        threads (int): Threads to run PSI-BLAST

    Returns:
        list: Complete command
    """
    dbpath = "/".join(db.split("/")[0:-1])
    dbname = db.split("/")[-1]
    cmd = ["/bin/sh", "-c", "Models/dopsiblast.sh" +\
        " -p " + psiblast +\
        " -e " + str(evalue) + " -i " + str(max_iter) +\
        " -d " + dbname + " -D " + dbpath  + \
        " -t " + str(threads)]

    return cmd


def splitfasta(fastafile:str) -> list: # , list
    """
    Generates a fasta file for every fasta sequences that the user introduces.
    They had to be separated by '"---"'.

    Args:
        fastafile (int): Path where the input fasta file is located 

    Returns:
        names (list): 
        paths (list): Contains the paths of the resulting fasta files
    """
    sep = "---"
    names = []
    paths = []
    output_path = "Results/temp/Fastafiles"
    os.mkdir(output_path)
    
    with open(fastafile, "r") as readfile:
        fastaseq = readfile.read().split(sep)

        for seq in fastaseq:
            name = seq.strip().split("\n")[0][1:].replace("|", "")
            path = output_path + "/" + name + ".fa"
            names.append(name)
            paths.append(path)
            
            with open(path, 'w') as result:
                result.write(seq.strip())

    return sorted(names), sorted(paths)

def runpsiblast(cmd:str):
    subprocess.run(cmd, shell=False, check=True)
