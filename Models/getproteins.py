#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
getproteins.py

Description
-----------
This program obtains the protein sequences and saves them and also classifies
them depending on if they are movement proteins, if they overlaps, and which 
protein overlaps with them. All this is indicated in every header of the fasta
sequences that are saved.

Autor: Esmeralda Garcia Legarda
Date: 27/09/2021
Version 0.1
"""
import os
import subprocess
import sys


def listfiles(directory:str, end:str) -> list:
    """List the files with located in a directory.

    Args:
        directory (str): directory

    Returns:
        list: paths of the files
    """

    list_files = []
    for file in os.listdir(directory):
        if file.endswith(end):
            list_files.append(directory + "/" + file)
    
    return sorted(list_files)


def readidproteins(file_name:str):
    """
    It reads the accession numbers.
    
    Parameters
    ----------
    file_name: str
        Name of the file (blast file, format 6)

    Returns
    -------
    list:
        List of accession numbers
    """
    cmd = "tail +2 " + file_name + "| awk -F ',' '{print $2}'"
    cmd = ["/bin/sh", "-c", cmd]
    id_proteins = subprocess.run(cmd,
        shell=False, check=True,
        stdout=subprocess.PIPE,
        universal_newlines=True,
        executable='/bin/bash').stdout.strip().split("\n")
    
    return id_proteins


def retrieve_seqsNUCL(num:str):
    """
    From an accession number it returns the corresponding sequence of 
    nucleotides (fasta).

    Parameters
    ----------
    num: str
        Accession number.

    Returns
    -------
    str:
        Sequence of nucleotides (fasta format).
    """

    cmd = "elink -db protein -id " + str(num) + " -target nuccore | efetch \
-format fasta"

    cmd = ["/bin/sh", "-c", cmd]

    fasta = subprocess.run(cmd,
        shell=False, check=True,
        stdout=subprocess.PIPE,
        universal_newlines=True,
        executable='/bin/bash').stdout

    return fasta

    
def retrieve_seqsPROT(num:str):
    """
    From an accession numbers it returns the corresponding sequence of 
    aminoacids (fasta).

    Parameters
    ----------
    num: str
        Accession number.

    Returns
    -------
    str:
        Sequence of aminoacids (fasta format).
    """
    cmd = "efetch -db protein -id " + str(num) + " -format fasta"

    cmd = ["/bin/sh", "-c", cmd]
    fasta = subprocess.run(cmd,
        shell=False, check=True,
        stdout=subprocess.PIPE,
        universal_newlines=True,
        executable='/bin/bash').stdout.split("\n")
    
    fasta = "".join(fasta[1:-1])
    
    return fasta

 
def processinfogenome(infogenome:list):
    info = ">" + "|".join(infogenome).replace(" ", "_")
    
    return info


def savefilefromlist(listinfo, name):
    
    with open(name, "w") as fastafile:
        for i in listinfo:
            fastafile.write(i + "\n")


def savedifferentfiles(nameoutput):
    patterns = ["|Y|Y|", "|Y|N|", "|N|Y|"]
    extensions = ["_query_overlap", "_query_nooverlap", "_other_overlapping"]

    for i in range(len(patterns)):

        cmd = "".join(["grep '", patterns[i], "' -A 1 ", nameoutput,
        " |  grep -v -- '^--$'", " > ", nameoutput[:-3], extensions[i], ".fa"])
        cmd = ["/bin/sh", "-c", cmd]
        subprocess.run(cmd,
            shell=False, check=False,
            stdout=subprocess.PIPE,
            universal_newlines=True,
            executable='/bin/bash')


def savefiles(fasta_seq, notfound, newpath, csvfile):

    nameoutput = "".join([newpath, "/database_", csvfile.split("/")[-1][:-12],
        ".fa"])
    savefilefromlist(fasta_seq, nameoutput)
    savedifferentfiles(nameoutput)
    
    nameoutput = "".join([newpath, "/notfound_", csvfile.split("/")[-1][:-12],
        ".txt"])
    savefilefromlist(notfound, nameoutput)
    

def getproteins():

    ###########################################################################
    #                           GET THE FASTA FILES                           #
    ###########################################################################
    
    # Initial settings
    newpath = "Results/Reports/data/builddb_last/"
    os.mkdir(newpath)
    directory_gnoverlap = "Results/Reports/data/genomesoverlap"
    filesgenomes = listfiles(directory_gnoverlap, "csv")
    notfound = []

    for csvfile in filesgenomes:

        list_fasta = []
        with open(csvfile, "r") as readingfile:
            readingfile.readline()
            for line in readingfile:
                infogenome = line.split(",")[1:]
                idprotein = infogenome[0]
                
                info = processinfogenome(infogenome)
                sequence = retrieve_seqsPROT(idprotein)

                if len(sequence) < 0:
                    notfound.append(idprotein)
                else:
                    list_fasta.append("".join([info, sequence]))

        savefiles(list_fasta, notfound, newpath, csvfile)


if __name__ == "__main__":
    getproteins()

