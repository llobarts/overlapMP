#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
readpsiblast.py

Description
-----------
This program reads the information obtained by executing PSI-BLAST. This
information is written in format 7. After reading the information, it saves the
number accession and the name of the species of each hit in two text files. 
Also it gives information about the convergence process (number of 
iterations until convergence) and the number of proteins until filtering, after 
filtering and the number of species with the corresponding evalue.

Author: Esmeralda Garcia Legarda
Date: 06/09/2021
Version 0.2
"""
import os
import pandas as pd
import re
import subprocess


def listblastfiles(directory:str) -> list:
    """List the blast files located in a directory.

    Args:
        directory (str): directory

    Returns:
        list: paths of the files
    """

    list_files = []
    for file in os.listdir(directory):
        if file.endswith(".blast"):
            list_files.append(file)

    return list_files


def listnumaccessfiles(directory:str) -> list:
    """List the blast files located in a directory.

    Args:
        directory (str): directory

    Returns:
        list: paths of the files
    """

    list_files = []
    for file in os.listdir(directory):
        if file.endswith("_info.txt"):
            list_files.append(file)

    return list_files


def read_accessnumandsp(file_name:str, identity:int=30, len_align:int=100) -> list:
    """
    From a blast file, it extracts the accession numbers.
    
    Parameters
    ----------
    file_name: str
        Name of the file (blast file, format 7)
    identity: int
        Minimum of percentage of identity
    len_align: int
        Minimum of alignment length

    Returns
    -------
    list:
        List of accession numbers
    """
    printing = '$2 "\t" $12 "\t" $13'
    cmd = "grep ^[^#] " + str(file_name) + ".blast | awk  -F '\t' '{ \
if ($3>" + str(identity) + " && $4>" + str(len_align) + ") \
{ print " + printing + "}}'"
    cmd = ["/bin/sh", "-c", cmd]
    accessnumbers = subprocess.run(cmd,
        shell=False, check=True,
        stdout=subprocess.PIPE,
        universal_newlines=True,
        executable='/bin/bash').stdout

    list_access = accessnumbers.rstrip().split("\n")
    list_access = list(set(list_access))

    return list_access


def save_species(name_file:str, numspecies:list):
    """
    It saves the species in a text file.

    Parameters
    ----------
    name_file: str
        Name of the file which is being processed (+ path)
    numspecies: list
        List of species
    """
    
    file_name = "".join([name_file, ".txt"])

    with open(file_name, "w") as file_writing:
    
        for i in numspecies:
            file_writing.write("".join([i, "\n"]))


def isvirus(info:str) -> bool:
    """
    It reads if the taxonomic name corresponds to a virus or not by detecting
    the presence of the word "virus".
    
    Parameters
    ----------
    taxname: str
        Taxonomic name

    Returns
    -------
    Bool
        True if it is a virus, False if it is not.

    """
    isvirus = False

    if len(re.findall("[Vv]irus", info)) >= 0:
        isvirus = True
    
    return isvirus


def save_numaccesionandsp(name_file:str, list_access:list):
    """
    It saves the accession numbers in a text file.

    Parameters
    ----------
    name_file: str
        Name of the file which is being processed
    list_access: list
        List of numbers of accession
    """
    name_file = "_".join(name_file.split("_")[:-2]) + "_info.txt"

    with open(name_file, "w") as file_num:
    
        for num in list_access:
            file_num.write(num)
            file_num.write("\n")


def get_stats(list_files:list, num_iter:int, evalues:list):
    """
    It saves a csv file which contains the number of proteins until filtering, 
    after filtering and the number of species with the corresponding evalue.

    Parameters
    ----------
    list_files: list
        List of files
    evalues: list
        List of used evalues
    prot_init: list
        List of quantities of proteins, considering all proteins
    prot_filter: list
        List of quantities of proteins, considering filtered proteins list
        (by alignment lenght and identity)
    prot_species: list
        List of quantities of species
    prot_virus: list
        List of quantities of viral species
    """
    dict = {"iterations":num_iter, "evalue":evalues}
    df_stats = pd.DataFrame.from_dict(dict, orient="index",columns=list_files)
    
    df_stats = df_stats.transpose()
    df_stats.to_csv("Results/Reports/data/stats.csv")


def converge(file_name:str, max_iter:int = 25) -> int:
    """
    It detects if the search had converged. If the convergence was not achieved,
    the maximum iteration limit, set by the user, can be indicated.

    Parameters
    ----------
    file_name: str
        Name of the file which is being processed
    max_iter: int
        Number of the maximum number of interation (in case convergence was not
        achieved)

    Returns
    -------
    int:
        Number of iterations
    """

    cmd = "cat " + file_name + ".blast | grep 'Search has CONVERGED!' | wc -l"
    cmd = ["/bin/sh", "-c", cmd]
    converge_file = int(subprocess.run(cmd,
            shell=False, check=True,
            stdout=subprocess.PIPE,
            universal_newlines=True,
            executable='/bin/bash').stdout.strip())

    if converge_file == 1:
        num_iter = num_iterations(file_name)

    else:
        num_iter = max_iter
    
    return num_iter


def num_iterations(file_name:str) -> int:
    """
    It reads the number of iterations used until converging.

    Parameters
    ----------
    file_name: str

    Returns
    -------
    int:
        Number of iterations until convergence
    
    """
    cmd = "cat " + file_name + ".blast | grep 'Iteration'"
    cmd = ["/bin/sh", "-c", cmd]
    BLASTiter = subprocess.run(cmd,
                shell=False, check=True,
                stdout=subprocess.PIPE,
                universal_newlines=True,
                executable='/bin/bash').stdout.strip()

    num_iter = int(BLASTiter.split(" ")[-1])

    return num_iter

def retrieve_accessnumandsp(all_name:str) -> list: # , list
    """
    It reads accession numbers taking into account the filters identity and
    length of the alignment and without, and it saves the first one list.

    Args:
        all_name (str): Name of the file (blast file, format 7)

    Returns:
        list: Accession numbers with filter
        list: Accession numbers without filter
    """
    list_access = read_accessnumandsp(all_name)
    save_numaccesionandsp(all_name, list_access)


def compareresults(list_directories:list) -> str:
    """
    It compares the number of number of accession in resulting of trying
    different evalues.

    Args:
        list_directories (list): folders to be compared.

    Returns:
        str: Path with the maximum number of number accessions or the 
        conservative
    """
    dict_directory = {}
    for directory in list_directories:
        count = 0
        
        for file_numacces in listnumaccessfiles(directory):
            with open(directory + "/" + file_numacces, 'r') as filep:
                num_lines = sum(1 for line in filep if line.rstrip())
                count += num_lines
        dict_directory[directory] = count

    max_path = max(dict_directory, key=dict_directory.get)

    return max_path


def createunion(list_directories:list, max_path:str) -> list:
    """
    It generates a list with all the accession numbers, doing an union from
    the different files.

    Args:
        list_directories (list): Directories with different files of accession
        numbers
        max_path (str): Path of the directory where there is the maximum number
        of accession numbers

    Returns:
        dict:  All accession numbers (key: file name of every familys)
    """
    list_max = {}
    for file_numacces in listnumaccessfiles(max_path):
        key = "".join(file_numacces.split("_")[0])
        with open(max_path + "/" + file_numacces, 'r') as file_max:
            list_max[key] = file_max.read().split("\n")[:-1]

    for directory in list_directories:
        print(max_path)
        for i, file_numacces in enumerate(listnumaccessfiles(directory)):
            key = "".join(file_numacces.split("_")[0])
            print(directory + "/" + file_numacces)
            with open(directory + "/" + file_numacces, 'r') as filep:
                list_num = filep.read().split("\n")[:-1]
            for n in list_num:
                if n not in list_max[key]:
                    list_max[key].append(n)
    
    return list_max


def saveunion(list_numaccess:dict):
    """It saves the files related to the union of accession numbers. 

    Args:
        list_numaccess (dict): [description]
    """
    for name_file in list_numaccess.keys():
        pathunion = "Results/temp/psiblastunion/" + name_file
        save_species(pathunion, list_numaccess[name_file])


def readpsiblast(list_directories:list):
    """
    It reads the results of the PSI-BLAST, calculates some stats to summarise
     all the results. Is saves stats, accession numbers, virus and other species
     in separated files.

     Args:
         list_directories (list): directories that contain different results of
         blast, probably because of the different used evalues.
     """

    # Initial setting
    evalues = []
    num_iter = []

    for directory in list_directories:
        
        # Setting directories and names
        num_dir = float(directory.split("_")[-2])
        list_files = listblastfiles(directory)

        for file_name in list_files:
            all_name = directory + "/" + file_name[:-6]

        #     # Retrieving accession numbers
        #     retrieve_accessnumandsp(all_name)
           
            # Getting stats
            evalues.append(num_dir)   
            num_itera = converge(all_name)
            num_iter.append(num_itera)

    # Get the directory with maximum number of accession numbers
    max_path = compareresults(list_directories)

    # Create a folder which generates the union of the different files
    os.mkdir("Results/temp/psiblastunion")
    list_numaccess = createunion(list_directories, max_path)
    saveunion(list_numaccess)

    # Summarizing values of the process and results
    get_stats(list_files, num_iter, evalues)

if __name__ == "__main__":
    readpsiblast(["Results/temp/psiblast_0.001_25"])