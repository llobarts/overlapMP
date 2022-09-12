#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
overlapping.py

This program allows you to 

Author: Esmeralda García Legarda
Date: 22/10/2021
"""
from Models import starting, askgenomefeatures , getproteins, readpsiblast, \
    reporting
import os


def main():

    ############################################################################
    #                    Process the args and default values                   #
    ############################################################################
    variables = starting.parsearguments()

    filepath = variables['fastafile']
    psiblast = variables['psiblastpath']
    db = variables['db']
    max_iter = variables['max_iter']
    evalue = variables['evalue']
    threads = variables['threads']
    nameproject = variables['nameproject']

    ############################################################################
    #            Create needed directories and process initial files           #
    ############################################################################
    
    # # Create temporary directory
    # os.mkdir("Results/temp")
    
    # # Create data directory
    # os.mkdir("Results/Reports/data")


    # # Splitting the fasta file into multiple fastafiles, one per family
    # namesfamilies, pathsfamilies = starting.splitfasta(filepath)


    ############################################################################
    #                               Run PSI-BLAST                              #
    ############################################################################
    
    # newpath = "Results/temp/psiblast_" + str(evalue) + "_" + str(max_iter)
    # os.mkdir(newpath)
    # path_results = []

    # cmd = starting.writecommand(psiblast, evalue, max_iter, db, threads)
    # starting.runpsiblast(cmd)
        
    # path_results.append(newpath)


    ############################################################################
    #                       Read the results of PSI-BLAST                      #
    ############################################################################

    # readpsiblast.readpsiblast(path_results)

    ############################################################################
    #                               Ask overlap                                #
    ############################################################################
    
    # try:
    #     os.mkdir("Results/Reports/data/genomesoverlap")
    # except FileExistsError:
    #     print("Carpeta genomesoverlap/ existe")
    
    # askgenomefeatures.askoverlap() # Generates csv with the information
    # askgenomefeatures.joincsv()


    # ############################################################################
    # #                              Get proteins                                #
    # ############################################################################

    try:
        os.mkdir("Results/Reports/data/builddb")
    except FileExistsError:
        print("Directory builddb/ already exists")
    try:
        os.mkdir("Results/Reports/data/builddb/infodb")
    except FileExistsError:
        print("Directory infodb/ already exists")

    getproteins.getproteins()


    # ############################################################################
    # #                                Get report                                #
    # ############################################################################
    
    reporting.generatereport(nameproject)
    reporting.openreportbrowser(nameproject)


if __name__ == "__main__": 
    main()
