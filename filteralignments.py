#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os

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


def read_aligntext(name):

    sequences = []
    with open(name, "r") as readenfile:
        seqs = readenfile.read().strip().split(">")[1:]
    
    for i in seqs:
        fastaseq = i.strip().split("\n")
        sequences.append([fastaseq[0], "".join(fastaseq[1:])])

    return sequences


def remove_empties(seqs, threshold):
    delete_elements = []
    size = len(seqs[0][-1])
    
    for i, seq in enumerate(seqs):
        if seq[-1].count("-") > threshold*size:
            delete_elements.append(i)
    
    for j in reversed(sorted(delete_elements)):
        del seqs[j]


def remove_gaps(seqs, threshold):
    size = len(seqs[0][-1])
    delete_elements = [0]*size
    num_seqs = len(seqs)

    for seq in seqs:
        for j, letter in enumerate(seq[-1]):
            if letter == "-":
               delete_elements[j] += 1

    gaps = [pos for pos, m in enumerate(delete_elements)  if m > threshold *num_seqs]

    for position in reversed(sorted(gaps)):
        for num in range(len(seqs)):
            seqs[num] = [seqs[num][0], seqs[num][-1][:position] + \
                seqs[num][-1][position+1:]]


def remove_equal_sequences(sequences):
    delete_elements = []

    for seq1 in range(len(sequences)):
        onecomp = sequences[seq1][-1]
        onecomp_sp = sequences[seq1][0].split("|")[2]
        
        for seq2 in range(len(sequences)):
            if  seq1 < seq2 and onecomp == sequences[seq2][-1] and\
                onecomp_sp == sequences[seq2][0].split("|")[2]:
                
                delete_elements.append(seq2)

    for j in reversed(sorted(list(set(delete_elements)))):
        del sequences[j]


def savealignment(sequences, name):

    with open(name, 'w') as towrite:
        for seq in sequences:
            towrite.write(">")
            towrite.write("\n".join(seq))
            towrite.write("\n")


def main():

    dir_alignments = "/Users/esmeralda/Documents/TFM/2_alignments/CLANS/results_6"
    filesalign = listfiles(dir_alignments, "fa")
    
    for alignm in filesalign:
        seqs = [] 
        seqs = read_aligntext(alignm)

        remove_equal_sequences(seqs)
        remove_gaps(seqs, 0.98)
        # remove_empties(seqs, 0.4)
        
        savealignment(seqs, alignm[:-3] + "_filtered.fa")
    

if __name__ == "__main__":
    main()
