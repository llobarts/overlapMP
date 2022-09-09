#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import squareform

def readinfofile(filename):
    listobjects = [[],[],[]]
    wordinit_list = ["<seq>", "<seqgroups>", "<pos>"]
    wordend_list = ["</seq>", "</seqgroups>", "</pos>"]

    with open(filename, 'r') as fileopen:

        for i in range(len(wordinit_list)):
            getinfochunk(wordinit_list[i], wordend_list[i], listobjects[i], 
                                                    fileopen)
    
    return listobjects


def getinfochunk(wordinit, wordend, listobjects, fileopen):
    end = True
    line = fileopen.readline().strip()

    while not line.startswith(wordend):
        if line.startswith(wordinit):
            end = False
            
        elif not end:
            listobjects.append(line.strip().split(" "))

        line = fileopen.readline().strip()


def getgroups(groups):
    sep_groups = []

    for i in groups:
        for l in i:
            if l.startswith("numbers"):
                sep_groups.append(l.replace("numbers=", "").split(";"))

    return sep_groups


def assigngroups(groups, size):
    
    numbers_group = []

    for row in range(size):
        for group in groups:
            if str(row) in group:
                numbers_group.append(groups.index(group) + 1)

        if len(numbers_group) == row:
            numbers_group.append(len(groups) + 1)

    return numbers_group


def joininfo(info, coord, numbers_group):

    info_df = pd.DataFrame([ x for x in info if info.index(x) % 2 == 0 ], 
                                                    columns = ["header"])

    info_df_aux = pd.DataFrame([ x for x in info if info.index(x) % 2 != 0 ], 
                                                    columns = ["sequence"])

    coords_df = pd.DataFrame(coord, columns = ["number", "x", "y", "z"])
    del coords_df["number"]

    group_df = pd.DataFrame(numbers_group, columns = ["group"])
    info_df = pd.concat([info_df, info_df_aux, coords_df, group_df], axis=1)
    
    return info_df


def savecoord(info_df, name):
    
    info_df[['x', 'y', 'z']].to_csv(name)


def calculate_centroids(info_df):
    centroids = {}
    for i in range(1, int(info_df['group'].max()) + 1):
        # REMOVE: JUST FOR ANALYSING TBSVP22
        if i != 12:
            info_df_group = info_df[info_df["group"] == i]
        else:
            info_df_group = info_df[info_df["group"] == i]
            info_df_group = info_df_group[info_df_group["header"].str.contains('TBSVP22')]
            
        centroids[i] = [np.mean([float(i) for i in list(info_df_group["x"])]), 
            np.mean([float(i) for i in list(info_df_group["y"])]), 
            np.mean([float(i) for i in list(info_df_group["z"])])]

    return centroids


def calculate_distances(centroids):
    matrix_distances = np.zeros(shape=(len(centroids),len(centroids)))
    
    for centroid_A in centroids.keys():
        for centroid_B in centroids.keys():
            matrix_distances[centroid_A - 1][centroid_B - 1] =  \
            np.sqrt((centroids[centroid_B][0] - centroids[centroid_A][0])**2 + \
            (centroids[centroid_B][1] - centroids[centroid_A][1])**2 + \
                (centroids[centroid_B][2] - centroids[centroid_A][2])**2)
    
    names = ['G1','G2','G3', 'G4', 'G5', 'G6', 'G7', 'G8', 'G9', 'G10', 'G11', 'G12']
    distances = pd.DataFrame(matrix_distances, columns = names, index = names)

    return distances

def plot_distances(matrix_distances):
    
    fig, ax = plt.subplots(figsize=(10, 8))
    # mask
    mask = np.triu(np.ones_like(matrix_distances, dtype=np.bool_))
    # adjust mask and df
    mask = mask[1:, :-1]
    corr = matrix_distances.iloc[1:,:-1].copy()
    # plot heatmap
    sns.heatmap(corr, 
                mask=mask, 
                annot=True, 
                fmt=".2f", 
                cmap='Reds',
                vmin=0, 
                vmax=matrix_distances.max().max(), 
                cbar_kws={"shrink": .8})
    # yticks
    plt.yticks(rotation=0)
    figure = plt.gcf()
    figure.set_size_inches(8, 6.8)
    plt.savefig("/Users/esmeralda/Documents/TFM/figures/heatmap.pdf", format="pdf", bbox_inches="tight")
    plt.show()


def extract_list_distances(matrix_distances):
    matrix_aux = pd.DataFrame(np.copy(matrix_distances))
    matrix_aux.values[np.triu_indices(len(matrix_distances))] = np.nan
    data = [ele for n in matrix_aux.columns for ele in matrix_aux[n].values.tolist()]

    return data


def plot_histogram(data):
    plt.style.use('ggplot')
    sns.kdeplot(data)
    plt.show()


def calculate_dispersion(info_df, centroids):

    dispersions = []
    
    for i in range(1, int(info_df['group'].max()) + 1):
        info_df_group = info_df[info_df["group"] == i]
        dispersion_sum = 0
        for element in list((info_df_group["x"]).index):

            dispersion_sum += (float(info_df_group["x"][element]) - centroids[i][0])**2 + \
                 (float(info_df_group["y"][element]) - centroids[i][1])**2 + \
                (float(info_df_group["z"][element]) - centroids[i][2])**2

        dispersions.append(np.sqrt(dispersion_sum/((len(info_df_group["x"]))**2)))

    return dispersions


def show_groupdata(groups_names, centroids, dispersions):
    
    names = ['G1','G2','G3','G4','G5','G6','G7','G8','G9','G10','G11','G12']
    groupdata = pd.DataFrame(centroids, index=["x","y","z"]).T
    groupdata['dispersion'] = dispersions
    groupdata["family"] = groups_names
    groupdata.insert(0, 'family', groupdata.pop('family'))
    groupdata.index = names
    print("INFORMATION GROUPS: CENTROIDS AND DISPERSION\n", groupdata)


def save_distances_to_edgelist(matrix_distances, name):
    edges = set()
    
    for i in list(matrix_distances.columns):
        for j in list(matrix_distances.columns):
            if i < j and matrix_distances[i][j] != 'NaN':
                edges.add(i+"\t"+j + "\t" + str(matrix_distances[i][j])+"\n")

    with open(name, 'w') as edgefile:
        for m in edges:
            edgefile.write(m)


def assingnamegroup(info_df):
    names = []

    for i in range(1, int(info_df['group'].max()) + 1):
        info_df_group = info_df[info_df["group"] == i]
        total_names = {}
        for element in list((info_df_group["header"]).index):
            namegroup = info_df_group["header"][element].split("|")[1]
            if namegroup not in total_names.keys():
                total_names[namegroup] = 1
            else:
                total_names[namegroup] += 1
        
        maxname = max(total_names, key=total_names.get)
        names.append("G" + str(i) + "-" + maxname)
    
    return names


def plotdendrogram(matrix_distances, groups_names):

    dists = squareform(matrix_distances.to_numpy())
    linkage_matrix = linkage(dists, "single")
    dendrogram(linkage_matrix, labels=groups_names)
    figure = plt.gcf()
    figure.set_size_inches(22, 11)
    plt.rcParams['lines.linewidth']=3
    plt.rcParams['axes.facecolor']='white'
    plt.rcParams['xtick.color']='black'
    plt.style.use("seaborn-whitegrid")
    plt.savefig("/Users/esmeralda/Documents/TFM/figures/dendro.pdf", format="pdf", bbox_inches="tight")
    plt.show()


def savegroups(info_df):
    
    for i in range(1, info_df["group"].max() + 1):
        df_group = info_df[info_df["group"] == i]
        with open("G"+ str(i) + '.fa', 'w') as groupfile:
            for index, row in df_group.iterrows():
                groupfile.write(row['header'] + "\n")
                groupfile.write(row['sequence'] + "\n")


def main():
    
    # Reading CLANS file 
    path = "/Users/esmeralda/Documents/TFM/2_alignments/CLANS/"
    os.chdir(path)
    name =  "all_clusters.clans"
    listobjects = readinfofile(name)

    info = listobjects[0]
    groups = listobjects[1]
    coord = listobjects[2]

    # Split information of CLANS file and crete a dataframe
    sep_groups = getgroups(groups)
    numbers_group = assigngroups(sep_groups, len(coord))
    info_df = joininfo(info, coord, numbers_group)
    groups_names = assingnamegroup(info_df)
    
    # Saving coordinates in csv
    # savecoord(info_df, "coord_clans.csv")

    # Calculating distances of centroids
    centroids = calculate_centroids(info_df)
    dispersions = calculate_dispersion(info_df, centroids)
    matrix_distances = calculate_distances(centroids)


    distances_list = extract_list_distances(matrix_distances)
    average = np.nanmean(distances_list)
    
    # Printing and ploting
    print(info_df)
    show_groupdata(groups_names, centroids, dispersions)
    print("DISTANCES BETWEEN CENTROIDS", matrix_distances)
    print("Average:", f'{average:.2f}')
    plot_distances(matrix_distances)
    plot_histogram(distances_list)
    plotdendrogram(matrix_distances, groups_names)

    # Save groups
    savegroups(info_df)


if __name__ == "__main__":
    main()