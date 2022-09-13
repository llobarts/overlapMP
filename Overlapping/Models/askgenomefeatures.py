import os
import pandas as pd
import subprocess

def idProteintoidGenome(id_protein:str):
    """It returns the ID of the Genome that corresponds to the id of an ID
    of a protein (accession numbers).

    Args:
        id_protein (str): ID of a protein (accession number)

    Returns:
        str: ID genome (accession number)
    """
    cmd =  "elink -db protein -id '" + id_protein + "' -target nuccore | efetch\
         -format acc"
    try:
        cmd = ["/bin/sh", "-c", cmd]
        id_genome = subprocess.run(cmd,
        shell=False, check=True,
        stdout=subprocess.PIPE,
        universal_newlines=True,
        executable='/bin/bash').stdout.strip()
    except:
        try:
            cmd = ["/bin/sh", "-c", cmd]
            id_genome = subprocess.run(cmd,
        shell=False, check=True,
        stdout=subprocess.PIPE,
        universal_newlines=True,
        executable='/bin/bash').stdout.strip()
        except:
            id_genome = ""
    
    if id_genome.split("\n"):
        id_genome = id_genome.split("\n")[0]
        
    return id_genome


def getinfogenome(id_genome:str, nameGFF:str):
    """With an ID of a genome, it retrieves the information about it, in 
    GFF format, and writes it in a file (nameGFF).

    Args:
        id_genome (str): ID of a genome
        nameGFF (str): processed name 

    Returns:
        bool: if the process was successfully.
    """
    information = False
    cmd =  "esearch -db nuccore -query '" + id_genome + "' | efetch -format ft" 
    cmd = ["/bin/sh", "-c", cmd]
    ft = subprocess.run(cmd,
        shell=False, check=False,
        stdout=subprocess.PIPE,
        universal_newlines=True,
        executable='/bin/bash').stdout.strip()

    if len(ft) > 0:
        information = True
        with open(nameGFF, "w") as filename:
            filename.write(ft)

    return information


def getinfoprotein(id_protein:str, ft:str) -> str:
    """It gets the information of the protein with the indicated ID.

    Args:
        id_protein (str): protein ID (accession number)
        ft (str): name of the file GFF

    Returns:
        str: information of the protein
    """
    infoprotein = ""
    num = 5
    limit = num*5

    cmd = "grep -B '" + str(num) + "' '" + id_protein + "' '" + ft + "'"
    cmd = ["/bin/sh", "-c", cmd]

    infoprotein = subprocess.run(cmd,
            shell=False, check=False,
            stdout=subprocess.PIPE,
            universal_newlines=True,
            executable='/bin/bash').stdout.strip()

    while not infoprotein.find("\tCDS") >= 0 and num < limit:
        cmd = "grep -B '" + str(num) + "' '" + id_protein + "' '" + ft + "'"
        cmd = ["/bin/sh", "-c", cmd]

        infoprotein = subprocess.run(cmd,
            shell=False, check=False,
            stdout=subprocess.PIPE,
            universal_newlines=True,
            executable='/bin/bash').stdout.strip()
        num += 2

    aux = False
    temp_list = infoprotein.split("\n")
    n = - 1
    while n > - len(temp_list) and not aux:
        if temp_list[n].find("\tCDS") >= 0:
            aux = True
        else:
            n -= 1

    new = temp_list[n:]
    newinfoprotein = "\n".join(new)

    return newinfoprotein


def processinfo(infoprotein:str) -> dict:
    """It process the information of a protein.

    Args:
        infoprotein (str): information from a GFF file.

    Returns:
        dict: contains features of the protein.
    """
    chain = "+"
    features = list(filter(bool, infoprotein.split("\t")))
    element = 0
    end = -1
    while element < len(features):
        f = features[element]
        if f.find("CDS") == 0:
            start = int(features[element - 2].replace("<", ""))
        
        if f.find("product") == 0:
            if features[element-1].find("CDS") >= 0 and end == -1:
                end = int(features[element-2].replace(">", ""))
            elif end == -1:
                end = int(features[element - 1].replace(">", ""))
            
            if infoprotein.find("prot_desc") < 0:
                functionname = features[element + 1].strip()

        if f.find("prot_desc") == 0:
            functionname = features[element+1].strip()

            if (not any([x == "product" for x in features])) and end == -1:
                end = int(features[element-2].replace(">", ""))
            elif end == -1:
                end = int(features[element-4].replace(">", ""))

        element += 1
    
    if end < start:
        chain = "-"
    featuresgenome = [min(start, end), max(start, end), chain,
                                                     functionname.strip()]
    return featuresgenome


def listoverlapping(positions:list, start:str, end:str):
    """Generates a list of positions that are overlapping with the gene of the
    protein of interest.

    Args:
        positions (list): total positions
        start (str): start of the gene of the protein of interest.
        end (str): end of the gene of the protein of interest.

    Returns:
        list: positions
    """
    list_overlaps = []
    for line in positions:
        start_init = int(line.strip().split("\t")[0].replace("<", ""))
        end_init = int(line.strip().split("\t")[1].replace(">", ""))
        start_other = min(start_init, end_init)
        end_other = max(start_init, end_init)

        same = (start == start_other) and (end == end_other)
        if not same:

            if (start_other <= start and\
                                end_other > start) or\
                                (start_other < end and\
                                end_other >= end) or\
                                (start_other <= start and\
                                end_other >= end) or\
                                (start_other > start and\
                                end_other < end):

                list_overlaps.append(line)

    return list_overlaps


def getpositions(nameft:str):
    """
    It obtains the positions of the CDSs.
    
    Args:
        nameft (str): name of the GFF file.

    Returns:
        list: positions of the CDSs.

    """
    cmd = """cat """ + '"' + nameft + '"' + """ | awk -F '\t' ' BEGIN { count=1 } 
    $3=="CDS" { print $0 ;  count=0 ; }
    $4~"product" || $4~"prot_desc" || $4~"transl_table" { count=1 ; }
    count==0 && $3!="CDS" {print $0}'"""
    cmd = ["/bin/sh", "-c", cmd]
    positions = subprocess.run(cmd,
        shell=False, check=True,
        stdout=subprocess.PIPE,
        universal_newlines=True,
        executable='/bin/bash').stdout.strip().split("\n")

    cont = 0
    infopositions = []

    while cont < len(positions):

        if cont == 0:
            division_start = positions[cont].strip().split("\t")[0]
        elif positions[cont].find("CDS") >= 0:
            if positions[cont-1].find("CDS") >= 0:
                division_end = positions[cont-1].strip().split("\t")[-2]
            else:
                division_end = positions[cont-1].strip().split("\t")[-1]
            
            infopositions.append(division_start + "\t" + division_end)
            division_start = positions[cont].strip().split("\t")[0]

        if cont == len(positions) - 1:
            division_end = positions[cont].strip().split("\t")[1]
            infopositions.append(division_start + "\t" + division_end)

        cont += 1
    
    return infopositions


def getinfopositions(i:str, nameGFF:str):
    """Obtains the information for a CDS.

    Args:
        i (str): CDS positions.
        nameGFF (str): file name

    Returns:
        list: features of a CDS
    """
    cmd = "grep -A 5 '" + "\t" + i.split("\t")[-1] + "' '" + nameGFF + "'"
    cmd = ["/bin/sh", "-c", cmd]

    infopositions = subprocess.run(cmd,
        shell=False, check=True,
        stdout=subprocess.PIPE,
        universal_newlines=True,
        executable='/bin/bash').stdout.strip()

    if infopositions.find("prot_desc") != -1 and \
        infopositions.find("product") != -1:
        if infopositions.find("product") > infopositions.find("prot_desc"):
            infopositions = infopositions.split("prot_desc")[1:]
        else:
            infopositions = infopositions.split("product")[1:]
    elif infopositions.find("product") != -1:
        infopositions = infopositions.split("product")[1:]
    else:
        infopositions = infopositions.split("prot_desc")[1:]
    
    infopositions = "|".join("".join(infopositions).split("|")[:2])

    return "\t".join([i, "CDS", "product",  infopositions])


def getidprotein(positions:str):
    """Asks the protein ID of the overlapping gene.

    Args:
        positions (str): information of the overlapping gene.

    Returns:
        str: "NA" or protein ID (accesion number).
    """
   
    list_protein = positions.split("|")

    if len(list_protein) > 1:
        id_protein = list_protein[1]
    elif positions.find("pseudo") and positions.find("RefSeq:"):
        id_protein = "_".join([positions.split("RefSeq:")[1], "pseudo"]) 
    else:
        id_protein = "NA"
    
    return id_protein


def lookforoverlap(nameGFF:str, start:str, end:str):
    """It looks for the overlap of the genes that codify proteins.

    Args:
        nameGFF (str): name of the GFF
        start (str): start of the gene of the protein of interest
        end (str): end of the gene of the protein of interest

    Returns:
        list: contains the information of overlap/s, if there is overlap.
    """
    infooverlap = []
    positions = getpositions(nameGFF)
    list_overlaps = listoverlapping(positions, start, end)
    
    if list_overlaps:

        for i in list_overlaps:
            info = getinfopositions(i, nameGFF)
            infooverlap.append(info)
    
    return infooverlap 


def getcoordinates(infolist_overlap:list, start, end):
    """It gets the coordinates of the overlapping gene ordered.

    Args:
        infolist_overlap (_type_): information of the overlap
        start (int): start of the gene
        end (int): end of the gene

    Returns:
        tuple: start and end, ordered.
    """
    o_start = max(infolist_overlap[0], start)
    o_end = min(infolist_overlap[1], end)
   
    return o_start, o_end


def saveinfospecies(proteins:list, family:str):
    """It saves the information of the species that have overlapping genes.

    Args:
        proteins (list): information of the proteins and its overlap.
        family (str): name of the family
    """
    namefile = "Results/Reports/data/genomesoverlap/" + family + "_overlap.csv"
    df = pd.DataFrame(proteins)
    df.to_csv(namefile)


def removeinfogenome(nameGFF:str):
    """It removes the GFF file.

    Args:
        nameGFF (str): name of the GFF file.
    """
    os.remove(nameGFF)


def askgenome(id_protein:str, nameGFF:str, familyandspecies:list, 
    genesgenomes:list):
    """This function does question if the protein of interest is overlapping
    or not.

    Args:
        id_protein (str): protein ID (accession number).
        nameGFF (str): name of the file in format GFF.
        familyandspecies (list): list with the information of family, species 
        and kingdom.
        genesgenomes (list): list of information with a structure.

    Returns:
        bool: if the process of asking has been succesfully. 
    """

    infoprotein = getinfoprotein(id_protein, nameGFF)
    
    if len(infoprotein) != 0:
        infolist = processinfo(infoprotein)
        infooverlapping = lookforoverlap(nameGFF, 
                        start=infolist[0], end=infolist[1])

        if infooverlapping:

            for infooverlap in infooverlapping:
                
                infolist_overlap = processinfo(infooverlap)
                id_overlapping = getidprotein(infooverlap)
                o_start, o_end = getcoordinates(infolist_overlap, 
                                        start=infolist[0],end=infolist[1])
                genesgenomes.append([id_protein] + familyandspecies + infolist +\
                                             ["Y", "Y", o_start, o_end])
                
                genesgenomes.append([id_overlapping] + familyandspecies + \
                    infolist_overlap + ["N", "Y", o_start, o_end])
        else:
            genesgenomes.append([id_protein] + familyandspecies + infolist + \
                                                    ["Y", "N", "0", "0"])
        obtained = True
    
    else:
        obtained = False

    return obtained


def savenotfound(notfound:list, family:str):
    """It saves in a file the accession numbers of the proteins of interest 
    that were not found in the database.

    Args:
        notfound (list): list of accession numbers.
        family (str): name of the family.
    """
    namefile = "Results/Reports/data/genomesoverlap/" + family + "_notfound.txt"
    with open(namefile, "w") as filetowrite:
        for species in notfound:
            filetowrite.write(species + "\n")


def read_sp_id(line:str):
    """Reads a line of the file of species, id and kingdom.

    Args:
        line (str): line of the read file

    Returns:
        tuple: species, id and kingdom
    """
    linesplit = line.split("\t")
    species = linesplit[1]
    id_protein = linesplit[0]
    kingdom = linesplit[2].strip()

    return species, id_protein, kingdom


def listcsvfiles(path:str):
    """Lists the csv files that are in a path.

    Args:
        path (str): path/folder.

    Returns:
        list: contains the paths (str) of the csv files.
    """
    listfiles = []
    for i in os.listdir(path):
        if i.endswith(".csv"):
            listfiles.append("Results/Reports/data/genomesoverlap/" + i)

    return listfiles


def joincsv():
    """The files the folder genomesoverlap are joined to an only file called
    allinfo.csv.
    """
    pathdir = "Results/Reports/data/genomesoverlap"
    pathscsv = " ".join(listcsvfiles(pathdir))

    cmd = "".join(["cat ", pathscsv," '>' ", pathdir, "/allinfo.csv"])
    cmd = ["/bin/sh", "-c", cmd]

    infopositions = subprocess.run(cmd,
        shell=False, check=True,
        stdout=subprocess.PIPE,
        universal_newlines=True,
        executable='/bin/bash').stdout.strip()


def askoverlap():
    """This is the main function. It reads the results of the union, searches
    the information and saves information about overlapping.
    """
    directory = "Results/temp/psiblastunion/"
    list_directories  = sorted(os.listdir(directory))

    for namefile in list_directories:
        family = namefile[:-4]
        genesgenomes = []
        notfound = []
        with open(directory + namefile, "r") as readfile:
            
            for line in readfile:
                species, id_protein, kingdom = read_sp_id(line)
                familyandspecies = [family, species, kingdom]
                
                nameGFF = "Results/Reports/data/genomesoverlap/" + \
                    species.replace(" ","_") + ".txt"
                id_genome = idProteintoidGenome(id_protein)
                
                if len(id_genome) > 0:
                    isinformation = getinfogenome(id_genome, nameGFF)

                    if isinformation:
                        obtained = askgenome(id_protein, nameGFF,
                        familyandspecies, genesgenomes)
                        
                        if not obtained:
                            notfound.append(id_protein)
                        removeinfogenome(nameGFF)
                    else:
                        notfound.append(id_protein)
                else:
                    notfound.append(id_protein)
        
        if len(genesgenomes) > 0:
            saveinfospecies(genesgenomes, family)
        savenotfound(notfound, family)
