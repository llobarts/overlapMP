# Overlap
This is the repository dedicated to the scripts and data from the TFM of Bioinformatics (UV) called **"STRUCTURAL AND EVOLUTIONARY CONSTRAINTS ON VIRAL MOVEMENT PROTEINS ENCODED BY OVERLAPPING GENES IN VIRAL GENOMES"**.
. 

Inside of the folder *Overlapping*, there is a set of useful scripts to retrieve, get information and classify the overlapping 
and non-overlapping sequences of proteins.

- *overlapping.py* is the main script.

- The folder *Example* contains an example file and the file with the sequences of movement proteins used for
the final project of the master.

- The folder *Models* includes scripts that run PSI-BLAST, read the information retrieved and give information about the search and
the overlap results.

- The folder *buildb* contains the results, for each family and all together. 

The folder *psiblast_results* stores the results and some information about the PSI-BLAST search for each family of movement proteins. The 
summary of the results of selection and some plots are generated in the folder *Selection*.

Other files are:

- filteralignments.py: it allow to filter the alignments by removing repeted sequences, columns with a percentage of gaps and/or 
sequences that represented a percentage below a limit (set in 0.98 and 0.4 in the present study).
- parse_clans.py: this script transforms the 3D topology of the clusters obtained by CLANS into a dendrogram and a heatmap.
- p19_alignment.fa: 33 sequences of p19 proteins. They have their corresponding sequence of movement protein in tombusMP_alignment.fa
- tombusMP_alignment.fa: 33 sequences of movement proteins. They have their corresponding sequence of p19 protein in p19_alignment.fa
- psiblast_report.html: report of the results of using different parameters in PSI-Blast and sequences for each family of 30k movement proteins.
- report_DMPfold.html: report of the DMPfold usage, another deep learning-based de novo protein modelling software.
- Dali_results.html: comparative analyses of the progressive use of a thicker multiple sequence alignment using DMPfold. 
- report_selection.html: report of the different methods used for the detection of selection in *Tombusvirus* movement protein. 

