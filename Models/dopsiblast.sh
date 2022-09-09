#!/bin/bash
# This program generates results of running PSI-BLAST
# Reading FLAGS
while getopts :p:e:i:d:D:t: flag
do
    case ${flag} in
        p) psiblast=${OPTARG};;
        e) evalue=${OPTARG};;
        i) iter=${OPTARG};;
        d) dbname=${OPTARG};;
        D) dbpath=${OPTARG};;
        t) threads=${OPTARG};;
    esac
done

# Initial settings
format='7 qaccver saccver pident length mismatch gapopen qstart qend sstart send stitle ssciname sskingdom'
currpath=$(pwd)
evalue_float=$(expr $evalue)
num_itera=$(expr $iter)

# Create directory
OUTPUT_DIR=${currpath}/Results/temp/psiblast_${evalue}_${iter}
cd $dbpath

for INPUT in $(ls ${currpath}/Results/temp/Fastafiles)
do
# Run PSI-BLAST
OUTPUT=${OUTPUT_DIR}/${INPUT%.*}_${evalue}_${iter}
code=$($psiblast -query ${currpath}/Results/temp/Fastafiles/$INPUT -db $dbname -evalue $evalue_float -num_iterations $num_itera -outfmt "$format" -out $OUTPUT.blast -out_pssm $OUTPUT.matrix.smp -seg yes -save_pssm_after_last_round -num_threads $threads)
eval $code
done

cd ${currpath}
