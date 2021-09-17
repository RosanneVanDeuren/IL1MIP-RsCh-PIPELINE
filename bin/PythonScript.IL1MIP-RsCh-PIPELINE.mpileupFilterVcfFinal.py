#!/usr/bin/env python

##Python script to filter QUAL-filtered vcf on mpileups - Nextflow.

# import relevant modules.
import sys
import pandas as pd
import numpy as np

# define function that extracts information from concludedMpileups.
def getConcludedMpileupInfo(row):
    chromosome = row[0]
    bpposition = row[1]
    ref = row[2]
    alt = row[3]
    zygosity = row[4]
    mpileupConcluded = row[5]
    mpileupPercentage = row[6]
    mpileupZygosity = row[7]
    mpileupMarker = row[8]
    mpileupNote = row[9]
    return chromosome, bpposition, ref, alt, zygosity, mpileupConcluded, mpileupPercentage, mpileupZygosity, mpileupMarker, mpileupNote
# define function that extracts information from QUALfilteredVcf.
def getInfoQUALfilteredVcf(line):
    cols = line.split()
    chromosome = cols[0]
    bpposition = cols[1]
    ref        = cols[3]
    alt        = cols[4]
    return chromosome, bpposition, ref, alt

# open concludedMpileups and store in dictionary.
concludedMpileups_dict = {}
concludedMpileups_df = pd.read_csv(str(sys.argv[1]), sep = "\t", header = 0)
for index, row in concludedMpileups_df.iterrows():
    chromosome, bpposition, ref, alt, zygosity, mpileupConcluded, mpileupPercentage, mpileupZygosity, mpileupMarker, mpileupNote = getConcludedMpileupInfo(row)
    concludedMpileupKey = "{0}:{1}:{2}:{3}".format(chromosome, bpposition, ref, alt)
    concludedMpileups_dict[concludedMpileupKey] = mpileupConcluded

# open QUALfiltered vcf file and write out only rare variants that pass in final vcf.
QUALfilteredVcf_open = open(str(sys.argv[2]), "r")
finalVcf_open = open(str(sys.argv[3]), "w")
for line in QUALfilteredVcf_open:
    if line.startswith('#'):
        finalVcf_open.write(line)
    else:
        chromosome, bpposition, ref, alt = getInfoQUALfilteredVcf(line)
        vcfKey = "{0}:{1}:{2}:{3}".format(chromosome, bpposition, ref, alt)
        if vcfKey in concludedMpileups_dict.keys():
            conclusion = concludedMpileups_dict[vcfKey]
            if conclusion == "TP":
                finalVcf_open.write(line)
        else:
            finalVcf_open.write(line)

# close QUALfilteredVcf and finalVcf.
QUALfilteredVcf_open.close()
finalVcf_open.close()
