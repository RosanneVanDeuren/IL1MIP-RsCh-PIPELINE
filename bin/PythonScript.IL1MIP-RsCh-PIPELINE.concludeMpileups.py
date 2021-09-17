#!/usr/bin/env python

##Python script to conclude mpileups - Nextflow.

# import relevant modules.
import sys
import pandas as pd
import numpy as np

# define function that extracts readMarkerInfo.
def getReadMarkerInfo(row):
    chromosome = row[0]
    bpPosition = row[1]
    reference = row[2]
    mpileupSeqMarker = row[3]
    return chromosome, bpPosition, reference, mpileupSeqMarker
# define function that extracts rarevariantInfo.
def getRarevariantInfo(row):
    chromosome = row[0]
    bpPosition = row[1]
    reference = row[2]
    alternative = row[3]
    zygosity = row[4]
    return chromosome, bpPosition, reference, alternative, zygosity
def getParsedMpileupInfo(row):
    chromosome = row[0]
    bpPosition = row[1]
    reference = row[2]
    alleleA = row[3]
    alleleT = row[4]
    alleleC = row[5]
    alleleG = row[6]
    insertion = row[7]
    deletion = row[8]
    return chromosome, bpPosition, reference, alleleA, alleleT, alleleC, alleleG, insertion, deletion
# define function that writes mpileupConcluded and mpileupNote.
def writeMpileupConclusionNote(altPercentage, zygosity, readmarker, rarevariantInfo_df, index):
    if (altPercentage < 25):
        rarevariantInfo_df.loc[index, "mpileupConcluded"] = "FP"
        rarevariantInfo_df.loc[index, "mpileupPercentage"] = str(altPercentage)
        rarevariantInfo_df.loc[index, "mpileupZygosity"] = "NA"
        rarevariantInfo_df.loc[index, "mpileupMarker"] = readmarker
        rarevariantInfo_df.loc[index, "mpileupNote"] = "LowAlternativePercentage"
    elif (25 <= altPercentage < 90):
        if ((zygosity == "0/1") or (zygosity == "1/0")):
            rarevariantInfo_df.loc[index, "mpileupConcluded"] = "TP"
            rarevariantInfo_df.loc[index, "mpileupPercentage"] = str(altPercentage)
            rarevariantInfo_df.loc[index, "mpileupZygosity"] = "OK"
            rarevariantInfo_df.loc[index, "mpileupMarker"] = readmarker
            rarevariantInfo_df.loc[index, "mpileupNote"] = "NA"
        else:
            rarevariantInfo_df.loc[index, "mpileupConcluded"] = "TP"
            rarevariantInfo_df.loc[index, "mpileupPercentage"] = str(altPercentage)
            rarevariantInfo_df.loc[index, "mpileupZygosity"] = "Uncertain"
            rarevariantInfo_df.loc[index, "mpileupMarker"] = readmarker
            rarevariantInfo_df.loc[index, "mpileupNote"] = "NA"
    else:
        if (zygosity == "1/1"):
            rarevariantInfo_df.loc[index, "mpileupConcluded"] = "TP"
            rarevariantInfo_df.loc[index, "mpileupPercentage"] = str(altPercentage)
            rarevariantInfo_df.loc[index, "mpileupZygosity"] = "OK"
            rarevariantInfo_df.loc[index, "mpileupMarker"] = readmarker
            rarevariantInfo_df.loc[index, "mpileupNote"] = "NA"
        else:
            rarevariantInfo_df.loc[index, "mpileupConcluded"] = "TP"
            rarevariantInfo_df.loc[index, "mpileupNote"] = str(altPercentage)
            rarevariantInfo_df.loc[index, "mpileupZygosity"] = "Uncertain"
            rarevariantInfo_df.loc[index, "mpileupMarker"] = readmarker
            rarevariantInfo_df.loc[index, "mpileupNote"] = "NA"
    return rarevariantInfo_df

# open parsedMpileups-file, calculate readCount and store information in dictionary.
mpileupVariantIndex_dict = {}
parsedMpileups_df = pd.read_csv(str(sys.argv[1]), sep = "\t", header = 0)
parsedMpileups_df["readCount"] = parsedMpileups_df["A"] + parsedMpileups_df["C"] + parsedMpileups_df["G"] + parsedMpileups_df["T"]
for index, row in parsedMpileups_df.iterrows():
    chromosome, bpPosition, reference, alleleA, alleleT, alleleC, alleleG, insertion, deletion = getParsedMpileupInfo(row)
    mpileupVariantIndexKey = "{0}:{1}:{2}".format(chromosome, bpPosition, reference)
    mpileupVariantIndex_dict[mpileupVariantIndexKey] = index

# open readMarkerInfo-file and store information in dictionary.
readMarkerInfo_dict = {}
readMarkerInfo_df = pd.read_csv(str(sys.argv[2]), sep = "\t", header = 0)
for index, row in readMarkerInfo_df.iterrows():
    chromosome, bpPosition, reference, mpileupSeqMarker = getReadMarkerInfo(row)
    readMarkerInfoKey = "{0}:{1}:{2}".format(chromosome, bpPosition, reference)
    readMarkerInfo_dict[readMarkerInfoKey] = mpileupSeqMarker

# open rarevariantInfo-file and extract corresponding mpileup-percentage.
rarevariantInfo_df = pd.read_csv(str(sys.argv[3]), sep = "\t", header = 0)
for index, row in rarevariantInfo_df.iterrows():
    chromosome, bpPosition, reference, alternative, zygosity = getRarevariantInfo(row)
    checkKey = "{0}:{1}:{2}".format(chromosome, bpPosition, reference[0])
    checkMpileupIndex = mpileupVariantIndex_dict[checkKey]
    # check for readMarker.
    if checkKey in readMarkerInfo_dict.keys():
        readmarker = readMarkerInfo_dict[checkKey]
    else:
        readmarker = "NA"
    # check for read count.
    if (parsedMpileups_df.loc[checkMpileupIndex, "readCount"] <= 20):
        rarevariantInfo_df.loc[index, "mpileupConcluded"] = "FP"
        rarevariantInfo_df.loc[index, "mpileupPercentage"] = "NA"
        rarevariantInfo_df.loc[index, "mpileupZygosity"] = "NA"
        rarevariantInfo_df.loc[index, "mpileupMarker"] = readmarker
        rarevariantInfo_df.loc[index, "mpileupNote"] = "LowReadCount"
    else:
        # substitutions.
        if (len(reference) == 1) and (len(alternative) == 1):
            altToCheck = alternative
            altPercentage  = round((float(parsedMpileups_df.loc[checkMpileupIndex, altToCheck]) / parsedMpileups_df.loc[checkMpileupIndex, "readCount"]) * 100, 2)
            rarevariantInfo_df = writeMpileupConclusionNote(altPercentage, zygosity, readmarker, rarevariantInfo_df, index)
        # deletions.
        elif (len(reference) > 1):
            altToCheck = reference[1:]
            if pd.isnull(parsedMpileups_df.loc[checkMpileupIndex, 'Deletion']):
                rarevariantInfo_df.loc[index, "mpileupConcluded"] = "FP"
                rarevariantInfo_df.loc[index, "mpileupPercentage"] = "NA"
                rarevariantInfo_df.loc[index, "mpileupZygosity"] = "NA"
                rarevariantInfo_df.loc[index, "mpileupMarker"] = readmarker
                rarevariantInfo_df.loc[index, "mpileupNote"] = "NoDeletions"
            else:
                deletions = parsedMpileups_df.loc[checkMpileupIndex, 'Deletion'].split('|')
                for option in deletions:
                    alleleCount = option.split(":")[0]
                    alleleDeletion = option.split(":")[1]
                    if alleleDeletion == altToCheck:
                        altPercentage  = round((float(alleleCount) / parsedMpileups_df.loc[checkMpileupIndex, "readCount"]) * 100, 2)
                        rarevariantInfo_df = writeMpileupConclusionNote(altPercentage, zygosity, readmarker, rarevariantInfo_df, index)
        # insertions.
        elif (len(alternative) > 1):
            altToCheck = alternative[1:]
            if pd.isnull(parsedMpileups_df.loc[checkMpileupIndex, 'Insertion']):
                rarevariantInfo_df.loc[index, "mpileupConcluded"] = "FP"
                rarevariantInfo_df.loc[index, "mpileupPercentage"] = "NA"
                rarevariantInfo_df.loc[index, "mpileupZygosity"] = "NA"
                rarevariantInfo_df.loc[index, "mpileupMarker"] = readmarker
                rarevariantInfo_df.loc[index, "mpileupNote"] = "NoInsertions"
            else:
                insertions = parsedMpileups_df.loc[checkMpileupIndex, 'Insertion'].split('|')
                for option in insertions:
                    alleleCount = option.split(":")[0]
                    alleleInsertion = option.split(":")[1]
                    if alleleInsertion == altToCheck:
                        altPercentage  = round((float(alleleCount) / parsedMpileups_df.loc[checkMpileupIndex, "readCount"]) * 100, 2)
                        rarevariantInfo_df = writeMpileupConclusionNote(altPercentage, zygosity, readmarker, rarevariantInfo_df, index)

# write mpileups concluded to file.
rarevariantInfo_df.to_csv(str(sys.argv[4]), sep = "\t", header = True, index = False)
