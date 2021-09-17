#!/usr/bin/env python

##Python script that computes summary statistics of parsed coverage depth - Nextflow.

# import relevant modules.
import sys
import pandas as pd
import numpy as np

# list with all genes as specified in MIP-probes.
geneList = ["AIM2_", "ATG16L1_", "ATG5_", "ATG7_", "BECN1_", "CASP1_", "CYBA_", "CYBB_",
           "IL18_", "IL18BP_", "IL18R1_", "IL18RAP_", "IL1A_", "IL1B_", "IL1F10_", "IL1R1_",
           "IL1R2_", "IL1RAP_", "IL1RAPL1_", "IL1RAPL2_", "IL1RL1_", "IL1RL2_", "IL1RN_", "IL33_",
           "IL36A_", "IL36B_", "IL36G_", "IL36RN_", "IL37_", "IRAK4_", "KIAA0226_", "MAP1LC3A_",
           "MYD88_", "NCF1_", "NCF2_", "NCF4_", "NLRP1_", "NLRP3_", "NOX4_", "PIK3C3_",
           "PYCARD_", "RAC1_", "RAC2_", "RAP1A_", "SIGIRR_", "TIMD4_", "TRAF6_", "UVRAG_"]

# open parsed coverage depth.
parsedCovDepth_df = pd.read_csv(str(sys.argv[1]), sep = "\t", header = 0)

# first compute panel summary statistics.
panelMean = round(parsedCovDepth_df['CoverageDepth'].mean(), 2)
panelMin = round(parsedCovDepth_df['CoverageDepth'].min(), 2)
panelMax = round(parsedCovDepth_df['CoverageDepth'].max(), 2)
panelMedian = round(parsedCovDepth_df['CoverageDepth'].median(), 2)

# create empty dataframe in which summary statistics will be stored.
sumstats_df = pd.DataFrame(columns=('sampleID', 'panel.mean', 'panel.min', 'panel.max', 'panel.median', 'chrX.genes.mean', 'chrX.genes.min', 'chrX.genes.max', 'chrX.genes.median', 'autosomal.genes.mean', 'autosomal.genes.min', 'autosomal.genes.max', 'autosomal.genes.median'))

# create lists to store summary statistics of chrX.genes and autosomal.genes.
chrX_meansList = []
chrX_minsList = []
chrX_maxsList = []
chrX_mediansList = []
autosomal_meansList = []
autosomal_minsList = []
autosomal_maxsList = []
autosomal_mediansList = []

# compute summary statistics per gene and store in empty dataframe.
ix = 0
for gene in geneList:
    subset_df = parsedCovDepth_df[parsedCovDepth_df['MIPs'].str.contains(gene)]
    meanVal = round(subset_df['CoverageDepth'].mean(), 2)
    minVal = round(subset_df['CoverageDepth'].min(), 2)
    maxVal = round(subset_df['CoverageDepth'].max(), 2)
    medianVal = round(subset_df['CoverageDepth'].median(), 2)
    sumstats_df.loc[ix, 'sampleID'] = str(sys.argv[2])
    sumstats_df.loc[ix, 'panel.mean'] = panelMean
    sumstats_df.loc[ix, 'panel.min'] = panelMin
    sumstats_df.loc[ix, 'panel.max'] = panelMax
    sumstats_df.loc[ix, 'panel.median'] = panelMedian
    sumstats_df.loc[ix, (gene[0:(len(gene)-1)] + '.mean')] = meanVal
    sumstats_df.loc[ix, (gene[0:(len(gene)-1)] + '.min')] = minVal
    sumstats_df.loc[ix, (gene[0:(len(gene)-1)] + '.max')] = maxVal
    sumstats_df.loc[ix, (gene[0:(len(gene)-1)] + '.median')] = medianVal
    #Check if chrX or autosomal gene and add to corresponding list.
    if ((gene == "CYBB_") or (gene == "IL1RAPL1_") or (gene == "IL1RAPL2_")):
        chrX_meansList.append(meanVal)
        chrX_minsList.append(minVal)
        chrX_maxsList.append(maxVal)
        chrX_mediansList.append(medianVal)
    else:
        autosomal_meansList.append(meanVal)
        autosomal_minsList.append(minVal)
        autosomal_maxsList.append(maxVal)
        autosomal_mediansList.append(medianVal)

# finally, add the summary statistics of chrX.genes and autosomal.genes.
sumstats_df.loc[ix, 'chrX.genes.mean'] = round(pd.Series(chrX_meansList).mean(), 2)
sumstats_df.loc[ix, 'chrX.genes.min'] = round(pd.Series(chrX_minsList).min(), 2)
sumstats_df.loc[ix, 'chrX.genes.max'] = round(pd.Series(chrX_maxsList).max(), 2)
sumstats_df.loc[ix, 'chrX.genes.median'] = round(pd.Series(chrX_mediansList).median(), 2)
sumstats_df.loc[ix, 'autosomal.genes.mean'] = round(pd.Series(autosomal_meansList).mean(), 2)
sumstats_df.loc[ix, 'autosomal.genes.min'] = round(pd.Series(autosomal_minsList).min(), 2)
sumstats_df.loc[ix, 'autosomal.genes.max'] = round(pd.Series(autosomal_maxsList).max(), 2)
sumstats_df.loc[ix, 'autosomal.genes.median'] = round(pd.Series(autosomal_mediansList).median(), 2)

# write summary statistics to file.
sumstats_df.to_csv(str(sys.argv[3]), sep = "\t", header = True, index = False)
