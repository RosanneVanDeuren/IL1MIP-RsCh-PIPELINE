#!/usr/bin/env python

##Python script to compute the average coverage depth over the entire panel and exclude sampleIDs with coverage below THRESHOLD - Nextflow.

# import relevant modules.
import sys
import pandas as pd
import numpy as np

# open parsedCoverage-file.
parsedCovDepth_df = pd.read_csv(str(sys.argv[1]), sep = "\t", header = 0)

# define sampleID and THRESHOLD.
sampleID = str(sys.argv[2])
THRESHOLD = int(sys.argv[3])

# compute average over the entire panel.
panelAverage = round(parsedCovDepth_df["CoverageDepth"].mean(), 2)

# if coverage above threshold return sampleID.
if panelAverage > THRESHOLD:
	print(sampleID)

# store average coverage over the entire panel in average-file.
panelAverageFile_open = open(str(sys.argv[4]), "w")
panelAverageFile_open.write("{0}\t{1}\n".format("sampleID", "panelAverage"))
panelAverageFile_open.write("{0}\t{1}\n".format(sampleID, panelAverage))
panelAverageFile_open.close()
