#!/usr/bin/env python

##Python script that filters reads on mappingQuality >= 60 && nsc && pprsF && pNMlt5 - Nextflow.

# import relevant modules.
import sys
import pysam

# define infile and outfile.
infile = pysam.AlignmentFile(str(sys.argv[1]), "rb")
outfile = pysam.AlignmentFile(str(sys.argv[2]), "wb", template=infile)

# process the reads.
for read in infile.fetch(until_eof=True, multiple_iterators=True):
    # filter on mappingQuality >= 60.
    if read.mapping_quality >= 60:
        # filter on softClipping.
        if not "S" in read.cigarstring:
            # filter on properlyPairedFlag.
            if read.flag & 2:
                # filter on NMlt5 or/else INDEL>4.
                if read.get_tag('NM') < 5:
                    outfile.write(read)
                else:
                    CIGARtuple = read.cigartuples
                    for tupl in range(0, len(CIGARtuple)):
                        if (CIGARtuple[tupl][0] == 1) or (CIGARtuple[tupl][0] == 2):
                            if (CIGARtuple[tupl][1] > 4):
                                outfile.write(read)

# close infile and outfile.
infile.close()
outfile.close()

# index filtered bamFile.
pysam.index(str(sys.argv[2]), str(sys.argv[2]) + ".bai")
