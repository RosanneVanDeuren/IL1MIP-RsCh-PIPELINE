#!/usr/bin/env python

##Python script to prefilter mpileups before parsing by perl - Nextflow.

# import relevant modules.
import sys

# define function that extracts information from mpileup file.
def getinfompileup(line):
    colmns = line.strip('\n').split('\t')
    pileup_chr = colmns[0]
    pileup_bp = int(colmns[1])
    pileup_ref = colmns[2]
    pileup_rc = int(colmns[3])
    pileup_seq = (colmns[4]).upper().replace(',', '.')
    pileup_bq = colmns[5]
    return(pileup_chr, pileup_bp, pileup_ref, pileup_rc, pileup_seq, pileup_bq)
# define function that writes to readmarker files.
def writeReadMarker(pileup_seq):
    if '^' in pileup_seq:
        readMarker_open.write('{0}\t{1}\t{2}\t{3}\n'.format(pileup_chr, pileup_bp, pileup_ref, "readStartMarker"))
    elif '$' in pileup_seq:
        readMarker_open.write('{0}\t{1}\t{2}\t{3}\n'.format(pileup_chr, pileup_bp, pileup_ref, "readEndMarker"))

# open mpileups.merged file.
pileup_open = open(str(sys.argv[1]), "r")
# create mpileups.merged.prepped to write to file.
pileupPrepped_open = open(str(sys.argv[2]), "w")
# create readmarker file to write to.
readMarker_open = open(str(sys.argv[3]), "w")
readMarker_open.write("{0}\t{1}\t{2}\t{3}\n".format("chromosome", "bpposition", "ref", "mpileupSeqMarker"))

# process mpileups.merged.
for line in pileup_open:
    pileup_chr, pileup_bp, pileup_ref, pileup_rc, pileup_seq, pileup_bq = getinfompileup(line)
    if '^' in pileup_seq or '$' in pileup_seq:
        writeReadMarker(pileup_seq)
        num = -1
        checkSeq = pileup_seq + 'X'
        NEWseq   = ""
        for base in checkSeq:
            num +=1
            if num+1 >= len(checkSeq):
                break
            elif checkSeq[num] == '^':
                NEWseq = NEWseq + '.'
                num +=2
            elif checkSeq[num+1] == '$':
                NEWseq = NEWseq + '.'
                num +=1
            else:
                NEWseq = NEWseq + checkSeq[num]
        pileupPrepped_open.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n'.format(pileup_chr, pileup_bp, pileup_ref, pileup_rc, NEWseq, pileup_bq))
    else:
        pileupPrepped_open.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n'.format(pileup_chr, pileup_bp, pileup_ref, pileup_rc, pileup_seq, pileup_bq))

# close files.        
pileupPrepped_open.close()
pileup_open.close()
readMarker_open.close()
