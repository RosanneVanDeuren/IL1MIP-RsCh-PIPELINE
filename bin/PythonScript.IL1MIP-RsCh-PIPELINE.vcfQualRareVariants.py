#!/usr/bin/env python

##Python script to filter rare variants on QUAL >= 1000 and extract those positions for samtools mpileup - Nextflow.

# import relevant modules.
import sys

# declare function that extracts information from vcf.
def getInfoVcf(line):
    cols = line.split()
    chromosome = cols[0]
    bpposition = cols[1]
    ref        = cols[3]
    alt        = cols[4]
    quality    = float(cols[5])
    zygosity   = cols[9].split(":")[0]
    return(chromosome, bpposition, ref, alt, quality, zygosity)
# define function that extracts information from dbSNPreferenceFile.
def getInfoDBsnp150(variant):
    cols = variant.split()
    CHR = cols[0]
    BP  = cols[1]
    REF = cols[3]
    ALT = cols[4]
    return(CHR, BP, REF, ALT)

# open dbSnp common variants file.
dbsnp150_open = open(str(sys.argv[2]), "r")

# store dbSnp150 common variants in list.
commonVarList = []
for variant in dbsnp150_open:
    if not variant.startswith("#CHROM"):
        CHR, BP, REF, ALT = getInfoDBsnp150(variant)
        dbSnpKey = "{0}{1}{2}{3}{4}{5}{6}".format(CHR, ":", BP, ":", REF, ":", ALT)
        commonVarList.append(dbSnpKey)
dbsnp150_open.close()

# open splitted vcf file.
vcf_open = open(str(sys.argv[1]), "r")
# create QUALfiltered vcf file to write to.
vcfQUALfiltered_open = open(str(sys.argv[3]), "w")
# create rareVariantLocs-file for samtools mpileup positions to write to.
rareVariantLocs_open = open(str(sys.argv[4]), "w")
# create rareVariantMpileupInfo-file for the background information we can use later.
rareVariantInfo_open = open(str(sys.argv[5]), "w")
rareVariantInfo_open.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format("chromosome", "bpposition", "ref", "alt", "zygosity"))

# process splitted vcf.
for line in vcf_open:
    if line.startswith('#'):
        vcfQUALfiltered_open.write(line)
    else:
        chromosome, bpposition, ref, alt, quality, zygosity = getInfoVcf(line)
        vcfKey = "{0}{1}{2}{3}{4}{5}{6}".format(chromosome, ":", bpposition, ":", ref, ":", alt)
	    rarevariantInfoKey = "{0}\t{1}\t{2}\t{3}\t{4}\n".format(chromosome, bpposition, ref, alt, zygosity)
        if vcfKey in commonVarList:
            vcfQUALfiltered_open.write(line)
        else:
            if quality >= 1000:
                vcfQUALfiltered_open.write(line)
                samtoolsKey = "{0}{1}{2}{3}{4}\n".format(chromosome, ":", bpposition, "-", bpposition)
                rareVariantLocs_open.write(samtoolsKey)
                rareVariantInfo_open.write(rarevariantInfoKey)

# close files.
vcf_open.close()
vcfQUALfiltered_open.close()
rareVariantLocs_open.close()
rareVariantInfo_open.close()
