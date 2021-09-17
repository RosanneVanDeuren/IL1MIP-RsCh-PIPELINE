#!/bin/bash

##Bashscript to perform samtools mpileup on rare variant locations - Nextflow.

locations=$(<"$3")
for pos in $locations
do
    snp=$pos
    echo $snp
    samtools mpileup -r $snp -d 10000 -q 30 -Q 20 -B -f $4 $2 > $snp.$1.mpileup.txt
done

cat *.$1.mpileup.txt > "$5"

rm chr*

echo All done.
