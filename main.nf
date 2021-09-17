#!/usr/bin/env nextflow

/*
 * Author: Rosanne C. van Deuren
 *
 *   This is the Nextflow main to run 'IL1MIP-RsCh-PIPELINE'.
 *
 *   IL1MIP-RsCh-PIPELINE is a pipeline designed to provide reliable common and 
 *   rare variant calls from Interleukin(IL)-1 Molecular Inversion Probe (MIP) 
 *   research data for the purpose of conducting cohort based analyses.
 *
 *   For more information or questions please contact: rc.vandeuren@outlook.com
 */

// Message.
println """
Running::

   ==    ==       ==   ==    ==  ==  ====        ====             ====             ====    ==  ====    ======  ==      ==  ==   ==  ======
   ||    ||     //||   ||\\\\//||  ||  ||  \\\\      ||  \\\\    **   //       *         ||  \\\\  ||  ||  \\\\  ||      ||      ||  |\\  ||  ||
   ||    ||       ||   ||    ||  ||  ||  //  --  ||  //   *    ||        *     --  ||  //  ||  ||  //  ||__    ||      ||  ||\\ ||  ||__
   ||    ||       ||   ||    ||  ||  ||==    --  ||==      **  ||        ***   --  ||==    ||  ||==    ||      ||      ||  || \\||  ||
   ||    ||       ||   ||    ||  ||  ||          ||  \\\\      *  \\\\       *  *      ||      ||  ||      ||      ||      ||  ||  \\|  ||
   ==    ======   ==   ==    ==  ==  ==          ==   ==   **     ====   *  *      ==      ==  ==      ======  ======  ==  ==   ==  ======

@author: Rosanne C. van Deuren
"""

// Input raw bamFiles.
inputBams_ch = Channel
    .fromPath( params.input + "*.bam")
    .map{ file -> tuple(file.simpleName, file) }

/*
 * -I.bamReadFilter-
 *
 *   Process bamReadFilter runs PythonScript.IL1MIP-RsCh-PIPELINE.bamReadFilter.py
 *   on raw bam-files containing IL-1 panel MIP research data aligned by BWA-MEM.
 *
 *   The bam-files are filtered on the following parameters:
 *     1. mapping quality >= 60
 *     2. no soft-clipped reads
 *     3. only properly paired reads as defined by the read flag
 *     4. discard reads with five or more variations from the reference
 *        with exception of large-base INDELs.
 */

// Run process bamReadFilter.
process bamReadFilter {

    publishDir params.outfilePath + params.publishDir_bamFolder, mode: params.publishDir_mode, overwrite: params.publishDir_overwrite
    tag "${sampleName}"

    input:
    tuple val(sampleName), file(bamFile) from inputBams_ch

    output:
    set val(sampleName), file("${sampleName}.sorted.hfix.rf.bam"), file("${sampleName}.sorted.hfix.rf.bam.bai") into filteredBams_ch

    """
    #!/bin/bash

    python ${params.scripts}${params.script_bamReadFilter} ${bamFile} "${sampleName}.sorted.hfix.rf.bam"
    """

}

// Duplicate filteredBams_ch for 1) bedtoolsCoverage, 2) variantCalling, 3) rareVariantMpileup.
filteredBams_ch.into{ filteredBams_bedtoolsCoverage_ch; filteredBams_variantCalling_ch; filteredBams_rareVariantMpileup_ch }

/*
 * -II.bedtoolsCoverage-
 * 
 *   Process bedtoolsCoverage computes the remaining coverage depth on the
 *   read-filtered bam-files.
 */

// Run process bedtoolsCoverage.
process bedtoolsCoverage {

    tag "${sampleName}"

    input:
    set val(sampleName), file(filteredBam), file(filteredBai) from filteredBams_bedtoolsCoverage_ch

    output:
    tuple val(sampleName), file("${sampleName}.rfBam.covDepth.txt") into rawCoverages_ch

    """
    #!/bin/bash

    bedtools coverage -a ${params.targetCallingfileInclPath} -b ${filteredBam} -d > "${filteredBam.simpleName}.rfBam.covDepth.txt"
    """

}

/*
 * -III.parseCovDepth-
 * 
 *   Process parseCovDepth runs PythonScript.IL1MIP-RsCh-PIPELINE.parseCovDepth.py
 *   which parses and sorts the raw coverage.
 */

// Run process parseCovDepth.
process parseCovDepth {

    tag "${sampleName}"

    input:
    tuple val(sampleName), file(rawCoverage) from rawCoverages_ch

    output:
    tuple val(sampleName), file("${sampleName}.rfBam.covDepth.parsed.txt") into parsedCoverages_ch

    """
    #!/bin/bash

    python ${params.scripts}${params.script_parseCovDepth} ${rawCoverage} "${sampleName}.tmpparsed.txt" "${sampleName}.rfBam.covDepth.parsed.txt"
    """

}

// Duplicate parsedCoverages_ch for 1) coverageFilter, 2) coverageSumstats.
parsedCoverages_ch.into{ parsedCoverages_coverageFilter_ch; parsedCoverages_coverageSumstats_ch }

/*
 * -IV.coverageFilter-
 *
 *   Process coverageFilter runs PythonScript.IL1MIP-RsCh-PIPELINE.coverageFilter.py
 *   in which the average coverage per sample over the entire panel is computed,
 *   after which this is combined for all samples and published in a file.
 *
 *   At the same time, only sampleIDs that pass the coverage threshold specified in 
 *   the nextflow.config parameter "coverageThreshold" are passed on to the next process.
 */

// Run process coverageFilter.
process coverageFilter {

    tag "${sampleName}"

    input:
    tuple val(sampleName), file(parsedCoverage) from parsedCoverages_coverageFilter_ch

    output:
    stdout into coveragePass_ch
    file "${sampleName}.rfBam.covDepth.parsed.average.txt" into averageCoverages

    """
    #!/bin/bash

    python ${params.scripts}${params.script_coverageFilter} ${parsedCoverage} ${sampleName} ${params.coverageThreshold} "${sampleName}.rfBam.covDepth.parsed.average.txt"
    """

}

// Store coveragePass in coveragePass_forVariantCalling_ch.
coveragePass_forVariantCalling_ch = coveragePass_ch
    .map { [it.replaceAll("\n", "")] }
    .flatten()

// Save average coverages in a file.
averageCoverages
    .collectFile(name: "IL1MIP-RsCh-PIPELINE_averagePanelCoverage_" + params.runID + "_" + params.runDATE + ".txt", newLine: params.collectAverages_newLine, keepHeader: params.collectAverage_keepHeader, storeDir: params.outfilePath + params.publishDir_covFolder )
    .subscribe{
        println "Average coverages are saved to file: $it"
    }

/*
 * -V.coverageSumstats-
 * 
 *   Process coverageSumstats runs PythonScript.IL1MIP-RsCh-PIPELINE.coverageSumstats.py
 *   in which summary statistics of the parsed and sorted coverage are computed for all samples
 *   and published into a combined file. 
 */

// Run process coverageSumstats.
process coverageSumstats {

    tag "$sampleName"

    input:
    tuple val(sampleName), file(parsedCoverage) from parsedCoverages_coverageSumstats_ch

    output:
    file "${sampleName}.rfBam.covDepth.parsed.sumstats.txt" into sumstatsCoverages

    """
    #!/bin/bash

    python ${params.scripts}${params.script_coverageSumstats} ${parsedCoverage} $sampleName "${sampleName}.rfBam.covDepth.parsed.sumstats.txt"
    """

}

// Save coverage summary statistics of all individuals in the cohort in a file.
sumstatsCoverages
    .collectFile(name: "IL1MIP-RsCh-PIPELINE_sumstatsCoverage_" + params.runID + ".txt", newLine: false, keepHeader: true, storeDir: params.outfilePath + params.publishDir_covFolder )
    .subscribe{
        println "Coverage summary statistics are saved to file: $it"
    }

/*
 * -VI.variantCalling-
 * 
 *   Process variantCalling calls variants from the read-filtered bam-files
 *   for sampleIDs that passed the coverageThreshold using GATK UnifiedGenotyper.
 *   This process uses "targetCallingFile_il1mips.bed" restricting variant calls to
 *   the regions specified in the bed-file, and a fasta reference of Hg19. These files,
 *   along with GATK variant calling parameters are specified in the nextflow.config.
 */

// Combine channels for process variantCalling.
inclBams_forVariantCalling_ch = coveragePass_forVariantCalling_ch.join(filteredBams_variantCalling_ch)

// Run process variantCalling.
process variantCalling {

    tag "${sampleName}"

    input:
    set val(sampleName), file(filteredBam), file(filteredBai) from inclBams_forVariantCalling_ch

    output:
    tuple val(sampleName), file("${sampleName}.rfBam.UnifiedGenotyperVariants.vcf") into rawVcfs_ch

    """
    #!/bin/bash

    java -jar /usr/GenomeAnalysisTK.jar -T UnifiedGenotyper -R ${params.fastaHg19fileInclPath} -I ${filteredBam} -o "${sampleName}.rfBam.UnifiedGenotyperVariants.vcf" -L ${params.targetCallingfileInclPath} -dcov ${params.gatk_dcov} -dt ${params.gatk_dt} -rf ${params.gatk_rf} -glm ${params.gatk_glm} -stand_call_conf ${params.gatk_stand_call_conf}
    """

}

/*
 * -VII.splitMultiAllelicSites-
 * 
 *   Process splitMultiAllelicSites splits multi-allelic variant sites from raw vcf-files
 *   using bcftools.
 */

// Run process splitMultiAllelicSites.
process splitMultiAllelicSites {

    tag "${sampleName}"

    input:
    tuple val(sampleName), file(rawVcf) from rawVcfs_ch

    output:
    tuple val(sampleName), file("${sampleName}.rfBam.UnifiedGenotyperVariants.splitted.vcf") into splittedVcfs_ch

    """
    #!/bin/bash

    bcftools norm -m - -f ${params.fastaHg19fileInclPath} ${rawVcf} -o "${sampleName}.rfBam.UnifiedGenotyperVariants.splitted.vcf"
    """

}

/*
 * -VIII.vcfQUALrarevariants-
 * 
 *   Process vcfQUALrarevariants runs PythonScript.IL1MIP-RsCh-PIPELINE.vcfQUALrarevariants.py
 *   which filters out rare variants with a QUAL-value <1000, and saves the rare variant positions
 *   that pass for samtools mpileup.
 *   Rare variants are defined as being absent in "dbSnp150common_il1mipTargetRegion.txt" (as specified
 *   in the nextflow.config), a reference file containing all common variant positions extracted from 
 *   dbSnp v150 within the IL-1 target region.
 */

// Run process vcfQUALrarevariants.
process vcfQUALrarevariants {

    tag "${sampleName}"

    input:
    tuple val(sampleName), file(splittedVcf) from splittedVcfs_ch

    output:
    tuple val(sampleName), file("${sampleName}.rfBam.UnifiedGenotyperVariants.splitted.QUALfiltered.vcf") into vcfsQUALfiltered_ch
    tuple val(sampleName), file("${sampleName}.rarevariantLocs.txt") into rareVariantLocs_ch
    tuple val(sampleName), file("${sampleName}.rarevariantInfo.txt") into rareVariantInfo_ch

    """
    #!/bin/bash

    python ${params.scripts}${params.script_vcfQUALrarevariants} ${splittedVcf} $params.dbSnp150commonInclPath "${sampleName}.rfBam.UnifiedGenotyperVariants.splitted.QUALfiltered.vcf" "${sampleName}.rarevariantLocs.txt" "${sampleName}.rarevariantInfo.txt"
    """

}

/*
 * -IX.rarevariantMpileup-
 * 
 *   Process rarevariantMpileup runs BashScript.IL1MIP-RsCh-PIPELINE.rarevariantMpileup.sh
 *   which computes mpileups for all rare variants from the QUAL-filtered vcf-file using
 *   samtools.
 */

// Combine channels for mpileup input.
filteredBams_rareVariantMpileup_ch
    .join(rareVariantLocs_ch)
    .map { it ->
        def sampleName = file(it[0]).getName();
        def filteredBam = file(it[1]);
        def filteredBai = file(it[2]);
        def locations = file(it[3]);

        [sampleName, filteredBam, filteredBai, locations]
    }
    .set { input_mpileups_ch }

// Run process rarevariantMpileup.
process rarevariantMpileup {

    tag "${sampleName}"

    input:
    set val(sampleName), file(filteredBam), file(filteredBai), file(locations) from input_mpileups_ch

    output:
    tuple val(sampleName), file("${sampleName}.rarevariant.mpileupsmerged.txt") into mergedMpileups_ch

    """
    #!/bin/bash

    bash ${params.scripts}${params.script_rarevariantMpileup} $sampleName ${filteredBam} ${locations} $params.fastaHg19fileInclPath "${sampleName}.rarevariant.mpileupsmerged.txt"
    """

}

/*
 * -X.preprocessMpileups-
 * 
 *   Process preprocessMpileups runs PythonScript.IL1MIP-RsCh-PIPELINE.preprocessMpileups.py
 *   which prepares the mpileups for parsing with perl library. 
 * 
 *   Preprocessing consists of removing strand-specific information (e.g. unifying forward and
 *   reverse alleles), and saving information regarding readMarkers.
 */

// Run process preprocessMpileups.
process preprocessMpileups {

    tag "${sampleName}"

    input:
    tuple val(sampleName), file(mpileup) from mergedMpileups_ch

    output:
    tuple val(sampleName), file("${sampleName}.rarevariant.mpileupsmerged.prepped.txt") into preppedMpileups_ch
    tuple val(sampleName), file("${sampleName}.rarevariant.mpileupsmerged.readMarkerInfo.txt") into readMarkerInfo_ch

    """
    #!/bin/bash

    python ${params.scripts}${params.script_preprocessMpileups} ${mpileup} "${sampleName}.rarevariant.mpileupsmerged.prepped.txt" "${sampleName}.rarevariant.mpileupsmerged.readMarkerInfo.txt"
    """

}

/*
 * -XI.parseMpileups-
 * 
 *   Process parseMpileups runs PerlScript.IL1MIP-RsCh-PIPELINE.pileup2baseindel.no.strand.pl
 *   which is a script from library pileup2base (https://github.com/riverlee/pileup2base.git)
 *   that parses the preprocessed mpileups so that variant allele frequencies can be
 *   calculated.
 */

// Run process parseMpileups.
process parseMpileups {

    tag "${sampleName}"

    input:
    tuple val(sampleName), file(preppedMpileup) from preppedMpileups_ch

    output:
    tuple val(sampleName), file("${preppedMpileup.simpleName}.rarevariant.mpileupsmerged.parsed.txt") into parsedMpileups_ch

    """
    #!/bin/bash

    perl ${params.scripts}${params.script_parseMpileups} "${preppedMpileup}" 1 "${sampleName}.rarevariant.mpileupsmerged.parsed.txt"
    """

}

/*
 * -XII.concludeMpileups-
 *
 *   Process concludeMpileups runs PythonScript.IL1MIP-RsCh-PIPELINE.concludeMpileups.py
 *   which uses the parsed mpileups to calculate the alternative allele frequency of rare
 *   variants, and additionally saves information regarding read-markers and zygosity 
 *   inconsistencies.
 *
 *   The following criteria define whether a rare variant is considered false positive:
 *     -Rare variant positions with a total read count < 20
 *     -Rare variant positions with a VAF < 25%
 *   This information is stored in the "mpileupConcluded" column in the published output-file.
 * 
 *   The following criteria define whether a rare variant is considered heterozygous/homozygous:
 *     -Rare variant positions with a 25% < VAF <= 90% are considered heterozygous
 *     -Rare variant positions with a VAF > 90% are considered homozygous
 *   This information is stored in the "zygosity" column in the published output-file.
 */

// Combine input channels.
parsedMpileups_ch
    .join(readMarkerInfo_ch)
    .join(rareVariantInfo_ch)
    .map { it ->
        def sampleName = file(it[0]).getName();
        def parsedMpileup = file(it[1]);
        def readMarkerInfo = file(it[2]);
        def rarevariantInfo = file(it[3]);

        [sampleName, parsedMpileup, readMarkerInfo, rarevariantInfo]
    }
    .set { parsedMpileups_readMarkerInfo_rarevariantInfo_ch }

// Run process concludeMpileups.
process concludeMpileups {

    publishDir params.outfilePath + params.publishDir_vcfFolder, mode: params.publishDir_mode, overwrite: params.publishDir_overwrite
    tag "${sampleName}"

    input:
    set val(sampleName), file(parsedMpileup), file(readMarkerInfo), file(rarevariantInfo) from parsedMpileups_readMarkerInfo_rarevariantInfo_ch

    output:
    tuple val(sampleName), file("${sampleName}.rarevariant.mpileupsmerged.concluded.txt") into concludedMpileups_ch

    """
    #!/bin/bash

    python ${params.scripts}${params.script_concludeMpileups} ${parsedMpileup} ${readMarkerInfo} ${rarevariantInfo} "${sampleName}.rarevariant.mpileupsmerged.concluded.txt"
    """

}

/*
 * -XIII.mpileupFilterVcfFinal-
 * 
 *   Process mpileupFilterVcfFinal runs PythonScript.IL1MIP-RsCh-PIPELINE.mpileupFilterVcfFinal.py
 *   which excludes rare variants from the QUAL-filtered vcf-files based on the concluded
 *   mpileup information and publishes the final filtered vcf-files.
 */

// Combine input channels.
concludedMpileups_ch
    .join(vcfsQUALfiltered_ch)
    .map { it ->
        def sampleName = file(it[0]).getName();
        def mpileupConcluded = file(it[1]);
        def vcfQUALfiltered = file(it[2]);

        [sampleName, mpileupConcluded, vcfQUALfiltered]
    }
    .set { mpileupsConcluded_vcfsQUALfiltered_ch }

// Run process mpileupFilterVcfFinal.
process mpileupFilterVcfFinal {

    publishDir params.outfilePath + params.publishDir_vcfFolder, mode: params.publishDir_mode, overwrite: params.publishDir_overwrite
    tag "${sampleName}"

    input:
    set val(sampleName), file(mpileupConcluded), file(vcfQUALfiltered) from mpileupsConcluded_vcfsQUALfiltered_ch

    output:
    tuple val(sampleName), file("${sampleName}.rfBam.UnifiedGenotyperVariants.final.vcf") into finalVcfs_ch

    """
    #!/bin/bash

    python ${params.scripts}${params.script_mpileupFilterVcfFinal} ${mpileupConcluded} ${vcfQUALfiltered} "${sampleName}.rfBam.UnifiedGenotyperVariants.final.vcf"
    """

}

