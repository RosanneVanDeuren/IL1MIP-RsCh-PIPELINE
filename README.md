# **IL1MIP-RsCh-PIPELINE**
A DNA analysis pipeline for research Molecular Inversion Probe sequencing Interleukin-1 panel data
#### **Description**
IL1MIP-RsCh-PIPELINE was designed to provide reliable common and rare variant calls from Interleukin(IL)-1 panel Molecular Inversion Probe (MIP)-sequencing data for the purpose of large-cohort analyses. In brief, this pipeline excludes poor quality reads from input bam-files, generates coverage (summary) statistics, calls variants with GATK UnifiedGenotyper, and finally excludes rare variants based on QUAL-parameter in the vcf-file and mpileup statistics.
&NewLine;

### **QUICKSTART**
Install Nextflow:
  - See installation instructions [link](https://www.nextflow.io/)

Install Singularity:
  - See installation instructions [link](https://sylabs.io/guides/3.0/user-guide/installation.html)

Download container images in `./containers/` folder:
  - bcftools 1.9: `singularity pull library://weizhu365/mocca-sv/bcftools_1-9:1.0.0`
  - bedtools 2.29: `singularity pull library://marialitovchenko/default/bedtools`
  - GATK 3.8.1: `singularity pull docker://broadinstitute/gatk3:3.8-1`
  - Perl 5.18.4: `singularity pull docker://perl:5.18.4-threaded-stretch`
  - Python 2.7 (includes modules pysam, numpy and pandas): ``
  - Samtools 1.11: `singularity pull library://daanjg98/rnaseq/samtools:1.11`
  
Specify required parameters in [nextflow.config](nextflow.config):
  - mainDir: Absolute path on your system to where you have installed IL1MIP-RsCh-PIPELINE
  - runID: Name of the run/cohort you wish to process
  - input: Relative path from `mainDir` to input bam-files.

Run IL1MIP-RsCh-PIPELINE with the following command:
  ```
  nextflow run main.nf -c nextflow.config
  ```

### **RUNNING**
First, make sure all components are installed properly (see QUICKSTART section).
###### **Input**


### **REQUIREMENTS**
- Nextflow 20.10.0 (or later) [link] (https://www.nextflow.io/)
- Java 8
- Singularity 3.5.1 (or later) [link] (https://sylabs.io/guides/3.0/user-guide/installation.html)
