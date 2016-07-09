NucProcess
----------

NucProcess is a Python program to perform single-cell Hi-C sequence processing
of paired read data. This software takes paired FASTQ sequence read files
(representing HiC data for only ONE cell), a reference genome sequence and
knowledge of experimental parameters (restriction enzyme type and the range of
size of DNA fragments sequenced in the library) to create processed, single-cell
Hi-C contact files. By default, the ouptput is generated in Nuc Chromoatin
Contact (NCC) format, which is a simple text based format described below. A
processing report document will also be generated for each run; in SVG format
that can be viewed in most web browsers.

To run NucProcess issue the 'nuc_process' command with the command line options
described below. The options -i (input FASTQ files) and -g (genome reference)
are mandatory, though its is usual to also use -re1 (primary restriction enzyme.
Default is MboI), -o (root name of output files) and -re2 (secondary restriction
enzyme in double-digest experiments).

IMPORTANT: If the files for the genome index and the corresponding restriction
enzyme (RE) cut location files have not been created they will be generated
automatically by NucProcess so long as the -f option is specfied. This option
specifies the location of complete chromosomal sequences in FASTA format,
typically by using a wildcard specification. Once these files are present the -f
option need not be specified, but it will not trigger the re-creation of the
files unless the -m or -x options are used (or if files are deleted).

NOTE: the names of the chromosomes used by NucProcess are determined by the file
names of the chromosome sequence files that were used to build the genome index
and RE files (and these must match). The chromosome sequence files should be
tagged with the chromosome name after an underscore but before the file
extension. For example "chr7" will be taken from a file named
"mm_ref_GRCm38.p2_chr7.fa".



Barcoded input
--------------

The splitFastqBarcodes.py script is provided to split FASTQ files that represent
many cells, each with a different barcode sequence, into separate paired read
files. The script can be run as follows, specifying the names of the two paired,
multiplex FASTQ read files after the script:

  python splitFastqBarcodes.py MULTIPLEXED_DATA_r_1.fq MULTIPLEXED_DATA_r_1.fq

This will generate paired FASTQ files of the form:

  MULTIPLEXED_DATA_r_1_CGC.fq  MULTIPLEXED_DATA_r_2_CGC.fq
  MULTIPLEXED_DATA_r_1_TAA.fq  MULTIPLEXED_DATA_r_2_TAA.fq

where the file names are tagged with the corresponding barcode sequence. These
demultiplexed files can then be used as input to NucProcess, specifying only one
barcode for each run.



Running NucProcess
------------------

Typical first-time use:

  nuc_process  -f /chromosomes/*.fa -o CELL_1 -v -a -k -re1 MboI -re2 AluI -s 150-2000 -n 12 -g /genome/GENOME_BUILD -i /data/SEQUENCING_DATA_r_?.fq

Typical use thereafter:

  nuc_process -o CELL_1 -v -a -k -re1 MboI -re2 AluI -s 150-2000 -n 8 -g /genome/GENOME_BUILD -i /data/SEQUENCING_DATA_r_?.fq


For the above commands:

-f /chromosomes/*.fa states that all FASTA files (ending in .fa) in the
/chromosomes/ directory will be used for creation of the genome index and RE cut
site files

-o spcifies CELL_1 will be used for naming the output. In this case the main
output contact file will be CELL_1.ncc

-v spcifies verbose output of processing progress

-a spcifies to generate ambigous contact files: CELL_1_ambig.ncc in this case

-k spcifies to keep all the intermediate processing files: filered NCC files,
clipped FASTQ files and the main Bowtie2 mapping SAM files 

-re1 is the primary restrition enzyme type at the ligation junction (see
enzymes.conf)

-re2 is the secondary restriction enzyme ued to release the fragments (option
not used for Nextera based protocol)

-s is the valid molecule/fragment size range, as used in the DNA sequencing

-n is the number of parallel CPU cores to use with Bowtie2

-g GENOME_BUILD is the root name for the Bowtie2 genome index without any file
extension and in this case would refer to files GENOME_BUILD.1.bt2,
GENOME_BUILD.rev.1.bt2 etc.

-i SEQUENCING_DATA_r_?.fq is a wildcard expression matching the two input FASTQ
files (though two separate file names, separated by a space can be specified).
In this case the expression would match SEQUENCING_DATA_r_1.fq and
SEQUENCING_DATA_r_2.fq - the paired sequence read files.



NCC data format
---------------

The main NCC output format for contact data takes the form of space-separate
lines, where each line represents a different chromosomal contact that pairs two
chromosomal regions. This format specifies more than just the chromosome contact
map: it also records the original read locations within the (ligated) primary RE
digest fragments, strand information, ambiguity information and information to
relate the data back to the original FASTQ input files. Accoringly this format
can be used for data filtering and validation.

The columns of NCC files correspond to: 

  Name of chromosome A
  First base position of sequence read A
  Last base position of sequence read A
  5' base position of primary RE fragment containing read A
  3' base position of primary RE fragment containing read A
  The strand of sequence read A
  Name of chromosome B
  First base position of sequence read B
  Last base position of sequence read B
  5' base position of primary RE fragment containing read B
  3' base position of primary RE fragment containing read B
  The strand of sequence read B
  The number of the ambiguity group to which the paired reads belong
  The ID number of the read pair in the original FASTQ files
  Whether read pairs are swapped relative to original FASTQ files

i.e:

  chr_a start_a end_a re1_a_start re1_a_end strand_a chr_b start_b end_b re1_b_start re1_b_end strand_b ambig_group pair_id swap_pair

For example two lines could be:

  chr10 100002828 100002899 100002733 100003107 + chr10 100015771 100015700 100015676 100016001 - 1498612 1534464 0
  chr10 100007729 100007658 100007551 100008354 - chr13 107185630 107185700 107184984 107185698 + 1009602 1032357 1



Command line options
--------------------

usage: nuc_process [-h] [-i FASTQ_FILE FASTQ_FILE] [-g GENOME_FILE]
                   [-re1 ENZYME] [-re2 ENZYME] [-s SIZE_RANGE] [-n CPU_COUNT]
                   [-r COUNT] [-o NCC_FILE] [-oa NCC_FILE] [-or REPORT_FILE]
                   [-b EXE_FILE] [-q SCHEME] [-m] [-p] [-x]
                   [-f FASTA_FILES [FASTA_FILES ...]] [-a] [-k] [-sam]
                   [-l SEQUENCE] [-z] [-v] [-u]

Chromatin contact paired-read Hi-C processing module for Nuc3D and NucTools

optional arguments:
  -h, --help            show this help message and exit
  -i FASTQ_FILE FASTQ_FILE
                        Input paired-read FASTQ files to process (accepts
                        wildcards that match two files)
  -g GENOME_FILE        Genome index file to map sequence reads to. A new
                        index will be created with this name if the index is
                        missing and genome FASTA files are specified
  -re1 ENZYME           Primary restriction enzyme (for ligation junctions).
                        Default: MboI. Available: AluI, BglII, MboI
  -re2 ENZYME           Secondary restruction enzyme (if used). Available:
                        AluI, BglII, MboI
  -s SIZE_RANGE         Allowed range of sequenced molecule sizes, e.g.
                        "150-1000", "100,800" or "200" (no maximum)
  -n CPU_COUNT          Number of CPU cores to use in parallel
  -r COUNT              Mimimum number of sequencing repeats required to
                        support a contact
  -o NCC_FILE           Optional output name for NCC format chromsome contact
                        file
  -oa NCC_FILE          Optional output name for ambiguous contact NCC file
  -or REPORT_FILE       Optional output name for HTML report file
  -b EXE_FILE           Path to bowtie2 (read aligner) executable (will be
                        searched for if not specified)
  -q SCHEME             Use a specific FASTQ quility scheme (normally not set
                        and deduced automatically). Available: integer,
                        phred33, phred64, solexa
  -m                    Force a re-mapping of genome restriction enzyme sites
                        (otherwise cached values will be used if present)
  -p                    The input data is population Hi-C; single-cell
                        processing steps are avoided
  -x, --reindex         Force a re-indexing of the genome (given appropriate
                        FASTA files)
  -f FASTA_FILES [FASTA_FILES ...]
                        Specify genome FASTA files for index bilding (accepts
                        wildcards)
  -a                    Whether to report ambiguously mapped contacts
  -k                    Keep any intermediate files (e.g. clipped FASTQ etc).
                        Note initial, primary SAM files are always kept.
  -sam                  Write paired contacts files to SAM format
  -l SEQUENCE           Seek a specific ligation junction sequence (otherwise
                        this is guessed from the primary restriction enzyme)
  -z                    GZIP compress any output FASTQ files
  -v, --verbose         Display verbose messages to report progress
  -u                    Whether to only accept uniquely mapping genome
                        positions and not attempt to resolve certain classes
                        of ambigous mapping.

Note enzymes.conf can be edited to add further restriction enzyme cut-site
defintions. For further help email tjs23@cam.ac.uk or wb104@cam.ac.uk
