NucProcess
----------

NucProcess is a Python program to perform single-cell Hi-C sequence processing
of paired read data. This software takes paired FASTQ sequence read files
(representing Hi-C data for only ONE cell), a reference genome sequence and
knowledge of experimental parameters (restriction enzyme type and the range of
size of DNA fragments sequenced in the library) to create processed, single-cell
Hi-C contact files. By default, the output is generated in Nuc Chromatin Contact
(NCC) format, which is a simple text based format described below. A processing
report document and a contact map will also be generated for each run; in SVG
format that can be viewed in most web browsers.

To run NucProcess issue the 'nuc_process' command with the command line options
described below. The options -i (input FASTQ files) and -g (genome reference)
are mandatory, though its is usual to also use -re1 (primary restriction enzyme.
Default is MboI), -o (root name of output files) and -re2 (secondary restriction
enzyme in double-digest experiments).

The 'nuc_contact_map' command takes the contact data from NCC format files to
make all-chromosome contact map graphics in SVG format. This is automatically
run on the main output of NucProcess, but can be run as required on any NCC
format file.

The 'nuc_contact_probability' command takes the contact data from one or more
NCC format files to create log plots of contact probability versus sequence
separation for intra chromosomal contacts.

IMPORTANT: If the files for the genome index and the corresponding restriction
enzyme (RE) cut location files have not been created they will be generated
automatically by NucProcess so long as the -f option is specified. This option
specifies the location of complete chromosomal sequences in FASTA format,
typically by using a wild-card specification. Once these files are present the -f
option need not be specified, but it will not trigger the re-creation of the
files unless the -m or -x options are used (or if files are deleted).

NOTE: the names of the chromosomes used by NucProcess are determined by the file
names of the chromosome sequence files that were used to build the genome index
and RE files (and these must match). The chromosome sequence files should be
tagged with the chromosome name after an underscore but before the file
extension. For example "chr7" will be taken from a file named
"mm_ref_GRCm38.p2_chr7.fa".


Citations
---------

If you use NucProcess in published work, please cite the following reference:

  Stevens et al. Nature. 2017 Apr 6;544(7648):59-64 [PMID:28289288]



Barcoded input
--------------

The splitFastqBarcodes.py script is provided to split FASTQ files that represent
many cells, each with a different barcode sequence, into separate paired read
files. The script can be run as follows, specifying the names of the two paired,
multiplex FASTQ read files after the script:

  python splitFastqBarcodes.py MULTIPLEXED_DATA_r_1.fq MULTIPLEXED_DATA_r_2.fq

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

-o specifies CELL_1 will be used for naming the output. In this case the main
output contact file will be CELL_1.ncc

-v specifies verbose output of processing progress

-a specifies to generate ambiguous contact files: CELL_1_ambig.ncc in this case

-k specifies to keep all the intermediate processing files: Filtered NCC files,
clipped FASTQ files and the main Bowtie2 mapping SAM files 

-re1 is the primary restriction enzyme type at the ligation junction (see
enzymes.conf)

-re2 is the secondary restriction enzyme used to release the fragments (option
not used for Nextera based protocol)

-s is the valid molecule/fragment size range, as used in the DNA sequencing

-n is the number of parallel CPU cores to use with Bowtie2

-g GENOME_BUILD is the root name for the Bowtie2 genome index without any file
extension and in this case would refer to files GENOME_BUILD.1.bt2,
GENOME_BUILD.rev.1.bt2 etc.

-i SEQUENCING_DATA_r_?.fq is a wild-card expression matching the two input FASTQ
files (though two separate file names, separated by a space can be specified).
In this case the expression would match SEQUENCING_DATA_r_1.fq and
SEQUENCING_DATA_r_2.fq - the paired sequence read files.


To generate contact map graphics from an NCC format file:

  nuc_contact_map -i CELL_1.ncc


This will generate the output graphics file CELL_1_contact_map.svg. However, the
output file name maye be specified via the -o option.



NCC data format
---------------

The main NCC output format for contact data takes the form of space-separate
lines, where each line represents a different chromosomal contact that pairs two
chromosomal regions. This format specifies more than just the chromosome contact
map: it also records the original read locations within the (ligated) primary RE
digest fragments, strand information, ambiguity information and information to
relate the data back to the original FASTQ input files. Accordingly this format
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



Command line options for nuc_process
------------------------------------

usage: nuc_process [-h] [-i FASTQ_FILE [FASTQ_FILE ...]] [-g GENOME_FILE]
                   [-g2 GENOME_FILE_2] [-re1 ENZYME] [-re2 ENZYME]
                   [-s SIZE_RANGE] [-n CPU_COUNT] [-r COUNT] [-o NCC_FILE]
                   [-oa NCC_FILE] [-or REPORT_FILE] [-b EXE_FILE] [-q SCHEME]
                   [-qm MIN_QUALITY] [-m] [-p]
                   [-pt PAIRED_READ_TAGS PAIRED_READ_TAGS] [-x]
                   [-f FASTA_FILES [FASTA_FILES ...]]
                   [-f2 FASTA_FILES_2 [FASTA_FILES_2 ...]] [-a] [-k] [-sam]
                   [-l SEQUENCE] [-z] [-v] [-hc HOM_CHROMO_TSV_FILE] [-u]
                   [-c GENOME_COPIES]

Chromatin contact paired-read Hi-C processing module for Nuc3D and NucTools

optional arguments:

  -h, --help            show this help message and exit
  -i FASTQ_FILE [FASTQ_FILE ...]
                        Input paired-read FASTQ files to process. Accepts
                        wildcards that match paired files. If more than two
                        files are input, processing will be run in batch mode
                        using the same parameters.
  -g GENOME_FILE        Location of genome index files to map sequence reads
                        to without any file extensions like ".1.b2" etc. A new
                        index will be created with the name if the index is
                        missing and genome FASTA files are specified
  -g2 GENOME_FILE_2     Location of secondary genome index files for hybrid
                        genomes. A new index will be created with the name if
                        the index is missing and genome FASTA files are
                        specified
  -re1 ENZYME           Primary restriction enzyme (for ligation junctions).
                        Default: MboI. Available: AluI, BglII, DpnII, HindIII,
                        MboI
  -re2 ENZYME           Secondary restriction enzyme (if used). Available:
                        AluI, BglII, DpnII, HindIII, MboI
  -s SIZE_RANGE         Allowed range of sequenced molecule sizes, e.g.
                        "150-1000", "100,800" or "200" (no maximum)
  -n CPU_COUNT          Number of CPU cores to use in parallel
  -r COUNT              Minimum number of sequencing repeats required to
                        support a contact
  -o NCC_FILE           Optional output name for NCC format chromosome contact
                        file. This option will be ignored if more than two
                        paired FASTA files are input (i.e. for batch mode);
                        automated naming will be used instead.
  -oa NCC_FILE          Optional output name for ambiguous contact NCC file.
                        This option will be ignored if more than two paired
                        FASTA files are input (i.e. for batch mode); automated
                        naming will be used instead.
  -or REPORT_FILE       Optional output name for SVG format report file. This
                        option will be ignored if more than two paired FASTA
                        files are input (i.e. for batch mode); automated
                        naming will be used instead.
  -b EXE_FILE           Path to bowtie2 (read aligner) executable (will be
                        searched for if not specified)
  -q SCHEME             Use a specific FASTQ quality scheme (normally not set
                        and deduced automatically). Available: phred33,
                        phred64, solexa
  -qm MIN_QUALITY       Minimum acceptable FASTQ quality score in range 0-40
                        for clipping 3' end of reads. Default: 10
  -m                    Force a re-mapping of genome restriction enzyme sites
                        (otherwise cached values will be used if present)
  -p                    The input data is multi-cell/population Hi-C; single-
                        cell processing steps are avoided
  -pt PAIRED_READ_TAGS PAIRED_READ_TAGS
                        When more than two FASTQ files are input (batch mode),
                        the subtrings/tags which differ between paired FASTQ
                        file paths. Default: r_1 r_2
  -x, --reindex         Force a re-indexing of the genome (given appropriate
                        FASTA files)
  -f FASTA_FILES [FASTA_FILES ...]
                        Specify genome FASTA files for genome index building
                        (accepts wildcards)
  -f2 FASTA_FILES_2 [FASTA_FILES_2 ...]
                        A second set of genome FASTA files for building a
                        second genome index when using hybrid strain cells
                        (accepts wildcards).
  -a                    Whether to report ambiguously mapped contacts
  -k                    Keep any intermediate files (e.g. clipped FASTQ etc).
  -sam                  Write paired contacts files to SAM format
  -l SEQUENCE           Seek a specific ligation junction sequence (otherwise
                        this is guessed from the primary restriction enzyme)
  -z                    GZIP compress any output FASTQ files
  -v, --verbose         Display verbose messages to report progress
  -hc HOM_CHROMO_TSV_FILE, --homologous_chromos HOM_CHROMO_TSV_FILE
                        File path specifying whitespace-separated pairs of
                        homologous chromosome names (to match first word in
                        header lines of genome sequence FASTQ files) and
                        corresponding pair name (e.g. "chrX"). Required for
                        hybrid strain analysis. See -g2 option.
  -u                    Whether to only accept uniquely mapping genome
                        positions and not attempt to resolve certain classes
                        of ambiguous mapping where a single perfect match is
                        found.
  -c GENOME_COPIES      Number of whole-genome copies, e.g. for S2 phase;
                        Default 1 unless homologous chromosomes (-hc) are
                        specified for hybrid genomes.

Note enzymes.conf can be edited to add further restriction enzyme cut-site
definitions. 



Command line options for nuc_contact_map
----------------------------------------

usage: nuc_contact_map [-h] [-i NCC_FILE [NCC_FILE ...]] [-o SVG_FILE_TAG]
                       [-w SVG_WIDTH] [-s BIN_SIZE] [-b] [-c RGB_COLOR]

Chromatin contact (NCC format) Hi-C contact map display module for Nuc3D and
NucTools

optional arguments:
  -h, --help    show this help message and exit
  -i NCC_FILES  Input NCC format chromatin contact file(s). Wildcards accepted
  -o SVG_FILE   Optional name tag to put at end of SVG format contact
                map file. Use "-" to print SVG to stdout rather than
                make a file. Default: "_contact_map"
  -w SVG_WIDTH  SVG document width
  -s BIN_SIZE   Sequence region size represented by each small square (the
                resolution) in megabases. Default is 5 kb
  -b            Specifies that the contact map should have a black background
                (default is white)
  -c RGB_COLOR  Optional main color for the contact points as a 24-bit
                hexidecimal RBG code e.g. "#0080FF" (with quotes)

For further help email tjs23@cam.ac.uk or wb104@cam.ac.uk


Command line options for nuc_contact_probability
------------------------------------------------

usage: nuc_contact_probability [-h] [-i NCC_FILE [NCC_FILE ...]] [-o SVG_FILE]
                               [-w SVG_WIDTH] [-s KB_BIN_SIZE]

Chromatin contact (NCC format) probability vs sequence separation graph module
for Nuc3D and NucTools

optional arguments:
  -h, --help            show this help message and exit
  -i NCC_FILE [NCC_FILE ...]
                        Input NCC format chromatin contact file(s). Wildcards
                        accepted
  -o SVG_FILE           Output SVF format file. Use "-" to print SVG to stdout
                        rather than make a file.
  -w SVG_WIDTH          SVG document width
  -s KB_BIN_SIZE        Sequence region size in kilobases for calculation of
                        contact probabilities. Default is 100 (kb)
