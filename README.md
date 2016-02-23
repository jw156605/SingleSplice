# SingleSplice
Algorithm for detecting alternative splicing in a population of single cells

System Requirements:
- R, Perl, gcc, and git
- Boost C++ library

Installation instructions:
- Check that your system meets the requirements
- Navigate to the desired installation directory
- Clone the SingleSplice repository:
  git clone https://github.com/jw156605/SingleSplice
- To install the Boost library, you can run the following commands:
  cd SingleSplice/
  wget https://sourceforge.net/projects/boost/files/boost/1.60.0/boost_1_60_0.tar.gz
  tar -xvf boost_1_60_0.tar.gz
- If you install in a directory other than SingleSplice/boost_1_60_0/, you need to modify the BOOST
  variable in the Makefile to point to the installation directory.
- Run the Makefile by simply typing: make

Sample input files for SingleSplice are provided. The FASTQ and BAM files are not included due to their large size.
To reproduce these results, follow the set of steps below:

- Download the FASTQ files from the 80 E18.5 cells in GEO record GSE52583.
- Align the reads to mm10 using your favorite RNA-seq aligner (MapSplice, TopHat, STAR) and convert to indexed BAM.
- Run diffsplice on the indexed BAM files:

cd diffsplice
bin/diffsplice -o <output_dir> -m full -s ../sample_input/sample_config.txt ../sample_input/sample_datafile.txt

- Align the reads to the ERCC transcript sequences using your favorite unspliced aligner (e.g., bowtie2).
- Count the numbers of mapped biological and spike-in reads and use this information to compute "cell size" values.
- Compute normalized coverage in units of reads per kilobase per million reads (RPKM) for the ERCC transcripts. Normalize to median cell size
as described in the paper.
- Run SingleSplice using the ASM abundance estimates contained in result/asm/. Note that diffsplice outputs a file for each chromosome, so
to run on all chromosomes, these files must be concatenated.
perl SingleSplice.pl -a sample_input/asm_all_chr.txt -p 1000 -s 80 -t sample_input/ERCC_rpkms_size_norm_median.csv -r sample_input/total_reads_DistalLungEpithelium.csv -g sample_input/groups.csv