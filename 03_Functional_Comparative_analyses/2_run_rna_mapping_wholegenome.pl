#!/usr/bin/perl
use strict;
use warnings;

#########################################################
#####
#
# Create scripts to submit GeMoMa annotation for GAGA assemblies with RNA-seq
#
#####

# Usage
# perl run_rna_mapping_wholegenome.pl


## Input variables

#my $idlist = $ARGV[0];
my $assemblydir = "/home/projects/ku_00039/people/joeviz/GAGA_annotations/Repeat_masked_assemblies/GAGA_annotations/"; # Folder with unzipped soft-masked genome assemblies classified in folders (as in the directory GAGA_annotations in ERDA)
#my $rnadir = "/home/projects/ku_00039/people/joeviz/RNA-seq/GAGA-IDs";
my $rnadir = "/home/projects/ku_00039/people/joeviz/RNA-seq_1June/GAGA-IDs";

#my $gemomainputdir = "/home/projects/ku_00039/people/joeviz/GAGA_annotations/Dataset/Whole_genomes/"; # Folder containing the folders with the annotations used in gemoma

my $gagaidlist = "/home/projects/ku_00039/people/joeviz/GAGA_species_list.txt";

##### Specify here the directory with the dfas folder to process all lgt candidates
system ("ls /home/projects/ku_00039/people/joeviz/LGT_RNA_mapping/dfast/*/genome.gff.fasta > inputgff.txt ");

## Start

# Read species list and prefix
my %shortname;
open(Gagafile, "<", $gagaidlist);
while(<Gagafile>){
	chomp;
	my $line = $_;
	next if ($line !~ /\S+/);
	$line =~ s/\n//g; $line =~ s/\r//g;
	my @subl = split (/\t/, $line);
	$shortname{$subl[0]} = $subl[3];
}
close Gagafile;
#

print "Processing GAGA IDs:\n";

my $numid = "0";
open(File, "<", "inputgff.txt");
while(<File>){
	chomp;
	my $lgtgenome = $_;
	next if ($lgtgenome !~ /\S+/);
	$lgtgenome =~ s/\s+//;


	my $id = "";
	if ($lgtgenome =~ /dfast\/(\S+)\/genome\.gff\.fasta/){
		$id = $1;
	}

	my $ingff = "/home/projects/ku_00039/people/joeviz/LGT_RNA_mapping/dfast/$id/genome.gff_mod.gff";
#	my $ingtf = "/home/projects/ku_00039/people/joeviz/LGT_RNA_mapping/dfast/$id/genome.gff_only.gtf";

#	system ("sed 's/ID=/gene_id /g' $ingff > $ingtf");


	my $gagaid = "";
	if ($lgtgenome =~ /(GAGA\-\d\d\d\d)/ || $lgtgenome =~ /(NCBI\-\d\d\d\d)/ || $lgtgenome =~ /(OUT\-\d\d\d\d)/ ){
		$gagaid = $1;
	}	


	## Find fasta
	system ("ls $assemblydir\/$gagaid\/$gagaid\*softMasked.fasta > inassembly.tmp 2> /dev/null");
	my $assembly = "";
	open (Tmpfile, "<", "inassembly.tmp");
	while (<Tmpfile>){
		chomp;
		my $nline = $_;
		next if ($nline !~ /\S+/);		
		if ($nline =~ /fasta/){
			$assembly = $nline;
		}
	}
	close Tmpfile;
	system ("rm inassembly.tmp");

	if ($assembly !~ /fasta/){
		print "ERROR: Skipping $gagaid Assembly could not be found, Maybe it is gzipped?, in $assemblydir\/$gagaid\/\n";
		next;
	}


	## Find RNA-seq reads
	system ("ls $rnadir\/$gagaid\/\*gz $rnadir\/$gagaid\/\*\/\*gz > inrna.tmp 2> /dev/null");
	my $paired = "";
	my %secondpair;
	my $unpaired = "";
	open (Tmpfile, "<", "inrna.tmp");
	while (<Tmpfile>){
		chomp;
		my $nline = $_;
		next if ($nline !~ /\S+/);		
		if ($nline =~ /(\S+)_1\.fq\.gz/ || $nline =~ /(\S+)_1\.fastq\.gz/ || $nline =~ /(\S+)\_R1\S*\.fastq\.gz/ || $nline =~ /(\S+)\_R1\S*\.fq\.gz/){
			$paired .= "$nline ";
		}
		elsif ($nline =~ /(\S+_)2(\.fq\.gz)/ || $nline =~ /(\S+_)2(\.fastq\.gz)/ || $nline =~ /(\S+\_R)2(\S*\.fastq\.gz)/ || $nline =~ /(\S+\_R)2(\S*\.fq\.gz)/){
			my $firstpair = "$1"."1$2";
			if ($paired =~ /$firstpair /){
				$secondpair{$firstpair} = $nline;
			}
			else {
				print "ERROR: No read pair could be found for $nline. Expected file $firstpair. Skipping this RNA-seq, check the files.\n";
				next;
			}
		}
		elsif ($nline =~ /(\S+)\.fq\.gz/ || $nline =~ /(\S+)\.fastq\.gz/){
			$unpaired .= "$nline ";
		}
	}
	close Tmpfile;
	system ("rm inrna.tmp");

	if ($paired !~ /gz/ && $unpaired !~ /gz/){
		print "ERROR, skipping $gagaid No RNA-seq data could be found in $rnadir\/$gagaid\/\n";
		next;
	}

	# Create script to process and map RNA-seq, and run GeMoMa

	$numid++;
	open ("Outscript", ">", "run_$numid\.sh");
	print Outscript "\#\!\/bin\/bash\n\n";
	print Outscript "mkdir -p $id\ncd $id\n\n";

	# Index genome
	print Outscript "mkdir -p $gagaid\_starindex\n";
	print Outscript "STAR --runThreadN 40 --runMode genomeGenerate --genomeSAindexNbases 13 --genomeDir $gagaid\_starindex --genomeFastaFiles $assembly\n\n";

	# Process RNA-seq

	my $bamfiles = "";


### Paired reads

	my @subpaired = split (/\s/, $paired);
	foreach my $fpair (@subpaired){
		my $spair = $secondpair{$fpair};
		my $rname = "";
		if ($fpair =~ /.*\/(\S+)\.f/){
			$rname = $1;
		} else {
			print "ERROR: Can't find read name in $fpair Skipping...\n";
			next;
		}


		# No trimming

		# Map RNA-seq
		print Outscript "STAR --genomeDir $gagaid\_starindex --runThreadN 40 --readFilesIn $fpair $spair --readFilesCommand zcat --outFileNamePrefix $rname --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 183273429364 --outSAMstrandField intronMotif --outSAMattrIHstart 0 --outFilterIntronMotifs RemoveNoncanonical --alignIntronMin 20 --outFilterMismatchNoverReadLmax 0.01 --outFilterScoreMinOverLread 0.9 --outFilterMatchNminOverLread 0.9\n";

		print Outscript "samtools index $rname"."Aligned.sortedByCoord.out.bam\n\n";

		my $bamfilename = "$rname"."Aligned.sortedByCoord.out.bam";

		$bamfiles .= "ERE.m=$rname"."Aligned.sortedByCoord.out.bam ";

# Quantification
#		print Outscript "featureCounts -a $ingff -p -T 40 -o countMatrix_$rname.txt $rname"."Aligned.sortedByCoord.out.bam\n\n";
		print Outscript "stringtie $bamfilename -G $ingff -e -B -p 40 > $rname\_abundance_fpkm.txt\n";
		print Outscript "mv e_data.ctab $rname\_abundance_counts.txt\n\n";

	}


### Unpaired reads

	my @subunpaired = split (/\s/, $unpaired);
	foreach my $unread (@subunpaired){
		my $rname = "";
		if ($unread =~ /.*\/(\S+)\.f/){
			$rname = "$1";
		} else {
			print "ERROR: Can't find read name in $unread Skipping...\n";
			next;
		}

		# Verify RNA-seq gzip integrity 


		# No trimming

		# Map RNA-seq
		print Outscript "STAR --genomeDir $gagaid\_starindex --runThreadN 40 --readFilesIn $unread --readFilesCommand zcat --outFileNamePrefix $rname --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 183273429364 --outSAMstrandField intronMotif --outSAMattrIHstart 0 --outFilterIntronMotifs RemoveNoncanonical --alignIntronMin 20 --outFilterMismatchNoverReadLmax 0.01 --outFilterScoreMinOverLread 0.9 --outFilterMatchNminOverLread 0.9\n";

		print Outscript "samtools index $rname"."Aligned.sortedByCoord.out.bam\n\n";

		my $bamfilename = "$rname"."Aligned.sortedByCoord.out.bam";

		$bamfiles .= "ERE.m=$rname"."Aligned.sortedByCoord.out.bam ";

		# Quantification
#		print Outscript "featureCounts -a $ingtf -T 40 -o countMatrix_$rname.txt $rname"."Aligned.sortedByCoord.out.bam\n\n";
		print Outscript "stringtie $bamfilename -G $ingff -e -B -p 40 > $rname\_abundance_fpkm.txt\n";
		print Outscript "mv e_data.ctab $rname\_abundance_counts.txt\n\n";
	}


# Quantification

#module load subread/2.0.2

### Assign reads to exons based on the GTF annotation:
# --a = the GTF/GFF file
# --F = specify the format of -a
# --p = data are paired-end (unset it if you have single-end)
# --T = set number of threads
# -P  = only consider pairs with ISIZE defined by -d & -D, default 50-600bp
# -o  = output file
#featureCounts -a $GTF -F GTF -p -T 8 -P -o countMatrix_all.txt *.bam	

#print Outscript "featureCounts -a $ingff -F GFF -p -T 40 -P -o countMatrix_all.txt *.bam";



# stringtie
#module load stringtie/2.1.5
#stringtie in.bam -G gff -e -B 

#	print Outscript "cat *counts.txt > RNAseq_mapping_count_summary.txt\n\n";
	print Outscript "cat *counts.txt | head -1 > RNAseq_mapping_count_summary.txt\n";
	print Outscript "grep -P \'\\d\' *counts.txt >> RNAseq_mapping_count_summary.txt\n";
#	print Outscript "cat *fpkm.txt | head -1 > RNAseq_mapping_fpkm_summary.txt\n";
#	print Outscript "grep -P \'\\d\' *fpkm.txt >> RNAseq_mapping_fpkm_summary.txt\n";
	print Outscript "grep -v \'\#\' *fpkm.txt >> RNAseq_mapping_fpkm_summary.txt\n";

	print "$gagaid OK\n";
}
close File;

print "\n";


## Create submission script

open ("Submitscript", ">", "submit_run_gemoma_all\.sh");
print Submitscript "#!/bin/sh
#PBS -W group_list=ku_00039 -A ku_00039
#PBS -m n
#PBS -l nodes=1:ppn=40
#PBS -l mem=180gb
#PBS -l walltime=24:00:00
#PBS -t 1-$numid

# Go to the directory from where the job was submitted (initial directory is \$HOME)
echo Working directory is \$PBS_O_WORKDIR
cd \$PBS_O_WORKDIR

### Here follows the user commands:
# Define number of processors
NPROCS=`wc -l < \$PBS_NODEFILE`
echo This job has allocated \$NPROCS nodes

# Load all required modules for the job

module load tools
module load ngs
module load anaconda3/4.4.0
module load perl
module load gcc
module load star/2.7.2b
module load samtools/1.12
module load ncbi-blast/2.2.31+
module load java/1.8.0
module load mmseqs2/release_12-113e3
module load subread/2.0.2
module load stringtie/2.1.5

# This is where the work is done

bash run_\$\{PBS_ARRAYID\}.sh


";



