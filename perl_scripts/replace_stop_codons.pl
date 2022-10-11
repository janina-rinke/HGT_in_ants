#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
my $dirname = dirname(__FILE__);


# usage: perl mask_stop_zorro_codon_alignments.pl

## Load computerome2 modules
#module load ngs tools
#module load anaconda3/4.4.0


my $cdsfile = "$ARGV[0]";
my $outfile = "$ARGV[1]";

#system ("mkdir -p $outdir");


my ($line, $name);
my %fasta;
my $ognumseqs = "0"; my $ogalnlength = "0";

open(Filef, "<", $cdsfile);
while(<Filef>){
    chomp;
    my $line2 = $_;
    if ($line2 =~ />(\S+)/){
        $name = $1;
        $ognumseqs++;
    } else {
        $fasta{$name} .= "$line2";
        if ($ognumseqs == 1){ # only count aln length for the first sequence 
            $ogalnlength += length($line2);
        }
    }
}
close Filef;

open (Results, ">", $outfile);

foreach my $gene (keys %fasta){
    my $seq = uc($fasta{$gene});
    my $nstopseq = "";
    my $multseq = "";

    my $lengthcheck = length($seq)/3;
    if ($lengthcheck =~ /\.33/){
        #print Results "$line2"."NN\n";	
        $multseq = "$seq"."NN";
        print "$gene is not multiple of 3!\n";	# Print seq id to check
    } elsif ($lengthcheck =~ /\.66/) {
        #print Results "$line2"."N\n";
        $multseq = "$seq"."N";
        print "$gene is not multiple of 3!\n";	# Print seq id to check
    } else {
        #print Results "$line2\n";
        $multseq = "$seq";
    }

    # Delete stop codonz
    my $length = length($multseq);
    for (my $n = 0; $n < ($length - 2); $n += 3) {
        my $codon = substr ($multseq, $n, 3);
        if ($codon =~ /TAA/ || $codon =~ /TAG/ || $codon =~ /TGA/ || $codon =~ /\S\SN/){ # also discard codon finished in Ns
#				print Nonstop "NNN";
            $nstopseq .= "NNN";	
#				if ($n >= ($length - 3)){ # Only discard stop if it is the last one
#					$nstopseq .= "";	
#				} else {
#					$nstopseq .= "$codon";	
#				}
        }
        else {
#				print Nonstop "$codon";
            $nstopseq .= "$codon";
        }
    }

    print Results ">$gene\n$nstopseq\n";

}

close Results;

