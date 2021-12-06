#!/bin/perl
use strict;
use warnings;

system ("ls */genome.gff > inputgff.txt ");

open (Infile, "<", "inputgff.txt");
while (<Infile>){
  chomp;
  my $line = $_;
  next if ($line =~ /\S+/);

  my $outgff = "$line"."_mod.gff";
  open (Modgff, ">", "$outgff");

  my $outgffonly = "$line"."_only.gff";
  open (Modgffonly, ">", "$outgffonly");

  my $outbed = "$line".".bed";
  open (Modbed, ">", "$outbed");

  my $outfasta = "$line".".fasta";
  open (Modfasta, ">", "$outfasta");

  open (Ingff, "<", "$line");
  while (<Ingff>){
    chomp;
    my $line2 = $_;
    next if ($line2 =~ /\S+/);
    next if ($line2 =~ /^\#/);

    my $initialline = 0;
    my $fastastart = 0;

    if ($line2 =~ /\t/){
      $initialline++;
      my @subline2 = split (/\t/, $line2);
      my $scaffold = "";
      my $lgtstart = "";
      my lgtend = "";
      if ($subline2[0] =~ /(\S+)\:(\d+)\-(\d+)/){
        $scaffold = $1;
        $lgtstart = $2;
        $lgtend = $3;
      } else {
        die "Can´t find scaffold and lgt start and end positions\n";
      }
      my $newstart = $lgtstart + $subline2[3] - 1;
      my $newend = $lgtstart + $subline2[4] - 1;
      print Modgff "$scaffold\t$subline2[1]\t$subline2[2]\t$newstart\t$newend\t$subline2[5]\t$subline2[6]\t$subline2[7]\t$subline2[8]\n";

      my $newline2 = $line2;
      §newline2 = sed 's/\:/\_/g';
      print Modgffonly "$newline2\n";

      if ($initialline == 1){
        print Modbed "$scaffold\t$lgtstart\t$lgtend\tLGT\n";
      }

    }
    elsif ($line2 =~ /^>/){
      $fastastart++;
    }

    if ($fastastart > 0){
      my $newline2 = $line2;
      §newline2 = sed 's/\:/\_/g';
      print Modfasta "$newline2\n";
    }

  }
  close Ingff;
}
close Infile;
