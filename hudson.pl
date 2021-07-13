#!/usr/bin/env perl

## author: reubwn July 2021

use strict;
use warnings;

use Getopt::Long;
use Sort::Naturally;
use List::Util qw(first);
use Data::Dumper;

my $usage = "
SYNOPSIS:
  Detect recombination breakpoints from WhatsHap phased blocks using Hudson's
  4-gamete test.

USAGE:
  hudson.pl -i <VCF> [-b <BLOCKS>] [-o <OUT>]

OPTIONS:
  -i|--vcf    [FILE] : VCF file with phased genotypes [required]
  -1|--sample1 [STR] : sample 1 ID [first in VCF]
  -2|--sample2 [STR] : sample 2 ID [second in VCF]
  -b|--blocks [FILE] : blocks info file (output using `whatshap stats`)
  -o|--out     [STR] : prefix for outfiles ['collated']
  -h|--help          : prints this help message
\n";

my ($vcf_file, $sample1_id, $sample2_id, $blocks_file, $help, $debug);
my $sample1_idx = 9;
my $sample2_idx = 10;
my $out_prefix = "hudson";

GetOptions (
  'i|vcf=s'     => \$vcf_file,
  '1|sample1:s' => \$sample1_id,
  '2|sample2:s' => \$sample2_id,
  'b|blocks:s'  => \$blocks_file,
  'o|out:s'     => \$out_prefix,
  'd|debug'     => \$debug,
  'h|help'      => \$help
);

die $usage if $help;
die $usage unless ($vcf_file);

my $VCF;
if ($vcf_file =~ m/.gz$/) {
  open ($VCF, "gunzip -c $vcf_file |") or die $!;
} else {
  open ($VCF, $vcf_file) or die $!;
}

my %seq_lengths;
my %block_lengths;
my %genotypes;
# my %blocks;

while (my $line = <$VCF>) {
  chomp $line;
  my @F = split (/\s+/, $line);

  if ( $line =~ /^##contig=<ID=(\w+),length=(\d+)>/ ) {
    print STDERR "Contig ID: $1, length: $2\n" if $debug;
    $seq_lengths{$1} = $2;
  } elsif ( $line =~ /^#CHROM/ ) {
    if ( ($sample1_id) && ($sample2_id) ) {
      $sample1_idx = first { $F[$_] eq $sample1_id } 0..$#F;
      $sample2_idx = first { $F[$_] eq $sample2_id } 0..$#F;
      print STDERR "Sample ID '$sample1_id' has index $sample1_idx in VCF\n" if $debug;
      print STDERR "Sample ID '$sample2_id' has index $sample2_idx in VCF\n" if $debug;
    } else {
      print STDERR "Samples are $F[$sample1_idx] and $F[$sample2_idx]\n" if $debug;
    }
  } elsif (scalar(@F)>=8) {
    ## only process SNVs with REF and ALT == 1 bp
    if ( (length($F[3])==1) && (length($F[4])==1) ) {
      print STDERR "\r[INFO] Processing SNV at position $F[0]:$F[1] $F[3]/$F[4]    "; $| = 1;

      my @ff1 = split (":", $F[$sample1_idx]);
      my @ff2 = split (":", $F[$sample2_idx]);

      ## FORMAT fields
      ## GT:DP:AD:AF:RAF:NVC:FLG:SB:SC:PF:PS
      ## 0  1  2  3  4   5   6   7  8  9  10
      ## FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
      ## FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
      ## FORMAT=<ID=AD,Number=1,Type=Integer,Description="Allele Depth">
      ## FORMAT=<ID=AF,Number=1,Type=Float,Description="Allele Fraction">
      ## FORMAT=<ID=RAF,Number=1,Type=Float,Description="Revised Allele Fraction">
      ## FORMAT=<ID=NVC,Number=1,Type=Integer,Description="Nearby Variation Count">
      ## FORMAT=<ID=FLG,Number=1,Type=Integer,Description="Flagged">
      ## FORMAT=<ID=SB,Number=1,Type=Float,Description="Strand Bias">
      ## FORMAT=<ID=SC,Number=1,Type=Float,Description="Score">
      ## FORMAT=<ID=PF,Number=1,Type=String,Description="Pass Filter">
      ## FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set identifier">
      if ( ($ff1[9] eq "PASS") && ($ff2[9] eq "PASS") ) { ## both PASS
        if ( ($ff1[0] =~ /\d\|\d/) && ($ff2[0] =~ /\d\|\d/) ) { ## both genotypes phased
          if ( $ff1[10] == $ff2[10] ) { ## both part of same phase set

          ## blocks data
          $genotypes{$F[0]}{$ff1[10]}{$F[1]}{GT1} = $ff1[0];
          $genotypes{$F[0]}{$ff1[10]}{$F[1]}{GT2} = $ff2[0];
          $genotypes{$F[0]}{$ff1[10]}{$F[1]}{REF} = $F[3];
          $genotypes{$F[0]}{$ff1[10]}{$F[1]}{ALT} = $F[4];

          }
        }
      }
    }
  } else {
    next;
  }
}
print STDERR "\n";

print Dumper(\%genotypes) if $debug;

## open outfile
open (my $OUT, ">".$out_prefix."_results.txt") or die $!;
print $OUT "#CHROM\tBLOCK\tSNV1_POS\tSNV2_POS\tGT:$sample1_id.1\tGT:$sample1_id.2\tGT:$sample2_id.1\tGT:$sample2_id.2\tNUM\tRECOMB\n";
## open fasta blocks file
open (my $FASTA, ">".$out_prefix."_blocks.fasta") or die $!;
## open haplotypes info file
open (my $HAP, ">".$out_prefix."_hap.txt") or die $!;
# print $HAP "#LOCATION\tGT:$sample1_id\tGT:$sample2_id\tEVENT\n";

foreach my $chrom (nsort keys %genotypes) {
  my %blocks = %{$genotypes{$chrom}};
  foreach my $block (sort {$a<=>$b} keys %blocks) {
    print $HAP "##block_id=$block\n";
    my %positions = %{$blocks{$block}};
    my (@prev_gt1, @prev_gt2, $prev_position);
    my (@hap1, @hap2, @hap3, @hap4);
    foreach my $curr_position (sort {$a<=>$b} keys %positions) {
      print $HAP "$chrom:$curr_position\t$positions{$curr_position}{GT1}\t$positions{$curr_position}{GT2}";
      my @curr_gt1 = split (/\|/, $positions{$curr_position}{GT1});
      my @curr_gt2 = split (/\|/, $positions{$curr_position}{GT2});
      $curr_gt1[0] == 0 ? push (@hap1, $positions{$curr_position}{REF}) : push (@hap1, $positions{$curr_position}{ALT});
      $curr_gt1[1] == 0 ? push (@hap2, $positions{$curr_position}{REF}) : push (@hap2, $positions{$curr_position}{ALT});
      $curr_gt2[0] == 0 ? push (@hap3, $positions{$curr_position}{REF}) : push (@hap3, $positions{$curr_position}{ALT});
      $curr_gt2[1] == 0 ? push (@hap4, $positions{$curr_position}{REF}) : push (@hap4, $positions{$curr_position}{ALT});

      ## evaluation block:
      if ( (@prev_gt1) && (@prev_gt2) && ($prev_position) ) {
        ## make haplotypes
        my $hap1 = join ("", $prev_gt1[0], $curr_gt1[0]);
        my $hap2 = join ("", $prev_gt1[1], $curr_gt1[1]);
        my $hap3 = join ("", $prev_gt2[0], $curr_gt2[0]);
        my $hap4 = join ("", $prev_gt2[1], $curr_gt2[1]);

        ## count unique haplotypes
        my %haplotypes;
        $haplotypes{$hap1}++;
        $haplotypes{$hap2}++;
        $haplotypes{$hap3}++;
        $haplotypes{$hap4}++;

        ## print to file
        print $OUT join ("\t", $chrom, $block, $prev_position, $curr_position, $hap1, $hap2, $hap3, $hap4, ($hap1 =~ tr/01//), scalar(keys %haplotypes));
        # scalar(keys %haplotypes) > 2 ? print $OUT "\tYES\n" : print $OUT "\tNO\n";
        if (scalar(keys %haplotypes) > 2) {
          print $OUT "\tYES\n";
          print $HAP "\t*\n";
        } else {
          print $OUT "\tNO\n";
          print $HAP "\t.\n";
        }

        ## record end position of block
        $block_lengths{"$chrom:$block"}{END} = $curr_position;

      } else {
        ## record start position of block
        $block_lengths{"$chrom:$block"}{START} = $curr_position;
        print $HAP "\t.\n";
      }
      @prev_gt1 = split (/\|/, $positions{$curr_position}{GT1});
      @prev_gt2 = split (/\|/, $positions{$curr_position}{GT2});
      $prev_position = $curr_position;
    }
    ## print haplotype blocks as fasta
    print $FASTA ">$sample1_id:$chrom:$block:1\n";
    print $FASTA join ("", @hap1) . "\n";
    print $FASTA ">$sample1_id:$chrom:$block:2\n";
    print $FASTA join ("", @hap2) . "\n";
    print $FASTA ">$sample2_id:$chrom:$block:1\n";
    print $FASTA join ("", @hap3) . "\n";
    print $FASTA ">$sample2_id:$chrom:$block:2\n";
    print $FASTA join ("", @hap4) . "\n";
    print $FASTA "\n";
  }
}
close $OUT;
close $FASTA;
close $HAP;

foreach (nsort keys %block_lengths) {
  print "$_ start = $block_lengths{$_}{START}\n";
  print "$_ end = $block_lengths{$_}{END}\n";
  my $len = $block_lengths{$_}{END} - $block_lengths{$_}{START};
  print "$_ length = $len\n";
}

print STDERR "[INFO] Finished " . `date`;

__END__
