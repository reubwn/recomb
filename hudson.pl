#!/usr/bin/env perl

## author: reubwn July 2021

use strict;
use warnings;

use Getopt::Long;
use Sort::Naturally;
use List::Util qw(first sum);
use Data::Dumper;

my $usage = "
SYNOPSIS:
  Detect recombination breakpoints from phased blocks using Hudson's 4-gamete test.

USAGE:
  hudson.pl -i <VCF> [-1 SAMPLE1 -2 SAMPLE2] [-o <OUT>]

OPTIONS:
  -i|--vcf    [FILE] : VCF file with phased genotypes [required]
  -1|--sample1 [STR] : sample 1 ID [first in VCF]
  -2|--sample2 [STR] : sample 2 ID [second in VCF]
  -o|--out     [STR] : prefix for outfiles ['recomb']
	-q|--quiet         : don't print progress
  -v|--verbose       : say more stuff
  -d|--debug         : say even more stuff
  -h|--help          : prints this help message
\n";

my ($vcf_file, $sample1_id, $sample2_id, $help, $quiet, $verbose, $debug);
my $sample1_idx = 9; ## default first and second samples in file
my $sample2_idx = 10;
my $out_prefix = "recomb";

GetOptions (
  'i|vcf=s'     => \$vcf_file,
  '1|sample1:s' => \$sample1_id,
  '2|sample2:s' => \$sample2_id,
  'o|out:s'     => \$out_prefix,
	'q|quiet'     => \$quiet,
  'v|verbose'   => \$verbose,
  'd|debug'     => \$debug,
  'h|help'      => \$help
);

die $usage if $help;
die $usage unless ($vcf_file);

if (system ("bcftools --version &> /dev/null") != 0) {
	die "\n[ERROR] Problem with bcftools!\n\n";
}

my $VCF;
## open with bcftools, should handle most file types
open ($VCF, "bcftools view $vcf_file |") or die $!;

my (%genotypes, %blocks, %seq_lengths, %supporting_sites);
my ($sum_of_blocks, $sum_of_recomb_tracts, $recombination_events) = (0,0,0,0);
print STDERR "\n";

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
      print STDERR "[INFO] Sample 1 ID: '$sample1_id' (position $sample1_idx in VCF)\n";
      print STDERR "[INFO] Sample 2 ID: '$sample2_id' (position $sample2_idx in VCF)\n";
    } else {
      print STDERR "[INFO] No sample IDs provided!\n";
      $sample1_id = $F[$sample1_idx];
      $sample2_id = $F[$sample2_idx];
      print STDERR "[INFO] Defaulting to 1st and 2nd samples in VCF: '$sample1_id', '$sample2_id'\n";
    }
  } elsif (scalar(@F)>=8) {
    ## only process SNVs with REF and ALT == 1 bp
    if ( (length($F[3])==1) && (length($F[4])==1) ) {
			
			unless ( $quiet ) {
				print STDERR "\r[INFO] Processing SNV at position $F[0]:$F[1] $F[3]/$F[4]    "; $| = 1;
			}

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
close $VCF;
print STDERR "\n[INFO] Done\n\n";

## add sample names to output filenames
$out_prefix .= ".$sample1_id.$sample2_id";

## open outfile
open (my $OUT, ">".$out_prefix.".hudson_results.txt") or die $!;
print $OUT "#CHROM:BLOCK\tSAMPLES\tSNP1_POS\tSNP2_POS\tDISTANCE\tGT:$sample1_id.1\tGT:$sample1_id.2\tGT:$sample2_id.1\tGT:$sample2_id.2\tHAP:$sample1_id.1\tHAP:$sample1_id.2\tHAP:$sample2_id.1\tHAP:$sample2_id.2\tNUM_HAPS\tRECOMB\n";
## open per-block outfile
open (my $BLOCKS, ">".$out_prefix.".block_results.txt") or die $!;
print $BLOCKS "#CHROM:BLOCK\tCOORDS\tSAMPLE1\tSAMPLE2\tLENGTH_SNPS\tLENGTH_BP\tNUM_RECOMBS\tMEAN_AFFECTED_SNPS\tMEAN_AFFECTED_BP\n";
## open fasta blocks file
open (my $FASTA, ">".$out_prefix.".blocks.fasta") or die $!;
## open haplotypes info file
open (my $HAP, ">".$out_prefix.".hap.txt") or die $!;
## open tracts info file
open (my $TRACTS, ">".$out_prefix.".tracts.txt") or die $!;
print $TRACTS "#CHROM:BLOCK\tNUM_RECOMBS\tRECOMB_POS\tWINDOW_START\tWINDOW_END\tSUPPORT\tTRACT_LENGTH\n";

foreach my $chrom (nsort keys %genotypes) {
  my %blocks = %{$genotypes{$chrom}};

  foreach my $block (sort {$a<=>$b} keys %blocks) {
    print $HAP "##block_id=$block\n";
    my %positions = %{$blocks{$block}};
    my (@ALL_positions, @RECOMB_indices, $prev_position);
    my (@prev_gt1, @prev_gt2, @prev_nuc1, @prev_nuc2);
    my (@hap1_nuc_array, @hap2_nuc_array, @hap3_nuc_array, @hap4_nuc_array);
    my (@hap1_geno_array, @hap2_geno_array, @hap3_geno_array, @hap4_geno_array);
    my (@affected_snps_per_block, @affected_bp_per_block);
    my $i = 0;

    foreach my $curr_position (sort {$a<=>$b} keys %positions) {
      print $HAP "$chrom:$block:$curr_position\t$positions{$curr_position}{GT1}\t$positions{$curr_position}{GT2}\t";
      my @curr_gt1 = split (/\|/, $positions{$curr_position}{GT1});
      my @curr_gt2 = split (/\|/, $positions{$curr_position}{GT2});
      ## construct haplotype nucleotide arrays
      my (@curr_nuc1, @curr_nuc2);
      $curr_gt1[0] == 0 ? push (@curr_nuc1, $positions{$curr_position}{REF}) : push (@curr_nuc1, $positions{$curr_position}{ALT});
      $curr_gt1[1] == 0 ? push (@curr_nuc1, $positions{$curr_position}{REF}) : push (@curr_nuc1, $positions{$curr_position}{ALT});
      $curr_gt2[0] == 0 ? push (@curr_nuc2, $positions{$curr_position}{REF}) : push (@curr_nuc2, $positions{$curr_position}{ALT});
      $curr_gt2[1] == 0 ? push (@curr_nuc2, $positions{$curr_position}{REF}) : push (@curr_nuc2, $positions{$curr_position}{ALT});
      ## print genotypes as nucleotides
      print $HAP join("\t", join("|",@curr_nuc1), join("|",@curr_nuc2));

      ## evaluation block:
      if ( (@prev_gt1) && (@prev_gt2) && ($prev_position) ) {
        ## make haplotypes of pairs of adjacent SNPs
        my $hap1_gt = join ('', $prev_gt1[0], $curr_gt1[0]);
        my $hap2_gt = join ('', $prev_gt1[1], $curr_gt1[1]);
        my $hap3_gt = join ('', $prev_gt2[0], $curr_gt2[0]);
        my $hap4_gt = join ('', $prev_gt2[1], $curr_gt2[1]);
        ## and as nucs
        my $hap1_nuc = join ('', $prev_nuc1[0], $curr_nuc1[0]);
        my $hap2_nuc = join ('', $prev_nuc1[1], $curr_nuc1[1]);
        my $hap3_nuc = join ('', $prev_nuc2[0], $curr_nuc2[0]);
        my $hap4_nuc = join ('', $prev_nuc2[1], $curr_nuc2[1]);

        ## count unique haplotypes
        my %haplotypes;
        $haplotypes{$hap1_gt}++;
        $haplotypes{$hap2_gt}++;
        $haplotypes{$hap3_gt}++;
        $haplotypes{$hap4_gt}++;

        ## print to file
        print $OUT join ("\t", "$chrom:$block", "$sample1_id,$sample2_id", $prev_position, $curr_position, ($curr_position-$prev_position), $hap1_gt, $hap2_gt, $hap3_gt, $hap4_gt, $hap1_nuc, $hap2_nuc, $hap3_nuc, $hap4_nuc, scalar(keys %haplotypes));

        ## test if the number of unique haplotypes is 2 or 4
        if ( scalar(keys %haplotypes) == 4 ) {
          ## recombination!
          print $OUT "\tR\n";
          print $HAP "\tR\n";
          $recombination_events++;
          push (@RECOMB_indices, $i); ## array of recombination event indices
        } elsif ( scalar(keys %haplotypes) == 2 ) {
          print $OUT "\t.\n";
          print $HAP "\t.\n";
        } else {
          print STDERR "[INFO] SNVs at positions $chrom:$block:$prev_position,$chrom:$block:$curr_position found to have intermediate number of haplotypes\n";
        }

      } else {
        print $HAP "\t.\n";
      }

      ## construct haplotype genotype arrays
      push (@hap1_geno_array, $curr_gt1[0]);
      push (@hap2_geno_array, $curr_gt1[1]);
      push (@hap3_geno_array, $curr_gt2[0]);
      push (@hap4_geno_array, $curr_gt2[1]);
      ## construct haplotype nucleotide arrays
      $curr_gt1[0] == 0 ? push (@hap1_nuc_array, $positions{$curr_position}{REF}) : push (@hap1_nuc_array, $positions{$curr_position}{ALT});
      $curr_gt1[1] == 0 ? push (@hap2_nuc_array, $positions{$curr_position}{REF}) : push (@hap2_nuc_array, $positions{$curr_position}{ALT});
      $curr_gt2[0] == 0 ? push (@hap3_nuc_array, $positions{$curr_position}{REF}) : push (@hap3_nuc_array, $positions{$curr_position}{ALT});
      $curr_gt2[1] == 0 ? push (@hap4_nuc_array, $positions{$curr_position}{REF}) : push (@hap4_nuc_array, $positions{$curr_position}{ALT});

      ## store GTs for next iteration
      @prev_gt1 = @curr_gt1;
      @prev_gt2 = @curr_gt2;
      @prev_nuc1 = @curr_nuc1;
      @prev_nuc2 = @curr_nuc2;
      $prev_position = $curr_position;

      ## record positions
      push (@ALL_positions, $curr_position);
      $i++;
    }

    ## test if genotypes downstream of inferred recombination support the switch
    ## LOGIC:$hap1_geno_up should match either $hap3_geno_up or $hap4_geno_up;
    ## then, if the switch is supported by downstream sites, $hap1_geno_down
    ## should match the alternative hap to whichever $hap1_geno_up matched to
    my @events = (0,@RECOMB_indices,$#hap1_geno_array+1); ## @events is an array with the indices: block start, any recombination events, block end (+1 needed to prevent index out of bounds for end of block)
    for my $i (1..$#events-1) {
      print $TRACTS join ("\t", "$chrom:$block", scalar(@RECOMB_indices), $ALL_positions[$events[$i]], $ALL_positions[$events[$i-1]], $ALL_positions[$events[$i+1]-1]) . "\t";
      print STDERR "$chrom:$block\tTaking window around event at index \#$events[$i] (zero-based), genomic positions ".($ALL_positions[$events[$i-1]])."<-[".($ALL_positions[$events[$i]])."]->".($ALL_positions[$events[$i+1]-1]).")\n" if $verbose;
      ## split haps into upstream and downstream segments
      my $hap1_geno_up = join ('', @hap1_geno_array[$events[$i-1]..$events[$i]-1]); ## upstream = slice genotype array from previous event, to index prior to current event
      my $hap1_geno_down = join ('', @hap1_geno_array[$events[$i]..$events[$i+1]-1]); ## downstream = slice genotype array from current event, to index prior to next event
      my $hap2_geno_up = join ('', @hap2_geno_array[$events[$i-1]..$events[$i]-1]);
      my $hap2_geno_down = join ('', @hap2_geno_array[$events[$i]..$events[$i+1]-1]);
      my $hap3_geno_up = join ('', @hap3_geno_array[$events[$i-1]..$events[$i]-1]);
      my $hap3_geno_down = join ('', @hap3_geno_array[$events[$i]..$events[$i+1]-1]);
      my $hap4_geno_up = join ('', @hap4_geno_array[$events[$i-1]..$events[$i]-1]);
      my $hap4_geno_down = join ('', @hap4_geno_array[$events[$i]..$events[$i+1]-1]);

      print STDERR "$chrom:$block\tUPSTREAM genomic coordinates: $ALL_positions[$events[$i-1]]<->$ALL_positions[$events[$i]-1]\n" if $verbose;
      print STDERR "$chrom:$block\tDOWNSTREAM genomic coordinates: $ALL_positions[$events[$i]]<->$ALL_positions[$events[$i+1]-1]\n" if $verbose;
      print STDERR "$chrom:$block\t$sample1_id hap1 upstream: $hap1_geno_up\n" if $verbose;
      print STDERR "$chrom:$block\t$sample1_id hap1 downstream: $hap1_geno_down\n" if $verbose;
      print STDERR "$chrom:$block\t$sample1_id hap2 upstream: $hap2_geno_up\n" if $verbose;
      print STDERR "$chrom:$block\t$sample1_id hap2 downstream: $hap2_geno_down\n" if $verbose;
      print STDERR "$chrom:$block\t$sample2_id hap3 upstream: $hap3_geno_up\n" if $verbose;
      print STDERR "$chrom:$block\t$sample2_id hap3 downstream: $hap3_geno_down\n" if $verbose;
      print STDERR "$chrom:$block\t$sample2_id hap4 upstream: $hap4_geno_up\n" if $verbose;
      print STDERR "$chrom:$block\t$sample2_id hap4 downstream: $hap4_geno_down\n" if $verbose;

      if ( $hap1_geno_up eq $hap3_geno_up ) {
        print STDERR "$chrom:$block\t$sample1_id hap1 upstream ($hap1_geno_up) matches $sample2_id hap3 upstream ($hap3_geno_up)\n" if $verbose;
        if ( $hap1_geno_down eq $hap4_geno_down ) {
          print STDERR "$chrom:$block\t$sample1_id hap1 downstream ($hap1_geno_down) matches $sample2_id hap4 downstream ($hap4_geno_down)\n" if $verbose;
          print STDERR "$chrom:$block\tNumber of supporting downstream SNVs = ".length($hap1_geno_down)."\n" if $verbose;
          $supporting_sites{length($hap1_geno_down)}++;
          ## get downstream tract length in num SNVs and base pairs
          my @slice = @ALL_positions[$events[$i]..$events[$i+1]-1];
          my $tract_length_bp = ($slice[$#slice]-$slice[0]);
          my $tract_length_snps = length($hap1_geno_down);
          # print STDERR "$chrom:$block\t@slice\n" if $verbose;
          print STDERR "$chrom:$block\tLength of downstream tract = $tract_length_bp bp\n" if $verbose;
          $sum_of_recomb_tracts += $tract_length_bp;
          print $TRACTS join ("\t", $tract_length_snps, $tract_length_bp)."\n";
          push (@affected_snps_per_block, $tract_length_snps);
          push (@affected_bp_per_block, $tract_length_bp);
        } else {
          print STDERR "$chrom:$block\t$sample1_id hap1 downstream ($hap1_geno_down) does not match $sample2_id hap4 downstream ($hap4_geno_down)\n" if $verbose;
          print STDERR "$chrom:$block\tSwitch is NOT supported at downstream sites\n" if $verbose;
          print $TRACTS join ("\t", "0", "-") . "\n";
        }
      } elsif ( $hap1_geno_up eq $hap4_geno_up ) {
        print STDERR "$chrom:$block\t$sample1_id hap1 upstream ($hap1_geno_up) matches $sample2_id hap4 upstream ($hap4_geno_up)\n" if $verbose;
        if ( $hap1_geno_down eq $hap3_geno_down ) {
          print STDERR "$chrom:$block\t$sample1_id hap1 downstream ($hap1_geno_down) matches $sample2_id hap3 downstream ($hap3_geno_down)\n" if $verbose;
          print STDERR "$chrom:$block\tNumber of supporting downstream SNVs = ".length($hap1_geno_down)."\n" if $verbose;
          $supporting_sites{length($hap1_geno_down)}++;
          ## get downstream tract length in num SNVs and base pairs
          my @slice = @ALL_positions[$events[$i]..$events[$i+1]-1];
          my $tract_length_bp = ($slice[$#slice]-$slice[0]);
          my $tract_length_snps = length($hap1_geno_down);
          # print STDERR "$chrom:$block\t@slice\n" if $verbose;
          print STDERR "$chrom:$block\tLength of downstream tract = $tract_length_bp bp\n" if $verbose;
          $sum_of_recomb_tracts += $tract_length_bp;
          print $TRACTS join ("\t", $tract_length_snps, $tract_length_bp)."\n";
          push (@affected_snps_per_block, $tract_length_snps);
          push (@affected_bp_per_block, $tract_length_bp);
        } else {
          print STDERR "$chrom:$block\t$sample1_id hap1 downstream ($hap1_geno_down) does not match $sample2_id hap3 downstream ($hap3_geno_down)\n" if $verbose;
          print STDERR "$chrom:$block\tSwitch is NOT supported at downstream sites\n" if $verbose;
          print $TRACTS join ("\t", "0", "-") . "\n";
        }
      }
    }

    ## record sum of block lengths
    $sum_of_blocks += ($ALL_positions[$#ALL_positions] - $ALL_positions[0]);

    ## print per-blocks info
    # CHROM:BLOCK BLOCK_COORDS SAMPLE1 SAMPLE2 LENGTH_SNPS LENGTH_BP NUM_RECOMBS
    print $BLOCKS join ("\t", "$chrom:$block","$chrom:$ALL_positions[0]-$ALL_positions[$#ALL_positions]",$sample1_id,$sample2_id,scalar(@ALL_positions),($ALL_positions[$#ALL_positions] - $ALL_positions[0]),scalar(@RECOMB_indices));
    scalar(@affected_snps_per_block) > 0 ? print $BLOCKS "\t".mean(@affected_snps_per_block) : print $BLOCKS "\t0";
    scalar(@affected_bp_per_block) > 0 ? print $BLOCKS "\t".mean(@affected_bp_per_block) : print $BLOCKS "\t0";
    print $BLOCKS "\n";

    ## print haplotype blocks as fasta
    print $FASTA ">$sample1_id:$chrom:$block:1\n";
    print $FASTA join ("", @hap1_nuc_array) . "\n";
    print $FASTA ">$sample1_id:$chrom:$block:2\n";
    print $FASTA join ("", @hap2_nuc_array) . "\n";
    print $FASTA ">$sample2_id:$chrom:$block:1\n";
    print $FASTA join ("", @hap3_nuc_array) . "\n";
    print $FASTA ">$sample2_id:$chrom:$block:2\n";
    print $FASTA join ("", @hap4_nuc_array) . "\n";
    print $FASTA "\n";
  }
}
close $OUT;
close $FASTA;
close $BLOCKS;
close $HAP;
close $TRACTS;

## print histogram of supporting SNVs
open (my $HIST, ">".$out_prefix.".hist.txt") or die $!;
my ($singleton_sites, $more_than_singleton_sites);
foreach (sort {$a<=>$b} keys %supporting_sites) {
  print $HIST "$_\t$supporting_sites{$_}\n";
  if ($_ == 1) {
    $singleton_sites = $supporting_sites{$_};
  } else {
    $more_than_singleton_sites += $supporting_sites{$_};
  }
}
close $HIST;

print STDERR "[INFO] Number of inferred recombination events = ".commify($recombination_events)."\n";
print STDERR "[INFO] Number of inferred recombination events with support from >1 downstream SNV = ".commify($more_than_singleton_sites)." (".commify(percentage($more_than_singleton_sites, $recombination_events))."\%)\n";
print STDERR "[INFO] Number of inferred recombination events with support from only one other SNV ('singleton' events) = ".commify($singleton_sites)." (".commify(percentage($singleton_sites, $recombination_events))."\%)\n";
print STDERR "[INFO] Total block length = ".commify($sum_of_blocks)." bp\n";
print STDERR "[INFO] Sum of recombinant tracts = ".commify($sum_of_recomb_tracts)." bp (".commify(percentage($sum_of_recomb_tracts, $sum_of_blocks))."\%)\n";
print STDERR "[INFO] Finished ".`date`;

## print to summary file
open (my $SUMMARY, ">".$out_prefix.".summary.txt") or die $!;
print $SUMMARY "Input VCF = $vcf_file\n";
print $SUMMARY "Sample IDs = $sample1_id, $sample2_id\n";
print $SUMMARY "Number of inferred recombination events = ".commify($recombination_events)."\n";
print $SUMMARY "Number of inferred recombination events with support from >1 downstream SNV = ".commify($more_than_singleton_sites)." (".commify(percentage($more_than_singleton_sites, $recombination_events))."\%)\n";
print $SUMMARY "Number of inferred recombination events with support from only one other SNV ('singleton' events) = ".commify($singleton_sites)." (".commify(percentage($singleton_sites, $recombination_events))."\%)\n";
print $SUMMARY "Total block length = ".commify($sum_of_blocks)." bp\n";
print $SUMMARY "Sum of recombinant tracts = ".commify($sum_of_recomb_tracts)." bp (".commify(percentage($sum_of_recomb_tracts, $sum_of_blocks))."\%)\n";
close $SUMMARY;

######### SUBS

sub commify {
    my $text = reverse $_[0];
    $text =~ s/(\d\d\d)(?=\d)(?!\d*\.)/$1,/g;
    return scalar reverse $text;
}

sub percentage {
    my $numerator = $_[0];
    my $denominator = $_[1];
    my $places = "\%.1f"; ## default is one decimal places
    if (exists $_[2]){$places = "\%.".$_[2]."f";};
    my $float = (($numerator / $denominator)*100);
    my $rounded = sprintf("$places",$float);
    return $rounded;
}

sub mean {
    return sum(@_)/@_;
}

__END__
