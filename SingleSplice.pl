#!/usr/bin/perl -w

######################################################################
# Program: SingleSplice.pl
#
# Author: Joshua Welch
#
# Date: 10/27/2014
#
# Description: Identifies differential isoform usage in single cell
# RNA-seq data.
######################################################################

use strict;
use Getopt::Long;
use Statistics::Basic qw(:all);
use IO::Handle qw( );
use Data::Dumper;

sub fitNoiseModel
{
    my $infile = $_[0];
    my $model = `Rscript fitNoiseModel.R $infile 1000`;
    chomp $model;
    print "Result of fitNoiseModel: $model\n";
    $model =~ /\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\-\d+\.\d+)/;
    my $a0 = $1;
    my $a1 = $2;
    my $b0 = $3;
    my $b1 = $4;
    return ($a0,$a1,$b0,$b1);
}

sub getReadDepth
{
    my $infile = $_[0];
    my @total_mapped;
    my @cell_scale_factors;

    #Lines 2, 3, and 9
    open(IN, $infile);
    my $line = <IN>;
    $line = <IN>;
    chomp $line;
    my @ERCC_mapped = split(/,/,$line);
    shift(@ERCC_mapped);
    $line = <IN>;
    chomp $line;
    my @cell_mapped = split(/,/,$line);
    shift(@cell_mapped);
    $line = <IN>;
    $line = <IN>;
    $line = <IN>;
    $line = <IN>;
    $line = <IN>;
    $line = <IN>;
    chomp $line;
    @cell_scale_factors = split(/,/,$line);
    shift(@cell_scale_factors);
    close IN;
    for (my $i = 0; $i < @ERCC_mapped; ++$i)
    {
	$total_mapped[$i] = $ERCC_mapped[$i]+$cell_mapped[$i]; 
    }
    return (\@total_mapped,\@cell_scale_factors);
}

sub testRatioChange
{
    my ($expr_file,$a0,$a1tilde,$b0,$b1,$num_perms) = @_;
    my $test_result = `./test_ratio_change $expr_file $a0 $a1tilde $b0 $b1 $num_perms`;
    my @results = split(/\s+/,$test_result);
    return (\@results);
}

sub testGroupChange
{
    my ($expr_file,$group_file) = @_;
    my $p_val = `Rscript kruskal-wallis.R $expr_file $group_file`;
}

sub printPaths
{
    my $path1 = $_[0];
    my $path2 = $_[1];
    my $file = $_[2];

    open (OUTF, ">$file");
    for (my $i = 0; $i < @{$path1}-1; $i++)
    {
	my $rounded = sprintf("%.4f",$path1->[$i]);
	print OUTF $rounded . ",";
    }
    my $rounded = sprintf("%.4f",$path1->[@{$path1}-1]);
    print OUTF "$rounded,\n";

    for (my $i = 0; $i < @{$path2}-1; $i++)
    {
	my $rounded = sprintf("%.4f",$path2->[$i]);
	print OUTF $rounded . ",";
    }
    $rounded = sprintf("%.4f",$path2->[@{$path2}-1]);
    print OUTF "$rounded,\n";
    close OUTF;
}

my ($asm_file, $num_perms, $num_samples, $tech_spike_ins, $read_depth_file, $group_file);
my $help = '';
if ( !GetOptions( 'a=s' => \$asm_file,
                  'p=s' => \$num_perms,
                  's=s' => \$num_samples,
                  't=s' => \$tech_spike_ins,
		  'r=s' => \$read_depth_file,
		  'g=s' => \$group_file
) )
{
  print <<END_USAGE;

Usage:
    -a  DiffSplice output file containing ASM abundance estimates
    -p  Number of permutations
    -s  Number of samples
    -t  Technical spike-ins file
    -r  File containing counts of total cellular and spike-in reads
    -g  CSV file listing a group for each sample
    
END_USAGE
  exit(1);
}

my ($a0,$a1tilde,$b0,$b1) = fitNoiseModel($tech_spike_ins);
my ($total_mapped,$cell_scale_factors) = getReadDepth($read_depth_file);

my $threshold = 10;
my $top_k_paths = 10;
my $sig = 0.05;

my $paths_read = 0;
my $paths_found = 0;
open(IN, $asm_file);
open(OUT, ">path_coverage.csv");
open(OUT2, ">ratio_tests.csv");
open(SOUT, ">significant_paths.csv");
print OUT2 "ASM,Path 1,Path 2,Mean 1,Mean 2,Expected Variance,Observed Variance,Variance P-value,Group Change P-value\n";
my $asm_counter = 0;
while (my $line = <IN>) {
    my @cols = split(/\s+/, $line);
    my $asm_id = $cols[0];
    my $asm_path_id = $cols[1];
    my $coords = $cols[2];
    $coords =~ /.*(chr\w+)\:(\d+)\-(\d+):[\-\+].*/;
    my $chr = $1;
    my $start = $2;
    my $end = $3;
    my $event_type = $cols[3];
    $asm_path_id =~ /.*\.p(\d+)/;
    my $num_paths = $1;

    ++$asm_counter;
    #Parse DiffSplice output to get
    #ASM path counts and ASM path fractions
    my @paths_above_thresh = ();
    my %means;
    
    my $pos_samples = 0;
    ++$paths_read;
    my @path_counts;
    for (my $j = 4; $j < $num_samples+4; $j++)
    {

	#print OUT1 $cols[$j] .",";
	my $coverage = $cols[$j];
	my $cell_num = $j-4;
	my $rpkm_size_norm = $coverage/$total_mapped->[$cell_num] * $cell_scale_factors->[$cell_num]*1000000000;
	if ($rpkm_size_norm >= $threshold)
	{
	    ++$pos_samples;
	}
	push(@path_counts,$rpkm_size_norm);
    }
    my $mean = mean(@path_counts);
    $mean =~ s/,//g; #remove annoying commas in numbers
    if ($pos_samples >= 10) {
        push(@paths_above_thresh,\@path_counts);
	$means{$mean} = scalar(@paths_above_thresh)-1;
    }

    for (my $i = 1; $i < $num_paths; $i++)
    {
	$pos_samples = 0;
	++$paths_read;
	if ($paths_read % 10 == 0)
	{
	    print "\r$paths_read ASM paths processed...";
	    my $fh = select(STDOUT);
	    $fh->flush();
	}
	$line = <IN>;
	chomp $line;
	@cols = split(/\s+/, $line);
	my @path_counts;
	for (my $j = 4; $j < $num_samples+4; ++$j)
	{
	    #print OUT1 $cols[$j] .",";
	    my $coverage = $cols[$j];
	    my $cell_num = $j-4;
	    my $rpkm_size_norm = $coverage/$total_mapped->[$cell_num] * $cell_scale_factors->[$cell_num] * 1000000000;
	    if ($rpkm_size_norm >= $threshold)
	    {
		++$pos_samples;
	    }
	    push(@path_counts,$rpkm_size_norm);
	}
	my $mean = mean(@path_counts);
	$mean =~ s/,//g; #remove annoying commas in numbers
	if ($mean >= $threshold && $pos_samples >= 10) {
	    push(@paths_above_thresh,\@path_counts);
	    $means{$mean} = scalar(@paths_above_thresh)-1;
	}
	#print OUT1 "\n";
	#print OUT2 "\n";
    }
    if (@paths_above_thresh > 1) {
      $paths_found += @paths_above_thresh;
      my @sorted_means = sort {$a <=> $b} keys %means;
      for (my $i = 0; $i < @paths_above_thresh; $i++)
      {
        print OUT "$asm_id,$chr,$start,$end,$event_type,";
        for (my $j = 0; $j < @{$paths_above_thresh[$i]} - 1; $j++)
        {
	  my $rounded = sprintf("%.4f",$paths_above_thresh[$i][$j]);
          print OUT $rounded . ",";
        }
	my $rounded = sprintf("%.4f",$paths_above_thresh[$i][@{$paths_above_thresh[$i]}-1]);
        print OUT "$rounded\n";
      }
      my %sig_inds;
      #Test ratios for each pair of ASM paths
      for (my $i = 0; $i < $top_k_paths && $i < @paths_above_thresh; $i++)
      {
	  my $ind1 = $means{$sorted_means[$i]};
	  for (my $j = 0; $j < $top_k_paths && $j < @paths_above_thresh; $j++)
	  {
	      if ($i < $j)
	      {
		  my $ind2 = $means{$sorted_means[$j]};
		  if (!defined $ind2)
		  {
		      print Dumper(%means);
		      print Dumper(@sorted_means);
		  }
		  printPaths($paths_above_thresh[$ind1],$paths_above_thresh[$ind2],"expr.csv");
		  my $test_result = testRatioChange("expr.csv",$a0,$a1tilde,$b0,$b1,$num_perms);
		  my $group_change = testGroupChange("expr.csv",$group_file);
		  print OUT2 "$asm_id,$ind1,$ind2," 
		      . $test_result->[0] .","
		      . $test_result->[1] .","
		      . $test_result->[2] .","
		      . $test_result->[3] .","
		      . $test_result->[4] .","
		      . $group_change . "\n";
		  if ($group_change <= $sig)
		  {
		      $sig_inds{$ind1} = 1;
		      $sig_inds{$ind2} = 1;
		      $asm_path_id =~ /(.*\.p)\d+/;
		      my $asm_path_minus_number = $1;
		      print SOUT "$asm_path_minus_number" . ($num_paths-$ind1) . ",$asm_path_minus_number" . ($num_paths-$ind2) . ",$event_type,";
		      for (my $j = 0; $j < @{$paths_above_thresh[$ind1]} - 1; $j++)
		      {
			  my $rounded = 0.5;
			  if ($paths_above_thresh[$ind2][$j] > 0)
			  {
			      $rounded = sprintf("%.4f",$paths_above_thresh[$ind1][$j]/($paths_above_thresh[$ind1][$j]+$paths_above_thresh[$ind2][$j]));
			  }
			  print SOUT $rounded . ",";
		      }
		      my $rounded = 0.5;
		      if ($paths_above_thresh[$ind2][@{$paths_above_thresh[$ind2]}-1] > 0)
		      {
			  $rounded = sprintf("%.4f",
				  $paths_above_thresh[$ind1][@{$paths_above_thresh[$ind1]}-1]/
				  ($paths_above_thresh[$ind1][@{$paths_above_thresh[$ind1]}-1]+$paths_above_thresh[$ind2][@{$paths_above_thresh[$ind2]}-1]));
		      }
		      print SOUT "$rounded\n";
		  }
	      }
	  }
      }
  }  
}

print "\n$paths_found ASM paths found after filtering for expression level.\n";
close OUT;
close OUT2;
close SOUT;

system("rm expr.csv");
system("rm path_coverage.csv");
system("rm significant_paths.csv");
