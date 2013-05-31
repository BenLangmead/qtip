#!/usr/bin/perl

##
# correct.pl
#
# Simple script that takes a .sam file and spits out one or both of: (a) a
# table with cumulative correct and incorrect alignment totals, stratified by
# mapping quality, and (b) a new .sam file with an additional extra field ZC:i
# that's set to 1 if the alignment is correct, 0 if incorrect.
#

use strict;
use warnings;
use Getopt::Long;

my $format = "art"; # could be "art" or "wgsim"
my $gap = 50;
my $checkChimeras = 1;

GetOptions (
	"format=s" => \$format,
	"gap=i" => \$gap,
	"no-chi" => sub { $checkChimeras = 0; }
) || die "Bad option";

$format eq "art" || $format eq "wgsim" || die "Bad format: '$format'";

scalar(@ARGV) == 3 || die "Expected 3 arguments";
my $sam_fn  = shift @ARGV;
my $roc_fn  = shift @ARGV;
my $samo_fn = shift @ARGV;

my %mapq = ();

##
# Parse ART-style name of simulated read to determine its point of origin.
# The standard format is <chr>:<off>:<strand>
#
sub _parse_art_sim_name($) {
	my @s = split(/:/, $_[0]);
	scalar(@s) >= 3 || die "Bad name: '$_[0]'";
	my $off = $s[-2];
	my $str = $s[-1];
	my $chr = join(":", @s[0..scalar(@s)-3]);
	length($chr) > 0 || die;
	return ($chr, $off, $str);
}

##
# Parse wgsim-style name of simulated read to determine its point of origin.
# The standard format is <chr>:<off>:<strand>
#
sub _parse_wgsim_sim_name($) {
	my $name = shift;
	my $mate = substr($name, -1);
	$name =~ s/\/[12]$//;
	my @s = split(/_/, $name);
	my $ii = pop @s;
	my $edits2 = pop @s;
	my $edits1 = pop @s;
	my $off2 = pop @s;
	my $off1 = pop @s;
	scalar(@s) >= 1 || die "Incorrectly formated read name: $name";
	my $rfname = join("_", @s);
	if($mate eq "1") {
		# I'm mate 1
		return ($rfname, $off1, "?");
	} else {
		# I'm mate 2
		return ($rfname, $off2, "?");
	}
}

sub parse_name($) {
	if($format eq "art") {
		return _parse_art_sim_name($_[0]);
	} else {
		return _parse_wgsim_sim_name($_[0]);
	}
}

open(my $sam_fh,    $sam_fn)   || die "Could not open '$sam_fn' for reading";
open(my $sam_ofh, ">$samo_fn") || die "Could not open '$samo_fn' for writing";
my $min_mapq = 999;
my $max_mapq = 0;
my $nchimeras = 0;
my $next_line = readline $sam_fh;
my $nal_ival = 20000;
my $nal = 0;
while(defined($next_line)) {
	my $cur_line = $next_line;
	$next_line = readline $sam_fh;
	next if $cur_line =~ /^\@/;
	$cur_line =~ /^([^\t]+)\t/;
	my $rdname = $1;
	defined($rdname) || die "Malformed line; no read name at the beginning:\n$cur_line";
	my @lines = ($cur_line);
	if($checkChimeras) {
		while(defined($next_line) && $next_line =~ /^$rdname/) {
			print STDERR "Chimeric alignment in '$sam_fn'; two SAM records for '$rdname'\n";
			push @lines, $next_line;
			$next_line = readline $sam_fh;
		}
	}
	if((++$nal % $nal_ival) == 0) {
		print STDERR "Processed $nal alignments from '$sam_fn' ...\n";
	}
	my $is_correct = 0;
	my $best_mapq = -1;
	scalar(@lines) >= 1 || die;
	$nchimeras++ if scalar(@lines) > 1;
	my $trimmed = 0;
	for my $line (@lines) {
		chomp($line);
		my @ts = split(/\t/, $line);
		my ($rdname, $flags, $cigar, $left, $mapq) = ($ts[0], $ts[1], $ts[5], $ts[3], $ts[4]);
		if(($flags & 4) != 0) {
			print {$sam_ofh} "$line\tZC:i:0\n";
			next;
		}
		# Consider soft and hard trimming, which moves offset up to after the
		# trim point
		my $ltrim = $1 if $cigar =~ (/^(\d+)[SH]/);
		$ltrim = $ltrim || 0;
		$trimmed = 1 if $cigar =~ (/[SH]/);
		$left -= $ltrim;
		$min_mapq = $mapq if $mapq < $min_mapq;
		$max_mapq = $mapq if $mapq > $max_mapq;
		$best_mapq = $mapq if $mapq > $best_mapq;
		my ($chr, $off, $fw) = parse_name($rdname);
		#print STDERR "chr=$chr, ts[2]=$ts[2], off=$off, left=$left\n";
		if($chr eq $ts[2] && abs($off - $left) <= $gap) {
			# Chromosome and position are good
			if($fw ne "?") {
				my $rdfw = (($flags & 16) == 0) ? "+" : "-";
				#print STDERR "chr=$chr, ts[2]=$ts[2], off=$off, left=$left, fw='$fw', rdfw='$rdfw'\n";
				if($fw eq $rdfw) {
					# Orientation is good
					$is_correct = 1;
				}
			} else {
				#print STDERR "chr=$chr, ts[2]=$ts[2], off=$off, left=$left\n";
				$is_correct = 1;
			}
		}
	}
	#print STDERR "is_correct=$is_correct\n";
	$mapq{$is_correct}{$best_mapq}++;
	$mapq{tot}{$best_mapq}++;
	$mapq{nchi}{$best_mapq}++ if scalar(@lines) > 1;
	$mapq{ntri}{$best_mapq}++ if $trimmed;
	my $first = 1;
	for my $line (@lines) {
		chomp($line);
		if($first) {
			print {$sam_ofh} "$line\tZC:i:$is_correct\n";
		} else {
			print {$sam_ofh} "$line\tZC:i:-1\n";
		}
		$first = 0;
	}
}
close($sam_fh);
close($sam_ofh);

# Print the ROC table
open(my $roc_ofh, ">$roc_fn") || die "Could not open '$roc_fn' for writing";
my $cum_correct = 0;
my $cum_incorrect = 0;
my $cum_nchi = 0;
my $cum_ntri = 0;
print {$roc_ofh} "mapq\tcorrect\tcorrect_cum\tincorrect\tincorrect_cum\tfrac_incorrect\tobs_mapq\tchimeras\tchimeras_cum\ttrimmed\ttrimmed_cum\n";
print STDERR "[Min, Max] MAPQ = [$min_mapq, $max_mapq] in '$sam_fn'\n";
for(my $m = $max_mapq; $m >= $min_mapq; $m--) {
	next unless $mapq{tot}{$m};
	my $ncor   = $mapq{1}{$m}    || 0;
	my $nincor = $mapq{0}{$m}    || 0;
	my $nchi   = $mapq{nchi}{$m} || 0;
	my $ntri   = $mapq{ntri}{$m} || 0;
	$cum_correct   += $ncor;
	$cum_incorrect += $nincor;
	$cum_nchi      += $nchi;
	$cum_ntri      += $ntri;
	my $frac  = $nincor / ($ncor + $nincor);
	my $omapq = (($frac > 0) ? (-10 * log($frac)/log(10)) : "Inf");
	print {$roc_ofh} join("\t", ($m, $ncor, $cum_correct, $nincor, $cum_incorrect, $frac, $omapq, $nchi, $cum_nchi, $ntri, $cum_ntri))."\n";
}
close($roc_ofh);
