#!/usr/bin/perl -w

##
# tabulate_sam.pl
#
# Turn a SAM file into a tab-delimited file appropriate for loading into R
# with read.table(filename, header=T, sep='\t', comment.char="", quote="", stringsAsFactors=F).
#

# Optional fields start at field 12 (1-based)

my $fn = shift @ARGV;
open(SAM, $fn) || die "Could not open '$fn'";
my %flags = ();
while(<SAM>) {
	chomp;
	next if substr($_, 0, 1) eq '@';
	my @ts = split(/\t/, $_, -1);
	for my $i (11..$#ts) {
		my @fl = split(/:/, $ts[$i]);
		scalar(@fl) > 2 || die "Bad flag: $ts[$i]";
		if($fl[1] eq "B") {
			# List of integers
			scalar(@fl) > 2 || die;
			my @cl = split(/,/, $fl[2], -1);
			for(my $j = 1; $j < scalar(@cl); $j++) {
				# Skip first elt, which is type specifier
				$flags{"$fl[0]:$fl[1]:$j"}++;
			}
		} else {
			# Just one value
			$flags{"$fl[0]:$fl[1]"}++;
		}
	}
}
close(SAM);
open(SAM, $fn) || die;
#print join("\t", ("qname", "flag", "rname", "pos", "mapq", "cigar", "rnext",
#                  "pnext", "tlen", "seq", "qual"));
print join("\t", ("qname", "flag", "rname", "pos", "mapq", "cigar", "tlen", "len"));
print "\t".join("\t", sort keys %flags)."\n";
while(<SAM>) {
	chomp;
	next if substr($_, 0, 1) eq '@';
	my @ts = split(/\t/, $_, -1);
	my %myflags = ();
	for my $i (11..$#ts) {
		my @fl = split(/:/, $ts[$i]);
		if($fl[1] eq "B") {
			# List of integers
			scalar(@fl) > 2 || die;
			my @cl = split(/,/, $fl[2], -1);
			for(my $j = 1; $j < scalar(@cl); $j++) {
				# Skip first elt, which is type specifier
				$myflags{"$fl[0]:$fl[1]:$j"} = $cl[$j];
			}
		} else {
			# Just one value
			$myflags{"$fl[0]:$fl[1]"} = join(":", @fl[2..$#fl]);
		}
	}
	print join("\t", ($ts[0], $ts[1], $ts[2], $ts[3], $ts[4], $ts[5], $ts[8], length($ts[9])));
	my @opts = ();
	for my $k (sort keys %flags) {
		push @opts, (defined($myflags{$k}) ? $myflags{$k} : "NA");
	}
	print "\t".join("\t", @opts)."\n";
}
