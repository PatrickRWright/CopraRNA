#!/usr/bin/env perl

# this file was taken from the domclust package and adapted 

# usage : fasta2genefile.pl <FASTAFILE>
# produces FASTAFILE.gene and FASTAFILE.tit from FASTA header information

$infile = $ARGV[0];
$outname = $infile if (! $outname);

# open streams or die trying
open(IN, $infile) || die("fasta2genefile.pl : cannot open/read from $infile");
open(GENE, ">$outname.gene") || die("fasta2genefile.pl : cannot  write to $outname.gene");
open(TIT, ">$outname.tit") || die("fasta2genefile.pl : cannot  write to $outname.tit");

# process input FASTA file
while(<IN>){
	chomp;
	# get FASTA id
	if (/^>\s*(\S.*)/) {
		# split fasta id in name (up to first whitespace) and rest (=title)
		($name0,$title0) = split(/\s+/, $1, 2);
		# if sequence was already parsed (last record) -> print
		if ($seq) {
			($sp, $gene) = split(/:/, $name);
			$seqlen = length($seq);
			print GENE "$sp $gene $seqlen\n";
			print TIT "$sp:$gene\t$title\n";
		}
		# set data for this (new) record
		$name = $name0;
		$title = $title0;
		# reset sequence
		$seq = '';
	} else {
		# parse sequence (concatenate subsequent lines)
		s/\s//g;
		$seq .= $_;
	}
}
# handle last sequence
if ($seq) {
	($sp, $gene) = split(/:/, $name);
	$seqlen = length($seq);
	print GENE "$sp $gene $seqlen\n";
	print TIT "$sp:$gene\t$title\n";
}
# close streams
close(TIT);
close(GENE);
close(IN);
