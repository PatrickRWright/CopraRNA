#!/usr/bin/env perl
$infile = $ARGV[0];
$outname = $infile if (! $outname);
open(IN, $infile) || die;
open(GENE, ">$outname.gene") || die;
open(TIT, ">$outname.tit") || die;
while(<IN>){
	chomp;
	if (/^>\s*(\S.*)/) {
		($name0,$title0) = split(/\s+/, $1, 2);
		if ($seq) {
			($sp, $gene) = split(/:/, $name);
			$seqlen = length($seq);
			print GENE "$sp $gene $seqlen\n";
			print TIT "$sp:$gene\t$title\n";
		}
		$name = $name0;
		$title = $title0;
		$seq = '';
	} else {
		s/\s//g;
		$seq .= $_;
	}
}
close(IN);
if ($seq) {
	($sp, $gene) = split(/:/, $name);
	$seqlen = length($seq);
	print GENE "$sp $gene $seqlen\n";
	print TIT "$sp:$gene\t$title\n";
}
