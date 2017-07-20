#!/usr/bin/perl -s

open(S, $titfile)  || die;
while (<S>) {
	chomp;
	($name, $descr) = split(/\t/);
	$tit{$name} = $descr;
}
close(S);

while (<>) {
	if (/^Cluster (.*)/) {
		$clustid = $1;
		$clustnum = 0;
		$out = '';
	} elsif (/^$/) {
		print "Cluster $clustid [$clustnum]\n";
		print "$out";
		print "===============\n";
	} elsif ($format eq 'tree') {
		chomp;
		$out .= $_;
		if (/\+\- (\S+)/) {
			$genename = $1;
			$genename =~ s/\(\d+\)//;
			$out .= "  $tit{$genename}";
		}
		$out .= "\n";
	} else {
		chomp;
		($name, $from, $to) = split;
		$genename = $name;
		$genename =~ s/\(\d+\)//;
		$out .= (join("\t", $name, $from, $to, $tit{$genename}) . "\n");
		$clustnum++;
	}
}
