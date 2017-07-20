#!/usr/bin/env perl

$EVAL_CUT = 0.001;
$infile = $ARGV[0];
$distconv = 1; ## edit 2.0.5.1 // this is no longer a parameter // needed to be done like this to change shebang

$BITUNIT = 1/3;    ## the default scoring system is in 1/3 bit units

if (! $skip_sort) {
	## eliminate comment lines (when using blastall -m 9)
	$infile = "egrep -v '^#' $infile | " .
	## excahnge name1 and name2, and also start and end positions
	q{awk ' BEGIN {OFS="\t"}
		$1<=$2 {print}
		$1>$2 {print $2,$1,$3,$4,$5,$6,$9,$10,$7,$8,$11,$12}' | } .
	## sort by name pair followed by E-value
	"sort -k 1,2 -k 11,11g | ";
}

open(IN, $infile) || die;

while (<IN>) {
	next if (/^#/);
	($qid, $sid, $ident, $alilen, $mismatch, $gap,
		$qstart, $qend, $sstart, $send, $eval, $score) = split;
	next if ($eval > $EVAL_CUT);
	if ($qid eq $prev_qid && $sid eq $prev_sid) {
		next;
	}
	$prev_qid = $qid; $prev_sid = $sid;

	if ( $evalcorr) {
		# corrected E-value
		$eval *= len
	}
	$score /= $BITUNIT;
	$dist = 100 - $ident;

        ## edit 2.0.1
        # fix for negative value in log(), low sequence identity causes this issue and 
        # terminates the program 
        $temp1 = $dist / 100;
        $temp2 = 1 - $temp1 - 0.2 * $temp1 * $temp1;
        if ($temp2 <= 0) { next; }

	if ( $distconv && $dist > 0 ) {
	    # Kimura's correction formula for protein sequence distances
		$dist /= 100;
		$dist = - log( 1 - $dist - 0.2 * $dist * $dist);
		$dist *= 100;
	}
	$dist = sprintf("%.0f", $dist);
	$score = sprintf("%.0f", $score);
	print "$qid $sid $qstart $qend $sstart $send $dist $score\n";
}
close(IN);
