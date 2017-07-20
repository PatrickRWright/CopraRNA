#!/usr/bin/perl -s

$tab1 = $ARGV[0];	# reference
$tab2 = $ARGV[1];	# target
$domcheck = 1;
$MAX_FPOS = 500;

if (! $OVLPRATIO1 && ! $OVLPRATIO2) {
	$OVLPRATIO1 = 0.5;
}

if ($self && ! $tab2) {
	$tab2 = $tab1;
}
if ($xref1 && -d $xref1) {
	$xrefdir = $xref1;
} elsif ($xref2 && -d $xref2) {
	$xrefdir = $xref2;
}
if ($xrefdir) {
	opendir(D, $xrefdir);
	while ($f = readdir(D)) {
#		($sp) = ($f =~ /^([^\.]+)\./);
		next if ($f =~ /^\./);
		$sp = $f;
		if ($sp) {
			next if (! open(F, "$xrefdir/$f"));
			while(<F>){
				($n1,$n2) = split;
				$n1 =~ tr/a-z/A-Z/;
				$n2 =~ tr/a-z/A-Z/;
				if ($xref1) {
					## Convert name1 into name2
					$XRef{$sp}->{$n1}->{$n2} = 1;
				} elsif ($xref2) {
					## Convert name1(=n2) into name1(=n1)
					$XRef{$sp}->{$n2}->{$n1} = 1;
				}
			}
		}
	}
}
if ($refspec & 1) {
	open(T, $tab1) || die "Can't open $tab1\n";
	while(<T>){
		if (/^([a-zA-Z]+):/){
			$spn1{$1} = 1;
		}
	}
	close(T);
}
if ($refspec & 2) {
	open(T, $tab2) || die "Can't open $tab2\n";
	while(<T>){
		if (/^([a-zA-Z]+):/){
			$spn2{$1} = 1;
		}
	}
	close(T);
}
if ($genetab) {
	open(G, $genetab) || die;
	while(<G>){
		($sp,$name,$aalen) = split;
		$Length{"$sp:$name"} = $aalen;
	}
	close(G);
}

$clnum = 0;
open(T2, $tab2) || die "Can't open $tab2\n";
while (<T2>) {
	chop;
	if(/^# SPEC=(.*)$/) {
		@spec2 = split(/,/, $1);
		$i = 0;
		foreach $sp (@spec2) {
			$spn2{$sp} = ++$i;
		}
	} elsif (/^Cluster\s*(\S+)/) {
		$clname = $1;
		$clnum++;
		$ClName[$clnum] = $clname;
print STDERR "." if ($clnum % 100 == 0);
		$Members = [];
		$split_cl = 0;
	} elsif (/^\s*$/) {
		if ($maxsize < @{$Members}) {
			$maxsize = @{$Members};
		}
		$Clust[$clnum] = $Members;
		$SplitNum[$clnum] = $cl_splitnum;
		$split_clnum++ if ($cl_splitnum);
		$cl_splitnum = 0;
	} else {
		my($name,$from,$to) = split;
		my($sp, $gname, $dom) = &convname($name);
		if (! $from && ! $to) {
			if (defined $Length{$name}) {
				$from = 1;
				$to = $Length{$name};
			} else {
			# segment data undefined
				if ($domcheck) {
					$domcheck = 0;
					print STDERR "turn off domcheck [$_]\n";
				}
			}
		}
		my $domdata = {name=>$name,dom=>$dom,from=>$from,to=>$to};
		$uname = &genename($sp,$gname,'u');
		push(@{$Dom{$uname,$clnum}}, $domdata);
		push(@{$ClNum{$uname}}, $clnum);
##		push(@{$Members},$name);
		push(@{$Members},$domdata);
		++$total_domnum;
		++$total_seqnum if ($dom == 1 || $dom eq '');
		++$total_split if ($dom == 1);
		$cl_splitnum++ if ($dom);
	}
}
print STDERR "Done\n";
close(T2);
$ClNum1 = $clnum;
$split_clnum++ if ($split_cl);

if (! $diff) {
	print "##reference cluster: $tab1\n";
	print "##target cluster: $tab2\n";
	print "##target clustnum: $clnum\n";
	print "##number of sequences: $total_seqnum\n";
	print "##number of domains: $total_domnum\n";
	print "##split sequences: $total_split\n";
	print "##split clusters: $split_clnum\n";
	print "##max cluster size: $maxsize\n";
	print "##\n";
}

$curr_clnum = 0;
open(T1, $tab1) || die "Can't open $tab1\n";
while (<T1>) {
	if(/^# SPEC=(.*)$/) {
		@spec1 = split(/,/, $1);
		$i = 0;
		foreach $sp (@spec1) {
			$spn1{$sp} = ++$i;
		}
	} elsif (/^Cluster\s*(\S+)/) {
		$clname = $1;
		$cnt = 0;
		$curr_clnum++;
		$domcnt = 0;
		@Genes = ();
	} elsif(/^\S+/) {
		my($name,$from,$to) = split;
		if (! $from && ! $to) {
			if (defined $Length{$name}) {
				$from = 1;
				$to = $Length{$name};
			} else {
			# segment data undefined
				if ($domcheck) {
					$domcheck = 0;
					print STDERR "turn off domcheck\n";
				}
			}
		}
		if ($name =~ /:/) {
			($sp,$name) = split(/:/,$name);
		}
		if ($name =~ /\(([0-9]+)\)$/) {
			$domnum = $1;
			$name =~ s/\([0-9]+\)$//;
			$domcnt++;
		}
##		$uname = uc($name);
		$uname = &genename($sp,$name,'u');
		push(@Genes,{
			sp=>$sp,name=>$name,uname=>$uname,dom=>$domnum,
			from=>$from,to=>$to
		});
		$cnt++;
	} else {
		undef %Cnt;
		undef %Jaccard;
		undef %Conv;
		undef %FoundGn;
		undef %FoundCl;
		undef %NotFound;
		$diffcnt = 0;
		## count intersection
#$dommatch = $domcntR = $domcntT = 0;
		foreach $gd (@Genes) {
			foreach $clnum (&getClnumList($gd->{sp},$gd->{name})) {
				foreach $d (@{$Dom{$gd->{uname},$clnum}}) {
				    if (&overlap($gd, $d)) {
#$dommatch++ if ( $gd->{dom} && $d->{dom} );
#$domcntR++ if ( $gd->{dom} ); #ref
#$domcntT++ if ( $d->{dom} ); #target
					$Cnt{$clnum}++;
					last;
				    }
				}
			}
		}
		## calculate Jaccard
		$gnum1 = @Genes;
		foreach $clnum (keys %Cnt) {
			$gnum2 = @{$Clust[$clnum]};
			if ($gnum1+$gnum2-$Cnt{$clnum}==0){
#				print "$gnum1,$gnum2,$Cnt{$clnum}\n";
			} else {
			$Jaccard{$clnum} =
				$Cnt{$clnum} / ($gnum1+$gnum2-$Cnt{$clnum});
			}
		}
		$j = 1;
		@d = sort {$Jaccard{$b}<=>$Jaccard{$a}}(keys %Jaccard);
#		if ($self) {
#			foreach $d (@d) {
#				if ($d == $clname) {
#				}
#			}
#			$truetop = shift @d;
#		}
		## the best compatible group
		$maxclnum = $d[0];
		if ($self && $maxclnum == $curr_clnum) {
			## for self match, eliminate identical group
			shift @d;
			$maxclnum = $d[0];
			next if (! $maxclnum);
		}
		## assign group id
		foreach $k (@d){
			$Conv{$k} = $j++;
		}
		foreach $gd (@Genes) {
			$min = 10000; $gene_clnum = 0;
			next if ( (defined %spn1 && ! $spn1{$gd->{sp}})
				 || (defined %spn2 && ! $spn2{$gd->{sp}}) );
			foreach $clnum (&getClnumList($gd->{sp},$gd->{name})) {
				next if ($self && $clnum==$curr_clnum);
				next if (! defined $Conv{$clnum});
				my $flag = 0;
				foreach $dm (@{$Dom{$gd->{uname},$clnum}}) {
					if (&overlap($gd,  $dm)) {
						$flag = 1; last;
					}
				}
				next if (! $flag);
				if ($min > $Conv{$clnum}) {
					$min = $Conv{$clnum};
					$gene_clnum = $clnum;
				}
			}
			$diffcnt++ if ($min != 1);
			if ($min != 10000) {
				$FoundCl{$gene_clnum} = 1;
				$gd->{clnum} = "$min";
				$gd->{orig_clnum} = "$gene_clnum";
			}
			$gene = &genename($gd->{sp}, $gd->{name}, 'u');
#                        if ($domnum = $Dom{$gene,$gene_clnum}->{dom}) {
#				$gd->{dom2} = $domnum;
#                        }
			$FoundGn{$gene} = 1;
			if ($xref1 || $xref2) {
				foreach $g (&getXRefList($gd->{sp},$gd->{name})) {
					$FoundGn{"$gd->{sp}:$g"} = 1;
				}
			}
		}
		foreach $clnum (keys %FoundCl) {
			$min = $Conv{$clnum};
			foreach $gn (@{$Clust[$clnum]}) {
				($sp,$gname,$dom) = &convname($gn->{name});
#				$gene = "$sp:$gname";
				$gene = &genename($sp,$gname,'u');
				next if ( (defined %spn1 && ! $spn1{$sp})
					 || (defined %spn2 && ! $spn2{$sp}) );
				if ($FoundGn{$gene}) {
				} else {
					$diffcnt++ if ($min==1);
					$NotFound{$gene} = {clnum=>$min};
##                                	if ($domnum = $Dom{$gene,$clnum}->{dom}) {
##						$NotFound{$gene}->{dom}=$domnum
##                               	}
				}
			}
		}
		next if ($MinSize && @Genes < $MinSize);
		$total_clnum++;
		next if ($diff && $diffcnt==0);
		$clustdiff++ if ($diff);
		&tableOut($clname, $ClName[$maxclnum], \@Genes, \%NotFound);
	}
}
close(T1);
if ($diff && $clustdiff) {
	print "Number of differences: $clustdiff/$total_clnum\n";
	exit(1);
}
print "test OK\n" if ($verbose);
exit(0);

sub tableOut {
	my($clname, $clname2, $Genes, $NotFound, $splitnum) = @_;
	my($setFlag);
	my($flag);
	print "#title: $Title\n" if ($Title);
	print "#subtitle: $SubTitle\n" if ($SubTitle);

	if (@spec) {
		print "#spec: ";
		for ($j = 0; $j < @spec; $j++) {
#			next if ($statout && $skipspec{$spec[$j]});
			next if (! $spn1{$spec[$j]} || ! $spn2{$spec[$j]});
			print "," if ($flag);
			print "$spec[$j]";
			$flag = 1;
		}
		print "\n";
	}
	print "#cluster: $clname\n";
	print "#maxmatch: $clname2\n";
	print "#split: $splitnum\n" if ($splitnum);
	foreach $gd (@{$Genes}) {
		next if ( (defined %spn1 && ! $spn1{$gd->{sp}})
			 || (defined %spn2 && ! $spn2{$gd->{sp}}) );
		if ($gd->{clnum} == 1) {
			$status = '*';
			foreach $dm (@{$Dom{$gd->{uname},$gd->{orig_clnum}}}){
				$ovlp = &overlap($gd, $dm, 1);
				last if ($ovlp >= 0.5);
			}
			if ($ovlp < 0.5) {
				$status = 'x';
			}
		} else {
			$status = '-';
			$gd->{clnum} = 0 if (! $gd->{clnum});
		}
		print "$gd->{sp}:$gd->{name} $status $gd->{clnum}\n";
	}
	my @NotFound =  keys %{$NotFound};
	if (@NotFound > $MAX_FPOS) {
	    print "#false_pos: ", 0+@NotFound, "\n";
	} else {
	    foreach $gene (@NotFound) {
		next if (! $gene);
		$clnum = $NotFound{$gene}->{clnum};
		$status = '+';
		if ($fullout || $clnum == 1) {
			print "$gene $status $clnum\n";
		}
 	    }
	}
	print "//\n";
}

sub getClnumList {
	my($sp, $gene) = @_;
	my(%tmpClNums);
	$gene = uc($gene);
	if ($xref1 || $xref2) {
		foreach $g (keys %{$XRef{$sp}->{$gene}}){
			foreach $cn (@{$ClNum{"$sp:$g"}}) {
				$tmpClNums{$cn} = 1;
			}
		}
		return keys %tmpClNums;
	} else {
		return @{$ClNum{"$sp:$gene"}};
	}
}
sub getXRefList {
	my($sp, $gene) = @_;
	$gene =~ tr/a-z/A-Z/;
	return keys %{$XRef{$sp}->{$gene}};
}
sub convname {
	my($name) = @_;
	my($sp,$name) = split(/:/, $name);
	my($dom);
	if ($name =~ /\(([0-9]+)\)$/) {
		$dom = $1;
		$name =~ s/\([0-9]+\)$//;
	}
	($sp, $name, $dom);
}
sub genename {
	my($sp,$name,$opt) = @_;
	if ($opt =~ /u/) {
		return $sp . ":" . uc($name);
	} else {
		return $sp . ":" . $name;
	}
}

sub overlap_SIMPLE {
	## obsolete
	my($seg1, $seg2) = @_;
	return 1 if (! $domcheck);
	if ($seg1->{from} <= $seg2->{to} && $seg2->{from} <= $seg1->{to}) {
		return 1;
	} else {
		return 0;
	}
}
sub overlap {
	my($seg1, $seg2, $covcheck) = @_;

	return 1 if (! $domcheck);
	my($maxfrom) = &max($seg1->{from}, $seg2->{from});
	my($minto) = &min($seg1->{to}, $seg2->{to});
	my($ovlplen) = $minto - $maxfrom;
	my($minlen) = &min( &seglen($seg1), &seglen($seg2) );
	my($maxlen) = &max( &seglen($seg1), &seglen($seg2) );


#print "$seg1->{from},$seg1->{to} $seg2->{from},$seg2->{to}; $minto $maxfrom $ovlplen $minlen\n";
	if ($OVLPRATIO1 && ! $OVLPRATIO2) {
		return ($ovlplen > $minlen * $OVLPRATIO1);
	} elsif (! $OVLPRATIO1 && $OVLPRATIO2) {
		return ($ovlplen > $maxlen * $OVLPRATIO2);
	} elsif ( ($ovlplen > $minlen * $OVLPRATIO1) ||
		($ovlplen > $maxlen * $OVLPRATIO2 && $ovlplen > $MINOV) ) {
		if ($covcheck) {
			return $ovlplen / $maxlen;
		} else {
			return 1;
		}
	} else {
		return 0;
	}
}
sub seglen {
	$_[0]->{to} - $_[0]->{from} + 1;
}
sub max {
	$_[0] > $_[1] ? $_[0] : $_[1];
}
sub min {
	$_[0] < $_[1] ? $_[0] : $_[1];
}
