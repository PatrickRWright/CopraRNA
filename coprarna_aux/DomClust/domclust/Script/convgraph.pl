#!/usr/bin/perl -s

@Color = (black, red, blue, green, magenta, cyan, orange, brown, gray, yellow, pink, darkgreen, navy, gold, salmon, midnightblue, purple, black);

if ($edgecolors) {
	@Color = split(/,/, $edgecolors);
	unshift @Color, 'black';
}

if ($nodecolor && open(F, $nodecolor)){
	while(<F>){
		($sp,$col)= split;
		$NodeColor{$sp} = $col;
	}
	close(F);
}
$leafshape = 'box' if (! $leafshape);

$FontSize = 24;
$FontName = "Helvetica";
$FontColor = "black";
$Size = "3.4,1.6";
$NodeSep= ".12";

print "digraph homgraph {\n";
print " size=\"$Size\" ratio=fill;\n";
print " node [fontsize=$FontSize,fontname=\"$FontName\",fontcolor=\"$FontColor\"]\n";
print "	nodesep=$NodeSep;\n";
if ($nodelabel eq 'none') {
	print "	ranksep=.3;\n";
} elsif ($nodelabel =~ /^sp/) {
	print "	ranksep=.1;\n";
} else {
	print "	ranksep=.7;\n";
}

while (<>) {
	if (/Cluster/) {
		$cln++;
		$flag = 1;
		next;
	} elsif (/^$/) {
		$flag = 0;
		next;
	}
	($id1, $id2, $opt) = split;
	if (! $opt){
		$Node{$id1} = $Node{$id2} = 1;
		$color = $Color[$cln % (0+@Color)];
		print qq{	"$id1" -> "$id2" [color=\"$color\",dir=none,style=bold];\n};
	} elsif ($NoPhysClust) {
		# do nothing
	} elsif ($opt eq 'R') {
		$PhysClust{$id1}->{$id2} = $PhysClust{$id2}->{$id1} = 1;
	} elsif ($opt eq 'L') {
		$PhysClust{$id1}->{$id2} = $PhysClust{$id2}->{$id1} = 1;
	}
}
&cluster;
foreach $n (keys %PhysClustNum) {
	push(@{$Clnum[$PhysClustNum{$n}]}, $n);
}
for ($i = 1; $i < @Clnum; $i++) {
	if (@{$Clnum[$i]} >= 2) {
		print "	subgraph cluster_$i {\n";
		print "		rankdir=BT;\n";
		print "		style=invis;\n";
		print "		rank=same;\n";
	}
	foreach $name (@{$Clnum[$i]}) {
		&print_node($name) if ($Node{$name});
	}
	if (@{$Clnum[$i]} >= 2) {
		print "	}\n";
	}
}
foreach $name (keys %Node) {
	&print_node($name);
}

print "}\n";


sub print_node {
	my($name) = @_;
	my($opt);
	return if ($Node{$name} == 2);
	$Node{$name} = 2;
	if ($name =~ /^[0-9]+/) {
		if ($intlabel) {
			$opt = "shape=circle,height=.2";
		} else {
			$opt = "shape=circle,height=.15,label=\"\"";
		}
	} else {
		my($sp,$name0) = split(/:/, $name);
		if ($nodelabel eq 'none') {
			$opt = "shape=$leafshape,label=\"\",height=.2,width=.2";
		} elsif ($nodelabel eq 'name') {
			$opt = "shape=$leafshape,label=\"$name0\",height=.2,width=.05";
		} elsif ($nodelabel =~ /^sp/) {
			$opt = "shape=$leafshape,label=\"$sp\",height=.2,width=.05";
		} else {
			$opt = "shape=$leafshape";
		}
		if (defined %NodeColor) {
			$opt .= "," if ($opt);
			$opt .= "style=filled,color=\"$NodeColor{$sp}\"";
		}
	}
	print qq{	"$name" [$opt]\n};
}

sub cluster {
	my($cln,$n1);
	foreach $n1 (keys %Node) {
		&cluster0($n1, ++$cln)
	}
}
sub cluster0 {
	my($n1,$cln) = @_;
	my($n2);
	return if ($PhysClustNum{$n1});
	$PhysClustNum{$n1} = $cln;
	foreach $n2 (keys %{$PhysClust{$n1}}) {
		&cluster0($n2,$cln);
	}
}
