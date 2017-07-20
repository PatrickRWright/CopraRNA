#!/usr/bin/env perl
use strict;
use File::Basename;
use POSIX qw(log10);

# --------------------------------------------------------------------------------------------------
# author:	steffen lott
# mail: 	steffen.lott@uni-freiburg.de
# date: 	16-01-2014
# version: 	1.0
# 
# description:
# 	the tool converts a standart output from CopraRNA into json. the output is stored in the same
# 	directory as the input file
# --------------------------------------------------------------------------------------------------

my @jsonOutArr = ();


# --------------------------------------------------------------------------------------------------
# -- read from command-line --
# CopraRNA output
# --------------------------------------------------------------------------------------------------
my $tableIn = $ARGV[0];





# --------------------------------------------------------------------------------------------------
# main structure
# --------------------------------------------------------------------------------------------------

# read table from file
my @jsonTable = readTable($tableIn);

# write output in json format -> copraRNA.json
writeJson($tableIn, @jsonOutArr);





# --------------------------------------------------------------------------------------------------
# functions 
# --------------------------------------------------------------------------------------------------

# read table
sub readTable{
	my $fileName   = $_[0];
	my $rowNum     = 0;
	my $sNum       = 1;
	my @groupNum   = ();
	my @filePara   = ();
	my @featLabels = ();
	my @rowArr     = ();
	my @jsonOut    = ();
	my $flagF      = 1;		# flagF=1 is needed only once
		
	open(INFILE, "<$fileName")  || die "File not found - \"CopraRNA-File\"!\n";
		while(<INFILE>){
			chomp($_);
			
			if($rowNum > 2){
				@rowArr  = split('\t', $_);
				createJsonEntry($flagF, $sNum, \@jsonOut, \@featLabels, \@filePara, \@rowArr, \@groupNum);
				$flagF   = 0;	
				$sNum++;
			}elsif($rowNum == 0){
				# every group has a specific value
				@groupNum   = split('\t', $_); 
			}elsif($rowNum == 1){
				# first value: max number of groups
				@filePara   = split('\t', $_);
			}elsif($rowNum == 2){
				# header line
				@featLabels = split('\t', $_); 
			}
			
			$rowNum++;		
		}
	close(INFILE);
	
	#return @jsonOut;
}



# write output in json format
sub writeJson{
	my($tableIn,@jsonArr1) = @_ ;
	my $path = dirname($tableIn);
	my $out  = "";
	
	if($path =~ /^\./){
		$out  = $path . "\/" . "copraRNA.json";
	}else{
		$out  = "./" . $path . "/" . "copraRNA.json";
	}
	
	open(FILE , ">$out")  || die "File can't be written - \"copraRNA.json - File\"!\n";
		print FILE "[" . "\n";
		for(my $i = 0 ; $i < @jsonOutArr; $i++){
			if($i == @jsonOutArr - 1){
				print FILE $jsonOutArr[$i] . "\n";
			}else{
				print FILE $jsonOutArr[$i] . ",\n";
			}
			
		}
		print FILE "]" . "\n";
	close(FILE);
}



# create json entry
sub createJsonEntry{
	my($flagF, $sNum,$jsonOutArr, $headerArr, $paraArr, $rowArr, $groupPara) = @_ ;
	my $jsonStr = "";
	my $flagS   = 1;
	my $featNum = 1;
	my $color   = "";
	my $bColor  = "";
	my $groupId = 0;
	my $groupStr= "";
	
	# compute group colors and build string. this string has the same value for every json entry
	for(my $i = 0 ; $i < @$groupPara ; $i++){
		$color      = getColor($paraArr->[1], ($i+1), 0, 1);
		my $rounded = sprintf("%.2f", $groupPara->[$i]);
		$groupStr  .= ($i + 1) . ":" . $rounded . ":" . $color;
		if( ($i+1) < @$groupPara){
			$groupStr .= ";";
		}
	}
	
	
	for(my $i = 11 ; $i < @$rowArr ; $i++){
		
		$jsonStr = "{";
		$jsonStr .= "\"flagF\":" . $flagF   . "," . "\"fMax\":"    .        (@$paraArr - 2)  . ","   ; 
		$jsonStr .=	"\"fNum\":"  . $featNum . "," . "\"feature\":" . "\"" . $headerArr->[$i] . "\"," ; 
		$jsonStr .= "\"flagS\":" . $flagS   . "," . "\"sMax\":"    .        $paraArr->[0]    . ","   ; 
		$jsonStr .= "\"sNum\":"  . $sNum    . "," . "\"source\":"  . "\"" . $rowArr->[4]     . "\"," ;
		
		if($rowArr->[$i] == 1){
			$groupId  = $paraArr->[($i - 9)];
			$color    = getColor($paraArr->[1], $groupId, $rowArr->[1], 0);
		}else{
			$color    = "#ffffff";
		}
		$jsonStr .= "\"color\":" . "\"" . $color . "\"" . "," . "\"value\":" . $rowArr->[1] . ",";
		$jsonStr .= "\"group\":" . "\"" . $groupStr . "\"";
		$jsonStr .= "}";
						
		push(@jsonOutArr, $jsonStr);
		$featNum++;
		$flagS = 0;
	}
}



# getColor computes for a specific groups and a specific p-value a specific color 
sub getColor{
	my $diffGroups = $_[0]; my $groupId = $_[1]; my $pVal = $_[2]; my $sFlag = $_[3];
	my $H = 0; my $S = 0; my $V = 1; my $hexStr = "";
	my @rgb = ();
		
	my $step = 360 / $diffGroups;
	$H = ($groupId * $step);
	
	if($sFlag == 0){
		if($pVal <= (1e-6)){
			$S =  1;
		}else{
			$S = -0.15 * log10($pVal);
		}
	}elsif($sFlag == 1){
		$S = 1;
	}		
	
			
	@rgb  = hsv2rgbTwo($H,$S,$V);
	
	my $hexR = sprintf("%x", $rgb[0]);
	my $hexG = sprintf("%x", $rgb[1]);
	my $hexB = sprintf("%x", $rgb[2]);
		
	if(length($hexR) == 1){
		$hexR = "0" . $hexR;
	}
	if(length($hexG) == 1){
		$hexG = "0" . $hexG;
	}
	if(length($hexB) == 1){
		$hexB = "0" . $hexB;
	}
	
	$hexStr = "#" . $hexR . $hexG . $hexB;
	return $hexStr;
}



# hsv to rgb converter 
sub hsv2rgbOne{
	my ( $h, $s, $v ) = @_;
	my @rgbTmp = ();
	my @rgb    = ();
	
	my $c = $v * $s;
	my $x = $c * (1 - abs( (($h / 60) % 2) - 1 ));
	my $m = $v - $c;
	
	if( ($h >= 0) && ($h < 60) ){
		$rgbTmp[0] = $c;	$rgbTmp[1] = $x;	$rgbTmp[2] = 0;
	}elsif( ($h >= 60) && ($h < 120) ){
		$rgbTmp[0] = $x;	$rgbTmp[1] = $c;	$rgbTmp[2] = 0;
	}elsif( ($h >= 120) && ($h < 180) ){
		$rgbTmp[0] = 0;		$rgbTmp[1] = $c;	$rgbTmp[2] = $x;
	}elsif( ($h >= 180) && ($h < 240) ){
		$rgbTmp[0] = 0;		$rgbTmp[1] = $x;	$rgbTmp[2] = $c;
	}elsif( ($h >= 240) && ($h < 300) ){
		$rgbTmp[0] = $x;	$rgbTmp[1] = 0;		$rgbTmp[2] = $c;
	}elsif( ($h >= 300) && ($h < 360) ){
		$rgbTmp[0] = $c;	$rgbTmp[1] = 0;		$rgbTmp[2] = $x;
	}
	
	$rgb[0] = int( ($rgbTmp[0] + $m) * 255);
	$rgb[1] = int( ($rgbTmp[1] + $m) * 255);
	$rgb[2] = int( ($rgbTmp[2] + $m) * 255);
	
	return @rgb;
}	





sub hsv2rgbTwo{
	my $H = $_[0];	my $S = $_[1];	my $V = $_[2];	
	my @rgb = ();
		
	my $hi = int( $H / 60 );
	my $f  = ( ($H / 60) - $hi);	
	
	my $p = $V * (1 - $S);
	my $q = $V * (1 - $S * $f);
	my $t = $V * (1 - $S * (1 - $f) );
	
	if( ($hi == 0) || ($hi == 6)){
		$rgb[0] = int($V*255); $rgb[1] = int($t*255); $rgb[2] = int($p*255);
	}elsif($hi == 1){
		$rgb[0] = int($q*255); $rgb[1] = int($V*255); $rgb[2] = int($p*255);
	}elsif($hi == 2){
		$rgb[0] = int($p*255); $rgb[1] = int($V*255); $rgb[2] = int($t*255);
	}elsif($hi == 3){
		$rgb[0] = int($p*255); $rgb[1] = int($q*255); $rgb[2] = int($V*255);
	}elsif($hi == 4){
		$rgb[0] = int($t*255); $rgb[1] = int($p*255); $rgb[2] = int($V*255);
	}elsif($hi == 5){
		$rgb[0] = int($V*255); $rgb[1] = int($p*255); $rgb[2] = int($q*255);
	}
	
	@rgb;
}


