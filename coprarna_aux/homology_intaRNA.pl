#!/usr/bin/env perl

use strict;
use warnings;
use Cwd 'abs_path'; 
use Parallel::ForkManager;

# get absolute path
my $ABS_PATH = abs_path($0); 
# remove script name at the end
# match all non slash characters at the end of the string
$ABS_PATH =~ s|[^/]+$||g; 
my $PATH_COPRA_SUBSCRIPTS = $ABS_PATH;

# files dedicated to capture output of subcalls for debugging
my $OUT_ERR = "CopraRNA2_subprocess.oe";

my $ncrnas = $ARGV[0]; # input_sRNA.fa
my $upfromstartpos = $ARGV[1]; # 200
my $down = $ARGV[2]; # 100
my $mrnapart = $ARGV[3]; # cds or 5utr or 3utr
my $core_count = $ARGV[4];
my $intarnaParamFile = $ARGV[5];
my $GenBankFiles = "";
my $orgcount = 0;

##########################################################################################
sub getOptionValue
##########################################################################################
# @param 1 option name from CopraRNA_option_file.txt file
# @return the respective value from CopraRNA_option_file.txt
##########################################################################################
{
	my $option = $_[0];
	my $value = `grep -m 1 -P '^\\s*$option:' CopraRNA_option_file.txt | sed 's/^\\s*$option://g'`;
	chomp $value;
	return( $value );
}
######################################################################## END OF SUBROUTINE


my $cores = getOptionValue("core count");

# check for verbose printing
my $verbose = getOptionValue("verbose");

# check for genomePath
my $genomePath = getOptionValue("genomePath");
# sanity check
if ($genomePath eq "") {
	$genomePath = ".";
}

# get amount of top predictions to return
my $topcount = getOptionValue("top count");
$topcount++; # need this to include the header

# check for websrv output printing
my $websrv = getOptionValue("websrv");

# check for enrichment on/off and count
my $enrich = getOptionValue("enrich"); 

# get window size option
my $winsize = getOptionValue("win size");

# get maximum base pair distance 
my $maxbpdist = getOptionValue("max bp dist");


####################################
# open error stream
####################################

open ERRORLOG, ">>err.log" or die("\nError: cannot open file err.log in homology_intaRNA.pl\n\n"); 




##########################################################################################
sub removeInvalidGenomeFiles
##########################################################################################
# @param 0 (optional) error message to be displayed for each removed invalid file
# @return the number of deleted files
##########################################################################################
{
    my $errorMsg = "";
    if (scalar(@_) > 0) {
    	$errorMsg = $_[0];
    }
    
    my $deletedFiles = 0;
    my $dgl = new Parallel::ForkManager($cores);

    # iterate all "gb" files
    foreach my $gbFile (<*.gb.gz>) {
        $dgl->start and next;
    	my $lastTwoLines = `zcat $gbFile | tail -n 2 | tr -d "\n"`;
    
    	# proper gb file ends with "//" line
    	unless ( $lastTwoLines =~ m/.*\/\/\s*$/ ) {
    		`rm -f $gbFile`; # remove link
    		`rm -f $genomePath/$gbFile`; # remove gb file to try download again later
    		# count deletion
    		$deletedFiles++;
    		# print error message if given
    		unless ($errorMsg eq "") {
    			print ERRORLOG "Genome file $gbFile : $errorMsg\n";
    		}
    	} # if invalid file
        $dgl->finish;
    } # for all gbFiles
    
    return( $deletedFiles );
$dgl->wait_all_children;

}
######################################################################## END OF SUBROUTINE





##########################################################################################
sub downloadGenomeAndLink
##########################################################################################
# @param 1 accession number of the genome to download if not already present
##########################################################################################
{
	my $accessionnumber = $_[0];

    my $refseqoutputfile = "$accessionnumber.gb.gz";
    
  	# download genome file if needed
	unless ( -e "$genomePath/$refseqoutputfile" ) {
	  	my $gbCall = $PATH_COPRA_SUBSCRIPTS  . "get_refseq_from_refid.pl -acc $accessionnumber -g $genomePath/$refseqoutputfile";
    	print $gbCall . "\n" if ($verbose); 
      	system $gbCall;
	}
  	# link genome file locally if not present
  	if ( (not -e "$refseqoutputfile") and (-e "$genomePath/$refseqoutputfile") ) {
	        system ("ln -s $genomePath/$refseqoutputfile .");
		# print("ln -s $genomePath/$refseqoutputfile .");
  	}
    
}
######################################################################## END OF SUBROUTINE

# finds replikons genomes 

####################################
print " create kegg 2 refseq mapping\n" if ($verbose);
####################################

my $keggtorefseqnewfile = $PATH_COPRA_SUBSCRIPTS . "kegg2refseqnew.csv";
# RefSeqID -> space separated RefSeqIDs // 'NC_005140' -> 'NC_005139 NC_005140 NC_005128'
my %refseqaffiliations = ();

# read kegg2refseqnew.csv
open(MYDATA, $keggtorefseqnewfile) or die("\nError: cannot open file $keggtorefseqnewfile in homology_intaRNA.pl\n\n");
    my @keggtorefseqnew = <MYDATA>;
close MYDATA;

# get the refseq affiliations
foreach(@keggtorefseqnew) {
    # split off quadruplecode (pseudokegg id)
    my @split = split("\t", $_);
    my $all_refseqs = $split[1];
    chomp $all_refseqs;
    # split up refseq ids
    my @split_refseqs = split(/\s/, $all_refseqs);
    foreach(@split_refseqs) {
        $refseqaffiliations{$_} = $all_refseqs;
        
    }
}

# add "ncRNA_" to fasta headers
system "sed 's/>/>ncRNA_/g' $ncrnas > ncrna.fa"; 

my $no_of_args = scalar(@ARGV); 

# assign correct refseq IDs for each sequence
for (my $i=6;$i<scalar(@ARGV);$i++) {
    # split up the refseq list for one organism
    my @split = split(/\s/, $refseqaffiliations{$ARGV[$i]});
    # get the first id entry
    my $first_refseq_id = $split[0];
    # override in ncrna.fa
    system "sed -i 's/$ARGV[$i]/$first_refseq_id/g' ncrna.fa";
}

# override $ncrnas variable
$ncrnas = "ncrna.fa";

# get Orgcount
$orgcount = (scalar(@ARGV) - 6);


## prepare input for combine_clusters.pl
## Download Refseq files by Refseq ID 
my $RefSeqIDs = `grep ">" input_sRNA.fa | tr '\n' ' ' | sed 's/>//g'`; 
my @split_RefIds = split(/\s+/, $RefSeqIDs);

my $dwn = new Parallel::ForkManager($cores);
foreach(@split_RefIds) {
    my $currRefSeqID = $_;
    my $presplitreplicons = $refseqaffiliations{$currRefSeqID};
    my @replikons = split(/\s/, $presplitreplicons);
  
    foreach my $accessionnumber (@replikons) {
	my $refseqoutputfile = "$accessionnumber.gb.gz";
        $GenBankFiles = $GenBankFiles . $accessionnumber . ".gb.gz" . ",";
        if ( (not -e "$refseqoutputfile") and (not -e "$genomePath/$refseqoutputfile") ){

	$dwn->start and next;
        my $dwnCall = "python ". $PATH_COPRA_SUBSCRIPTS  . "download_genomes.py -g $accessionnumber -p $genomePath/";
        system $dwnCall;
        $dwn->finish;    }
    }
           
	# remove trailing comma
    chop $GenBankFiles;
	# add space separator to begin next block of genome files
    $GenBankFiles = $GenBankFiles . " ";

    
}
$dwn->wait_all_children; 

## RefSeq correct download check for 2nd try
removeInvalidGenomeFiles();


## list of available refseq files
my @gbFiles = ();

my @totalrefseqFiles = split(/\s|,/, $GenBankFiles);
my $consistencyswitch = 1;

my $limitloops = 0;
my $sleeptimer = 30;

while($consistencyswitch) {
    @gbFiles = <*.gb.gz>;
    foreach my $gbFile (@totalrefseqFiles) {
	chomp $gbFile;
        if(grep( /^$gbFile$/, @gbFiles )) { 
            $consistencyswitch = 0;
        } else {
             $limitloops++;
             $consistencyswitch = 1;
             if($limitloops > 1000) { 
                $consistencyswitch = 0;
				die( "Could not download all RefSeq *gb files... Restart your job.\n", 1 ); 
             }
	# extract genome ID from genome file name
             my $accNr = $gbFile;
	     if ($gbFile =~ m/^(.+)\.gb(\.gz)?$/) {
		$accNr = $1;
	     }
 	     # print("downloadGenomeAndLink\n");
	     downloadGenomeAndLink( $accNr );
             last;
        }
	
    }
}


### end availability check


### refseq correct DL check kill job 
if ( 0 < removeInvalidGenomeFiles("Genome file did not download correctly. This is probably due to a connectivity issue with the NCBI servers. Please retry later..") ) {
	die( "Not all RefSeq *gb files downloaded correctly. Restart your job.\n", 1 );
}

###########################################
print "check sanity of CDS features\n" if ($verbose);
###########################################
my $pm = new Parallel::ForkManager($cores);
foreach my $gbFile (@gbFiles) {

	## fixing issue with CONTIG and ORIGIN both in gb file (can't parse without this) 
	# check if gbFile has to be corrected
	my $gbContainsCONTIG = `zgrep -m 1 -c -P "^CONTIG" $gbFile`;
	if ($gbContainsCONTIG > 0) {
		# remove "CONTIG" string
    	system "zcat $gbFile | sed '/^CONTIG/d' | gzip -9 > tmp.gz; mv -f tmp.gz $gbFile"; ## d stands for delete
	}

    $pm->start and next; 
	# check for CDS features
    system $PATH_COPRA_SUBSCRIPTS . "check_for_gene_CDS_features.pl $gbFile >> gene_CDS_exception.txt";
    $pm->finish;
}
$pm->wait_all_children;
{ 
    open(MYDATA, "gene_CDS_exception.txt") or die("\nError: cannot open file gene_CDS_exception.txt at homology_intaRNA.pl\n\n");
        my @exception_lines = <MYDATA>;
    close MYDATA;
    
    if (scalar(@exception_lines) >= 1) {
        my $exceptionRefSeqs = "";
        foreach(@exception_lines) {
            my @split = split(/\s+/,$_);
            $exceptionRefSeqs = $exceptionRefSeqs . $split[-1] . " ";
        }
        print ERRORLOG "Error: gene but no CDS features present in $exceptionRefSeqs.\n This is most likely connected to currently corrupted RefSeq record(s) at the NCBI.\nPlease resubmit your job without the currently errorous organism(s) or wait some time with your resubmission.\nUsually the files are fixed within ~1 week.\n"; 
    }
}
## end CDS gene exception check


unless (-e "cluster.tab") { # only do if cluster.tab has not been imported
###########################################
print "get cluster.tab with DomClust\n" if ($verbose);
###########################################

    ### get AA fasta for homolog clustering

    foreach my $gbFile (@gbFiles) {
        system $PATH_COPRA_SUBSCRIPTS . "get_CDS_from_gbk.pl $gbFile >> all.fas"; 
    }

    # prep for DomClust
	print "formatdb for all.fas\n" if ($verbose);
    system "formatdb -i all.fas" unless (-e "all.fas.blast.gz"); 
    # blast sequences
	print "blastall for all.fas\n" if ($verbose);
    system "blastall -a $cores -p blastp -d all.fas -e 0.001 -i all.fas -Y 1e9 -v 30000 -b 30000 -m 8 2>> $OUT_ERR | gzip -9 > all.fas.blast.gz " unless (-e "all.fas.blast.gz"); # change the -a parameter to qdjust core usage 
    # remove empty error file
    system $PATH_COPRA_SUBSCRIPTS . "blast2homfile.pl all.fas.blast.gz > all.fas.hom"; 
    system $PATH_COPRA_SUBSCRIPTS . "fasta2genefile.pl all.fas";
    # DomClust
	print "domclust for all.fas.*\n" if ($verbose);
    my $domclustExitStatus = system "domclust all.fas.hom all.fas.gene -HO -S -c60 -p0.5 -V0.6 -C80 -o5 > cluster.tab 2>> ".$OUT_ERR;
    $domclustExitStatus /= 256; # get original exit value
    # ensure domclust went fine
    if ($domclustExitStatus != 0) {
    	# restart domclust with --nobreak option
	    my $domclustExitStatus = system "domclust all.fas.hom all.fas.gene -HO -S -c60 -p0.5 -V0.6 -C80 -o5 --nobreak > cluster.tab 2>> $OUT_ERR"; 
	    $domclustExitStatus /= 256; # get original exit value
	    # check if second run was successful
	    if ($domclustExitStatus != 0) {
	    	die("\nERROR: 'domclust' returned with non-zero exit status $domclustExitStatus.\n\n");
	    }
    }

    
    system "grep '>' all.fas | uniq -d > duplicated_CDS.txt";
    if (-s "duplicated_CDS.txt") {
        print ERRORLOG "duplicated CDS for some genes. Please check locus tags:\n";
		my $fileContent = do{local(@ARGV,$/)="duplicated_CDS.txt";<>};
		print ERRORLOG $fileContent . "\n";
    }
}

# 16s sequence parsing 
# system $PATH_COPRA_SUBSCRIPTS . "parse_16s_from_gbk.pl  $GenBankFiles > 16s_sequences.fa" unless (-e "16s_sequences.fa");
system $PATH_COPRA_SUBSCRIPTS . "parse_16s_from_gbk3.pl $core_count $GenBankFiles > 16s_sequences.fa" unless (-e "16s_sequences.fa");
# print $PATH_COPRA_SUBSCRIPTS . "parse_16s_from_gbk3.pl $core_count $GenBankFiles > 16s_sequences.fa\n" if ($verbose);
# check 16s
open(MYDATA, "16s_sequences.fa") or die("\nError: cannot open file 16s_sequences.fa in homology_intaRNA.pl\n\n");
    my @sixteenSseqs = <MYDATA>;
close MYDATA;

my $sixteenScounter = 0;
my $temp_16s_ID = ""; 
foreach (@sixteenSseqs) {
    if ($_ =~ m/>/) {
        $temp_16s_ID = $_; 
        chomp $temp_16s_ID;
        $sixteenScounter++;
    #} else {
       # if ($_ =~ m/N/) { print ERRORLOG "\nError: 'N' characters present in 16s_sequences.fa. Remove $temp_16s_ID from the input for the job to execute correctly.\n"; }
    }
}

if ($sixteenScounter ne $orgcount) {
    my $no16sOrgs = `(grep ">" 16s_sequences.fa && grep ">" input_sRNA.fa) | sort | uniq -u | tr '\n' ' '`; 
    chomp $no16sOrgs;
    print ERRORLOG "\nError: wrong number of sequences in 16s_sequences.fa.\nOne (or more) of your entered organisms does not contain a correctly annotated 16s rRNA sequence and needs to be removed.\nPlease remove $no16sOrgs\n";
}

## prepare single organism whole genome target predictions 
system "echo $GenBankFiles > merged_refseq_ids.txt"; # need this for iterative region plot construction

# my $prepare_intarna_out_call = $PATH_COPRA_SUBSCRIPTS . "prepare_intarna_out.pl $ncrnas $upfromstartpos $down $mrnapart $GenBankFiles";
my $prepare_intarna_out_call = $PATH_COPRA_SUBSCRIPTS . "prepare_intarna_out.pl $ncrnas $upfromstartpos $down $mrnapart $core_count $intarnaParamFile $GenBankFiles";
print $prepare_intarna_out_call . "\n" if ($verbose);
system $prepare_intarna_out_call;
## end

# do CopraRNA combination 
print "\n" . $PATH_COPRA_SUBSCRIPTS . "combine_clusters.pl $orgcount\n\n" if ($verbose);
system $PATH_COPRA_SUBSCRIPTS . "combine_clusters.pl $orgcount";

# make annotations
my $annotateCall = undef;
  $annotateCall = $PATH_COPRA_SUBSCRIPTS . "annotate_raw_output.pl CopraRNA2_prep_sorted.csv opt_tags.clustered $GenBankFiles > CopraRNA2_prep_anno.csv";

print "\n$annotateCall\n" if ($verbose);
system $annotateCall; 

# get additional homologs in cluster.tab
my $parseHomologsCall = undef;
	$parseHomologsCall = $PATH_COPRA_SUBSCRIPTS . "parse_homologs_from_domclust_table.pl CopraRNA2_prep_anno.csv cluster.tab > CopraRNA2_prep_anno_addhomologs.csv";

print "$parseHomologsCall\n" if ($verbose);
system $parseHomologsCall; 

 
# just for formatting
system "awk -F',' '{ print \$1 }' CopraRNA2_prep_anno_addhomologs.csv > CopraRNA2_pvalues.txt";

# just for formatting
system "R --slave -f $PATH_COPRA_SUBSCRIPTS/calc_padj.R --args CopraRNA2_pvalues.txt";
system "paste padj.csv CopraRNA2_prep_anno_addhomologs.csv -d ',' > CopraRNA2_prep_anno_addhomologs_padj.csv"; 


# make consistent names
system $PATH_COPRA_SUBSCRIPTS . "get_amount_sampled_values_and_add_to_table.pl CopraRNA2_prep_anno_addhomologs_padj.csv 1 > CopraRNA2_prep_anno_addhomologs_padj_amountsamp.csv"; 

# get ooi refseq id
my @split = split(/\s/, $refseqaffiliations{$ARGV[6]});
# get the first id entry
my $ooi_refseq_id = $split[0];




######################################################
print "compute phylogenetic distances to the ooi UTRs\n" if ($verbose);
######################################################
system "R --slave -f " . $PATH_COPRA_SUBSCRIPTS . "copraRNA2_phylogenetic_sorting.r 2>> $OUT_ERR 1>&2"; 
# perform actual CopraRNA 2 p-value combination
system "R --slave -f " . $PATH_COPRA_SUBSCRIPTS . "join_pvals_coprarna_2.r 2>> $OUT_ERR 1>&2"; 



# truncate final output // 
system "head -n $topcount CopraRNA_result_all.csv > CopraRNA_result.csv"; 

#######################################################
print "find conserved sites, plot CopraRNA evo heatmap, jalview files for selection\n" if ($verbose);
#######################################################
print "copraRNA2_find_conserved_sites.r\n" if ($verbose);
system "R --slave -f " . $PATH_COPRA_SUBSCRIPTS . "copraRNA2_find_conserved_sites.r 2>> $OUT_ERR 1>&2";
print "copraRNA2_conservation_heatmaps.r\n" if ($verbose);
system "R --slave -f " . $PATH_COPRA_SUBSCRIPTS . "copraRNA2_conservation_heatmaps.r 2>> $OUT_ERR 1>&2"; 


# check for run fail CopraRNA
open(MYDATA, "CopraRNA_result.csv") or die("\nError: cannot open file CopraRNA_result.csv at homology_intaRNA.pl\n\n");
    my @CopraRNA_out_lines = <MYDATA>;
close MYDATA;

if (scalar(@CopraRNA_out_lines) <= 1) { 
    print ERRORLOG "Error: No predictions in CopraRNA_result.csv. CopraRNA run failed.\n"; 
}

# trim off last column (initial_sorting) if CopraRNA 2 prediction mode
system "awk -F',' '{ print \$NF }' CopraRNA_result.csv > CopraRNA_result.map_evo_align" if ($websrv);
system "awk -F, -vOFS=, '{NF-=1;print}' CopraRNA_result.csv > CopraRNA_result_temp.csv";
system "mv CopraRNA_result_temp.csv CopraRNA_result.csv";
system "awk -F, -vOFS=, '{NF-=1;print}' CopraRNA_result_all.csv > CopraRNA_result_all_temp.csv";
system "mv CopraRNA_result_all_temp.csv CopraRNA_result_all.csv";
# change header
system "sed -i 's/,Additional.homologs,/,Additional homologs,/g' CopraRNA_result.csv";
system "sed -i 's/,Amount.sampled/,Amount sampled/g' CopraRNA_result.csv";
system "sed -i 's/p.value/p-value/g' CopraRNA_result.csv";
system "sed -i 's/,Additional.homologs,/,Additional homologs,/g' CopraRNA_result_all.csv";
system "sed -i 's/,Amount.sampled/,Amount sampled/g' CopraRNA_result_all.csv";
system "sed -i 's/p.value/p-value/g' CopraRNA_result_all.csv";


if ($websrv) { 
	#######################################################
	print "generate webserver output\n" if ($verbose); 
	#######################################################

    my $allrefs = $refseqaffiliations{$ARGV[6]};
    my @splitallrefs = split(/\s/,$allrefs);

    my $themainrefid = $splitallrefs[0]; # organism of interest RefSeq ID
    my $orgofintTargets = $themainrefid . "_upfromstartpos_" . $upfromstartpos . "_down_" . $down . ".fa";
    my $orgofintsRNA = "ncRNA_" . $themainrefid . ".fa";
	system "R --slave -f " . $PATH_COPRA_SUBSCRIPTS . "prepare_webserver_output.r";

    system "cp $orgofintTargets target_sequences_orgofint.fa";
}

system $PATH_COPRA_SUBSCRIPTS . "print_archive_README.pl > README.txt";

if ($enrich) { 
	#######################################################
	print "create DAVID enrichment table\n" if ($verbose); 
	#######################################################

    ## this has all been changed to python in version 2.0.3.1 because the DAVID-WS perl client was flawed
    system $PATH_COPRA_SUBSCRIPTS . "DAVIDWebService_CopraRNA.py CopraRNA_result_all.csv $enrich > DAVID_enrichment_temp.txt"; 
    system "grep -P 'termName\\s=|categoryName\\s=|score\\s=|listHits\\s=|percent\\s=|ease\\s=|geneIds\\s=|listTotals\\s=|popHits\\s=|popTotals\\s=|foldEnrichment\\s=|bonferroni\\s=|benjamini\\s=|afdr\\s=' DAVID_enrichment_temp.txt | sed 's/^[ ]*//g' | sed 's/ = /=/g' | sed 's/, /,/g' > DAVID_enrichment_grepped_temp.txt"; ##  only removing obsolete spaces and keeping others
	# ensure there is some enrichment output
	if( -s "DAVID_enrichment_grepped_temp.txt" ) {
        system $PATH_COPRA_SUBSCRIPTS . "make_enrichment_table_from_py_output.pl DAVID_enrichment_grepped_temp.txt > termClusterReport.txt"; 
    	if ( -s "termClusterReport.txt" ) {
	        open(MYDATA, "termClusterReport.txt");
            my @enrichment_lines = <MYDATA>;
       		close MYDATA;
            ## enrichment visualization
            system "cp $PATH_COPRA_SUBSCRIPTS" . "copra_heatmap.html ."; 
            system "R --slave -f " . $PATH_COPRA_SUBSCRIPTS . "extract_functional_enriched.R --args CopraRNA_result_all.csv termClusterReport.txt enrichment.txt";
            system $PATH_COPRA_SUBSCRIPTS . "make_heatmap_json.pl enrichment.txt"; 
            system "cp $PATH_COPRA_SUBSCRIPTS" . "index-thumb.html ."; 
            system "cp $PATH_COPRA_SUBSCRIPTS" . "index-pdf.html ."; 
            system "phantomjs " . $PATH_COPRA_SUBSCRIPTS . "rasterize.js " . "./index-thumb.html enriched_heatmap_big.png"; 
            system "phantomjs " . $PATH_COPRA_SUBSCRIPTS . "rasterize.js " . "./index-pdf.html enriched_heatmap_big.pdf"; 
            system "rm index-thumb.html"; 
            system "rm index-pdf.html"; 
            ## end add enrichment vis
        } else {
			system "echo -e 'If you are reading this, then your prediction did not return an enrichment, your organism of interest is not in the DAVID database\nor the DAVID webservice is/was termporarily down. You can either rerun your CopraRNA\nprediction or create your enrichment manually at the DAVID homepage.' > termClusterReport.txt";
		}
    } else {
		system "echo -e 'If you are reading this, then your prediction did not return an enrichment, your organism of interest is not in the DAVID database\nor the DAVID webservice is/was termporarily down. You can either rerun your CopraRNA\nprediction or create your enrichment manually at the DAVID homepage.' > termClusterReport.txt";
	}

    ##### end DAVID enrichment

}


# close error stream and be done
close ERRORLOG;
exit(0);

