#!/usr/bin/env bash

##########################################################

# default: remove temporary files at end of call
KEEPTMPFILES=
# comment next line to remove temporary files at end of cal
KEEPTMPFILES="--noclean"

TESTDIR=$PWD
TESTFILE=$0

##########################################################

CALLARGS=" -verbose -srnaseq $TESTDIR/test1.fa -genomePath $TESTDIR/$TESTFILE-genomes -ntup 200 -ntdown 100 -region 5utr -enrich 200 -topcount 200 -cores 4 $KEEPTMPFILES"

##########################################################

TESTTMPDIR=$TESTFILE-tmp
echo "create and enter test directory '$TESTDIR/$TESTTMPDIR'"
rm -rf $TESTTMPDIR
mkdir -p $TESTTMPDIR
cd $TESTTMPDIR

echo "run test call : $CALLARGS"
../../CopraRNA2.pl $CALLARGS
CALLEXITCODE=$?
echo "exit code of call = $CALLEXITCODE"

# leave test directory
cd $TESTDIR


if [ "$CALLEXITCODE" -eq "0" ]; then
	DIFFOUT=$TESTFILE-result.diff
	diff $TESTTMPDIR/CopraRNA_result_all.csv $TESTFILE-results/CopraRNA_result_all.csv > $DIFFOUT;
	if [ -s "$DIFFOUT" ]; then
		echo "TEST FAILURE: check $DIFFOUT for details!";
		exit -1;
	else
		echo "TEST SUCCESSFUL!";
	fi
fi

# exit test with test's exit code
exit $CALLEXITCODE

