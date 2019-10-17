#!/usr/bin/env bash

##########################################################

# default: remove temporary files at end of call
KEEPTMPFILES=
# comment next line to remove temporary files at end of cal
KEEPTMPFILES="--noclean"

TESTDIR=$PWD

##########################################################

CALLARGS=" -verbose -srnaseq $TESTDIR/test1.fa -genomePath $TESTDIR/genomes -ntup 200 -ntdown 100 -region 5utr -enrich 200 -topcount 200 -cores 4 $KEEPTMPFILES"

##########################################################

TESTTMPDIR=test1-tmp
echo "create and enter test directory '$TESTDIR/$TESTTMPDIR'"
rm -rf $TESTTMPDIR
mkdir -p $TESTTMPDIR
cd $TESTTMPDIR

echo "run test call : $CALLARGS"
../../CopraRNA2.pl $CALLARGS

# leave test directory
cd $TESTDIR

echo "end of test: exit code = $?"

# exit test with test's exit code
exit $?

