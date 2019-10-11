# standard test call 
CopraRNA2.pl -srnaseq sRNAs.fa -ntup 200 -ntdown 100 -region 5utr -enrich 200 -topcount 200 -cores 4

# standard test call noclean
CopraRNA2.pl -srnaseq sRNAs.fa -ntup 200 -ntdown 100 -region 5utr -enrich 200 -topcount 200 -cores 4 --noclean
