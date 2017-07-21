
# omit incompatible / screwed up genomes
#sed -i '/NC_022785/d' kegg2refseqnew.csv
#sed -i '/NC_022785/d' CopraRNA_available_organisms.tmp

# make help list for webserver
echo "## last updated: \c" > Date.txt; date >> Date.txt
cat Date.txt ../head.txt > fullhead.txt
cat fullhead.txt CopraRNA_available_organisms.tmp > CopraRNA_available_organisms.txt
rm CopraRNA_available_organisms.tmp fullhead.txt Date.txt

