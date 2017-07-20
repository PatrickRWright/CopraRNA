NAME=myco
BLASTOPT=-Y 1e9 -v 30000 -b 30000 -m 8
OSNAME=$(shell uname -s)

all:	$(NAME).hom $(NAME).gene

$(NAME).pin:	$(NAME)
	/usr/local/ncbiblast/2.2.15/bin/formatdb -i $<

$(NAME).blast:	$(NAME) $(NAME).pin
	/usr/local/ncbiblast/2.2.15/bin/blastall -p blastp -d $(NAME) -a 4 -i $(NAME) $(BLASTOPT) -o $@

$(NAME).hom:	$(NAME).blast
	./blast2homfile.pl -distconv $< > $@

$(NAME).gene:	$(NAME)
	./fasta2genefile.pl $<

test:	$(NAME).gene $(NAME).hom
	../bin-$(OSNAME)/domclust $(NAME).hom $(NAME).gene | ../Script/addtit.pl -titfile=$(NAME).tit

cleantmp:
	-rm *.blast *.pin *.psq *.phr *.log 2> /dev/null

clean:
	-rm *.gene *.hom *.blast *.pin *.psq *.phr *.tit *.log 2> /dev/null
