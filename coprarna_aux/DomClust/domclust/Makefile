OSNAME := $(shell uname -s)
BINDIR=bin-$(OSNAME)

sinclude ../conf/$(OSNAME).conf

## You might need to edit the following lines
CC=gcc
CFLAGS=-g -O2 $(CFLAGS0)
#####

## Probably you need not edit the following lines

OBJS1 = domclust.o util.o graph.o memalloc.o namehash.o hash.o \
	cluster.o clustout.o vararray.o plist.o bin.o seqreg.o \
	spec.o sptree.o domain.o readfile.o dump.o neighbor.o \
	postproc.o mergeclust.o newali.o

SELOBJS = select.o memalloc.o plist.o hash.o namehash.o

OBJS = $(OBJS1)

TESTOPT=-S60 -C80 -V.6 -ao.8

all: $(BINDIR)
	echo $(BINDIR)
	cd $(BINDIR); \
	$(MAKE) -f ../Makefile VPATH=.. OSNAME=$(OSNAME) all_target

all_target:	domclust

$(BINDIR):
	@if [ ! -d $(BINDIR) ]; then mkdir $(BINDIR); fi

domclust:	$(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -o $@ -lm

sptree:	sptree.c plist.c namehash.c memalloc.c hash.c spec.c
	$(CC) -g sptree.c plist.c namehash.c memalloc.c hash.c spec.c -DDEBUGMAIN_SPTREE -o $@

clean:
	-rm -rf $(BINDIR)

select:	$(SELOBJS)
	$(CC) $(CFLAGS) $(SELOBJS) -o $@

test:	all
	cd tst;\
	../$(BINDIR)/domclust $(TESTOPT) test.hom test.gene | \
		../Script/cmpr.pl -verbose -OVLPRATIO2=.5 -diff test.out -

testout:	all
	cd tst;\
	../$(BINDIR)/domclust $(TESTOPT) test.hom test.gene > test.out

tarfile:
	tar cvf domclust.tar *.c *.h Makefile Script conf tst
