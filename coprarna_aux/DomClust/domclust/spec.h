/*
 * DomClust: Hierarchical Clustering for Orthologous Domain Classification
 * Copyright (c) 2000-2006, Ikuo Uchiyama
 * All rights reserved.
 */

#ifndef _SPEC_H_
#define _SPEC_H_

#define SPFLAGSIZ 80
#define MAXSP (SPFLAGSIZ * sizeof(char) * 8)
#define MAXTAXNUM 20


int SPnum;
/*
int bitcnt[256] = {
	0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4,1,2,2,3,2,3,3,4,2,3,3,4,3,4,
	4,5,1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,
	4,5,5,6,1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,
	4,5,4,5,5,6,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,
	4,5,5,6,5,6,6,7,1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,
	4,5,3,4,4,5,4,5,5,6,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,
	4,5,5,6,4,5,5,6,5,6,6,7,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,
	4,5,4,5,5,6,4,5,5,6,5,6,6,7,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
	4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8
};
double wbitcnt[SPFLAGSIZ][256];
*/
NameHash *SpHash;
char *SPnames[MAXSP];
double SPweights[MAXSP];
typedef unsigned char specFlag[SPFLAGSIZ];
typedef unsigned char *specFlagP;

#define FLAG_TAXOR 1
#define FLAG_TAXAND 2

typedef struct SPTreeNode {
	int spid;
	int parent;
	int child, sibling;
	double weight;
	specFlag spflag;
	char flag;
} SPTreeNode;

typedef struct SPTree {
	SPTreeNode node[MAXSP*2];
	int nodenum;
} SPTree;
SPTree spTree;

struct {
	specFlag inGroup, outGroup;
	specFlag spMask;
} SPflags;

char *getSPname();
double spFlagCntW(), spFlagANDcntW(), getSPweight(), sptree_MatchFlagsCntW(),
	sptree_spFlagCountTaxOrW();
int parse_spinfo(char *);
int readSPfile(char *);
#endif
