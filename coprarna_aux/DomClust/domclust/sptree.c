/*
 * DomClust: Hierarchical Clustering for Orthologous Domain Classification
 * Copyright (c) 2000-2006, Ikuo Uchiyama
 * All rights reserved.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "domclust.h"
#include "namehash.h"
#include "spec.h"

static struct {
	int prev_child;
	int count;
} node_buf[MAXSP*2];
static int sp_treenode_idx[MAXSP];
int sptree_create_newnode(int spid, int bufidx, int parent, char flag);

readSPfile(char *spfile)
{
	FILE *fp;
	static char buf[BUFSIZ];
	char *patbuf, *p;
	int i;
	char flag = 1;
	int allocsize = 2000;

	if ((fp = fopen(spfile, "r")) == NULL) {
		fprintf(stderr, "Can't open file\n");
		exit(1);
	}
	if ((patbuf = (char *)malloc(allocsize))==NULL) {
		fprintf(stderr, "malloc failed: patbuf\n");
		exit(1);
	}
	i = 0;
	while(fgets(buf, BUFSIZ, fp)!=NULL) {
		for (p = buf; *p; p++) {
			if (*p == '#'){
				flag = 0;
			} else if (*p == '\n') {
				flag = 1;
				continue;
			}
			if (! flag) {
				continue;
			}
			patbuf[i] = *p;
			if (++i >= allocsize) {
				allocsize *= 1.5;
				if ((patbuf = (char *)realloc(patbuf, allocsize))==NULL) {
					fprintf(stderr, "malloc failed: patbuf\n");
					exit(1);
				}
			}
		}
	}
	patbuf[i] = '\0';
	parse_spinfo(patbuf);
	fclose(fp);
}
parse_spinfo(char *str)
{
	char *p;
	double weight;
	int i, l, id;
	char *sp, spname_buf[10];
	enum {NAME,OTHER} status = OTHER;
	pList *pl;
	int intnode_buf[MAXSP], *parentp;
	int bufidx, parent;
	char buf[BUFSIZ];

	if (! SpHash) {
		SpHash = initNames(MAXSP*2);
	}

/***
      (  A  (  B  C) (  D  E  (  F  G)  H)  I)
idx   0  1  2  3  4  5  6  7  8  9 10  11  12
par  -1  0  0  2  2  0  5  5  5  8  8   5   0
spid -1  A -1  B  C -1  D  E -1  F  G   H   I
1/w0  1  4  4  2  2  4  4  4  4  2  2   4   4

parentp(id,cnt):  (0,4)->(2,2)->(5,3)->(8,2)
***/

	pl = create_pList();
	/* root */
	spTree.node[0].spid = -1;
	spTree.node[0].parent = -1;
	spTree.node[0].weight = 1;

	parentp = &intnode_buf[0];
	node_buf[0].prev_child = -1;
	*parentp = 0;

	bufidx = 1;

	for (p = str; *p; p++) {
		if (isalnum(*p)) {
			if (status == NAME) {
				l++;
			} else {
				status = NAME;
				sp = p;
				l = 1;
			}
		} else if (status == NAME) {
			strncpy(spname_buf,sp,l);
			spname_buf[l] = 0;
			id = getSPid(spname_buf);
			/* create new node */
			sptree_create_newnode(id, bufidx, *parentp, (char)0);
			bufidx++;
			status = OTHER;
		}
		if (status == OTHER) {
			if (*p == '(' || *p == '{') {
				char flag = 0;
				if (*p == '{') {
					flag = FLAG_TAXOR;
				}
				sptree_create_newnode(-1, bufidx, *parentp,flag);

				/* save current parent and
					create a new internal node */
				pushList(pl, parentp);
				parentp++;
				*parentp = bufidx;

				bufidx++;
			} else if (*p == ')' || *p == '}') {
				if ((parentp = (int *) popList(pl)) == NULL) {
					parse_spinfo_error();
				}
			} else if (*p == ':') {
				weight = strtod(p, &p);
				if (p == NULL) {
					parse_spinfo_error();
				}
			} else {
			}
		}
	}
	spTree.nodenum = bufidx;
	for (i = 0; i < bufidx; i++) {
		if ((parent = spTree.node[i].parent) == -1) continue;
		if (parent == 0) continue;
		if (spTree.node[parent].flag & FLAG_TAXOR) {
			/* taxonomic OR: do not count each child */
			spTree.node[i].weight = 0;
		} else {
			spTree.node[i].weight = spTree.node[parent].weight
						/ node_buf[parent].count;
		}
		SPweights[spTree.node[i].spid] = spTree.node[i].weight;
	}
/*
	sptree_set_spflag();
*/
	free_pList(pl);
}
sptree_create_newnode(int spid, int bufidx, int parent, char flag)
{
	int prev_sibling;
	spTree.node[bufidx].spid = spid;
	spTree.node[bufidx].parent = parent;
	spTree.node[bufidx].weight = 1;
	spTree.node[bufidx].child = -1;
	spTree.node[bufidx].sibling = -1;
	spTree.node[bufidx].flag = flag;
	node_buf[bufidx].count = 0;
	node_buf[bufidx].prev_child = -1;
	node_buf[parent].count++;
	if ((prev_sibling = node_buf[parent].prev_child) >= 0) {
		spTree.node[prev_sibling].sibling = bufidx;
	} else {
		spTree.node[parent].child = bufidx;
	}
	node_buf[parent].prev_child = bufidx;
	if (spid>=0) {
		setSPflag(spTree.node[bufidx].spflag, spid);
		sp_treenode_idx[spid] = bufidx;
	}
}
sptree_add_species()
{
	int i;
	int spnum = getSPnum();
/* add species that are not found in sptreefile: called from preproc */
	for (i = 0; i < spnum; i++) {
		if (sp_treenode_idx[i] == 0) {
			sptree_create_newnode(i, spTree.nodenum++, 0, (char)0);
		}
	}
	sptree_set_spflag();
}
sptree_set_spflag()
{
	int i, j;
	for (i = spTree.nodenum-1; i >= 0; i--) {
		clearSPflag( spTree.node[i].spflag );
		if (spTree.node[i].child < 0) {
			setSPflag(spTree.node[i].spflag, spTree.node[i].spid);
		} else {
			for (j = spTree.node[i].child; j >= 0;
				j = spTree.node[j].sibling) {
				spFlagOR( spTree.node[i].spflag,
					spTree.node[j].spflag,
					spTree.node[i].spflag);
			}
		}
	}
}
/* return the deepest node that includes the set of species specified by spflag */
sptree_MatchFlags(specFlag *spflag)
{
	register int i;
	for (i = spTree.nodenum-1; i >= 0; i--) {
		if (spFlagInclude(spTree.node[i].spflag, spflag)) {
			return i;
		}
	}
	return -1;
}
double sptree_spFlagCountTaxOrW(specFlag spflag)
{
	register int i;
	double cnt = 0;
	for (i = 0; i < spTree.nodenum; i++) {
		if (spTree.node[i].flag & FLAG_TAXOR) {
			if (spFlagANDcnt(spTree.node[i].spflag, spflag)) {
				cnt += spTree.node[i].weight;
			}
		}
	}
	return cnt;
}

double sptree_MatchFlagsCntW(specFlag spflag1, specFlag spflag2, int nodenum)
{
	register int i, j;
	int cnt1, cnt2;
	static specFlag tmp_spflag1, tmp_spflag2, tmp_spflag3, tmp_spflag;
	double cnt;
	
	cnt = 0.0;
	clearSPflag(tmp_spflag);
	clearSPflag(tmp_spflag1);
	clearSPflag(tmp_spflag2);
	clearSPflag(tmp_spflag3);

	for (i = spTree.node[nodenum].child; i >= 0;
			i = spTree.node[i].sibling) {
		spFlagAND(spflag1, spTree.node[i].spflag, tmp_spflag1);
		spFlagAND(spflag2, spTree.node[i].spflag, tmp_spflag2);
		spFlagOR(tmp_spflag1, tmp_spflag2, tmp_spflag3);
		spFlagOR(tmp_spflag, tmp_spflag3, tmp_spflag);
	}
	cnt = spFlagCntW(tmp_spflag3);
	return cnt;
}
sptree_MatchFlagsInv(specFlag *spflag)
{
	register int i;
	for (i = spTree.nodenum-1; i >= 0; i--) {
		if (spFlagInclude(spflag, spTree.node[i].spflag)) {
			return i;
		}
	}
	return -1;
}
print_sptree()
{
	int i;
printf("%d\n",spTree.nodenum);
	for (i = 0; i < spTree.nodenum; i++) {
		printf("%3d %3d %3d %3d %s %3d %1d %7.3lf ",
			i,spTree.node[i].parent,
			spTree.node[i].child, spTree.node[i].sibling,
			(spTree.node[i].spid >= 0 ?
				SPnames[spTree.node[i].spid] : "---"),
			spTree.node[i].spid, spTree.node[i].flag,
			spTree.node[i].weight*100);
		print_specFlag(spTree.node[i].spflag);
	}
}
parse_spinfo_error()
{
	fprintf(stderr, "parse error\n");
	exit(1);
}

#ifdef DEBUGMAIN_SPTREE
main(int argc,char **argv)
{
	int c;
	specFlag spf;
	if (argc <=1){
		fprintf(stderr, "Usage: sptree spfile\n");
		exit(1);
	}
	readSPfile(argv[1]);
	sptree_set_spflag();
	print_sptree();

/***
	setSPflag(spf,0);
	addSPflag(spf,1);
	addSPflag(spf,2);
	addSPflag(spf,3);
	addSPflag(spf,5);
	addSPflag(spf,6);
	print_specFlag(spf);
	c =sptree_spFlagCountWithTaxOr(spf);
	printf("%d\n",c);
**/
}
#endif
