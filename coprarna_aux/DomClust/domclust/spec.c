/*
 * DomClust: Hierarchical Clustering for Orthologous Domain Classification
 * Copyright (c) 2000-2006, Ikuo Uchiyama
 * All rights reserved.
 */

#include <stdio.h>
#include <string.h>
#include "domclust.h"
#include "spec.h"

static int spflagsiz = SPFLAGSIZ;
static unsigned char bitmask[8] = {1,2,4,8,16,32,64,128};

static int bitcnt[256] = {
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
static double wt_bitcnt[SPFLAGSIZ][256];
static int taxFlag[MAXTAXNUM][SPFLAGSIZ][256];
int _getSPid(char *sp, char add);

getSPid(char *sp)
{
	return _getSPid(sp, 1);
}
_getSPid(char *sp, char add)
{
	int id;
	if (! SpHash) {
		SpHash = initNames(MAXSP*2);
	}
	if ( (id = getNameID(SpHash, sp)) < 0 ) {
	    if (add) {
		/* Not found in the hash. Add this name and assign a new ID */
		id = SPnum++;
		if (SPnum >= MAXSP) {
			fprintf(stderr,
				"Too many species. Raise SPFLAGSIZ value. (%d,%d)\n",SPnum,MAXSP);
			exit(1);
		}
		setSPnum(SPnum);
		SPnames[id] = addName(SpHash, sp, id);
		SPweights[id] = 1.0;
	    } else {
		return -1;
	    }
	}
	return id;
}
char *getSPname(int id)
{
	if (id < SPnum) {
		return SPnames[id];
	} else {
		return NULL;
	}
}
setSPnum(int spnum)
{
	SPnum = spnum;
	spflagsiz = spnum / 8 + 1;
}
getSPnum()
{
	return SPnum;
}
setSPweight(int id, double w)
{
	SPweights[id] = w;
}
double getSPweight(int id)
{
	return SPweights[id];
}

setSpMask(char *str)
{
	setSPflagBySPname(str, SPflags.spMask);
}
setOutGroup(char *str, char rev)
{
	if (rev) {
		setSPflagBySPname(str, SPflags.inGroup);
		/* reverse the inGroup flags to obtain the outGroup flags */
		spFlagCOMPL(SPflags.inGroup, SPflags.outGroup);
	} else {
		setSPflagBySPname(str, SPflags.outGroup);
		spFlagCOMPL(SPflags.outGroup, SPflags.inGroup);
	}
}
setSPflagBySPname(char *str, specFlagP readGrpP)
{
	char *sp;
	int id;
	static char *Delim = ":, ";

	if (! str) return 0;
	if ( (sp = strtok(str, Delim)) == NULL) {
		return 0;
	}
	if ( (id = _getSPid(sp, 0)) < 0 ) {
		fprintf(stderr, "Not found: %s\n",sp);
	} else {
		addSPflag(readGrpP, id);
	}
	while ( (sp = strtok(NULL, Delim)) != NULL ) {
		if ( (id = _getSPid(sp, 0)) < 0 ) {
			fprintf(stderr, "Not found: %s\n",sp);
		} else {
			addSPflag(readGrpP, id);
		}
	}
	return 0;
}
spFlagInGroup(specFlag spflag, specFlag newflag)
{
	spFlagAND(spflag, SPflags.inGroup, newflag);
}
spFlagOutGroup(specFlag spflag, specFlag newflag)
{
	spFlagAND(spflag, SPflags.outGroup, newflag);
}

/** preproc_for_SpecInfo():
	this routine should be called before the clustering process
	but after reading the homology data (or reading all of
	the species names) **/
preproc_for_SpecInfo()
{
	sptree_add_species();
	create_wt_bitcnt();
}
create_wt_bitcnt()
{
	int i, j, jj, k;
	for (i = 0; i < spflagsiz; i++) {
		for (j = 0; j < 256; j++) {
			wt_bitcnt[i][j] = 0;
			jj = j;
			for (k = 0; k < 8 && jj>0; k++) {
				if (jj % 2 == 1) {
					wt_bitcnt[i][j] += SPweights[i*8+k];
				}
				jj /= 2;
			}
		}
	}
}
clearSPflag(specFlag spflag)
{
	bzero(spflag, SPFLAGSIZ);
}
setSPflag(specFlag spflag, int spnum)
{
	clearSPflag(spflag);
	spflag[spnum / 8] = bitmask[spnum % 8];
}
addSPflag(specFlag spflag, int spnum)
{
	spflag[spnum / 8] |= bitmask[spnum % 8];
}
copySPFlag(specFlag spflag1, specFlag spflag2)
{
	int i;
	for (i = 0; i < spflagsiz; i++) {
		spflag2[i] = spflag1[i];
	}
}

spFlagOR(specFlag spflag1, specFlag spflag2, specFlag newflag)
{
	register int i;
	for (i = 0; i < spflagsiz; i++) {
		newflag[i] = spflag1[i] | spflag2[i];
	}
}
spFlagAND(specFlag spflag1, specFlag spflag2, specFlag newflag)
{
	register int i;
	for (i = 0; i < spflagsiz; i++) {
		newflag[i] = spflag1[i] & spflag2[i];
	}
}
spFlagCOMPL(specFlag spflag, specFlag newflag)
{
	int i;
	for (i = 0; i < spflagsiz; i++) {
		newflag[i] = ~ spflag[i];
	}
}
spFlagANDcnt(specFlag spflag1, specFlag spflag2)
{
	register int i;
	register int cnt = 0;
	for (i = 0; i < spflagsiz; i++) {
		cnt += bitcnt[spflag1[i] & spflag2[i]];
	}
	return cnt;
}
double spFlagANDcntW(specFlag spflag1, specFlag spflag2)
{
	register int i;
	double cnt = 0;
	specFlag tmpflag;
	spFlagAND(spflag1, spflag2, tmpflag);
	return spFlagCntW(tmpflag);
}
spFlagCnt(specFlag spflag)
{
	register int i;
	register int cnt = 0;
	for (i = 0; i < spflagsiz; i++) {
		cnt += bitcnt[spflag[i]];
	}
	return cnt;
}
double spFlagCntW(specFlag spflag)
{
	register int i;
	register double cnt = 0;
	double cnt2 = 0;
	for (i = 0; i < spflagsiz; i++) {
		cnt += wt_bitcnt[i][spflag[i]];
	}
	cnt2 = sptree_spFlagCountTaxOrW(spflag);
	return cnt + cnt2;
}
/* spflag1 >= (includes) spflag2 */
spFlagInclude(specFlag spflag1, specFlag spflag2)
{
	return (spFlagCnt(spflag2) == spFlagANDcnt(spflag1, spflag2));
}

print_specFlag(specFlag spflag)
{
	int i;
	for (i = 0; i < SPnum; i++) {
		printf("%d", (spflag[i / 8] & bitmask[i % 8]) != 0);
	}
	putchar('\n');
}
dump_specFlag(FILE *ofp, specFlag spflag)
{
	int i;
	for (i = 0; i < spflagsiz; i++) {
		fprintf(ofp, " %d", spflag[i]);
	}
}
restore_specFlag(char *str, specFlag spflag)
{
	int i;
	for (i = 0; i < spflagsiz; i++) {
		spflag[i] = (unsigned char) strtol(str, &str, 10);
		if (! str) break;
	}
}
