/*
 * DomClust: Hierarchical Clustering for Orthologous Domain Classification
 * Copyright (c) 2000-2006, Ikuo Uchiyama
 * All rights reserved.
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "seqreg.h"
#include "domclust.h"

#ifdef DEBUGMAIN
main(int argc, char *argv)
{
	Region reg1, reg2, mreg;
	reg1.from = 100; reg1.to = 200;
	reg2.from = 200; reg2.to = 350;
	meanSeqReg(&reg1, 2, &reg2, 1, &mreg);
	printf("%d,%d\n",mreg.from, mreg.to);
}
#endif

setSeqReg(Region *reg, SeqPos from, SeqPos to)
{
	reg->from = from;
	reg->to = to;
}
/* convert from Region to StrRegion (for stroage) */
Reg2StrReg(StrRegion *reg, Region *reg2)
{
	reg->from = reg2->from; reg->to = reg2->to;
}
/* convert from StrRegion to Region (for calculation) */
StrReg2Reg(Region *reg, StrRegion *reg2)
{
	reg->from = reg2->from; reg->to = reg2->to;
}
copyReg(Region *reg, Region *reg2)
{
	reg->from = reg2->from;
	reg->to = reg2->to;
}

wmeanSeqReg(Region *reg1, Region *reg2, Region *newreg)
{
	meanSeqReg(reg1, 1, reg2, 1, newreg);
}
meanValue(int val1, int cnt1, int val2, int cnt2)
{
	return ((val1 * cnt1) + (val2 * cnt2)) / (cnt1 + cnt2);
}
meanSeqReg(Region *reg1, int cnt1, Region *reg2, int cnt2, Region *newreg)
{
	/* weighted mean */
	if (definedPos(reg1->from) && definedPos(reg2->from)) {
		newreg->from = (reg1->from * cnt1 + reg2->from * cnt2)
				/ (cnt1+cnt2);
	} else if (definedPos(reg1->from)) {
		newreg->from = reg1->from;
	} else if (definedPos(reg2->from)) {
		newreg->from = reg2->from;
	} else {
		newreg->from = INFPOS;
	}
	if (definedPos(reg1->to) && definedPos(reg2->to)) {
		newreg->to = (reg1->to * cnt1 + reg2->to * cnt2) / (cnt1+cnt2);
	} else if (definedPos(reg1->to)) {
		newreg->to = reg1->to;
	} else if (definedPos(reg2->to)) {
		newreg->to = reg2->to;
	} else {
		newreg->to = SUPPOS;
	}
}
mergeReg(Region *reg1, Region *reg2, Region *newreg)
{
	newreg->from = (reg1->from < reg2->from) ? reg1->from : reg2->from;
	newreg->to = (reg1->to > reg2->to) ? reg1->to : reg2->to;
}
ovlpReg(Region *reg1, Region *reg2, Region *andreg, Region *orreg, int minov)
{
	char flag = 0;
	if (reg1->to < reg2->from + minov) {
		return -1;
	} else if (reg2->to < reg1->from + minov) {
		return 1;
	} else {
		if (reg1->from < reg2->from) {
			if (andreg) andreg->from = reg2->from;
			if (orreg) orreg->from = reg1->from;
		} else {
			if (andreg) andreg->from = reg1->from;
			if (orreg) orreg->from = reg2->from;
		}
		if (reg1->to < reg2->to) {
			if (andreg) andreg->to = reg1->to;
			if (orreg) orreg->to = reg2->to;
		} else {
			if (andreg) andreg->to = reg2->to;
			if (orreg) orreg->to = reg1->to;
		}
	}
	return 0;
}
ovlpReg2(Region *reg1, Region *reg2, Region *andreg, Region *orreg,
			int minov, float ratio)
{
	if (regLen(reg1) < minov) {
		minov = regLen(reg1) * ratio;
	}
	if (regLen(reg2) < minov) {
		minov = regLen(reg2) * ratio;
	}
	return ovlpReg(reg1, reg2, andreg, orreg, minov);
}
ovlpReg3(Region *reg1, Region *reg2, Region *andreg, Region *orreg, float ratio)
{
	int minlen, len1 = regLen(reg1), len2 = regLen(reg2);
	if (len1 < len2) {
		minlen = len1 * ratio;
	} else {
		minlen = len2 * ratio;
	}
	return ovlpReg(reg1, reg2, andreg, orreg, minlen);
}
#define EQ(a,b) (abs((a)-(b))<=(MinOV))
#define LT(a,b) ((a) < (b) - (MinOV))
/* compare overlapping patterns of segment pairs
		{reg11,reg12} and {reg21,reg22} */
cmprOvlpStatusReg(Region *reg11, Region *reg12, Region *reg21, Region *reg22,
	int minov)
{
	int status1, status2;
	int dist11, dist12, dist13, dist21, dist22, dist23;
	status1 = getOvlpStatusReg(reg11, reg12, minov,
				&dist11, &dist12, &dist13);
	status2 = getOvlpStatusReg(reg21, reg22, minov,
				&dist21, &dist22, &dist23);
	if ( cmprOvlpStat(status1,status2) ) {
		if (distcheck(dist11,dist21) && distcheck(dist12,dist22)
			&& distcheck(dist13,dist23)) {
			/* equivalent */
			return 1;
		}
	}
	/* not equivalent */
	return 0;
}
getOvlpStatusReg(Region *reg1, Region *reg2, int MinOV,
		int *dist1, int *dist2, int *dist3)
{
	int status, prevpos;
	if (LT(reg1->to,reg2->from)) {
		/* T1 < F2 (seg1 << seg2) */
		status = 1;
	} else if ( EQ(reg1->to, reg2->from) ) {
		/* T1 == F2 (seg1 < seg2) */
		status = 2;
	} else if (LT(reg1->from,reg2->from) && LT(reg1->to,reg2->to)) {
		/* F1 < F2 and F2 < T1 < T2 (seg1 <= seg2) */
		status = 3;
	} else if (EQ(reg1->from, reg2->from)  && LT(reg1->to,reg2->to)) {
		/* F1 == F2 and T1 < T2 */
		status = 4;
	} else if (LT(reg1->from,reg2->from) && EQ(reg1->to,reg2->to)) {
		/* F1 < F2 and T1 == T2 */
		status = 5;
	} else if (LT(reg2->from,reg1->from) && LT(reg1->to,reg2->to)) {
		/* F2 < F1 and T1 < F2 (seg2 inc. seg1) */
		status = 6;
	} else if (EQ(reg1->from,reg2->from) && EQ(reg1->to,reg2->to)) {
		/* F1 == F2 and T1 == T2 (seg1==seg2) */
		status = 7;
	} else if (LT(reg1->from,reg2->from) && LT(reg2->to,reg1->to)) {
		/* F1 < F2 and T2 < F1 (seg1 inc. seg2) */
		status = 8;
	} else if (LT(reg2->from,reg1->from) && EQ(reg2->to,reg1->to)) {
		/* F2 < F1 and T2 == T1 */
		status = 9;
	} else if (EQ(reg2->from,reg1->from) && LT(reg2->to,reg1->to)) {
		/* F2 == F1 and T2 < T1 */
		status = 10;
	} else if (LT(reg2->from,reg1->from) && LT(reg2->to,reg1->to)) {
		/* F2 < F1 < T2 and T2 < T1 */
		status = 11;
	} else if (LT(reg2->from,reg1->from) && EQ(reg2->to,reg1->from)) {
		/* T2 == F1 (seg2 < seg1) */
		status = 12;
	} else if (LT(reg2->to,reg1->from)) {
		/* T2 < F1 (seg2 << seg1)*/
		status = 13;
	} else {
		status = 14;
		fprintf(stderr, "error: %d,%d,%d,%d\n",
			reg1->from,reg1->to,reg2->from, reg2->to);
	}
	if (reg1->from < reg2->from) {
		if (reg2->from < reg1->to) {
			*dist1 = reg2->from - reg1->from;
			if (reg1->to < reg2->to) {
				*dist2 = reg1->to - reg2->from;
				*dist3 = reg2->to - reg1->to;
			} else {
				*dist2 = reg2->to - reg2->from;
				*dist3 = reg1->to - reg2->to;
			}
		} else {
			*dist1 = reg1->to - reg1->from;
			*dist2 = reg2->from - reg1->to;
			*dist3 = reg2->to - reg2->from;
		}
	} else {
		if (reg1->from < reg2->to) {
			*dist1 = reg1->from - reg2->from;
			if (reg2->to < reg1->to) {
				*dist2 = reg2->to - reg1->from;
				*dist3 = reg1->to - reg2->to;
			} else {
				*dist2 = reg1->to - reg1->from;
				*dist3 = reg2->to - reg1->to;
			}
		} else {
			*dist1 = reg2->to - reg2->from;
			*dist2 = reg1->from - reg2->to;
			*dist3 = reg1->to - reg1->from;
		}
	}
	return status;
}
cmprOvlpStat(int stat1, int stat2)
{
	int equiv[13][13] = {
		{1,1,0,0,0,0,0,0,0,0,0,0,0},
		{1,1,1,0,0,0,0,0,0,0,0,0,0},
		{0,1,1,1,1,0,0,0,0,0,0,0,0},
		{0,0,1,1,0,1,1,0,0,0,0,0,0},
		{0,0,1,0,1,0,1,1,0,0,0,0,0},
		{0,0,0,1,0,1,0,0,1,0,0,0,0},
		{0,0,0,1,1,0,1,0,1,1,0,0,0},
		{0,0,0,0,1,0,0,1,0,1,0,0,0},
		{0,0,0,0,0,1,1,0,1,0,1,0,0},
		{0,0,0,0,0,0,1,1,0,1,1,0,0},
		{0,0,0,0,0,0,0,0,1,1,1,1,0},
		{0,0,0,0,0,0,0,0,0,0,1,1,1},
		{0,0,0,0,0,0,0,0,0,0,0,1,1},
	};
	return (equiv[stat1-1][stat2-1]);
}
distcheck(int dist1, int dist2)
{
	double ratio = 0.5;
	int shortdiff = 20;
	if (dist1 == dist2) {
		return 1;
	} else if ( (abs(dist1 - dist2) < shortdiff) ) {
		return 1;
	} else if (dist1 > dist2) {
		if ((double) dist2 / dist1 > ratio) {
			return 1;
		}
	} else if (dist2 > 0) {
		if ((double) dist1 / dist2 > ratio) {
			return 1;
		}
	}
	return 0;
}
double ovlpRatio(Region *reg1, Region *reg2, Region *ANDreg)
{
	Region tmpANDreg, ORreg;
	int len1 = regLen(reg1), len2 = regLen(reg2);
	int smallLen = (len1 < len2 ? len1 : len2);
	int ovlp = 0;
	
	if (ANDreg == NULL) {
		ANDreg = &tmpANDreg;
		ovlp = ovlpReg(reg1, reg2, ANDreg, &ORreg, 0);
	}
	if (ovlp == 0) {
		return (double) regLen(ANDreg) / smallLen;
	} else {
		return 0.0;
	}
}
double gapRatio(Region *ali1, Region *ali2)
{
	int diag1 = ali1->from - ali2->from;
	int diag2 = ali1->to - ali2->to;
	int len1 = regLen(ali1), len2 = regLen(ali2);
	int smallLen = (len1 < len2 ? len1 : len2);
	int largeLen = (len1 > len2 ? len1 : len2);
	return (double) (fabs(diag1 - diag2) - fabs(len1 - len2)) / smallLen;
}
definedReg(Region *reg)
{
	return(definedPos(reg->from) && definedPos(reg->to));
}

/* reg is split by reg2 */
/*    L   |<---reg2--->| R */
/*     |<---reg-->| */

splitReg(Region *reg, Region *reg2, Region *L, Region *M, Region *R, int ovlen)
{
	SeqPos maxfrom, minto;

	if (reg2->from - reg->from > ovlen) {
		L->from = reg->from;
		L->to = (reg->to < reg2->from) ? reg->to : reg2->from;
	} else {
		L->from = L->to = SUPPOS;
	}
	maxfrom = (reg->from < reg2->from) ? reg2->from : reg->from;
	minto = (reg->to < reg2->to) ? reg->to : reg2->to;
	if (minto - maxfrom > ovlen) {
		M->from = maxfrom; M->to = minto;
	} else {
		M->from = M->to = SUPPOS;
	}
	if (reg->to - reg2->to > ovlen) {
		R->from = (reg->from > reg2->to) ? reg->from : reg2->to;
		R->to = reg->to;
	} else {
		R->from = R->to = SUPPOS;
	}
}
shiftReg(Region *reg, SeqPos offset)
{
	reg->from += offset;
	reg->to += offset;
}
resetReg(Region *reg)
{
	reg->from = INFPOS;
	reg->to = SUPPOS;
}

transformSeqPos(SeqPos pos, Region *ali1, Region *ali2)
{
	if (regLen(ali1) <= 1) {
		printf("%d,%d\n",ali1->from,ali1->to);
		fprintf(stderr, "zero length ???:(%d,%d)\n",ali1->from,ali1->to);
		abort();
	}
	if (pos < ali1->from) {
		return pos - ali1->from + ali2->from;
	} else if (pos > ali1->to) {
		return pos - ali1->to + ali2->to;
	} else {
		return (pos - ali1->from) *
			(ali2->to - ali2->from) / (ali1->to - ali1->from)
			+ ali2->from;
	}
}
/* transformation of reg on a coordinate ali1 to a coordinate ali2 */
Region *transformReg2(Region *reg, Region *reg2, Region *ali1, Region *ali2)
{
	if (definedPos(reg->from)) {
		reg2->from = transformSeqPos(reg->from, ali1, ali2);
	} else {
		reg2->from = reg->from;
	}
	if (definedPos(reg->to)) {
		reg2->to = transformSeqPos(reg->to, ali1, ali2);
	} else {
		reg2->to = reg->to;
	}
	return reg2;
}
Region *transformReg(Region *reg, Region *ali1, Region *ali2)
{
	return transformReg2(reg, reg, ali1, ali2);
}
/* mapping subregion reg on ali1 to ali2: discarding outside of the alignment */
mapReg2(Region *reg, Region *newreg, Region *ali1, Region *ali2)
{
	Region tmp;
	int ov;
	transformReg2(reg, newreg, ali1, ali2);
	copyReg(&tmp, newreg);
	/* discard the region outside the destination region (ali2) */
	ov = ovlpReg(&tmp, ali2, newreg, NULL, 1);
	return ov;
}
mapReg(Region *reg, Region *ali1, Region *ali2)
{
	mapReg2(reg, reg, ali1, ali2);
}
printReg(Region *reg)
{
	if (reg) {
		printf("%d %d\n",reg->from, reg->to);
	} else {
		printf("(null)\n");
	}
}
