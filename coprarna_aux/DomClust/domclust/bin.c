/*
 * DomClust: Hierarchical Clustering for Orthologous Domain Classification
 * Copyright (c) 2000-2006, Ikuo Uchiyama
 * All rights reserved.
 */
#include <stdio.h>
#include <math.h>
#include "memalloc.h"
#include "plist.h"
#include "bin.h"
#include "domclust.h"

#ifdef DEBUG
main()
{
	Bin *bin = createBin(0,50,10);
	addBin(bin, 33.4, "ABC");
	addBin(bin, 5.4, "EFG");
	addBin(bin, 24.2, "HIJ");
	addBin(bin, 15.2, "KL");
	addBin(bin, 45.2, "MN");
	printf(">%s\n", shiftMinData(bin));
	printf(">%s\n", shiftMinData(bin));
	printf(">%s\n", shiftMinData(bin));
	printf(">%s\n", shiftMinData(bin));
	printf(">%s\n", shiftMinData(bin));
}
#endif

Bin *createBin(int min, int max, int scale, int logscale)
{
	Bin *bin;
	int i;

	if (min > max) {
		fprintf(stderr, "internal error: createBin: min>max???\n");
		exit(1);
	}
	if ((bin = (Bin *) malloc(sizeof(Bin))) == NULL) {
		fprintf(stderr, "Can't allocate bin\n");
		exit(1);
	}

	if (logscale) {
		if (min <= 0) min = 1;
		bin->minval = min;
		bin->size = bin->maxidx = (int)(log((double)max/min)*scale);
	} else {
		bin->minval = min;
		bin->size = bin->maxidx = (max-min)*scale;
	}
	bin->minidx = 0;
	bin->scale = scale;
	bin->logscale = logscale;

	if ((bin->bin = (pList *) create_pList_array(bin->size + 1)) == NULL) {
		fprintf(stderr, "size: %d (%d,%d)\n",bin->size, min,max);
		allocError("bin->bin");
	}
	init_pList(&(bin->uflw));
	init_pList(&(bin->oflw));
	return bin;
}

addBin(Bin *bin, double value, char *datum)
{
	addBinRetElem(bin, value, datum);
}

int getBinIdx(Bin *bin, double value)
{
	int idx;
	if (bin->logscale) {
		if ( (value - (double) bin->minval) < 1 ) {
			idx = 0;
		} else {
			idx = (int) (log(value / (double) bin->minval)
					* bin->scale);
		}
	} else {
		idx = (int) ((value - (double) bin->minval) * bin->scale);
	}
	return idx;
}
listElem *addBinRetElem(Bin *bin, double value, char *datum, pList *tmp)
{
	int idx;
	idx = getBinIdx(bin, value);


	if (idx < 0) {
		pushList(&(bin->uflw), datum);
		fprintf(stderr, "Warning: score bucket underflows: %d,%lf\n", idx, value);
		return getListElemIdx(&(bin->uflw), -1);
	} else if (idx > bin->size) {
		pushList(&(bin->oflw), datum);
		fprintf(stderr, "Warning: score bucket overflows: %d,%lf\n", idx,value);
		return getListElemIdx(&(bin->oflw), -1);
	} else {
		pushList(&(bin->bin[idx]), datum);
		return getListElemIdx(&(bin->bin[idx]), -1);
	}
}

char *shiftMinData(Bin *bin)
{
	char *ret = NULL;
	while ( (ret = shiftList(&(bin->bin[bin->minidx]))) == NULL ) {
		if (++bin->minidx >= bin->maxidx) {
			break;
		}
	}
	return ret;
}
char *popMaxData(Bin *bin)
{
	char *ret = NULL;
	while ( (ret = popList(&(bin->bin[bin->maxidx]))) == NULL ) {
		if (--bin->maxidx <= bin->minidx) {
			break;
		}
	}
	return ret;
}
char *getBestData(Bin *bin, char maxflag)
{
	return ((maxflag) ? popMaxData(bin) : shiftMinData(bin));
}

moveBinData(Bin *bin, int idx1, int idx2, listElem *e)
{
}

delBinData(Bin *bin, double value, listElem *e)
{
	int idx = getBinIdx(bin, value);
	if (idx < 0) {
		delelemList(&(bin->uflw), e);
	} else if (idx > bin->size) {
		delelemList(&(bin->oflw), e);
	} else {
		delelemList(bin->bin + idx, e);
	}
	if (idx == bin->minidx) {
		while (numelemList(bin->bin + bin->minidx) == 0) {
			bin->minidx++;
		}
	}
	if (idx == bin->maxidx) {
		while (numelemList(bin->bin + bin->maxidx) == 0) {
			bin->maxidx--;
		}
	}
}

freeBin(Bin *bin)
{
	free_pList_array(bin->bin, bin->size + 1);
	free(bin);
}
