/*
 * DomClust: Hierarchical Clustering for Orthologous Domain Classification
 * Copyright (c) 2000-2006, Ikuo Uchiyama
 * All rights reserved.
 */
#ifndef _BIN_H_
#define MAXBINS 5
typedef struct {
        int minval;
        int maxidx;
        int minidx;
	int size;
        int scale;
	int logscale;
        pList *bin;
	pList *bins[MAXBINS];
        pList uflw, oflw;
} Bin;

Bin *createBin();
listElem *addBinRetElem();
char *shiftMinData(), *popMaxData();
char *getBestData(Bin *bin, char maxflag);
#define _BIN_H_
#endif
