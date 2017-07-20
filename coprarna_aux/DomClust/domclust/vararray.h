/*
 * vararray.h
 * Copyright (c) 2000-2006, Ikuo Uchiyama
 * All rights reserved.
 */

#ifndef _VARLIST_H_
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
        char *array;
        int num;
	int datasize;
        int allocsize;
} varArray;

typedef struct {
	varArray *array;
	int idx;
} arrayIter;

varArray *createVarArray(), *create2DVarArray();
arrayIter *createArrayIter();
char *getArrayItem(), *get2DArrayItem(), *getArrayIter();
char *getArrayIdx(), *get2DArrayIdx();
/*
char *getArrayIdx(varArray *vararray, const void *key,
		int (*func)(const void *, const void *));
*/

#define _VARLIST_H_
#endif
