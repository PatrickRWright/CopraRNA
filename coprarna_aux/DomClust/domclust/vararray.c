/*
 * vararray: variable size array
 * Copyright (c) 2000-2006, Ikuo Uchiyama
 * All rights reserved.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "vararray.h"

#ifdef ARRAYDEBUG
typedef struct {
        char *name;
        int id;
} NameID;
namecmp(const NameID *a, const NameID *b)
{
        return strcmp(a->name, b->name);
}

main()
{
	int idx;
	NameID data[10] = {{"ABC",1},{"NHK",2},{"TBS",3},{"NTV",4}};
	NameID *nid;
	varArray *namearray;
	int i;
	create2DVarArray(100, 100, sizeof(NameID));
	namearray = createVarArray(100, sizeof(NameID));
	for (i = 0; i < 4; i++) {
		addArray(namearray, &data[i]);
	}
	sortArray(namearray, namecmp);
	nid = getArrayIdx(namearray, &data[1], namecmp);
	printf("%s %d\n",nid->name, nid->id);
}
#endif

varArray *createVarArray(int initsize, int datasize)
{
	varArray *vararray;
	if ((vararray = (varArray *) malloc(sizeof(varArray))) == NULL) {
		fprintf(stderr, "Can't alloc vararray\n"); exit(1);
	}
	vararray->array = NULL;
	vararray->datasize = datasize;
	reallocArray(vararray, initsize);
	vararray->num = 0;
	return vararray;
}
varArray *create2DVarArray(int initsize1, int initsize2, int datasize)
{
	int i;
	varArray *vararray2 =
			createVarArray(initsize1, sizeof(varArray *));
	varArray *var;
	for (i = 0; i < initsize1; i++) {
		var = createVarArray(initsize2, datasize);
		addArray(vararray2, &var);
	}
	return vararray2;
}

add2DArray(varArray *vararray, int idx, char *datum)
{
	addArray(((varArray **)vararray->array)[idx], datum);
}
addArray(varArray *vararray, char *datum)
{
	setArray(vararray, vararray->num, datum);
	vararray->num++;
}
setArray(varArray *vararray, int idx, char *datum)
{
	while (idx >= vararray->allocsize) {
		reallocArray(vararray, (int) (vararray->allocsize * 1.8));
	}

	memcpy(vararray->array + vararray->datasize * idx,
		datum, vararray->datasize);
}
set2DArrayElem(varArray *vararray, int idx, varArray *newarray)
{
	varArray *var = ((varArray **)vararray->array)[idx];
	if (idx < vararray->num) {
		if (var) {
			freeArray(var);
		}
		memcpy((varArray **)vararray->array + idx, &newarray,
				vararray->datasize);
	}
}
set2DArray(varArray *vararray, int idx1, int idx2, char *datum)
{
	if (idx1 < vararray->num) {
		setArray(((varArray **)vararray->array)[idx1], idx2, datum);
	}
}

reallocArray(varArray *vararray, int size)
{
	if ((vararray->array = realloc(vararray->array,
				size * vararray->datasize)) == NULL) {
		fprintf(stderr, "Can't alloc array\n"); exit(1);
	}
	vararray->allocsize = size;
}
freeArray(varArray *vararray)
{
	free(vararray->array);
	free(vararray);
}
free2DArray(varArray *vararray)
{
	int idx;
	for (idx = 0; idx < vararray->num; idx++) {
		freeArray(((varArray **)vararray->array)[idx]);
	}
}
free2DArrayElem(varArray *vararray, int idx)
{
	freeArray(((varArray **)vararray->array)[idx]);
}
clearArray(varArray *vararray)
{
/*
	reallocArray(vararray, 0);
*/
	vararray->num = 0;
}
clear2DArrayElem(varArray *vararray, int idx)
{
	clearArray(((varArray **)vararray->array)[idx]);
}

printArray(varArray *vararray, int (*printfunc)(char *))
{
	int i;
	for (i = 0; i < vararray->num; i++) {
		printfunc(getArrayItem(vararray, i));
	}
}

char *getArrayItem(varArray *vararray, int idx)
{
	if (idx >= 0 && idx < vararray->num) {
		return vararray->array + vararray->datasize * idx;
	} else {
		return NULL;
	}
}
char *get2DArrayItem(varArray *vararray, int idx1, int idx2)
{
	if (idx1 >= 0 && idx1 < vararray->num) {
		return getArrayItem(((varArray **)vararray->array)[idx1], idx2);
	} else {
		return NULL;
	}
}

sortArray(varArray *vararray, int (*func)(const void *, const void *))
{
	qsort(vararray->array, vararray->num, vararray->datasize, func);
}
sort2DArray(varArray *vararray, int (*func)(const void *, const void *))
{
	int i;
	varArray **var = (varArray **) vararray->array;
	for (i = 0; i < vararray->num; i++) {
		qsort(var[i]->array, var[i]->num, var[i]->datasize, func);
	}
}

char *getArrayIdx(varArray *vararray, const void *key,
		int (*func)(const void *, const void *))
{
        return bsearch(key, vararray->array, vararray->num,
                                vararray->datasize, func);
}
char *get2DArrayIdx(varArray *vararray, int idx, const void *key,
		int (*func)(const void *, const void *))
{
	varArray **var = (varArray **) vararray->array;
        return bsearch(key, var[idx]->array, var[idx]->num,
                                var[idx]->datasize, func);
}

mergeArray(varArray *vararray1, varArray *vararray2,
		int (*func)(const void *, const void *))
{
	int i = 0, j = 0;
	char *p1, *p2;
	int num = 0;
	int datasize = vararray1->datasize;
	varArray *newarray;
	if (vararray2->datasize != datasize) {
		fprintf(stderr, "merge failed\n");
		return -1;
	}
	newarray = createVarArray(vararray1->allocsize, datasize);
	p1 = vararray1->array;
	p2 = vararray2->array;
	while (i < vararray1->num && j < vararray2->num) {
		if (func(p1, p2) == 0) {
			addArray(newarray, p1);
			p1 += datasize; p2 += datasize;
			i++; j++;
		} else if (func(p1, p2) < 0) {
			addArray(newarray, p1);
			p1 += datasize;
			i++;
		} else if (func(p1, p2) > 0) {
			addArray(newarray, p2);
			p2 += datasize;
			j++;
		}
	}
	while (i < vararray1->num) {
		if (func(p1, p2) < 0) {
			addArray(newarray, p1);
			p1 += datasize;
		}
		i++;
	}
	while (j < vararray2->num) {
		if (func(p1, p2) > 0) {
			addArray(newarray, p2);
			p2 += datasize;
		}
		j++;
	}
	freeArray(vararray1);
	freeArray(vararray2);
	return 0;
}
arraySize(varArray *array)
{
	return array->num;
}

arrayIter *createArrayIter(varArray *vararray)
{
	arrayIter *iter;
	if ((iter = (arrayIter *)malloc(sizeof(arrayIter)))== NULL) {
		fprintf(stderr, "Can't alloc arrayiter\n");
		exit(1);
	}
	setArrayIter(iter, vararray);
	return iter;
}
setArrayIter(arrayIter *iter, varArray *vararray)
{
	iter->array = vararray;
	iter->idx = 0;
}
initArrayIter(arrayIter *iter)
{
	iter->idx = 0;
}
char *getArrayIter(arrayIter *iter)
{
	return getArrayItem(iter->array, iter->idx++);
}
