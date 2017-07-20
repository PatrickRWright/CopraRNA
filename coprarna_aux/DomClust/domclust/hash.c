/*
 * hash.c: hash data structure
 * Copyright (c) 2000-2006, Ikuo Uchiyama
 * All rights reserved.
 */

#include <stdio.h>
#include <stdlib.h>
#include "hash.h"
#define PRIME 2047232111


static HENTRY **hashtable;
static HENTRY *allocrec(), *getrec();
static HashDataBlock *init_datablock();
int Hashsiz;

#define BLOCKSIZ 2500
#define MAXBLOCK 8000


/*
#define DEBUGMAIN
*/
#ifdef DEBUGMAIN
main(int argc, char **argv)
{
	Hash *h;
	HENTRY ent, *e;
	HashIter iter;
	char p[4][100];
	strcpy(p[0], "AAA");
	ent.datum = p[0];

	h = Hcreate(1001);
	ent.key = "ABCDEFGHIJK";
	Hsearch(h, &ent, ENTER);

	ent.key = "ABCDEFGHIJL";
	strcpy(p[1], "AAb");
	ent.datum = p[1];
	Hsearch(h, &ent, ENTER);

	ent.datum = p[1];
	strcpy(p[1], "BBB");
	printf("%s\n", ent.datum);
	ent.key = "ABCDEFGHIJ";
	Hsearch(h, &ent, FIND);
	printf("%s\n", ent.datum);
	ent.key = "ABCDEFGHIJK";
	Hsearch(h, &ent, FIND);
	printf("%s\n", ent.datum);
	initHashIter(&iter, h);
	while (e = nextHashIter(&iter)) {
		printf("%s,%s\n", e->key, e->datum);
	}
}
#endif

Hash *Hcreate(int hashsize)
{
	Hash *hashp;
	if ((hashp = (Hash *) malloc(sizeof(Hash))) == NULL) {
		fprintf(stderr, "Can't alloc hash\n");
		return 0;
	}
	if ((hashp->table = (HENTRY **) calloc(hashsize, sizeof(HENTRY *))) == NULL) {
		fprintf(stderr, "Can't alloc table\n");
		return 0;
	}
	hashp->hashdatablock = init_datablock();
	hashp->hashsize = hashsize;
	return hashp;
}

Hsearch(Hash *hashp, HENTRY *entry, ACTION action)
{
	HENTRY *p, *prevp;
	int i = 0;
	int hashsize = hashp->hashsize;
	int hashidx;
	register unsigned int keyi = 0;
	char *kp;
	int len = strlen(entry->key);


	for (kp = entry->key; *kp; kp++) {
		keyi *= 113; keyi += *kp; /* keyi %= PRIME; */
	}
	hashidx = hashf(keyi,i,hashsize);

	p = hashp->table[hashidx];
	while (p != NULL) {
		if (strcmp(p->key, entry->key) == 0) {
			entry->datum = p->datum;
			return 1;
		}
		if (++i >= hashsize) {
			fprintf(stderr, "Hash table overflows: %d\n", i);
			exit(1);
		}
		hashidx = hashf(keyi,i,hashsize);
		p = hashp->table[hashidx];
	}
	if (action == ENTER) {
		if ((hashp->table[hashidx] = allocrec(hashp->hashdatablock)) == NULL) {
			fprintf(stderr, "Can't alloc record\n");
			return -1;
		}
		hashp->table[hashidx]->key = entry->key;
		if (entry->datum)
			hashp->table[hashidx]->datum = entry->datum;
	}
	return 0;
}
HIsearch(Hash *hashp, HENTRY *entry, ACTION action)
{
	HENTRY *p;
	int i = 0;
	int hashsize = hashp->hashsize;
	int hashidx;
	unsigned int keyi = (unsigned int) entry->key;

	hashidx = hashf(keyi,i,hashsize);

	p = hashp->table[hashidx];
	while (p != NULL) {
		if (p->key == entry->key) {
			entry->datum = p->datum;
			return 1;
		}
		if (++i >= hashsize) {
			fprintf(stderr, "Hash table overflows: %d\n", i);
			exit(1);
		}
		hashidx = hashf(keyi,i,hashsize);
		p = hashp->table[hashidx];
	}
	if (action == ENTER) {
		if ((hashp->table[hashidx] = allocrec(hashp->hashdatablock)) == NULL) {
			fprintf(stderr, "Can't alloc record\n");
			return -1;
		}
		hashp->table[hashidx]->key = entry->key;
		if (entry->datum)
			hashp->table[hashidx]->datum = entry->datum;
	}
	return 0;
}

Hdestroy(Hash *hashp)
{
	free(hashp->table);
	destroy_datablock(hashp->hashdatablock);
	free(hashp);
}

initHashIter(HashIter *iter, Hash *hashp)
{
	iter->hash = hashp;
	iter->idx = 0;
}
HENTRY *nextHashIter(HashIter *iter)
{
	return getrec(iter->hash->hashdatablock, iter->idx++);
}
numelemHash(Hash *hash)
{
	HashDataBlock *dblock = hash->hashdatablock;
	if (dblock->blksiz < 0) {
		return 0;
	} else {
		return dblock->blksiz * BLOCKSIZ + dblock->recnum + 1;
	}
}


static HashDataBlock *init_datablock()
{
	HashDataBlock *dblock;
	if ((dblock = (HashDataBlock *) malloc(sizeof(HashDataBlock))) == NULL) {
		fprintf(stderr, "Can't alloc hash datablock\n");
		exit(1);
	}
	if ((dblock->datablock = (HENTRY **)
			calloc(BLOCKSIZ, sizeof(HENTRY *))) == NULL) {
		fprintf(stderr, "Can't alloc hash datablock\n");
		exit(1);
	}
	dblock->blksiz = -1;
	dblock->recnum = BLOCKSIZ;
	return dblock;
}
destroy_datablock(HashDataBlock *dblock)
{
	int i;
	for (i = 0; i < dblock->blksiz; i++) {
		free(dblock->datablock[i]);
	}
	free(dblock->datablock);
	free(dblock);
}
static HENTRY *allocrec(HashDataBlock *dblock)
{
	if (++(dblock->recnum) >= BLOCKSIZ) {
		if (++(dblock->blksiz) >= MAXBLOCK) {
			fprintf(stderr, "Hash Block overflows\n");
			exit(1);
		}
		if ((dblock->datablock[dblock->blksiz] = (HENTRY *) malloc(sizeof(HENTRY) * BLOCKSIZ)) == NULL) {
			fprintf(stderr, "Can't alloc memory (Hash)\n");
			exit(1);
		}
		dblock->recnum = 0;
	}
	return &(dblock->datablock[dblock->blksiz][dblock->recnum]);
}

static HENTRY *getrec(HashDataBlock *dblock, int idx)
{
	if (idx > dblock->blksiz * BLOCKSIZ + dblock->recnum) {
		return NULL;
	}
	return &(dblock->datablock[idx / BLOCKSIZ][idx % BLOCKSIZ]);
}
