/*
 * DomClust: Hierarchical Clustering for Orthologous Domain Classification
 * Copyright (c) 2000-2006, Ikuo Uchiyama
 * All rights reserved.
 */

#ifndef _NAMEHASH_H_
#include "hash.h"

#define INIT_NAMLSTSIZ 80000
#define NAMEREG_BLKSIZ 100000

#define hashf(x,i,hashsize) \
        (((x) % hashsize + (1 + (x) % (hashsize-2)) * (i)) % hashsize)

typedef struct {
        char *name;
        int id;
} NameID;

typedef struct {
        Alloc_Object *nameobj, *datumobj;
        Hash *hash;
} NameHash;
NameHash *initNames();
char *addName();
char *getName();

#define _NAMEHASH_H_
#endif
