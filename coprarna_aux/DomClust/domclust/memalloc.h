/*
 * memalloc.h
 * Copyright (c) 2000-2006, Ikuo Uchiyama
 * All rights reserved.
 */

#ifndef MEMALLOC_H
#include <sys/types.h>
#include "plist.h"

typedef struct RegList {
	char *ptr;
	struct RegList *next, *prev;
} RegList;

typedef struct Object {
	RegList *first;
	RegList *last;
	size_t blksiz;
	size_t recno;
	size_t size;
	pList *free; 
} Alloc_Object;

Alloc_Object *init_alloc_object(), *init_alloc_object_with_freelist();
char *alloc_object();
char *get_objdata_idx();
char *memalloc();
char *memalloc_size();

#define MEMALLOC_H
#endif /* ifndef MEMALLOC_H */
