/*
 * memalloc.c: a general memory allocation library
 * Copyright (c) 2000-2006, Ikuo Uchiyama
 * All rights reserved.
 */

#include <stdio.h>
#include <stdlib.h>
#include "memalloc.h"

Alloc_Object *init_alloc_object(size, blksiz)
	size_t size, blksiz;
{
	Alloc_Object *obj;

	if ((obj = (Alloc_Object *) malloc(sizeof(Alloc_Object))) == NULL) {
		return NULL;
	}
	obj->blksiz = blksiz;
	obj->recno = obj->blksiz;
	obj->size = size;
	obj->first = obj->last = NULL;
	obj->free = NULL;
	return obj;
}
Alloc_Object *init_alloc_object_with_freelist(size, blksiz)
	size_t size, blksiz;
{
	Alloc_Object *obj;
	obj = init_alloc_object(size, blksiz);
	obj->free = create_pList();
	return obj;
}

char *memalloc_size(obj, size)
	Alloc_Object *obj;
	int size;
{
	if (!obj) return(NULL);
	if (size == 1 && obj->free && numelemList(obj->free)) {
		char *ret = shiftList(obj->free);
		if (ret) return ret;
	}
	if (obj->recno + size > obj->blksiz) {
		RegList *reglist;
		char *alloc_obj;
		if ((reglist = (RegList *)malloc(sizeof(RegList))) == NULL) {
			return NULL;
		}
		if ((alloc_obj = (char *)malloc(obj->size * obj->blksiz)) == NULL) {
			return NULL;
		}
		reglist->next = NULL;
		reglist->ptr = alloc_obj;
		if (obj->first == NULL) {
			reglist->prev = NULL;
			obj->first = obj->last = reglist;
		} else {
			reglist->prev = obj->last;
			obj->last->next = reglist;
			obj->last = reglist;
		}
		obj->recno = size;		/* next point */
		return obj->last->ptr;
	} else {
		(obj->recno)+=size;		/* next point */
		return obj->last->ptr + obj->size * (obj->recno - size);
	}
}
char *memalloc(obj)
	Alloc_Object *obj;
{
	return memalloc_size(obj, 1);
}

char *memalloc_zero(obj)
	Alloc_Object *obj;
{
	char *p;
	p = memalloc(obj);
	bzero(p, obj->size);
	return p;
}

char *get_objdata_idx(obj, idx)
	Alloc_Object *obj;
	int idx;
{
	register RegList *reglist = obj->first;
	register int blkno = idx / obj->blksiz;
	int recno = idx % obj->blksiz;
	
	while (reglist && blkno--) {
		reglist = reglist->next;
	}
	if (reglist == NULL) return NULL;
	return reglist->ptr + obj->size * recno;
}

get_allocobj_num(obj)
	Alloc_Object *obj;
{
	RegList *p;
	int blknum = 0;
	for (p = obj->first; p; p = p->next) {
		blknum++;
	}
	return (blknum - 1) * obj->blksiz + obj->recno;
}
char **get_allocobj_array(obj)
	Alloc_Object *obj;
{
	int i, j;
	int allocnum = get_allocobj_num(obj);
	char **array;
	RegList *reg;

	if ((array = (char **) malloc(allocnum * sizeof(char *))) == NULL) {
		fprintf(stderr, "Can't alloc memory (%d)\n", allocnum);
		exit(1);
	}
	reg = obj->first;
	for (i = 0, j = 0; i < allocnum; i++, j++) {
		if (j >= obj->blksiz) {
			reg = reg->next;
			j = 0;
		}
		array[i] = reg->ptr + obj->size * j;
	}
	return array;
}

free_object(obj)
	Alloc_Object *obj;
{
	free_object_alldata(obj);
	free(obj);
	return 0;
}
clear_object(obj)
	Alloc_Object *obj;
{
	free_object_data(obj);
	obj->recno = obj->blksiz;
	obj->first = obj->last = NULL;
	return 0;
}
free_object_alldata(obj)
	Alloc_Object *obj;
{
	RegList *p, *next_p;
	for (p = obj->first; p; p = next_p) {
		free(p->ptr);
		next_p = p->next;
		free(p);
	}
	return 0;
}

free_lastobject(obj)
	Alloc_Object *obj;
{
	RegList *prev_p;
	if (obj->recno) {
		(obj->recno)--;
	} else {
		prev_p = obj->last->prev;
		free(obj->last->ptr);
		free(obj->last);
		obj->last = prev_p;
		obj->recno = obj->blksiz;
	}
}
free_object_data(obj, ptr)
	Alloc_Object *obj;
	char *ptr;
{
	if (ptr == 0){
		fprintf(stderr, "ERROR: attempting to free a null pointer\n");
		exit(10);
	}
	pushList(obj->free, ptr);
}

isAllocObj(obj, ptr)
	Alloc_Object *obj;
	char *ptr;
{
	RegList *p;
	for (p = obj->first; p; p = p->next) {
		if ((unsigned) ptr >= (unsigned) p->ptr &&
		   (unsigned) ptr < (unsigned) p->ptr + obj->blksiz * obj->size) {
			return 1;
		}
	}
	return 0;
}
