/*
 * plist.h
 * Copyright (c) 2000-2006, Ikuo Uchiyama
 * All rights reserved.
 */

#ifndef _PLIST_H_ 
typedef struct listElem {
	struct listElem *next;
	struct listElem *prev;
	char *datum;
} listElem;

typedef struct {
	listElem *first;
	listElem *last;
	int numelem;
} pList;

typedef struct {
	pList *plist;
	listElem *ptr;
	int dir;
} listIter;
typedef enum {LIST_FIND, LIST_INSERT} findElemMode;

listElem *allocListElem();
pList *create_pList(), *create_pList0(), *create_pList_array();
pList *dup_pList();
char *popList(), *shiftList();
listElem *getListElemIdx();
listIter *createListIter();
char *getListIter();
char *getListIdx();
#define _PLIST_H_
#endif
