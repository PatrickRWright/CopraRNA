/*
 * plist.c: linked list data structure
 * Copyright (c) 2000-2006, Ikuo Uchiyama
 * All rights reserved.
 */

#include <stdio.h>
#include <stdlib.h>
#include<string.h>
#include "memalloc.h"
#include "plist.h"

#define LISTBLKSIZ 10000
#define FREED ((char *) -1)

/********* FOR DEBUG **************/
/* #define _DEBUGMAIN_PLIST */
#ifdef _DEBUGMAIN_PLIST
#include <string.h>
char *test[5] = {"DEF", "ABC","JK","GHI","LM"};

print_datum(char *datum)
{
	printf(" %s", datum);
}
mystrcmp(char **str1, char **str2)
{
	return strcmp(*str1, *str2);
}

main()
{
	int i;
	pList *plist = create_pList();
	listElem *e;
	for (i = 0; i <= 4; i++) {
		printf("Push: %s\n", test[i]);
		pushList(plist,test[i]);
		printStrList(plist);
	}
	e = getListElemIdx(plist, 3);
	sortList(plist, mystrcmp);
	printf("Sort:\n");
	printStrList(plist);
	e = getListElemIdx(plist, 3);
	printf("Del: %s\n", e->datum);
	delelemList(plist, e);
	printStrList(plist);
	for (i = 1; i <= 3; i++) {
		printf("Pop: %s\n", popList(plist));
		printStrList(plist);
		printf("Unshift: %s\n", test[i-1]);
		unshiftList(plist,test[i-1]);
		printStrList(plist);
	}
	clearList(plist);
	printf("Clear\n"); printStrList(plist);
	for (i = 0; i <= 2; i++) {
		pushList(plist,test[i]);
		printf("Push: %s\n", test[i]); printStrList(plist);
	}
	for (i = 0; i <= 3; i++) {
		char *s = shiftList(plist);
		if (s) {
			printf("Shift: %s\n",s);
		}
		printStrList(plist);
	}
}
printStrList(pList *plist)
{
	printf("[");
	execList(plist, print_datum);
	printf(" ]\n");
}
#endif
/***********************/

/* global memory allocator */
static Alloc_Object *plistElemObj;
/* global free list */
static listElem *_freelist;


pList *create_pList()
{
	return create_pList0(LISTBLKSIZ, 1);
}
pList *create_pList_array(int nel)
{
	return create_pList0(LISTBLKSIZ, nel);
}
pList *create_pList0(int blksiz, int nel)
{
	pList *plist;
	Alloc_Object *eobj = NULL;
	int i;

	if (! plistElemObj) {
		plistElemObj = 
			init_alloc_object(sizeof(listElem), blksiz);
	}

	if ((plist = (pList *) malloc(sizeof(pList)*nel)) == NULL) {
		return NULL;
	}
	for (i = 0; i < nel; i++) {
		init_pList(&plist[i], blksiz);
	}
	return plist;
}

init_pList(pList *plist)
{
	plist->first = plist->last = NULL;
	plist->numelem = 0;
}
pList *init_pList_array(pList *plist, register int nel)
{
	while (--nel>=0) {
		init_pList(plist + nel);
	}
}

pushList(pList *plist, char *datum)
{
	listElem *elem;
	elem= allocListElem(plist);
	elem->datum = datum;
	elem->next = NULL;
	if (plist->first == NULL) {
		plist->first = plist->last = elem;
		elem->prev = NULL;
	} else {
		plist->last->next = elem;
		elem->prev = plist->last;
		plist->last = elem;
	}
	plist->numelem++;
}

char *popList(pList *plist)
{
	if (plist && plist->last) {
		listElem *elem = plist->last;
		char *ret = elem->datum;
		if (plist->last->prev) {
			plist->last->prev->next = NULL;
		} else {
			plist->first = NULL;
		}
		plist->last = plist->last->prev;
		freeElem(elem);
		plist->numelem--;
		return ret;
	} else {
		return NULL;
	}
}

unshiftList(pList *plist, char *datum)
{
	listElem *elem = allocListElem(plist);
	elem->datum = datum;
	elem->prev = NULL;
	if (plist->last == NULL) {
		plist->first = plist->last = elem;
		elem->next = NULL;
	} else {
		plist->first->prev = elem;
		elem->next = plist->first;
		plist->first = elem;
	}
	plist->numelem++;
}

char *shiftList(pList *plist)
{
	if (plist->first) {
		listElem *elem = plist->first;
		char *ret = elem->datum;
		if (plist->first->next) {
			plist->first->next->prev = NULL;
		} else {
			plist->last = NULL;
		}
		plist->first = plist->first->next;
		freeElem(elem);
		plist->numelem--;
		return ret;
	} else {
		return NULL;
	}
}
pushListAll(pList *list1, pList *list2)
{
	listIter iter;
	char *datum;
	setListIter(&iter, list2, 1);
	while (datum = getListIter(&iter)) {
		pushList(list1, datum);
	}
}
concatList(pList *list1, pList *list2)
{
	if (numelemList(list1)==0) {
		list1->first = list2->first;
		list1->last = list2->last;
	} else {
		list1->last->next = list2->first;
		list2->first->prev = list1->last;
		list1->last = list2->last;
	}
	list1->numelem += list2->numelem;
	list2->first = list2->last = NULL;
	list2->numelem = 0;
}

char *getListIdx(pList *plist, int idx)
{
	listElem *e = getListElemIdx(plist, idx);
	if ( ! e)  {
		return NULL;
	} else {
		return e->datum;
	}
}
sortList(pList *plist, int (*func)(const void *, const void *))
{
	char **array;
	char *datum;
	int i = 0;
	int numelem = plist->numelem;
	if (numelem <=1){
		return;
	}
	if ((array = (char **) malloc(sizeof(char*) * numelem)) == NULL) {
		fprintf(stderr, "Can't allocate memory\n");
		exit(1);
	}
	while ( datum = shiftList(plist) ) {
		array[i++] = datum;
	}
	qsort(array, numelem, sizeof(char *), func);
	for (i = 0; i < numelem; i++) {
		pushList(plist, array[i]);
	}
	free(array);
}

addElemSortedList(pList *plist, char *datum,
	int (*func)(const void *, const void *))
{
	/* insert mode */
	findElemSortedList(plist, datum, func, NULL, LIST_INSERT);
}
findElemSortedList(pList *plist, char *datum,
	int (*func)(const void *, const void *), char **ret,
	findElemMode mode)
{
	listElem *e;
	int cmp;
	listElem *new;

	if (! plist->first) {
		if (mode == LIST_INSERT) {
			new = allocListElem(plist);
			new->datum = datum;
			plist->first = plist->last = new;
			plist->numelem++;
		}
		return 0;
	}
	for (e = plist->first; e; e = e->next) {
		cmp = func(e->datum, datum);
		if (cmp == 0) {
			if (ret) {
				/* return the found object */
				*ret = e->datum;
			}
			return 1;
		} else if (cmp >= 0) {
			if (mode == LIST_INSERT) {
				new = allocListElem(plist);
				new->datum = datum;
				new->prev = e->prev;
				new->next = e;
				if (e->prev) {
					e->prev->next = new;
				} else {
					plist->first = new;
				}
				e->prev = new;
				plist->numelem++;
			}
			return 0;
		}
	}
	/* last */
	if (mode == LIST_INSERT) {
		new = allocListElem(plist);
		new->datum = datum;
		new->prev = plist->last;
		new->next = NULL;
		plist->last->next = new;
		plist->last = new;
		plist->numelem++;
	}
	return 0;
}
addElemSortedListInv(pList *plist, char *datum,
	int (*func)(const void *, const void *))
{
	/* insert mode */
	findElemSortedListInv(plist, datum, func, NULL, LIST_INSERT);
}
findElemSortedListInv(pList *plist, char *datum,
	int (*func)(const void *, const void *), char **ret,
	findElemMode mode)
{
	listElem *e;
	int cmp;
	listElem *new;

	if (! plist->first) {
		if (mode == LIST_INSERT) {
			new = allocListElem(plist);
			new->datum = datum;
			plist->first = plist->last = new;
			plist->numelem++;
		}
		return 0;
	}
	for (e = plist->last; e; e = e->prev) {
		cmp = func(e->datum, datum);
		if (cmp == 0) {
			if (ret) {
				/* return the found object */
				*ret = e->datum;
			}
			return 1;
		} else if (cmp <= 0) {
			if (mode == LIST_INSERT) {
				new = allocListElem(plist);
				new->datum = datum;
				new->next = e->next;
				new->prev = e;
				if (e->next) {
					e->next->prev = new;
				} else {
					plist->last = new;
				}
				e->next = new;
				plist->numelem++;
			}
			return 0;
		}
	}
	/* last */
	if (mode == LIST_INSERT) {
		new = allocListElem(plist);
		new->datum = datum;
		new->next = plist->first;
		new->prev = NULL;
		plist->first->prev = new;
		plist->first = new;
		plist->numelem++;
	}
	return 0;
}
addElemPrevList(pList *plist, listElem *e, char *datum)
{
	listElem *new = allocListElem(plist);
	new->datum = datum;
	new->next = e;
	new->prev = e->prev;
	if (e->prev) {
		e->prev->next = new;
	} else {
		plist->first = new;
	}
	e->prev = new;
	plist->numelem++;
}
addElemNextList(pList *plist, listElem *e, char *datum)
{
	listElem *new = allocListElem(plist);
	new->datum = datum;

	new->prev = e;
	new->next = e->next;
	if (e->next) {
		e->next->prev = new;
	} else {
		plist->last = new;
	}
	e->next = new;
	plist->numelem++;
}
delelemListIdx(pList *plist, int idx)
{
	listElem *e = getListElemIdx(plist, idx);
	delelemList(plist, e);
}
mergelist(pList *plist1, pList *plist2,
		int (*func)(const char *datum1, const char *datum2))
/* plist1 and plist2 must be ordered */
{
	listElem *e1 = plist1->first, *e2 = plist2->first;
	int cmp;
	while (e1 && e2) {
		cmp = func(e1->datum,e2->datum);
		if (cmp == 0) {
			e1 = e1->next;
			e2 = e2->next;
		} else if (cmp < 0) {
			e1 = e1->next;
		} else {
			addElemPrevList(plist1, e1, e2->datum);
			e2 = e2->next;
		}
	}
	while (e2) {
		pushList(plist1, e2->datum);
		e2 = e2->next;
	}
}

execList(pList *plist, int (*func)(const char *datum))
{
	listElem *e;
	if (! plist) return 0;
	for (e = plist->first; e; e = e->next) {
		func(e->datum);
	}
}
printList(pList *plist, int (*func)(const char *datum))
{
	listElem *e;
	if (! plist) return 0;
	printf("[ ");
	for (e = plist->first; e; e = e->next) {
		func(e->datum);
		if (e->next) putchar(',');
	}
	printf(" ]");
}
numelemList(pList *plist)
{
	if (! plist) {
		return 0;
	}
	return plist->numelem;
}

freeElem(listElem *e)
{
	e->next = _freelist;
	e->datum = (char *) FREED;
	_freelist = e;
}
freeList(pList *plist)
{
	if (plist->first) {
		plist->last->next = _freelist;
		_freelist = plist->first;
		init_pList(plist);
	}
}
free_pList(pList *plist)
{
	freeList(plist);
	free(plist);
}
free_pList_array(pList *plist, int nel)
{
	int i;
	for (i = 0; i < nel; i++) {
		freeList(&plist[i]);
	}
	free(plist);
}

clearList_array(pList *plist, int nel)
{
	while (--nel>=0) {
		clearList(plist + nel);
	}
}
clearList(pList *plist)
{
	freeList(plist);
	plist->first = plist->last = NULL;
	plist->numelem = 0;
}
delelemList(pList *plist, listElem *elem)
{
	if (! elem) {
		return 0;
	} else if (elem->datum == FREED) {
fprintf(stderr, "delelem: warning: already freed: %d\n", elem);
		return -1;
	}
	if (elem->next) {
		elem->next->prev = elem->prev;
	} else {
		plist->last = elem->prev;
	}
	if (elem->prev) {
		elem->prev->next = elem->next;
	} else {
		plist->first = elem->next;
	}

	freeElem(elem);
	plist->numelem--;
	return 0;
}
listElem *getListElemIdx(pList *plist, register int idx)
{
	register listElem *e;
	if (plist->numelem == 0) {
		return NULL;
	}
	if (idx >= 0) {
		for (e = plist->first; e; e = e->next) {
			if (idx-- <= 0) break;
		}
	} else {
		for (e = plist->last; e; e = e->prev) {
			if (++idx >= 0) break;
		}
	}
	return e;
}

listElem *allocListElem(pList *plist)
{
	listElem *ret;
	if (_freelist) {
		/* reuse a free object if available */
		ret = _freelist;
		_freelist = _freelist->next;
	} else {
		ret = ((listElem *) memalloc(plistElemObj));
		if (! ret) {
			fprintf(stderr, "Can't alloc memory for plist\n");
			exit(1);
		}
	}
	ret->next = ret->prev = NULL;
	return ret;
}

listIter *createListIter(pList *plist, int dir)
{
	listIter *iter;
	if ((iter = (listIter*) malloc(sizeof(listIter))) == NULL) {
		fprintf(stderr, "Can't alloc memory for listiter\n");
		exit(1);
	}
	setListIter(iter, plist, dir);
	return iter;
}
setListIter(listIter *iter, pList *plist, int dir)
{
	iter->dir = dir;
	iter->plist = plist;
	if (dir < 0) {
		iter->ptr = plist->last;
	} else {
		iter->ptr = plist->first;
	}
}
setListIterRev(listIter *iter)
{
	iter->dir = - iter->dir;
	if (iter->ptr) {
		getListIter(iter);
	} else if (iter->dir > 0) {
		iter->ptr = iter->plist->first;
	} else {
		iter->ptr = iter->plist->last;
	}
}
copyListIter(listIter *iter1, listIter *iter2)
{
	iter2->plist = iter1->plist;
	iter2->ptr = iter1->ptr;
	iter2->dir = iter1->dir;
}
char *getListIter(listIter *itr)
{
	char *p;
	if (! itr || ! itr->ptr) {
		return NULL;
	}
	p = itr->ptr->datum;
	if (itr->dir < 0) {
		itr->ptr = itr->ptr->prev;
	} else {
		itr->ptr = itr->ptr->next;
	}
	return p;
}
freeListIter(listIter *itr)
{
	free(itr);
}
delCurrElem(listIter *itr)
{
	/** itr->ptr points the next element to be got !! **/
	if (itr->dir > 0) {
		if (! itr->ptr) {
			delelemList(itr->plist, itr->plist->last);
		} else if (itr->ptr->prev) {
			delelemList(itr->plist, itr->ptr->prev);
		}
	} else {
		if (! itr->ptr) {
			delelemList(itr->plist, itr->plist->first);
		} else if (itr->ptr->next) {
			delelemList(itr->plist, itr->ptr->next);
		}
	}
}
insAfterCurrElem(listIter *itr, char *datum)
{
	/** itr->ptr points the next element of the target one !! **/
	if (itr->dir > 0) {
		if (! itr->ptr) {
			pushList(itr->plist, datum);
		} else {
			addElemPrevList(itr->plist, itr->ptr, datum);
		}
	} else {
		if (! itr->ptr) {
			unshiftList(itr->plist, datum);
		} else {
			addElemNextList(itr->plist, itr->ptr, datum);
		}
	}
}
insBeforeCurrElem(listIter *itr, char *datum)
{
	listElem *curr;
	/** itr->ptr points the next element of the target one !! **/
	if (itr->dir > 0) {
		if (! itr->ptr) {
			curr = itr->plist->last;
		} else {
			curr = itr->ptr->prev;
		}
		if (! curr) {
			unshiftList(itr->plist, datum);
		} else {
			addElemPrevList(itr->plist, curr, datum);
		}
	} else {
		if (! itr->ptr) {
			curr = itr->plist->first;
		} else {
			curr = itr->ptr->next;
		}
		if (! curr) {
			pushList(itr->plist, datum);
		} else {
			addElemNextList(itr->plist, curr, datum);
		}
	}
}
newrefList(pList *p1, pList *p2)
{
	if (p2->first) {
		p2->last->next = _freelist; 
		_freelist = p2->first;
	}
	p2->first = p1->first;
	p2->last = p1->last;
	p2->numelem = p1->numelem;
}
moveList(pList *p1, pList *p2)
{
	newrefList(p1, p2);
	p1->first = p1->last = NULL;
	p1->numelem = 0;
}
pList *dup_pList(pList *plist)
{
	listElem *e;
	pList *newlist;
	newlist = create_pList();
	for (e = plist->first; e; e = e->next) {
		pushList(newlist, e->datum);
	}
	return newlist;
}

checkElem(pList *plist, listElem *e0)
{
	listElem *e;
	if (! plist) return 0;
	for (e = plist->first; e; e = e->next) {
		if (e0 == e) {
			return 1;
		}
	}
	return 0;
}
findList(pList *plist, char *ptr)
{
	listElem *e;
	if (! plist) return 0;
	for (e = plist->first; e; e = e->next) {
		if (e->datum == ptr) {
			return 1;
		}
	}
	return 0;
}

isFree(listElem *e0)
{
	if (e0->datum == FREED) {
		return 1;
	} else {
		return 0;
	}
}
isFreeDeep(listElem *e0)
{
	listElem *e;
	if (e0->datum == FREED) {
		return 1;
	}
	for (e = _freelist; e; e = e->next) {
		if (e == e0) {
			return 2;
		}
	} 
	return 0;
}


