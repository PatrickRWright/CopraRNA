/*
 * DomClust: Hierarchical Clustering for Orthologous Domain Classification
 * Copyright (c) 2000-2006, Ikuo Uchiyama
 * All rights reserved.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "memalloc.h"
#include "vararray.h"
#include "namehash.h"
#include "hash.h"

NameHash *initNames0();

#ifdef DEBUG
char *addName();
main()
{
	int idx;
	int nid;
	NameHash *nhash = initNames(1000);
	addName(nhash, "HIJ",1);
	addName(nhash, "ABC",2);
	addName(nhash, "XYZ",3);
	addName(nhash, "EFG",4);
	nid = getNameID(nhash, "ABC");
	printf("%d\n",nid);
}
#endif

NameHash *initNames(int estim_size)
{
	return initNames0( estim_size * 17 - 11, estim_size * 6, estim_size );
}
NameHash *initNames0(int hashsize, int name_blksize, int datum_blksize)
{
	NameHash *nhash;

	if ((nhash = (NameHash *) malloc(sizeof(NameHash))) == NULL) {
		fprintf(stderr, "Can't alloc memory\n");
		exit(1);
	}
	nhash->nameobj = init_alloc_object(sizeof(char), name_blksize);
	nhash->datumobj = init_alloc_object(sizeof(NameID), datum_blksize);
	nhash->hash = Hcreate(hashsize);
	return nhash;
}

char *addName(NameHash *nhash, char *name, int id)
{
	char *ptr;
	NameID *tmpnameid;
	HENTRY hent;

	hent.key = name;
/*
	if (Hsearch(nhash->hash, &hent, FIND) == 1) {
		return NULL;
	}
*/

 	ptr = memalloc_size(nhash->nameobj, strlen(name)+1);
 	tmpnameid = (NameID *) memalloc(nhash->datumobj);
	strcpy(ptr, name);
	tmpnameid->name = ptr;
	tmpnameid->id = id;
	hent.key = ptr;
	hent.datum = (char*) tmpnameid;
	Hsearch(nhash->hash, &hent, ENTER);
	return ptr;
}

NameID *getNameIDStr(NameHash *nhash, char *name)
{
	int idx;
	NameID *nid;
	HENTRY hent;
/*
	tmpnid.name = name;
*/
	hent.key = name;
	if (! Hsearch(nhash->hash, &hent, FIND)) {
		return NULL;
	}
	nid = (NameID *) hent.datum;
	return nid;
}
char *getName(NameHash *nhash, char *name)
{
	NameID *nid = getNameIDStr(nhash, name);
	if (nid == NULL) {
		return NULL;
	} else {
		return nid->name;
	}
}
getNameID(NameHash *nhash, char *name)
{
	NameID *nid = getNameIDStr(nhash, name);
	if (nid == NULL) {
		return -1;
	} else {
		return nid->id;
	}
}
resetNameID(NameHash *nhash, char *name,  int id)
{
	NameID *nid = getNameIDStr(nhash, name);
	if (nid == NULL) {
		addName(nhash, name, id);
	} else {
		nid->id= id;
	}
}

isName(NameHash *nhash, char *ptr)
{
	if (! ptr) {
		return 0;
	}
	return isAllocObj(nhash->nameobj, ptr);
}
