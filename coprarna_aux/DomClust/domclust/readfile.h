/*
 * DomClust: Hierarchical Clustering for Orthologous Domain Classification
 * Copyright (c) 2000-2006, Ikuo Uchiyama
 * All rights reserved.
 */

#ifndef _READFILE_H_
#define _READFILE_H_
/*
#define HOMFILE_NAMELEN 64
#define SELFILE_NAMELEN 26
*/
/*
#define HOMFILE_NAMELEN 32
#define SELFILE_NAMELEN 32
*/
#define HOMFILE_NAMELEN 30
#define SELFILE_NAMELEN 30
#define BIN_MAGIC 197
#define MIN_HOMDATA 7

typedef struct {
	char name1[HOMFILE_NAMELEN];
	char name2[HOMFILE_NAMELEN];
	unsigned short from1, to1 /*, len1 */;
	unsigned short from2, to2 /*, len2 */;
	float ident;
	float eval;
	float score;
	float pam;
} HomData;

typedef struct {
	char name1[SELFILE_NAMELEN];
	char name2[SELFILE_NAMELEN];
#ifdef LARGE
	int from1, to1, from2, to2;
#else
	unsigned short from1, to1, from2, to2;
#endif
	float dist;
	float score;
#ifdef LARGE
	int dir;
#endif
} SelData;
typedef struct {
	FILE *fp;
	int (*read)();
	int rsize;
} SelFile;
int read_seldata_ascii();
int read_seldata_bin();
#endif
