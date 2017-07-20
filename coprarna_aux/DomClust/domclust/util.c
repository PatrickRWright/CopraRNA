/*
 * DomClust: Hierarchical Clustering for Orthologous Domain Classification
 * Copyright (c) 2000-2006, Ikuo Uchiyama
 * All rights reserved.
 */

#include <stdio.h>
#include "domclust.h"
allocError(char *msg)
{
	fprintf(stderr, "Can't alloc memory: %s\n", msg);
	exit(1);
}

cmpr_int(int *i, int *j)
{
	return *i - *j;
}
cmpr_pos(SeqPos *i, SeqPos *j)
{
	return *i - *j;
}
