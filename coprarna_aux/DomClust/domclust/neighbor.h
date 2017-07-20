/*
 * DomClust: Hierarchical Clustering for Orthologous Domain Classification
 * Copyright (c) 2000-2006, Ikuo Uchiyama
 * All rights reserved.
 */

#ifndef _NEIGHBOR_H
#define _NEIGHBOR_H
typedef struct {
/*
        Node *node1, *node2;
*/
        char *node1, *node2;
        Count count;
        signed char dir;
        char status;
} Neighbor;
pList *nodeNeighbor_toposort();
pList *nodeNeighbor_gapsearch();
Node *getNeighborNode();
typedef enum {Ignore, AddCount} addElemMode;
typedef enum {chkIncluded, chkOverlap, chkCount} chkClustOvlpMode;
#endif
