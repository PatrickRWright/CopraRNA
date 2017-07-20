/*
 * DomClust: Hierarchical Clustering for Orthologous Domain Classification
 * Copyright (c) 2000-2006, Ikuo Uchiyama
 * All rights reserved.
 */
#ifndef _DOMCLUST_
#define _DOMCLUST_
#include "memalloc.h"
#include "vararray.h"
#include "plist.h"
#include "bin.h"
#include "seqreg.h"
#include "namehash.h"
#include "spec.h"

#define DOMCLUST_VERSION "1.0.1"

#define MAXDOM 256
#define todigit(c) (c - '0')
#define BETTER(a,b) (Better((double) a,(double) b,Opt.sim))
#define DISTRATIO(a,b) ( DistRatio((double) a,(double) b, Opt.sim) )
#define MEASURE(e) (Opt.sim ? (e)->score : (e)->dist)

#define BADVALUE 99999999
#define DISTSCALE 20
#define MAXDIST 270
#define UNDEFINED 99999999
#define MISSRATIO 0.95
#define MINSCORE 60
#define MAXSCORE 50000
#define MAXUCHR 255
#define MAXNBR 25000

#define INIT_NODENUM 200000
#define NODE_BLKSIZ 40000
#define EDGE_BLKSIZ 120000
#define INT_BLKSIZ 100000
#define AVE_REL 50

#define NAMELEN 200
#define SPNAMELEN 10

#define otherNodeNoCheck(e,n) (e->node1->id == n ? e->node2 : e->node1)
#define otherNodeCheck(e,n) \
  ( e->node1->id == n ? e->node2 : ((e->node2->id == n) ? e->node1 : NULL) )
#define otherNode(e,n) otherNodeCheck(e,n)
#define otherNodeID(e,n) otherNodeNoCheck(e,n)->id
#define otherNodeAli(e,n) (e->node1->id == n ? e->ali2 : e->ali1);

#define NODE_CLEARED ((NodeFlag) 1)
#define NODE_INT1 ((NodeFlag) 2)
#define NODE_INT2 ((NodeFlag) 4)
#define NODE_INTERNAL (NODE_INT1|NODE_INT2)
#define NODE_DELETED ((NodeFlag) 8)
#define NODE_MERGED ((NodeFlag) 16)	/* the cluster is included in the adjacent cluster */
#define NODE_MERGED2 ((NodeFlag) 32) /* sequences once split are merged during the clustering */
#define NODE_TMPMARK ((NodeFlag) 128)
#define NODE_TMPMARK2 ((NodeFlag) 64)
#define NODE_DUPMARK ((NodeFlag) 256)

#define EDGE_CLEARED 1
#define EDGE_SELECTED 2
#define EDGE_TMPMARK 128
#define EDGE_TMPMARK2 64

enum outStyle {
	NORMALOUT, TREEOUT, NEWICK, NEWICK_DIST, HIER_DETAILOUT,
		CLUSTTAB, GRAPHOUT, NEIGHBOR, DOMAINOUT, TABOUT, DUMP
};
typedef enum {NORMAL_EDGE, MULTI_EDGE, SELF_EDGE, SELF_MATCH} EdgeType;
typedef enum {B_BREAK, B_ALIREG} BoundaryOutput;
typedef enum {D_MATCH, D_LEFT1, D_LEFT2, D_RIGHT1, D_RIGHT2,
                D_MATCH1, D_MATCH2, D_THIRD} DomIdx;
typedef enum {DBG_basic=1, DBG_cluster=2, DBG_nbrconst=4, DBG_nbrrestore=8}
	DBG_Mode;

typedef unsigned char DomID;
typedef unsigned char ConnCount;
typedef unsigned short Count;
typedef unsigned int NodeID;
typedef unsigned int EdgeID;
typedef unsigned int ClustID;
/*
typedef double Dist; 
*/
typedef float Dist; 
typedef unsigned short NodeFlag; 
typedef unsigned char EdgeFlag; 

/* SeqRegion */
/*
#define regLen(reg) ((reg)->to - (reg)->from + 1)
#define checkOvlpReg(reg1,reg2) \
        ((reg1)->from < (reg2)->to && (reg2)-> from < (reg1)->to)

typedef struct {
	SeqPos from, to;
} Region;
*/

/* Graph */
typedef struct Node {
	char *name;
	NodeID id;
	NodeFlag flag;
	signed char dir;
	pList *domains;
	SeqPos from, len;
	SeqPos meanlen;
	SeqPos minlen;
	SeqPos totlen;
	struct Node *parent, *parentL, *parentR;
	pList *left, *right;
	Region consreg;		/* conserved region within the seq */

	Region brk, newreg;

	struct Node *parentM;
	Region brk2, newreg2;	/* for self match */

	struct Edge *child;
	Count cnt;
	Dist meandist;
	specFlag spflag;
} Node;
typedef struct {
	Alloc_Object *nodeobj, *intobj, *lenobj, *idobj;
	int leafnum;
	int total_nodenum;	/* leaf nodes + internal nodes */
	int nodenum;		/* extra nodes may be added */
} NodeSet;
typedef struct {
	NodeSet *nodes;
	int i;
	int dir;
	int end;
} NodeSetIter;

typedef struct Edge {
	Node *node1, *node2;
	EdgeID id;
	StrRegion *ali1, *ali2;
	Dist dist, score;
	listElem *binelem;
	listElem *ordelem1, *ordelem2;
#ifndef NONEIGHBOR
	pList *left, *right;
#endif
	EdgeFlag flag;
	signed char dir;	/* reverse match on nucleotide comparison */
	ConnCount connect;	/* number of true connection */
} Edge;
typedef struct {
	Alloc_Object *edgeobj, *idobj, *lenobj, *regobj;
	int edgenum;
	int *bestdist;
	Edge **bestedge;
	varArray *tmplist;
	pList *ordlist;
	int *nodeindex;
	Bin *bin;

/*
	int (*better)(double a,double b);
	char sim;
*/
} EdgeSet;

typedef enum {
	n3flag_skip = 1
} Node3Flag;
typedef struct {
	NodeID n3;
	Edge *e1, *e2;
	Region *ali13, *ali31, *ali23, *ali32;
	Region ali_e11, ali_e12, ali_e21, ali_e22;
	char flag;
} Node3List;

typedef struct {
	NodeID id1, id2;
	Edge *edge;
	listElem *ordelem;
} pairList;

typedef struct {
	NodeID id;
	ClustID clid;
} Cluster;

typedef struct {
	Node *root;
	Node *leaf;
	int num;
	SeqPos from;
	SeqPos to;
	char mark;
} Domain;

typedef struct {
	NodeSet *nodes;
	EdgeSet *edges;
	NameHash *nhash;
} SimGraph;

typedef struct {
	Region break1, break2, boundary1, boundary2;
	Region seq, aliM;
	Region consreg;
	SeqPos meanlen;
	SeqPos minlen;
/*
	SeqPos maxlen;
*/
} NewSeqInfo;
typedef struct {
	Node *root;
	pList *clusters;
	pList *members;
	specFlag spflag;
} ClusterInfo;

typedef struct {
	int cutoff, cutoff2, sumcut;
	int distscale;
	int logdistscale;
	Dist missdist, missscore;
	float missratio, cutoff2ratio;
	int DEBUG, DEBUG_val;
	char *DEBUG_ent;
	int VERBOSE;
	int verbose_step1, verbose_step2;
	int minsp, minent;	/* minimum num. spececies/genes */
	int delete_small;	/* delete small grp. before dom. split */
	float ovlpratio, ovlpratio2;
	float coverage, covfilt;
	int mincutcnt;
	int minovlp, minovlp2;
	int minlen, minlen2;
	int min_alilen;
	enum outStyle outstyle;
	int sim;
	int nobreak;
	double phylocutratio, phylocutratio2;
	double distdiffcut;
	int neighbor;
	int dpWin;
/*
	double chkInclClst, chkOvlpClst;
*/
	double adjInclRatio, adjOvlpRatio;
	char chkSpFlg;
	BoundaryOutput domBoundary;
	int domcut;
	float chkConnect;
	Dist iniGap, extGap, skipGap;
	double nbrScoreRatio,nbrScoreDecay,nbrScoreLim;
	double nbrConnRatio;
	int nbrConnGap;
	int phyloMode;
	int skipErrEnt;
	int revMatch;
	char *dumpfile;
	int treecheck;
	int noextend; /*tmp*/

	char *outgroupStr;
	char *spmaskStr;
	char outgroupRev;	/* specifying ingroup */
	char outgroupMode;

	char *genes;
	char NEW;
} Options;

Options Opt;

NodeSet *createNodes();
Node *getNode();
Node *dupNode(), *dupNodeID();
/*
Node *addNode(NodeSet *nodes, char *name, int cnt, int len,
	int minlen,int maxlen, Edge *child, char flag, specFlag spflag);
*/
Node *addNode(NodeSet *nodes, char *name, int cnt, int len,
	Region *consreg, int totlen,
	Edge *child, NodeFlag flag, specFlag spflag);
EdgeSet *createEdges();
Edge *getBestEdge();
Edge *getEdgeByNode(), *getEdgeByNode0();

Edge *addEdge(EdgeSet *, Node*, Node*, Region *, Region *, Dist,
		ConnCount, signed char);
Edge *addEdgeWithScore(EdgeSet *, Node*, Node*, Region *, Region *, Dist, Dist,
		ConnCount, signed char);
Node *getRootNext();

int calNewDist(int nn, Count cnt1, Count cnt2, Dist dist1, Dist dist2,
                Dist score1, Dist score2, ConnCount conn1, ConnCount conn2,
                Dist *newdist, Dist *newscore, ConnCount *newconn);

int restoreGraph(char *, SimGraph *, int curr_argc, char **cur_argv);
pairList *getPairByNode(), *getPairByNode0();
Cluster *slink();
Region *transformReg();
char *addName();
char *getName();
void printNode();
void printEdge();
void printEdge2();
int deleted();
int isNotRoot();
int nodeMarked();
int edgeMarked();
int cmpr_nodeid();
Dist getEdgeScoreSum();
double DistRatio();
#endif
