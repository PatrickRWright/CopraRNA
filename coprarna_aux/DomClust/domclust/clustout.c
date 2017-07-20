/*
 * DomClust: Hierarchical Clustering for Orthologous Domain Classification
 * Copyright (c) 2000-2006, Ikuo Uchiyama
 * All rights reserved.
 */
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "domclust.h"
#include "util.h"
#include "neighbor.h"
#include "spec.h"

#define MAXLINE 4000
#define WITH_POSITION 1

static char head[MAXLINE];
static Node *rootnode;

int outputHier0(Node *node, int lev, signed char dir);
int clearNodeOutMarks(pList *nodelist, NodeFlag flag);
int outputNewick0(Node *node, Dist prevdist);
int printDomInfo0(Domain *dom, char flag);
static pList **clustTab;


outputClusterInfo(ClusterInfo *cInfo, int clustnum, NodeSet *nodes)
{
	Domain *dom;
	int i, clnum = 1, clsize, spcnt;
	for (i = 0; i < clustnum; i++) {
		clsize = numelemList(cInfo[i].members);
		if (clsize == 0) {
			continue;
		}
		spcnt = spFlagCntW(cInfo[i].spflag);
		if (spcnt < Opt.minsp || clsize < Opt.minent) {
			continue;
		}
		if (Opt.outstyle != TABOUT) {
			printf("Cluster %d\n", clnum);
		}
		switch (Opt.outstyle) {
		case TREEOUT:
		case GRAPHOUT:
		case NEWICK:
		case HIER_DETAILOUT:
			outputClusterInfoMulti(&cInfo[i], clustnum, nodes);
			break;
		default:
			while (dom = (Domain *) shiftList(cInfo[i].members)) {
				if (Opt.outstyle == TABOUT) {
					printf("%d ", clnum);
				}
				printDomInfo(dom);
				putchar('\n');
			}
		}
		if (Opt.outstyle != TABOUT) {
			putchar('\n');
		}
		clnum++;
	}
}
outputClusterInfoMulti(ClusterInfo *cInfo, int clustnum, NodeSet *nodes)
{
	listIter iter;
	Node *node;
	int tn = 0;
	pList *nodelist;
	int testflag = 1;
	int clnum;

	clnum = numelemList(cInfo->clusters);

	setListIter(&iter, cInfo->clusters, 1);
	if (testflag) {
		nodelist = create_pList();
		while (node = (Node *) getListIter(&iter)) {
			testTreeOverlap(node, nodelist);
		}
	}
	setListIter(&iter, cInfo->clusters, 1);
	while (node = (Node *) getListIter(&iter)) {
		outputCluster_sub(nodes, node, cInfo->root, 0);
		if (++tn < clnum){
			if (Opt.outstyle != NEWICK &&
					Opt.outstyle != NEWICK_DIST) {
				printf("--");
			}
			printf("\n");
		}
	}
	if (testflag) {
		clearNodeOutMarks(nodelist, NODE_TMPMARK);
		clearNodeOutMarks(nodelist, NODE_DUPMARK);
		free_pList(nodelist);
	}
}

outputClusters(pList *roots, NodeSet *nodes)
{
	int i;
	NodeSetIter nodeiter;
	Node *node;
	int clustnum = 0;
	int spcnt;

	if (Opt.outstyle == DOMAINOUT) {
		return;
	}
	if (Opt.outstyle == CLUSTTAB) {
		outputSpNamesForClustTab();
		alloc_clustTabList();
	}

	getRootInit(&nodeiter, nodes, -1, 1);
	while (node = getRootNext(&nodeiter)) {
		if (Opt.minent > 1 && ! isIntNode(node)) {
			/* leaf == singleton */
			return 0;
		}
		spcnt = spFlagCntW(node->spflag);
		if (spcnt < Opt.minsp || node->cnt<Opt.minent) {
			continue;
		}
		if (Opt.outstyle != NEIGHBOR) {
			++clustnum;
			if (Opt.outstyle != TABOUT && Opt.outstyle != CLUSTTAB){
				printf("Cluster %d", clustnum);
				putchar('\n');
			}
		}
		outputCluster_sub(nodes, node, node, clustnum);
		if (Opt.outstyle != TABOUT) {
			putchar('\n');
		}
	}
	if (Opt.outstyle == CLUSTTAB) {
		free_clustTabList();
	}
}
clearNodeOutMarks(pList *nodelist, NodeFlag flag)
{
	listIter iter;
	Node *n;
	setListIter(&iter, nodelist, 1);
	while (n = (Node *) getListIter(&iter)) {
		unsetFlagNode(n, flag);
	}
}
testTreeOverlap(Node *node, pList *nodelist)
{
	Node *child1 = NULL, *child2 = NULL;
	
	if (testFlagNode(node, NODE_TMPMARK)) {
		/* an overlapping node */
		setFlagNode(node, NODE_DUPMARK);
	} else {
		setFlagNode(node, NODE_TMPMARK);
		pushList(nodelist, node);
	}
	getChilds(node, &child1, &child2);
	if (child1) {
		testTreeOverlap(child1, nodelist);
	}
	if (child2) {
		testTreeOverlap(child2, nodelist);
	}
}
outputCluster_sub(NodeSet *nodes, Node *node, Node *root, int clustnum)
{
	rootnode = root;
	switch (Opt.outstyle) {
	case TREEOUT:
		head[0] = '\0';
		outputTree(node,0,1);
		break;
	case GRAPHOUT:
		outputGraph(node);
		break;
#ifdef NEIGHBOROUT
	case NEIGHBOR:
		outputNeighborNode(nodes, node);
		return 0;
		break;
#endif
	case NEWICK:
	case NEWICK_DIST:
		outputNewick(node);
		break;
	case TABOUT:
		outputSimpleTab(node, clustnum);
		break;
	case NORMALOUT:
		outputNormal(node);
		break;
	case CLUSTTAB:
		outputClustTab(node, clustnum);
		break;
	case HIER_DETAILOUT:
/*
	case NORMALOLD:
*/
	default:
		outputHier(node);
		break;
	}
}

outputNormal(Node *node)
{
	pList *members = create_pList();
	Domain *dom;
	collectCluster(node, members);
	while (dom = (Domain *) shiftList(members)) {
		printDomInfo(dom);
		printf("\n");
	}
	freeList(members);
}
alloc_clustTabList()
{
	int i;
	if ( (clustTab = (pList **) malloc(sizeof(pList*) * SPnum))==NULL ){
		allocError("clusters");
	}
	for (i = 0; i < SPnum; i++) {
		if ( (clustTab[i] = create_pList())==NULL ) {
			allocError("clustTab");
		}
	}
}
free_clustTabList()
{
	int i;
	for (i = 0; i < SPnum; i++) {
		free(clustTab[i]);
	}
	free(clustTab);
}
outputSpNamesForClustTab(Node *node, int clustnum)
{
	int i;
	char *spname;
	printf("#id");
	for (i = 0; i < SPnum; i++) {
		spname = getSPname(i);
		printf("\t%s",spname);
	}
	printf("\n");
}
outputClustTab(Node *node, int clustnum)
{
	pList *members = create_pList();
	Domain *dom;
	Node *lnode;
	char spname[SPNAMELEN];
	char *p, *q;
	int i, j, spid;

	collectCluster(node, members);
	while (dom = (Domain *) shiftList(members)) {
		lnode = dom->leaf;
		for (p = lnode->name, q = spname; *p && *p != ':'; p++, q++) {
			*q = *p;
		}
		*q = '\0';
		spid = getSPid(spname);
		pushList(clustTab[spid], dom);
	}
	printf("%d",clustnum);
	for (i = 0; i < SPnum; i++) {
		j = 0;
		printf("\t");
		while (dom = (Domain *) shiftList(clustTab[i])) {
			if (j++ > 0) printf(" ");
			printDomName(dom);
		}
	}
	freeList(members);
}
outputSimpleTab(Node *node, int clustnum)
{
	pList *members = create_pList();
	Domain *dom;
	collectCluster(node, members);
	while (dom = (Domain *) shiftList(members)) {
		printf("%d ", clustnum);
		printDomInfo(dom);
		printf("\n");
	}
	freeList(members);
}
printNodeInfo(Node *node, Node *root)
{
	Domain *dom;
	int domn;
	printf("%s", node->name);
	if (node->domains) {
		domn = getDomainNoMark(node->domains, root, &dom);
		if (domn) {
			if (numelemList(node->domains) > 1) {
				printf("(%d)", domn);
			}
			printf(" %d %d",(int)dom->from, (int)dom->to);
		} else {
			printf(" %d %d",1,(int)node->len);
		}
	} else {
		printf(" %d %d",1,(int)node->len);
	}
}
printDomName(Domain *dom)
{
	printDomInfo0(dom, 0);
}
printDomInfo(Domain *dom)
{
	printDomInfo0(dom, WITH_POSITION);
}
printDomInfo0(Domain *dom, char flag)
{
	Node *node = dom->leaf, *root = dom->root;
	int domn;
	printf("%s", node->name);
	if (node->domains) {
		if (dom->num) {
			if (Opt.outstyle == TABOUT) {
				printf(" %d", dom->num);
			} else if (numelemList(node->domains) > 1) {
				printf("(%d)", dom->num);
			}
			if (flag & WITH_POSITION) {
			    printf(" %d %d",(int)dom->from, (int)dom->to);
			}
		} else {
			domn = getDomainNoMark(node->domains, root, &dom);
			if (domn) {
				if (numelemList(node->domains) > 1) {
					printf("(%d)", domn);
				}
				if (flag & WITH_POSITION) {
				    printf(" %d %d",(int)dom->from, (int)dom->to);
				}
			} else if (flag & WITH_POSITION) {
				printf(" %d %d",1,(int)node->len);
			}
		}
	} else if (flag & WITH_POSITION) {
		printf(" %d %d",1,(int)node->len);
	}
}

/* collect cluster members by traversing the tree */
collectCluster(Node *root, pList *nodelist)
{
	collectCluster_sub(root, root, nodelist);
	clearDomainMark(nodelist);
}
collectCluster_sub(Node *root, Node *node, pList *nodelist)
{
	Edge *child = node->child;
	Domain *dom;
	int domn;
	if (isFlankNode1(node)) {
		collectCluster_sub(root, child->node1, nodelist);
	} else if (isFlankNode2(node)) {
		collectCluster_sub(root, child->node2, nodelist);
	} else if (isIntNode(node)) {
		collectCluster_sub(root, child->node1, nodelist);
		collectCluster_sub(root, child->node2, nodelist);
	} else {
		/* leaf node */
		if (node->domains) {
			domn = getDomain(node->domains, root, &dom);
			if (domn) {
				pushList(nodelist, dom);
			}
		}
	}
}

outputHier(Node *node)
{
	rootnode = node;
	outputHier0(node, 1, (signed char) 1);
}
outputHier0(Node *node, int lev, signed char dir)
{
	int i;
	Edge *child;
	char indent[200];
	Node *maxnode;
	int maxcnt;
	Domain *dom;

	indent[0] = '\0';
	if (Opt.outstyle == HIER_DETAILOUT) {
		for (i = 0; i < lev; i++) {
			strcat(indent, "  ");
		}
	}
	if (isFlankNode1(node)) {
		if (Opt.outstyle == HIER_DETAILOUT) {
			printf("%s - (%s,%d) %d/%d (flanking node)\n",
				indent, node->name,node->id,
				nodeLen(node),node->len);
		}
		outputHier0(node->child->node1, lev, dir);
	} else if (isFlankNode2(node)) {
		if (Opt.outstyle == HIER_DETAILOUT) {
			printf("%s - (%s,%d) %d/%d (flanking node)\n",
				indent, node->name,node->id,
				nodeLen(node),node->len);
		}
		outputHier0(node->child->node2, lev, dir*node->child->dir);
	} else if (isIntNode(node)) {
		child = node->child;
		if (Opt.outstyle == HIER_DETAILOUT) {
			printf("%s", indent);
			printf("[%d](%s,%d) [%d:%d] %d/%d %.2f %d  ",
				lev,node->name,node->id,
				node->newreg.from,node->newreg.to,
				nodeLen(node), node->len,
				MEASURE(child),child->connect);
			print_specFlag(node->spflag);

			if (node->left) {
				printf("%s", indent);
				printf("    Left: ");
/*
				printList(node->left, printNode);
				printNeighborList(node->left, node);
*/
				maxcnt = getMaxNeighbor(node->left,node,1,
					&maxnode,NULL);
				if (maxcnt && maxcnt >=
					(double) node->cnt * Opt.nbrConnRatio) {
					printNode(maxnode);
					printf(" [%d]",maxcnt);
				}
				putchar('\n');
			}
			if (node->right) {
				printf("%s", indent);
				printf("    Right: ");
/*
				printList(node->right, printNode);
				printNeighborList(node->right, node);
*/
				maxcnt = getMaxNeighbor(node->right,node,1,
					&maxnode,NULL);
				if (maxcnt && maxcnt >=
					(double) node->cnt * Opt.nbrConnRatio) {
					printNode(maxnode);
					printf(" [%d]",maxcnt);
				}
				putchar('\n');
			}
		}
		outputHier0(child->node1, lev+1, dir);
		outputHier0(child->node2, lev+1, dir*node->child->dir);
	} else {
		printf("%s", indent);
		printNodeInfo(node,rootnode);
		if (Opt.revMatch) {
			printf(" %d", dir);
		}
		putchar('\n');
	}
}
outputGraph(Node *node)
{
	Node *maxnode = NULL;
	enum {Leaf,NonLeaf} type = NonLeaf;
	int maxcnt;

	if (isFlankNode1(node)) {
		if (isIntNode(node->child->node1)) {
			printf("%d %d\n", node->id, node->child->node1->id);
		} else {
			printf("%d %s\n", node->id, node->child->node1->name);
		}
		outputGraph(node->child->node1);
	} else if (isFlankNode2(node)) {
		if (isIntNode(node->child->node2)) {
			printf("%d %d\n", node->id, node->child->node2->id);
		} else {
			printf("%d %s\n", node->id, node->child->node2->name);
		}
		outputGraph(node->child->node2);
	} else if (isIntNode(node)) {
		if (isIntNode(node->child->node1)) {
			printf("%d %d\n", node->id, node->child->node1->id);
		} else {
			printf("%d %s\n", node->id, node->child->node1->name);
		}
		if (isIntNode(node->child->node2)) {
			printf("%d %d\n", node->id, node->child->node2->id);
		} else {
			printf("%d %s\n", node->id, node->child->node2->name);
		}
		outputGraph(node->child->node1);
		outputGraph(node->child->node2);
	} else {
		type = Leaf;
	}

	if (node->left) {
		maxcnt = getMaxNeighbor(node->left,node,0,&maxnode,NULL);
		if (maxcnt && maxcnt >=
				rint(node->cnt * Opt.nbrConnRatio)) {
			if (type == NonLeaf) {
				printf("%d %d L\n", node->id, maxnode->id);
			} else {
				if (isLeaf(maxnode) && ! isRoot(maxnode)) {
					printf("%s %s L\n",node->name,maxnode->name);
				}
			}
		}
	}
	if (node->right) {
		maxcnt = getMaxNeighbor(node->right,node,0,&maxnode,NULL);
		if (maxcnt && maxcnt >=
				rint(node->cnt * Opt.nbrConnRatio)) {
			if (type == NonLeaf) {
				printf("%d %d R\n", node->id, maxnode->id);
			} else {
				if (isLeaf(maxnode) && ! isRoot(maxnode)) {
					printf("%s %s R\n",node->name,maxnode->name);
				}
			}
		}
	}
}

outputTree(Node *node, int lev, int dir)
{

	Edge *child;
	Domain *dom = NULL;
	char markchar = ( (node->flag & NODE_DUPMARK) ? '*' : '+');

	if (! isIntNode(node)) {
		printf("%s%c%s", head, markchar, "- ");
		if (node->domains) {
			if (numelemList(node->domains) > 1) {
				getDomainNoMark(node->domains, rootnode, &dom);
			}
			if (dom) {
				printDomInfo(dom);
			} else {
				printf("%s %d %d", node->name, 1, (int) node->len);
			}
		} else {
			printf("%s", node->name);
		}
		putchar('\n');
	} else {
		if (isFlankNode1(node)) {
			outputTree(node->child->node1, lev, dir);
		} else if (isFlankNode2(node)) {
			outputTree(node->child->node2, lev, dir);
		} else {
			child = node->child;
			if (child) {
				head[lev*2] = '\0';
				strcat(head, (lev > 0 && dir == -1 ? "| " : "  "));
				outputTree(child->node1, lev+1, 1);
	
				head[lev*2] = '\0';
				printf("%s%c%s%.1f\n", head, markchar,
						"-| ", MEASURE(child));
	
				strcat(head, (lev > 0 && dir == 1 ? "| " : "  "));
				outputTree(child->node2, lev+1, -1);
			} else {
				printf("%s%c%s[%s]", head, markchar,
					"- ", node->name);
				printf(" %d %d",1,(int)node->len);
				putchar('\n');
			}
		}
	}
}

static char namebuf[NAMELEN];
outputNewick(Node *node)
{
	outputNewick0(node, ((Dist) -1.0) );
	printf(";");
}
outputNewick0(Node *node, Dist prevdist)
{
	Dist dist = 0;
	Dist currdist = 0;

	if (isFlankNode1(node)) {
		outputNewick0(node->child->node1, prevdist);
	} else if (isFlankNode2(node)) {
		outputNewick0(node->child->node2, prevdist);
	} else {
		if (isIntNode(node)) {
			currdist = node->child->dist;
			printf("(");
			outputNewick0(node->child->node1, currdist);
			printf(", ");
			outputNewick0(node->child->node2, currdist);
			printf(")");
			if (Opt.outstyle == NEWICK_DIST) {
				if (! Opt.sim && prevdist >= 0) {
					dist = prevdist - currdist;
					printf(":%.1f ",  (float)dist);
				}
			}
		} else {
			Domain *dom = NULL;
			int domn;
			if (numelemList(node->domains) > 1) {
				domn = getDomainNoMark(
					node->domains, rootnode, &dom);
				sprintf(namebuf, "%s#%d", node->name, domn);
			} else {
				strcpy(namebuf, node->name);
			}
			convname(namebuf);
			if (Opt.outstyle == NEWICK_DIST) {
				printf("%s", namebuf);
				if (! Opt.sim && prevdist >= 0) {
					dist = prevdist;
					printf(":%.1f ",  (float)dist);
				}
			} else {
				printf("%s", namebuf);
			}
		}
	}
}
convname(char *name)
{
	register char *p;
	for (p = name; *p; p++) {
		if (*p == ':') {
			*p = '_';
		}
	}
}
