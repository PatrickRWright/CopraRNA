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

#define ONE .99999999

cmpr_node_by_size(Node **a, Node **b)
{
	return (*b)->cnt - (*a)->cnt;
}
pList *createRootNodeList(NodeSet *nodes);

outputResults(SimGraph *SimG, int argc, char **argv)
{
	NodeSet *nodes = SimG->nodes;
	EdgeSet *edges = SimG->edges;
	pList *clustRoots;
	ClusterInfo *cInfo;
	int i, clustnum;

	setTotalNodeNum(nodes);

	if (Opt.outstyle == DUMP) {
		dumpGraph(SimG, argc, argv);
		exit(0);
	}

	if (Opt.spmaskStr) {
		checkSkipMaskAll(nodes);
	}
	/** phylogenetic tree cutting procedure **/
        phyloCutAll(nodes, edges);

	outgroupCheckAll(nodes);


	/** reexamine neighboring clusters **/
	if (Opt.adjInclRatio >= ONE) {
		checkIncludedClustersAll(nodes);
	}

	if (Opt.delete_small) {
		/** delete small clusters before domain definition **/
		checkClusterSize(nodes);
	}
	/** domain definition **/
	checkDomains(nodes);

	clustRoots = createRootNodeList(nodes);
	sortList(clustRoots, cmpr_node_by_size);

	if (Opt.adjOvlpRatio ||
			(Opt.mincutcnt > 1) || 
			(Opt.adjInclRatio && (Opt.adjInclRatio < ONE)) ) {
		clustnum = checkOverlapClustersAll(clustRoots,
			&cInfo, nodes);
	        outputClusterInfo(cInfo, clustnum, nodes);
	} else {
	        outputClusters(clustRoots, nodes);
	}
}

checkConnect(Node *node) {
	float cutcnt;
	if (Opt.chkConnect > 1) {
		/** count **/
		cutcnt = Opt.chkConnect;
	} else {
		/** ratio **/
		cutcnt = ((float)node->child->node1->cnt +
			node->child->node2->cnt) * Opt.chkConnect;
	}
	if (node->child->connect < cutcnt) {
		/* cut */
		return 0;
	} else {
		return 1;
	}
}

#define DEBUG_ON(x) if(x){ tmpflag=1; } else { tmpflag=0; }
int tmpflag;

phyloCutAll(NodeSet *nodes, EdgeSet *edges)
{
	int i;
	Node *node;
	NodeSetIter iter;
	int cut;

	/* decreasing order:
		note that new root nodes may be created after the cut */
	getRootInit(&iter, nodes, -1, 0);
	while (node = getRootNext(&iter)) {
		/** cut non-recursively
		   further cut should be done in the succeeding loops **/
		cut = phyloCut(node);
	}
}

phyloCut_disttest(Node *node, Node *parent)
{

	/* test mode */
	if (node->child && parent->child) {
		double distr = DISTRATIO( MEASURE(node->child),
				MEASURE(parent->child) );
		if (distr <= Opt.distdiffcut){
			return 1;
		}
	}
	return 0;
}
postPhyloCut(Node *node, Node *parent)
{
	int stat;

if (nodeDeleted0(node)){
	printf("phylocut: ????\n");
	if (node->child) {
		printf("%d\n",isRoot(node));
		printf("%d\n",isRoot(node->child->node1));
		printf("%d\n",isRoot(node->child->node2));
	}
}
	if (! parent) {
		return 0;
	}
	/** When there is a parent, we have entered this node after
		cutting this parent node. So we must cancel
		the domain boundaries created during the cluster
		merging process.
	**/

	restoreNeighborList(parent, node);
	restoreBreak(node, parent);

	if ( (stat=checkUnvisitedParent(node)) == 0) {
		/** do not proceed further
			if there remains an unvisited parent **/
		return 0;
	}


	/* ### There is no unvisited parent */

	/* delete the flanking root node if exist,
		which is probably retained as a long domain */
	if (node->parentL && isRoot(node->parentL)) {
		deleteNode(node->parentL);
		restoreBreak(node, node->parentL);
	}
	if (node->parentR && isRoot(node->parentR)) {
		deleteNode(node->parentR);
		restoreBreak(node, node->parentR);
	}

	/* if the parent is not a root then make this node as a root */
	if (! node->parent || isNotRoot(node->parent)) {
		makeRootNode(node);
	}

	return 0;
}
phyloCut(Node *node)
{
	int cut = 0;
	if (node->child) {
		if (isFlankNode(node)) {
			Node *nbrNode = NULL, *child = NULL;
			/* nbrNode = central node (newnode0) */
			if (isFlankNode1(node)) {
				child = node->child->node1;
			} else if (isFlankNode2(node)) {
				child = node->child->node2;
			} else {
				fprintf(stderr, "????\n");
			}
			/* find the neighboring (central) node */
			if (node == child->parentL) {
				/* right neighbor */
				nbrNode = (Node *) getNeighborNode(
							node, 1, NULL);
			} else if (node == child->parentR) {
				/* left neighbor */
				nbrNode = (Node *) getNeighborNode(
							node, -1, NULL);
			} else {
				/** ???? never come here **/
				/** node == child->parent || child->parentM */
				fprintf(stderr, "ERROR?? a central node is treated as a flanking node\n");
				assert(0);
				makeRootNode(child);
				deleteRootNode(node);
				return 0;
			}

			if (nbrNode && nbrNode->flag & NODE_DELETED) {
				/* the central node is already cut */
				/* never come here ? */
				fprintf(stderr, "????\n");
				deleteNode(node);
				return postPhyloCut(child, node);
			} else {
				if (checkNodeLen(node) < 2) {
					/* short segment */
					deleteNode(node);
					checkIncludedCluster(node);
				} else {
					/* having sufficient length */
					/* make this node as a root */
					makeRootNode(node);
				}
				return 0;
			}
			/* never come here */
		}

	/** Test for the phylogenetic tree cutting procedure **/
		cut = duplicationCheckAll(node);
		if (cut == 1) {
			/** immediately cut this node **/
		} else if ( cut == 2 ) {
			if ((Opt.chkConnect && ! checkConnect(node)) ||
				(Opt.sumcut &&
			  	Opt.sumcut > getEdgeScoreSum(node->child))) {

				/* weak cut condition (connect/sumcut) */
			} else if (Opt.distdiffcut) {
				/* weak cut condition (distdiffcut) */
				int cut1, cut2; 
				cut1 = phyloCut_disttest(node->child->node1, node);
				cut2 = phyloCut_disttest(node->child->node2, node);
				if (cut1 || cut2) {
					/* cut */
				} else {
					cut = 0;
				}
			} else {
				cut = 0;
			}
		}
		if (cut) {
			cutNode(node);
		}
	}
	return cut;
}
cutNode(Node *node)
{
	if (Opt.DEBUG & DBG_nbrrestore) {
		printf("Cut: "); printNode(node); putchar('\n');
		print_specFlag(node->child->node1->spflag);
		print_specFlag(node->child->node2->spflag);
	}

	deleteRootNode(node);
	postPhyloCut(node->child->node1, node);
	if (node->child->node1 != node->child->node2) {
		/* NOT self_match */
		postPhyloCut(node->child->node2, node);
	}
}

/* skip nodes that contain only masked species */
checkSkipMaskAll(NodeSet *nodes)
{
	Node *node;
	NodeSetIter iter;

	getRootInit(&iter, nodes, -1, 0);
	while (node = getRootNext(&iter)) {
		if (! nodeDeleted0(node)) {
			checkSkipMask(node);
		}
	}
}
checkSkipMask(Node *node)
{
	Node *child1, *child2;
	int maskchk = spFlagMaskCheck(node->spflag, SPflags.spMask);
	if (maskchk == 2) {
		/** all organisms are matched with specified ones **/
		return 0;
	} else if (maskchk == 0) {
		/** no organism is matched with specified ones **/
		deleteRootNode(node);
		return 0;
	}
	if (node->child) {
		child1 = node->child->node1;
		child2 = node->child->node2;
		if (! isFlankNode(node)) {
			if (spFlagMaskCheck(child1->spflag, SPflags.spMask) == 0) {
				/* child1 should be deleted */
				unsetFlagNode(node, NODE_INT1);
			} else if (spFlagMaskCheck(child2->spflag, SPflags.spMask) == 0) {
				/* child2 should be deleted */
				unsetFlagNode(node, NODE_INT2);
			}
		}
	}
	return 0;
}

/* delete node and recover neighbor relationships */
deleteNode(Node *node)
{
	listIter *iterL, *iterR;
	Neighbor *nbrL, *nbrR;
	Node *nodeL, *nodeR;
	int cnt;
	char status;

	if (nodeDeleted0(node)) {
		return 0;
	}

if (Opt.DEBUG & DBG_nbrrestore) {
	printf("DelFlank: ");
	printNode(node);
	printf("\n");
}

	/***
		 +<-------newlink------->+
		 |     nbrL       nbrR   |
		nodeL --+-- node --+-- nodeR
			   delete
	**/

	deleteRootNode(node);
	if (node->left && node->right) {
		iterL = createListIter(node->left, 1);
		iterR = createListIter(node->right, 1);
		while (nbrL = (Neighbor *)getListIter(iterL)) {
			nodeL = (Node *) (((Node *) nbrL->node1 == node)
					? nbrL->node2 : nbrL->node1);

			setListIter(iterR, node->right, 1);

			while (nbrR = (Neighbor *)getListIter(iterR)) {
				nodeR = (Node *) (
					    ((Node *) nbrR->node1 == node) ?
					    nbrR->node2 : nbrR->node1);


				cnt = (nbrL->count < nbrR->count) ?
					nbrL->count : nbrR->count;
			/***
				   nbrL       nbrR
				L -{1}- curr -{-1}- R
				1         1        -1
			**/
				if (nbrL->status == 0 && nbrR->status == 0){
					status = 0;
				} else if (nbrL->status == 2 || nbrR->status == 2){
					status = 2;
				} else {
					status = 1;
				}

if (Opt.DEBUG & DBG_nbrrestore) {
	printf("newnode: ");
	printNode(nodeL); printNode(nodeR);
	printf(" %d\n", status);
}
				addNeighborRestore(nodeL,nodeR,
					nbrL->dir,nbrR->dir,cnt, status);
			}
		}
	}
	return 1;
}
pList *createRootNodeList(NodeSet *nodes)
{
	NodeSetIter iter;
	Node *node;
	pList *rootNodes = create_pList();

	getRootInit(&iter, nodes, -1, 1);
	while (node = getRootNext(&iter)) {
		pushList(rootNodes, node);
	}
	return rootNodes;
}
outgroupCheckAll(NodeSet *nodes)
{
	NodeSet iter;
	Node *node;

	getRootInit(&iter, nodes, -1, 0);
	while (node = getRootNext(&iter)) {
		if (outgroupCheck(node->spflag) == 2) {
		/* the set contains only the outgroup organisms */
			deleteRootNode(node);
		}
	}
}
outgroupCheck(specFlag spflag)
{
	return spFlagMaskCheck(spflag, SPflags.outGroup);
}
spFlagMaskCheck(specFlag spflag, specFlag spFlagPat)
{
	specFlag matched;
	int matched_spcnt, spcnt;

	spFlagAND(spflag, spFlagPat, matched);
	matched_spcnt = spFlagCnt(matched);
	if (matched_spcnt) {
		spcnt = spFlagCnt(spflag);
		if (matched_spcnt == spcnt) {
		/* all organisms in the set are matched */
			return 2;
		} else {
		/* there exists at least one matched organism */
			return 1;
		}
	}
	return 0;
}
duplicationCheckAll(Node *node)
{
	if (Opt.outgroupMode == 1) {
		int ret = outgroupCheck(node->spflag);
	/* cut if there is an outgroup organism in the set */
		if (ret == 1) return ret;
	}
	return duplicationCheck(node);
}
/* check the phylogenetic tree cutting criterion */
duplicationCheck(Node *node)
{
	return duplicationCheck0(node->spflag, node->child->node1->spflag,
		node->child->node2->spflag);
}
duplicationCheck0(specFlag spflag0, specFlag spflag1, specFlag spflag2)
{
	int mchcnt, spcnt1, spcnt2, mchcnt2;
	double wt_mchcnt, wt_spcnt1, wt_spcnt2;

	mchcnt = spFlagANDcnt(spflag1, spflag2);
	spcnt1 = spFlagCnt(spflag1);
	spcnt2 = spFlagCnt(spflag2);

	if (mchcnt * spcnt1 * spcnt2 == 1) {
		return 0;
	}

	wt_mchcnt = spFlagANDcntW(spflag1, spflag2);
	wt_spcnt1 = spFlagCntW(spflag1);
	wt_spcnt2 = spFlagCntW(spflag2);

	/** cut node when two subgroups share organisms with a ratio of
		more than Opt.phylocutratio: do not cut when the
		ratio is just Opt.phylocutratio */
	if ((wt_mchcnt > (double) wt_spcnt1 * Opt.phylocutratio
		 || wt_mchcnt > (double) wt_spcnt2 * Opt.phylocutratio)
			 && (mchcnt * spcnt1 * spcnt2 != 1) ) {
		/** Immediately cut this node **/
		return 1;
	}
	if (Opt.treecheck) {
		/* Mapping the nodes of gene tree onto the species tree */
		int spnode0, spnode1, spnode2;
		spnode0 = sptree_MatchFlags(spflag0);
		spnode1 = sptree_MatchFlags(spflag1);
		spnode2 = sptree_MatchFlags(spflag2);
		if (spnode0 == spnode1 || spnode0 == spnode2) {
			/* duplication !! */
			/* sum of the weights of the duplicated species */
			wt_mchcnt =
			    sptree_MatchFlagsCntW(spflag1,spflag2,spnode0);

			if ((wt_mchcnt > (double) wt_spcnt1 * Opt.phylocutratio
			 || wt_mchcnt > (double) wt_spcnt2 * Opt.phylocutratio)
		 		&& (mchcnt * spcnt1 * spcnt2 != 1) ) {
				/** Cut this node **/
				return 1;
			}
		}
	}
	if ((wt_mchcnt > (double) wt_spcnt1 * Opt.phylocutratio2
		 || wt_mchcnt > (double) wt_spcnt2 * Opt.phylocutratio2)
	) {
		return 2;
	}
	/** Do not cut this node **/
	return 0;
}

restoreBreak(Node *node, Node *parent)
{
	if (! parent) return;

	if (parent == node->parentL) {
		node->brk.from = INFPOS;
		node->parentL = NULL;
	} else if (parent == node->parentR) {
		if (node->parentM) {
			node->brk2.to = SUPPOS;
		} else {
			node->brk.to = SUPPOS;
		}
		node->parentR = NULL;
	} else if (parent == node->parentM) {
		node->brk.to = SUPPOS;
		node->brk2.from = INFPOS;
		node->parentM = NULL;
	} else if (parent == node->parent) {
		if(definedPos(node->brk.from)) {
			node->brk.from = INFPOS;
		} else if (definedPos(node->brk.to)){
			node->brk.to = SUPPOS;
		}
		if (parent->child->node1 == parent->child->node2) {
			if(definedPos(node->brk2.from)) {
				node->brk2.from = INFPOS;
			} else if(definedPos(node->brk2.to)) {
				node->brk2.to = SUPPOS;
			}
		}
		node->parent = NULL;
	} else {
		/*** error **/
		if (node->parent != NULL &&
	node->parentR != NULL && node->parentL != NULL && node->parentM!=NULL) {
		fprintf(stderr, ">>ERROR %d: %d,%d,%d\n",
			parent->id,node->parent,node->parentL,node->parentR);
			printf("??? %s,%d %d,%d,%d,%d,%d\n",node->name,node->id,node->parent,node->parentL,node->parentR,parent,parent->id);
			printf("??? %d %d,%d\n",parent->child->node1, parent->child->node2);
			printNode(node); putchar(' ');
			if (node->parent){
				printNode(node->parent); putchar('\n');
			}
			if (node->parentL){
				putchar('L');
				printNode(node->parentL); putchar('\n');
			}
			if (node->parentR){
				putchar('R');
				printNode(node->parentR); putchar('\n');
			}
		}
	}
}
checkUnvisitedParent(Node *node)
{
	if ( (node->parent == NULL || nodeVisited(node->parent)) &&
		(node->parentL == NULL || nodeVisited(node->parentL)) &&
		(node->parentR == NULL || nodeVisited(node->parentR)) ) {

		/** there remains no unvisited parent **/
		return 1;
	}
	/** there remains unvisited parent **/
	return 0;
}

checkClusterSize(NodeSet *nodes)
{
	NodeSetIter iter;
	Node *node;
	double spcnt;
	int i;

	getRootInit(&iter, nodes, -1, 0);
	while (node = getRootNext(&iter)) {
		if (node->child) {
			spcnt = spFlagCntW(node->spflag);
			if (spcnt < Opt.minsp || node->cnt < Opt.minent) {
				deleteNode(node);
			}
		}
	}
}
