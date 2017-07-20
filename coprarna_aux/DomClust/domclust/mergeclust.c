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
#include "hash.h"
#include "plist.h"

int concatNeighbor(Node *nodeS, Node *nodeL, char dir);

cmpr_node_by_id(Domain **dom1, Domain **dom2)
{
	return (*dom1)->leaf->id - (*dom2)->leaf->id;
}
cmpr_node_by_id_pos(Domain **dom1, Domain **dom2)
{
	if ((*dom1)->leaf->id == (*dom2)->leaf->id) {
		return (*dom1)->from - (*dom2)->from;
	} else {
		return (*dom1)->leaf->id - (*dom2)->leaf->id;
	}
}

cinfoHsearch(Hash *cinfoHash, int nodeid, ClusterInfo **ci)
{
	HENTRY hent;
	hent.key = (char *) nodeid;
	if (HIsearch(cinfoHash, &hent, FIND)){
		*ci = (ClusterInfo *) hent.datum;
		return 1;
	} else {
		return 0;
	}
}
createClusterInfo(pList *nodes, ClusterInfo **ret_cinfo, Hash **cinfoHash)
{
	ClusterInfo *cinfo, *ci;
	Node *node;
	listIter iter;
	int clustnum = numelemList(nodes);
	int i;
	HENTRY hent;

	if ( (cinfo = (ClusterInfo *)
		malloc(clustnum * sizeof(ClusterInfo)) ) == NULL) {
		fprintf(stderr, "Can't allocate memory\n");
		exit(1);
	}
	setListIter(&iter, nodes, 1);
	i = 0;
	while (node = (Node *) getListIter(&iter)) {
		if (isRoot(node)) {
			cinfo[i].root = node;
			cinfo[i].clusters = create_pList();
			pushList(cinfo[i].clusters, node);
			cinfo[i].members = create_pList();
			collectCluster(node, cinfo[i].members);
			copySPFlag(node->spflag, cinfo[i].spflag);
			i++;
		}
	}
	clustnum = i;
	if (cinfoHash) {
		*cinfoHash = Hcreate(clustnum * 50);
		for (i = 0; i < clustnum; i++) {
			hent.key = (char *) cinfo[i].root->id;
			hent.datum = (char *) &(cinfo[i]);
			HIsearch(*cinfoHash, &hent, ENTER);
		}
	}
	*ret_cinfo = cinfo;
	return clustnum;
}
cinfoMeanLen(ClusterInfo *cinfo)
{
	listIter iter;
	Domain *dom;
	int len = 0, cnt = 0;

	setListIter(&iter, cinfo->members, 1);
	while (dom = (Domain *) getListIter(&iter)) {
		len += (dom->to - dom->from + 1);
		cnt++;
	}
	if (len == 0) return 0;
	return (len / cnt);
}
printCinfo(ClusterInfo *cinfo)
{
	listIter iter;
	Domain *dom;
	setListIter(&iter, cinfo->members, 1);
	printf("Root:"); printNode(cinfo->root); printf("\n");
	while (dom = (Domain *) getListIter(&iter)) {
		printf("%s %d %d %d\n",dom->leaf->name,dom->num,dom->from,dom->to);
	}
}
renumCinfo(ClusterInfo *cinfo, int clustnum)
{
	int i;
	pList *alldom = create_pList();
	Domain *dom, *prevdom;
	int domn;
	for (i = 0; i < clustnum; i++) {
		pushListAll(alldom, cinfo[i].members);
	}
	sortList(alldom, cmpr_node_by_id_pos);
	prevdom = NULL;
	while (dom = (Domain *) shiftList(alldom)) {
		if (! prevdom || dom->leaf->id != prevdom->leaf->id) {
			domn = 1;
		} else {
			domn++;
		}
		dom->num= domn;
		prevdom = dom;
	}
	free_pList(alldom);
}
/** merge cluster1 into cluster2: cinfo1=>null, cinfo2=>merged_list **/
mergeClusters(ClusterInfo *cinfo1, ClusterInfo *cinfo2, NodeSet *nodes)
{
	listIter iter1, iter2;
	pList *merged_list;
	Domain *d1, *d2;
	int mergecnt = 0, cnt1 = 0, cnt2 = 0;
	int mergeflag = 0;
	Node *newnode;


	sortList(cinfo1->members, cmpr_node_by_id);
	sortList(cinfo2->members, cmpr_node_by_id);
	merged_list = create_pList();

	/* Check merge condition */

	setListIter(&iter1, cinfo1->members, 1);
	d1 = (Domain *) getListIter(&iter1);
	setListIter(&iter2, cinfo2->members, 1);
	d2 = (Domain *) getListIter(&iter2);

	while (d1 || d2) {
		if (! d2 || (d1 && cmpr_node_by_id(&d1,&d2) < 0)) {
			d1 = (Domain *) getListIter(&iter1);
			cnt1++;
		} else if (! d1 || (d2 && cmpr_node_by_id(&d2,&d1) < 0)) {
			d2 = (Domain *) getListIter(&iter2);
			cnt2++;
		} else {
			mergecnt++;
			cnt1++; cnt2++;
			d1 = (Domain *) getListIter(&iter1);
			d2 = (Domain *) getListIter(&iter2);
		}
	}

	if (Opt.adjOvlpRatio &&
			mergecnt >= rint(max(cnt1, cnt2) * Opt.adjOvlpRatio)) {
		mergeflag = 1;
	} else if (Opt.adjInclRatio || Opt.mincutcnt) {
		if ( cnt1 <= cnt2 && (
				(mergecnt >= rint(cnt1 * Opt.adjInclRatio)) ||
				(cnt1 - mergecnt < Opt.mincutcnt)) ) {
			if (cinfoMeanLen(cinfo1) <= Opt.minlen2) {
				mergeflag = 1;
			}
		} else if ( cnt2 <= cnt1 && (
				(mergecnt >= rint(cnt2 * Opt.adjInclRatio)) ||
				(cnt2 - mergecnt < Opt.mincutcnt)) ) {
			if (cinfoMeanLen(cinfo2) <= Opt.minlen2) {
				mergeflag = 1;
			}
		}
	}
	if (! mergeflag) {
		/* NG: do not update members */
		return 0;
	}

	/* OK: update members */

	setListIter(&iter1, cinfo1->members, 1);
	d1 = (Domain *) getListIter(&iter1);
	setListIter(&iter2, cinfo2->members, 1);
	d2 = (Domain *) getListIter(&iter2);

	while (d1 || d2) {
		if (! d2 || (d1 && cmpr_node_by_id(&d1,&d2) < 0)) {
			pushList(merged_list, d1);
			d1 = (Domain *) getListIter(&iter1);
		} else if (! d1 || (d2 && cmpr_node_by_id(&d2,&d1) < 0)) {
			pushList(merged_list, d2);
			d2 = (Domain *) getListIter(&iter2);
		} else {
			/** concatenate domains **/
			if (d1->to <= d2->from) {
				d1->to = d2->to;
				deleteDomain(d2);
			} else if (d1->from >= d2->to) {
				d1->from = d2->from;
				deleteDomain(d2);
			} else {
				fprintf(stderr, "Warning: overlapping domain: %d,%d: [%d,%d],[%d,%d]\n",
					cinfo1->root->id,cinfo2->root->id,
					d1->from,d1->to,d2->from,d2->to);
			}
			pushList(merged_list, d1);
			d1 = (Domain *) getListIter(&iter1);
			d2 = (Domain *) getListIter(&iter2);
		}
	}

	if (mergeflag) {
		newnode = dupNode(nodes, cinfo2->root);

		/* OK: update members */
		freeList(cinfo1->members);
		freeList(cinfo2->members);
		cinfo2->members = merged_list;

		concatList(cinfo2->clusters, cinfo1->clusters);
		free(cinfo1->clusters);
		cinfo1->clusters = NULL;
		cinfo1->root = newnode;

		mergeNeighborList(cinfo1->root, cinfo2->root, newnode, 1);
		mergeNeighborList(cinfo1->root, cinfo2->root, newnode, -1);
		spFlagOR(cinfo1->spflag, cinfo2->spflag, cinfo2->spflag);
		cinfo2->root = newnode;

		setListIter(&iter1, merged_list, 1);
		while (d1 = (Domain *) getListIter(&iter1)) {
			d1->root = newnode;
		}

		return 1;
	}
}

int checkOverlapClustersAll(pList *nodes, ClusterInfo **ret_cinfo,
		NodeSet *nodeset)
{
	int mergecnt;
	ClusterInfo *clust;
	Node *node;
	listIter iter;
	ClusterInfo *cinfo;
	int clustnum;
	int i, j;
	Hash *cinfoHash;

	clustnum = createClusterInfo(nodes, &cinfo, &cinfoHash);

	/** ascending order of size **/
	for (i = clustnum-1; i >= 0; i--) {
		checkOverlapCluster(&cinfo[i], cinfoHash, nodeset);
	}
	Hdestroy(cinfoHash);
	*ret_cinfo = cinfo;
	renumCinfo(cinfo, clustnum);	/* renum domains */
	return clustnum;
}
checkOverlapCluster(ClusterInfo *cinfo, Hash *cinfoHash, NodeSet *nodeset)
{
	Neighbor *nbr;
	pList *nbrlist = create_pList();
	pList *merged_list;
	Node *nbrnode;
	Node *node = cinfo->root;
	Node *newnode;
	ClusterInfo *nbrcinfo;

	if (node->left) {
		getLargeNeighbors(node->left,Opt.adjOvlpRatio,Opt.adjInclRatio,
				Opt.mincutcnt, node, 1,nbrlist);
	}
	if (node->right) {
		getLargeNeighbors(node->right,Opt.adjOvlpRatio,Opt.adjInclRatio,
				Opt.mincutcnt, node, 1,nbrlist);
	}
	while (nbr = (Neighbor *) popList(nbrlist)) {
		nbrnode = ((Node*)nbr->node1 == node) ?
			(Node*)nbr->node2 : (Node*)nbr->node1;
		cinfoHsearch(cinfoHash, nbrnode->id, &nbrcinfo);
		while (! nbrcinfo->clusters && nbrcinfo->root) {
			int origid = nbrcinfo->root->id;
			cinfoHsearch(cinfoHash, origid, &nbrcinfo);
			if(origid==nbrcinfo->root->id) break;
		}
		if (! nbrcinfo->clusters) {
			continue;
		}
		if (nbrcinfo->root->id == node->id) {
			/** already merged ? **/
			continue;
		}
		if (mergeClusters(cinfo, nbrcinfo, nodeset)) {
			HENTRY hent;
			Node *newnode = cinfo->root;
			hent.key = (char *) (newnode->id);
			hent.datum = (char *) nbrcinfo;
			HIsearch(cinfoHash, &hent, ENTER);

			break;
		}
	}
	free_pList(nbrlist);
}

checkIncludedClustersAll(NodeSet *nset)
{
	int mergecnt;
	Node *node;
	listIter iter;
	int i;

	do {
		/** ascending order of size **/
		mergecnt = 0;
		for (i = nset->nodenum-1; i>=nset->leafnum; i--){
			node=getNode(nset,i);
			if (isRoot(node) && node->child) {
				mergecnt += checkIncludedCluster(node);
			}
		}
	} while (mergecnt);
}
checkIncludedCluster(Node *node)
{
	Node *nbrnode;
	Neighbor *nbr = NULL;
	int mincnt;
	int mergecnt = 0;
	pList *nbrlist = create_pList();

	mincnt = node->cnt;

	if (node->left) {
		getLargeNeighbors(node->left, 0.0, 1.0, 0,
				node,1,nbrlist);
	}
	while (nbr = (Neighbor *) popList(nbrlist)) {
		nbrnode = ((Node *) nbr->node1 == node) ?
			(Node *) nbr->node2 : (Node *) nbr->node1;
		if (node->cnt != nbr->count) continue;
		if (concatNeighbor(node, nbrnode, 'L')) {
			mergecnt++;
		}
	}
	if (node->right) {
		getLargeNeighbors(node->right, 0.0, 1.0, 0,
				node,1,nbrlist);
	}
	while (nbr = (Neighbor *) popList(nbrlist)) {
		nbrnode = ((Node *) nbr->node1 == node) ?
			(Node *) nbr->node2 : (Node *) nbr->node1;
		if (node->cnt != nbr->count) continue;
		if (concatNeighbor(node, nbrnode, 'R')) {
			mergecnt++;
		}
	}
	free_pList(nbrlist);
	return mergecnt;
}
concatNeighbor(Node *nodeS, Node *nodeL, char dir)
{
	if(nodeS->cnt > nodeL->cnt) {
		/* swap*/
		Node *tmp = nodeS;
		nodeS = nodeL; nodeL = tmp;
		dir = (dir == 'R') ? 'L' : 'R';
	}
	if (checkNodeLen(nodeS)==2) {
		return 0;
	}

	if (Opt.DEBUG & DBG_nbrrestore) {
		printf("DEL: ");
		printNode(nodeS); printf("=%c=>",dir);
		printNode(nodeL); putchar('\n');
	}

	nodeL->len += nodeS->len;
	nodeS->len = 0;
	/* mark the smaller side node as 'MERGED' */
	setFlagNode(nodeS, NODE_MERGED);

	deleteNode(nodeS);

	/* mark the smaller node as 'DELETED' if it is on the left side */
	if (dir == 'R') {
		setFlagNode(nodeS, NODE_DELETED);
	}
	return 1;
}

