/*
 * DomClust: Hierarchical Clustering for Orthologous Domain Classification
 * Copyright (c) 2000-2006, Ikuo Uchiyama
 * All rights reserved.
 */
#include <stdio.h>
#include "domclust.h"
#include "util.h"

static Alloc_Object *domobj;
static int checkNodeDomains(Node *node, SeqPos from, SeqPos to,
		SeqPos savefrom, SeqPos saveto,
		Region *anchor, SeqPos seqlen);
static int getDomainMark(pList *domains, Node *rootnode, Domain **dom, char mark);

createDomainObj()
{
	domobj = init_alloc_object(sizeof(Domain), NODE_BLKSIZ);
}
createDomain(Node *node)
{
	Domain *dom;
	dom = (Domain *) memalloc(domobj);
	dom->root = NULL;
	dom->leaf = node;
	dom->from = 1;
	dom->to = node->len;
	dom->mark = 0;
}
deleteDomain(Domain *deldom)
{
	Node *leaf = deldom->leaf;
	Domain *dom;
	listIter iter;
	if (leaf==NULL) return;
	setListIter(&iter, leaf->domains, 1);
	while (dom = (Domain *) getListIter(&iter)) {
		if (dom == deldom) {
			delCurrElem(&iter);
		}
	}
}
getDomainNoMark(pList *domains, Node *rootnode, Domain **dom)
{
	return getDomainMark(domains, rootnode, dom, 0);
}
getDomain(pList *domains, Node *rootnode, Domain **dom)
{
	return getDomainMark(domains, rootnode, dom, 1);
}
getDomainMark(pList *domains, Node *rootnode, Domain **dom, char mark)
{
	listIter iter;
	Domain *dm;
	int domn = 0, i = 1;
	if (domains == NULL) {
		return 0;
	}
	*dom = NULL;
	setListIter(&iter, domains, 1);
	while (dm = (Domain *) getListIter(&iter)) {
		if (dm->root == rootnode && dm->mark == 0) {
			*dom = dm;
			dm->mark = mark;
			dm->num = domn = i;
			break;
		}
		i++;
	}
	/** domn might be 0 when the domain has been discarded **/

	return domn;
}
clearDomainMark(pList *domains)
{
	listIter iter;
	Domain *dm;
	if (domains == NULL) {
		return 0;
	}
	setListIter(&iter, domains, 1);
	while (dm = (Domain *) getListIter(&iter)) {
		dm->mark = 0;
	}
}
clearDomainMarkAll(NodeSet *nodes)
{
	int i;
	Node *node;
	for (i = 0; i < nodes->leafnum; i++) {
		node = getNode(nodes, i);
		clearDomainMark(node->domains);
	}
}

static Node *leafnode; 
static int domid;
static int save_from;

int tmpflag;

/* For each sequence, determine the domain boundaries */
checkDomains(NodeSet *nodes)
{
	Node *node;
	int i;

	createDomainObj();
	for (i = 0; i < nodes->leafnum; i++) {
		node = getNode(nodes, i);
		leafnode = node;
		domid = 1;
		leafnode->domains = (pList *) create_pList();
		save_from = 0;
		checkNodeDomains(node, (SeqPos) 1, node->len,
				(SeqPos) 1, node->len, NULL, node->len);
	}
}

/* from, to: anchor points on on the current sequence */
/* savefrom, saveto: saved domain boundary on the current sequence */
checkNodeDomains(Node *node, SeqPos from, SeqPos to,
		SeqPos savefrom, SeqPos saveto, Region *anchor, SeqPos seqlen)
{
	char flag = 0;
	SeqPos from1, to1;
	SeqPos savefrom1, saveto1;
	SeqPos len = to - from + 1;
	Domain *dom;
	Edge *e, *child;
	SeqPos alilen;
	Region *reg, thisNode, currReg;

	/* currReg: coordinate on the examined sequence */
	/* thisNode: coordinate on this node */
	/*     they are aligned to each other */
	if (anchor != NULL) {
		copyReg(&thisNode, anchor);
	} else {
		thisNode.from = 1; thisNode.to = node->len;
	}
	currReg.from = from; currReg.to = to;
	
	if (isRoot(node)) {
		int tmplen;
		if (Opt.outstyle == DOMAINOUT) {
			printf("%s(%d) %d %d %d %d\n",
				leafnode->name, domid++,
				node->id, from, to, node->len);
		}
		if (save_from > 0) {
			savefrom = save_from;
		}
		savefrom = max(savefrom, 1);
		saveto   = min(saveto, seqlen);
		tmplen = saveto - savefrom + 1;
		if (tmplen < Opt.minlen && tmplen < seqlen * Opt.ovlpratio2) {
			/** discard short domains **/
			return 0;
		}

		dom = (Domain *) memalloc(domobj);
		if (save_from > 0) {
			save_from = 0;
		}
		dom->from = savefrom;
		dom->to = saveto;
		dom->root = node;
		dom->leaf = leafnode;
		dom->mark = 0;
		pushList(leafnode->domains,dom);
	} else if (node->flag & NODE_MERGED) {
		if (node->flag & NODE_DELETED) {
			/* merge this node into the right (next) node */
			if (save_from == 0) {
				save_from = savefrom;
			}
		} else {
			/* merge this node into the left (previous) node */
			dom = (Domain *) getListIdx(leafnode->domains, -1);
			if (dom) {
				/** merged region **/
				dom->to = saveto;
			}
		}
	} else if (nodeDeleted0(node)) {
		/* node that have been already cut--do nothing */
	} else {
		/***** ------ParentL------- *****/
		if (node->parentL) {
			if (definedPos(node->brk.from)) {
				saveto1 = transformSeqPos(node->brk.from-1, &thisNode, &currReg);
			} else if (definedPos(node->brk.to)) {
				saveto1 = transformSeqPos(node->brk.to, &thisNode, &currReg);
			} else{
				saveto1 = saveto;
			}
			from1 = from;
			to1 = transformSeqPos(node->newreg.from-1, &thisNode, &currReg);

			checkNodeDomains(node->parentL, from1, to1,
				savefrom, saveto1, NULL, seqlen);
			flag = 1;
		}

		/***** ------ParentM------- *****/
		if (node->parent) {
			if (Opt.domBoundary == B_BREAK) {
				/* break point */
				reg = &(node->brk);
			} else {
				/* aligned region */
				reg = &(node->newreg);
			}
			if (definedPos(node->brk.from)) {
				savefrom1 = transformSeqPos(node->brk.from, &thisNode, &currReg);
			} else {
				savefrom1 = savefrom;
			}

			if (definedPos(node->brk.to)) {
				saveto1 = transformSeqPos(node->brk.to, &thisNode, &currReg);
			} else {
				saveto1 = saveto;
			}

			/* mapping: node->newreg <-> parent->consreg */
			from1 = transformSeqPos(node->newreg.from, &thisNode, &currReg);
			to1 = transformSeqPos(node->newreg.to, &thisNode, &currReg);

			checkNodeDomains(node->parent, from1, to1,
					savefrom1, saveto1,
					&(node->parent->consreg), seqlen);


			if (definedPos(node->newreg2.from)
					&& definedPos(node->newreg2.to)) {
				int from2, to2;
				Region *reg2;

				if (Opt.domBoundary == 1) {
					/* break point */
					reg2 = &(node->brk2);
				} else 	{
					/* aligned region */
					reg2 = &(node->newreg2);
				}
				
				from2 = transformSeqPos(node->newreg2.from, &thisNode, &currReg);
				to2 = transformSeqPos(node->newreg2.to, &thisNode, &currReg);

if (node->parentM){
				checkNodeDomains(node->parentM, to1+1, from2-1,
					to1+1, from2-1,
					NULL, seqlen);
}
				checkNodeDomains(node->parent, from2, to2,
					from2, to2,
					&(node->parent->consreg), seqlen);
			}
			flag = 1;
		}

		/***** ------ParentR------- *****/
		if (node->parentR) {
			if (definedPos(node->brk.to)) {
				savefrom1 = transformSeqPos(node->brk.to+1,&thisNode,&currReg);
			} else {
				savefrom1 = savefrom;
			}
			from1 = transformSeqPos(node->newreg.to+1,&thisNode,&currReg);
			to1 = to;


			checkNodeDomains(node->parentR, from1, to1,
					savefrom1, saveto, NULL, seqlen);
			flag = 1;
		}
		if (! flag) {
			/* never come here */
			fprintf(stderr, "No parent at unrooted node: %s (%d)\n",
				node->name,node->id);
		}
	}
}

printDomain(Domain *dom)
{
	printNode(dom->leaf);
	printf(" (%d) %d %d\n",dom->num, dom->from, dom->to);
}
