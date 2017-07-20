/*
 * DomClust: Hierarchical Clustering for Orthologous Domain Classification
 * Copyright (c) 2000-2006, Ikuo Uchiyama
 * All rights reserved.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "domclust.h"
#include "memalloc.h"
#include "plist.h"
#include "neighbor.h"
#include "spec.h"

dumpGraph(SimGraph *SimG, int argc, char **argv)
{
	int i;
	Node *node;
	char *redcheck;
	char *sp;
	FILE *ofp;
	int spnum = getSPnum();
	if (Opt.dumpfile) {
		if ( (ofp = fopen(Opt.dumpfile, "w")) == NULL) {
			fprintf(stderr, "Can't open dump file\n");
			exit(1);
		}
	} else {
		ofp = stdout;
	}
	fprintf(ofp, "## dumpfile for domclust\n");
	fprintf(ofp, "AC %d\n", argc);
	for (i = 1; i < argc; i++) {
		if (strcmp(argv[i], "-o10")==0) {
		} else {
			fprintf(ofp, "AV %s\n", argv[i]);
		}
	}
	fprintf(ofp, "//\n");
	for (i = 0; i < spnum; i++) {
		sp = getSPname(i);
		fprintf(ofp, "SP %s %lf\n", sp, getSPweight(i));
	}
	fprintf(ofp, "//\n");
	fprintf(ofp, "G %d %d", SimG->nodes->nodenum, SimG->nodes->leafnum);
	fprintf(ofp, " %d\n", SimG->edges->edgenum);
	for (i = 0; i < SimG->nodes->nodenum; i++) {
		node = getNode(SimG->nodes, i);
		dumpNode1(ofp, node);
	}
	fprintf(ofp, "//\n");
	for (i = 0; i < SimG->nodes->nodenum; i++) {
		node = getNode(SimG->nodes, i);
		dumpNode2(ofp, node);
	}
}
dumpNode1(FILE *ofp, Node *node)
{
	fprintf(ofp, "N %s %d ", node->name, node->id);
	fprintf(ofp, "%d %d ", node->flag, node->dir);
	fprintf(ofp, "%d %d ", node->from, node->len);
	fprintf(ofp, "%d %d ", node->consreg.from, node->consreg.to);
	fprintf(ofp, "%d ", node->totlen);
	fprintf(ofp, "%d %d ", node->brk.from, node->brk.to);
	fprintf(ofp, "%d %d ", node->newreg.from, node->newreg.to);
	fprintf(ofp, "%d\n", node->cnt);
}
dumpNode2(FILE *ofp, Node *node)
{
	listIter iter;
	Domain *dom;
	Neighbor *nbr;
	fprintf(ofp, "N2 %d\n", node->id);
	fprintf(ofp, "S ");
	dump_specFlag(ofp, node->spflag);
	fprintf(ofp,"\n");
	if (node->parent) {
		fprintf(ofp, "P %d\n", node->parent->id);
	}
	if (node->parentL) {
		fprintf(ofp, "PL %d\n", node->parentL->id);
	}
	if (node->parentM) {
		fprintf(ofp, "PL %d\n", node->parentM->id);
	}
	if (node->parentR) {
		fprintf(ofp, "PR %d\n", node->parentR->id);
	}
	if (node->child) {
		dumpEdge1(ofp, node->child);
	}
	if (node->left) {
		dumpNeighborList(ofp, node->left, "NB", node, -1);
	}
	if (node->right) {
		dumpNeighborList(ofp, node->right, "NB", node, 1);
	}
}
dumpEdge1(FILE *ofp, Edge *edge)
{
	fprintf(ofp, "E %d %d %d %lf %lf",
			edge->id, edge->node1->id, edge->node2->id,
			edge->dist, edge->score);
	fprintf(ofp, " %d %d %d %d", edge->ali1->from, edge->ali1->to,
			edge->ali2->from, edge->ali2->to);
	fprintf(ofp, " %d %d %d\n", edge->flag, edge->dir, edge->connect);
}
restoreGraph(char *filename, SimGraph *SimG,
		int curr_argc, char **curr_argv)
{
	FILE *fp;
	char buf[BUFSIZ];
	int ln = 0, statflag = 0;
	int nodenum,leafnum,edgenum,spnum;
	int p1,p2,p3;
	char name[NAMELEN], *currname, header[10];
	int id,flag,dir,from,len, cnt,status;
	int consfrom, consto;
	int totlen;
	int id1, id2,from1,to1,from2,to2,connect;
	double score,dist,weight;
	NameHash *nhash;
	NodeSet *nodes;
	EdgeSet *edges;
	Node *node, *node1, *node2;
	Edge *edge;
	specFlag spflag;
	Region ali1, ali2;
	Neighbor *nbr;
	int argc = 0, agn = 0;
	char **argv;
	Region consreg;

	nodes = createNodes();
	if ((fp = fopen(filename, "r")) == NULL) {
		fprintf(stderr, "Can't open file\n");
		exit(1);
	}
	while (fgets(buf, BUFSIZ, fp) != NULL) {
		ln++;
		if (buf[0] == '#') {
			continue;
		} else if (buf[0] == '/') {
			if (statflag == 0) {
				getargs(agn+1, argv);
				getargs(curr_argc, curr_argv);
			}
			statflag++;
		} else if (statflag == 0) {
			if (strncmp(buf, "AC ", 3) == 0) {
				if (sscanf(buf, "AC %d", &argc) < 1) {
					formatError(ln, "argc", buf);
				}
				agn = 0;
				if ((argv = (char **)malloc(sizeof(char*)*argc))
						== NULL) {
					allocError("restore");
				}
			} else if (strncmp(buf, "AV ", 3) == 0) {
				buf[strlen(buf)-1] = '\0';
				if (++agn >= argc) {
					formatError(ln, "Too many args", buf);
				}
				argv[agn] = strdup(&buf[3]);
			}
		} else if (statflag == 1) {
			sscanf(buf, "SP %s%lf", name, &weight);
			id = getSPid(name);
			setSPweight(id, weight);
		} else if (statflag == 2) {
			if (buf[0] == 'G') {
				sscanf(buf, "G %d%d%d%d", &nodenum,&leafnum,
							&edgenum);
 				nhash = initNames(nodenum);
				edges = createEdges(nodenum);
				nodes->leafnum = leafnum;
				continue;
			} else {
				sscanf(buf, "N %s%d%d%d%d%d%d%d%d%d%d%d%d%d",
					name,&id,&flag,&dir,&from,&len,
					&consfrom, &consto, &totlen,
					&from1,&to1,&from2,&to2,&cnt);
				currname = addName(nhash, name, id);
				consreg.from=consfrom; consreg.to=consto;
				node = addNode(nodes, currname, cnt, len,
					&consreg, totlen,
					NULL, 0, spflag);
				node->flag = (NodeFlag)flag;
				node->dir = dir;
				setSeqReg(&node->brk, (SeqPos)from1, (SeqPos)to1);
				setSeqReg(&node->newreg, (SeqPos)from2, (SeqPos)to2);
				if (node->id != id){
					fprintf(stderr, "Warning: node id mismatch\n");
					node->id = id;
				}
			}
		} else if (strncmp(buf, "N2 ", 3) == 0) {
			sscanf(buf, "N2 %d", &id);
			node = getNode(nodes, id);
		} else if (buf[0] == 'S') {
			restore_specFlag(buf+2, node->spflag);
		} else if (buf[0] == 'P') {
			if (sscanf(buf, "%s%d", header,&id) < 2) {
				formatError(ln, "E", buf);
			}
			if (header[1]=='R') {
				node->parentR = getNode(nodes, id);
			} else if (header[1]=='L') {
				node->parentL = getNode(nodes, id);
			} else if (header[1]=='M') {
				node->parentM = getNode(nodes, id);
			} else {
				node->parent = getNode(nodes, id);
			}
		} else if (strncmp(buf, "E ", 2) == 0) {
			if (sscanf(buf, "E %d%d%d%lf%lf%d%d%d%d%d%d%d",
				&id, &id1, &id2, &dist, &score,
				&from1,&to1,&from2,&to2,
				&flag,&dir,&connect) < 12) {
				formatError(ln, "E", buf);
			}
			node1 = getNode(nodes, id1);
			node2 = getNode(nodes, id2);
			setSeqReg(&ali1, (SeqPos)from1, (SeqPos)to1);
			setSeqReg(&ali2, (SeqPos)from2, (SeqPos)to2);
			edge = addEdgeWithScore(edges, node1, node2,
				&ali1, &ali2,
				(Dist) dist, (Dist) score, (ConnCount) connect,
				(signed char) dir);
			edge->id = id;	/* reset edge ID */
			node->child = edge;
		} else if (strncmp(buf, "NB ", 3) == 0) {
			int dir1, dir2;
			if (sscanf(buf, "%s%d%d%d%d%d%d",
					header,&id1,&id2,&cnt,
					&dir1,&dir2,&status) < 7) {
				formatError(ln, "NB", buf);
				exit(1);
			}
			node1 = getNode(nodes, id1);
			node2 = getNode(nodes, id2);
			addNeighborRestore(node1, node2, dir1, dir2,
				cnt, status);
		}
	}
	fclose(fp);
	SimG->nodes = nodes;
	SimG->edges = edges;
	SimG->nhash = nhash;
}

formatError(int ln, char *str, char *buf)
{
	fprintf(stderr, "restore failed: line %d: %s\n>>> %s", ln, str, buf);
	exit(1);
}
