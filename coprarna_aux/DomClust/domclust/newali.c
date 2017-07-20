/*
 * DomClust: Hierarchical Clustering for Orthologous Domain Classification
 * Copyright (c) 2000-2006, Ikuo Uchiyama
 * All rights reserved.
 */

#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include "domclust.h"
#include "vararray.h"
#include "util.h"
#include "seqreg.h"

#define MAXNEWDOM 5

post_calNewAli(Node *n1, Node *n2, Node *n3, Node *n3_2, Node **newnodes,
	EdgeSet *edges, Edge *best, Edge *e1, Edge *e2, Edge **newedges,
	Region *newali, Region *newali3,
	char *hitflag, pList *newEdges)
{
	int nn;
	Node *newn;
	Dist newdist, newscore;
	ConnCount newconn;
	signed char edir1, edir2, newedir;
	Count cnt1 = n1->cnt;
	Count cnt2 = n2->cnt;
	Dist score1, score2, dist1, dist2;
	ConnCount conn1, conn2;
	Edge *newedge;
	int dir;

	if (e1) {
		score1 = e1->score;
		dist1 = e1->dist;
		edir1 = e1->dir;
		conn1 = e1->connect;
	} else {
		score1 = Opt.missscore;
		dist1 = Opt.missdist;
		edir1 = 0;
		conn1 = 0;
	}
	if (e2) {
		score2 = e2->score;
		dist2 = e2->dist;
		edir2 = e2->dir;
		conn2 = e2->connect;
	} else {
		score2 = Opt.missscore;
		dist2 = Opt.missdist;
		edir2 = 0;
		conn2 = 0;
	}

	/** create new edges (similarities) for each new node (domain) **/
	    for (nn = 0; nn < MAXNEWDOM; nn++) {
		newedges[nn] = NULL;
		if (! (newn = newnodes[nn])) {
			continue;
		}

		if (Opt.DEBUG & DBG_cluster) {
			if (e1 && e2) {
				printf("%d %d\n",n3->id, n3_2->id);
			} else if (e1) {
				printf("%d -\n",n3->id);
			} else {
				printf("- %d\n",n3->id);
			}
		}
		if (Opt.DEBUG & DBG_basic) {
			printf("NEWDOM: %d: OVLP-", nn);
		}
		if (hitflag[nn] == 0) {
			/* no overlap */
			if (Opt.DEBUG & DBG_basic) {
				printf("NO overlap -- skip\n");
			}
			continue;
		}
		if (Opt.DEBUG & DBG_basic) {
			printf("YES: %d,%d\n",newali[nn].from,newali[nn].to);
		}

/*
		assert (newali[nn].to <= newnodes[nn]->len);
*/


	/** update the distance and the number of connections **/
		calNewDist(nn, cnt1, cnt2, dist1, dist2, score1, score2,
			conn1, conn2, &newdist, &newscore, &newconn);

		if (Opt.revMatch) {
			/* check direction incongruence for DNA sequence
				matches  (edir1 = bestdir * edir2) */
			if ((int)best->dir * (int)edir1 * (int)edir2 < 0) {
				fprintf(stderr, "contradictory direction %d,%d,%d %s,%s,%s\n",
					(int)best->dir,(int)edir1,(int)edir2,
					n1->name,n2->name,n3->name);
			}
		}
		if (edir1) {
			newedir = edir1;
		} else {
			newedir = best->dir * edir2;
		}

		if (Opt.DEBUG & DBG_cluster) {
			printf("OK: %s(%d)={%s,%s}[%d]-%s,%f,%f  [%d,%d]\n",
		  newn->name,newn->id,n1->name,n2->name,nn,n3->name, newdist,
				newscore, newali3[nn].from, newali3[nn].to);
			printf("-------------\n");
		}

		if (Opt.DEBUG & DBG_cluster) {
			printf("NewAlign: %d\n", nn);
			printf("newali : "); printReg(&newali[nn]);
			printf("newali3: "); printReg(&newali3[nn]);
			printf("---\n");
		}

		/** create a new edge **/
		newedge = addEdgeWithScore(edges, newn, n3,
				&newali[nn], &newali3[nn], newdist, newscore,
				newconn, newedir);
		if (Opt.DEBUG & DBG_basic) {
			printf("NewEdge:"); printEdge(newedge);
		}

		newedges[nn] = newedge;
		pushList(newEdges, newedge);
	    }
        if (Opt.neighbor) {
                /** add neighboring links to the new edges **/
                addLinksForNewEdge(newedges, e1, e2);
        }
}

/* calculate new alignment positions using three pair-alignment positions
	between three sequences and break points on the new sequence */
calNewAli(
    /* aliXY: aligment region on sequence X against sequence Y */
	Region *ali12,		/* region on seq1 aligned against seq2 */
	Region *ali21,		/* region on seq2 aligned against seq1 */
	Region *ali13,		/* region on seq1 aligned against seq3 */
	Region *ali31,		/* region on seq3 aligned against seq1 */
	Region *ali23,		/* region on seq2 aligned against seq3 */
	Region *ali32,		/* region on seq3 aligned against seq2 */
	Region *brk1, Region *brk2,	/* break points on seq1 and seq2 */
	Region *newseq,
	Region *newreg,		/* new sequence by merging seq1 and seq2 */
	Region *newconsreg,	/* aligned reg. (seq1 vs seq2) on newreg */
	Region *newali,		/* output: aligned reg on newreg against seq3 */
	Region *newali3,	/* output: aligned reg on seq3 against newreg */
	char *hitflag,
	Count cnt1, Count cnt2, Count cnt3, Node *n1, Node *n2, Node *n3,
	Edge *e1, Edge *e2)
{
	Region brk3, tmp1, tmp2, tmp3;
	Region tmp_newali;
	Region tmp13, tmp23, tmp31, tmp32, tmpdmy;
	int i, j;
	int ali_consistence = 0;

	if (Opt.DEBUG & DBG_cluster){
		printReg1(ali12); printReg1(ali21);
		printReg1(newreg); printReg1(newconsreg);
		printReg1(ali13); printReg1(ali31);
		printReg1(ali23); printReg1(ali32);
		printReg1(brk1); printReg1(brk2);
		printf("----\n");
	}

	/* calculate brk3 (break point on seq3) from brk1 and/or brk2 */

	if (ali13 && ali23) {
		ali_consistence = alignmentConsistency(ali12,ali21,
				ali13,ali31,ali23,ali32);
		transformReg2(brk1, &tmp1, ali13, ali31);
		transformReg2(brk2, &tmp2, ali23, ali32);
		if (ali_consistence) {
			/* the break points on seq1,2 are mapped onto seq3
				and are averaged */
			meanSeqReg(&tmp1, cnt1, &tmp2, cnt2, &brk3);
		} else if (BETTER(MEASURE(e1),MEASURE(e2))) {
			/* e1 precedes e2 */
			brk3.from = definedPos(tmp1.from) ? tmp1.from : tmp2.from;
			brk3.to = definedPos(tmp1.to) ? tmp1.to : tmp2.to;
		} else {
			/* e2 precedes e1 */
			brk3.from = definedPos(tmp2.from) ? tmp2.from : tmp1.from;
			brk3.to = definedPos(tmp2.to) ? tmp2.to : tmp1.to;
		}
	} else if (ali13) {
		transformReg2(brk1, &brk3, ali13, ali31);
	} else if (ali23) {
		transformReg2(brk2, &brk3, ali23, ali32);
	}

	if (Opt.DEBUG & DBG_cluster){
		printReg1(&brk3);
		printf("----\n");
	}

	if (hitflag[0]) {
		int ov1,ov2;

	/* calculate alignment region on newseq against seq3 (newali[0]) */
		if (ali13 && ali23) {
		/* alignment regions in seq1,2 are mapped onto newseq */

			transformReg2(ali13, &tmp13, ali12, newconsreg);
			transformReg2(ali23, &tmp23, ali21, newconsreg);

			ov1 = ovlpReg(ali31, &brk3, &tmp31, NULL, 1);
			ov2 = ovlpReg(ali32, &brk3, &tmp32, NULL, 1);

			if (ov1 != 0 || ov2 != 0) {
				/** alignment inconsistency: no overlap **/
				hitflag[0] = 0;
				/* break the reg0 block */
				goto check_reg0_end;
			}

			if (! ali_consistence) {
				/* non-order preserving */
				/* internal repeat ? */
				/* choose alignment with better score */

				if (Opt.DEBUG & DBG_cluster){
					printf("internal repeat ??\n");
					printReg1(&tmp31); printReg1(&tmp32);
					printReg1(ali31); printReg1(ali32);
					putchar('\n');
				}

				if (BETTER(MEASURE(e1),MEASURE(e2))) {
					copyReg(&tmp_newali, &tmp13);
					copyReg(&newali3[0], &tmp31);
				} else {
					copyReg(&tmp_newali, &tmp23);
					copyReg(&newali3[0], &tmp32);
				}
			} else {
				Region tmpmeanM, tmpmean3;
				/* the alignments are consistent */
				/* union of the two segments */
				mergeReg(&tmp13, &tmp23, &tmp_newali);
				mergeReg(&tmp31, &tmp32, &newali3[0]);


			/**** a heuristic to avoid extremely erroneous case ***/
				meanSeqReg(&tmp13, cnt1, &tmp23, cnt2, &tmpmeanM);
				meanSeqReg(&tmp31, cnt1, &tmp32, cnt2, &tmpmean3);
				if (abs(tmp31.from-tmp32.from)<=10) {
					tmp_newali.from = tmpmeanM.from + newali3[0].from - tmpmean3.from;
				}
				if (abs(tmp31.to-tmp32.to)<=10) {
					tmp_newali.to = tmpmeanM.to + newali3[0].to - tmpmean3.to;
				}
				
			}

		} else if (ali13) {
			transformReg2(ali13, &tmp_newali, ali12, newconsreg);
			ovlpReg(ali31, &brk3, &newali3[0], NULL, 1);
		} else if (ali23) {
			transformReg2(ali23, &tmp_newali, ali21, newconsreg);
			ovlpReg(ali32, &brk3, &newali3[0], NULL, 1);
		}

		/** tmp_newali ^ newreg ==> newali[0] **/
		if (ovlpReg(&tmp_newali, newreg, &newali[0], NULL, 1) != 0) {
			/** alignment inconsistency: no overlap **/
			if (Opt.DEBUG & DBG_basic){
				printf("No overlap -- discarded\n");
				printReg1(&tmp_newali); printReg1(newreg);
			}
			hitflag[0] = 0;
			goto check_reg0_end;	/* break the reg0 block */
		}

		if (Opt.DEBUG & DBG_basic){
			printf("---- new alignments ---\n");
			printReg2(newreg, "newreg");
			printReg2(&newali[0], "newali");
			printReg2(&newali3[0], "newali3");
			printf("----\n");
		}

	}

	/** check overlap lengths **/
    check_reg0_end:
	if (hitflag[0]) {
		if (Opt.DEBUG & DBG_basic){
			printReg1(&newali[0]); printReg1(newreg);
			printReg1(ali13); printReg1(ali23);
		}

		/** length check */
		if ((ali31 && ! regLenCheck(&newali3[0], newreg, ali31)) || 
		    (ali32 && ! regLenCheck(&newali3[0], newreg, ali32)) ) {
			hitflag[0] = 0;
		}
		if ((ali13 && ! regLenCheck(&newali[0], newreg, ali13)) || 
		    (ali23 && ! regLenCheck(&newali[0], newreg, ali23)) ) {
			hitflag[0] = 0;
		}
	}
	if (hitflag[1]) {
		setSeqReg(&newali3[1], ali31->from, min(ali31->to,brk3.from-1));
		if ( ! regLenCheck(&newali3[1], newreg, ali31) ) {
			hitflag[1] = 0;
		} else {
			setSeqReg(&newali[1], ali13->from, min(ali13->to,brk1->from-1));
			if ( ! regLenCheck(&newali[1], newreg, ali13) ) {
				hitflag[1] = 0;
			}
		}
	}
	if (hitflag[2]) {
		setSeqReg(&newali3[2], ali32->from, min(ali32->to,brk3.from-1));
		if ( ! regLenCheck(&newali3[2], newreg, ali32) ) {
			hitflag[2] = 0;
		} else {
			setSeqReg(&newali[2], ali23->from, min(ali23->to,brk2->from-1));
			if ( ! regLenCheck(&newali[2], newreg, ali23) ) {
				hitflag[2] = 0;
			}
		}
	}
	if (hitflag[3]) {
		setSeqReg(&newali3[3], max(ali31->from,brk3.to+1), ali31->to);
		if ( ! regLenCheck(&newali3[3], newreg, ali31) ) {
			hitflag[3] = 0;
		} else {
			setSeqReg(&newali[3],
			 max(ali13->from - brk1->to, 1), ali13->to - brk1->to);
			if ( ! regLenCheck(&newali[3], newreg, ali13) ) {
				hitflag[3] = 0;
			}
		}
	}
	if (hitflag[4]) {
		setSeqReg(&newali3[4], max(ali32->from, brk3.to+1), ali32->to);
		if ( ! regLenCheck(&newali3[4], newreg, ali32) ) {
			hitflag[4] = 0;
		} else {
			setSeqReg(&newali[4],
			  max(ali23->from - brk2->to, 1), ali23->to - brk2->to);
			if ( ! regLenCheck(&newali[4], newreg, ali23) ) {
				hitflag[4] = 0;
			}
		}
	}
	for (i = 0; i < MAXNEWDOM; i++) {
		if (hitflag[i]) {
			if (newali[i].from < 0) {
				/** should not come here **/
				if (newali[i].to < 0) {
					newali[i].to = 1;
					hitflag[i] = 0;
				}
				printf("Warning: minus position: %d, %d\n",
					newali[i].from, i);
				newali[i].from = 1;
			}
		}
	}

}

regLenCheck(Region *ali, Region *reg1, Region *reg2)
{
	int len = regLen(ali);
	int reglen1 = regLen(reg1);
	int reglen2 = regLen(reg2);

	return (len >= overlapCutoff(reglen1,reglen2));
}

/** averaging the distances, scores and the number of connections **/
calNewDist(int nn, Count cnt1, Count cnt2, Dist dist1, Dist dist2,
		Dist score1, Dist score2, ConnCount conn1, ConnCount conn2,
		Dist *newdist, Dist *newscore, ConnCount *newconn)
{
	int flag1 = (nn == 0 || nn == 1 || nn == 3);
	int flag2 = (nn == 0 || nn == 2 || nn == 4);

	*newdist = (Dist) (dist1 * cnt1 * flag1 + dist2 * cnt2 * flag2)
			/ (Dist) (cnt1 * flag1 + cnt2 * flag2);
	*newscore =(Dist) (score1 * cnt1 * flag1 + score2 * cnt2 * flag2)
			/ (Dist) (cnt1 * flag1 + cnt2 * flag2);
	*newconn = ((int) conn1 + conn2 > MAXUCHR)
			? MAXUCHR : conn1 + conn2;
}

