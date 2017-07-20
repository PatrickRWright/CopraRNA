/*
 * DomClust: Hierarchical Clustering for Orthologous Domain Classification
 * Copyright (c) 2000-2006, Ikuo Uchiyama
 * All rights reserved.
 */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <ctype.h>
#include "domclust.h"
#include "readfile.h"


char filename[100] = "stdin";
char genefilename[100];
char *sptreefile;
char *restorefile;
SimGraph SimG;
int getargs(int argc, char **argv);

/** setting default parameters */
defaultOpt() {
	Opt.distscale = DISTSCALE;
	Opt.logdistscale = 0;
	Opt.DEBUG = 0;
	Opt.VERBOSE = 0;
	Opt.verbose_step1 = 0;
	Opt.verbose_step2 = 0;
	Opt.minsp = 2;
	Opt.minent = 2;
	Opt.delete_small = 1;
	Opt.outstyle = NORMALOUT;

	Opt.minlen = 40;
	Opt.minlen2 = 400; /* used in checkIncludedCluster() */

	/* used for overlap condition defined in overlapCheck() */
	/* ($ovlen > Opt.ovlpratio * $len1 && $ovlen > Opt.minlen)
	    || (($ovlen > Opt.ovlpratio2 * $len2) && $ovlen > Opt.minovlp)
		($len1 < $len2) */
	Opt.minovlp = 50;
	Opt.ovlpratio = 0.6;
	Opt.ovlpratio2 = 0.3;

	Opt.sim = 1;	/* 1..similarity, 0..distance; default=similarity */
	Opt.phylocutratio = 0.5;
	Opt.phylocutratio2 = 0.2;	/* checking distdiffcut */
	Opt.mincutcnt = 1;	/* minimum number of seq. for dom. cut */

	Opt.min_alilen = 250;

	Opt.distdiffcut = 0;
	Opt.domBoundary = B_BREAK;
	Opt.neighbor = 0;
	Opt.dpWin = 5;
	Opt.iniGap = 20.0;
	Opt.extGap = 0.2;
	Opt.skipGap = 0.1;

	Opt.nbrScoreRatio = 0.3;
	Opt.nbrScoreDecay = 1;
	Opt.nbrConnRatio = 0.5;

	Opt.adjInclRatio = 1.0;
	Opt.adjOvlpRatio = 0.0;

	Opt.noextend=0;
}

main(int argc, char **argv)
{
	FILE *fp;
	Dist dist, score;
	char name1[NAMELEN], name2[NAMELEN], lenstr[NAMELEN];
	int id1, id2, domnum;
	static char buf[BUFSIZ];
	int dir = 0;
	Node *node, *node1, *node2;
	Edge *edge;
	Cluster *sclust;
	int i;
	int from1, to1, from2, to2;
	Region ali1, ali2;
	int tmplen;
	int errflg;
	int scannum;
	SelFile selfile;
	SelData seldata;
	int minlen, alilen;
	int cnt = 0;

	defaultOpt();
	getargs(argc, argv);
	if (restorefile) {
		restoreGraph(restorefile, &SimG, argc, argv);
	}
	if (Opt.sim) {
		if (! Opt.cutoff) {
			Opt.cutoff = MINSCORE;
		}
		if (! Opt.missscore) {
			if (Opt.missratio) {
				Opt.missscore = Opt.cutoff * Opt.missratio;
			} else {
				Opt.missscore = Opt.cutoff * MISSRATIO;
			}
		} 
	} else {
		if (! Opt.cutoff) {
			Opt.cutoff = MAXDIST;
		}
		if (! Opt.missdist) {
			if (Opt.missratio) {
				Opt.missdist = Opt.cutoff / Opt.missratio;
			} else {
				Opt.missdist = Opt.cutoff / MISSRATIO;
			}
		}
	}
	if (sptreefile) {
		if (*sptreefile == '(') {
			parse_spinfo(sptreefile);
		} else {
			readSPfile(sptreefile);
		}
	} else if (Opt.treecheck) {
		fprintf(stderr, "No tree is specified "
			"-- the treecheck option is canceled\n");
		Opt.treecheck = 0;
	}
	if (restorefile) {
		preproc_for_SpecInfo();
		setOutGroup(Opt.outgroupStr, Opt.outgroupRev);
		setSpMask(Opt.spmaskStr);
		outputResults(&SimG, argc, argv);
		exit(0);
	}

	SimG.nodes = createNodes();
	SimG.nhash = initNames(INIT_NODENUM);

	if (! *genefilename) {
		FILE *fp;
		strcpy(genefilename, filename);
		read_genefile(filename, &SimG, &fp);
		open_seldata(filename, &selfile, fp);
	} else {
		read_genefile(genefilename, &SimG, NULL);
		open_seldata(filename, &selfile, NULL);
	}


	if (! SimG.nodes->nodenum) {
		fprintf(stderr, "No node read (Is the gene file correct?)\n");
	}

	SimG.edges = createEdges(SimG.nodes->nodenum);

	while ((scannum = read_seldata(&selfile, &seldata)) > 0) {
		if (scannum < MIN_HOMDATA) {
			/** read error **/
			continue;
		}
		if (Opt.VERBOSE) {
			if (++cnt % Opt.verbose_step1 == 0) {
				fprintf(stderr, "read %d\n", cnt);
			}
		}
		if (Opt.sim) {
			if (seldata.score < Opt.missscore) continue;
		} else {
			if (seldata.dist > Opt.missdist) continue;
		}

		id1 = getNameID(SimG.nhash, seldata.name1);
		id2 = getNameID(SimG.nhash, seldata.name2);
		errflg = 0;
		if ((node1 = getNode(SimG.nodes, id1)) == NULL) {
			fprintf(stderr, "not found: %s\n", seldata.name1);
			errflg = 1;
		}
		if ((node2 = getNode(SimG.nodes, id2)) == NULL) {
			fprintf(stderr, "not found: %s\n", seldata.name2);
			errflg = 1;
		}
		if (errflg) {
			if (Opt.skipErrEnt) {
				continue;
			} else {
				exit(1);
			}
		}

		if (node1->len < node2->len) {
			minlen = node1->len;
			alilen = (int)seldata.to1 - seldata.from1 + 1;
		} else {
			minlen = node2->len;
			alilen = (int)seldata.to2 - seldata.from2 + 1;
		}
		if (Opt.covfilt) {
			if ((float) alilen / minlen < Opt.covfilt) {
				continue;
			}
		}

		dist = (Dist) seldata.dist; score = (Dist) seldata.score;
#ifdef LARGE
		dir = seldata.dir;
#endif
		if (! Opt.sim) {
			if (alilen < Opt.min_alilen) {
			    minlen = (minlen < Opt.min_alilen)
				? minlen : Opt.min_alilen;
			    dist = (dist * alilen + (Dist) Opt.missdist *
					(minlen - alilen))/minlen;
			}
		}
		setSeqReg(&ali1, (SeqPos)seldata.from1, (SeqPos)seldata.to1);
		setSeqReg(&ali2, (SeqPos)seldata.from2, (SeqPos)seldata.to2);

		addEdgeWithScore(SimG.edges, node1, node2, &ali1, &ali2,
			(Dist) dist, (Dist) score, (ConnCount) 1, (char)dir);

	}
	if (SimG.edges->edgenum == 0) {
		fprintf(stderr, "No data\n");
		exit(0);
	}
	if (Opt.VERBOSE) {
		fprintf(stderr, "Done (%d)\n", cnt);
	}
	if (SimG.edges->edgenum == 0) {
		fprintf(stderr, "No edge read (Is the homology file correct?)\n");
		exit(0);
	}
	setLeafNodeNum(SimG.nodes);

	if (Opt.VERBOSE) fprintf(stderr, "Sorting\n");

	createOrdList(&SimG);

	if (Opt.VERBOSE) fprintf(stderr, "Ceating Index\n");
	createIndex(&SimG);
	if (Opt.VERBOSE) fprintf(stderr, "Done\n");

	preproc_for_SpecInfo();
	setOutGroup(Opt.outgroupStr, Opt.outgroupRev);
	setSpMask(Opt.spmaskStr);
	domCluster(SimG.edges, SimG.nodes, &Opt);
	outputResults(&SimG, argc, argv);
	exit(0);
}

getargs(int argc, char **argv)
{
	int i, nn = 1;
	char *p;
	if (argc == 1 && isatty(0) == 1) {
		usage();
		fprintf(stderr, "  -h for more help\n");
		exit(0);
	}
	for (i = 1; i < argc; i++) {
		p = argv[i];
		if (*p == '-') {
			switch (*++p) {
			case 'a':
			/* parameters for merging adjacent clusters */
				p++;
				if (*p == 'i') {
					Opt.adjInclRatio = atof(++p);
				} else if (*p == 'o') {
					Opt.adjOvlpRatio = atof(++p);
				}
				break;
			case 'c':
			/* cutoff score or distance */
				if (isdigit(*++p)) {
					Opt.cutoff = atoi(p);
				} else if (*p == 'd') {
					Opt.distdiffcut = atof(++p);
				}
				break;
			case 'C':
			/* cutoff for splitting domains */
				if (*(++p)=='r') {
					Opt.cutoff2ratio = atof(++p);
				} else {
					Opt.cutoff2 = atoi(p);
				}
				break;
			case 's':
				Opt.sumcut = atoi(++p);
				break;
			case 'd':
				Opt.sim = 0;
				if (*(++p)) {
					Opt.cutoff = atoi(p);
				}
				break;
			case 'l':
				Opt.logdistscale = 1;
				break;
			case 'm':
			/* set missscore or missdist */
				if (*(++p)=='r') {
					Opt.missratio = atof(++p);
				} else {
					if (Opt.sim) {
						Opt.missscore = (Dist) atof(p);
					} else {
						Opt.missdist = (Dist) atof(p);
					}
				}
				break;
			case 'V':
			/* minimum alignment coverage for domain split */
				Opt.coverage = atof(++p);
				break;
			case 'v':
				Opt.VERBOSE = 1;
				int verbose_step1 = atoi(++p);
				if (verbose_step1) {
					Opt.verbose_step1 = verbose_step1;
				} else if (! Opt.verbose_step1) {
					Opt.verbose_step1 = 100000;
				}
				if (! Opt.verbose_step2) {
					Opt.verbose_step2 = 5000;
				}
				break;
			case 'D':
			/* debug flag */
				++p;
				if (*p=='\0') {
				} else if (isdigit(*p)) {
					Opt.DEBUG = Opt.DEBUG_val = atoi(p);
				} else {
					Opt.DEBUG_ent = p;
				}
				if (! Opt.DEBUG) {
					Opt.DEBUG = Opt.DEBUG_val = 1;
				}
				break;
			case 'o':
			/* output format */
				Opt.outstyle = atoi(++p);
				break;
			case 'n':
			/* minimum size of clusters to be output */
				p++;
				if (isdigit(*p)) {
					Opt.minsp = atoi(p);
					if (Opt.minent<Opt.minsp){
						Opt.minent = Opt.minsp;
					}
				} else if (*p == 's') {
					Opt.minsp = atoi(++p);
					if (Opt.minent<Opt.minsp){
						Opt.minent = Opt.minsp;
					}
				} else if (*p == 'e') {
					Opt.minent = atoi(++p);
					if (Opt.minent<Opt.minsp){
						Opt.minsp = Opt.minent;
					}
				} else if (*p == 'n') {
					Opt.delete_small = 0;
				}
				break;
			case 'S':
			/* use similarity measure */
				Opt.sim = 1;
				if (*(++p)) {
					Opt.cutoff = atoi(p);
				}
				break;
			case 'p':
			/* ratio for phylogenetic tree cutting */
				++p;
				if (*p == 'p') {
					Opt.phylocutratio2 = atof(++p);
				} else {
					Opt.phylocutratio = atof(p);
				}
				break;
			case 'H':
			/* homologous rather than orthologous clustering */
				Opt.phylocutratio = 2.0;
				break;
			case 'R':
			/* restore from dump file */
				restorefile = ++p;
				break;
			case 't':
			/* specifying related genomes with a tree file */
				if (*(p+1)) {
					sptreefile = ++p;
				} else if (*argv[i+1]!='-'){
					sptreefile = argv[++i];
				}
				break;
			case 'g':
				Opt.genes = ++p;
				break;
			case 'O':
			case '-':
			/* additional options */
			/* many of them are experimental */
				++p;
				if (strncmp(p, "minovlp=", 8)==0) {
					Opt.minovlp = atoi(p+8);
				} else if (strncmp(p, "ovlpratio=", 10)==0) {
					Opt.ovlpratio = atof(p+10);
				} else if (strncmp(p, "ovlpratio2=", 11)==0) {
					Opt.ovlpratio2 = atof(p+11);
				} else if (strncmp(p, "minlen=", 7)==0) {
					Opt.minlen = atoi(p+7);
				} else if (strncmp(p, "minlen2=", 8)==0) {
					Opt.minlen2 = atoi(p+8);
				} else if (strncmp(p, "min_alilen=", 11)==0) {
					Opt.min_alilen = atof(p+11);
				} else if (strncmp(p, "covfilt=", 8)==0) {
					Opt.covfilt = atof(p+8);
				} else if (strncmp(p, "distscale=",10)==0) {
					Opt.distscale = atoi(p+10);
				} else if (strncmp(p, "dpwin=", 6)==0) {
					Opt.dpWin = atoi(p+6);
				} else if (strncmp(p, "nbrScoreRatio=", 14)==0) {
					Opt.nbrScoreRatio = atof(p+14);
				} else if (strncmp(p, "nbrScoreDecay=", 14)==0) {
					Opt.nbrScoreDecay = atof(p+14);
				} else if (strncmp(p, "nbrScoreLim=", 12)==0) {
					Opt.nbrScoreLim = atof(p+12);
				} else if (strncmp(p, "iniGap=", 7)==0) {
					Opt.iniGap = (Dist)atof(p+7);
				} else if (strncmp(p, "extGap=", 7)==0) {
					Opt.extGap = (Dist)atof(p+7);
				} else if (strncmp(p, "skipGap=", 8)==0) {
					Opt.skipGap = (Dist)atof(p+8);
				} else if (strncmp(p, "chkInclClst=", 12)==0) {
					Opt.adjInclRatio = atof(p+12);
				} else if (strncmp(p, "adjInclRatio=", 13)==0) {
					Opt.adjInclRatio = atof(p+13);
				} else if (strncmp(p, "chkOvlpClst=", 12)==0) {
					Opt.adjOvlpRatio = atof(p+12);
				} else if (strncmp(p, "adjOvlpRatio=", 13)==0) {
					Opt.adjOvlpRatio = atof(p+13);
				} else if (strncmp(p, "nochkInclClst", 13)==0) {
					Opt.adjInclRatio = 0;
				} else if (strncmp(p, "chkConnect", 10)==0) {
					if (p[10] == '=') {
						Opt.chkConnect = atof(p+11);
					} else {
						Opt.chkConnect = 2.0;
					}
				} else if (strncmp(p, "phyloMode=", 10)==0) {
					Opt.phyloMode = atoi(p+10);
				} else if (strncmp(p, "treecheck", 9)==0) {
					Opt.treecheck = 1;
				} else if (strncmp(p, "skipErrEnt", 10)==0) {
					Opt.skipErrEnt = 1;
				} else if (strncmp(p, "revMatch", 8)==0) {
					Opt.revMatch = 1;
				} else if (strncmp(p, "domBoundary", 11)==0) {
					Opt.domBoundary = B_ALIREG;
				} else if (strncmp(p, "mincutcnt=", 10)==0) {
					Opt.mincutcnt = atoi(p+10);
				} else if (strncmp(p, "domcut", 6)==0) {
					Opt.domcut = 1;
				} else if (strncmp(p, "nobreak", 7)==0) {
					Opt.nobreak = 1;
				} else if (strncmp(p, "nbrConnRatio=", 13)==0) {
					Opt.nbrConnRatio = atof(p+13);
				} else if (strncmp(p, "outgroup=", 9)==0) {
					Opt.outgroupStr = p+9;
					Opt.outgroupRev = 0;
					if (! Opt.outgroupMode)
						Opt.outgroupMode = 1;
				} else if (strncmp(p, "ingroup=", 8)==0) {
					Opt.outgroupStr = p+8;
					Opt.outgroupRev = 1;
					if (! Opt.outgroupMode)
						Opt.outgroupMode = 1;
				} else if (strncmp(p, "spmask=", 7)==0) {
					Opt.spmaskStr = p+7;
				} else if (strncmp(p, "outgroupMode=",13)==0){
					Opt.outgroupMode = atoi(p+13);
				} else if (strncmp(p, "noextend", 8)==0) {
					Opt.noextend = 1;
				} else if (strncmp(p, "verbose_step1=", 14)==0) {
					Opt.verbose_step1 = atoi(p+14);
				} else if (strncmp(p, "verbose_step2=", 14)==0) {
					Opt.verbose_step2 = atoi(p+14);
				} else if (strncmp(p, "NEW", 3)==0) {
					/** version flag: temporal option **/
					Opt.NEW = 1;
				}
				break;
			case 'h':
				help();
				exit(0);
			case '\0':
				if (nn == 1) {
					strcpy(filename, "stdin");
				} else {
					strcpy(genefilename, "stdin");
				}
				break;
			default:
				break;
			}
		} else {
			switch (nn) {
			case 1:
				strcpy(filename, p);
				break;
			case 2:
				strcpy(genefilename, p);
				break;
			}
			nn++;
		}
	}
	if (! Opt.cutoff2) {
		if (Opt.cutoff2ratio) {
			Opt.cutoff2 = Opt.cutoff * Opt.cutoff2ratio;
		} else if (Opt.sim) {
			Opt.cutoff2 = Opt.cutoff;
		}
	}
}
usage()
{
	fprintf(stderr, "DomClust ver.%s\n", DOMCLUST_VERSION);
	fprintf(stderr, "Usage: domclust [options] homfile genetab\n");
}
help()
{
	usage();
	fprintf(stderr, "\
    -S     use similarity as a measure of relatedness [on]\n\
    -d     use distance (or disimilarity) as a measure of relatedness\n\
    -c#    cutoff score/distance (can also be spcified as -S# or -d#) [60]\n\
    -m#    score/distance for missing relationships (m<c)\n\
    -mr#   specify a missing score as a ratio to c (0<mr<1) [0.95]\n\
    -C#    cutoff score for domain split (c<=C)\n\
    -V#    alignment coverage for domain split (0<=V<=1)\n\
    -n#    minimum # of organisms in clusters to be output [2]\n\
    -ne#   minimum # of entries in clusters to be output [2]\n\
    -p#    ratio of phylogenetic pattern overlap for tree cutting [0.5]\n\
    -H     homology clustering (i.e. skip the tree cutting)\n\
    -ai#   member overlap for absorbing adjacent small clusters (0<=ai<=1)\n\
    -ao#   member overlap for merging adjacent clusters (0<=ao<=ai)\n\
    -t<fn> use a tree file for weighting related genomes\n\
    -R<fn> restore from dump file\n\
    -o#    output format\n\
");
}
