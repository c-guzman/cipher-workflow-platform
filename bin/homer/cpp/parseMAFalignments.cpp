
// Copyright 2009 - 2015 Christopher Benner <cbenner@salk.edu>
// 
// This file is part of HOMER
//
// HOMER is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// HOMER is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

#include "SeqTag.h"

#define MAF_IO_BUFFER 1000000
#define MAX_SPECIES 200
#define WHITE_SPACE 7

void printCMD() {
	fprintf(stderr, "\n\tparseMAFalignments <peak file> <MAF directory> [options]\n");
	fprintf(stderr, "\n\tOptions:\n");
	fprintf(stderr, "\t\t-ref <genome> (reference genome, def: hg19)\n");
	fprintf(stderr, "\t\t-cmp <genome> (genome to compare, def: panTro4)\n");
	fprintf(stderr, "\n");
	exit(0);
}

void split(char* string, char** cols, int &numCols, char delim) {
	cols[0] = string;
	numCols=1;
	char delim2 = 0;
	if (delim == WHITE_SPACE) {
		delim = '\t';
		delim2 = 32;
	}
	int len = strlen(string);
	for (int i=0;i<len;i++) {
		if (string[i] == delim || string[i] == delim2) {
			string[i] = '\0';
			while (i<len && (string[i] == delim || string[i] == delim2)) {
				i++;
			}
			if (i>=len) break;
			cols[numCols] = &(string[i+1]);
			numCols++;
		} else if (string[i] == '\n') {
			string[i] = '\0';
		} else if (string[i] == '\r') {
			string[i] = '\0';
		}
	}
}


int main(int argc, char** argv) {


	char* mafDirectory = NULL;
	char* peakFile = NULL;
	char* refGenome = "hg19";
	char* cmpGenome = "panTro4";
	char* outputFileName = NULL;

	if (argc < 3) printCMD();
	peakFile = argv[1];
	mafDirectory = argv[2];

	for (int i=3;i<argc;i++) {
		if (strcmp(argv[i],"-ref")==0) {
			refGenome = argv[++i];
		} else if (strcmp(argv[i],"-cmp")==0) {
			cmpGenome = argv[++i];
		} else {
			printCMDdenovo();
		}
	}

	fprintf(stderr, "\n");
	fprintf(stderr, "\tAnalyzing peak file: %s\n", peakFile);
	fprintf(stderr, "\tLooking for gzipped maf files in directory: %s\n", mafDirectory);
	fprintf(stderr, "\tReference Genome: %s\n", targetGenome);
	fprintf(stderr, "\tComparison Genome: %s\n", cmpGenome);
	fprintf(stderr, "\n");


	PeakLibrary* peaks = new PeakLibrary(peakfile,PEAK_READ_MODE_NORMAL);

	FILE* outputFp = stdout;
	if (outputFileName != NULL) {
		outputFp = fopen(outputFileName, "w");
	}

	fprintf(outputFp, "Peak\tchr\tstart\tend\tstrand\tsize\talignSize\talignMatch\talignMisMatch\tmisR\tmisC\tmisNA\n");

	char* openStr = new char[100000];
	char* buf = new char[MAF_IO_BUFFER];
	char** cols = new char[1000];
	int numCols=0;
	char** cols2 = new char[1000];
	int numCols2=0;


	char** aSeq = new char*[MAX_SPECIES];	
	char** aOrg = new char*[MAX_SPECIES];	
	char** aChr = new char*[MAX_SPECIES];	
	int* aPos = new int[MAX_SPECIES];	
	int* aLen = new int[MAX_SPECIES];	
	int* aCLen = new int[MAX_SPECIES];	
	int* aSLen = new int[MAX_SPECIES];	
	char* aStrand = new char[MAX_SPECIES];	
	for (int i=0;i<MAX_SPECIES;i++) {
		aSeq[i] = new char[MAF_IO_BUFFER];
		aOrg[i] = new char[10000];
		aChr[i] = new char[10000];
	}

	char** chr = peaks->chrs->keys();
	qsort(chr,peaks->chrs->total,sizeof(char*),&chrcmp);
	for (int i=0;i<peaks->chrs->total;i++) {
		ChrPeaks* cp = (ChrPeaks*)peaks->chrs->search(chr[i]);

		sprintf(openStr, "zcat %s/%s.maf.gz", mafDirectory, chr[i]);
		FILE* fp = popen(openStr,"r");
		if (fp == NULL) {
			fprintf(stderr, "\t! Could not open/find file for %s: %s\n", chr[i], openStr);
			continue;
		}

		fprintf(stderr, "\tAnalyzing %s\n", chr[i]);

		int* asize = new int[cp->numPeaks];
		int* amatch = new int[cp->numPeaks];
		int* amis = new int[cp->numPeaks];
		int* amisR = new int[cp->numPeaks];
		int* amisT = new int[cp->numPeaks];
		int* amisNA = new int[cp->numPeaks];
		for (int j=0;j<cp->numPeaks;j++) {
			asize[j]=0;
			amatch[j]=0;
			amis[j]=0;
			amisR[j]=0;
			amisT[j]=0;
			amisNA[j]=0;
		}

		int peakIndex = 0;
		int numSpecies = 0;

		while (fgets(buf, BUFFER, inputfp) != NULL) {
			if (buf[0] == 'a') {
				// new alignment
				peakIndex = analyzeAlignment()xxxxx;

				numSpecies = 0;
			} else if (buf[0] == 'i') {
				continue;
			} else if (buf[0] == 'e') {
				continue;
			} else if (buf[0] == 's') {
				split(buf,cols,numCols,WHITE_SPACE);
				if (numCols < 7) {
					fprintf(stderr, "!!! Something might be wrong, a line has %d columns!\n", numCols);
				}
				split(cols[1],cols2,numCols2,WHITE_SPACE);
				int position = 0;
				int length = 0;
				sscanf(cols[2],"%d",&position);
				sscanf(cols[3],"%d",&length);
	
				strcpy(aSeq[numSpecies],cols[6]);
				strcpy(aOrg[numSpecies],cols2[0]);
				strcpy(aChr[numSpecies],cols2[1]);
				aPos[numSpecies]=position;
				aLen[numSpecies]=length;
				aSLen[numSpecies]=strlen(cols[6]);
				aStrand[numSpecies]=cols[4][0];
				
				numSpecies++;
			}
		}
		pclose(fp);

		if (peakIndex < cp->numPeaks && numSpecies > 0) {
			peakIndex = analyzeAlignment();
		}
	
		delete []asize;
		delete []amatch;
		delete []amis;
		delete []amisR;
		delete []amisT;
		delete []amisNA;
		delete [](chr[i]);
	}
	delete []chr;

	return 0;
}

int analyzeAlignment() {
}
