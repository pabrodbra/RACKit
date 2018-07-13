#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <ctype.h>

#define MAX_ID 1000
#define MAX_BUFFER 50000

struct ParsedBlastInfo{
	// Sequences Info
	char seqInfo[MAX_ID];
	// BLAST score or E-Value
	float score;
	// Identity
	uint64_t identity;
	// Length
	uint64_t length;
	// Similarity
	uint64_t similarity;
	// Initial Gaps
	uint64_t i_gaps;
	// Extended Gaps
	uint64_t e_gaps;
	// Strand
	char strand[2];
	// Read Start
	uint64_t r_start;
	// Read End
	uint64_t r_end;
	// Genome Start
	uint64_t g_start;
	// Genome End
	uint64_t g_end;
};

struct PBISeqInfo{
	// Read ID
	char readID[MAX_ID];
	// Genome ID
	char genomeID[MAX_ID];
	// Genome Length
	uint64_t genLength;
	// BLAST k Result
	uint64_t k_BLAST;
};

/*******************************/

// Print the error message 's' and exit(-1)
void terror(const char *s);
// Parse Line to ParsedBlastInfo structure
int parsePBI(char * buffer, struct ParsedBlastInfo* pbi);
// Parse SeqInfo to PBISeqInfo
int parseSeqInfo(char* seqInfo, struct PBISeqInfo* psi);
// Print ParsedBlastInfo structure
void printParsedBlastInfo(struct ParsedBlastInfo* pbi);
// Print PBISeqInfo structure
void printPBISeqInfo(struct PBISeqInfo* psi);

/*******************************/

// gcc -g -D_FILE_OFFSET_BITS=64 filterParsedBlast.c -o ../bin/filterParsedBlast
int main(int argc, char** argv){
	FILE * fi, * fo;
	char buffer[MAX_BUFFER];
	struct ParsedBlastInfo pbi, pbi_current_best;
	struct PBISeqInfo psi, psi_current_best;
	int scan;

	// ./filterParsedBlast /home/pablorod/data/metagenomics/REVCO/parsed_format_gastro_reads_fix /home/pablorod/data/metagenomics/REVCO/FILTER_pb_format_gastro_reads 0
	if (argc < 4)//				0					1							2						3	
			terror("USE: ./filterParsedBlast <Fixed_Parsed_Blast> <OUTPUT_Fixed_Parsed_Blast> <FilterType(0=MultipleHits;1=UniqueHit>");
	
	if((fi=fopen(argv[1],"rt"))==NULL){
		terror("Can't open INPUT file");
	}

	if((fo=fopen(argv[2],"wt"))==NULL){
		terror("Can't open OUTPUT file");
	}

	char filterFlag = argv[3][0];
	char pbi_loaded = '0';

	while(!feof(fi)){
		// Scan line by line and load to PBI structure
		fgets(buffer, MAX_BUFFER, fi);
		scan = parsePBI(buffer, &pbi);
		if(scan != 12) terror("PBI structure loaded incorrectly");

		// Parse seqInfo to PSI structure
		scan = parseSeqInfo(pbi.seqInfo, &psi);
		if(scan != 4) terror("PSI structure loaded incorrectly");

		// Check if it's a new read
		if(strcmp(psi.readID,psi_current_best.readID)==0){
			pbi_loaded = '1'; // 1 => Different pbi for same read
		}else{
			pbi_loaded = '0'; // 0 => Different read, new pbi
			psi_current_best = psi;
			pbi_current_best = pbi;
		}

		// Apply filter
		if(pbi_loaded == '0'){
			// Unique Hits + Multiple Hits | Print first/top PBI
			fprintf(fo, "%s;%f;%"PRIu64";%"PRIu64";%"PRIu64";%"PRIu64";%"PRIu64";%s;%"PRIu64";%"PRIu64";%"PRIu64";%"PRIu64"\n",
		pbi.seqInfo, pbi.score, pbi.identity, pbi.length, pbi.similarity, pbi.i_gaps, pbi.e_gaps, pbi.strand, pbi.r_start, pbi.r_end, pbi.g_start, pbi.g_end);
		}

		if(filterFlag == '0'){
			// Multiple Hits
			if(pbi_loaded == '1' && pbi.identity == pbi_current_best.identity && pbi.length == pbi_current_best.length && pbi.score == pbi_current_best.score){
				fprintf(fo, "%s;%f;%"PRIu64";%"PRIu64";%"PRIu64";%"PRIu64";%"PRIu64";%s;%"PRIu64";%"PRIu64";%"PRIu64";%"PRIu64"\n",
		pbi.seqInfo, pbi.score, pbi.identity, pbi.length, pbi.similarity, pbi.i_gaps, pbi.e_gaps, pbi.strand, pbi.r_start, pbi.r_end, pbi.g_start, pbi.g_end);
			}
		}
	}

	fclose(fi);
	fclose(fo);
	printf("filterParsedBlast succesfully executed...\n");
	return(0);
}

/*******************************/

//Print the error message 's' and exit(-1)
void terror(const char *s) {
	printf("ERROR::**** %s ****\n", s);
	exit(-1);
}

// Parse Line to ParsedBlastInfo structure
int parsePBI(char * buffer, struct ParsedBlastInfo* pbi){
	int scan;

	scan=sscanf(buffer,"%[^;];%f;%"PRIu64";%"PRIu64";%"PRIu64";%"PRIu64";%"PRIu64";%[^;];%"PRIu64";%"PRIu64";%"PRIu64";%"PRIu64"\n",
				&pbi->seqInfo, &pbi->score, &pbi->identity, &pbi->length, &pbi->similarity, &pbi->i_gaps, &pbi->e_gaps, &pbi->strand, &pbi->r_start, &pbi->r_end, &pbi->g_start, &pbi->g_end);

	//DEBUG | printParsedBlastInfo(pbi);

	return scan;
}

// Parse SeqInfo to PBISeqInfo
int parseSeqInfo(char* seqInfo, struct PBISeqInfo* psi){
	int scan;

	scan = sscanf(seqInfo, ">%[^>]>%[^ ] %"PRIu64" %"PRIu64,
		&psi->readID, &psi->genomeID, &psi->genLength, &psi->k_BLAST);

	//DEBUG | printPBISeqInfo(psi);
	return scan;
}
// Print ParsedBlastInfo structure
void printParsedBlastInfo(struct ParsedBlastInfo* pbi){
	printf("%s\n%f\n%"PRIu64"\n%"PRIu64"\n%"PRIu64"\n%"PRIu64"\n%"PRIu64"\n%s\n%"PRIu64"\n%"PRIu64"\n%"PRIu64"\n%"PRIu64"\n",
		pbi->seqInfo, pbi->score, pbi->identity, pbi->length, pbi->similarity, pbi->i_gaps, pbi->e_gaps, pbi->strand, pbi->r_start, pbi->r_end, pbi->g_start, pbi->g_end);
}

// Print PBISeqInfo structure
void printPBISeqInfo(struct PBISeqInfo* psi){
	printf("ReadID: %s | GenomeID: %s | GenLenth: %"PRIu64" | K Blast: %"PRIu64"\n",
		psi->readID, psi->genomeID, psi->genLength, psi->k_BLAST);
}