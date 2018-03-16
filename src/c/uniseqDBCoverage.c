#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <ctype.h>

#define MAX_ID 1000
#define MAX_BUFFER 50000
#define PNG_DIM 1024

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

// TBD

struct TaxoItem{
	// TaxoID
	int64_t taxoID;
	// Taxo Name
	char taxoName[MAX_ID];
	// Reference Length to UniSeq
	uint64_t refLength;
	// Reference to Next
	struct TaxoItem * next;
};

struct TaxoFile{
	struct TaxoItem** items;
	uint64_t num_items;
};

/*******************************/

//Ascii to UInt64
uint64_t asciiToUint64(const char *text);
// Print the error message 's' and exit(-1)
void terror(const char *s);
// Retrieve Output of CMD (Depends on MAX_BUFFER)
char * getStdoutFromCommand(char* cmd);
// Get lines in FILE
uint64_t getLinesInFile(char* filename);
// Check if character valid nucleotide annotation
int validc(char c);
// Get amount of nucleotide in FASTA
int64_t getLengthFasta(char* file);
// Get amount of nucleotide per sequence in FASTA
int64_t* getLengthFastaSequences(char* file);
// Parse Line to ParsedBlastInfo structure
int parsePBI(char * buffer, struct ParsedBlastInfo* pbi);
// Parse SeqInfo to PBISeqInfo
int parseSeqInfo(char* seqInfo, struct PBISeqInfo* psi);
// Print ParsedBlastInfo structure
void printParsedBlastInfo(struct ParsedBlastInfo* pbi);
// Print PBISeqInfo structure
void printPBISeqInfo(struct PBISeqInfo* psi);
// Create TaxoFile
void makeTaxoFile(char * filename, struct TaxoFile * tf);
// Initialize TaxoFile
void initTaxoFile(struct TaxoFile* tf, uint64_t size);
// Hash Sequence ID
uint64_t hashSequenceID(char* genomeID, uint64_t hashsize);
// Retrieve UniSeq Reference Length of the Genome ID in PBISeqInfo
int64_t retrieveRefLength(struct TaxoFile*, struct PBISeqInfo);

/*******************************/

// gcc -g -D_FILE_OFFSET_BITS=64 uniseqDBCoverage.c -o ../bin/uniseqDBCoverage
int main(int argc, char** argv){
	FILE * fi, *fi_c, * fo;
	struct ParsedBlastInfo pbi;
	struct PBISeqInfo psi;
	struct TaxoFile *tf;
	int64_t DB_LENGTH;
	char * DB_seq, * DB_seq_c;
	char buffer[MAX_BUFFER];
	int64_t scan, temp;
	int64_t current_reflength, current_start, current_end, coverage_counter = 0, coverage_counter_c = 0;
	double coverage;

	int i = 0;
	int64_t ii = 0;
	int DEBUG = 1;

	printf("Starting...\n");

	// ./uniseqDBCoverage /home/pablorod/data/metagenomics/REVCO/parsed_format_gastro_reads_fix /home/pablorod/data/metagenomics/REVCO/gastrointestinal_tract.fasta.taxo /home/pablorod/data/metagenomics/REVCO/UNIDB_gastrointestinal_tract.fasta /home/pablorod/data/metagenomics/REVCO/CoverageREADS
	if (argc < 8)//				0					1							2						3		4			5		6			7
			terror("USE: ./uniseqDBCoverage <Fixed_Parsed_Blast_Reads> <Fixed_Parsed_Blast_Contigs> <DB.taxo> <uniseqDB> <Output> <n_reads> <n_contigs>");
	// Load <Fixed_Parsed_Blast>
	if((fi=fopen(argv[1],"rt"))==NULL){
		terror("Can't open INPUT file");
	}
	if((fi_c=fopen(argv[2],"rt"))==NULL){
		terror("Can't open INPUT file");
	}
	// Load <DB.taxo> + Make TaxoFile
	if((tf = (struct TaxoFile*) calloc(1, sizeof(struct TaxoFile))) == NULL)
		terror("No memory for TaxoFile");

	printf("Making TaxoFile...\n");
	makeTaxoFile(argv[3], tf);

	//DEBUG printf("D1: TF Size: %"PRIu64"\n", tf->num_items);
	
	// Load <uniseqDB> + Obtain DB length¡
	printf("Loading UniSeq...\n");
	DB_LENGTH = getLengthFasta(argv[4]);

	if ((DB_seq = (char*) calloc(DB_LENGTH, sizeof(char))) == NULL)
		terror("No memory for Seq");

	if ((DB_seq_c = (char*) calloc(DB_LENGTH, sizeof(char))) == NULL)
		terror("No memory for Seq_c");

	//DEBUG printf("D2: FASTA Length: %"PRId64"\n", DB_LENGTH);

	char ** error_check;
	double reads_match_ratio = 0.0, contigs_match_ratio = 0.0;
	uint64_t reads_match_counter = 0, contigs_match_counter = 0;

	// Read Fixed_Parsed_Blast
	printf("Reading Fixed_Parsed_Blast_Reads...\n");
	while(!feof(fi)){
		// Scan line by line and load to PBI structure
		fgets(buffer, MAX_BUFFER, fi);
		scan = parsePBI(buffer, &pbi);
		if(scan != 12) terror("PBI structure loaded incorrectly");

		// Parse seqInfo to PSI structure
		scan = parseSeqInfo(pbi.seqInfo, &psi);
		if(scan != 4) terror("PSI structure loaded incorrectly");

		reads_match_counter++;
		// Obtain RefLength of GenomeID from Taxo file
		current_reflength = retrieveRefLength(tf, psi);
		//printf("current RefLength: %"PRId64"\n", current_reflength);

		if(current_reflength==-1){
			printf("--------------NOT FOUND--------------\n");
		}else{
			// Make changes to DB_seq
			current_start = current_reflength + pbi.g_start;
			current_end = current_reflength + pbi.g_end;

			if(current_end < current_start){
				temp = current_start;
				current_start = current_end;
				current_end = temp;
			}
			for(ii = current_start; ii<current_end; ii++){
				DB_seq[ii]='1';
			}
		}
	}

	reads_match_ratio = (double)reads_match_counter/strtod(argv[6], error_check);

	printf("Reading Fixed_Parsed_Blast_Contigs...\n");
	while(!feof(fi_c)){
		// Scan line by line and load to PBI structure
		fgets(buffer, MAX_BUFFER, fi_c);
		scan = parsePBI(buffer, &pbi);
		if(scan != 12) terror("PBI structure loaded incorrectly");

		// Parse seqInfo to PSI structure
		scan = parseSeqInfo(pbi.seqInfo, &psi);
		if(scan != 4) terror("PSI structure loaded incorrectly");

		contigs_match_counter++;

		// Obtain RefLength of GenomeID from Taxo file
		current_reflength = retrieveRefLength(tf, psi);
		//printf("current RefLength: %"PRId64"\n", current_reflength);

		if(current_reflength==-1){
			printf("--------------NOT FOUND--------------\n");
		}else{
			// Make changes to DB_seq
			current_start = current_reflength + pbi.g_start;
			current_end = current_reflength + pbi.g_end;

			if(current_end < current_start){
				temp = current_start;
				current_start = current_end;
				current_end = temp;
			}
			for(ii = current_start; ii<current_end; ii++){
				DB_seq_c[ii]='1';
			}
		}
	}

	contigs_match_ratio = (double)contigs_match_counter/strtod(argv[7], error_check);

	// Load <Output>
	if((fo=fopen(argv[5],"wt"))==NULL){
		terror("Can't open OUTPUT file");
	}
	FILE * fseq, *fseq_c, * finfo, * fratio;


	char name[250];
	strcpy(name, argv[5]);

	if((fseq=fopen(strcat(name,"READS.seq"),"wt"))==NULL){
		terror("Can't open FSEQ file");
	}
	strcpy(name, argv[5]);
	if((fseq_c=fopen(strcat(name,"CONTIGS.seq"),"wt"))==NULL){
		terror("Can't open FSEQ file");
	}
	strcpy(name, argv[5]);
	if((finfo=fopen(strcat(name,".info"),"wt"))==NULL){
		terror("Can't open INFO file");
	}
	strcpy(name, argv[5]);
	if((fratio=fopen(strcat(name,".ratio"),"wt"))==NULL){
		terror("Can't open RATIO file");
	}

	fprintf(fratio, "SequenceSet;TotalMatches;NumSequences;MatchesPerSequence\n");
	fprintf(fratio, "Reads;%"PRIu64";%"PRIu64";%lf\n", reads_match_counter, asciiToUint64(argv[6]), reads_match_ratio);
	fprintf(fratio, "Contigs;%"PRIu64";%"PRIu64";%lf\n", contigs_match_counter, asciiToUint64(argv[7]), contigs_match_ratio);

	uint64_t max_hits_per_pixel = (uint64_t)(DB_LENGTH / (PNG_DIM * PNG_DIM)) + 1;
	uint64_t pixel_counter = 0, hits_at_pixel = 0, hits_at_pixel_c = 0;
	uint64_t nucleotide_counter = 0;

	uint64_t read_hits = 0, contig_hits = 0, both_yes = 0, both_no = 0;

	// ---
	fprintf(fseq, "Pixel;Hits<%"PRIu64"\n", max_hits_per_pixel);
	fprintf(fseq_c, "Pixel;Hits<%"PRIu64"\n", max_hits_per_pixel);

	// Calculate DB Coverage
	printf("Calculating coverage...\n");

	for(ii = 0; ii<DB_LENGTH; ii++){
		if(DB_seq[ii]=='1'){
			coverage_counter++;
			hits_at_pixel++;
			read_hits++;
		}
		if(DB_seq_c[ii]=='1'){
			coverage_counter_c++;
			hits_at_pixel_c++;
			contig_hits++;
		}
		if(DB_seq[ii] == DB_seq_c[ii]){
			if(DB_seq[ii] == '1')
				both_yes++;
		}


		if (nucleotide_counter == max_hits_per_pixel){
			fprintf(fseq, "%"PRIu64";%"PRIu64"\n", pixel_counter, hits_at_pixel);
			fprintf(fseq_c, "%"PRIu64";%"PRIu64"\n", pixel_counter, hits_at_pixel_c);
			hits_at_pixel = 0;
			hits_at_pixel_c = 0;
			nucleotide_counter = 0;
			pixel_counter++;
		}
		nucleotide_counter++;
	}
	coverage = 100 * (double)coverage_counter/(double)DB_LENGTH;
	both_no = DB_LENGTH - read_hits - contig_hits - both_yes;
	// ---

	fprintf(finfo, "ReadHits\tContigHits\tBothYes\tBothNo\n");
	fprintf(finfo, "%"PRIu64"\t%"PRIu64"\t%"PRIu64"\t%"PRIu64"\n", read_hits, contig_hits, both_yes, both_no);
	fclose(fseq);
	fclose(fseq_c);
	fclose(finfo);

	printf("UniSeq Length\tCovered Bases\tCoverage Percentage\nREADS\t%"PRId64"\t%lf\n",
		coverage_counter, coverage);
	fprintf(fo, "UniSeq Length\tCovered Bases\tCoverage Percentage\nREADS\t%"PRId64"\t%lf\n",
		coverage_counter, coverage);

	coverage = 100 * (double)coverage_counter_c/(double)DB_LENGTH;
	printf("CONTIGS\t%"PRId64"\t%lf\n", 
		coverage_counter_c, coverage);
	fprintf(fo, "CONTIGS\t%"PRId64"\t%lf\n", 
		coverage_counter_c, coverage);

	// Free Memory
	free(DB_seq); free(DB_seq_c);

	struct TaxoItem * temp_ti, * holder;
	for(i=0; i<tf->num_items; i++){
		temp_ti = tf->items[i];
		while(temp_ti->next != NULL){
			holder = temp_ti->next->next;
			free(temp_ti->next);
			temp_ti->next = holder;
		}
		free(temp_ti);
	}
	free(tf->items);
	free(tf);

	// Close INPUT and OUTPUT files
	printf("UniseqDBCoverage succesfully executed...\n");
	fclose(fi);
	fclose(fo);
}

/*******************************/

//Print the error message 's' and exit(-1)
void terror(const char *s) {
	printf("ERROR::**** %s ****\n", s);
	exit(-1);
}
// 
uint64_t asciiToUint64(const char *text){
	uint64_t number=0;

	for(;*text;text++){
		char digit=*text-'0';           
		number=(number*10)+digit;
	}
	return number;
}

// Get lines in FILE
uint64_t getLinesInFile(char * filename){
	FILE * f = fopen(filename, "rt");
	if(f==NULL) terror("Could not open file");
	uint64_t lines = 0;
	while(!feof(f)){
		if(getc(f) == '\n') lines++;
	}
	fclose(f);
	return lines;
}

// Check if character valid nucleotide annotation
int validc(char c){
	//if( c=='A' || c== 'C' || c== 'G' || c== 'T')return 1;
	if( c=='A' || c== 'C' || c== 'G' || c== 'T' || c=='N' || c=='U' )return 1;
	//fprintf(stderr,"Warning: -%c- in seq",c);
	return 0;
}

// Get amount of nucleotide in FASTA
int64_t getLengthFasta(char* file){
	int64_t *V;
	int64_t length=0;
	int i=0;

	V=getLengthFastaSequences(file);
	for(i=1;i<V[0];i++)length += V[i];

	return length;
}

// Get amount of nucleotide per sequence in FASTA
int64_t* getLengthFastaSequences(char* file){
	int DEBUG=0;

	FILE* fe;
	if((fe=fopen(file,"r"))==NULL){
		printf("***ERROR Opening seq input file %s\n",file);
		exit(-1);
	}
	int64_t n=0;
	char c;

	// Creates an array of 50 cromosomas
	int64_t nSeq[5000];
	int64_t *lengthsArray;
	int i;

	int nSequence=1;
	// Skip first line and copy
	c=getc(fe);

	// Loop for each sequence in fasta. Skip ID then get length of current sequence
	while(!feof(fe)){

		// Skip first line (ID)
		while(c!='\n') {
			c=getc(fe);
			c=toupper(c);
			if(DEBUG)printf("%c",c);
		}

		// Get length of current sequence
		while(!feof(fe) && c!= '>') {
			c=getc(fe);
			c=toupper(c);
			if(validc(c)){
				n++;
			}
		}

		// DEBUG
		if(DEBUG)printf("\nSequence: %d\t Longitud: %d\n",nSequence,n);


		nSeq[nSequence++]=n;
		n=0;
	}

	lengthsArray = (int64_t*)malloc(sizeof(int64_t)*nSequence);
	lengthsArray[0]=nSequence;
	for(i=1;i<nSequence;i++){
		lengthsArray[i]=nSeq[i];
	}

	fclose(fe);
	return lengthsArray;
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
	printf("ReadID: %s | GenomeID: %s | GenLenth: %d | K Blast: %d\n",
		psi->readID, psi->genomeID, psi->genLength, psi->k_BLAST);
}

// Hash Sequence ID
uint64_t hashSequenceID(char* genomeID, uint64_t hashsize){
	char *c = genomeID;
	uint64_t hash = 0;

	while(*c != '\0'){
		hash = (uint64_t)(*c + 31*hash);
		c++;
	}

	return hash % hashsize;
}

// Retrieve UniSeq Reference Length of the Genome ID in PBISeqInfo
int64_t retrieveRefLength(struct TaxoFile * tf, struct PBISeqInfo psi){
	int64_t reference = -1, current_hash;
	struct TaxoItem * ti;
	int FOUND = 0;

	current_hash = hashSequenceID(psi.genomeID, tf->num_items);

	ti = tf->items[current_hash];
	while(FOUND == 0 && ti != NULL){
		/*printf("G1: |%s|%d|\n", psi.genomeID, hashSequenceID(psi.genomeID, tf->num_items));
		printf("G2: |%s|%d|\n", ti->taxoName, hashSequenceID(ti->taxoName, tf->num_items));*/

		if(strcmp(psi.genomeID, ti->taxoName)==0){
			FOUND = 1;
			reference = ti->refLength;
		}
		ti = ti->next;
	}
	return reference;
}

// Initialize TaxoFile
void initTaxoFile(struct TaxoFile *tf, uint64_t size){
	int i;

	for(i = 0; i < size; i++){
		tf->items[i]->taxoID = -1;
		//tf.items[i]->refLength=-1;
		//tf.items[i]->taxoName==NULL;
		//tf->items[i]->next==NULL;
		//DEBUG printf("base taxoID: %"PRId64"\n", tf.items[i]->taxoID);
	}
	tf->num_items = size;
}

// Create TaxoFile
void makeTaxoFile(char* filename, struct TaxoFile * tf){
	struct TaxoItem* current_item;
	struct TaxoItem new_ti;
	uint64_t temp_hash, tf_size;
	FILE * fi;
	char buffer[MAX_BUFFER], c;
	int scan, tempi;
	
	tf_size = getLinesInFile(filename);

	//printf("DEBUG: lines - %"PRIu64"\n", tf_size);
	/*if( (current_item = (struct TaxoItem*) calloc(1, sizeof(struct TaxoItem)) )==NULL)
        	terror("No memory for TaxoFile::items");
*/
	// Initialize new TaxoFile
	if ( (tf->items = (struct TaxoItem**) calloc(tf_size, sizeof(struct TaxoItem*)) )==NULL ){
        terror("No memory for TaxoFile");
    }

    for(tempi=0; tempi<tf_size;tempi++){
        if( (tf->items[tempi] = (struct TaxoItem*) calloc(1, sizeof(struct TaxoItem)) )==NULL)
        	terror("No memory for TaxoFile::items");
    }

	if((fi=fopen(filename,"r"))==NULL){
		terror("***ERROR Opening Taxo file");
	}

	initTaxoFile(tf, tf_size);
    // Read line by line and create new TaxoItem
    int counter = 0;
	while(!feof(fi)){

		fgets(buffer, MAX_BUFFER, fi);

		// 0       36381   0       ADMP01000001    Bacteroides xylanisolvens SD CC 2a contig00137, whole genome shotgun sequence. [Bacteroides xylanisolvens SD CC 2a]
		scan=sscanf(buffer,"%"PRId64"\t%*"PRIu64"\t%*"PRIu64"\t%"PRIu64"\t%s\t%*s\n",
				&new_ti.taxoID, &new_ti.refLength, new_ti.taxoName);

		//printf(buffer);
		//printf("scan : %d | TaxoID: %"PRId64" | Reference Length: %"PRIu64" | TaxoName: %s\n", scan, new_ti.taxoID, new_ti.refLength, new_ti.taxoName);

		// If read properly, then add TaxoItem to TaxoFile
		if(scan==3){
			temp_hash = hashSequenceID(new_ti.taxoName, tf_size);
						
			if(tf->items[temp_hash]->taxoID == -1){
				tf->items[temp_hash]->taxoID = new_ti.taxoID;
				tf->items[temp_hash]->refLength = new_ti.refLength;
				strcpy(tf->items[temp_hash]->taxoName, new_ti.taxoName);
				tf->items[temp_hash]->next = NULL;
			}
			else{
				current_item = tf->items[temp_hash];

				struct TaxoItem * aux_ti;
				if( (aux_ti = (struct TaxoItem*) calloc(1, sizeof(struct TaxoItem)) )==NULL)
        			terror("No memory for TaxoFile::items");

				aux_ti->taxoID = new_ti.taxoID;
				aux_ti->refLength = new_ti.refLength;
				strcpy(aux_ti->taxoName, new_ti.taxoName);
				aux_ti->next = current_item;

				tf->items[temp_hash] = aux_ti;
			}
			counter++;
		}
	}
	//printf("Counter: %d\n", counter);
}