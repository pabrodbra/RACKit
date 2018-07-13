#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*******************************/

//Print the error message 's' and exit(-1)
void terror(char *s);
// Check if character valid nucleotide annotation
int validc(char c);
// Create new FASTA with 1 Sequence made by all the Sequences from the Input File
void mergeMultiFasta(char* i_filename, char* o_filename);

/*******************************/

int main(int argc, char** argv){
	if (argc < 3)
			terror("USE: ./mergeMultiFasta <multifasta> <output>");

	mergeMultiFasta(argv[1], argv[2]);
	printf("mergeMultiFasta succesfully executed...\n");
}

/*******************************/

//Print the error message 's' and exit(-1)
void terror(char *s) {
	printf("ERROR::**** %s ****\n", s);
	exit(-1);
}

// Check if character valid nucleotide annotation
int validc(char c){
	//if( c=='A' || c== 'C' || c== 'G' || c== 'T')return 1;
	if( c=='A' || c== 'C' || c== 'G' || c== 'T' || c=='N' || c=='U' )return 1;
	//fprintf(stderr,"Warning: -%c- in seq",c);
	return 0;
}

// Create new FASTA with 1 Sequence made by all the Sequences from the Input File
void mergeMultiFasta(char* i_filename, char* o_filename){
	FILE * fi,* fo;
	char id_buffer[1000];
	char c;

	// Open Files
	if((fi=fopen(i_filename,"rt"))==NULL){
		terror("Can't open INPUT file");
	}
	if((fo=fopen(o_filename,"wt"))==NULL){
		terror("Can't open OUTPUT file");
	}

	// Output 1st Line
	fprintf(fo, ">0 ORIGINAL: %s\n", i_filename);
	// Read Sequence by Sequence and Output it
	c=getc(fi);
	while(!feof(fi)){
		// Skip ID
		if(c=='>'){
			while(c!='\n'){c=getc(fi);}
			c=getc(fi);
		}
		// Read Seq and Copy to output
		if(validc(c)) fputc(c,fo);

		c=getc(fi);
	}

	fclose(fi);
	fclose(fo);
}