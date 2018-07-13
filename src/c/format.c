/*********

File		format.c
Author		EPW <estebanpw@uma.es>
Description	Formats a .fasta database adding a small ID value

USAGE		<genomes_database.fasta>	The DB file with all genomes concatenated in fasta format
		<genomes_formatted.fasta>	The output DB file

**********/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#define MAXLINE	5000




int main(int argc, char ** av){
	if(argc != 3){ printf("USE: format <genomeDB> <genomeDBout>\n"); exit(-1);}
	FILE * f, * g;
	 
	f = fopen(av[1], "rt");
	if(f == NULL){
		printf("Error opening genome DB file\n");
		exit(-1);
	}

	g = fopen(av[2], "wt");
	if(g == NULL){
		printf("Error opening output file\n");
		exit(-1);
	}
	
	char cr;
	uint64_t gc = 0;
	while(!feof(f)){
		cr = getc(f);
		if(cr == '>'){
			//Entering new genome, write ref
			while(cr != ' '){
				fputc(cr, g);
				cr = getc(f);
			}
			fprintf(g, ".%d\n", gc);
			gc++;

			while(cr != '\n'){
				cr = getc(f);
			}
		}else{
			fputc(cr, g);
		}
	}
	
	fclose(f);
	fclose(g);
	printf("format succesfully executed...\n");
	return 0;
}



