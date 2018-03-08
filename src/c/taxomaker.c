/*********

File		taxomaker.c
Author		EPW <estebanpw@uma.es>
Description	Produces a custom bounded taxonomical file reading the DB. It is suggested to only set the ask option to one if using a relatively small DB, and if using a larger one, use a spreadsheet.

USAGE		<genomes_database.fasta>	The DB file with all genomes concatenated in fasta format
			<OPTION: 1/0>	1 to ask for every genome to be included as substrain from the previous one, 0 to just write them all different species
		
OUTPUT
			Writes a new file with the name of the DB and ".taxo" with the custom boundaries. 


**********/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#define MAXLINE	5000

//goes [a, b)

uint64_t asciiToUint64(const char *text){
	uint64_t number=0;

	for(;*text;text++){
		char digit=*text-'0';           
		number=(number*10)+digit;
	}
	return number;
}


int main(int argc, char ** av){
	if(argc != 3){ printf("USE: taxomaker genomeDB Ask?[1/0]\n"); exit(-1);}
	FILE * f, * g;
	 
	f = fopen(av[1], "rt");
	if(f == NULL){
		printf("Error opening genome DB file\n");
		exit(-1);
	}
	char ask = av[2][0];
	char decision;
	char gi[MAXLINE]="";
	char name[MAXLINE] = "";
	char show[MAXLINE] = "";
	char decision_str[MAXLINE] = "";
	strcat(name, av[1]);
	strcat(name, ".taxo");
	printf("Opening :%s\n", name);
	g = fopen(name, "wt");
	uint64_t str_idx, buf_idx, gLevel1=0, gLevel2=0, len=0, ref_tot=0;
	uint64_t position;
	char c;

	fgets(show, MAXLINE, f);
	while(!feof(f)){
		while(show[0]!='>' && !feof(f)){
			fgets(show, MAXLINE, f);
		}
		if(show[0]=='>'){
		
			position = ftello(f);
			c=getc(f);
			len = 0;
			while(c!='>' && !feof(f)){
				c=getc(f);
				if(c!='\n') len++;
			}
			fseek(f, position, SEEK_SET);
		
			str_idx = 1;
			buf_idx = 0;
			while(show[str_idx] != ' ') gi[buf_idx++] = show[str_idx++];
			str_idx++;
			gi[buf_idx] = '\0';
			buf_idx = 0;
			while(show[str_idx] != '\n') name[buf_idx++] = show[str_idx++];
			name[buf_idx] = '\0';
			
			if(ask == '1'){
				printf("Is %s %s substrain from previous? [t/f]", gi, name);
				fflush(stdin);
				scanf("%s", decision_str);
				decision = decision_str[0];
				if(decision == 't'){
					gLevel2++;
					fprintf(g, "%"PRIu64"\t%"PRIu64"\t%"PRIu64"\t%s\t%s\n", gLevel1, gLevel2, len, gi, name);
				}else{
					
					gLevel2 = 0;
					fprintf(g, "%"PRIu64"\t%"PRIu64"\t%"PRIu64"\t%s\t%s\n", gLevel1, gLevel2, len, gi, name);
					gLevel1++;
				}
			}else{
				//printf("%s %s\n", gi, name);
				
				fprintf(g, "%"PRIu64"\t%"PRIu64"\t%"PRIu64"\t%"PRIu64"\t%s\t%s\n", gLevel1, gLevel2, len, ref_tot, gi, name);
				gLevel1++;
				ref_tot += len;
			}
			
			


			fgets(show, MAXLINE, f);
		}
	}
	


	fclose(f);
	fclose(g);

	return 0;
}



