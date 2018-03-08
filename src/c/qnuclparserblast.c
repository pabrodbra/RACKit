/*********

File		qnuclparserblast.c
Author		EPW <estebanpw@uma.es>
Description	Converts blastn's defaults output to the used format in the workflow

INPUT	
		<blastfile.in> 	Alignments file produced by BLASTn
		<file.out>		Location to store the parsed file of alignments
OUTPUT
		<file.out>		Parsed file of alignments
		
NOTES
		Please forgive the messy code.

**********/


#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define MAX_SEQ 512
#define MAXLINE 1000



int main(int argc, char **av){
	FILE * f, * g;
	if(argc != 3){
		printf("USE: ./qnuclparserblast <blastfile.in> <file.out>\n");
		exit(-1);
	}
	
	f = fopen(av[1], "rt");
	if(f==NULL){
		printf("Could not open input file\n");
		exit(-1);
	}
	g = fopen(av[2], "wt");
	if(g==NULL){
		printf("Could not open output file\n");
		exit(-1);
	}
	
	int descNum = 0, nextMeta = 0, fragmentNum = 0;
	char c;
	char line[MAXLINE];
	
	char tab[2];
	tab[0]='\t';
	tab[1]='\0';
	
	int globalQuery1=0;
	int globalQuery2=0;
	int globalSubject1=0;
	int globalSubject2=0;
	
	
	fgets(line, MAXLINE, f);
	while(!feof(f)){

		char sequenceID[MAX_SEQ]="";
		int sequenceIDlength = 0;
		char genomeID[MAX_SEQ]="";
		int genomeIDlength = 0;
		while(strncmp(line, "Query= ", 6) != 0) fgets(line, MAXLINE, f);
		int guide=7;
		while(line[guide]!='\n'){
			if(line[guide]!='\n') sequenceID[sequenceIDlength++]=line[guide];
			guide++;
		}

		int justJumped = 0;
		while(nextMeta == 0 && !feof(f)){
			
			if(justJumped==0){
				fgets(line, MAXLINE, f);	
			}
			justJumped=0;
			if(line[0] == '>'){
				int fakeLength = 1;
				while(line[fakeLength] != ' '){
					genomeID[genomeIDlength]= line[fakeLength];
					genomeIDlength++;
					fakeLength++;
				}
				while(strncmp(line, "Length = ", 5) != 0){
					fgets(line,MAXLINE,f);
				}
				int k=7;
				char fullGenomeLength[64]="";
				while(line[k] != '\r' && line[k] != '\n'){
					if(line[k] != '\r' && line[k] != '\n')fullGenomeLength[k-7] = line[k];
					k++;
				}
				genomeID[genomeIDlength++] = '\0';
				fullGenomeLength[k++] = '\0';
				char header[MAXLINE]="";
				sequenceID[sequenceIDlength] = ' ';
				char mayor[2];
				mayor[0]='>';
				mayor[1]='\0';
				strcat(header, mayor);
				strcat(header, sequenceID);
				strcat(header, mayor);
				strcat(header, genomeID);
				strcat(header, " ");
				strcat(header, fullGenomeLength);
				strcat(header, "\n");				

				fwrite(header, sizeof(char), strlen(header), g);
				genomeIDlength = 0;
				genomeID[0]='\0';
				
				fgets(line, MAXLINE, f);
				int fragmentMoment=0;
				int igaps=0, egaps=0;
				char wrline[MAXLINE]="";
				while(line[0] != '>' && strncmp(line, "Query= ", 6) != 0 && !feof(f)){
					if(strncmp(line, " Score", 4) == 0){
						char expectFloat[64]="";
						char score[64]="";
						fragmentNum++;
						if(fragmentMoment != 0){
							
							char gapsbuffer[256]="", finalWRLINE[MAXLINE]="";
							int findStrand=0;
							while(wrline[findStrand] != '+'){
								finalWRLINE[findStrand]=wrline[findStrand];
								findStrand++;
							}
							sprintf(gapsbuffer, "%d", igaps);
							strcat(finalWRLINE, gapsbuffer);
							strcat(finalWRLINE, tab);
							sprintf(gapsbuffer, "%d", egaps);
							strcat(finalWRLINE, gapsbuffer);
							strcat(finalWRLINE, tab);
							int newlenWR=strlen(finalWRLINE);
							while(findStrand < strlen(wrline)){
								finalWRLINE[newlenWR++] = wrline[findStrand++];
							}
							
							
							fprintf(g, "%d\t%s\t%d\t%d\t%d\t%d\n",fragmentNum-1, finalWRLINE, globalQuery1, globalQuery2, globalSubject1, globalSubject2);
							igaps=0;
							egaps=0;
							finalWRLINE[0] ='\0';
							wrline[0]='\0';
							globalQuery1=0;
							globalQuery2=0;
							globalSubject1=0;
							globalSubject2=0;
							
						}
						fragmentMoment = 0;
						int scoreCount=0;
						int counter= 0;
						while(line[scoreCount] != '=') scoreCount++;
						scoreCount++;
						while(line[scoreCount] == ' ') scoreCount++;
						while(line[scoreCount] != ' '){
							score[counter++] = line[scoreCount++];
						}
						int lnk=strlen(line)-1;
						int copyFloat=0;
						while(line[lnk] != ' ') lnk--;
						while(line[lnk] != '\n' && line[lnk] != '\n'){
							lnk++;
							if(line[lnk] != '\n' && line[lnk] != '\n')expectFloat[copyFloat++] = line[lnk];
						}
						//strcat(wrline, score);
						//strcat(wrline, tab);
						//strcat(wrline, expectFloat);
						//strcat(wrline, tab);
						strcat(wrline, "0\t");
					}
					if(strncmp(line," Identities", 9) == 0){
						char identities[64]="";
						char fragLength[64]="";
						char percentage[64] ="";
						char gaps[64]="";
						int lkn = 0;
						while(line[lkn] != ' ') lkn++;
						lkn++;
						while(line[lkn] != ' ') lkn++;
						lkn++;
						while(line[lkn] != ' ') lkn++;
						lkn++;
						int firstNumber=lkn;
						while(line[firstNumber] != '/'){
							identities[firstNumber-lkn] = line[firstNumber];
							firstNumber++;
						}
						firstNumber++;
						int secondNumber=firstNumber;
						while(line[secondNumber] != ' '){
							fragLength[secondNumber-firstNumber] = line[secondNumber];
							secondNumber++;
						}
						secondNumber += 2;
						int thirdNumber = secondNumber;
						while(line[thirdNumber] != '%'){
							percentage[thirdNumber-secondNumber] = line[thirdNumber];
							thirdNumber++;
						}
						while(line[thirdNumber] != '=')thirdNumber++;
						thirdNumber++;
						int fourth=0;
						while(line[thirdNumber] != '/'){ //Basic gaps are no longer used
							gaps[fourth] = line[thirdNumber];
							thirdNumber++;
							fourth++;
						}
						
						strcat(wrline, identities);
						strcat(wrline, tab);
						strcat(wrline, fragLength);
						strcat(wrline, tab);
						strcat(wrline, percentage);
						strcat(wrline, tab);
					}
					
					if(strncmp(line, " Strand ", 5) == 0){
						int firstNumber = 0;
						char strand1;
						char strand2;
						while(line[firstNumber] != '=') firstNumber++;
						firstNumber+=1;
						if(line[firstNumber] == 'P') strand1 = '+'; else strand1 = '-';
						while(line[firstNumber] != '/') firstNumber++;
						firstNumber+=1;
						if(line[firstNumber] == 'P') strand2 = '+'; else strand2 = '-';
						
						char fakestrand[2];
						fakestrand[0] = strand1;
						fakestrand[1] = '\0';
						
						strcat(wrline, fakestrand);
						fakestrand[0] = strand2;
						fakestrand[1] = '\0';
						strcat(wrline, fakestrand);
					}
					if((strncmp(line, "Query: ", 5) == 0) || (strncmp(line, "Sbjct: ", 5) == 0)){
						char query1[64]="";
						char query2[64]="";
						char sbj1[64]="";
						char sbj2[64]="";
						
						if(strncmp(line, "Query: ", 5) == 0){
							
							int firstNumber = 7;
							int k=0;
							while(line[firstNumber] != ' '){
								query1[k++] = line[firstNumber];
								firstNumber++;
							}
							while(line[firstNumber] == ' ') firstNumber++;
							while(line[firstNumber] != ' '){
								if(line[firstNumber] == '-' && line[firstNumber-1] != '-') igaps++;
								if(line[firstNumber] == '-' && line[firstNumber-1] == '-') egaps++;
								firstNumber++;
							}
							firstNumber++;
							int secondNumber = firstNumber;
							while(line[secondNumber] != '\n' && line[secondNumber] != '\n'){
								query2[secondNumber-firstNumber] = line[secondNumber];
								secondNumber++;
							}
							if(fragmentMoment==0){
								globalQuery1 = atoi(query1);
								globalQuery2 = atoi(query2);
							}else{
								globalQuery2 = atoi(query2);
							}
						}
						if(strncmp(line, "Sbjct: ", 5) == 0){
							int firstNumber = 7;
							int k=0;
							while(line[firstNumber] != ' '){
								sbj1[k++] = line[firstNumber];
								firstNumber++;
							}
							while(line[firstNumber] == ' ') firstNumber++;
							while(line[firstNumber] != ' '){
								if(line[firstNumber] == '-' && line[firstNumber-1] != '-') igaps++;
								if(line[firstNumber] == '-' && line[firstNumber-1] == '-') egaps++;
								firstNumber++;
							}
							firstNumber++;
							int secondNumber = firstNumber;
							while(line[secondNumber] != '\n' && line[secondNumber] != '\n'){
								sbj2[secondNumber-firstNumber] = line[secondNumber];
								secondNumber++;
							}
							if(fragmentMoment==0){
								globalSubject1 = atoi(sbj1);
								globalSubject2 = atoi(sbj2);
								
							}else{
								globalSubject2 = atoi(sbj2);
							}	
							fragmentMoment++;
						}
						
							
					}
					
					fgets(line, MAXLINE, f);
					justJumped = 1;
				}
				
				
				char gapsbuffer[256]="", finalWRLINE[MAXLINE]="";
				int findStrand=0;
				while(wrline[findStrand] != '+'){
					finalWRLINE[findStrand]=wrline[findStrand];
					findStrand++;
				}
				sprintf(gapsbuffer, "%d", igaps);
				strcat(finalWRLINE, gapsbuffer);
				strcat(finalWRLINE, tab);
				sprintf(gapsbuffer, "%d", egaps);
				strcat(finalWRLINE, gapsbuffer);
				strcat(finalWRLINE, tab);
				int newlenWR=strlen(finalWRLINE);
				while(findStrand < strlen(wrline)){
					finalWRLINE[newlenWR++] = wrline[findStrand++];
				}
				
				fprintf(g, "%d\t%s\t%d\t%d\t%d\t%d\n",fragmentNum, finalWRLINE, globalQuery1, globalQuery2, globalSubject1, globalSubject2);
				wrline[0] = '\0';
				igaps=0;
				egaps=0;
				finalWRLINE[0]='\0';
				globalQuery1=0;
				globalQuery2=0;
				globalSubject1=0;
				globalSubject2=0;
				fragmentNum = 0;
				
			}
			else if(strncmp(line, "***** No hits found *****", 10) == 0){
			}else if(strncmp(line, "Query= ", 6) == 0){
				nextMeta = 1;
			}
		}
		nextMeta = 0;
	}
	fclose(f);
	fclose(g);
	return 0;
}
