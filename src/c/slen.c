/* slen : histogram of sequences lengths

    ./slen fileIN  0/1: verbose
   -------------------------------------------------*/
#include "stdio.h"
#include "stdlib.h"
#include "ctype.h"

#define MAXLS  1000000
#define maxH   5000
#define MAXLID   200

struct Sequence
{ 
  char id[MAXLID+1];
  char *data;
};



int loadSeqDB(FILE *, struct Sequence *, int);
int** LoadScores(char *, int , char*, char*, int*, int*); 

void terror(char* m) {
    printf("\nERR**  %s ***\n",m); 
    exit(-1);
}



/* ===============================M A I N==============*/
void main(int ac,char **av){
  int HIS[maxH];
  FILE   *fdb;
  char fV=0;
  int fA=0, fC=0, fG=0, fT=0, fN=0, fNo=0;
  /* from varias.c---------------------------*/
  struct Sequence SSH;
  int Nseq=0, n0, totlen=0;
  int MaxL=0, i, van=0;
  int kk, kkk;

/* -----Load sequences to memory ----------------------------*/
  if (ac!=3)  
     terror("Use changeReadHeader slen file.fasta -v=0/1 verbose");
  fV = atoi(av[2]);

  if (( SSH.data = (char*) malloc(MAXLS*sizeof(char)))==NULL)
	terror("memory for sequence");


  if ((fdb=fopen(av[1],"rt"))==NULL)
     terror("file open error on sequences file");

   for (i=0;i<maxH;i++) HIS[i]=0;
   n0 = loadSeqDB(fdb, &SSH, 1);

   while (n0) {

      if (fV && (Nseq%100000)==0) printf("\n %7d: %40s   [%8d]",Nseq,SSH.id, n0);
      Nseq++;
      totlen+=n0;
      if (MaxL<n0) MaxL=n0;
for (kk=0; kk<n0;kk++) {
  if (SSH.data[kk]=='A') fA++; else
  if (SSH.data[kk]=='C') fC++; else
  if (SSH.data[kk]=='G') fG++; else
  if (SSH.data[kk]=='T') fT++; else
  if (SSH.data[kk]=='N') fN++; else fNo++;
}
      i=n0/100;
      if (i>=maxH) i=maxH-1;
      HIS[i]++;
      if (feof(fdb)) break;
      n0 = loadSeqDB(fdb, &SSH, 0);
  }

  fclose(fdb);
  printf("\nTot Seq = %d  MaxLen=%d  L.Prom=%d  No.Bases=%d\n",Nseq,MaxL,totlen/Nseq,totlen);
  printf("\nHISTOGRAMA-------------------\n");
  for (i=0;i<maxH;i++)
    { van+=HIS[i];
      if (HIS[i]) printf("[%5d]   %6d   %7d  %7.2f\n",i*100, HIS[i], van, van/(float)Nseq*100.);
    }

fprintf(stdout,"fA=%d\nfC=%d\nfG=%d\nfT=%d\nfN=%d\nfNo=%d\n",fA,fC,fG,fT,fN,fNo);

}

// FUNCTIONS=============================================
// fAst=1:first sequence  / =0 : other seqs 
int loadSeqDB(FILE *f, struct Sequence *sX, int fAst) {
  char c;
  int i,lon=0,lee=1,k=0;

  if (fAst) while((c=getc(f))!='>' && !feof(f)); // start seq
  if (feof(f)) return -1;
  while((c=getc(f))==' ');
  while(k< MAXLID && c!='\n' && c!=' ') {
    sX->id[k++] = c;
    c=getc(f);
  }
  sX->id[k]=0x00;
  while(c!='\n') c=getc(f);
  c=getc(f);
  while(c!='>' && !feof(f)) {
      c=toupper(c);
      if (c>64 && c<123)  sX->data[lon++]=c;
      c=getc(f);
      if (lon==MAXLS-1) terror("Sequence too long");
    }
  sX->data[lon]=0x00;
  return lon;
}
