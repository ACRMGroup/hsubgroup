#include <stdio.h>
#include <ctype.h>
#include "bioplib/general.h"
#include "bioplib/macros.h"

#define MAXPOS 21
#define MAXBUFF 256
#define MAXWORD 8
#define MAXSUBCLASS 100
typedef struct _fullmatrix {
   char        type[MAXWORD];
   int         index;
   char        title[MAXBUFF];
   REAL        scores[MAXPOS][26];
} FULLMATRIX;

FULLMATRIX *ReadFullMatrix(FILE *fp, int *nEntries)
{
   char       buffer[MAXBUFF],
              *chp;
   int        entryCount = 0;
   BOOL       inData = FALSE;
   static FULLMATRIX fullMatrix[MAXSUBCLASS];
   
   while(fgets(buffer, MAXBUFF, fp))
   {
      if(buffer[0] == '#')
      {
         continue; /* Skip comment lines */
      }
      else if(buffer[0] == '>')
      {
         inData = TRUE;
         chp = buffer+1;
         sscanf(chp, "%s %d", fullMatrix[entryCount].type, &(fullMatrix[entryCount].index));
      }
      else if (inData && (buffer[0] == '"'))
      {
         chp = buffer+1;
         TERMAT(chp, '"');
         strncpy(fullMatrix[entryCount].title, chp, MAXBUFF);
      }
      else if ((buffer[0] == '/') && (buffer[1] == '/'))
      {
         entryCount++;
         inData = FALSE;
      }
      else if (inData)
      {
         char aa    = buffer[0];
         int  aaidx = (int)aa - 65;
         char word[MAXWORD];
         int  valCount = 0;
         
         chp        = buffer+2;
         while((chp = blGetWord(chp, word, MAXWORD))!=NULL)
         {
            fullMatrix[entryCount].scores[valCount++][aaidx] = atof(word);
         }
      }
   }

   *nEntries = entryCount;
   return(fullMatrix);
}

int main(int argc, char **argv)
{
   FILE *fp;
   if((fp=fopen("../data/human_full.dat", "r"))!=NULL)
   {
      FULLMATRIX *fullMatrix;
      int        nEntries, i;
      
      fullMatrix = ReadFullMatrix(fp, &nEntries);
      for(i=0; i<nEntries; i++)
      {
         int j, k;
         fprintf(stdout, "Class: %s Index: %d\n", fullMatrix[i].type, fullMatrix[i].index);
         fprintf(stdout, "Title: %s\n", fullMatrix[i].title);
         for(j=0; j<26; j++)
         {
            if((j==(int)'B'-65) ||
               (j==(int)'J'-65) ||
               (j==(int)'O'-65) ||
               (j==(int)'U'-65) ||
               (j==(int)'X'-65) ||
               (j==(int)'Z'-65)) continue;
            fprintf(stdout, "%c", (char)(j+65));
            for(k=0; k<MAXPOS; k++)
            {
               fprintf(stdout, "%6.3f", fullMatrix[i].scores[j][k]);
            }
            fprintf(stdout, "\n");
         }
      }
      fclose(fp);
   }
   else
   {
      fprintf(stderr, "Can't open ../data/human_full.dat\n");
      return(1);
   }
   return(0);
}
