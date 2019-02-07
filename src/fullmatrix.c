
#include <stdio.h>
#include <ctype.h>
#include "bioplib/general.h"
#include "bioplib/macros.h"
#include "subgroup.h"

#define ILETTER(x) ((int)(x)-65)
#define LETTER(x) ((char)((x)+65))


static void PopulateTopScores(FMSUBGROUPINFO *subGroupInfo);



int ReadFullMatrix(FILE *fp, FMSUBGROUPINFO *fullMatrix)
{
   char       buffer[MAXBUFF],
              *chp;
   int        entryCount = 0;
   BOOL       inData = FALSE;
   
   while(fgets(buffer, MAXBUFF, fp))
   {
      if(buffer[0] == '#')
      {
         continue; /* Skip comment lines */
      }
      else if(buffer[0] == '>')
      {
         if(entryCount >= MAXSUBTYPES)
         {
            fprintf(stderr,"Error: too many subtypes in data file. Increase MAXSUBTYPES\n");
            exit(1);
         }
         
         inData = TRUE;
         chp = buffer+1;
         sscanf(chp, "%s %d", fullMatrix[entryCount].type, &(fullMatrix[entryCount].index));

         switch(fullMatrix[entryCount].type[0])
         {
         case 'L':
            fullMatrix[entryCount].chainType = CHAINTYPE_LAMBDA;
            break;
         case 'K':
            fullMatrix[entryCount].chainType = CHAINTYPE_KAPPA;
            break;
         case 'H':
            fullMatrix[entryCount].chainType = CHAINTYPE_HEAVY;
            break;
         default:
            return(0);
         }
         
      }
      else if (inData && (buffer[0] == '"'))
      {
         chp = buffer+1;
         TERMAT(chp, '"');
         strncpy(fullMatrix[entryCount].name, chp, MAXBUFF);
      }
      else if ((buffer[0] == '/') && (buffer[1] == '/'))
      {
         PopulateTopScores(&(fullMatrix[entryCount]));

         entryCount++;
         inData = FALSE;
      }
      else if (inData)
      {
         char aa    = buffer[0];
         int  aaidx = ILETTER(aa);
         char word[MAXWORD];
         int  valCount = 0;
         
         chp        = buffer+2;
         while((chp = blGetWord(chp, word, MAXWORD))!=NULL)
         {
            fullMatrix[entryCount].scores[valCount++][aaidx] = atof(word);
         }
      }
   }

   return(entryCount);
}

static void PopulateTopScores(FMSUBGROUPINFO *subGroupInfo)
{
   int aaNum, posNum;
   REAL maxScore = 0.0;
   

   for(posNum=0; posNum<MAXREFSEQLEN; posNum++)
   {
      maxScore = 0.0;
      
      for(aaNum=0; aaNum<26; aaNum++)
      {
         if((aaNum==ILETTER('B')) ||
            (aaNum==ILETTER('J')) ||
            (aaNum==ILETTER('O')) ||
            (aaNum==ILETTER('U')) ||
            (aaNum==ILETTER('X')) ||
            (aaNum==ILETTER('Z'))) continue;
         if(subGroupInfo->scores[posNum][aaNum] > maxScore)
            maxScore = subGroupInfo->scores[posNum][aaNum];
      }
      subGroupInfo->topScores[posNum] = maxScore;
   }
}


REAL CalcFullScore(FMSUBGROUPINFO subGroupInfo, char *sequence, int offset, int offsetType)
{
   REAL score    = 0.0,
        scoreMax = 0.0;
   int i;
   
   if(offsetType == OFFSETTRUNCATION)
   {
      for(i = 0; i < (MAXREFSEQLEN-offset); i++) 
      {
         score += subGroupInfo.scores[i+offset][ILETTER(sequence[i])];
                                                
         /* Calculate the best score we could get against the most common 
            residues
         */
         scoreMax += subGroupInfo.topScores[i+offset];
      }
   }
   else  /* OFFSETEXTENSION                                             */
   {
      /* Return score of zero if we don't have enough residues in the
         test sequence
      */
      if(strlen(sequence) < (MAXREFSEQLEN + offset))
         return(0.0);
      
      for(i = 0; i < MAXREFSEQLEN; i++) 
      {
         score += subGroupInfo.scores[i][ILETTER(sequence[i+offset])];
         
         /* Calculate the best score we could get against the most common 
            residues
         */
         scoreMax += subGroupInfo.topScores[i];
      }
   }

   return((score*100.0)/scoreMax);
}


#ifdef TEST_READFMSUBGROUPINFO
int main(int argc, char **argv)
{
   FILE *fp;
   if((fp=fopen("../data/human_full.dat", "r"))!=NULL)
   {
      FMSUBGROUPINFO fullMatrix[MAXSUBTYPES];
      int        nEntries, entryNum;
      
      nEntries = ReadFullMatrix(fp, fullMatrix);
      for(entryNum=0; entryNum<nEntries; entryNum++)
      {
         int aaNum, position;
         fprintf(stdout, "Class: %s Index: %d\n", fullMatrix[entryNum].type, fullMatrix[entryNum].index);
         fprintf(stdout, "Name: %s\n", fullMatrix[entryNum].name);
         for(aaNum=0; aaNum<26; aaNum++)
         {
            if((aaNum==ILETTER('B')) ||
               (aaNum==ILETTER('J')) ||
               (aaNum==ILETTER('O')) ||
               (aaNum==ILETTER('U')) ||
               (aaNum==ILETTER('X')) ||
               (aaNum==ILETTER('Z'))) continue;
                fprintf(stdout, "%c", LETTER(aaNum));
            for(position=0; position<MAXREFSEQLEN; position++)
            {
               fprintf(stdout, "%6.3f", fullMatrix[entryNum].scores[position][aaNum]);
            }
            fprintf(stdout, "\n");
         }

         fprintf(stdout, "Top Scores:\n ");
         for(position=0; position<MAXREFSEQLEN; position++)
         {
            fprintf(stdout, "%6.3f", fullMatrix[entryNum].topScores[position]);
         }
         fprintf(stdout, "\n");
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
#endif
