/*************************************************************************

   Program:    hsubgroup
   File:       hsubgroup.c
   
   Version:    V3.0
   Date:       12.02.19
   Function:   Assign human subgroups from antibody sequences in PIR file
   
   Copyright:  (c) Dr. Andrew C. R. Martin / UCL 1997-2019
   Author:     Dr. Andrew C. R. Martin
   Address:    Biomolecular Structure & Modelling Unit,
               Department of Biochemistry & Molecular Biology,
               University College,
               Gower Street,
               London.
               WC1E 6BT.
   EMail:      andrew@bioinf.org.uk
               
**************************************************************************

   This program is not in the public domain, but it may be copied
   according to the conditions laid out in the accompanying file
   COPYING.DOC

   The code may be modified as required, but any modifications must be
   documented so that the person responsible can be identified. If someone
   else breaks this code, I don't want to be blamed for code that does not
   work! 

   The code may not be sold commercially or included as part of a 
   commercial product except as described in the file COPYING.DOC.

**************************************************************************

   Description:
   ============

**************************************************************************

   Usage:
   ======

**************************************************************************

   Revision History:
   =================
   V1.0  16.06.97   Original
   V2.0  01.08.18   Complete rewrite of underlying code
   V2.1  27.08.18   Allows data to be read from a file and records best
                    and second-best scores
   V2.2  08.01.19   Fixes problem with DOS files
   V2.3  05.02.19   Added info to verbose output on the second best match
   V3.0  12.02.19   Added support for full matrices

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <ctype.h>
#include <math.h>
#include "bioplib/general.h"
#include "bioplib/macros.h"
#include "subgroup.h"


/************************************************************************/
/* Defines and macros
*/
#define ILETTER(x) ((int)(x)-65)
#define LETTER(x) ((char)((x)+65))

/************************************************************************/
/* Prototypes
*/
static void PopulateTopScores(FMSUBGROUPINFO *subGroupInfo);

/************************************************************************/
/*>int ReadFullMatrix(FILE *fp, FMSUBGROUPINFO *fullMatrix)
   --------------------------------------------------------
*//**
   \param[in]   *fp          File pointer for matrix file
   \param[out]  *fullMatrix  The populated matrix
   \return                   The number of matrix entries

   Reads a full-matrix representation of residue frequencies

-  12.02.19 Original   By: ACRM
*/
int ReadFullMatrix(FILE *fp, FMSUBGROUPINFO *fullMatrix)
{
   char       buffer[MAXBUFF],
              *chp;
   int        entryCount = 0;
   BOOL       inData = FALSE;
   

   for(entryCount=0; entryCount<MAXSUBTYPES; entryCount++)
   {
      int i, j;
      for(i=0; i<MAXREFSEQLEN; i++)
      {
         for(j=0; j<26; j++)
         {
            fullMatrix[entryCount].scores[i][j] = 0.0;
         }
         fullMatrix[entryCount].topScores[i] = 0.0;
      }
   }

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
            fprintf(stderr,"Error: too many subtypes in data file. \
Increase MAXSUBTYPES\n");
            exit(1);
         }
         
         inData = TRUE;
         chp = buffer+1;
         sscanf(chp, "%s %d",
                fullMatrix[entryCount].type,
                &(fullMatrix[entryCount].index));

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


/************************************************************************/
/*>static void PopulateTopScores(FMSUBGROUPINFO *subGroupInfo)
   -----------------------------------------------------------
*//**
   \param[in,out] *subGroupInfo  The subgroup information matrix

   Populates the top score information in the subgroup full matrix.

-  12.02.19 Original   By: ACRM
*/
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


/************************************************************************/
/*>REAL CalcFullScore(FMSUBGROUPINFO subGroupInfo, char *sequence, 
                      int offset, int offsetType, BOOL includeX)
   ---------------------------------------------------------------
*//**
   \param[in]   subGroupInfo   The full matrix information for a subgroup
   \param[in]   sequence       The sequence to test
   \param[in]   offset         Offset into the sequence
   \param[in]   offsetType     Truncation or extension
   \param[in]   includeX       Include X characters in calculation
   \return                     The score for this sequence

   Calculates the score for a sequence against a sub group matrix

-  12.02.19 Original   By: ACRM
-  13.02.19 Added X checking
*/
REAL CalcFullScore(FMSUBGROUPINFO subGroupInfo, char *sequence,
                   int offset, int offsetType, BOOL includeX)
{
   REAL score    = 0.0,
        scoreMax = 0.0;
   int i;
   
   if(offsetType == OFFSETTRUNCATION)
   {
      for(i = 0; i < (MAXREFSEQLEN-offset); i++) 
      {
         if((sequence[i] != 'X') || includeX)
         {
            score += subGroupInfo.scores[i+offset][ILETTER(sequence[i])];
                                                
            /* Calculate the best score we could get against the most 
               common residues
            */
            scoreMax += subGroupInfo.topScores[i+offset];
         }
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
         if((sequence[i+offset] != 'X') || includeX)
         {
            score += subGroupInfo.scores[i][ILETTER(sequence[i+offset])];
         
            /* Calculate the best score we could get against the most 
               common residues
            */
            scoreMax += subGroupInfo.topScores[i];
         }
      }
   }

   return((score*100.0)/scoreMax);
}

/************************************************************************/
void fmTakeLogs(FMSUBGROUPINFO *subGroupInfo, int nSubGroups)
{
   int sgNum, i, j;
   for(sgNum=0; sgNum<nSubGroups; sgNum++)
   {
      for(i=0; i<MAXREFSEQLEN; i++)
      {
         for(j=0; j<26; j++)
         {
            subGroupInfo[sgNum].scores[i][j] = 
               log(subGroupInfo[sgNum].scores[i][j] + 1);
         }
         subGroupInfo[sgNum].topScores[i] = 
            log(subGroupInfo[sgNum].topScores[i] + 1);
      }
   }
}


/************************************************************************/
#ifdef TEST_READFMSUBGROUPINFO
/*>int main(int argc, char **argv)
   -------------------------------
*//**
   Test code for reading the matrix

-  12.02.19 Original   By: ACRM
*/
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
         fprintf(stdout, "Class: %s Index: %d\n",
                 fullMatrix[entryNum].type, fullMatrix[entryNum].index);
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
               fprintf(stdout, "%6.3f",
                       fullMatrix[entryNum].scores[position][aaNum]);
            }
            fprintf(stdout, "\n");
         }

         fprintf(stdout, "Top Scores:\n ");
         for(position=0; position<MAXREFSEQLEN; position++)
         {
            fprintf(stdout, "%6.3f",
                    fullMatrix[entryNum].topScores[position]);
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
