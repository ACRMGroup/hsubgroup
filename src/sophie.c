/*************************************************************************

   Program:    hsubgroup
   File:       sophie.c
   
   Version:    V2.3
   Date:       05.02.19
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

   This code is based loosely on the original by Sophie Deret 
   (deret@ceylan.necker.fr) from part of her Subim program which you can
   download from http://www.bioinf.org.uk/abs/subim.tar.gz

   Sophie kindly allowed us to integrate her code into our programs and
   it is provided free and with no warranty. This is a complete
   rewrite of her code using a structure to maintain all the info
   related to each subgroup instead of diverse arrays. Some redundant
   steps were eliminated compared with the original. The only thing
   remaining is the algorithm (as described in Deret, S. et al (1995),
   "SUBIM: a Program for Analysing the Kabat Database and Determining 
   the Variability Subgroup of a new Immunoglobulin Sequence", Comput.
   Appl. Biosci., 11:435-439) and the consensus sequence data and
   scores.

**************************************************************************

   Usage:
   ======

**************************************************************************

   Revision History:
   =================
   V1.0  16.06.97   Original
   V2.0  01.08.18   Complete rewrite
   V2.1  27.11.18   Allows data to be read from a file and records best
                    and second-best scores
   V2.2  08.01.19   Fixed problem with DOS data files and corrected 
                    code to deal with data files having other than 13
                    subtypes!
   V2.3  05.02.19   Added info to verbose output on the second best match

*************************************************************************/

/* Includes
*/
#include <stdio.h>
#include <string.h>
#include <stddef.h>
#include <stdlib.h>
#include <ctype.h>
#include "bioplib/macros.h"
#include "bioplib/general.h"

/************************************************************************/
/* Defines and macros
*/
#define MAXSUBTYPES     100  /* The max number of subtypes              */
#define MAXREFSEQLEN     21  /* The length of the reference sequences   */
#define MAXTRUNCATION     6  /* Amount we can Nter truncate a sequence  */
#define MAXEXTENSION     20  /* Amount we can Nter extend a sequence    */
#define MAXBUFF         320  /* General purpose buffer                  */
#define OFFSETTRUNCATION  1  /* Using offsets for Nter truncation       */
#define OFFSETEXTENSION   2  /* Using offsets for Nter extension        */

/* Used to store info on a subgroup                                     */
typedef struct
{
   REAL topScores[MAXREFSEQLEN],
        secondScores[MAXREFSEQLEN];
   int  chainType,
        subGroup;
   char name[MAXBUFF],
        topSeq[MAXBUFF],
        secondSeq[MAXBUFF];
} SUBGROUPINFO;


/************************************************************************/
/* Globals
*/
static BOOL sVerbose = FALSE;


/************************************************************************/
/* Prototypes
*/
#include "sophie.h"
static REAL CalcScore(SUBGROUPINFO subGroupInfo, char *sequence, 
                      int offset, int offsetType);
static void InitSubgroupInfo(SUBGROUPINFO *subGroupInfo, int chainType, 
                             int subGroup,
                             char *name, 
                             char *topSeq, 
                             REAL tv0,  REAL tv1,  REAL tv2,  REAL tv3, 
                             REAL tv4,  REAL tv5,  REAL tv6,  REAL tv7, 
                             REAL tv8,  REAL tv9,  REAL tv10, REAL tv11,
                             REAL tv12, REAL tv13, REAL tv14, REAL tv15,
                             REAL tv16, REAL tv17, REAL tv18, REAL tv19,
                             REAL tv20, 
                             char *secondSeq, 
                             REAL sv0,  REAL sv1,  REAL sv2,  REAL sv3,
                             REAL sv4,  REAL sv5,  REAL sv6,  REAL sv7,
                             REAL sv8,  REAL sv9,  REAL sv10, REAL sv11,
                             REAL sv12, REAL sv13, REAL sv14, REAL sv15,
                             REAL sv16, REAL sv17, REAL sv18, REAL sv19,
                             REAL sv20);
static int InitializeAllSubgroups(SUBGROUPINFO *subGroupInfo);
static int ReadSubgroupData(FILE *fp, SUBGROUPINFO *subGroupInfo);


/************************************************************************/
/*>static void InitSubgroupInfo(SUBGROUPINFO *subGroupInfo, 
                                int chainType, int subGroup,
                                char *name, 
                                char *topSeq, 
                                REAL tv0, REAL tv1,  REAL tv2,  REAL tv3, 
                                REAL tv4, REAL tv5,  REAL tv6,  REAL tv7, 
                                REAL tv8, REAL tv9,  REAL tv10, REAL tv11,
                                REAL tv12,REAL tv13, REAL tv14, REAL tv15,
                                REAL tv16,REAL tv17, REAL tv18, REAL tv19,
                                REAL tv20, 
                                char *secondSeq, 
                                REAL sv0, REAL sv1,  REAL sv2,  REAL sv3,
                                REAL sv4, REAL sv5,  REAL sv6,  REAL sv7,
                                REAL sv8, REAL sv9,  REAL sv10, REAL sv11,
                                REAL sv12,REAL sv13, REAL sv14, REAL sv15,
                                REAL sv16,REAL sv17, REAL sv18, REAL sv19,
                                REAL sv20)
   -----------------------------------------------------------------------
*//**
   \param[out]   *subGroupInfo  Pointer to structure to be populated
   \param[in]    chainType      CHAINTYPE_HEAVY
                                CHAINTYPE_KAPPA
                                CHAINTYPE_LAMBDA
   \param[in]    subGroup       Subgroup number
   \param[in]    name           Text version fo chain type and subgroup
   \param[in]    topSeq         The reference sequence of top amino acid
                                at each position
   \param[in]    tv0-tv20       The scores for the 21 positions in topSeq
   \param[in]    secondSeq      The reference sequence for the second
                                amino acid at each position
   \param[in]    sv0-sv20       The scores for the 21 positions in 
                                secondSeq

   Initialize a subGroupInfo structure

-  01.08.18  Original   By: ACRM
*/
static void InitSubgroupInfo(SUBGROUPINFO *subGroupInfo, int chainType, 
                             int subGroup,
                             char *name, 
                             char *topSeq, 
                             REAL tv0,  REAL tv1,  REAL tv2,  REAL tv3, 
                             REAL tv4,  REAL tv5,  REAL tv6,  REAL tv7, 
                             REAL tv8,  REAL tv9,  REAL tv10, REAL tv11,
                             REAL tv12, REAL tv13, REAL tv14, REAL tv15,
                             REAL tv16, REAL tv17, REAL tv18, REAL tv19,
                             REAL tv20, 
                             char *secondSeq, 
                             REAL sv0,  REAL sv1,  REAL sv2,  REAL sv3,
                             REAL sv4,  REAL sv5,  REAL sv6,  REAL sv7,
                             REAL sv8,  REAL sv9,  REAL sv10, REAL sv11,
                             REAL sv12, REAL sv13, REAL sv14, REAL sv15,
                             REAL sv16, REAL sv17, REAL sv18, REAL sv19,
                             REAL sv20)
{
   strncpy(subGroupInfo->name,      name,      MAXBUFF);
   strncpy(subGroupInfo->topSeq,    topSeq,    MAXBUFF);
   strncpy(subGroupInfo->secondSeq, secondSeq, MAXBUFF);

   subGroupInfo->chainType        = chainType;
   subGroupInfo->subGroup         = subGroup;

   subGroupInfo->topScores[0]     = tv0;
   subGroupInfo->topScores[1]     = tv1;
   subGroupInfo->topScores[2]     = tv2;
   subGroupInfo->topScores[3]     = tv3;
   subGroupInfo->topScores[4]     = tv4;
   subGroupInfo->topScores[5]     = tv5;
   subGroupInfo->topScores[6]     = tv6;
   subGroupInfo->topScores[7]     = tv7;
   subGroupInfo->topScores[8]     = tv8;
   subGroupInfo->topScores[9]     = tv9;
   subGroupInfo->topScores[10]    = tv10;
   subGroupInfo->topScores[11]    = tv11;
   subGroupInfo->topScores[12]    = tv12;
   subGroupInfo->topScores[13]    = tv13;
   subGroupInfo->topScores[14]    = tv14;
   subGroupInfo->topScores[15]    = tv15;
   subGroupInfo->topScores[16]    = tv16;
   subGroupInfo->topScores[17]    = tv17;
   subGroupInfo->topScores[18]    = tv18;
   subGroupInfo->topScores[19]    = tv19;
   subGroupInfo->topScores[20]    = tv20;

   subGroupInfo->secondScores[0]  = sv0;
   subGroupInfo->secondScores[1]  = sv1;
   subGroupInfo->secondScores[2]  = sv2;
   subGroupInfo->secondScores[3]  = sv3;
   subGroupInfo->secondScores[4]  = sv4;
   subGroupInfo->secondScores[5]  = sv5;
   subGroupInfo->secondScores[6]  = sv6;
   subGroupInfo->secondScores[7]  = sv7;
   subGroupInfo->secondScores[8]  = sv8;
   subGroupInfo->secondScores[9]  = sv9;
   subGroupInfo->secondScores[10] = sv10;
   subGroupInfo->secondScores[11] = sv11;
   subGroupInfo->secondScores[12] = sv12;
   subGroupInfo->secondScores[13] = sv13;
   subGroupInfo->secondScores[14] = sv14;
   subGroupInfo->secondScores[15] = sv15;
   subGroupInfo->secondScores[16] = sv16;
   subGroupInfo->secondScores[17] = sv17;
   subGroupInfo->secondScores[18] = sv18;
   subGroupInfo->secondScores[19] = sv19;
   subGroupInfo->secondScores[20] = sv20;
}


/************************************************************************/
/*>static int InitializeAllSubgroups(SUBGROUPINFO *subGroupInfo)
   -------------------------------------------------------------
*//**
   \param[out]  *subGroupInfo   Array of SUBGROUPINFO structures to
                                be initialized

   Initializes all the subgroup information

-  01.08.18  Original   By: ACRM
-  08.01.19  Now returns the number of subgroups done
*/
static int InitializeAllSubgroups(SUBGROUPINFO *subGroupInfo)
{
   int nSubGroups = 0;
   
   InitSubgroupInfo(&subGroupInfo[0],  CHAINTYPE_KAPPA,  1, 
                    "Human Kappa Light chain subgroup I",    
                    "XDIQMTQSPSSLSASVGDRVT", 
                    0.016,0.748,0.759,0.710,0.721,0.819,0.803,
                    0.819,0.814,0.776,0.579,0.841,0.879,0.672,
                    0.798,0.699,0.819,0.639,0.743,0.683,0.803,
                    "ZBVZLMZAATTVPLTPRESAI",
                    0.005,0.022,0.016,0.022,0.087,0.005,0.038,
                    0.005,0.005,0.011,0.289,0.033,0.022,0.093,
                    0.022,0.120,0.005,0.109,0.027,0.114,0.032);
   nSubGroups++;

   InitSubgroupInfo(&subGroupInfo[1],  CHAINTYPE_KAPPA,  2, 
                    "Human Kappa Light chain subgroup II",   
                    "DIVMTQSPLSLPVTPGEPASI", 
                    0.980,0.921,0.921,0.902,0.921,0.902,0.568,
                    0.843,0.862,0.843,0.784,0.588,0.784,0.666,
                    0.608,0.607,0.392,0.686,0.686,0.686,0.686,
                    "DVILTQTPLSSSGTLVQPSAI",  
                    0.000,0.019,0.019,0.039,0.000,0.000,0.235,
                    0.000,0.000,0.000,0.039,0.235,0.019,0.000,
                    0.078,0.039,0.274,0.000,0.000,0.000,0.000);
   nSubGroups++;

   InitSubgroupInfo(&subGroupInfo[2],  CHAINTYPE_KAPPA,  3, 
                    "Human Kappa Light chain subgroup III",  
                    "EIVLTQSPGTLSLSPGERATL", 
                    0.801,0.801,0.880,0.722,0.900,0.861,0.894,
                    0.894,0.477,0.834,0.847,0.841,0.680,0.814,
                    0.841,0.961,0.859,0.867,0.828,0.851,0.914,
                    "DTLMRZVPASMCVTVGZKVAI",
                    0.039,0.013,0.006,0.198,0.006,0.039,0.006,
                    0.000,0.285,0.013,0.006,0.006,0.165,0.019,
                    0.006,0.000,0.047,0.023,0.094,0.023,0.007);
   nSubGroups++;

   InitSubgroupInfo(&subGroupInfo[3],  CHAINTYPE_KAPPA,  4, 
                    "Human Kappa Light chain subgroup IV",   
                    "DIVMTQSPDSLAVSLGERATI", 
                    0.944,0.888,0.944,0.833,0.888,0.944,0.888,
                    0.944,0.555,0.722,0.888,0.888,0.888,0.722,
                    0.555,0.722,0.388,0.555,0.777,0.777,0.611,
                    "DLVLSQSPBTLAVSPGDQATV",  
                    0.000,0.055,0.000,0.111,0.055,0.000,0.000,
                    0.000,0.166,0.111,0.000,0.000,0.000,0.000,
                    0.166,0.000,0.222,0.111,0.000,0.000,0.111);
   nSubGroups++;
   
   InitSubgroupInfo(&subGroupInfo[4],  CHAINTYPE_LAMBDA, 1, 
                    "Human Lambda Light chain subgroup I",   
                    "ZSVLTQPPSVSGAPGQRVTIS", 
                    0.585,0.878,0.878,0.926,0.951,0.902,0.926,
                    0.902,0.926,0.487,0.951,0.512,0.512,0.902,
                    0.902,0.853,0.488,0.731,0.756,0.829,0.829,
                    "QALLTZPSSASATSGEKASLT",  
                    0.292,0.024,0.024,0.000,0.000,0.024,0.000,
                    0.024,0.000,0.390,0.000,0.512,0.390,0.024,
                    0.000,0.024,0.268,0.170,0.048,0.024,0.024);
   nSubGroups++;

   InitSubgroupInfo(&subGroupInfo[5],  CHAINTYPE_LAMBDA, 2, 
                    "Human Lambda Light chain subgroup II",  
                    "ZSALTQPASVSGSPGQSITIS", 
                    0.736,0.789,0.736,0.947,0.894,0.947,0.894,
                    0.605,0.973,0.894,0.947,0.763,0.973,0.763,
                    1.000,0.763,0.921,0.579,0.815,0.579,0.605,
                    "HVIVAZSPRATATLGATVKVT",  
                    0.184,0.158,0.131,0.026,0.052,0.026,0.052,
                    0.289,0.026,0.052,0.026,0.184,0.026,0.184,
                    0.000,0.131,0.052,0.368,0.157,0.157,0.184);
   nSubGroups++;

   InitSubgroupInfo(&subGroupInfo[6],  CHAINTYPE_LAMBDA, 3, 
                    "Human Lambda Light chain subgroup III", 
                    "XSYELTQPPSVSVSPGQTARI", 
                    0.032,0.290,0.919,0.467,0.887,0.854,0.887,
                    0.822,0.919,0.887,0.790,0.887,0.854,0.629,
                    0.806,0.822,0.742,0.822,0.806,0.354,0.806,
                    "XFFVVSZASVLFLAALZPVSA",  
                    0.000,0.032,0.032,0.209,0.048,0.032,0.064,
                    0.048,0.016,0.016,0.08,0.016,0.032,0.193,
                    0.048,0.016,0.048,0.016,0.016,0.290,0.016);
   nSubGroups++;

   InitSubgroupInfo(&subGroupInfo[7],  CHAINTYPE_LAMBDA, 4, 
                    "Human Lambda Light chain subgroup IV",  
                    "SELTQDPAVSVALGQTVRITC", 
                    0.777,0.666,0.888,0.777,0.888,0.444,0.777,
                    0.444,0.888,0.777,0.777,0.666,0.555,0.888,
                    0.666,0.777,0.555,0.555,0.888,0.777,0.777,
                    "SALVQPASVZGSPGZSASIGC",  
                    0.000,0.111,0.000,0.111,0.000,0.333,0.111,
                    0.333,0.000,0.111,0.111,0.222,0.333,0.000,
                    0.222,0.111,0.222,0.111,0.000,0.111,0.000);
   nSubGroups++;
   
   InitSubgroupInfo(&subGroupInfo[8],  CHAINTYPE_LAMBDA, 5, 
                    "Human Lambda Light chain subgroup V",   
                    "ZSALTQPPSASGSPGQSVTIS", 
                    1.000,1.000,1.000,1.000,1.000,1.000,1.000,
                    1.000,1.000,1.000,1.000,1.000,1.000,0.666,
                    1.000,1.000,1.000,1.000,1.000,1.000,1.000,
                    "ZSALTQPPSASGSLGQSVTIS",  
                    0.000,0.000,0.000,0.000,0.000,0.000,0.000,
                    0.000,0.000,0.000,0.000,0.000,0.000,0.333,
                    0.000,0.000,0.000,0.000,0.000,0.000,0.000);
   nSubGroups++;

   InitSubgroupInfo(&subGroupInfo[9],  CHAINTYPE_LAMBDA, 6, 
                    "Human Lambda Light chain subgroup VI",  
                    "NFMLTQPHSVSESPGKTVTIS", 
                    0.642,0.857,0.928,1.000,0.857,0.928,1.000,
                    0.714,1.000,0.928,1.000,0.785,0.857,0.928,
                    0.857,0.785,0.785,0.785,0.642,0.642,1.000,
                    "DLILIEPLSLSDSPEQKIIFS",  
                    0.357,0.142,0.071,0.000,0.071,0.071,0.000,
                    0.071,0.000,0.071,0.000,0.071,0.000,0.000,
                    0.071,0.071,0.071,0.071,0.142,0.071,0.000);
   nSubGroups++;
   
   InitSubgroupInfo(&subGroupInfo[10], CHAINTYPE_HEAVY,  1, 
                    "Human Heavy chain subgroup I",          
                    "XQVQLVQSGAEVKKPGASVKV", 
                    0.008,0.443,0.756,0.747,0.782,0.634,0.539,
                    0.704,0.686,0.643,0.686,0.669,0.591,0.686,
                    0.704,0.686,0.339,0.669,0.486,0.591,0.460,
                    "XZMHVLASASDLNRLPETLRI",  
                    0.000,0.200,0.017,0.026,0.008,0.034,0.130,
                    0.000,0.017,0.026,0.017,0.043,0.052,0.008,
                    0.008,0.008,0.182,0.008,0.173,0.104,0.200);
   nSubGroups++;

   InitSubgroupInfo(&subGroupInfo[11], CHAINTYPE_HEAVY,  2, 
                    "Human Heavy chain subgroup II",         
                    "QVQLQESGPGLVKPSQTLSLT", 
                    0.594,0.554,0.524,0.722,0.564,0.514,0.643,
                    0.702,0.623,0.613,0.702,0.653,0.623,0.673,
                    0.584,0.376,0.673,0.693,0.603,0.702,0.712,
                    "ZLTVRQWSAALLRATEAFTVI",  
                    0.099,0.128,0.099,0.009,0.059,0.188,0.059,
                    0.009,0.049,0.049,0.000,0.059,0.089,0.009,
                    0.099,0.267,0.019,0.009,0.089,0.009,0.019);
   nSubGroups++;

   InitSubgroupInfo(&subGroupInfo[12], CHAINTYPE_HEAVY,  3, 
                    "Human Heavy chain subgroup III",        
                    "XEVQLVESGGGLVQPGGSLRL", 
                    0.004,0.712,0.828,0.768,0.847,0.643,0.768,
                    0.810,0.842,0.879,0.699,0.699,0.800,0.546,
                    0.754,0.736,0.620,0.703,0.750,0.662,0.703,
                    "XQMHALQXTADVIKAERFMRV", 
                    0.000,0.069,0.027,0.032,0.009,0.199,0.046,
                    0.031,0.007,0.004,0.092,0.115,0.041,0.115,
                    0.009,0.009,0.106,0.009,0.004,0.069,0.046);
   nSubGroups++;

   return(nSubGroups);
}


/************************************************************************/
/*static REAL CalcScore(SUBGROUPINFO subGroupInfo, char *sequence, 
                        int offset, int offsetType)
  ------------------------------------------------------------
*//**
   \param[in]  subGroupInfo - Information for the subgroup we are 
                              looking at
   \param[in]  sequence     - The sequence we are looking at
   \param[in]  offset       - Offset into the reference sequence or
                              sequence
   \param[in]  offsetType   - OFFSETTRUNCATION - account for the sequence
                                                 being Nter truncated
                              OFFSETEXTENSION  - account for the sequence
                                                 being Nter extended

   Calculates the score for the test sequence against a specified
   subgroup

-  16.06.97 Original from Sophie's code
-  01.08.18 Complete rewrite
*/
static REAL CalcScore(SUBGROUPINFO subGroupInfo, char *sequence, 
                      int offset, int offsetType)
{
   REAL score    = 0.0,
        scoreMax = 0.0;
   int i;
   
   if(offsetType == OFFSETTRUNCATION)
   {
      for(i = 0; i < (MAXREFSEQLEN-offset); i++) 
      {
         /* Test if the residue matches the most common one at this 
            position 
         */
         if(sequence[i] == subGroupInfo.topSeq[i+offset])
         {
            score += subGroupInfo.topScores[i+offset];
         }
         else /* If not, test if it matches the next most common        */
         {
            if(sequence[i] == subGroupInfo.secondSeq[i+offset])
            {
               score += subGroupInfo.secondScores[i+offset];
            }
         }
         
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
         /* Test if the residue matches the most common one at this 
            position 
         */
         if(sequence[i+offset] == subGroupInfo.topSeq[i])
         {
            score += subGroupInfo.topScores[i];
         }
         else /* If not, test if it matches the next most common        */
         {
            if(sequence[i+offset] == subGroupInfo.secondSeq[i])
            {
               score += subGroupInfo.secondScores[i];
            }
         }
         
         /* Calculate the best score we could get against the most common 
            residues
         */
         scoreMax += subGroupInfo.topScores[i];
      }
   }

   return((score*100.0)/scoreMax);
}


/************************************************************************/
/*>static int ReadSubgroupData(FILE *fp, SUBGROUPINFO *subGroupInfo)
   -----------------------------------------------------------------
*//**
   \param[in]    *fp            File pointer
   \param[out]   *subGroupInfo  The sub group information
   \return                      Number of subgroups (0 for error)

   Reads the subgroup information file which contains records of the
   form

>KAPPA,  1,
"Human Kappa Light chain subgroup I",
XDIQMTQSPSSLSASVGDRVT,
0.016,0.748,0.759,0.710,0.721,0.819,0.803,0.819,0.814,0.776,0.579,0.841,0.879,0.672,0.798,0.699,0.819,0.639,0.743,0.683,0.803,
ZBVZLMZAATTVPLTPRESAI,
0.005,0.022,0.016,0.022,0.087,0.005,0.038,0.005,0.005,0.011,0.289,0.033,0.022,0.093,0.022,0.120,0.005,0.109,0.027,0.114,0.032

- 27.11.18 Original   By: ACRM
//

*/
static int ReadSubgroupData(FILE *fp, SUBGROUPINFO *subGroupInfo)
{
   int  subGroupCount = 0,
        dataNum       = 0,
        chainType     = 0,
        chainTypeNum  = 0,
        nSubGroups    = 0;
   char buffer[MAXBUFF],
        label[MAXBUFF],
        seq1[MAXBUFF],
        seq2[MAXBUFF],
        *chp;
   REAL freq1[MAXREFSEQLEN],
        freq2[MAXREFSEQLEN];
   
   while(fgets(buffer, MAXBUFF, fp))
   {
      TERMINATE(buffer);            /* Terminate normally               */
      TERMAT(buffer, '\r');         /* Terminate at DOS cursor return   */
      TERMAT(buffer, '#');          /* Remove comments                  */
      KILLTRAILSPACES(buffer);      /* Trailing spaces                  */
      KILLLEADSPACES(chp, buffer);  /* Leading spaces                   */
      if(strlen(chp))               /* Anything left?                   */
      {
         if(chp[0] == '>')          /* Start of new entry               */
         {
            dataNum = 0;
            chp++;
            KILLLEADSPACES(chp, chp);
         }
         
         do
         {
            char word[MAXBUFF];
            chp=blGetWord(chp, word, MAXBUFF);

            /* Test for the end of a block */
            if((word[0] == '/') && (word[1] == '/'))
            {
               if(dataNum != 47)
               {
                  char chainTypeLabel[16];
                  switch(chainType)
                  {
                  case CHAINTYPE_HEAVY:
                     strcpy(chainTypeLabel, "HEAVY");
                     break;
                  case CHAINTYPE_KAPPA:
                     strcpy(chainTypeLabel, "KAPPA");
                     break;
                  case CHAINTYPE_LAMBDA:
                     strcpy(chainTypeLabel, "LAMBDA");
                     break;
                  default:
                     strcpy(chainTypeLabel, "????");
                     break;
                  }
                  
                  fprintf(stderr,"Datafile invalid at %s %d. \
Got %d fields instead of 46\n",
                          chainTypeLabel, chainTypeNum, dataNum);
                  return(0);
               }
               InitSubgroupInfo(&subGroupInfo[subGroupCount++],
                                chainType, chainTypeNum,
                                label,
                                seq1,
                                freq1[0],  freq1[1],  freq1[2],
                                freq1[3],  freq1[4],  freq1[5],
                                freq1[6],  freq1[7],  freq1[8],
                                freq1[9],  freq1[10], freq1[11], 
                                freq1[12], freq1[13], freq1[14],
                                freq1[15], freq1[16], freq1[17],
                                freq1[18], freq1[19], freq1[20],
                                seq2,
                                freq2[0],  freq2[1],  freq2[2],
                                freq2[3],  freq2[4],  freq2[5],
                                freq2[6],  freq2[7],  freq2[8],
                                freq2[9],  freq2[10], freq2[11], 
                                freq2[12], freq2[13], freq2[14],
                                freq2[15], freq2[16], freq2[17],
                                freq2[18], freq2[19], freq2[20]);
               if(++nSubGroups > MAXSUBTYPES)
               {
                  fprintf(stderr, "Data file contains too many types. \
Increase MAXSUBTYPES.\n");
                  exit(1);
               }
            }
            
            if(dataNum == 0)
            {
               switch(word[0])
               {
               case 'L':
                  chainType = CHAINTYPE_LAMBDA;
                  break;
               case 'K':
                  chainType = CHAINTYPE_KAPPA;
                  break;
               case 'H':
                  chainType = CHAINTYPE_HEAVY;
                  break;
               default:
                  return(0);
               }
            }
            else if(dataNum == 1)
            {
               sscanf(word, "%d", &chainTypeNum);
            }
            else if(dataNum == 2)
            {
               strncpy(label, word, MAXBUFF);
            }
            else if(dataNum == 3)
            {
               strncpy(seq1, word, MAXBUFF);
            }
            else if((dataNum >= 4) && (dataNum <= 24))
            {
               sscanf(word, "%lf", &(freq1[dataNum-4]));
            }
            else if(dataNum == 25)
            {
               strncpy(seq2, word, MAXBUFF);
            }
            else if((dataNum >= 26) && (dataNum <= 46))
            {
               sscanf(word, "%lf", &(freq2[dataNum-26]));
            }
            
            dataNum++;
         } while(chp!=NULL);
      }
   }
   return(nSubGroups);
}


/************************************************************************/
/*>void FindSubgroupSetVerbose(BOOL verbose)
   -----------------------------------------
*//*
   \param[in]    verbose    Verbose setting

   Sets the static 'verbose' variable

-  27.11.18  Original   By: ACRM
*/
void FindSubgroupSetVerbose(BOOL verbose)
{
   sVerbose = verbose;
}


/************************************************************************/
/*>BOOL FindHumanSubgroup(FILE *fp, char *sequence, int *chainType, 
                          int *subGroup)
   ----------------------------------------------------------------
*//**
   \param[in]   fp           - file of residue subgroup specifications
                               (NULL - use default hardcoded values)
   \param[in]   sequence     - the sequence of interest
   \param[out]  chainType    - chain type: CHAINTYPE_HEAVY
                                           CHAINTYPE_KAPPA
                                           CHAINTYPE_LAMBDA
   \param[out]  subGroup     - subgroup
   \return                   - Success in reading data file

   Assigns the subgroup information for a sequence

-  16.06.97 Original from Sophie's code
-  01.08.18 Complete rewrite
-  27.11.18 Now returns BOOL and can read file of residue frequencies
            Also deals with verbose printing
*/
BOOL FindHumanSubgroup(FILE *fp, char *sequence, int *chainType,
                       int *subGroup)
{
   static SUBGROUPINFO subGroupInfo[MAXSUBTYPES];
   static int          sInitialized            = 0,
                       sNSubGroups             = 0;
   int                 bestSubGroupCount       = -1,
                       secondBestSubGroupCount = -1;
   REAL                val                     = 0.0,
                       maxVal                  = 0.0,
                       secondMaxVal            = 0.0;
   int                 subGroupCount,
                       offset;
#ifdef DEBUG
   int                 bestOffset              = 0;
#endif
   
   if(!sInitialized)
   {
      sInitialized = 1;
      if(fp != NULL)
      {
         sNSubGroups = ReadSubgroupData(fp, subGroupInfo);
         if(!sNSubGroups)
            return(FALSE);
      }
      else
      {
         sNSubGroups = InitializeAllSubgroups(subGroupInfo);
      }
   }
   
   /* For each sub-group                                                */
   for(subGroupCount = 0; subGroupCount < sNSubGroups; subGroupCount++) 
   { 
      /* Shift along the reference sequence to account for N-terminal
         truncation of the test sequence
      */
      for(offset = 0; offset < MAXTRUNCATION; offset++)
      {
         val = CalcScore(subGroupInfo[subGroupCount], sequence, 
                         offset, OFFSETTRUNCATION);
         if(val > maxVal) 
         {
            maxVal            = val;
            bestSubGroupCount = subGroupCount;
#ifdef DEBUG
            bestOffset        = offset;
#endif
         }
         else if((val < maxVal) && (val > secondMaxVal))
         {
            secondMaxVal            = val;
            secondBestSubGroupCount = subGroupCount;
         }
         
      }

      /* Shift along the test sequence to account for N-terminal 
         extension of the test sequence
      */
      for(offset = 0; offset < MAXEXTENSION; offset++)
      {
         val = CalcScore(subGroupInfo[subGroupCount], sequence, 
                         offset, OFFSETEXTENSION);
         if(val > maxVal) 
         {
            maxVal            = val;
            bestSubGroupCount = subGroupCount;
#ifdef DEBUG
            bestOffset        = -offset;
#endif
         }
         else if((val < maxVal) && (val > secondMaxVal))
         {
            secondMaxVal            = val;
            secondBestSubGroupCount = subGroupCount;
         }
      }
   }

   /* Print the winning name                                            */
   printf("%s",subGroupInfo[bestSubGroupCount].name);
   if(sVerbose)
   {
      printf(" [%f %f] ", maxVal, secondMaxVal);
      printf("(%s)",subGroupInfo[secondBestSubGroupCount].name);
   }
   printf("\n");
#ifdef DEBUG
   printf("Offset: %d (%s)\n", bestOffset, 
          ((bestOffset<0)?"extension":"truncation"));
#endif

   /* Set the chain type and sub group                                  */
   *chainType = subGroupInfo[bestSubGroupCount].chainType;
   *subGroup  = subGroupInfo[bestSubGroupCount].subGroup;

   return(TRUE);
}


#ifdef DEMO
/************************************************************************/
/* Prototypes
*/
int main(void);
static void RequestSequence(char *sequence);

int main(void)
{
   char sequence[200];
   int  chainType, subGroup;
   
   RequestSequence(sequence);
   FindHumanSubgroup(sequence, &chainType, &subGroup);
   
   return(0);
}


/************************************************************************/
static void RequestSequence(char *sequence)
{
   int i;
   
   printf(" Enter your sequence (one letter code) : \n");
   scanf("%s",sequence);
   i = 0;
   while (sequence[i] != '\0') 
   {
      sequence[i] = toupper(sequence[i]);
      i++;
   }
}

#endif
