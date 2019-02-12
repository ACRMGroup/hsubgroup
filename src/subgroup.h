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
#define MAXWORD           8  /* Small buffer for parsing files          */
#define CHAINTYPE_HEAVY   0
#define CHAINTYPE_KAPPA   1
#define CHAINTYPE_LAMBDA  2

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


typedef struct {
   REAL        scores[MAXREFSEQLEN][26];
   REAL        topScores[MAXREFSEQLEN];
   int         index, chainType, subGroup;
   char        type[MAXWORD];
   char        name[MAXBUFF];
} FMSUBGROUPINFO;


/************************************************************************/
/* Prototypes
*/
BOOL FindHumanSubgroup(FILE *fp, BOOL fullMatrix, char *testSequence,
                       int *chainType, int *subGroup);
void FindSubgroupSetVerbose(BOOL verbose);

/* Not for end-user use                                                 */
int ReadFullMatrix(FILE *fp, FMSUBGROUPINFO *fullMatrix);
REAL CalcFullScore(FMSUBGROUPINFO subGroupInfo, char *sequence,
                   int offset, int offsetType);
