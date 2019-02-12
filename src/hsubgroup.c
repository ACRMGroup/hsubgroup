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
#include <stdlib.h>
#include "bioplib/MathType.h"
#include "bioplib/seq.h"
#include "bioplib/general.h"
#include "subgroup.h"

/************************************************************************/
/* Defines and macros
*/
#define MAXSEQ 8


/************************************************************************/
/* Prototypes
*/
int main(int argc, char **argv);
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile,
                  char *dataFile, BOOL *verbose, BOOL *fullMatrix);
void Usage(void);


/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
   Main program to read a PIR file and make human antibody subgroup
   assignments for each chain.

   12.06.97 Original   By: ACRM
   16.06.97 Fixed memory leak --- wasn't freeing sequence data
   26.11.18 Added data file and verbose options
   12.02.18 Added full matrix support
*/
int main(int argc, char **argv)
{
   FILE *in     = stdin,
        *out    = stdout,
        *fpData = NULL;
   char infile[MAXBUFF],
        outfile[MAXBUFF],
        dataFile[MAXBUFF],
        *seqs[MAXSEQ];
   int  nchain, i,
        class, subGroup;
   BOOL punct, error, verbose, fullMatrix;

   dataFile[0] = '\0';
   if(ParseCmdLine(argc, argv, infile, outfile, dataFile, &verbose,
                   &fullMatrix))
   {
      FindSubgroupSetVerbose(verbose);
      
      if(dataFile[0] != '\0')
      {
         if((fpData=fopen(dataFile, "r"))==NULL)
         {
            fprintf(stderr, "hsubgroup Error: Unable to open data \
file (%s)\n", dataFile);
            return(1);
         }
      }
      
      if(blOpenStdFiles(infile, outfile, &in, &out))
      {
         while((nchain=blReadPIR(in,FALSE,seqs,MAXSEQ,NULL,
                                 &punct,&error)))
         {
            for(i=0; i<nchain; i++)
            {
               BOOL ok = FindHumanSubgroup(fpData, fullMatrix, seqs[i],
                                           &class, &subGroup);
               free(seqs[i]);
               if(!ok)
               {
                  fprintf(stderr, "hsubgroup Error: Unable to read data \
from data file (%s)\n", dataFile);
                  return(1);
               }
            }
         }
      }
   }
   else
   {
      Usage();
   }
   
   return(0);
}

/************************************************************************/
/*>BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile,
                     char *dataFile, BOOL *verbose, BOOL *fullMatrix)
   ---------------------------------------------------------------------
   Input:   int    argc         Argument count
            char   **argv       Argument array
   Output:  char   *infile      Input file (or blank string)
            char   *outfile     Output file (or blank string)
            char   *dataFile    Optional data file (or blank string)
            BOOL   *verbose     Verbose output from subgroup code
            BOOL   *fullMatrix  Data file is a full scoring matrix
   Returns: BOOL                Success?

   Parse the command line
   
   12.06.97 Original    By: ACRM
   26.11.18 Added data file and verbose options
   05.02.19 Added -f
*/
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile,
                  char *dataFile, BOOL *verbose, BOOL *fullMatrix)
{
   argc--;
   argv++;

   infile[0] = outfile[0]  = dataFile[0] = '\0';
   *verbose  = *fullMatrix = FALSE;
   
   while(argc)
   {
      if(argv[0][0] == '-')
      {
         switch(argv[0][1])
         {
         case 'd':
            argc--; argv++;
            if(!argc)
               return(FALSE);
            strcpy(dataFile, argv[0]);
            break;
         case 'v':
            *verbose = TRUE;
            break;
         case 'f':
            *fullMatrix = TRUE;
            break;
         default:
            return(FALSE);
            break;
         }
      }
      else
      {
         /* Check that there are only 1 or 2 arguments left             */
         if(argc > 2)
            return(FALSE);
         
         /* Copy the first to infile                                    */
         strcpy(infile, argv[0]);
         
         /* If there's another, copy it to outfile                      */
         argc--;
         argv++;
         if(argc)
            strcpy(outfile, argv[0]);
            
         return(TRUE);
      }
      argc--;
      argv++;
   }
   
   return(TRUE);
}


/************************************************************************/
/*>void Usage(void)
   ----------------
   Prints a usage message

   12.06.97 Original    By: ACRM
   27.11.18 V2.1
   08.01.19 V2.2
   05.02.19 V2.3
   05.02.19 V2.4
   12.02.19 V3.0
*/
void Usage(void)
{
   fprintf(stderr,"\nhsubgroup V3.0 (c) 1997-2019, Andrew C.R. Martin, \
UCL\n");
   fprintf(stderr,"Original subgroup assignment code (c) Sophie Deret, \
Necker Entants Malade, Paris\n");
   fprintf(stderr,"   Used with permission\n");
   
   fprintf(stderr,"\nUsage: hsubgroup [-d datafile [-f]][-v] [in.pir \
[out.txt]]\n");
   fprintf(stderr,"       -d Specify data file\n");
   fprintf(stderr,"       -f Data file is a full matrix\n");
   fprintf(stderr,"       -v Verbose - shows best and 2nd best scores\n");
   fprintf(stderr,"          and the second best match\n");
   fprintf(stderr,"\nAssigns sub-group information for antibody \
sequences\n\n");
}
