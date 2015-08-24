/*************************************************************************

   Program:    hsubgroup
   File:       hsubgroup.c
   
   Version:    V1.0
   Date:       16.06.97
   Function:   Assign human subgroups from antibody sequences in PIR file
   
   Copyright:  (c) Dr. Andrew C. R. Martin / UCL 1997-2015
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
   V1.0    16.06.97  Original   By: ACRM
   V1.1    24.08.15  Updated for new bioplib

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <stdlib.h>
#include "bioplib/MathType.h"
#include "bioplib/seq.h"
#include "bioplib/general.h"


/************************************************************************/
/* Defines and macros
*/
#define MAXBUFF 160
#define MAXSEQ 8

/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/
int main(int argc, char **argv);
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile);
void Usage(void);

/* External */
extern void det_sgpe(char *tseq, long *class, long *sgpe);

/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
   Main program to read a PIR file and make human antibody subgroup
   assignments for each chain.

   12.06.97 Original   By: ACRM
   16.06.96 Fixed memory leak --- wasn't freeing sequence data
*/
int main(int argc, char **argv)
{
   FILE *in  = stdin,
        *out = stdout;
   char infile[MAXBUFF],
        outfile[MAXBUFF],
        *seqs[MAXSEQ];
   int  nchain, i;
   long class, sgpe;
   BOOL punct, error;
   
   if(ParseCmdLine(argc, argv, infile, outfile))
   {
      if(blOpenStdFiles(infile, outfile, &in, &out))
      {
         while((nchain=blReadPIR(in,FALSE,seqs,MAXSEQ,NULL,&punct,&error)))
         {
            for(i=0; i<nchain; i++)
            {
               det_sgpe(seqs[i], &class, &sgpe);
               free(seqs[i]);
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
/*>BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile)
   ---------------------------------------------------------------------
   Input:   int    argc         Argument count
            char   **argv       Argument array
   Output:  char   *infile      Input file (or blank string)
            char   *outfile     Output file (or blank string)
   Returns: BOOL                Success?

   Parse the command line
   
   12.06.97 Original    By: ACRM
*/
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile)
{
   argc--;
   argv++;

   infile[0] = outfile[0] = '\0';
   
   while(argc)
   {
      if(argv[0][0] == '-')
      {
         switch(argv[0][1])
         {
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
*/
void Usage(void)
{
   fprintf(stderr,"\nhsubgroup V1.1 (c) 1997-2015, Andrew C.R. Martin, \
UCL\n");
   fprintf(stderr,"Subgroup assignment code (c) Sophie Deret, \
Necker Entants Malade, Paris\n");
   fprintf(stderr,"\nUsage: sophie [in.pir [out.txt]]\n");
   fprintf(stderr,"Assigns sub-group information for antibody \
sequences\n\n");
}
