#include  <stdio.h>
#include  <string.h>
#include  <stddef.h>
#include  <stdlib.h>

#define MAX_NB_SEQ 26
#define MAX_SEQ_GEN 21

char *tab_seq_gen[MAX_NB_SEQ] = {
"XDIQMTQSPSSLSASVGDRVT", "ZBVZLMZAATTVPLTPRESAI" ,
"DIVMTQSPLSLPVTPGEPASI" , "DVILTQTPLSSSGTLVQPSAI" ,
"EIVLTQSPGTLSLSPGERATL" , "DTLMRZVPASMCVTVGZKVAI" ,
"DIVMTQSPDSLAVSLGERATI" , "DLVLSQSPBTLAVSPGDQATV" ,
"ZSVLTQPPSVSGAPGQRVTIS" , "QALLTZPSSASATSGEKASLT" ,
"ZSALTQPASVSGSPGQSITIS" , "HVIVAZSPRATATLGATVKVT" ,
"XSYELTQPPSVSVSPGQTARI" , "XFFVVSZASVLFLAALZPVSA" ,
"SELTQDPAVSVALGQTVRITC" , "SALVQPASVZGSPGZSASIGC" ,
"ZSALTQPPSASGSPGQSVTIS" , "ZSALTQPPSASGSLGQSVTIS" ,
"NFMLTQPHSVSESPGKTVTIS" , "DLILIEPLSLSDSPEQKIIFS",
"XQVQLVQSGAEVKKPGASVKV" , "XZMHVLASASDLNRLPETLRI" ,
"QVQLQESGPGLVKPSQTLSLT" , "ZLTVRQWSAALLRATEAFTVI" ,
"XEVQLVESGGGLVQPGGSLRL" , "XQMHALQXTADVIKAERFMRV"
};

float tab_stat_gen[26][21] = {
0.016,0.748,0.759,0.710,0.721,0.819,0.803,0.819,0.814,0.776,0.579,
0.841,0.879,0.672,0.798,0.699,0.819,0.639,0.743,0.683,0.803, /* HKL1 */
0.005,0.022,0.016,0.022,0.087,0.005,0.038,0.005,0.005,0.011,0.289,
0.033,0.022,0.093,0.022,0.120,0.005,0.109,0.027,0.114,0.032, /* HKL1B */
0.980,0.921,0.921,0.902,0.921,0.902,0.568,0.843,0.862,0.843,0.784,
0.588,0.784,0.666,0.608,0.607,0.392,0.686,0.686,0.686,0.686, /* HKL2 */
0.000,0.019,0.019,0.039,0.000,0.000,0.235,0.000,0.000,0.000,0.039,
0.235,0.019,0.000,0.078,0.039,0.274,0.000,0.000,0.000,0.000, /* HKL2b */
0.801,0.801,0.880,0.722,0.900,0.861,0.894,0.894,0.477,0.834,0.847,
0.841,0.680,0.814,0.841,0.961,0.859,0.867,0.828,0.851,0.914, /* HKL3 */
0.039,0.013,0.006,0.198,0.006,0.039,0.006,0.000,0.285,0.013,0.006,
0.006,0.165,0.019,0.006,0.000,0.047,0.023,0.094,0.023,0.007, /* HKL3B */
0.944,0.888,0.944,0.833,0.888,0.944,0.888,0.944,0.555,0.722,0.888,
0.888,0.888,0.722,0.555,0.722,0.388,0.555,0.777,0.777,0.611, /* HKL4 */
0.000,0.055,0.000,0.111,0.055,0.000,0.000,0.000,0.166,0.111,0.000,
0.000,0.000,0.000,0.166,0.000,0.222,0.111,0.000,0.000,0.111, /* HKL4B */
0.585,0.878,0.878,0.926,0.951,0.902,0.926,0.902,0.926,0.487,0.951,
0.512,0.512,0.902,0.902,0.853,0.488,0.731,0.756,0.829,0.829, /* HLL1 */
0.292,0.024,0.024,0.000,0.000,0.024,0.000,0.024,0.000,0.390,0.000,
0.512,0.390,0.024,0.000,0.024,0.268,0.170,0.048,0.024,0.024, /* HLL1B */
0.736,0.789,0.736,0.947,0.894,0.947,0.894,0.605,0.973,0.894,0.947,
0.763,0.973,0.763,1.00,0.763,0.921,0.579,0.815,0.579,0.605, /* HLL2 */
0.184,0.158,0.131,0.026,0.052,0.026,0.052,0.289,0.026,0.052,0.026,
0.184,0.026,0.184,0.000,0.131,0.052,0.368,0.157,0.157,0.184, /* HLL2B */
0.032,0.290,0.919,0.467,0.887,0.854,0.887,0.822,0.919,0.887,0.790,
0.887,0.854,0.629,0.806,0.822,0.742,0.822,0.806,0.354,0.806, /* HLL3 */
0.000,0.032,0.032,0.209,0.048,0.032,0.064,0.048,0.016,0.016,0.08,
0.016,0.032,0.193,0.048,0.016,0.048,0.016,0.016,0.290,0.016, /* HLL3B */
0.777,0.666,0.888,0.777,0.888,0.444,0.777,0.444,0.888,0.777,0.777,
0.666,0.555,0.888,0.666,0.777,0.555,0.555,0.888,0.777,0.777, /* HLL4 */
0.000,0.111,0.000,0.111,0.000,0.333,0.111,0.333,0.000,0.111,0.111,
0.222,0.333,0.000,0.222,0.111,0.222,0.111,0.000,0.111,0.000, /* HLL4B */
1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,
1.000,1.000,0.666,1.000,1.000,1.000,1.000,1.000,1.000,1.000, /* HLL5 */
0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,
0.000,0.000,0.333,0.000,0.000,0.000,0.000,0.000,0.000,0.000, /* HLL5B */
0.642,0.857,0.928,1.000,0.857,0.928,1.000,0.714,1.000,0.928,1.000,
0.785,0.857,0.928,0.857,0.785,0.785,0.785,0.642,0.642,1.000, /* HLL6 */
0.357,0.142,0.071,0.000,0.071,0.071,0.000,0.071,0.000,0.071,0.000,
0.071,0.000,0.000,0.071,0.071,0.071,0.071,0.142,0.071,0.000, /* HLL6B */
0.008,0.443,0.756,0.747,0.782,0.634,0.539,0.704,0.686,0.643,0.686,
0.669,0.591,0.686,0.704,0.686,0.339,0.669,0.486,0.591,0.460, /* HHC1 */
0.000,0.200,0.017,0.026,0.008,0.034,0.130,0.000,0.017,0.026,0.017,
0.043,0.052,0.008,0.008,0.008,0.182,0.008,0.173,0.104,0.200, /* HHC1B */
0.594,0.554,0.524,0.722,0.564,0.514,0.643,0.702,0.623,0.613,0.702,
0.653,0.623,0.673,0.584,0.376,0.673,0.693,0.603,0.702,0.712, /* HHC2 */
0.099,0.128,0.099,0.009,0.059,0.188,0.059,0.009,0.049,0.049,0.000,
0.059,0.089,0.009,0.099,0.267,0.019,0.009,0.089,0.009,0.019, /* HHC2B */
0.004,0.712,0.828,0.768,0.847,0.643,0.768,0.810,0.842,0.879,0.699,
0.699,0.800,0.546,0.754,0.736,0.620,0.703,0.750,0.662,0.703, /* HHC3 */
0.000,0.069,0.027,0.032,0.009,0.199,0.046,0.031,0.007,0.004,0.092,
0.115,0.041,0.115,0.009,0.009,0.106,0.009,0.004,0.069,0.046, /* HHC3B */
};

char *tab_res_gen[MAX_NB_SEQ/2] = {
"Human Kappa Light chain subgroup I  ",
"Human Kappa Light chain subgroup II ",
"Human Kappa Light chain subgroup III",
"Human Kappa Light chain subgroup IV ",
"Human Lambda Light chain subgroup I  ",
"Human Lambda Light chain subgroup II ",
"Human Lambda Light chain subgroup III",
"Human Lambda Light chain subgroup IV ",
"Human Lambda Light chain subgroup V  ",
"Human Lambda Light chain subgroup VI ",
"Human Heavy chain subgroup I  ",
"Human Heavy chain subgroup II ",
"Human Heavy chain subgroup III"
};

main ()
{
	char tseq[200];
	long *class;
	long *sgpe;

	saisie_seq(tseq);
	det_sgpe(tseq,&class,&sgpe);
}


saisie_seq(tseq)
        char tseq[200];

{
        long i,j=0;
        long choix;

        printf(" Enter your sequence (one letter code) : \n");
        scanf("%s",tseq);
        i = 0;
        while (tseq[i] != '\0') {
                tseq[i] = toupper(tseq[i]);
                i++;
        }
}

float calc_stat(sgp,tseq,deb)
        long sgp,deb;
        char tseq[200];

{
        float val,valmax;
        long i;

        val = 0;
        valmax = 0;
        for (i=0;i<MAX_SEQ_GEN-deb;i++) {
                if (tseq[i] == tab_seq_gen[sgp][i+deb]) {
                        val += tab_stat_gen[sgp][i+deb];
                }
                else {
                    if (tseq[i] == tab_seq_gen[sgp+1][i+deb]) {
                        val += tab_stat_gen[sgp+1][i+deb];
                    }
                }
                valmax += tab_stat_gen[sgp][i+deb];
        }
        return((val*100)/valmax);
}

det_sgpe(tseq,class,sgpe)
        char tseq[200];
        long *class;
        long *sgpe;

{
        long sgp1,sgp2;
        float val,max1,max2;
        long sgp,i,j;
        long sgp_deb,sgp_fin;

        max1 = 0;
        max2 = 0;
        sgp1 = -1;
        sgp2 = -1;
        sgp_deb = 0;
        sgp_fin = MAX_NB_SEQ;

        for (sgp=sgp_deb;sgp<sgp_fin;sgp+=2) { /* sous groupe interessant */
            for (i=0;i<6;i++) { /* decalage */
                val = calc_stat(sgp,tseq,i);
                if (val > max2) {
                    if (val > max1) {
                        max2 = max1;
                        max1 = val;
                        sgp2 = sgp1;
                        sgp1 = sgp;
                    }
                    else {
                        max2 = val;
                        sgp2 = sgp;
                    }
                }
            }
        }
        printf("\n%s \n ",tab_res_gen[sgp1/2]);
       if(strstr(tab_res_gen[sgp1/2], "Kappa Light chain subgroup I  ")
                 != NULL) {
                *sgpe = 1;
                *class = 1;
        }
        if(strstr(tab_res_gen[sgp1/2], "Kappa Light chain subgroup II ")
                != NULL) {
                *class = 1;
                *sgpe = 2;
        }
        if(strstr(tab_res_gen[sgp1/2], "Kappa Light chain subgroup III")
                != NULL) {
                *class = 1;
                *sgpe = 3;
        }
        if(strstr(tab_res_gen[sgp1/2], "Kappa Light chain subgroup IV ")
                != NULL) {
                *class = 1;
                *sgpe = 4;
        }
        if(strstr(tab_res_gen[sgp1/2], "Lambda Light chain subgroup I  ")
                != NULL) {
                *class = 2;
                *sgpe = 1;
        }
       if(strstr(tab_res_gen[sgp1/2], "Lambda Light chain subgroup II ")
                != NULL) {
                *class = 2;
                *sgpe = 2;
        }
        if(strstr(tab_res_gen[sgp1/2], "Lambda Light chain subgroup III")
                != NULL) {
                *class = 2;
                *sgpe = 3;
        }
        if(strstr(tab_res_gen[sgp1/2], "Lambda Light chain subgroup IV ")
                != NULL) {
                *class = 2;
                *sgpe = 4;
        }
        if(strstr(tab_res_gen[sgp1/2], "Lambda Light chain subgroup V  ")
                != NULL) {
                *class = 2;
                *sgpe = 5;
        }
        if(strstr(tab_res_gen[sgp1/2], "Lambda Light chain subgroup VI ")
                != NULL) {
                *class = 2;
                *sgpe = 6;
        }
        if(strstr(tab_res_gen[sgp1/2], "Heavy chain subgroup I  ") != NULL) {
                *class = 0;
                *sgpe = 1;
        }
        if(strstr(tab_res_gen[sgp1/2], "Heavy chain subgroup II ") != NULL) {
                *class = 0;
                *sgpe = 2;
        }
        if(strstr(tab_res_gen[sgp1/2], "Heavy chain subgroup III ") != NULL) {
                *class = 0;
                *sgpe = 3;
        }
}
