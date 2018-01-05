
#include <stdint.h>
#include <stdio.h>

#include "tldevel.h"

#include "delve_struct.h"

char* glocal_viterbi_log_Y(struct hmm* hmm,char* a, char* b, int n,int m)
{
        int i,j,c;
	
        unsigned short int tmp_tb;
        float** M = hmm->fM;
        float** X = hmm->fX;
        float** Y = hmm->fY;
        char* seqa = a -1;
        char* seqb = b -1;
        float* prob = hmm->prob;
        unsigned short int** t = hmm->traceback;
        char* alignment = NULL;


        MMALLOC(alignment, sizeof( char)* (n+m+2));
	
        t[0][2] = 4 << 6;
	
        c = 0;
        M[0][0] = prob2scaledprob(0.0f);
        X[0][0] = prob2scaledprob(0.0f);
        Y[0][0] = prob2scaledprob(0.0f);
	
        M[0][1] = prob2scaledprob(0.0f);
        X[0][1] = prob2scaledprob(0.0f);
        Y[0][1] = prob2scaledprob(0.0f);
        tmp_tb = (4 << 6);
        t[0][1] = tmp_tb;
	
        M[0][2] = prob2scaledprob(0.0f);
        X[0][2] = prob2scaledprob(0.0f);
        Y[0][2] = prob[   (EGAP << 12) | (seqb[1] << 8) | (EGAP << 4) | seqb[2]];
        tmp_tb = (4 << 6);
        t[0][2] = tmp_tb;
	
        for(j = 3; j <= m-2;j++){
                M[0][j] = prob2scaledprob(0.0f);
                X[0][j] = prob2scaledprob(0.0f);
                Y[0][j] = prob[   (EGAP << 12) | (seqb[j-1] << 8) | (EGAP << 4) | seqb[j]] + Y[0][j-1];
                tmp_tb = (4 << 6);
                t[0][j] = tmp_tb;
        }
	
        M[0][m-1] = prob2scaledprob(0.0f);
        X[0][m-1] = prob2scaledprob(0.0f);
        Y[0][m-1] = prob2scaledprob(0.0f);
        t[0][m-1] = 0;
	
        M[0][m] = prob2scaledprob(0.0f);
        X[0][m] = prob2scaledprob(0.0f);
        Y[0][m] = prob2scaledprob(0.0f);
        t[0][m] = 0;
	
        i = 1;
	
        M[i][0] = prob2scaledprob(0.0f);
        X[i][0] = prob2scaledprob(0.0f);
        Y[i][0] = prob2scaledprob(0.0f);
        t[i][0] = 0;
        M[i][1] = prob2scaledprob(0.0f);
        X[i][1] = prob2scaledprob(0.0f);
        Y[i][1] = prob2scaledprob(0.0f);
	
        t[i][1] = 0;
	
	
        for(j = 2; j <= m-2;j++){
                M[i][j] = prob[              (EGAP << 12) | (seqb[j-1] << 8) | (seqa[i] << 4) | (seqb[j])] + Y[i-1][j-1];
                tmp_tb = 4;
                t[i][j] = tmp_tb;
		
                X[i][j] = prob[              (EGAP << 12)  | (seqb[j] << 8) | (seqa[i] << 4) | IGAP] + Y[i-1][j];
                tmp_tb = 4 << 3;
                t[i][j] |= tmp_tb;
		
                tmp_tb = 4 << 6;
                Y[i][j] =                            prob[           (IGAP << 12) | (seqb[j-1] << 8) | (IGAP << 4) | seqb[j]] + Y[i][j-1];
		
                if(prob[(seqa[i] << 12) | (seqb[j-1] << 8) | (IGAP << 4) | seqb[j]] + M[i][j-1] > Y[i][j]){
                        Y[i][j] = prob[(seqa[i] << 12) | (seqb[j-1] << 8) | (IGAP << 4) | seqb[j]] + M[i][j-1];
                        tmp_tb = 1 << 6;
                }
                if(prob[(seqa[i] << 12) |               (IGAP << 8)  | (IGAP << 4) | seqb[j]] + X[i][j-1] > Y[i][j]){
                        Y[i][j] = prob[(seqa[i] << 12) |               (IGAP << 8)  | (IGAP << 4) | seqb[j]] + X[i][j-1];
                        tmp_tb = 2 << 6;
			
                }
                t[i][j] |= tmp_tb;
		
        }
        M[i][m-1] = prob2scaledprob(0.0f);
        X[i][m-1] = prob2scaledprob(0.0f);
        Y[i][m-1] = prob2scaledprob(0.0f);
        t[i][m-1] = 0;
	
        M[i][m] = prob2scaledprob(0.0f);
        X[i][m] = prob2scaledprob(0.0f);
        Y[i][m] = prob2scaledprob(0.0f);
        t[i][m] = 0;
	
	
        for(i = 2; i < n;i++){
                M[i][0] = prob2scaledprob(0.0f);
                X[i][0] = prob2scaledprob(0.0f);
                Y[i][0] = prob2scaledprob(0.0f);
                t[i][0] = 0;
                M[i][1] = prob2scaledprob(0.0f);
                X[i][1] = prob2scaledprob(0.0f);
                Y[i][1] = prob2scaledprob(0.0f);
		
                t[i][1] = 0;
		
		
                for(j = 2; j <= m-2;j++){
			
                        tmp_tb = 1;
                        M[i][j] =                             prob[(seqa[i-1] << 12) | (seqb[j-1] << 8) | (seqa[i] << 4) | (seqb[j])] + M[i-1][j-1];
			
                        if(prob[(seqa[i-1] << 12) |               (IGAP << 8) | (seqa[i] << 4) | (seqb[j])] + X[i-1][j-1] > M[i][j]){
                                M[i][j] = prob[(seqa[i-1] << 12) |               (IGAP << 8) | (seqa[i] << 4) | (seqb[j])] + X[i-1][j-1];
                                tmp_tb = 2;
                        }
			
                        if(prob[              (IGAP << 12) | (seqb[j-1] << 8) | (seqa[i] << 4) | (seqb[j])] + Y[i-1][j-1] > M[i][j]){
                                M[i][j] = prob[              (IGAP << 12) | (seqb[j-1] << 8) | (seqa[i] << 4) | (seqb[j])] + Y[i-1][j-1];
                                tmp_tb = 4;
                        }
                        t[i][j] = tmp_tb;
			
			
                        tmp_tb = 2 << 3;
                        X[i][j] =                            prob[(seqa[i-1] << 12) |            (IGAP << 8) | (seqa[i] << 4) | IGAP] + X[i-1][j] ;
                        if(prob[(seqa[i-1] << 12) | (seqb[j] << 8) | (seqa[i] << 4) | IGAP] + M[i-1][j] > X[i][j]){
                                X[i][j] = prob[(seqa[i-1] << 12) | (seqb[j] << 8) | (seqa[i] << 4) | IGAP] + M[i-1][j];
                                tmp_tb = 1 << 3;
                        }
                        if(prob[              (IGAP << 12)  | (seqb[j] << 8) | (seqa[i] << 4) | IGAP] + Y[i-1][j] > X[i][j]){
                                X[i][j] = prob[              (IGAP << 12)  | (seqb[j] << 8) | (seqa[i] << 4) | IGAP] + Y[i-1][j];
                                tmp_tb = 4 << 3;
                        }
                        t[i][j] |= tmp_tb;
			
                        tmp_tb = 4 << 6;
                        Y[i][j] =                            prob[           (IGAP << 12) | (seqb[j-1] << 8) | (IGAP << 4) | seqb[j]] + Y[i][j-1];
			
                        if(prob[(seqa[i] << 12) | (seqb[j-1] << 8) | (IGAP << 4) | seqb[j]] + M[i][j-1] > Y[i][j]){
                                Y[i][j] = prob[(seqa[i] << 12) | (seqb[j-1] << 8) | (IGAP << 4) | seqb[j]] + M[i][j-1];
                                tmp_tb = 1 << 6;
                        }
                        if(prob[(seqa[i] << 12) |               (IGAP << 8)  | (IGAP << 4) | seqb[j]] + X[i][j-1] > Y[i][j]){
                                Y[i][j] = prob[(seqa[i] << 12) |               (IGAP << 8)  | (IGAP << 4) | seqb[j]] + X[i][j-1];
                                tmp_tb = 2 << 6;
				
                        }
                        t[i][j] |= tmp_tb;
			
                }
                M[i][m-1] = prob2scaledprob(0.0f);
                X[i][m-1] = prob2scaledprob(0.0f);
                Y[i][m-1] = prob2scaledprob(0.0f);
                t[i][m-1] = 0;
		
                M[i][m] = prob2scaledprob(0.0f);
                X[i][m] = prob2scaledprob(0.0f);
                Y[i][m] = prob2scaledprob(0.0f);
                t[i][m] = 0;
        }
	
	
        M[n][0] = prob2scaledprob(0.0f);
        X[n][0] = prob2scaledprob(0.0f);
        Y[n][0] = prob2scaledprob(0.0f);
        t[n][0] = 0;
	
        M[n][1] = prob2scaledprob(0.0f);
        X[n][1] = prob2scaledprob(0.0f);
        Y[n][1] = prob2scaledprob(0.0f);
        t[n][1] = 0;
	
        for(j = 2; j <= m-2;j++){
                tmp_tb = 1;
                M[n][j] =                             prob[(seqa[n-1] << 12) | (seqb[j-1] << 8) | (seqa[n] << 4) | (seqb[j])] + M[n-1][j-1];
		
                if(prob[(seqa[n-1] << 12) |               (IGAP << 8) | (seqa[n] << 4) | (seqb[j])] + X[n-1][j-1] > M [n][j]){
                        M[n][j] = prob[(seqa[n-1] << 12) |               (IGAP << 8) | (seqa[n] << 4) | (seqb[j])] + X[n-1][j-1];
                        tmp_tb  = 2;
                }
		
                if(prob[              (IGAP << 12) | (seqb[j-1] << 8) | (seqa[i] << 4) | (seqb[j])] + Y[n-1][j-1] > M[n][j]){
                        M[n][j] = prob[              (IGAP << 12) | (seqb[j-1] << 8) | (seqa[i] << 4) | (seqb[j])] + Y[n-1][j-1];
                        tmp_tb = 4;
                }
                t[n][j] = tmp_tb;
		
                tmp_tb = 2 << 3;
                X[n][j] =                            prob[(seqa[n-1] << 12) |            (IGAP << 8) | (seqa[i] << 4) | IGAP] + X[n-1][j] ;
                if(prob[(seqa[n-1] << 12) | (seqb[j] << 8) | (seqa[i] << 4) | IGAP] + M[n-1][j] > X[n][j]){
                        X[n][j] = prob[(seqa[n-1] << 12) | (seqb[j] << 8) | (seqa[i] << 4) | IGAP] + M[n-1][j];
                        tmp_tb = 1 << 3;
                }
                if(prob[              (IGAP << 12)  | (seqb[j] << 8) | (seqa[i] << 4) | IGAP] + Y[n-1][j] > X[n][j]){
                        X[n][j] = prob[              (IGAP << 12)  | (seqb[j] << 8) | (seqa[i] << 4) | IGAP] + Y[n-1][j];
                        tmp_tb = 4 << 3;
                }
                t[n][j] |= tmp_tb;
		
                tmp_tb = 4 << 6;
                Y[n][j] =                            prob[           (EGAP << 12) | (seqb[j-1] << 8) | (EGAP << 4) | seqb[j]] + Y[n][j-1];
                if(prob[(seqa[n] << 12) | (seqb[j-1] << 8) | (EGAP << 4) | seqb[j]] + M[n][j-1] > Y[n][j]){
                        Y[n][j] = prob[(seqa[n] << 12) | (seqb[j-1] << 8) | (EGAP << 4) | seqb[j]] + M[n][j-1];
                        tmp_tb = 1 << 6;
                }
                if(prob[(seqa[n] << 12) |               (IGAP << 8)  | (EGAP << 4) | seqb[j]] + X[n][j-1] > Y[n][j]){
                        Y[n][j] = prob[(seqa[n] << 12) |               (IGAP << 8)  | (EGAP << 4) | seqb[j]] + X[n][j-1];
                        tmp_tb = 2 << 6;
			
                }
                t[n][j] |= tmp_tb;
		
        }
        M[n][m-1] = prob2scaledprob(0.0f);
        X[n][m-1] = prob2scaledprob(0.0f);
	
        tmp_tb = 4 << 6;
        Y[n][m-1] =                            prob[           (EGAP << 12) | (seqb[m-2] << 8) | (EGAP << 4) | seqb[m-1]] + Y[n][m-2];
	
        if(prob[(seqa[n] << 12) | (seqb[m-2] << 8) | (EGAP << 4) | seqb[m-1]] + M[n][m-2] > Y[n][m-1]){
                Y[n][m-1] = prob[(seqa[n] << 12) | (seqb[m-2] << 8) | (EGAP << 4) | seqb[m-1]] + M[n][m-2];
                tmp_tb = 1 << 6;
        }
	
        if(prob[(seqa[n] << 12) |               (IGAP << 8)  | (EGAP << 4) | seqb[m-1]] + X[n][m-2] > Y[n][m-1]){
                Y[n][m-1] = prob[(seqa[n] << 12) |               (IGAP << 8)  | (EGAP << 4) | seqb[m-1]] + X[n][m-2];
                tmp_tb = 2 << 6;
		
        }
	
        t[n][m-1] = tmp_tb;
	
        M[n][m] = prob2scaledprob(0.0f);
        X[n][m] = prob2scaledprob(0.0f);
	
        Y[n][m] =                            prob[           (EGAP << 12) | (seqb[m-1] << 8) | (EGAP << 4) | seqb[m]] + Y[n][m-1];
        t[n][m] = 4 << 6;
        i = n;
        j = m;
        if(M[i][j] >= Y[i][j]){
                if(M[i][j] >= X[i][j]){
                        tmp_tb = 1;
                }else{
                        tmp_tb = 2;
                }
        }else{
                if(Y[i][j] > X[i][j]){
                        tmp_tb = 4;
                }else{
                        tmp_tb = 2;
                }
		
        }
        tmp_tb = 4;
        float id = 0.0f;
        c = 1;
	
        while(1){
                switch (tmp_tb) {
                case 1:
                        if(seqa[i] == NCODE){
                                if(seqb[j] == NCODE){
                                        alignment[c] = (4<< 4) |4;
                                }else{
                                        alignment[c] = (4<< 4) | (seqb[j]);
                                }
                        }else {
                                if(seqb[j] == NCODE){
                                        alignment[c] = ((seqa[i] & 3)<< 4) | 4;
                                }else{
                                        alignment[c] = ((seqa[i] & 3)<< 4) | (seqb[j]);
                                }
                        }
				
                        if(seqa[i] == seqb[j]){
                                id += 1.0f;
                        }
				
                        tmp_tb = (t[i][j]  ) & 0x7;
                        i--;
                        j--;
                        break;
                case 2:
                        if(seqa[i] == NCODE){
                                alignment[c] = (4<< 4) | 5;
                        }else {
                                alignment[c] = ((seqa[i] & 3)<< 4) | 5;
                        }
				
                        tmp_tb = (t[i][j] >> 3 ) & 0x7;
                        i--;
                        break;
                case 4:
                        if(seqb[j] == NCODE){
                                alignment[c] = (5 << 4) |4 ;
                        }else{
                                alignment[c] = (5 << 4) | (seqb[j]) ;
                        }
                        tmp_tb = (t[i][j] >> 6 ) & 0x7;
                        j--;
                        break;
				
                default:
				
                        fprintf(stderr,"ERROR:code: %d not found\n",tmp_tb);
                        for(i = 0; i < n;i++){
                                fprintf(stderr,"%c","ACGTN"[(int)a[i]]);
                        }
                        fprintf(stderr,"\n");
                        for(i = 0; i < m;i++){
                                fprintf(stderr,"%c","ACGTN"[(int)b[i]]);
                        }
                        fprintf(stderr,"\n");
                        return NULL;
                        //exit(0);
                        break;
                }
                c++;
                if(i == 0 && j == 0){
                        break;
                }
        }
        alignment[0] = c-1;
	
        return alignment;
ERROR:
        return NULL;
}

int glocal_forward_log_Y(struct hmm* hmm, char* a,  char* b, int n,int m,int frame )
{
        //init - terminal gap penalty = 0 (p = 1)
        ASSERT(n > 10, "read too short ");
        ASSERT(m > 10, " genome seq too short");
	
        int i,j;
        float** M = 0;
        float** X = 0;
        float** Y = 0;
	
        char* seqa = a -1;
        char* seqb = b -1;
        float* prob = hmm->prob;
	
        for(j = 0; j < n;j++){
                ASSERT(a[j] < 5,"problem in a");
        }
	
        for(j = 0; j < m;j++){
                ASSERT(b[j] < 5,"problem in b: %d",b[j]  );
        }
	
        if(frame == -1){
                M = hmm->fM;
                //float** C = hmm->fC;
                X = hmm->fX;
                Y = hmm->fY;
        }else{
                M = hmm->tfM[frame];
                //float** C = hmm->fC;
                X = hmm->tfX[frame];
                Y = hmm->tfY[frame];
        }
        //c = 0;
        M[0][0] = prob2scaledprob(0.0f);
        X[0][0] = prob2scaledprob(0.0f);
        Y[0][0] = prob2scaledprob(0.0f);
	
        M[0][1] = prob2scaledprob(0.0f);
        X[0][1] = prob2scaledprob(0.0f);
        Y[0][1] = prob2scaledprob(0.0f);
	
        M[0][2] = prob2scaledprob(0.0f);
        X[0][2] = prob2scaledprob(0.0f);
        Y[0][2] = prob[   (EGAP << 12) | (seqb[1] << 8) | (EGAP << 4) | seqb[2]];
        //fprintf(stderr,"%d start %d\n",i,c);
        for(j = 3; j <= m-2;j++){
                M[0][j] = prob2scaledprob(0.0f);
                X[0][j] = prob2scaledprob(0.0f);
                Y[0][j] = prob[   (EGAP << 12) | (seqb[j-1] << 8) | (EGAP << 4) | seqb[j]] + Y[0][j-1];
        }
	
	
        M[0][m-1] = prob2scaledprob(0.0f);
        X[0][m-1] = prob2scaledprob(0.0f);
        Y[0][m-1] = prob2scaledprob(0.0f);
	
        M[0][m] = prob2scaledprob(0.0f);
        X[0][m] = prob2scaledprob(0.0f);
        Y[0][m] = prob2scaledprob(0.0f);
	
        i = 1;
        M[i][0] = prob2scaledprob(0.0f);
        X[i][0] = prob2scaledprob(0.0f);
        Y[i][0] = prob2scaledprob(0.0f);
	
        M[i][1] = prob2scaledprob(0.0f);
        X[i][1] = prob2scaledprob(0.0f);
        Y[i][1] = prob2scaledprob(0.0f);
	
	
        for(j = 2; j <= m-2;j++){
                //M[i][j] =                             prob[(seqa[i-1] << 12) | (seqb[j-1] << 8) | (seqa[i] << 4) | (seqb[j])] + M[i-1][j-1];
                //M[i][j] = logsum(M[i][j], prob[(seqa[i-1] << 12) |               (IGAP << 8) | (seqa[i] << 4) | (seqb[j])] + X[i-1][j-1]);
                M[i][j] = prob[              (EGAP << 12) | (seqb[j-1] << 8) | (seqa[i] << 4) | (seqb[j])] + Y[i-1][j-1];
		
                //X[i][j] =                            prob[(seqa[i-1] << 12) |            (IGAP << 8) | (seqa[i] << 4) | IGAP] + X[i-1][j] ;
                //X[i][j] = logsum(X[i][j], prob[(seqa[i-1] << 12) | (seqb[j] << 8) | (seqa[i] << 4) | IGAP] + M[i-1][j]);
                X[i][j] = prob[              (EGAP << 12)  | (seqb[j] << 8) | (seqa[i] << 4) | IGAP] + Y[i-1][j];
		
                Y[i][j] =                            prob[           (IGAP << 12) | (seqb[j-1] << 8) | (IGAP << 4) | seqb[j]] + Y[i][j-1];
		
                //fprintf(stdout,"seqa[i]:%d seqb[j-1]:%d seqb[j]:%d	prob:%f\n",seqa[i],seqb[j-1],seqb[j],prob[(seqa[i] << 12) | (seqb[j-1] << 8) | (IGAP << 4) | seqb[j]] + M[i][j-1]);
		
                Y[i][j] = logsum(Y[i][j] ,prob[(seqa[i] << 12) | (seqb[j-1] << 8) | (IGAP << 4) | seqb[j]] + M[i][j-1]);
                Y[i][j] = logsum(Y[i][j] ,prob[(seqa[i] << 12) |               (IGAP << 8)  | (IGAP << 4) | seqb[j]] + X[i][j-1]);
		
        }
        M[i][m-1] = prob2scaledprob(0.0f);
        X[i][m-1] = prob2scaledprob(0.0f);
        Y[i][m-1] = prob2scaledprob(0.0f);
	
        M[i][m] = prob2scaledprob(0.0f);
        X[i][m] = prob2scaledprob(0.0f);
        //X[i][m] = prob[              (IGAP << 12)  | (seqb[m] << 8) | (seqa[i] << 4) | IGAP] + Y[i-1][m];
        Y[i][m] = prob2scaledprob(0.0f);
	
	
        for(i = 2; i < n;i++){//= 4;i++){
		
		
                M[i][0] = prob2scaledprob(0.0f);
                X[i][0] = prob2scaledprob(0.0f);
                Y[i][0] = prob2scaledprob(0.0f);
		
                M[i][1] = prob2scaledprob(0.0f);
                X[i][1] = prob2scaledprob(0.0f);
                Y[i][1] = prob2scaledprob(0.0f);
		
		
                for(j = 2; j <= m-2;j++){
			
                        M[i][j] =                             prob[(seqa[i-1] << 12) | (seqb[j-1] << 8) | (seqa[i] << 4) | (seqb[j])] + M[i-1][j-1];
                        M[i][j] = logsum(M[i][j], prob[(seqa[i-1] << 12) |               (IGAP << 8) | (seqa[i] << 4) | (seqb[j])] + X[i-1][j-1]);
                        M[i][j] = logsum(M[i][j], prob[              (IGAP << 12) | (seqb[j-1] << 8) | (seqa[i] << 4) | (seqb[j])] + Y[i-1][j-1]);
			
                        X[i][j] =                            prob[(seqa[i-1] << 12) |            (IGAP << 8) | (seqa[i] << 4) | IGAP] + X[i-1][j] ;
                        X[i][j] = logsum(X[i][j], prob[(seqa[i-1] << 12) | (seqb[j] << 8) | (seqa[i] << 4) | IGAP] + M[i-1][j]);
                        X[i][j] = logsum(X[i][j], prob[              (IGAP << 12)  | (seqb[j] << 8) | (seqa[i] << 4) | IGAP] + Y[i-1][j]);
			
                        Y[i][j] =                            prob[           (IGAP << 12) | (seqb[j-1] << 8) | (IGAP << 4) | seqb[j]] + Y[i][j-1];
                        Y[i][j] = logsum(Y[i][j] ,prob[(seqa[i] << 12) | (seqb[j-1] << 8) | (IGAP << 4) | seqb[j]] + M[i][j-1]);
                        Y[i][j] = logsum(Y[i][j] ,prob[(seqa[i] << 12) |               (IGAP << 8)  | (IGAP << 4) | seqb[j]] + X[i][j-1]);
			
                }
                M[i][m-1] = prob2scaledprob(0.0f);
                X[i][m-1] = prob2scaledprob(0.0f);
                Y[i][m-1] = prob2scaledprob(0.0f);
		
                M[i][m] = prob2scaledprob(0.0f);
                X[i][m] = prob2scaledprob(0.0f);
                //X[i][m] =                            prob[(seqa[i-1] << 12) |            (IGAP << 8) | (seqa[i] << 4) | IGAP] + X[i-1][m] ;
                Y[i][m] = prob2scaledprob(0.0f);
        }
	
	
        M[n][0] = prob2scaledprob(0.0f);
        X[n][0] = prob2scaledprob(0.0f);
        Y[n][0] = prob2scaledprob(0.0f);
	
        M[n][1] = prob2scaledprob(0.0f);
        X[n][1] = prob2scaledprob(0.0f);
        Y[n][1] = prob2scaledprob(0.0f);
	
        for(j = 2; j <= m-2;j++){
                M[n][j] =                             prob[(seqa[n-1] << 12) | (seqb[j-1] << 8) | (seqa[n] << 4) | (seqb[j])] + M[n-1][j-1];
                M[n][j] = logsum(M[n][j], prob[(seqa[n-1] << 12) |               (IGAP << 8) | (seqa[n] << 4) | (seqb[j])] + X[n-1][j-1]);
                M[n][j] = logsum(M[n][j], prob[              (IGAP << 12) | (seqb[j-1] << 8) | (seqa[i] << 4) | (seqb[j])] + Y[n-1][j-1]);
		
                X[n][j] =                            prob[(seqa[n-1] << 12) |            (IGAP << 8) | (seqa[i] << 4) | IGAP] + X[n-1][j] ;
                X[n][j] = logsum(X[n][j], prob[(seqa[n-1] << 12) | (seqb[j] << 8) | (seqa[i] << 4) | IGAP] + M[n-1][j]);
                X[n][j] = logsum(X[n][j], prob[              (IGAP << 12)  | (seqb[j] << 8) | (seqa[i] << 4) | IGAP] + Y[n-1][j]);
		
                Y[n][j] =                            prob[           (EGAP << 12) | (seqb[j-1] << 8) | (EGAP << 4) | seqb[j]] + Y[n][j-1];
                Y[n][j] = logsum(Y[n][j] ,prob[(seqa[n] << 12) | (seqb[j-1] << 8) | (EGAP << 4) | seqb[j]] + M[n][j-1]);
                Y[n][j] = logsum(Y[n][j] ,prob[(seqa[n] << 12) |               (IGAP << 8)  | (EGAP << 4) | seqb[j]] + X[n][j-1]);
		
        }
        M[n][m-1] = prob2scaledprob(0.0f);
        X[n][m-1] = prob2scaledprob(0.0f);
        Y[n][m-1] =                            prob[           (EGAP << 12) | (seqb[m-2] << 8) | (EGAP << 4) | seqb[m-1]] + Y[n][m-2];
        Y[n][m-1] = logsum(Y[n][m-1] ,prob[(seqa[n] << 12) | (seqb[m-2] << 8) | (EGAP << 4) | seqb[m-1]] + M[n][m-2]);
        Y[n][m-1] = logsum(Y[n][m-1] ,prob[(seqa[n] << 12) |               (IGAP << 8)  | (EGAP << 4) | seqb[m-1]] + X[n][m-2]);
	
	
        M[n][m] = prob2scaledprob(0.0f);
        X[n][m] = prob2scaledprob(0.0f);
	
        Y[n][m] =                            prob[           (EGAP << 12) | (seqb[m-1] << 8) | (EGAP << 4) | seqb[m]] + Y[n][m-1];
	
        hmm->score =Y[n][m];
        return OK;
ERROR:
        return FAIL;
}

int glocal_backward_log_Y(struct hmm* hmm,char* a, char* b, int n,int m)
{
        int i,j;//c;
	
        float** M = hmm->bM;
        float** X = hmm->bX;
        float** Y = hmm->bY;
        float* prob = hmm->prob;
        char* seqa = a -1;
        char* seqb = b -1;
        //float score;
	
        M[n][m] = prob2scaledprob(0.0f);
        X[n][m] = prob2scaledprob(0.0f);
        Y[n][m] = prob2scaledprob(1.0f);
	
	
        M[n][m-1] = prob2scaledprob(0.0f);
        X[n][m-1] = prob2scaledprob(0.0f);
        Y[n][m-1] = prob[            (EGAP << 12) | (seqb[m-1] << 8) | (EGAP << 4) | seqb[m]] + Y[n][m];
	
	
	
        for(j = m-2; j > 1;j--){
                M[n][j] = prob[(seqa[n] << 12) | (seqb[j] << 8) | (EGAP << 4) | seqb[j+1]] + Y[n][j+1];
                X[n][j] =  prob[(seqa[n] << 12) |            (IGAP << 8) | (EGAP << 4) | seqb[j+1]] + Y[n][j+1];
                Y[n][j] =  prob[            (EGAP << 12) | (seqb[j] << 8) | (EGAP << 4) | seqb[j+1]] + Y[n][j+1];
        }
	
        M[n][1] =  prob2scaledprob(0.0f);
        X[n][1] =  prob2scaledprob(0.0f);
        Y[n][1] =  prob2scaledprob(0.0f);
	
        M[n][0] = prob2scaledprob(0.0f);
        X[n][0] = prob2scaledprob(0.0f);
        Y[n][0] = prob2scaledprob(0.0f);
        //was here ......
	
        for(i = n-1; i >= 1;i--){
                M[i][m] = prob2scaledprob(0.0f);
                X[i][m] = prob2scaledprob(0.0f);
                Y[i][m] = prob2scaledprob(0.0f);
		
		
                M[i][m-1] = prob2scaledprob(0.0f);
                X[i][m-1] = prob2scaledprob(0.0f);
                Y[i][m-1] = prob2scaledprob(0.0f);
		
                for(j = m-2; j > 1;j--){
                        M[i][j] =                             prob[(seqa[i] << 12) | (seqb[j] << 8) | (seqa[i+1] << 4) | seqb[j+1]] + M[i+1][j+1];
                        M[i][j] = logsum(M[i][j], prob[(seqa[i] << 12) | (seqb[j] << 8) |                 (IGAP << 4) | seqb[j+1]] + Y[i][j+1]);
                        M[i][j] = logsum(M[i][j], prob[(seqa[i] << 12) | (seqb[j] << 8) | (seqa[i+1] << 4)  |                IGAP] + X[i+1][j]);
			
                        X[i][j] =                             prob[(seqa[i] << 12) | (IGAP << 8) | (seqa[i+1] << 4) |               IGAP ] + X[i+1][j];
                        X[i][j] =  logsum(X[i][j], prob[(seqa[i] << 12) | (IGAP << 8) | (seqa[i+1] << 4) | seqb[j+1]] + M[i+1][j+1]);
                        X[i][j] =  logsum(X[i][j], prob[(seqa[i] << 12) | (IGAP << 8) |                (IGAP << 4) | seqb[j+1]] + Y[i][j+1]);
			
                        Y[i][j] =                             prob[(IGAP << 12) | (seqb[j] << 8) |                (IGAP << 4) | seqb[j+1]] + Y[i][j+1];
                        Y[i][j] = logsum(Y[i][j], prob[(IGAP << 12) | (seqb[j] << 8) | (seqa[i+1] << 4) | seqb[j+1]] + M[i+1][j+1]);
                        Y[i][j] = logsum(Y[i][j], prob[(IGAP << 12) | (seqb[j] << 8) | (seqa[i+1] << 4) | IGAP  ] + X[i+1][j]);
			
			
			
                }
		
                M[i][1] =  prob2scaledprob(0.0f);
                X[i][1] =  prob2scaledprob(0.0f);
                Y[i][1] =  prob2scaledprob(0.0f);
		
                M[i][0] = prob2scaledprob(0.0f);
                X[i][0] = prob2scaledprob(0.0f);
                Y[i][0] = prob2scaledprob(0.0f);
		
        }
	
	
	
        M[0][m] = prob2scaledprob(0.0f);
        X[0][m] = prob2scaledprob(0.0f);
        Y[0][m] = prob2scaledprob(0.0f);
	
	
        M[0][m-1] = prob2scaledprob(0.0f);
        X[0][m-1] = prob2scaledprob(0.0f);
        Y[0][m-1] = prob2scaledprob(0.0f);
	
	
	
        for(j = m-2; j > 1;j--){
                M[0][j] =  prob2scaledprob(0.0f);
                X[0][j] = prob2scaledprob(0.0f);
                Y[0][j] =                             prob[ (EGAP << 12) | (seqb[j] << 8) |                (EGAP << 4) | seqb[j+1]] + Y[0][j+1];
                Y[0][j] = logsum(Y[0][j], prob[(EGAP << 12) | (seqb[j] << 8) | (seqa[1] << 4) | seqb[j+1]] + M[1][j+1]);
                Y[0][j] = logsum(Y[0][j], prob[(EGAP << 12) | (seqb[j] << 8) | (seqa[1] << 4) | IGAP  ] + X[1][j]);
        }
	
        M[0][1] =  prob2scaledprob(0.0f);
        X[0][1] =  prob2scaledprob(0.0f);
        Y[0][1] =                               prob[(EGAP << 12) | (seqb[1] << 8) |             (EGAP << 4) | seqb[2]] + Y[0][2];
	
        M[0][0] = prob2scaledprob(0.0f);
        X[0][0] = prob2scaledprob(0.0f);
        Y[0][0] = prob[(EGAP << 12) | (seqb[1] << 8) |             (EGAP << 4) | seqb[2]] + Y[0][1];
	
        return OK;
}


int get_prob_log_Y(struct hmm* hmm,char* a, char* b, int n,int m,float frac,int frame )
{
        int i,j,g;
        char* seqa = a -1;
        char* seqb = b -1;
        float* prob = hmm->prob;
        float* tprob = hmm->tprob;
	
        float** M = 0;
        float** X = 0;
        float** Y = 0;
	
        if(frame == -1){
                M = hmm->fM;
                X = hmm->fX;
                Y = hmm->fY;
        }else{
                M = hmm->fM;
                X = hmm->fX;
                Y = hmm->fY;
                hmm->score =  hmm->tfY[frame][n][m];
		
                hmm->fM = hmm->tfM[frame];
                hmm->fX = hmm->tfX[frame];
                hmm->fY = hmm->tfY[frame];
        }
        
        tprob[   (EGAP << 12) | (seqb[1] << 8) | (EGAP << 4) | seqb[2]] = logsum(tprob[   (EGAP << 12) | (seqb[1] << 8) | (EGAP << 4) | seqb[2]], hmm->fY[0][2] + hmm->bY[0][2]  - hmm->score+ frac );
	
        for(j = 3; j <= m-2;j++){
                //M[0][j] = prob2scaledprob(0.0f);
                //X[0][j] = prob2scaledprob(0.0f);
                //Y[0][j] = prob[c][   (6 <<9) | (seqb[j-1] << 6) | (6 << 3) | seqb[j]] + Y[0][j-1];
                g =  (EGAP << 12) | (seqb[j-1] << 8) | (EGAP << 4) | seqb[j];
		
                tprob[g]  = logsum(tprob[g],  prob[g] + hmm->fY[0][j-1] + hmm->bY[0][j]  - hmm->score+ frac );
		
                //Y[i][j] = logsum(Y[i][j], prob[(seqa[i] << 9) | (seqb[j-1] << 6) | (5 << 3) | seqb[j]] );
        }
	
        i = 1;
	
        for(j = 2; j <= m-2;j++){
                //M[i][j] =                             prob[c][(seqa[i-1] << 9) | (seqb[j-1] << 6) | (seqa[i] << 3) | (seqb[j])] + M[i-1][j-1];
                //M[i][j] = logsum(M[i][j], prob[c][(seqa[i-1] << 9) |               (5 << 6) | (seqa[i] << 3) | (seqb[j])] + X[i-1][j-1]);
		
                //M[i][j] = logsum(M[i][j], prob[c][              (5 << 9) | (seqb[j-1] << 6) | (seqa[i] << 3) | (seqb[j])] + Y[i-1][j-1]);
                g =           (EGAP << 12) | (seqb[j-1] << 8) | (seqa[i] << 4) | (seqb[j]);
                tprob[g] = logsum(tprob[g], prob[g] +  hmm->fY[i-1][j-1]+ hmm->bM[i][j]  - hmm->score + frac);
		
                //X[i][j] =                            prob[c][(seqa[i-1] << 9) |            (5 << 6) | (seqa[i] << 3) | 5] + X[i-1][j] ;
		
                //X[i][j] = logsum(X[i][j], prob[c][              (5 << 9)  | (seqb[j] << 6) | (seqa[i] << 3) | 5] + Y[i-1][j]);
                g =               (EGAP << 12)  | (seqb[j] << 8) | (seqa[i] << 4) | IGAP;
                tprob[g] = logsum(tprob[g], prob[g] + hmm->fY[i-1][j] + hmm->bX[i][j] - hmm->score + frac);
		
		
                //Y[i][j] =                            prob[c][           (5 << 9) | (seqb[j-1] << 6) | (5 << 3) | seqb[j]] + Y[i][j-1];
                g =           (IGAP << 12) | (seqb[j-1] << 8) | (IGAP << 4) | seqb[j];
                tprob[g] = logsum(tprob[g], prob[g] + hmm->fY[i][j-1] + hmm->bY[i][j] - hmm->score + frac);
		
                //Y[i][j] = logsum(Y[i][j] ,prob[c][(seqa[i] << 9) | (seqb[j-1] << 6) | (5 << 3) | seqb[j]] + M[i][j-1]);
                g = (seqa[i] << 12) | (seqb[j-1] << 8) | (IGAP << 4) | seqb[j];
                tprob[g] = logsum(tprob[g], prob[g] + hmm->fM[i][j-1] + hmm->bY[i][j] - hmm->score + frac);
		
                //Y[i][j] = logsum(Y[i][j] ,prob[c][(seqa[i] << 9) |               (5 << 6)  | (5 << 3) | seqb[j]] + X[i][j-1]);
                g = (seqa[i] << 12) |               (IGAP << 8)  | (IGAP << 4) | seqb[j];
                tprob[g] = logsum(tprob[g], prob[g] + hmm->fX[i][j-1] + hmm->bY[i][j] - hmm->score + frac);
		
		
        }
	
	
	
        for(i = 2; i < n;i++){//= 4;i++){
                for(j = 2; j <= m-2;j++){
			
                        //fprintf(stderr,"%f	%f	%f\n",hmm->fM[i-1][j-1] ,M[i-1][j-1] , hmm->fM[i-1][j-1] -M[i-1][j-1]);
			
                        //M[i][j] =                             prob[c][(seqa[i-1] << 9) | (seqb[j-1] << 6) | (seqa[i] << 3) | (seqb[j])] + M[i-1][j-1];
                        g = (seqa[i-1] << 12) | (seqb[j-1] << 8) | (seqa[i] << 4) | (seqb[j]);
                        tprob[g] = logsum(tprob[g], prob[g] + hmm->fM[i-1][j-1] + hmm->bM[i][j]  - hmm->score + frac);
                        //M[i][j] = logsum(M[i][j], prob[c][(seqa[i-1] << 9) |               (5 << 6) | (seqa[i] << 3) | (seqb[j])] + X[i-1][j-1]);
                        g = (seqa[i-1] << 12) |               (IGAP << 8) | (seqa[i] << 4) | (seqb[j]);
                        tprob[g] = logsum(tprob[g], prob[g] + hmm->fX[i-1][j-1] + hmm->bM[i][j]  - hmm->score + frac);
                        //M[i][j] = logsum(M[i][j], prob[c][              (5 << 9) | (seqb[j-1] << 6) | (seqa[i] << 3) | (seqb[j])] + Y[i-1][j-1]);
                        g =           (IGAP << 12) | (seqb[j-1] << 8) | (seqa[i] << 4) | (seqb[j]);
                        tprob[g] = logsum(tprob[g], prob[g] +  hmm->fY[i-1][j-1]+ hmm->bM[i][j]  - hmm->score + frac);
			
                        //X[i][j] =                            prob[c][(seqa[i-1] << 9) |            (5 << 6) | (seqa[i] << 3) | 5] + X[i-1][j] ;
                        g = (seqa[i-1] << 12) |            (IGAP << 8) | (seqa[i] << 4) | IGAP;
                        tprob[g] = logsum(tprob[g], prob[g] + hmm->fX[i-1][j] + hmm->bX[i][j] - hmm->score + frac);
                        //X[i][j] = logsum(X[i][j], prob[c][(seqa[i-1] << 9) | (seqb[j] << 6) | (seqa[i] << 3) | 5] + M[i-1][j]);
                        g = (seqa[i-1] << 12) | (seqb[j] << 8) | (seqa[i] << 4) | IGAP;
                        tprob[g] = logsum(tprob[g], prob[g] + hmm->fM[i-1][j] + hmm->bX[i][j] - hmm->score + frac);
                        //X[i][j] = logsum(X[i][j], prob[c][              (5 << 9)  | (seqb[j] << 6) | (seqa[i] << 3) | 5] + Y[i-1][j]);
                        g =               (IGAP << 12)  | (seqb[j] << 8) | (seqa[i] << 4) | IGAP;
                        tprob[g] = logsum(tprob[g], prob[g] + hmm->fY[i-1][j] + hmm->bX[i][j] - hmm->score + frac);
			
			
                        //Y[i][j] =                            prob[c][           (5 << 9) | (seqb[j-1] << 6) | (5 << 3) | seqb[j]] + Y[i][j-1];
                        g =           (IGAP << 12) | (seqb[j-1] << 8) | (IGAP << 4) | seqb[j];
                        tprob[g] = logsum(tprob[g], prob[g] + hmm->fY[i][j-1] + hmm->bY[i][j] - hmm->score + frac);
			
                        //Y[i][j] = logsum(Y[i][j] ,prob[c][(seqa[i] << 9) | (seqb[j-1] << 6) | (5 << 3) | seqb[j]] + M[i][j-1]);
                        g = (seqa[i] << 12) | (seqb[j-1] << 8) | (IGAP << 4) | seqb[j];
                        tprob[g] = logsum(tprob[g], prob[g] + hmm->fM[i][j-1] + hmm->bY[i][j] - hmm->score + frac);
			
                        //Y[i][j] = logsum(Y[i][j] ,prob[c][(seqa[i] << 9) |               (5 << 6)  | (5 << 3) | seqb[j]] + X[i][j-1]);
                        g = (seqa[i] << 12) |               (IGAP << 8)  | (IGAP << 4) | seqb[j];
                        tprob[g] = logsum(tprob[g], prob[g] + hmm->fX[i][j-1] + hmm->bY[i][j] - hmm->score + frac);
			
			
                }
        }
	
        for(j = 2; j <= m-2;j++){
                //M[n][j] =                             prob[c][(seqa[n-1] << 9) | (seqb[j-1] << 6) | (seqa[n] << 3) | (seqb[j])] + M[n-1][j-1];
                g = (seqa[n-1] << 12) | (seqb[j-1] << 8) | (seqa[n] << 4) | (seqb[j]);
                tprob[g] = logsum(tprob[g], prob[g] + hmm->fM[n-1][j-1] + hmm->bM[n][j]  - hmm->score + frac);
		
                //M[n][j] = logsum(M[n][j], prob[c][(seqa[n-1] << 9) |               (5 << 6) | (seqa[n] << 3) | (seqb[j])] + X[n-1][j-1]);
                g = (seqa[n-1] << 12) |               (IGAP << 8) | (seqa[n] << 4) | (seqb[j]);
                tprob[g] = logsum(tprob[g], prob[g] + hmm->fX[n-1][j-1] + hmm->bM[n][j]  - hmm->score + frac);
		
                //M[n][j] = logsum(M[n][j], prob[c][              (5 << 9) | (seqb[j-1] << 6) | (seqa[i] << 3) | (seqb[j])] + Y[n-1][j-1]);
                g =               (IGAP << 12) | (seqb[j-1] << 8) | (seqa[n] << 4) | (seqb[j]);
                tprob[g] = logsum(tprob[g], prob[g] + hmm->fY[n-1][j-1] + hmm->bM[n][j]  - hmm->score + frac);
		
		
                //X[n][j] =                            prob[c][(seqa[n-1] << 9) |            (5 << 6) | (seqa[i] << 3) | 5] + X[n-1][j] ;
                g = (seqa[n-1] << 12) |            (IGAP << 8) | (seqa[n] << 4) | IGAP;
                tprob[g] = logsum(tprob[g], prob[g] + hmm->fX[n-1][j] + hmm->bX[n][j] - hmm->score + frac);
		
                //X[n][j] = logsum(X[n][j], prob[c][(seqa[n-1] << 9) | (seqb[j] << 6) | (seqa[i] << 3) | 5] + M[n-1][j]);
                g = (seqa[n-1] << 12) | (seqb[j] << 8) | (seqa[n] << 4) | IGAP;
                tprob[g] = logsum(tprob[g], prob[g] + hmm->fM[n-1][j] + hmm->bX[n][j] - hmm->score + frac);
		
                //X[n][j] = logsum(X[n][j], prob[c][              (5 << 9)  | (seqb[j] << 6) | (seqa[i] << 3) | 5] + Y[n-1][j]);
                g =                (IGAP << 12)  | (seqb[j] << 8) | (seqa[n] << 4) | IGAP;
                tprob[g] = logsum(tprob[g], prob[g] + hmm->fY[n-1][j] + hmm->bX[n][j] - hmm->score + frac);
		
                //Y[n][j] =                            prob[c][           (6 << 9) | (seqb[j-1] << 6) | (6 << 3) | seqb[j]] + Y[n][j-1];
                g =            (EGAP << 12) | (seqb[j-1] << 8) | (EGAP << 4) | seqb[j];
                tprob[g] = logsum(tprob[g], prob[g] + hmm->fY[n][j-1] + hmm->bY[n][j] - hmm->score+ frac );
		
                //Y[n][j] = logsum(Y[n][j] ,prob[c][(seqa[n] << 9) | (seqb[j-1] << 6) | (5 << 3) | seqb[j]] + M[n][j-1]);
                g = (seqa[n] << 12) | (seqb[j-1] << 8) | (EGAP << 4) | seqb[j];
                tprob[g] = logsum(tprob[g], prob[g] + hmm->fM[n][j-1] + hmm->bY[n][j] - hmm->score + frac);
		
                //Y[n][j] = logsum(Y[n][j] ,prob[c][(seqa[n] << 9) |               (5 << 6)  | (5 << 3) | seqb[j]] + X[n][j-1]);
                g = (seqa[n] << 12) |               (IGAP << 8)  | (EGAP << 4) | seqb[j];
                tprob[g] = logsum(tprob[g], prob[g] + hmm->fX[n][j-1] + hmm->bY[n][j] - hmm->score + frac);
		
        }
	
        g =           (EGAP << 12) | (seqb[m-2] << 8) | (EGAP << 4) | seqb[m-1];
        tprob[g] = logsum(tprob[g], prob[g] + hmm->fY[n][m-2] + hmm->bY[n][m-1]  - hmm->score+ frac );
	
        //Y[n][m-1] = logsum(Y[n][m-1] ,prob[c][(seqa[n] << 9) | (seqb[m-2] << 6) | (5 << 3) | seqb[m-1]] + M[n][m-2]);
        g = (seqa[n] << 12) | (seqb[m-2] << 8) | (EGAP << 4) | seqb[m-1];
        tprob[g] = logsum(tprob[g], prob[g] + hmm->fM[n][m-2] + hmm->bY[n][m-1]  - hmm->score + frac);
	
        //Y[n][m-1] = logsum(Y[n][m-1] ,prob[c][(seqa[n] << 9) |               (5 << 6)  | (5 << 3) | seqb[m-1]] + X[n][m-2]);
        g = (seqa[n] << 12) |               (IGAP << 8)  | (EGAP << 4) | seqb[m-1];
        tprob[g] = logsum(tprob[g], prob[g] + hmm->fX[n][m-2] + hmm->bY[n][m-1]  - hmm->score + frac);
	
        g =           (EGAP << 12) | (seqb[m-1] << 8) | (EGAP << 4) | seqb[m];
        tprob[g] = logsum(tprob[g], prob[g] + hmm->fY[n][m-1] + hmm->bY[n][m] - hmm->score+ frac );
	
	
        if(frame != -1){
                hmm->fM = M;// hmm->tfM[frame];
                hmm->fX = X;//hmm->tfX[frame];
                hmm->fY = Y;//hmm->tfY[frame];
        }
	
        return OK;
}





int random_model_read_score(struct hmm* hmm, char* seq,int len,float frac)
{
        int i,code;

        float score = prob2scaledprob(1.0);

        for(i = 0; i < len-1;i++){
                if(seq[i] == NCODE){
                        code = 20;
                }else{
                        code = (seq[i] & 0x3) *5;
                }
                if(seq[i+1] == NCODE){
                        code += 4;
                }else{
                        code+= (seq[i+1] & 0x3);
                }
                score += hmm->random_read_model[code];
                //hmm->random_score += hmm->random_read_model[code];
                hmm->random_read_model_t[code] =  logsum(hmm->random_read_model_t[code], frac );
        }

        hmm->unaligned_read_score = score;
	
        return OK;
}


int random_model_genome_score(struct hmm* hmm, char* seq,int len,float frac)
{
        int i,code;

        float score = prob2scaledprob(1.0);
        for(i = 0; i < len-1;i++){
                if(seq[i] == NCODE){
                        code = 20;
				
                }else{
                        code = (seq[i] & 0x3) *5;
                }
                if(seq[i+1] == NCODE){
                        code += 4;
                }else{
                        code+= (seq[i+1] & 0x3);
                }
                score += hmm->random_genome_model[code];
                hmm->random_genome_model_t[code] =  logsum(hmm->random_genome_model_t[code],  frac);

        }
	

        hmm->unaligned_genome_score = score;
	
        return OK;
}
