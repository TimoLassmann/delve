
#include "tldevel.h"

#include "delve_struct.h"
#include "hmm.h"

struct hmm* init_hmm(int x, int y, int z)
{
	struct hmm* hmm = NULL;
	int i,j,c;
	int r1,r2,g1,g2,key;
	float* sum_prob = 0;
	float sum;
	x = x + 2;
	y = y + 2;
	
	MMALLOC(hmm, sizeof(struct hmm));


	hmm->x = x;
	hmm->y = y;
	hmm->num_hits = z;
	
	hmm->fM = NULL;// malloc(sizeof(float*)*x);
	hmm->fX = NULL;//malloc(sizeof(float*)*x);
	hmm->fY = NULL;//malloc(sizeof(float*)*x);
	hmm->bM = NULL;//= malloc(sizeof(float*)*x);
	hmm->bX = NULL;//= malloc(sizeof(float*)*x);
	hmm->bY = NULL;//= malloc(sizeof(float*)*x);
	
	hmm->traceback = NULL;
	hmm->tprob = NULL;
	hmm->prob = NULL;
	hmm->print_prob = NULL;

	hmm->unaligned_scores = NULL;
	hmm->alignment_scores = NULL;
	hmm->tmp_scores = NULL;

	MMALLOC(hmm->unaligned_scores,sizeof(float) * (z+1));
	MMALLOC(hmm->alignment_scores,sizeof(float) * (z+1));
	MMALLOC(hmm->tmp_scores,sizeof(float) * (z+1));
	
	
	
	MMALLOC(hmm->fM, sizeof(float*)*x);
	MMALLOC(hmm->fX, sizeof(float*)*x);
	MMALLOC(hmm->fY, sizeof(float*)*x);
	
	MMALLOC(hmm->bM, sizeof(float*)*x);
	MMALLOC(hmm->bX, sizeof(float*)*x);
	MMALLOC(hmm->bY, sizeof(float*)*x);
	
	//hmm->fM = malloc(sizeof(float*)*x);
	//hmm->fX = malloc(sizeof(float*)*x);
	//hmm->fY = malloc(sizeof(float*)*x);
	hmm->tfM = NULL;//malloc(sizeof(float**)*LIST_STORE_SIZE);
	hmm->tfX = NULL;//malloc(sizeof(float**)*LIST_STORE_SIZE);
	hmm->tfY = NULL;//malloc(sizeof(float**)*LIST_STORE_SIZE);
	
	MMALLOC(hmm->tfM, sizeof(float**)* hmm->num_hits);
	MMALLOC(hmm->tfX, sizeof(float**)* hmm->num_hits);
	MMALLOC(hmm->tfY, sizeof(float**)* hmm->num_hits);
	
	
	//hmm->tfM = malloc(sizeof(float**)*LIST_STORE_SIZE);
	//hmm->tfX = malloc(sizeof(float**)*LIST_STORE_SIZE);
	//hmm->tfY = malloc(sizeof(float**)*LIST_STORE_SIZE);
	
	MMALLOC(hmm->traceback , sizeof(unsigned short int*)*x);
	
	//hmm->traceback = malloc(sizeof(unsigned short int*)*x);
	
	//hmm->prob = malloc(sizeof(float) * 65536);
	//hmm->tprob = malloc(sizeof(float) * 65536);
	//hmm->print_prob = malloc(sizeof(float) * 65536);
	MMALLOC(hmm->prob, sizeof(float) * 65536);
	MMALLOC(hmm->tprob, sizeof(float) * 65536);
	MMALLOC(hmm->print_prob, sizeof(float) * 65536);

	for(i = 0; i < 65536;i++){
		hmm->prob[i] = 0.0;
		hmm->tprob[i] = 0.0;
		hmm->print_prob[i] = 0.0;
	}
	
	
	for(i = 0; i < hmm->num_hits;i++){
		hmm->tfM[i] = NULL;
		hmm->tfX[i] = NULL;
		hmm->tfY[i] = NULL;
		MMALLOC(hmm->tfM[i], sizeof(float*)*x);
		MMALLOC(hmm->tfX[i], sizeof(float*)*x);
		MMALLOC(hmm->tfY[i], sizeof(float*)*x);
		
		for(j = 0;j < x;j++){
			hmm->tfM[i][j] = NULL;
			hmm->tfX[i][j] = NULL;
			hmm->tfY[i][j] = NULL;
			MMALLOC(hmm->tfM[i][j],sizeof(float*)*y);
			MMALLOC(hmm->tfX[i][j],sizeof(float*)*y);
			MMALLOC(hmm->tfY[i][j],sizeof(float*)*y);
			
			for(c = 0; c < y;c++){
				hmm->tfM[i][j][c] = 0;
				hmm->tfX[i][j][c] = 0;
				hmm->tfY[i][j][c] = 0;
			}
		}
	}
	
	for(i = 0;i < x;i++){
		
		hmm->fM[i] = NULL;//malloc(sizeof(float*)*y);
		hmm->fX[i] = NULL;//malloc(sizeof(float*)*y);
		hmm->fY[i] = NULL;//malloc(sizeof(float*)*y);
		
		hmm->bM[i] = NULL;//malloc(sizeof(float*)*y);
		hmm->bX[i] = NULL;//malloc(sizeof(float*)*y);
		hmm->bY[i] = NULL;//malloc(sizeof(float*)*y);
		
		hmm->traceback[i] = NULL;//malloc(sizeof(unsigned short int*)*y);
		
		
		MMALLOC(hmm->fM[i],sizeof(float*)*y);
		MMALLOC(hmm->fX[i],sizeof(float*)*y);
		MMALLOC(hmm->fY[i],sizeof(float*)*y);
		
		MMALLOC(hmm->bM[i],sizeof(float*)*y);
		MMALLOC(hmm->bX[i],sizeof(float*)*y);
		MMALLOC(hmm->bY[i],sizeof(float*)*y);
		
		MMALLOC(hmm->traceback[i],sizeof(unsigned short int*)*y);
		for(j = 0; j < y;j++){
			hmm->fM[i][j] = 0;
			hmm->fX[i][j] = 0;
			hmm->fY[i][j] = 0;
			
			hmm->bM[i][j] = 0;
			hmm->bX[i][j] = 0;
			hmm->bY[i][j] = 0;
			hmm->traceback[i][j] = 0;
		}
	}
	sum_prob = NULL;
	MMALLOC(sum_prob,sizeof(float) * 256 );
	
	hmm->random_genome_model = NULL;//malloc(sizeof(float)* 25);
	hmm->random_read_model = NULL;//malloc(sizeof(float)* 25);
	hmm->random_genome_model_t = NULL;//malloc(sizeof(float)* 25);
	hmm->random_read_model_t = NULL;//malloc(sizeof(float)* 25);
	
	MMALLOC(hmm->random_genome_model,sizeof(float)* 25);
	MMALLOC(hmm->random_read_model,sizeof(float)* 25);
	MMALLOC(hmm->random_genome_model_t,sizeof(float)* 25);
	MMALLOC(hmm->random_read_model_t,sizeof(float)* 25);
	
	for(i = 0; i < 256;i++){
		sum_prob[i] = 0;
	}
	r1 = -1;
	for(j = 0; j < 25;j++){
		if(j % 5 == 0){
			r1++;
		}
		hmm->random_genome_model[j] = prob2scaledprob(0.00001f);
		hmm->random_read_model[j] = prob2scaledprob(0.00001f);
		hmm->random_genome_model_t[j] = prob2scaledprob(0.0f);
		hmm->random_read_model_t[j] = prob2scaledprob(0.0f);
		
		sum_prob[r1] = logsum(sum_prob[r1], hmm->random_genome_model[j]);
	}
	
	r1 = -1;
	for(j = 0; j < 25;j++){
		if(j % 5 == 0){
			r1++;
		}
		hmm->random_genome_model[j]  = hmm->random_genome_model[j]  - sum_prob[r1] ;
		hmm->random_read_model[j]  = hmm->random_read_model[j]  - sum_prob[r1] ;
	}
	for(i = 0; i < 256;i++){
		sum_prob[i] = 0;
	}
	for(i = 0; i < 65536;i++){
		r1 = (i >> 12) & 15;
		g1= (i >> 8) & 15;
		r2 = (i >> 4) & 15;
		g2 = i &15;
	       
		hmm->prob[i] = 1.0f;
		hmm->print_prob[i] = 0.0f;
		
		sum_prob[(r1  << 4) | g1 ]  += hmm->prob[i];
	}
	
	for(i = 0; i < 65536;i++){
		r1 = (i >> 12) & 15;
		g1= (i >> 8) & 15;
		r2 = (i >> 4) & 15;
		g2 = i & 15;
		if(r1>= 15 || g1 >= 15  || r2 >= 15 || g2 >= 15){
			hmm->prob[i] = prob2scaledprob(0.0f);
		}else{
			hmm->prob[i] = prob2scaledprob( hmm->prob[i] /  sum_prob[(r1 << 4) | g1 ] );
		}
	}
	free(sum_prob);
	
	for(r1 = 0;r1 < 5;r1++){
		for(g1 = 0;g1 < 5;g1++){
			//substitutions....
			// A A
			// A A
			for(r2 = 0;r2 < 5;r2++){
				sum = prob2scaledprob(0.0f);
				for(g2 = 0;g2 < 5;g2++){
					if(r2 == g2){
						//key = (r1 << 12) | (g1 << 8) | (r2 << 4) | g2;
						sum = logsum(sum, prob2scaledprob( 0.99)  );
					}else{
						sum =  logsum(sum, prob2scaledprob( 0.025)  );
					}
				}
				for(g2 = 0;g2 < 5;g2++){
					key = (r1 << 12) | (g1 << 8) | (r2 << 4) | g2;
					if(r2 == g2){
						//key = (r1 << 12) | (g1 << 8) | (r2 << 4) | g2;
						hmm->prob[key] = prob2scaledprob( 0.99)  - sum;
						//sum = logsum(sum, prob2scaledprob( 0.99)  );
					}else{
						hmm->prob[key] = prob2scaledprob( 0.025) - sum;
					}
				}
			}
			//gap open in read;
			// A -
			// A A
			r2 = 13;
			for(g2 = 0;g2 < 5;g2++){
				key = (r1 << 12) | (g1 << 8) | (r2 << 4) | g2;
				hmm->prob[key] = prob2scaledprob( 0.01) ;
			}
			//gap open in genome;
			// A A
			// A -
			for(r2 = 0;r2 < 5;r2++){
				key = (r1 << 12) | (g1 << 8) | (r2 << 4) | 13;
				hmm->prob[key] = prob2scaledprob( 0.01) ;
			}
			//read end.... gap open in read;
			// A .
			// A A
			r2 = 14;
			
			for(g2 = 0;g2 < 5;g2++){
				hmm->prob[key] = prob2scaledprob( 0.02) ;
			}
		}
	}
	r1 = 14;
	for(g1 = 0;g1 < 5;g1++){
		//start of read....
		// . A
		// A A
		for(r2 = 0;r2 < 5;r2++){
			for(g2 = 0;g2 < 5;g2++){
				key = (r1 << 12) | (g1 << 8) | (r2 << 4) | g2;
				hmm->prob[key] = prob2scaledprob( 0.02) ;
				
			}
		}
		// 5' gaps / 3' gaps
		// .  .
		// A A
		r2 = 14;
		for(g2 = 0;g2 < 5;g2++){
			key = (r1 << 12) | (g1 << 8) | (r2 << 4) | g2;
			hmm->prob[key] = prob2scaledprob( 0.02) ;
			//hmm->prob[key] = hmm->tprob[key]-sum + main_probs[9] -  main_probs[0];
		}
	}
	
	r1 = 13;
	for(g1 = 0;g1 < 5;g1++){
		//gap close
		//  - A
		// A A
		for(r2 = 0;r2 < 5;r2++){
			for(g2 = 0;g2 < 5;g2++){
				key = (r1 << 12) | (g1 << 8) | (r2 << 4) | g2;
				hmm->prob[key] = prob2scaledprob( 0.01) ;
			}
		}
		// gap extension
		// - -
		// A A
		r2 = 13;
		for(g2 = 0;g2 < 5;g2++){
			key = (r1 << 12) | (g1 << 8) | (r2 << 4) | g2;
			hmm->prob[key] = prob2scaledprob( 0.01) ;
		}
	}
	
	for(r1 = 0;r1 < 5;r1++){
		g1 = 13;
		//gap close genome
		// A A
		// - A
		
		for(r2 = 0;r2 < 5;r2++){
			for(g2 = 0;g2 < 5;g2++){
				key = (r1 << 12) | (g1 << 8) | (r2 << 4) | g2;
				hmm->prob[key] = prob2scaledprob( 0.01) ;
			}
		}
		//gap extension genome;
		// A A
		// -  -
		g2 = 13;
		for(r2 = 0;r2 < 5;r2++){
			key = (r1 << 12) | (g1 << 8) | (r2 << 4) | g2;
			hmm->prob[key] = prob2scaledprob( 0.01) ;
		}
	}
	return hmm;
ERROR:
	free_hmm(hmm);
	return NULL;
}

void free_hmm(struct hmm* hmm)
{
	int i,j;
	int x;//,y;
	
	
	//y = hmm->y;
	
	if(hmm){
		x = hmm->x;
		if(hmm->tfM){
			for(i = 0; i < hmm->num_hits;i++){
				if(hmm->tfM[i]){
					for(j = 0;j < x;j++){
						MFREE(hmm->tfM[i][j]);
					}
					MFREE(hmm->tfM[i]);
				}
			}
			MFREE(hmm->tfM);
		}
		if(hmm->tfX){
			for(i = 0; i <  hmm->num_hits;i++){
				if(hmm->tfX[i]){
					for(j = 0;j < x;j++){
						MFREE(hmm->tfX[i][j]);
					}
					MFREE(hmm->tfX[i]);
				}
			}
			MFREE(hmm->tfX);
		}
		if(hmm->tfY){
			for(i = 0; i <  hmm->num_hits;i++){
				if(hmm->tfY[i]){
					for(j = 0;j < x;j++){
						MFREE(hmm->tfY[i][j]);
					}
					MFREE(hmm->tfY[i]);
				}
			}
			MFREE(hmm->tfY);
		}
		if(hmm->fM){
			for(i = 0;i < x;i++){
				MFREE(hmm->fM[i]);
			}
			MFREE(hmm->fM);
		}
		if(hmm->fX){
			for(i = 0;i < x;i++){
				MFREE(hmm->fX[i]);
			}
			MFREE(hmm->fX);
		}
		if(hmm->fY){
			for(i = 0;i < x;i++){
				MFREE(hmm->fY[i]);
			}
			MFREE(hmm->fY);
		}
		if(hmm->bM){
			for(i = 0;i < x;i++){
				MFREE(hmm->bM[i]);
			}
			MFREE(hmm->bM);
		}
		if(hmm->bX){
			for(i = 0;i < x;i++){
				MFREE(hmm->bX[i]);
			}
			MFREE(hmm->bX);
		}
		if(hmm->bY){
			for(i = 0;i < x;i++){
				MFREE(hmm->bY[i]);
			}
			MFREE(hmm->bY);
		}
		if(hmm->traceback){
			for(i = 0;i < x;i++){
				MFREE(hmm->traceback[i]);
			}
			MFREE(hmm->traceback);
		}
		if(hmm->random_genome_model){
			MFREE(hmm->random_genome_model);
		}
		
		if(hmm->random_read_model){
			MFREE(hmm->random_read_model);
		}
		
		if(hmm->random_genome_model_t){
			MFREE(hmm->random_genome_model_t);
		}
		
		if(hmm->random_read_model_t){
			MFREE(hmm->random_read_model_t);
		}
		
		if(hmm->unaligned_scores){
			MFREE(hmm->unaligned_scores);
		}
		
		if(hmm->alignment_scores){
			MFREE(hmm->alignment_scores);
		}

		if(hmm->tmp_scores){
			MFREE(hmm->tmp_scores);
		}
		
		if(hmm->prob){
			MFREE(hmm->prob);
		}
		if(hmm->tprob){
			MFREE(hmm->tprob);
		}
		
		if(hmm->print_prob){
			MFREE(hmm->print_prob);
		}
		MFREE(hmm);
	}
}


int init_thread_hmms(struct shared_data* bsd)
{
	int i,c;
	int num_threads = 0;
	struct hmm* master_hmm = NULL;
	struct hmm* thread_hmm = NULL;
	num_threads = bsd->num_threads;

	master_hmm = bsd->master_hmm;

	for(c = 0; c < num_threads;c++){
		thread_hmm = bsd->thread_hmm[c];
		for(i = 0; i < 65536;i++){
			thread_hmm->prob[i] = master_hmm->prob[i];
			thread_hmm->tprob[i] = prob2scaledprob(0.0f);
		}
		for(i = 0; i < 25;i++){
			
			thread_hmm->random_genome_model[i] = master_hmm->random_genome_model[i];
			thread_hmm->random_genome_model_t[i] = prob2scaledprob(0.0f);
		        thread_hmm->random_read_model[i] = master_hmm->random_read_model[i];
		        thread_hmm->random_read_model_t[i] = prob2scaledprob(0.0f);
			
		}
	}
	return OK;
}

int init_shared_data_hmms(struct shared_data* bsd)
{
	struct hmm* tmp_hmm = NULL;
	int max_len = 0;
	int num_hits = 0;
	int i;
	
	ASSERT(bsd != NULL,"bsd is not allocated.");
	ASSERT(bsd->master_hmm == NULL,"master hmm already allocated.");

	
	/* hmms need to have space for longest read sequence plus the flanking genomic regions..  */
	max_len = bsd->max_seq_len + ALIGNMENT_FLANKING_LENGTH + ALIGNMENT_FLANKING_LENGTH; 
	num_hits = bsd->param->num_maxhits;
	
	
	
	RUNP(tmp_hmm = init_hmm(max_len ,max_len,num_hits));
	RUN(add_pseudo_count(tmp_hmm, bsd->pseudo_counts));
	RUN(re_estimate(tmp_hmm));
	bsd->master_hmm = tmp_hmm;
	tmp_hmm = NULL;

	MMALLOC(bsd->thread_hmm,sizeof(struct hmm*) * bsd->num_threads);
	for (i = 0; i < bsd->num_threads; i++) {
		bsd->thread_hmm[i] = NULL;
		RUNP(tmp_hmm = init_hmm(max_len ,max_len,num_hits));
		bsd->thread_hmm[i]  = tmp_hmm;
		tmp_hmm = NULL;
	}
	
	
	
	return OK;
ERROR:
	if(tmp_hmm){
		free_hmm(tmp_hmm);
	}
	return FAIL;
}




int entangle_hmms(struct shared_data* bsd)
{
	int i,j;
	int num_threads = 0;
	struct hmm* master = NULL;
	struct hmm* thread_hmm = NULL;
	num_threads = bsd->num_threads;

	master = bsd->master_hmm;
	
	for(i = 0; i < num_threads;i++){
		thread_hmm = bsd->thread_hmm[i];
		for(j = 0; j < 65536;j++){
			master->tprob[j] = logsum(master->tprob[j] , thread_hmm->tprob[j]);
		}
		for (j = 0; j < 25; j++) {
			master->random_read_model_t[j] = logsum(master->random_read_model_t[j],thread_hmm->random_read_model_t[j]);
			master->random_genome_model_t[j] = logsum(master->random_genome_model_t[j],thread_hmm->random_genome_model_t[j]);
		}
	}
	return OK;
}

int add_pseudo_count(struct hmm* hmm, float total)
{
	int r1,g1,r2,g2,key;
	for(r1 = 0;r1 < 5;r1++){
		for(g1 = 0;g1 < 5;g1++){
			//substitutions....
			// A A
			// A A
			for(r2 = 0;r2 < 5;r2++){
				for(g2 = 0;g2 < 5;g2++){
					key = (r1 << 12) | (g1 << 8) | (r2 << 4) | g2;
					if(r2 == g2){
						hmm->tprob[key]  = prob2scaledprob(total);
					}else{
						hmm->tprob[key]  = prob2scaledprob(total * 0.02);
					}
					DPRINTF2("%d %d: %f\n",r2,g2, hmm->tprob[key]);
				}
			}
			//gap open in read;
			// A -
			// A A
			r2 = 13;
			for(g2 = 0;g2 < 5;g2++){
				key = (r1 << 12) | (g1 << 8) | (r2 << 4) | g2;
				hmm->tprob[key]  = prob2scaledprob(total * 0.01);
			}
			
			//gap open in genome;
			// A A
			// A -
			for(r2 = 0;r2 < 5;r2++){
				key = (r1 << 12) | (g1 << 8) | (r2 << 4) | 13;
				hmm->tprob[key]  = prob2scaledprob(total * 0.01);
				
			}
			
			//read end.... gap open in read;
			// A .
			// A A
			r2 = 14;
			for(g2 = 0;g2 < 5;g2++){
				key = (r1 << 12) | (g1 << 8) | (r2 << 4) | g2;
				hmm->tprob[key]  = prob2scaledprob(total * 0.02);
			}
		}
	}
	r1 = 14;
	for(g1 = 0;g1 < 5;g1++){
		//start of read....
		// . A
		// A A
		for(r2 = 0;r2 < 5;r2++){
			
			for(g2 = 0;g2 < 5;g2++){
				key = (r1 << 12) | (g1 << 8) | (r2 << 4) | g2;
				if(r2 == g2){
					hmm->tprob[key]  = prob2scaledprob(total);
				}else{
					hmm->tprob[key]  = prob2scaledprob(total * 0.02);
				}
			}
		}
		// 5' gaps / 3' gaps
		// .  .
		// A A
		r2 = 14;
		for(g2 = 0;g2 < 5;g2++){
			key = (r1 << 12) | (g1 << 8) | (r2 << 4) | g2;
			hmm->tprob[key]  = prob2scaledprob(total * 0.02);
		}
	}
	
	r1 = 13;
	for(g1 = 0;g1 < 5;g1++){
		//gap close
		//  - A
		// A A
		for(r2 = 0;r2 < 5;r2++){
			for(g2 = 0;g2 < 5;g2++){
				key = (r1 << 12) | (g1 << 8) | (r2 << 4) | g2;
				if(r2 == g2){
					hmm->tprob[key]  = prob2scaledprob(total);
				}else{
					hmm->tprob[key]  = prob2scaledprob(total * 0.02);
				}
			}
		}
		// gap extension
		// - -
		// A A
		r2 = 13;
		for(g2 = 0;g2 < 5;g2++){
			key = (r1 << 12) | (g1 << 8) | (r2 << 4) | g2;
			hmm->tprob[key]  = prob2scaledprob(total * 0.01);
		}
	}
	
	for(r1 = 0;r1 < 5;r1++){
		g1 = 13;
		//gap close genome
		// A A
		// - A
		
		for(r2 = 0;r2 < 5;r2++){
			for(g2 = 0;g2 < 5;g2++){
				key = (r1 << 12) | (g1 << 8) | (r2 << 4) | g2;
				if(r2 == g2){
					hmm->tprob[key]  = prob2scaledprob(total);
				}else{
					hmm->tprob[key]  = prob2scaledprob(total * 0.02);
				}
			}
		}
		//gap extension genome;
		// A A
		// -  -
		g2 = 13;
		//added just now in 2016 - correct?
		//sum = prob2scaledprob(0.0f);
		for(r2 = 0;r2 < 5;r2++){
			key = (r1 << 12) | (g1 << 8) | (r2 << 4) | g2;
			hmm->tprob[key]  = prob2scaledprob(total * 0.01);
			//sum = logsum(sum, hmm->tprob[key] );
			
		}
	}
	return OK;
}

int re_estimate(struct hmm* hmm)
{
	float sum = prob2scaledprob(0.0f);
	//sum_prob = malloc(sizeof(float)*4096);
	int r1,g1,r2,g2,key;
	for(r1 = 0;r1 < 5;r1++){
		for(g1 = 0;g1 < 5;g1++){
			
			sum = prob2scaledprob(0.0f);
			//substitutions....
			// A A
			// A A
			for(r2 = 0;r2 < 5;r2++){
				for(g2 = 0;g2 < 5;g2++){
					key = (r1 << 12) | (g1 << 8) | (r2 << 4) | g2;
					sum = logsum(sum,hmm->tprob[key] );					
				}
			}
			//gap open in read;
			// A -
			// A A
			r2 = 13;
			for(g2 = 0;g2 < 5;g2++){
				key = (r1 << 12) | (g1 << 8) | (r2 << 4) | g2;
				sum = logsum(sum,hmm->tprob[key] );
				//hmm->tprob[key]  = prob2scaledprob(total * 0.005);
			}
			
			//gap open in genome;
			// A A
			// A -
			for(r2 = 0;r2 < 5;r2++){
				key = (r1 << 12) | (g1 << 8) | (r2 << 4) | 13;
				sum = logsum(sum,hmm->tprob[key] );
				//hmm->tprob[key]  = prob2scaledprob(total * 0.005);
			}
			
			//read end.... gap open in read;
			// A .
			// A A
			r2 = 14;
			for(g2 = 0;g2 < 5;g2++){
				key = (r1 << 12) | (g1 << 8) | (r2 << 4) | g2;
				sum = logsum(sum,hmm->tprob[key] );
				//hmm->tprob[key]  = prob2scaledprob(total * 0.02);
			}
			
			// set pro..
			
			for(r2 = 0;r2 < 5;r2++){
				for(g2 = 0;g2 < 5;g2++){
					key = (r1 << 12) | (g1 << 8) | (r2 << 4) | g2;
					hmm->prob[key]  =hmm->tprob[key] -sum;
				}
			}
			//gap open in read;
			// A -
			// A A
			r2 = 13;
			for(g2 = 0;g2 < 5;g2++){
				key = (r1 << 12) | (g1 << 8) | (r2 << 4) | g2;
				hmm->prob[key]  =hmm->tprob[key] -sum;
				//hmm->tprob[key]  = prob2scaledprob(total * 0.005);
			}
			
			//gap open in genome;
			// A A
			// A -
			for(r2 = 0;r2 < 5;r2++){
				key = (r1 << 12) | (g1 << 8) | (r2 << 4) | 13;
				hmm->prob[key]  =hmm->tprob[key] -sum;
				//hmm->tprob[key]  = prob2scaledprob(total * 0.005);
			}
			
			//read end.... gap open in read;
			// A .
			// A A
			r2 = 14;
			for(g2 = 0;g2 < 5;g2++){
				key = (r1 << 12) | (g1 << 8) | (r2 << 4) | g2;
				hmm->prob[key]  =hmm->tprob[key] -sum;
				//hmm->tprob[key]  = prob2scaledprob(total * 0.02);
			}
			
		}
	}
	r1 = 14;
	for(g1 = 0;g1 < 5;g1++){
		//start of read....
		// . A
		// A A
		
		sum = prob2scaledprob(0.0);
		
		for(r2 = 0;r2 < 5;r2++){
			
			for(g2 = 0;g2 < 5;g2++){
				key = (r1 << 12) | (g1 << 8) | (r2 << 4) | g2;
				sum = logsum(sum,hmm->tprob[key] );
			}
		}
		// 5' gaps / 3' gaps
		// .  .
		// A A
		r2 = 14;
		for(g2 = 0;g2 < 5;g2++){
			key = (r1 << 12) | (g1 << 8) | (r2 << 4) | g2;
			sum = logsum(sum,hmm->tprob[key] );
		}
		
		for(r2 = 0;r2 < 5;r2++){
			
			for(g2 = 0;g2 < 5;g2++){
				key = (r1 << 12) | (g1 << 8) | (r2 << 4) | g2;
				hmm->prob[key]  =hmm->tprob[key] -sum;
			}
		}
		// 5' gaps / 3' gaps
		// .  .
		// A A
		r2 = 14;
		for(g2 = 0;g2 < 5;g2++){
			key = (r1 << 12) | (g1 << 8) | (r2 << 4) | g2;
			hmm->prob[key]  =hmm->tprob[key] -sum;
		}
		
		
	}
	
	r1 = 13;
	for(g1 = 0;g1 < 5;g1++){
		
		sum = prob2scaledprob(0.0);
		//gap close
		//  - A
		// A A
		for(r2 = 0;r2 < 5;r2++){
			for(g2 = 0;g2 < 5;g2++){
				key = (r1 << 12) | (g1 << 8) | (r2 << 4) | g2;
				sum = logsum(sum,hmm->tprob[key] );
			}
		}
		// gap extension
		// - -
		// A A
		r2 = 13;
		for(g2 = 0;g2 < 5;g2++){
			key = (r1 << 12) | (g1 << 8) | (r2 << 4) | g2;
			sum = logsum(sum,hmm->tprob[key] );
		}
		
		for(r2 = 0;r2 < 5;r2++){
			for(g2 = 0;g2 < 5;g2++){
				key = (r1 << 12) | (g1 << 8) | (r2 << 4) | g2;
				hmm->prob[key]  =hmm->tprob[key] -sum;
			}
		}
		// gap extension
		// - -
		// A A
		r2 = 13;
		for(g2 = 0;g2 < 5;g2++){
			key = (r1 << 12) | (g1 << 8) | (r2 << 4) | g2;
			hmm->prob[key]  =hmm->tprob[key] -sum;
		}
		
	}
	g1 = 13;

	for(r1 = 0;r1 < 5;r1++){
		
		sum = prob2scaledprob(0.0);
				//gap close genome
		// A A
		// - A
		
		for(r2 = 0;r2 < 5;r2++){
			for(g2 = 0;g2 < 5;g2++){
				key = (r1 << 12) | (g1 << 8) | (r2 << 4) | g2;
				sum = logsum(sum,hmm->tprob[key] );
			}
		}
		//gap extension genome;
		// A A
		// -  -
		g2 = 13;
		//added just now in 2016 - correct?
		//sum = prob2scaledprob(0.0f);
		for(r2 = 0;r2 < 5;r2++){
			key = (r1 << 12) | (g1 << 8) | (r2 << 4) | g2;
			sum = logsum(sum,hmm->tprob[key] );
			//sum = logsum(sum, hmm->tprob[key] );
			
		}
		
		for(r2 = 0;r2 < 5;r2++){
			for(g2 = 0;g2 < 5;g2++){
				key = (r1 << 12) | (g1 << 8) | (r2 << 4) | g2;
				hmm->prob[key]  =hmm->tprob[key] -sum;
			}
		}
		//gap extension genome;
		// A A
		// -  -
		g2 = 13;
		//added just now in 2016 - correct?
		//sum = prob2scaledprob(0.0f);
		for(r2 = 0;r2 < 5;r2++){
			key = (r1 << 12) | (g1 << 8) | (r2 << 4) | g2;
			hmm->prob[key]  =hmm->tprob[key] -sum;
			//sum = logsum(sum, hmm->tprob[key] );
			
		}
		
	}
	//free(sum_prob);
	//return hmm;
	return OK;
//ERROR:
	
//	return FAIL;
}


int re_estimate_random(struct hmm* hmm)
{
	int j;
	float* sum_prob =  NULL;
	
	MMALLOC(sum_prob, sizeof(float)*256);
	//m_prob = malloc();
	
	int r1;
	
	r1 = -1;
	for(j = 0; j < 25;j++){
		if(j % 5 == 0){
			r1++;
			sum_prob[r1] = prob2scaledprob(0.0f);
		}
		sum_prob[r1] = logsum(sum_prob[r1], hmm->random_genome_model_t[j]);
		
	}
	for(j = 0; j < 5;j++){
		if(sum_prob[j] == prob2scaledprob(0.0f)){
			sum_prob[j] = prob2scaledprob(1.0f);
		}
	}
	r1 = -1;
	for(j = 0; j < 25;j++){
		if(j % 5 == 0){
			r1++;
			//sum_prob[r1] = prob2scaledprob(0.0f);
		}
		hmm->random_genome_model[j]  = hmm->random_genome_model_t[j]  - sum_prob[r1];
	}
	
	r1 = -1;
	for(j = 0; j < 25;j++){
		if(j % 5 == 0){
			r1++;
			sum_prob[r1] = prob2scaledprob(0.0f);
		}
		sum_prob[r1] = logsum(sum_prob[r1], hmm->random_read_model_t[j]);
		
	}
	for(j = 0; j < 5;j++){
		if(sum_prob[j] == prob2scaledprob(0.0f)){
			sum_prob[j] = prob2scaledprob(1.0f);
		}
	}
	r1 = -1;
	for(j = 0; j < 25;j++){
		if(j % 5 == 0){
			r1++;
			//sum_prob[r1] = prob2scaledprob(0.0f);
		}
		hmm->random_read_model[j]  = hmm->random_read_model_t[j]  - sum_prob[r1];
	}
	MFREE(sum_prob);
	return OK;
ERROR:
	return FAIL;
}
