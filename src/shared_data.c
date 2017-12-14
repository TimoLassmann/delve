

/* other libraries  */
//#include <stdio.h>
//#include <stdint.h>

/* my libraries */
#include "tldevel.h"
#include "thr_pool.h"
#include "htsglue.h"
#include "pwrite.h"
#include "rtr.h"

/* delve headers */
#include "delve_struct.h"
#include "genome_sequences.h"
#include "hmm.h"
#include "shared_data.h"


struct shared_data* init_shared_data(struct parameters* param, int buffer_size)
{
        struct shared_data* bsd = NULL;
        int i;
	
        MMALLOC(bsd,sizeof(struct shared_data));
        bsd->param = NULL;
        bsd->gc = NULL;
        bsd->sb_file = NULL;
        bsd->rtree = NULL;
        bsd->pool = NULL;
        bsd->index = NULL;
        bsd->master_hmm = NULL;
        bsd->thread_hmm = NULL;
        bsd->g_int_working = NULL;
        bsd->thread_forward = NULL;
		
        bsd->free = free_shared_data;
        bsd->buffer_size = 0;
        bsd->num_threads = 0;
        bsd->num_maxhits = 0;
        bsd->max_seq_len = 0;
        //bsd->pseudo_counts = param->pseudocounts;
        bsd->iteration = 0;
        /* assignment  */
        bsd->free = free_shared_data;
        bsd->buffer_size = buffer_size;
        bsd->num_maxhits = param->num_maxhits;
        bsd->param = param;
        bsd->num_threads = param->num_threads;
        bsd->statmode = 0;
        /* allocation  */
        if((bsd->pool = thr_pool_create(bsd->num_threads+1, bsd->num_threads+1, 0, 0)) == NULL) ERROR_MSG("Creating pool thread failed.");
        RUNP(bsd->gc = init_genome_sequences(bsd->buffer_size, bsd->num_maxhits));
        RUNP(bsd->index = get_faidx(bsd->param->genome));

        //RUNP(bsd->pw = init_pwrite_main(param->out_file,param->num_threads,BUFFER_P_WRITE_SIZE));
        bsd->pw = NULL;
        MMALLOC(bsd->thread_forward,sizeof(double)* param->num_threads);
	
        MMALLOC(bsd->g_int_working, sizeof(struct genome_interval*) * param->num_threads);
        for (i = 0; i < param->num_threads; i++) {
                bsd->g_int_working[i] = NULL;
                bsd->g_int_working[i] = init_genome_interval(0,0,0);
        }
        return bsd;
ERROR:
        if(bsd){
                bsd->free(bsd);
        }
        return NULL;
}


void free_shared_data(struct shared_data* bsd)
{
        int i;
        if(bsd){
                if(bsd->g_int_working){
                        for (i = 0; i < bsd-> num_threads; i++) {
                                free_genome_interval(bsd->g_int_working[i]);
                        }
                        MFREE(bsd->g_int_working);
                }
                if(bsd->thread_forward){
                        MFREE(bsd->thread_forward);
                }
                if(bsd->master_hmm){
                        free_hmm(bsd->master_hmm);
                }

                if(bsd->thread_hmm){
                        for (i = 0; i < bsd->num_threads; i++) {
                                free_hmm(bsd->thread_hmm[i]);
					 
                        }
                        MFREE(bsd->thread_hmm);
                }

                if(bsd->rtree){		
                        bsd->rtree->free(bsd->rtree);
                }
                if(bsd->pool){
                        thr_pool_destroy(bsd->pool);
                }
                if(bsd->gc){
                        free_genome_sequences(bsd->gc ,bsd->buffer_size,bsd->num_maxhits);
                }
                if(bsd->index){
                        free_faidx(bsd->index);
                }
                if(bsd->pw){
                        bsd->pw->free(bsd->pw);
                }
                MFREE(bsd);
        }
}

