

#include <stdio.h>
#include <stdint.h>

/* my libs  */
#include "tldevel.h"
#include "rtr.h"
#include "htsglue.h"

#include "delve_struct.h"

int resolve_delve_genome_region_data(struct  rtree_interval* a,struct rtree_interval* b);



struct rtr_data* build_rtree(struct sam_bam_file* sb_file )
{
	struct rtr_data* rtree = NULL;
	struct delve_genome_region_data* dgd = NULL;
	int64_t* val = NULL;
	int i,j;
	int32_t id = 1;
	
	MMALLOC(val,sizeof(int64_t)*2);
	RUNP(rtree = init_rtr_data(1 , 5 ,sb_file->buffer_size ));

	rtree->resolve_dup = resolve_delve_genome_region_data;
	
	while(1){
		RUN(read_SAMBAM_chunk(sb_file,1,0));
	        if(!sb_file->num_read){
			break;
		}
		LOG_MSG("read: %d",sb_file->num_read);
		for (i = 0; i < sb_file->num_read; i++) {
			for (j=0; j < sb_file->buffer[i]->num_hits ; j++) {
				MMALLOC(dgd,sizeof(struct delve_genome_region_data));
				dgd->count = 1;
				dgd->alignment_weigth = 0.0f;
				dgd->prior = prob2scaledprob(1.0f);
				dgd->prior_e = prob2scaledprob(0.0f);
				
				
				val[0] = sb_file->buffer[i]->start[j];
				val[1] = sb_file->buffer[i]->stop[j];
				rtree->insert(rtree,val,dgd , id,1,1);
				dgd = NULL;
				id++;
			}
		}
        }
	MFREE(val);
	RUN(rtree->flatten_rtree(rtree));

	
	RUN(rtree->print_rtree(rtree,rtree->root));
	for(i = 0; i < rtree->stats_num_interval;i++){
		fprintf(stdout,"%d %d\n", rtree->flat_interval[i]->count,i);
	}
	
	return rtree;
ERROR:
	MFREE(val);
	if(rtree){
		rtree->free(rtree);
	}
	return NULL;
}

int free_delve_region_data_from_tree(struct rtr_data* rtree)
{
	int i;

	ASSERT(rtree != NULL,"No tree");
	LOG_MSG("cleaning tree: %d entries.",rtree->stats_num_interval);
	for(i = 0; i < rtree->stats_num_interval;i++){
		if(rtree->flat_interval[i]->data){
			MFREE(rtree->flat_interval[i]->data);
		}
	}
	return OK;
ERROR:
	return FAIL;
}


int set_sequence_weigth(struct shared_data* bsd)
{
	int64_t* val = NULL;
	struct sam_bam_entry** buffer = NULL;
	struct rtr_data* rtree = NULL;
	struct genome_sequences** gc = NULL;
	struct delve_genome_region_data* dgd = NULL;
	int i,j;
	int num_seq; 
	float weigth = 0.0;
	float tmp;
	
	ASSERT(bsd != NULL,"No shared data found.");
	num_seq = bsd->sb_file->num_read;
	buffer = bsd->sb_file->buffer;
	rtree = bsd->rtree;
	gc = bsd->gc;
	
	MMALLOC(val,sizeof(int64_t)*2);
	
	for (i = 0; i < num_seq;i++) {
		weigth = 0.0;
		for (j = 0; j < buffer[i]->num_hits ; j++) {
			val[0] = buffer[i]->start[j];
			val[1] = buffer[i]->stop[j];
			// int query(struct rtr_data* rtrd , int64_t* val,int32_t* identifier,int32_t* count)
			RUNP(dgd = rtree->query(rtree,val));
			
			tmp =  log((float) dgd->count) / (float) dgd->count;
			if(tmp > weigth){
				weigth = tmp;
			}
			gc[i]->alignment_weigth[j] =  prob2scaledprob(weigth);
		}
	}
	MFREE(val);
	return OK;
ERROR:
	MFREE(val);
	return FAIL;
}

int resolve_delve_genome_region_data(struct  rtree_interval* a,struct rtree_interval* b)
{
	struct delve_genome_region_data* org = NULL;
	struct delve_genome_region_data* new = NULL;

	ASSERT(a != NULL,"org interval is NULL");	
	ASSERT(b != NULL,"new interval is NULL");

	org = (struct delve_genome_region_data*) a->data;
	new = (struct delve_genome_region_data*) b->data;

	org->count = org->count + new->count;
	
	a->data = (void*) org;	

	MFREE(b->data);
	b->data = NULL;
	return OK;
ERROR:
	return FAIL;
}
