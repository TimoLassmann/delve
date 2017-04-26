

#include "tldevel.h"
#include "htsglue.h"
#include "rtr.h"


#include "delve_struct.h"


int add_genome_pseudocounts(struct shared_data* bsd)
{
	struct delve_genome_region_data* dgd = NULL;
	struct rtr_data* rtree = NULL;
	float pseudocounts = 0.0f;
	int i;

	ASSERT(bsd != NULL,"No shared data.");

	rtree = bsd->rtree;

	pseudocounts = prob2scaledprob((float) bsd->param->genome_pseudocounts);
	
	
	for (i = 0; i < rtree->stats_num_interval; i++) {
		dgd = (struct delve_genome_region_data*) rtree->flat_interval[i]->data;
		dgd->prior_e =  pseudocounts;
		
	}
	
	return OK;
ERROR:
	return FAIL;
}

int re_estimate_genome_priors(struct shared_data* bsd)
{
	struct delve_genome_region_data* dgd = NULL;
	struct rtr_data* rtree = NULL;
	float sum = 0.0;
	int i;

	ASSERT(bsd != NULL,"No shared data.");

	rtree = bsd->rtree;

	
	sum = prob2scaledprob(0.0);
	for (i = 0; i < rtree->stats_num_interval; i++) {
		dgd = (struct delve_genome_region_data*) rtree->flat_interval[i]->data;
		sum = logsum(sum,dgd->prior_e);
		
	}

	for (i = 0; i < rtree->stats_num_interval; i++) {
		dgd = (struct delve_genome_region_data*) rtree->flat_interval[i]->data;
		dgd->prior = dgd->prior_e -sum;
//		fprintf(stdout,"%d:%d %f %f %f\n",i, dgd->count,scaledprob2prob(dgd->prior),dgd->prior_e,sum);
		
	}
	
	return OK;
ERROR:
	return FAIL;
}	



int copy_genome_priors_to_gc(struct shared_data* bsd)
{
	int64_t* val = NULL;
	struct sam_bam_entry** buffer = NULL;
	struct rtr_data* rtree = NULL;
	struct genome_sequences** gc = NULL;
	struct delve_genome_region_data* dgd = NULL;
	int i,j;
	int num_seq; 
	
	ASSERT(bsd != NULL,"No shared data found.");
	ASSERT(bsd->sb_file != NULL,"No open file found");

	num_seq = bsd->sb_file->num_read;
	buffer = bsd->sb_file->buffer;
	rtree = bsd->rtree;
	gc = bsd->gc;

	MMALLOC(val,sizeof(int64_t)*2);
	
	for(i = 0; i < num_seq; i++) {
	        for (j = 0; j < buffer[i]->num_hits ; j++) {
			val[0] = buffer[i]->start[j];
			val[1] = buffer[i]->stop[j];
        		RUNP(dgd = rtree->query(rtree,val));
        		gc[i]->prior[j] = dgd->prior;
			gc[i]->prior_e[j] = prob2scaledprob(0.0f);
			
		}
	}
	MFREE(val);
	return OK;
ERROR:
	MFREE(val);
	return FAIL;
}

int entangle_genome_priors(struct shared_data* bsd)
{
	int64_t* val = NULL;
	struct sam_bam_entry** buffer = NULL;
	struct rtr_data* rtree = NULL;
	struct genome_sequences** gc = NULL;
	struct delve_genome_region_data* dgd = NULL;
	int i,j;
	int num_seq; 
	
	ASSERT(bsd != NULL,"No shared data found.");
	num_seq = bsd->sb_file->num_read;
	buffer = bsd->sb_file->buffer;
	rtree = bsd->rtree;
	gc = bsd->gc;

	MMALLOC(val,sizeof(int64_t)*2);
	
	for(i = 0; i < num_seq; i++) {
	        for (j = 0; j < buffer[i]->num_hits ; j++) {
			val[0] = buffer[i]->start[j];
			val[1] = buffer[i]->stop[j];
        		RUNP(dgd = rtree->query(rtree,val));
        	        dgd->prior_e = logsum(dgd->prior_e,gc[i]->prior_e[j]);
			
		}
	}
	MFREE(val);
	return OK;
ERROR:
	MFREE(val);
	return FAIL;
}
