#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <string.h>
#include <inttypes.h>
#include  <ctype.h>

#include "tldevel.h"
#include "rbtree.h"
#include "htsglue.h"
#include "thr_pool.h"
#include "rtr.h"
#include "pwrite.h"
#include "delve.h"

#define OPT_NTHREAD 1
#define OPT_OUT 2

struct parameters{
	char** infiles;
	char* genome;
	char* aln_infile;
	char* hmm_file;
	char* out_file;
	int num_infiles;
	int num_threads;
	int num_maxhits;
};

struct shared_data{
	struct parameters* param;
	struct genome_sequences** gc;
	struct sam_bam_file* sb_file;
	struct rtr_data* rtree;
	struct thr_pool* pool;
	struct pwrite_main* pw;
	struct hmm** thread_hmm;
	struct hmm* master_hmm;
	struct genome_interval** g_int_working;	
	double* thread_forward;
	faidx_t*  index;
	void (*free)(struct shared_data* bsd);
	FILE* fptr_out;
	int buffer_size;
	int num_threads;
	int num_maxhits;
	int max_seq_len; 
	float pseudo_counts;
};

struct thread_data{
	struct shared_data* bsd;
	int thread_id;
};

static pthread_mutex_t avail_mtx;

/* related to shared memory */
struct shared_data* init_shared_data(struct parameters* param, int buffer_size);
int init_shared_data_hmms(struct shared_data* bsd);
void free_shared_data(struct shared_data* bsd);

/* related to thread data  */
struct thread_data** init_thread_data(struct shared_data* bsd,int num_threads);
void free_thread_data(struct thread_data** td,int num_threads);

struct rtr_data* build_rtree(struct sam_bam_file* sb_file );
int set_sequence_weigth(struct shared_data* bsd);

int add_genome_sequences(struct shared_data* bsd);

//int add_genome_sequences(struct sam_bam_file* sb_file,char* genome);
//int remove_genome_sequences(struct sam_bam_file* sb_file);
int init_thread_hmms(struct shared_data* bsd);

int get_start_stop_for_threads(const int num_threads, const int num_seq, const int id, int *start,int *stop);

int run_estimate_random_models(struct shared_data* bsd);
int run_estimate_sequence_model(struct shared_data* bsd);
int run_score_alignments(struct shared_data* bsd);

int entangle_hmms(struct shared_data* bsd);

int run_delve(struct parameters* param);

struct genome_sequences** init_genome_sequences(int num, int num_maxhits);
int clear_genome_sequences(struct genome_sequences** gc, int num,int num_maxhits);
void free_genome_sequences(struct genome_sequences** gc, int num,int num_maxhits);

int reverse(uint8_t* p,int len);

int convert_buffer_ACGT_to_0123(struct sam_bam_file* sb_file);
int write_sam_header(struct sam_bam_file* sb_file,FILE* out);

int main (int argc,char *argv[]) 
{		
	struct parameters* param = NULL;
	int i,c;
	tlog.echo_build_config();
     	MMALLOC(param, sizeof(struct parameters));
	param->infiles = NULL;
	param->genome = NULL;
	param->aln_infile = NULL;
	param->hmm_file = NULL;
	param->out_file = NULL;
	param->num_infiles = 0;
	param->num_threads = 4;
	param->num_maxhits = 10;
	
	while (1){	
		static struct option long_options[] ={
			{"t",required_argument,0,OPT_NTHREAD},
			{"out",required_argument,0,OPT_OUT},
			{"help",0,0,'h'},
			{0, 0, 0, 0}
		};
		
		int option_index = 0;
		c = getopt_long_only (argc, argv,"h",long_options, &option_index);
		
		if (c == -1){
			break;
		}
		
		switch(c) {
			case OPT_NTHREAD:
				param->num_threads = atoi(optarg);
				break;
			case OPT_OUT:
				param->out_file = optarg;
				break;
		case 'h':
			fprintf(stdout,"GAGA\n");
			MFREE(param);
			exit(EXIT_SUCCESS);
			break;
			default:
				fprintf(stdout,"GAGASFGA\n");
				ERROR_MSG("not recognized");
			break;
		}
	}

	MMALLOC(param->infiles,sizeof(char*)* (argc-optind));
	
	c = 0;
	while (optind < argc){
		param->infiles[c] = argv[optind++];
		param->num_infiles++;
		c++;
	}

	tlog.log_message("Starting run");
	
	if(!param->num_infiles){
		ERROR_MSG("ks tools requires at least one input\n");
	}
	
	//Sanity checks
	//check if input files exist:
	for(i = 0; i < param->num_infiles;i++){
		if(my_file_exists(param->infiles[i]) == 0){
			ERROR_MSG("file \"%s\" does not exist.",param->infiles[i] );
		}
	}
	param->genome = param->infiles[0];
	param->aln_infile = param->infiles[1];
		
	init_logsum();
	
	RUN(run_delve(param));
	
	if(param){
		//for(i = 0; i < param->num_infiles;i++){
		//	MFREE(param->infiles[i]);
		//}
		MFREE(param->infiles);
		//MFREE(param->out_prefix);
		MFREE(param);
	}
	return EXIT_SUCCESS;
ERROR:
	fprintf(stdout,"\n  Try run with  --help.\n\n");
	if(param){
		MFREE(param->infiles);
		MFREE(param);
	}
	return EXIT_FAILURE;
}

int run_delve(struct parameters* param)
{
	struct shared_data* bsd = NULL;
	
	int buffer_size = MAXNUMQUERY;	
	
	init_logsum();

	RUNP(bsd = init_shared_data(param,buffer_size));	
	/* 1. build rtree
	   requires reading all alignments in file. This means that I need to open the 
	   file, loop, close and then re-open.... 
	   Additional : this means that max_seq_len will be recorded.. 
	*/
        LOG_MSG("Quantifying read depth at pu1tatitve mapping locations.");
	RUNP(bsd->sb_file = open_SAMBAMfile(bsd->param->aln_infile,bsd->buffer_size,bsd->num_maxhits , 0,0));
	RUNP(bsd->rtree = build_rtree(bsd->sb_file));


	
	/* assign max_seq_len */
	RUN(get_max_seq_len_from_sb_buffer(bsd->sb_file,&bsd->max_seq_len));
	
	RUN(close_SAMBAMfile(bsd->sb_file));

	/* init all HMMs; now that I have the max_seq_len... */
        /* allocate HMMs... */
	RUN(init_shared_data_hmms(bsd));
	
	/* 2. build model 
	   Here I read in the first X (MAXNUMQUERY sequences and re-use them for model building 
	   etc. I.e. I open, read and only after all models are trained I close....
	*/
	/* open file again - read all alignments (param "1","0") */
	LOG_MSG("Opening alignment file...");
	RUNP(bsd->sb_file = open_SAMBAMfile(bsd->param->aln_infile,bsd->buffer_size,bsd->num_maxhits,0,0));
     	RUN(read_SAMBAM_chunk(bsd->sb_file,1,0));
	RUN(convert_buffer_ACGT_to_0123(bsd->sb_file));
	/* Retrieve all genome sequences "hit" by reads in the first SAM/BAM chunk
	   Note: add_genome_sequences, specific to delve (not a library function) will convert the sequences to 0123... 
	 */
	RUN(add_genome_sequences(bsd));

	/* set weigth of allocated sequences */
	RUN(set_sequence_weigth(bsd));

	/* Estimate random models  */
	LOG_MSG("Estimating Random model (%d)",bsd->sb_file->num_read);
	RUN(run_estimate_random_models(bsd));
        LOG_MSG("done");

	/* Estmate Sequence Model  */
	LOG_MSG("Estimating Sequence model...");
	RUN(run_estimate_sequence_model(bsd));
        RUN(clear_genome_sequences(bsd->gc,bsd->buffer_size, bsd->num_maxhits));
     	RUN(close_SAMBAMfile(bsd->sb_file));
	LOG_MSG("done");	

	/* "gently" add in estimation of genome priors using sumuklated annealing..  */

	
	
	
		
	/* Generating alignments  */
	RUNP(bsd->sb_file = open_SAMBAMfile(bsd->param->aln_infile,bsd->buffer_size,bsd->num_maxhits,0,0));

	LOG_MSG("Generating alignments...");
	RUN(write_sam_header(bsd->sb_file, bsd->pw->out_ptr ));
	
	while(1){
		RUN(read_SAMBAM_chunk(bsd->sb_file,1,0));
		if(!bsd->sb_file->num_read){
			break;
		}
	        //LOG_MSG("read:%d",bsd->sb_file->num_read);
		RUN(add_genome_sequences(bsd));
		RUN(convert_buffer_ACGT_to_0123(bsd->sb_file));
		RUN(run_score_alignments(bsd));	       
		RUN(clear_genome_sequences(bsd->gc,bsd->buffer_size, bsd->num_maxhits));
        }
	
	RUN(close_SAMBAMfile(bsd->sb_file));
	LOG_MSG("Done.");
	bsd->free(bsd);
	return OK;
ERROR:
	if(bsd){
		bsd->free(bsd);
	}
	return FAIL;
}


int run_estimate_random_models(struct shared_data* bsd)
{
	struct thread_data** td = NULL;
	int i;
	int num_threads;
	int status;
	num_threads = bsd->param->num_threads;

	/* I think I should add pseudocounts  */

	RUN(add_pseudo_count(bsd->master_hmm, bsd->pseudo_counts));
	RUN(re_estimate(bsd->master_hmm)); /*  */
	
	/* copy parameters from master hmm into copies used by threads...  */
	RUN(init_thread_hmms(bsd));

	/* initialize datastructs to pass bsd (shared) and thread_id's (private)  */
	RUNP(td = init_thread_data(bsd,num_threads));

	/* kick off jobs  */
	for(i = 0; i < num_threads;i++){
        	if((status = thr_pool_queue(bsd->pool,do_baum_welch_thread_random_model,td[i])) == -1) fprintf(stderr,"Adding job to queue failed.");	
	}
	/* wait for all jobs to finish */
	thr_pool_wait(bsd->pool);

	/* entangle HMMs - copy estimated counts from thread hmm copies back into the master HMM. */
	RUN(entangle_hmms(bsd));
	
	/* re-estimate parameters.. */
	RUN(re_estimate_random(bsd->master_hmm));
	
	
	free_thread_data(td,num_threads);
	return OK;
ERROR:
	free_thread_data(td,num_threads);
	return FAIL;
}

int run_estimate_sequence_model(struct shared_data* bsd)
{
	struct thread_data** td = NULL;
	int i,iter;
	int num_threads;
	int status;
	int iterations;

	iterations = 3;
	num_threads = bsd->param->num_threads;
	/* initialize datastructs to pass bsd (shared) and thread_id's (private)  */
	RUNP(td = init_thread_data(bsd,num_threads));
	/* I think I should add pseudocounts  */
	for(iter = 0 ; iter< iterations;iter++){
		LOG_MSG("Iteration %d.",iter);
		RUN(add_pseudo_count(bsd->master_hmm, bsd->pseudo_counts));
	 	/* copy parameters from master hmm into copies used by threads...  */
		RUN(init_thread_hmms(bsd));

		for(i = 0; i < num_threads;i++){
			bsd->thread_forward[i] = 0.0; // (i.e P = 1.0);
			LOG_MSG("%d score %f.",i,bsd->thread_forward[i]);
		}
		/* kick off jobs  */
		for(i = 0; i < num_threads;i++){
			if((status = thr_pool_queue(bsd->pool,do_baum_welch_thread,td[i])) == -1) fprintf(stderr,"Adding job to queue failed.");	
		}
		/* wait for all jobs to finish */
		thr_pool_wait(bsd->pool);
		
		for(i = 0; i < num_threads;i++){
			fprintf(stderr,"%d score %f\n",i,bsd->thread_forward[i]);
//		        LOG_MSG("%d score: %f\n", bsd->sb_file->buffer[i]->fscore);
		}
		/* entangle HMMs - copy estimated counts from thread hmm copies back into the master HMM. */
		RUN(entangle_hmms(bsd));
		/* re-estimate parameters.. */
		RUN(re_estimate(bsd->master_hmm));
		/* re-estimate random model more... */
		RUN(re_estimate_random(bsd->master_hmm));
	}
	free_thread_data(td,num_threads);
	return OK;
ERROR:
	free_thread_data(td,num_threads);
	return FAIL;
}


int run_score_alignments(struct shared_data* bsd)
{
	struct thread_data** td = NULL;
	int i;
	int num_threads;
	int status;
	num_threads = bsd->param->num_threads;
	/* I think I should add pseudocounts  */
	
	/* copy parameters from master hmm into copies used by threads...  */
	RUN(init_thread_hmms(bsd));

	/* initialize datastructs to pass bsd (shared) and thread_id's (private)  */
	RUNP(td = init_thread_data(bsd,num_threads));
	if((status = thr_pool_queue(bsd->pool, bsd->pw->write_thread_function, bsd->pw)) == -1) fprintf(stderr,"Adding job to queue failed.");
	/* Important! otherwise a deadlock is possible/likely... */
	RUN(bsd->pw->write_wait(bsd->pw));
	
	/* kick off jobs  */
	for(i = 0; i < num_threads;i++){ // numthreads includes the one writer thread ... need to take away one
        	if((status = thr_pool_queue(bsd->pool,do_score_alignments_thread_hmm ,td[i])) == -1) fprintf(stderr,"Adding job to queue failed.");
	}
        /* wait for all jobs to finish */
	thr_pool_wait(bsd->pool);
	
	
	free_thread_data(td,num_threads);
	return OK;
ERROR:
	free_thread_data(td,num_threads);
	return FAIL;
}

struct thread_data** init_thread_data(struct shared_data* bsd,int num_threads)
{
	int i;
	struct thread_data** td = NULL;
		
	MMALLOC(td,sizeof(struct thread_data*) * num_threads);

	for(i = 0; i < num_threads;i++){
		td[i] = NULL;
		MMALLOC(td[i],sizeof(struct thread_data));
		td[i]->thread_id = i;
		td[i]->bsd = bsd;
		/* kick off jobs? */
	}
	return td;
ERROR:
	return NULL;
}


void free_thread_data(struct thread_data** td, int num_threads)
{
	int i;
	if(td){
		for(i = 0; i < num_threads;i++){
			if(td[i]){
				MFREE(td[i]);
			}
		}
		MFREE(td);
	}
}


struct rtr_data* build_rtree(struct sam_bam_file* sb_file )
{
	struct rtr_data* rtree = NULL;
	int64_t* val = NULL;
	int i,j;
	int32_t id = 1;
	
	MMALLOC(val,sizeof(int64_t)*2);
	RUNP(rtree = init_rtr_data(1 , 5 ,sb_file->buffer_size ));
	while(1){
		RUN(read_SAMBAM_chunk(sb_file,1,0));
		//DPRINTF2("read:%d",sb_file->num_read );
		//RUN(read_sam_bam_chunk(infile,data->si, buffer,1,&numseq),"Read sam/bam chunk in thread %d failed", data->threadID);
		//if((status = read_sam_bam_chunk(infile,data->si, buffer,1,&numseq)) != kslOK)  exit(status);
		
		if(!sb_file->num_read){
			break;
		}
		LOG_MSG("read: %d",sb_file->num_read);
		for (i = 0; i < sb_file->num_read; i++) {
			for (j=0; j < sb_file->buffer[i]->num_hits ; j++) {
				val[0] = sb_file->buffer[i]->start[j];
				val[1] = sb_file->buffer[i]->stop[j];
				rtree->insert(rtree,val,id,1,1);
				id++;
			}
		}
        }
	MFREE(val);
	RUN(rtree->flatten_rtree(rtree));	
	return rtree;
ERROR:
	MFREE(val);
	if(rtree){
		rtree->free(rtree);
	}
	return NULL;
}


int set_sequence_weigth(struct shared_data* bsd)
{
	int64_t* val = NULL;
	struct sam_bam_entry** buffer = NULL;
	struct rtr_data* rtree = NULL;
	struct genome_sequences** gc = NULL;
	
	int i,j;
	int num_seq; 
	int32_t count,id,sum;
	float weigth = 0.0;
	float tmp;
	
	ASSERT(bsd != NULL,"No shared data found.");

	num_seq = bsd->sb_file->num_read;
	buffer = bsd->sb_file->buffer;
	rtree = bsd->rtree;
	gc = bsd->gc;
	sum = 0;
	for (i = 0; i < rtree->stats_num_interval; i++) {
		sum = sum + rtree->flat_interval[i]->count; 
	}
	LOG_MSG("total positions hit: %d\n",sum);

	MMALLOC(val,sizeof(int64_t)*2);
	
	for (i = 0; i < num_seq; i++) {
		weigth = 0.0;
		for (j = 0; j < buffer[i]->num_hits ; j++) {
			val[0] = buffer[i]->start[j];
			val[1] = buffer[i]->stop[j];
			// int query(struct rtr_data* rtrd , int64_t* val,int32_t* identifier,int32_t* count)
			RUN(rtree->query(rtree,val,&id,&count));
			tmp =  log((float) count) / (float) count;
			if(tmp > weigth){
				weigth = tmp;
			}
			//if(i % 1000 == 99){
				//	fprintf(stdout,"%lld-%lld id:%d, count:%d (%d) w:%f %f\n",val[0],val[1],id,count,sum,weigth , prob2scaledprob(weigth));
			//}
			gc[i]->priors[j] =  prob2scaledprob(weigth);
		}
	}
	MFREE(val);
	return OK;
ERROR:
	
	MFREE(val);
	return FAIL;
}

/* int run_pHMM(struct hmm* localhmm,struct sam_bam_file* sb_file ,struct genome_sequences** gc, faidx_t*  index,int num_threads ,int size, int mode, FILE* fout)
{
	struct thread_data* thread_data_array = NULL;

	
	
	pthread_t threads[num_threads];
	pthread_attr_t attr;
	int t;
	int interval = 0;
	int rc;


	interval =  (int)((double)size /(double)num_threads);

	MMALLOC(thread_data_array, sizeof(struct thread_data)* num_threads);
	
	for(t = 0;t < num_threads - 1;t++) {
		thread_data_array[t].sb_file = sb_file;
		thread_data_array[t].index = index;
		thread_data_array[t].gc   = gc;
	//	thread_data_array[t].db = db;
		thread_data_array[t].hmm = init_hmm(100,100,10);
		thread_data_array[t].start = t*interval;
		thread_data_array[t].end = t*interval + interval;
		thread_data_array[t].fout = fout;
		
		RUNP(thread_data_array[t].g_int = init_genome_interval(0,0,0));
	
	}
	thread_data_array[num_threads - 1].gc = gc;
//	thread_data_array[num_threads - 1].db = db;
	thread_data_array[num_threads - 1].sb_file = sb_file;
	thread_data_array[num_threads - 1].index = index;
	thread_data_array[num_threads - 1].hmm = init_hmm(100, 100,10);
	thread_data_array[num_threads - 1].start = t*interval;
	thread_data_array[num_threads - 1].end = size;
	thread_data_array[num_threads - 1].fout = fout;
	RUNP(thread_data_array[num_threads-1].g_int = init_genome_interval(0,0,0));
	
	
	init_thread_hmms(thread_data_array,localhmm,num_threads);
	
	rc = pthread_attr_init(&attr);
	if(rc){
		fprintf(stderr,"ERROR; return code from pthread_attr_init() is %d\n", rc);
		exit(-1);
	}
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	
	for(t = 0;t < num_threads;t++) {
		switch (mode) {
			case RUN_FULL_HMM:
				rc = pthread_create(&threads[t], &attr, do_baum_welch_thread, (void *) &thread_data_array[t]);
				break;
			case RUN_RANDOM_HMM:
				rc = pthread_create(&threads[t], &attr, do_baum_welch_thread_random_model, (void *) &thread_data_array[t]);
				break;
			case RUN_SCORE_HMM:
				rc = pthread_create(&threads[t], &attr, do_score_alignments_thread_hmm, (void *) &thread_data_array[t]);
				break;
			default:
				fprintf(stderr,"ERROR: hmm mode not found...%d (%d	%d	%d)\n",  mode, RUN_FULL_HMM, RUN_RANDOM_HMM,RUN_SCORE_HMM );
				break;
		}
		if (rc) {
			fprintf(stderr,"ERROR; return code from pthread_create() is %d\n", rc);
			exit(-1);
		}
	}
	
	pthread_attr_destroy(&attr);
	
	for (t = 0;t < num_threads;t++){
		rc = pthread_join(threads[t], NULL);
		if (rc){
			fprintf(stderr,"ERROR; return code from pthread_join()is %d\n", rc);
			exit(-1);
		}
	}
	
	if(mode != RUN_SCORE_HMM){
		entangle_hmms(thread_data_array,localhmm,num_threads);
	}
	
	for(t = 0;t < num_threads;t++) {
		free_hmm(thread_data_array[t].hmm);
		free_genome_interval(thread_data_array[t].g_int);
	}
	
	free(thread_data_array);
	return OK;
ERROR:
	return FAIL;
}
*/





void* do_baum_welch_thread_random_model(void *threadarg)
{
	struct thread_data *data;
	struct genome_sequences** gc = NULL;
	struct sam_bam_entry** buffer = NULL;
	struct hmm* hmm = NULL;
	int thread_id = -1;
	int num_threads = 0;
	int num_sequences = 0;
	int i,hit;
	int start,stop;
	
	data = (struct thread_data *) threadarg;

	thread_id = data->thread_id;

	num_sequences = data->bsd->sb_file->num_read;
	num_threads = data->bsd->num_threads;
	
	hmm = data->bsd->thread_hmm[thread_id];
	buffer = data->bsd->sb_file->buffer;
	gc = data->bsd->gc;
	start = 0;
	stop = 0;
	
	RUN(get_start_stop_for_threads(num_threads, num_sequences,thread_id,&start,&stop));
	
	//struct genome_interval* g_int = NULL;

	//struct read_info** ri = data->ri;

	//struct db* db = data->db;
	//int j;
	//int strand = 0;
	
//	char* seq = 0;
	char* genomic_sequence = NULL;
	int len;
	
//	seq = malloc(sizeof(char) * 100);
	//g_int = data->g_int;
	for(i = start; i < stop;i++){
       		for(hit = 0; hit < buffer[i]->num_hits;hit++){
			
			//DPRINTF2("%s: %d ->%d", data->sb_file->buffer[i]->name,data->sb_file->buffer[i]->start[hit],data->sb_file->buffer[i]->stop[hit] );
			
			//RUN(get_chr_start_stop(data->sb_file,i,hit,g_int),"get_chr_start_stop failed");
			//g_int->start -=ALIGNMENT_FLANKING_LENGTH;
			//g_int->stop += ALIGNMENT_FLANKING_LENGTH;
			//DPRINTF2("%s:%d-%d \n",g_int->chromosome,g_int->start,g_int->stop);
			genomic_sequence = gc[i]->genomic_sequences[hit];//  data->sb_file->buffer[i]->genomic_sequences[hit];
			len = gc[i]->g_len[hit];//  data->sb_file->buffer[i]->g_len[hit];
			hmm = random_model_calc(hmm,buffer[i]->sequence,  genomic_sequence,buffer[i]->len,  len, prob2scaledprob(1.0),0);
			

		}
		//DPRINTF1("modeling read sequence %d : %s",i,data->sb_file->buffer[i]->sequence);
		hmm = random_model_calc(hmm,buffer[i]->sequence,  genomic_sequence,buffer[i]->len,  len, prob2scaledprob(1.0),1);
		//		ri[i]->random_read_score = hmm->random_score;
	}
//	free(seq);
//	free_genome_interval(g_int);

	return NULL;
ERROR:
	return NULL;
}

void* do_baum_welch_thread(void *threadarg)
{
	struct thread_data *data;
	struct genome_sequences** gc = NULL;
	struct sam_bam_entry** buffer = NULL;
	struct hmm* hmm = NULL;
	int thread_id = -1;
	int num_threads = 0;
	int num_sequences = 0;
	int i;
	int start;
	int stop;
	int max_num_hits  = 0;
	
	data = (struct thread_data *) threadarg;

	thread_id = data->thread_id;

	num_sequences = data->bsd->sb_file->num_read;
	num_threads = data->bsd->num_threads;
	max_num_hits = data->bsd->sb_file->max_num_hits;
	
	hmm = data->bsd->thread_hmm[thread_id];
	buffer = data->bsd->sb_file->buffer;
	gc = data->bsd->gc;
	float scores[max_num_hits+1];// = malloc(sizeof(float)* (LIST_STORE_SIZE+1));
	float genome_scores[max_num_hits+1];// = malloc(sizeof(float)*  (LIST_STORE_SIZE+1));


	start = 0;
	stop = 0;
	RUN(get_start_stop_for_threads(num_threads, num_sequences,thread_id,&start,&stop));

	int j,len;
	//int strand = 0;
	char* genomic_sequence = NULL;
	//char* seq = 0;
	
	double sum = 0;
	int c;
	float max,max2,unaligned;
       
	LOG_MSG("thread %d : %d -%d.",thread_id,start,stop);

//	g_int = data->g_int;
	//seq = malloc(sizeof(char) * 100);
	for(i = start; i < stop;i++){
		if(buffer[i]->num_hits){
			//if(ri[i]->identity[0] >= 0.0f   ){
			//hit = 0;
			max = -SCALEINFTY;
			
			sum = prob2scaledprob(0.0);
			for(c = 0; c < buffer[i]->num_hits;c++){
				genomic_sequence = gc[i]->genomic_sequences[c];
				len = gc[i]->g_len[c];

				
				hmm = glocal_forward_log_Y(hmm,buffer[i]->sequence,  genomic_sequence,buffer[i]->len,  len ,c);
				if(!hmm){
					
					DPRINTF2("%s (%d): %d ->%d", buffer[i]->name,i,buffer[i]->start[c],buffer[i]->stop[c] );
					for(len = 0;len <  gc[i]->g_len[c];len++){
						fprintf(stdout,"%d ", gc[i]->genomic_sequences[c][len]);
					}
					fprintf(stdout,"\n");
					//return FAIL;
				}
				
				
				hmm = random_model_calc(hmm,buffer[i]->sequence,  genomic_sequence,buffer[i]->len,  len, prob2scaledprob(1.0),0);
				
				scores[c] = hmm->score;// + ri[i]->priors[hit];
				genome_scores[c] = hmm->random_score ;// ri[i]->priors[c];
			}
			
			hmm = random_model_calc(hmm,buffer[i]->sequence,  genomic_sequence,buffer[i]->len,  len , prob2scaledprob(1.0),1);
			unaligned = hmm->random_score;// + prob2scaledprob(data->param->unaligned_prior);
			
			for(c = 0; c < buffer[i]->num_hits;c++){
				unaligned += genome_scores[c];// + prob2scaledprob(1.0 - scaledprob2prob(genome_scores[c]));
			}
			
			for(c = 0; c < buffer[i]->num_hits;c++){
				sum = prob2scaledprob(1.0f);
				for(j = 0; j < buffer[i]->num_hits;j++){
					if(c == j){
						//read is aligned at this seed position
						sum += scores[j];// + genome_scores[j];
					}else{
						sum += genome_scores[j];// + prob2scaledprob(1.0 - scaledprob2prob(genome_scores[j])) ;
					}
				}
				scores[c] = sum;
			}
			
			// CRITICAL I THINK...
			//if(ri[i]->priors[0] == 0){
				sum = prob2scaledprob(0.0);
			//}else{
			//	sum = unaligned;
			//}
			max2= prob2scaledprob(0.0);
			max = prob2scaledprob(0.0);
			for(c = 0; c < buffer[i]->num_hits;c++){
				//for(c = 0; c <  LIST_STORE_SIZE+1;c++){
				//	if(ri[i]->identity[c] >= 0.0f){
				sum = logsum(sum,scores[c]);
				if(scores[c] >= max){
					max2 = max;
					max = scores[c];
				}else{
					if(scores[c] >= max2){
						max2 = scores[c];
					}
				}
			}
			/* add max prob to thread score for bookkeeping. */
			//fprintf(stdout,"score: best: %f\n",scores[max - sum);
			//fprintf(stdout,"thread %d: %f adding %f\n",thread_id,  data->bsd->thread_forward[thread_id],scores[max] - sum);
			for(c = 0; c < buffer[i]->num_hits;c++){
				data->bsd->thread_forward[thread_id] += scores[c] - sum;
			}
			//hit = 0;
			//if(i < 5){
			
			///}
			max =  prob2scaledprob(scaledprob2prob(max - sum ) -scaledprob2prob(max2 - sum ) );

			
		
			for(c = 0; c < buffer[i]->num_hits;c++){
				//DPRINTF2("%s:%d-%d \n",g_int->chromosome,g_int->start,g_int->stop);
				genomic_sequence = gc[i]->genomic_sequences[c];//  data->sb_file->buffer[i]->genomic_sequences[c];
				len = gc[i]->g_len[c];//   data->sb_file->buffer[i]->g_len[c];
				
				DPRINTF2("len: %d and %d",buffer[i]->len,  len);
				hmm = glocal_backward_log_Y(hmm,buffer[i]->sequence,  genomic_sequence,buffer[i]->len,  len);
				
				hmm = get_prob_log_Y(hmm,buffer[i]->sequence,  genomic_sequence,buffer[i]->len,  len, scores[c] - sum + gc[i]->priors[c]  ,c);// max ,c);//  hmm->score + hmm->score - sum + ri[i]->priors[hit]  );
				//fprintf(stderr,"%s:%d	%f	%f\n",data->sb_file->buffer[i]->name,c,scores[c] - sum,scaledprob2prob(scores[c] - sum) );
			}
			
		}
	}
	//free_genome_interval(g_int);
	//free(genome_scores);
	//free(scores);
	//free(seq);
	return NULL;
ERROR:
	return NULL;
}

void* do_score_alignments_thread_hmm(void *threadarg)
{
	struct thread_data *data;
	struct genome_sequences** gc = NULL;
	struct sam_bam_entry** buffer = NULL;
	struct hmm* hmm = NULL;
	struct genome_interval* g_int = NULL;
	struct pwrite_main* pw = NULL; 
	int thread_id = -1;
	int num_threads = 0;
	int num_sequences = 0;
	int max_num_hits = 0;
	int i;
	int start,stop;
	
	data = (struct thread_data *) threadarg;

	thread_id = data->thread_id;

	num_sequences = data->bsd->sb_file->num_read;
	num_threads = data->bsd->num_threads;
	max_num_hits = data->bsd->sb_file->max_num_hits;
	
	hmm = data->bsd->thread_hmm[thread_id];
	pw = data->bsd->pw;
	buffer = data->bsd->sb_file->buffer;
	g_int = data->bsd->g_int_working[thread_id];
	gc = data->bsd->gc;
	float genome_scores[max_num_hits+1];
	float scores[max_num_hits+1];
	

	start = 0;
	stop = 0;
	
	RUN(get_start_stop_for_threads(num_threads, num_sequences,thread_id,&start,&stop));

	
	int j,c,k,f;
	//int hit = 0;
	//char* seq = 0;
	char* aln = 0;
	int num_scores = 0;
	
	//unsigned long int*  top_scores = 0;
	//float* scores = 0;
	float sum;//genome_sum;
	float unaligned;
	float max;
	int flag;
	
	char* genomic_sequence = NULL;
	int len;
	//top_scores = malloc(sizeof(unsigned long int)* (LIST_STORE_SIZE+1));
	LOG_MSG("Looking at %d-%d.",start,stop);
	for(i = start;i < stop;i++){
	
		//hit = 0;
		num_scores = 0;
		if(!buffer[i]->num_hits){
		//	unaligned_to_sam(ri[i]);
		}else{
		//	hit = 0;
			sum = prob2scaledprob(0.0);
			
			for(c = 0; c < buffer[i]->num_hits;c++){
				
				//DPRINTF2("%s: %d ->%d", data->sb_file->buffer[i]->name,data->sb_file->buffer[i]->start[hit],data->sb_file->buffer[i]->stop[hit] );
				
				genomic_sequence = gc[i]->genomic_sequences[c];// data->sb_file->buffer[i]->genomic_sequences[c];
				len = gc[i]->g_len[c];//    data->sb_file->buffer[i]->g_len[c];
				hmm = glocal_forward_log_Y(hmm,buffer[i]->sequence,  genomic_sequence,buffer[i]->len,  len ,c);
				if(!hmm){
					
					fprintf(stdout,"%s (%d): hit:%d %" PRId64 " ->%" PRId64 "\n", buffer[i]->name,i,c,buffer[i]->start[c],buffer[i]->stop[c] );
				
					for(len = 0;len < buffer[i]->len;len++){
						fprintf(stdout,"%d ", buffer[i]->sequence[len]);
					}
					fprintf(stdout,"\n");
					for(len = 0;len < gc[i]->g_len[c];len++){
						fprintf(stdout,"%d ", gc[i]->genomic_sequences[c][len]);
					}
					fprintf(stdout,"\n");
					//return FAIL;
				}
				
				hmm = random_model_calc(hmm,buffer[i]->sequence,  genomic_sequence,buffer[i]->len,  len, prob2scaledprob(1.0),0);
				
				scores[c] = hmm->score;// + ri[i]->priors[hit];
				genome_scores[c] = hmm->random_score ;// ri[i]->priors[c];
			
			}
			hmm = random_model_calc(hmm,buffer[i]->sequence,  genomic_sequence,buffer[i]->len,  len , prob2scaledprob(1.0),1);
			unaligned = hmm->random_score;// + prob2scaledprob(data->param->unaligned_prior);
			
			for(c = 0; c < buffer[i]->num_hits;c++){
				unaligned += genome_scores[c];// + prob2scaledprob(1.0 - scaledprob2prob(genome_scores[c]));
			}
			
			for(c = 0; c < buffer[i]->num_hits;c++){
				sum = prob2scaledprob(1.0f);
				for(j = 0; j < buffer[i]->num_hits;j++){
				//for(j = 0; j < ri[i]->nhits;j++){
					if(c == j){
						//read is aligned at this seed position
						sum += scores[j];// + genome_scores[j];
					}else{
						sum += genome_scores[j];// + prob2scaledprob(1.0 - scaledprob2prob(genome_scores[j])) ;
					}
				}
				scores[c] = sum;
			}
			
			sum = unaligned;
			//sum = scores[0];
			for(c = 0; c < buffer[i]->num_hits;c++){
				
				sum = logsum(sum,scores[c]);
			}
			unaligned = unaligned - sum;

			max = prob2scaledprob(0.0);
			
			k = 1;
			
			for(c = 0; c < buffer[i]->num_hits;c++){
				scores[c] = scores[c] - sum;
				if(scores[c] > max){
					max = scores[c];
					k = 1;
				}else if(scores[c] == max){
					k++;
				}
			}
			f = 0;
			if(k > 1){
				k = 0;
				
				//k = random_int(k) + 1;
				f = 1;
			}else{
				k = 0;
			}
			
			if(unaligned  > max){
			//	unaligned_to_sam(ri[i]);
			}else{
				
				//if
		//		hit = 0;
				//c =0;
				
				
				//fprintf(stderr,"\n%s	un	%f	%f	%f\n", data->sb_file->buffer[i]->name, scaledprob2prob(unaligned),unaligned,unaligned) ;
				for(j = 0; j < buffer[i]->num_hits;j++){
					flag = 0;
					if(scores[j]  == max){
						if(k){
							if(k != f ){
								flag |= 0x100;
							}
							f++;
						}
					}else if(scaledprob2prob(scores[j]) > (1.0 / (double)LIST_STORE_SIZE)){
						flag |= 0x100;
					}else{
						flag = -1;
					}
					if(flag != -1){
						RUN(get_chr_start_stop(data->bsd->sb_file->si,g_int,buffer[i]->start[j],buffer[i]->stop[j]));
						g_int->start -=ALIGNMENT_FLANKING_LENGTH;
						g_int->stop += ALIGNMENT_FLANKING_LENGTH;
						
						
						//fprintf(stderr,"%s	%d	%f	%f	%d\n",data->sb_file->buffer[i]->name, j,scaledprob2prob(scores[j] ),scores[j] ,flag) ;
						genomic_sequence = gc[i]->genomic_sequences[j];//  data->sb_file->buffer[i]->genomic_sequences[j];
						len = gc[i]->g_len[j];//   data->sb_file->buffer[i]->g_len[j];
						//DPRINTF2("%s ( %d)\n",genomic_sequence,len);
						
						
						/*if(g_int->strand){
						  RUN(reverse_complement_sequence(data->sb_file->buffer[i]->sequence , data->sb_file->buffer[i]->len),"revcomp failed.");
						  RUN(reverse(data->sb_file->buffer[i]->basequal, data->sb_file->buffer[i]->len),"rev failed.");
						  }*/
						
						
						aln = glocal_viterbi_log_Y(hmm,buffer[i]->sequence,  genomic_sequence,buffer[i]->len,  len);
						
						if(!aln){
							fprintf(stderr,"seq: %d hitNO:%d\n",i,j);
							fprintf(stderr,"%s aligned against:\n%s %" PRId64 " -%" PRId64 " %d\n",buffer[i]->name, g_int->chromosome,g_int->g_start,g_int->g_stop,g_int->strand);
						}
						
						
						if(g_int->strand){
							
							aln = reverse_path(aln);
							RUN(reverse_complement_sequence(buffer[i]->sequence , buffer[i]->len));//,"revcomp failed.");
							RUN(reverse(buffer[i]->base_qual, buffer[i]->len));//,"rev failed.");
							flag |= 16;
						}
						if(num_scores == 1){
							if( scaledprob2prob( scores[j]) >= buffer[i]->qual){
								scores[j] = prob2scaledprob(buffer[i]->qual);
							}
						}
						RUN(align_to_sam(pw, g_int, buffer[i], thread_id, aln, flag,scaledprob2prob( scores[j])));
						
						if(g_int->strand){
							
							RUN(reverse_complement_sequence(buffer[i]->sequence , buffer[i]->len));//,"revcomp failed.");
							RUN(reverse(buffer[i]->base_qual, buffer[i]->len));//,"rev failed.");
						}
						MFREE(aln);
						
					}
				}
			}
		}
	}

	pw->flush(pw,thread_id);	
	
	return NULL;
ERROR:
	return NULL;
}


int get_start_stop_for_threads(const int num_threads, const int num_seq, const int id, int *start,int *stop)
{

	int interval = 0;

	if(num_threads > num_seq){
		if(id > num_seq){
			*start = 0;
			*stop = 0;
		}else{
			*start = id;
			*stop = id+1;
		}
		return OK;
	}
	
	interval = num_seq / num_threads;
	*start = id* interval;
	*stop =  id* interval + interval;
	LOG_MSG("Trying to spread work: %d total_threads %d.",id, num_threads);
	if(id +1 == num_threads){
		*stop = num_seq;
	}
		
	return OK;
}




char* glocal_viterbi_log_Y(struct hmm* hmm,char* a, char* b, int n,int m)
{
	int i,j,c;
	
	unsigned short int tmp_tb;
	float** M = hmm->fM;
	//float** C = hmm->fC;
	float** X = hmm->fX;
	float** Y = hmm->fY;
	char* seqa = a -1;
	char* seqb = b -1;
	float* prob = hmm->prob;
	unsigned short int** t = hmm->traceback;
	char* alignment =  malloc(sizeof( char)* (n+m+2));
	
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
}

struct hmm* glocal_forward_log_Y(struct hmm* hmm, char* a,  char* b, int n,int m,int frame )
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
	return hmm;
ERROR:
	return NULL;
}

struct hmm* glocal_backward_log_Y(struct hmm* hmm,char* a, char* b, int n,int m)
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
	
	return hmm;
}



struct hmm* get_prob_log_Y(struct hmm* hmm,char* a, char* b, int n,int m,float frac,int frame )
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
	
	return hmm;
}



struct hmm* random_model_calc(struct hmm* hmm, char* rseq, char* gseq,int rlen, int glen,float frac,int mode)
{
	int i,code;
	hmm->random_score = prob2scaledprob(1.0);
	if(mode == 0){
		for(i = 0; i < glen-1;i++){
			if(gseq[i] == NCODE){
				code = 20;
				
			}else{
				code = (gseq[i] & 0x3) *5;
			}
			if(gseq[i+1] == NCODE){
				code += 4;
			}else{
				code+= (gseq[i+1] & 0x3);
			}
			hmm->random_score += hmm->random_genome_model[code];
			hmm->random_genome_model_t[code] =  logsum(hmm->random_genome_model_t[code],  frac);
		}
	}
	if(mode == 1){
		for(i = 0; i < rlen-1;i++){
			if(rseq[i] == NCODE){
				code = 20;
				
			}else{
				code = (rseq[i] & 0x3) *5;
			}
			if(rseq[i+1] == NCODE){
				code += 4;
			}else{
				code+= (rseq[i+1] & 0x3);
			}
			hmm->random_score += hmm->random_read_model[code];
			hmm->random_read_model_t[code] =  logsum(hmm->random_read_model_t[code], frac );
		}
	}
	
	return hmm;
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

int align_to_sam(struct pwrite_main* pw,struct genome_interval* g_int,struct sam_bam_entry* entry,int id, char* aln,unsigned flag,float score)
{
	char nuc[] = "ACGTN";
	char cigarline[128];
	char mdline[128];
	char pseq[500];
	//unsigned char* qual = 0;
	
	//unsigned char noqual[] = "*";
	int mismatches = 0;
	int i,j,c,g,t;
	unsigned int add = 0;
	c = 0;
	int state = 0;
	mdline[0] = 0;
	c = aln[0];
	j = 0;
	t = 1;
	
	
	while((aln[c]  >> 4) == 5 ){
	//	fprintf(stderr,"%d	%d	len: %d\n",((aln[c]   >> 4)& 0xF ),(aln[c]  & 0xF),entry->len);
		c--;
	}
	for(i = c;i >= 1;i--){
	//	fprintf(stderr,"%d	%d	%d\n",((aln[i]   >> 4)& 0xF ),(aln[i]  & 0xF),t);
		if((aln[i]  >> 4) != 5 && (aln[i]  & 0xF) !=5){
			if((aln[i]  >> 4) != (aln[i]  & 0xF)){
				g = (int)strlen(mdline);
				sprintf(mdline+g,"%d%c",j,nuc[(aln[i]  & 0xF)]);
				if((aln[i]  & 0xF) > 4 ){
					
					fprintf(stderr,"ERROR:%d	%d\n",aln[i] ,aln[i]  & 0xF) ;
				}
	//			fprintf(stderr,"%d	%c\n",j,nuc[(aln[i]  & 0xF)]);
				j = 0;
			}else{
				j++;
			}
			
			if(t == entry->len){
				i = -1;
				break;
			}
			t++;
			state = 0;
		}else if((aln[i]  >> 4) != 5 && (aln[i]  & 0xF) == 5){
			
			
			
			if(t == entry->len){
				i = -1;
				break;
			}
			t++;
			
	//		fprintf(stderr,"%c	\n",nuc[(aln[i]  & 0xF)]);
			//j = 0;
			state = 0;
		}else if((aln[i]  >> 4) == 5  && (aln[i]  & 0xF) != 5){
			if(state == 0){
				g = (int)strlen(mdline);
				sprintf(mdline+g,"%d^",j);
	//			fprintf(stderr,"%d	\n",j);
				state = 1;
			}
			g = (int)strlen(mdline);
			sprintf(mdline+g,"%c",nuc[(aln[i]  & 0xF )]);
			j = 0;
			//state = 0;
		}
	}
	
	if(j){
		g = (int)strlen(mdline);
		sprintf(mdline+g,"%d",j);
	//		fprintf(stderr,"%d	\n",j);
	}
	
	g = (int)strlen(mdline);
	if(!isdigit((int)  mdline[g-1])){
		sprintf(mdline+g,"%d",0);
	}
	
	//fprintf(stderr,"LINE:%s\n\n",mdline);
	
	
	i =  aln[0];
	t = 0;
	cigarline[0] = 0;
	while((aln[i]  >> 4) == 5){
		i--;
		add++;
	}
	
	while(i >= 1){
		j = 0;
		if((aln[i]  >> 4) != 5 && (aln[i]  & 0xF) != 5){
			while ( (aln[i]  >> 4) != 5 && (aln[i]  & 0xF) != 5 && i >= 1) {
				if((aln[i]  >> 4) != (aln[i]  & 0xF) ){
					mismatches++;
				}
				i--;
				j++;
				t++;
				
			}
			c = (int)strlen(cigarline);
			sprintf(cigarline+c,"%dM",j);
			if(t ==entry->len){
				break;
			}
			
		}else if((aln[i]  >> 4) != 5 && (aln[i]  & 0xF) == 5 ){
			while ((aln[i]  >> 4) != 5 && (aln[i]  & 0xF) == 5&& i >= 1) {
				i--;
				j++;
				t++;
				mismatches++;
			}
			c = (int)strlen(cigarline);
			sprintf(cigarline+c,"%dI",j);
			if(t ==entry->len){
				break;
			}
		}else if((aln[i]  >> 4) == 5 && (aln[i]  & 0xF) != 5 ){
			while ((aln[i]  >> 4) == 5 && (aln[i]  & 0xF) != 5&& i >= 1) {
				i--;
				j++;
				mismatches++;
			}
			c = (int)strlen(cigarline);
			sprintf(cigarline+c,"%dD",j);
		}
		//fprintf(stderr,"%d	%s	%d	%d	%d\n",i,cigarline,aln[i]>>4, aln[i] & 0xF,add);
	}
	//fprintf(stderr,"CIGAR ::: %s\n",cigarline);
	
	for(i = 0;i < entry->len;i++){
		if(entry->sequence[i] == 12){
			pseq[i] = 'N';
		}else{
			pseq[i] =  "ACGTN"[(int)entry->sequence[i]];//   nuc[(int)ri->seq[i] & 3 ];
		}
	}
	pseq[entry->len]  = 0;
	
	if(1.0f - score < 0.0001f){
		score = 0.0001f;
	}else{
		score = 1.0f - score;
	}
	
	//get_chr_name (seq_info,(unsigned int) pos),
	//get_chr_pos (seq_info,(unsigned int) pos),
	
	if(flag & 0x100){
		RUN(pw->write(pw,id,"%s\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%s\tNM:i:%d\tMD:Z:%s\n",
			      entry->name,//1 QNAME Query NAME of the read or the read pair
			      flag,//2 FLAG bitwise FLAG (pairing, strand, mate strand, etc.)
			      g_int->chromosome,
			      g_int->start+add +1,
			      //get_chr_name (db,(unsigned int) pos + add + 1), //3 RNAME Reference sequence NAME
			      //get_chr_pos (db,(unsigned int) pos+add +1) ,//4 POS 1-based leftmost POSition of clipped alignment
			      (int)((-10.0f * log10(score ))+0.5f),//5 MAPQ MAPping Quality (Phred-scaled)
			      cigarline,//6 CIGAR extended CIGAR string (operations: MIDNSHP)
			      "*",//7 MRNM Mate Reference NaMe (‘=’ if same as RNAME)
			      0,//8 MPOS 1-based leftmost Mate POSition
			      0,//9 ISIZE inferred Insert SIZE
			      "*" ,//10 SEQ query SEQuence on the same strand as the reference
			      "*",//11 QUAL query QUALity (ASCII-33=Phred base quality)
			      mismatches,
			      mdline
			    ));
	}else{
		RUN(pw->write(pw,id,"%s\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%s\tNM:i:%d\tMD:Z:%s\n",
			      entry->name,//1 QNAME Query NAME of the read or the read pair
			      flag,//2 FLAG bitwise FLAG (pairing, strand, mate strand, etc.)
			      g_int->chromosome,
			      g_int->start+add +1,
			      (int)((-10.0f * log10(score))+0.5f),//5 MAPQ MAPping Quality (Phred-scaled)
			      cigarline,//6 CIGAR extended CIGAR string (operations: MIDNSHP)
			      "*",//7 MRNM Mate Reference NaMe (‘=’ if same as RNAME)
			      0,//8 MPOS 1-based leftmost Mate POSition
			      0,//9 ISIZE inferred Insert SIZE
			      pseq ,//10 SEQ query SEQuence on the same strand as the reference
			      entry->base_qual,//11 QUAL query QUALity (ASCII-33=Phred base quality)
			      mismatches,
			      mdline
			    ));
	}

	
	return OK;
ERROR:
	return FAIL;
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


char* reverse_path(char* org )
{
	char* new = NULL;
	
	int i,c;
	int nuc[5] = {3,2,1,0,4};
	MMALLOC(new, sizeof(char)* (org[0]+1));
	
	new[0] = org[0];
	c = 1;
	for (i = org[0]; i > 0; i--){
		if(( org[i] >> 4) < 5){
			if((org[i] & 0xF) < 5){
				new[c] = (nuc[org[i] >> 4] << 4) |  nuc[org[i] & 0xF];
			}else{
				new[c] = (nuc[org[i] >> 4] << 4) |  (org[i] & 0xF);
			}
		}else{
			if((org[i] & 0xF) < 5){
				new[c] = (org[i] & 0xF0) |  nuc[org[i] & 0xF];
			}else{
				new[c] =org[i];
			}
		}
		//new[c] = org[i];
		c++;
	}
	MFREE(org);
	return new;
ERROR:
	return NULL;
}

int reverse(uint8_t* p,int len)
{
	int c, i, j;

	if(p[0] == '*'){
		return OK;
	}
	
	for (i = 0, j = len - 1; i < j; i++, j--)
	{
		c = p[i];
		p[i] = p[j];
		p[j] = c;
	}
	return OK;
}

int reverse_complement_sequence(char* p,int len)
{
	int c, i, j;
	
	for (i = 0, j = len - 1; i < j; i++, j--)
	{
		c = p[i];
		p[i] = p[j];
		p[j] = c;
	}
	for (i = 0; i <  len; i++){
		
		switch (p[i]) {
			case 0:
				
				p[i] = 3;
				break;
			case 1:
				p[i] = 2;
				break;
			case 2:
				p[i] = 1;
				break;
			case 3:
				p[i] = 0;
				break;
			default:
				break;
		}
		
	}
	return OK;
}

int convert_buffer_ACGT_to_0123(struct sam_bam_file* sb_file)

{
	int i;
	for(i = 0; i < sb_file->num_read;i++){
		RUN(ACGT_to_0123(sb_file->buffer[i]->sequence,&sb_file->buffer[i]->len));
	}
	return OK;
ERROR:
	return FAIL;
}


int ACGT_to_0123(char* seq,int* len)
{
	int i = 0;
	int c = 0;
	while(seq[i] !=0){
		switch(seq[i]){
			case 'A':
			case 'a':
				seq[c] = 0;
				c++;
				break;
			case 'C':
			case 'c':
				seq[c] = 1;
				c++;
				break;
			case 'G':
			case 'g':
				seq[c] = 2;
				c++;
				break;
			case 'T':
			case 't':
				seq[c] = 3;
				c++;
				break;
			default:
				seq[c] = NCODE;
				c++;
				//ERROR_MSG("Non ACGT letter in sequence:%s.",seq);
				break;
		}
		i++;
	}
	seq[c] = 0;
	*len = c;
	return OK;
}


int add_genome_sequences(struct shared_data* bsd)
{
//	faidx_t*  index = NULL;
	struct sam_bam_file* sb_file = NULL;
	struct genome_sequences** gc = NULL;
	faidx_t*  index = NULL;
	
	struct genome_interval* g_int = NULL;

	int i,j,len,c;

	
	sb_file = bsd->sb_file;
	gc = bsd->gc;
	index = bsd->index;
	
	//int len_b,len_c;
//	RUNP(index = get_faidx(genome));
	RUNP(g_int = init_genome_interval(0,0,0));
	for(i = 0; i < sb_file->num_read;i++){
		for(j = 0; j < sb_file->buffer[i]->num_hits;j++){
			DPRINTF2("NEWHIT:");
			DPRINTF2("%s: %d ->%d", sb_file->buffer[i]->name,sb_file->buffer[i]->start[j],sb_file->buffer[i]->stop[j] );
			
			RUN(get_chr_start_stop(sb_file->si  ,g_int,sb_file->buffer[i]->start[j], sb_file->buffer[i]->stop[j]));
			g_int->start -=ALIGNMENT_FLANKING_LENGTH;
			g_int->stop += ALIGNMENT_FLANKING_LENGTH;
			DPRINTF2("%s:%d-%d ",g_int->chromosome,g_int->start,g_int->stop);
//			sb_file->buffer[i]->genomic_sequences[j] = get_sequence(index,g_int, &len_a);

			gc[i]->genomic_sequences[j] = get_sequence(index,g_int);
			
			len = (int) strlen(gc[i]->genomic_sequences[j]);
			//len_b = len;
//			sb_file->buffer[i]->g_len[j] = len;
			gc[i]->g_len[j] = len;
			//DPRINTF2("%s	(len:%d)",sb_file->buffer[i]->genomic_sequences[j],sb_file->buffer[i]->g_len[j]);

			RUN(ACGT_to_0123(gc[i]->genomic_sequences[j] ,&len));
			
			
			//len_c = len;
			
//			ASSERT(len_a == len_b,"%d %d %d %d %d\n",i,j, len_a,len_b,len_c);
//			ASSERT(len_b == len_c,"%d %d %d %d %d\n",i,j, len_a,len_b,len_c);
//			ASSERT(len_a == len_c,"%d %d %d %d %d\n",i,j, len_a,len_b,len_c);
			for(c =0; c < len;c++){
				ASSERT(gc[i]->genomic_sequences[j][c] < 5, "%d %d pos:%d  nuc:%d  ", i,j, c, gc[i]->genomic_sequences[j][c] ); 
			}		
		}
	}
	free_genome_interval(g_int);
//	free_faidx(index);
	return OK;
ERROR:
	if(g_int){
		free_genome_interval(g_int);
	}
//	if(index){
//		free_faidx(index);
//	}
	return FAIL;
}

/*
int remove_genome_sequences(struct sam_bam_file* sb_file)
{
	int i,j;
	for(i = 0; i < sb_file->num_read;i++){
		for(j = 0; j < sb_file->buffer[i]->num_hits;j++){
			
			MFREE(sb_file->buffer[i]->genomic_sequences[j]);//
			sb_file->buffer[i]->genomic_sequences[j]= NULL;
			sb_file->buffer[i]->g_len[j]  = 0;
		}
		sb_file->buffer[i]->num_hits = 0;
	}
	
	return OK;
ERROR:
	return FAIL;
}
*/
		
struct genome_sequences** init_genome_sequences(int num, int num_maxhits)
{
	struct genome_sequences** gc = NULL;
	int i,j;

	MMALLOC(gc,sizeof(struct genome_sequences*) * num);
	for(i = 0; i < num;i++){
		gc[i] = NULL;
		MMALLOC(gc[i], sizeof(struct genome_sequences));
		gc[i]->genomic_sequences = NULL;
		gc[i]->priors = NULL;
		gc[i]->g_len = NULL;
		MMALLOC(gc[i]->g_len,sizeof(int) * num_maxhits);
		MMALLOC(gc[i]->genomic_sequences,sizeof(char*) * num_maxhits);
		MMALLOC(gc[i]->priors,sizeof(float) * num_maxhits);
		for (j = 0; j < num_maxhits; j++) {
			gc[i]->g_len[j] = 0;
			gc[i]->genomic_sequences[j] = NULL;
			gc[i]->priors[j] = 0.0f;
			
		}
	}
	return gc;
ERROR:
	free_genome_sequences(gc,num,num_maxhits);			    
	return NULL;
	 
}


int clear_genome_sequences(struct genome_sequences** gc, int num,int num_maxhits)
{
	int i,j;
	for(i = 0;i< num;i++){
		if(gc[i]){
			for(j = 0; j < num_maxhits;j++){
				if(gc[i]->genomic_sequences[j]){
					MFREE(gc[i]->genomic_sequences[j]);
					gc[i]->genomic_sequences[j] = NULL;
				}
				gc[i]->g_len[j] = 0;
				       
			}
		}
				      
	}
	
	
	return OK;		
}


void free_genome_sequences(struct genome_sequences** gc, int num,int num_maxhits)
{
	int i,j;
	if(gc){
		for(i = 0;i< num;i++){
			if(gc[i]){
				for(j = 0; j < num_maxhits;j++){
					if(gc[i]->genomic_sequences[j]){
						MFREE(gc[i]->genomic_sequences[j]);
					}
				}
				MFREE(gc[i]->priors);
				MFREE(gc[i]->genomic_sequences);
				MFREE(gc[i]->g_len);
				MFREE(gc[i]);
			}
				      
		}
		MFREE(gc);
	}
}



int write_sam_header(struct sam_bam_file* sb_file, FILE* out)
{
	int i;
	fprintf(out,"@HD\tVN:1.3\tSO:coordinate\n");
	for(i = 0; i < sb_file->header->n_targets;i++){
		fprintf(out,"@SQ\tSN:%s\tLN:%d\n", sb_file->header->target_name[i],(int)sb_file->header->target_len[i]);
	}
	return OK;
}
	

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
	bsd->pseudo_counts = 1.0f;
	/* assignment  */
	bsd->free = free_shared_data;
	bsd->buffer_size = buffer_size;
	bsd->num_maxhits = param->num_maxhits;
	bsd->param = param;
	bsd->num_threads = param->num_threads;

	/* initialize mutex  */
	pthread_mutex_init(&avail_mtx,NULL);
	
	/* allocation  */
	if((bsd->pool = thr_pool_create(bsd->num_threads+1, bsd->num_threads+1, 0, 0))  ==NULL) ERROR_MSG("Creating pool thread failed.");
	RUNP(bsd->gc = init_genome_sequences(bsd->buffer_size, bsd->num_maxhits));
	RUNP(bsd->index = get_faidx(bsd->param->genome));
	RUNP(bsd->pw = init_pwrite_main(param->out_file,param->num_threads,BUFFER_P_WRITE_SIZE));

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
};

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
		pthread_mutex_destroy(&avail_mtx);
		MFREE(bsd);
	}
}
