#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <string.h>

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif


#include "tldevel.h"
#include "rbtree.h"
#include "htsglue.h"
#include  <ctype.h>
#include "rtr.h"
#include "delve.h"

#define OPT_NTHREAD 1
#define OPT_OUT 2

struct parameters{
	char** infiles;
	char* genome;
	char* aln_infile;
	char* hmm_file;
	char* out_prefix;
	int num_infiles;
	int num_threads;
};

struct rtr_data* build_rtree(struct sam_bam_file* sb_file );

int add_genome_sequences(struct sam_bam_file* sb_file,struct genome_sequences** gc, char* genome);

//int add_genome_sequences(struct sam_bam_file* sb_file,char* genome);
//int remove_genome_sequences(struct sam_bam_file* sb_file);

int run_delve(struct parameters* param);


struct genome_sequences** init_genome_sequences(int num, int num_maxhits);
int clear_genome_sequences(struct genome_sequences** gc, int num,int num_maxhits);
void free_genome_sequences(struct genome_sequences** gc, int num,int num_maxhits);


int main (int argc,char *argv[]) 
{
	struct parameters* param = NULL;
	int i,c;
	
	//const char description[] = "Aligns HMMs to new sequences.";
	//const char usage[] = "<genome> <seed aln> ";


	tlog.echo_build_config();
       
	MMALLOC(param, sizeof(struct parameters));

	param->infiles = NULL;
	param->genome = NULL;
	param->aln_infile = NULL;
	param->hmm_file = NULL;
	param->out_prefix = NULL;
	param->num_infiles = 0;
	param->num_threads = 4;
	
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
				param->out_prefix = optarg;
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
		for(i = 0; i < param->num_infiles;i++){
			MFREE(param->infiles[i]);
		}
		MFREE(param->infiles);
		MFREE(param->out_prefix);
		MFREE(param);
	}
	return EXIT_SUCCESS;
ERROR:
	fprintf(stdout,"\n  Try run with  --help.\n\n");
	if(param){
		for(i = 0; i < param->num_infiles;i++){
			MFREE(param->infiles[i]);
		}
		MFREE(param->infiles);
		MFREE(param->out_prefix);
		MFREE(param);
	}
	return EXIT_FAILURE;
}




int run_delve(struct parameters* param)
{
	struct sam_bam_file* sb_file = NULL;
	struct rtr_data* rtree = NULL;
	struct hmm* hmm = NULL;

	struct genome_sequences** gc = NULL; 
	int i,j,iter;


	int buffer_size = 1000000;
	init_logsum();
	
	faidx_t*  index = NULL;
	RUNP(sb_file = open_SAMBAMfile(param->aln_infile,buffer_size,0,0));

	/* Alloc space for genomic sequences.. */

	RUNP(gc = init_genome_sequences(buffer_size,sb_file->max_num_hits));
	
	RUNP(index = get_faidx(param->genome));
	RUNP(hmm = init_hmm(100, 100));
	RUN(add_pseudo_count(hmm, 100.0f));
	RUN(re_estimate(hmm));

	tlog.log_message("Opening alignment file...");
	
	RUN(read_SAMBAM_chunk(sb_file,1,0));
	
	RUN(add_genome_sequences(sb_file,gc,param->genome));
	
	//for(i = 0; i < sb_file->num_read;i++){
	//	RUN(ACGT_to_0123(sb_file->buffer[i]->sequence,&sb_file->buffer[i]->len  ));
	//	for(j = 0; j < sb_file->buffer[i]->num_hits;j++){
	//		DPRINTF2("%d %d %p",i,j,sb_file->buffer[i]->genomic_sequences[j] );
	//	}
	//}
	tlog.log_message("Quantifying read depth at putatitve mapping locations.");
	RUNP(rtree =  build_rtree(sb_file) );
	exit(0);
	tlog.log_message("Estimating Random model...");
	
	
	RUN(run_pHMM(hmm, sb_file ,gc,index,param->num_threads ,sb_file->num_read,RUN_RANDOM_HMM,stdout));
	RUN(re_estimate_random(hmm));
	
	tlog.log_message("done");

	
	for(i = 0; i < sb_file->num_read;i++){
		//RUN(ACGT_to_0123(sb_file->buffer[i]->sequence,&sb_file->buffer[i]->len  ),"ACGT to 0123 failed");
		for(j = 0; j < sb_file->buffer[i]->num_hits;j++){
//			DPRINTF2("%d %d %p",i,j,sb_file->buffer[i]->genomic_sequences[j] );
		}
	}
	//total_numseq = 0;
	tlog.log_message("Estimating Sequence model...");
	for(iter = 0; iter < 10;iter++){
		tlog.log_message("%d iteration...",iter);
		RUN(run_pHMM(hmm,sb_file ,gc,  index,param->num_threads,sb_file->num_read,RUN_FULL_HMM,stdout));

	
		RUN(re_estimate(hmm));
		RUN(add_pseudo_count(hmm, 100.0f));
		RUN(re_estimate_random(hmm));
	}
	tlog.log_message("done");
	RUN(clear_genome_sequences(gc,buffer_size,sb_file->max_num_hits  ));
//	RUN(remove_genome_sequences(sb_file));
	//Start from beginning....
	RUN(close_SAMBAMfile(sb_file));
	RUNP(sb_file = open_SAMBAMfile(param->aln_infile,1000000,0,0));
	
	RUN(echo_header(sb_file));
	tlog.log_message("Generating alignments...");

	while(1){
		RUN(read_SAMBAM_chunk(sb_file,1,0));
		//DPRINTF2("read:%d",sb_file->num_read );
		//RUN(read_sam_bam_chunk(infile,data->si, buffer,1,&numseq),"Read sam/bam chunk in thread %d failed", data->threadID);
		//if((status = read_sam_bam_chunk(infile,data->si, buffer,1,&numseq)) != kslOK)  exit(status);
		if(!sb_file->num_read){
			break;
		}
		tlog.log_message("read:%d",sb_file->num_read);
		RUN(add_genome_sequences(sb_file,gc,param->genome));
		for(i = 0; i < sb_file->num_read;i++){
			RUN(ACGT_to_0123(sb_file->buffer[i]->sequence,&sb_file->buffer[i]->len  ));
		}
		
		RUN(run_pHMM(hmm,sb_file ,gc,  index,param->num_threads,sb_file->num_read,RUN_SCORE_HMM,stdout));
		RUN(clear_genome_sequences(gc,buffer_size,sb_file->max_num_hits  ));
        }
	tlog.log_message("done");
	
	//rewind(file);
	
	//param->num_query = 1000000;
	/*numseq = read_sam_chunk(ri,db,file);
	for(i = 0;i < numseq;i++){
		if(ri[i]->strand[0] != 0){
			ri[i]->seq = reverse_complement(ri[i]->seq,ri[i]->len);
		}
		ri[i]->seq = transform_16to5_sequence(ri[i]->seq,ri[i]->len);
	}
	
	hmm = run_pHMM(hmm,ri, db,num_threads,numseq,RUN_RANDOM_HMM,out );
	hmm = re_estimate_random(hmm);
	fprintf(stderr,"	done.\n");
	*/
	free_faidx(index);
	free_hmm(hmm);
	RUN(close_SAMBAMfile(sb_file));

	
	return OK;
ERROR:
	return FAIL;
}


struct rtr_data* build_rtree(struct sam_bam_file* sb_file )
{
	struct rtr_data* rtree = NULL;

	int64_t* val = NULL;
	int i,j;
	int32_t id = 0;

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
		tlog.log_message("read:%d",sb_file->num_read);
		for (i = 0; i < sb_file->num_read; i++) {
			for (j=0; j < sb_file->buffer[i]->num_hits ; j++) {
				val[0] = sb_file->buffer[i]->start[j];
				val[1] = sb_file->buffer[i]->stop[j];
				rtree->insert(rtree,val,id,1,1);
			}
			
		}
        }
	MFREE(val);
	RUN(rtree->flatten_rtree(rtree));
	RUN(rtree->print_rtree(rtree, rtree->root));
	for (i = 0; i < rtree->stats_num_interval; i++) {
//		fprintf(stdout,"%d\n", rtree->flat_interval[i]->count);
	}
	return rtree;
ERROR:
	MFREE(val);
	if(rtree){
		rtree->free(rtree);
	}
	return NULL;
}

int run_pHMM(struct hmm* localhmm,struct sam_bam_file* sb_file ,struct genome_sequences** gc, faidx_t*  index,int num_threads ,int size, int mode, FILE* fout)
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
		thread_data_array[t].hmm = init_hmm(100,100);
		thread_data_array[t].start = t*interval;
		thread_data_array[t].end = t*interval + interval;
		thread_data_array[t].fout = fout;
	}
	thread_data_array[num_threads - 1].gc = gc;
//	thread_data_array[num_threads - 1].db = db;
	thread_data_array[num_threads - 1].sb_file = sb_file;
	thread_data_array[num_threads - 1].index = index;
	thread_data_array[num_threads - 1].hmm = init_hmm(100, 100);
	thread_data_array[num_threads - 1].start = t*interval;
	thread_data_array[num_threads - 1].end = size;
	thread_data_array[num_threads - 1].fout = fout;
	
	
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
	}
	
	free(thread_data_array);
	return OK;
ERROR:
	return FAIL;
}

void* do_baum_welch_thread(void *threadarg)
{
	struct thread_data *data;
	data = (struct thread_data *) threadarg;
	
	//struct read_info** ri = data->ri;
	struct hmm* hmm = data->hmm;
	//struct db* db = data->db;
	const int start = data->start;
	const int end = data->end;
	int i,j,len;
	//int strand = 0;
	char* genomic_sequence = NULL;
	char* seq = 0;
	
	double sum = 0;
	int c;
	float* scores = malloc(sizeof(float)* (LIST_STORE_SIZE+1));
	float* genome_scores = malloc(sizeof(float)*  (LIST_STORE_SIZE+1));
	float max,max2,unaligned;
	struct genome_interval* g_int = NULL;
	RUNP(g_int = init_genome_interval(0,0,0));
	
	seq = malloc(sizeof(char) * 100);
	c = LIST_STORE_SIZE;
	for(i = start; i < end;i++){
		if(data->sb_file->buffer[i]->num_hits){
			//if(ri[i]->identity[0] >= 0.0f   ){
			//hit = 0;
			max = -SCALEINFTY;
			
			sum = prob2scaledprob(0.0);
			for(c = 0; c < data->sb_file->buffer[i]->num_hits;c++){
				
				//DPRINTF2("%s: %d ->%d", data->sb_file->buffer[i]->name,data->sb_file->buffer[i]->start[hit],data->sb_file->buffer[i]->stop[hit] );
				
			//DPRINTF2("%s:%d-%d \n",g_int->chromosome,g_int->start,g_int->stop);
				//		DPRINTF2("%d hit %d	%p\n",i,c, data->sb_file->buffer[i]->genomic_sequences[c]);
				//genomic_sequence = data->sb_file->buffer[i]->sequence[c];// get_sequence(data->index,g_int);
				//len =   data->sb_file->buffer[i]->len[c];

				genomic_sequence = data->gc[i]->genomic_sequences[c];
				len = data->gc[i]->g_len[c];

				
				hmm = glocal_forward_log_Y(hmm,data->sb_file->buffer[i]->sequence,  genomic_sequence,data->sb_file->buffer[i]->len,  len ,c);
				if(!hmm){
					
					DPRINTF2("%s (%d): %d ->%d", data->sb_file->buffer[i]->name,i,data->sb_file->buffer[i]->start[c],data->sb_file->buffer[i]->stop[c] );
					for(len = 0;len <  data->gc[i]->g_len[c];len++){
						fprintf(stdout,"%d ", data->gc[i]->genomic_sequences[c][len]);
					}
					fprintf(stdout,"\n");
					//return FAIL;
				}
				
				
				hmm = random_model_calc(hmm,data->sb_file->buffer[i]->sequence,  genomic_sequence,data->sb_file->buffer[i]->len,  len, prob2scaledprob(1.0),0);
				
				scores[c] = hmm->score;// + ri[i]->priors[hit];
				genome_scores[c] = hmm->random_score ;// ri[i]->priors[c];
			}
			
			hmm = random_model_calc(hmm,data->sb_file->buffer[i]->sequence,  genomic_sequence,data->sb_file->buffer[i]->len,  len , prob2scaledprob(1.0),1);
			unaligned = hmm->random_score;// + prob2scaledprob(data->param->unaligned_prior);
			
			for(c = 0; c < data->sb_file->buffer[i]->num_hits;c++){
				unaligned += genome_scores[c];// + prob2scaledprob(1.0 - scaledprob2prob(genome_scores[c]));
			}
			
			for(c = 0; c < data->sb_file->buffer[i]->num_hits;c++){
				sum = prob2scaledprob(1.0f);
				for(j = 0; j < data->sb_file->buffer[i]->num_hits;j++){
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
			for(c = 0; c < data->sb_file->buffer[i]->num_hits;c++){
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
			
			//hit = 0;
			//if(i < 5){
			
			///}
			max =  prob2scaledprob(scaledprob2prob(max - sum ) -scaledprob2prob(max2 - sum ) );
			
		
			for(c = 0; c < data->sb_file->buffer[i]->num_hits;c++){
				//DPRINTF2("%s:%d-%d \n",g_int->chromosome,g_int->start,g_int->stop);
				genomic_sequence = data->gc[i]->genomic_sequences[c];//  data->sb_file->buffer[i]->genomic_sequences[c];
				len = data->gc[i]->g_len[c];//   data->sb_file->buffer[i]->g_len[c];
				
				DPRINTF2("len: %d and %d",data->sb_file->buffer[i]->len,  len);
				hmm = glocal_backward_log_Y(hmm,data->sb_file->buffer[i]->sequence,  genomic_sequence,data->sb_file->buffer[i]->len,  len);
				
				hmm = get_prob_log_Y(hmm,data->sb_file->buffer[i]->sequence,  genomic_sequence,data->sb_file->buffer[i]->len,  len, scores[c] - sum  ,c);// max ,c);//  hmm->score + hmm->score - sum + ri[i]->priors[hit]  );
				//fprintf(stderr,"%s:%d	%f	%f\n",data->sb_file->buffer[i]->name,c,scores[c] - sum,scaledprob2prob(scores[c] - sum) );
			}
			
		}
	}
	free(genome_scores);
	free(scores);
	free(seq);
	pthread_exit((void *) 0);
ERROR:
	pthread_exit((void *) 0);
}



void* do_baum_welch_thread_random_model(void *threadarg)
{
	struct thread_data *data;
	data = (struct thread_data *) threadarg;
	struct genome_interval* g_int = NULL;

	//struct read_info** ri = data->ri;
	struct hmm* hmm = data->hmm;
	//struct db* db = data->db;
	const int start = data->start;
	const int end = data->end;
	int i,hit;
	int j;
	//int strand = 0;
	
	char* seq = 0;
	char* genomic_sequence = NULL;
	int len;
	
	seq = malloc(sizeof(char) * 100);
	
	RUNP(g_int = init_genome_interval(0,0,0));
	
	for(i = start; i < end;i++){
		
		for(hit = 0; hit < data->sb_file->buffer[i]->num_hits;hit++){
			
			//DPRINTF2("%s: %d ->%d", data->sb_file->buffer[i]->name,data->sb_file->buffer[i]->start[hit],data->sb_file->buffer[i]->stop[hit] );
			
			//RUN(get_chr_start_stop(data->sb_file,i,hit,g_int),"get_chr_start_stop failed");
			//g_int->start -=ALIGNMENT_FLANKING_LENGTH;
			//g_int->stop += ALIGNMENT_FLANKING_LENGTH;
			//DPRINTF2("%s:%d-%d \n",g_int->chromosome,g_int->start,g_int->stop);
			genomic_sequence = data->gc[i]->genomic_sequences[hit];//  data->sb_file->buffer[i]->genomic_sequences[hit];
			len = data->gc[i]->g_len[hit];//  data->sb_file->buffer[i]->g_len[hit];
			hmm = random_model_calc(hmm,data->sb_file->buffer[i]->sequence,  genomic_sequence,data->sb_file->buffer[i]->len,  len, prob2scaledprob(1.0),0);
			

		}
		hmm = random_model_calc(hmm,data->sb_file->buffer[i]->sequence,  genomic_sequence,data->sb_file->buffer[i]->len,  len, prob2scaledprob(1.0),1);
		//		ri[i]->random_read_score = hmm->random_score;
	}
	free(seq);
	free_genome_interval(g_int);

	pthread_exit((void *) 0);
ERROR:
	pthread_exit((void *) 0);
}



void* do_score_alignments_thread_hmm(void *threadarg)
{
	struct thread_data *data;
	data = (struct thread_data *) threadarg;
	struct genome_interval* g_int = NULL;
	//struct read_info** ri = data->ri;
	struct hmm* hmm = data->hmm;
	//struct db* db = data->db;
	const int start = data->start;
	const int end = data->end;
	
	int i,j,c,k,f;
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
	
	
	//float* random_scores = malloc(sizeof(float)*  (LIST_STORE_SIZE+1));
	//float* genome_prior = malloc(sizeof(float)*  (LIST_STORE_SIZE+1));
	//float* genome_scores =  malloc(sizeof(float)*  (LIST_STORE_SIZE+1));
	RUNP(g_int = init_genome_interval(0,0,0));
	//seq = malloc(sizeof(char) * 100);
	//scores = malloc(sizeof(float)* LIST_STORE_SIZE);
	
	float genome_scores[LIST_STORE_SIZE+1];
	float scores[LIST_STORE_SIZE];
	//top_scores = malloc(sizeof(unsigned long int)* (LIST_STORE_SIZE+1));
	
	for(i = start;i < end;i++){
		//hit = 0;
		num_scores = 0;
		if(!data->sb_file->buffer[i]->num_hits){
		//	unaligned_to_sam(ri[i]);
		}else{
		//	hit = 0;
			sum = prob2scaledprob(0.0);
			
			for(c = 0; c < data->sb_file->buffer[i]->num_hits;c++){
				
				//DPRINTF2("%s: %d ->%d", data->sb_file->buffer[i]->name,data->sb_file->buffer[i]->start[hit],data->sb_file->buffer[i]->stop[hit] );
				
				genomic_sequence = data->gc[i]->genomic_sequences[c];// data->sb_file->buffer[i]->genomic_sequences[c];
				len = data->gc[i]->g_len[c];//    data->sb_file->buffer[i]->g_len[c];
				hmm = glocal_forward_log_Y(hmm,data->sb_file->buffer[i]->sequence,  genomic_sequence,data->sb_file->buffer[i]->len,  len ,c);
				if(!hmm){
					
					fprintf(stdout,"%s (%d): hit:%d %lld ->%lld\n", data->sb_file->buffer[i]->name,i,c,data->sb_file->buffer[i]->start[c],data->sb_file->buffer[i]->stop[c] );
				
					for(len = 0;len <  data->sb_file->buffer[i]->len;len++){
						fprintf(stdout,"%d ", data->sb_file->buffer[i]->sequence[len]);
					}
					fprintf(stdout,"\n");
					for(len = 0;len < data->gc[i]->g_len[c];len++){
						fprintf(stdout,"%d ", data->gc[i]->genomic_sequences[c][len]);
					}
					fprintf(stdout,"\n");
					//return FAIL;
				}
				
				hmm = random_model_calc(hmm,data->sb_file->buffer[i]->sequence,  genomic_sequence,data->sb_file->buffer[i]->len,  len, prob2scaledprob(1.0),0);
				
				scores[c] = hmm->score;// + ri[i]->priors[hit];
				genome_scores[c] = hmm->random_score ;// ri[i]->priors[c];
			
			}
			hmm = random_model_calc(hmm,data->sb_file->buffer[i]->sequence,  genomic_sequence,data->sb_file->buffer[i]->len,  len , prob2scaledprob(1.0),1);
			unaligned = hmm->random_score;// + prob2scaledprob(data->param->unaligned_prior);
			
			for(c = 0; c < data->sb_file->buffer[i]->num_hits;c++){
				unaligned += genome_scores[c];// + prob2scaledprob(1.0 - scaledprob2prob(genome_scores[c]));
			}
			
			for(c = 0; c < data->sb_file->buffer[i]->num_hits;c++){
				sum = prob2scaledprob(1.0f);
				for(j = 0; j < data->sb_file->buffer[i]->num_hits;j++){
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
			for(c = 0; c < data->sb_file->buffer[i]->num_hits;c++){
				
				sum = logsum(sum,scores[c]);
			}
			unaligned = unaligned - sum;
			
			
			
			
			max = prob2scaledprob(0.0);
			
			k = 1;
			
			for(c = 0; c < data->sb_file->buffer[i]->num_hits;c++){
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
			/*if(i < 10){
			 
				fprintf(stderr,"%s	un	-1	%f	%f	%f\n", ri[i]->name, scaledprob2prob(unaligned),unaligned,unaligned) ;
				for(c = 0; c < ri[i]->nhits;c++){
			 fprintf(stderr,"%d	%d	%f	%f	%f\n",i, c,scaledprob2prob(scores[c] ),scores[c] ,scores[c]) ;
				}
				fprintf(stderr,"\n");
				fprintf(stderr,"Which is best??? : %f %f\n",max, unaligned);
			 }*/
			
			if(unaligned  > max){
			//	unaligned_to_sam(ri[i]);
			}else{
				
				//if
		//		hit = 0;
				//c =0;
				
				
				//fprintf(stderr,"\n%s	un	%f	%f	%f\n", data->sb_file->buffer[i]->name, scaledprob2prob(unaligned),unaligned,unaligned) ;
				for(j = 0; j < data->sb_file->buffer[i]->num_hits;j++){
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
						RUN(get_chr_start_stop(data->sb_file->si,g_int,data->sb_file->buffer[i]->start[j], data->sb_file->buffer[i]->stop[j]));
						g_int->start -=ALIGNMENT_FLANKING_LENGTH;
						g_int->stop += ALIGNMENT_FLANKING_LENGTH;
						
						
						//fprintf(stderr,"%s	%d	%f	%f	%d\n",data->sb_file->buffer[i]->name, j,scaledprob2prob(scores[j] ),scores[j] ,flag) ;
						genomic_sequence = data->gc[i]->genomic_sequences[j];//  data->sb_file->buffer[i]->genomic_sequences[j];
						len = data->gc[i]->g_len[j];//   data->sb_file->buffer[i]->g_len[j];
						//DPRINTF2("%s ( %d)\n",genomic_sequence,len);
						
						
						/*if(g_int->strand){
						  RUN(reverse_complement_sequence(data->sb_file->buffer[i]->sequence , data->sb_file->buffer[i]->len),"revcomp failed.");
						  RUN(reverse(data->sb_file->buffer[i]->basequal, data->sb_file->buffer[i]->len),"rev failed.");
						  }*/
						
						
						aln = glocal_viterbi_log_Y(hmm,data->sb_file->buffer[i]->sequence,  genomic_sequence,data->sb_file->buffer[i]->len,  len);
						
						if(!aln){
							fprintf(stderr,"seq: %d hitNO:%d\n",i,j);
							fprintf(stderr,"%s aligned against:\n%s %lld -%lld %d\n",data->sb_file->buffer[i]->name, g_int->chromosome,g_int->g_start,g_int->g_stop,g_int->strand);
						}
						
						
						if(g_int->strand){
							
							aln = reverse_path(aln);
							RUN(reverse_complement_sequence(data->sb_file->buffer[i]->sequence , data->sb_file->buffer[i]->len));//,"revcomp failed.");
							RUN(reverse(data->sb_file->buffer[i]->qual, data->sb_file->buffer[i]->len));//,"rev failed.");
							flag |= 16;
						}
						if(num_scores == 1){
							if( scaledprob2prob( scores[j]) >= data->sb_file->buffer[i]->qual){
								scores[j] = prob2scaledprob(data->sb_file->buffer[i]->qual);
							}
						}
						
						align_to_sam(g_int, data->sb_file->buffer[i], aln, flag,scaledprob2prob( scores[j]));
						
						if(g_int->strand){
							
							RUN(reverse_complement_sequence(data->sb_file->buffer[i]->sequence , data->sb_file->buffer[i]->len));//,"revcomp failed.");
							RUN(reverse(data->sb_file->buffer[i]->qual, data->sb_file->buffer[i]->len));//,"rev failed.");
						}
						
					}
					/*align_to_sam(
						     aln,
						     struct sam_bam_entry* sb_ptr
						     seed_start,
						     ri[i]->strand[j] | flag,
						     scaledprob2prob( scores[j]),
						     ri[i],
						     0
						     );*/
					

					//if(i < 10){
					//	fprintf(stderr,"%f %d (flag...)\n",scaledprob2prob(scores[j]),flag );
					//}
					/*
					
					if(flag != -1){
						RUN(get_chr_start_stop(data->sb_file,i,hit,g_int),"get_chr_start_stop failed");
						g_int->start -=ALIGNMENT_FLANKING_LENGTH;
						g_int->stop += ALIGNMENT_FLANKING_LENGTH;
						//DPRINTF2("%s:%d-%d \n",g_int->chromosome,g_int->start,g_int->stop);
						genomic_sequence = get_sequence(data->index,g_int);
						len = (int) strlen(genomic_sequence);
						RUN(ACGT_to_0123(genomic_sequence,&len),"ACGT_to_0123 failed");
						aln = glocal_viterbi_log_Y(hmm,data->sb_file->buffer[i]->sequence,  genomic_sequence,data->sb_file->buffer[i]->len,  len);
						aln = reverse_path(aln);
						
						
						
						if(ri[i]->strand[j] != 0){
							//reverse_complement2(ri[i]->seq,ri[i]->len);
							if(ri[i]->hits[j] > ALIGNMENT_FLANKING_LENGTH){
								strand = get_sequence(db,seq,ri[i]->hits[j]-ALIGNMENT_FLANKING_LENGTH, ri[i]->len+2*ALIGNMENT_FLANKING_LENGTH,1);
								seed_start = ri[i]->hits[j]-ALIGNMENT_FLANKING_LENGTH;
							}else{
								strand = get_sequence(db,seq,0, ri[i]->len+2*ALIGNMENT_FLANKING_LENGTH,1);
								seed_start =  0;
							}
							seq = transform_16to5_sequence(seq , ri[i]->len+2*ALIGNMENT_FLANKING_LENGTH);
							//reverse_complement2(ri[i]->seq,ri[i]->len);
							aln = glocal_viterbi_log_Y(hmm,ri[i]->seq, seq, ri[i]->len, ri[i]->len+2*ALIGNMENT_FLANKING_LENGTH);
							aln = reverse_path(aln);
							ri[i]->seq = reverse_complement(ri[i]->seq,ri[i]->len);
							
							
							if(num_scores == 1){
								if( scaledprob2prob( scores[j]) < ri[i]->mapq){
									ri[i]->mapq = scaledprob2prob( scores[j]) ;
								}
								align_to_sam(
									     aln,
									     db,
									     seed_start,
									     ri[i]->strand[j] | flag,
									     ri[i]->mapq,
									     ri[i],
									     0
									     );
							}else{
								
								align_to_sam(
									     aln,
									     db,
									     seed_start,
									     ri[i]->strand[j] | flag,
									     scaledprob2prob( scores[j]),
									     ri[i],
									     0
									     );
							}
							free(aln);
							ri[i]->seq = reverse_complement(ri[i]->seq,ri[i]->len);
						}else{
							if(ri[i]->hits[j] > ALIGNMENT_FLANKING_LENGTH){
								strand = get_sequence(db,seq,ri[i]->hits[j]-ALIGNMENT_FLANKING_LENGTH, ri[i]->len+2*ALIGNMENT_FLANKING_LENGTH,0);
								seed_start = ri[i]->hits[j]-ALIGNMENT_FLANKING_LENGTH;
							}else{
								strand = get_sequence(db,seq,0, ri[i]->len+2*ALIGNMENT_FLANKING_LENGTH,0);
								seed_start = 0;
							}
							
							
							//strand = get_sequence(db,seq,ri[i]->hits[j]-10, ri[i]->len+20,0);
							seq = transform_16to5_sequence(seq , ri[i]->len+2*ALIGNMENT_FLANKING_LENGTH);
							
							
							aln = glocal_viterbi_log_Y(hmm,ri[i]->seq, seq, ri[i]->len, ri[i]->len+2*ALIGNMENT_FLANKING_LENGTH );
							if(num_scores == 1){
								if( scaledprob2prob( scores[j]) < ri[i]->mapq){
									ri[i]->mapq = scaledprob2prob( scores[j]) ;
								}
								align_to_sam(
									     aln,
									     db,
									     seed_start,
									     ri[i]->strand[j] | flag,
									     ri[i]->mapq,
									     ri[i],
									     0
									     );
							}else{
								align_to_sam(
									     aln,
									     db,
									     seed_start,
									     ri[i]->strand[j] | flag,
									     scaledprob2prob(scores[j] ),
									     ri[i],
									     0
									     );
							}
							free(aln);
						}
					}*/
				}
				
				
				//for(j = 0; j < LIST_STORE_SIZE+1;j++){
				//	top_scores[j] = 0ul;
				//}
			}
		}
	}
	//free(random_scores);
	//free(genome_scores);
	//free(genome_prior);
	//free(scores);
	free_genome_interval(g_int);
	//free(top_scores);
	//free(seq);
	pthread_exit((void *) 0);
ERROR:
	pthread_exit((void *) 0);
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
					fprintf(stderr,"%c","ACGTN"[a[i]]);
				}
				fprintf(stderr,"\n");
				for(i = 0; i < m;i++){
					fprintf(stderr,"%c","ACGTN"[b[i]]);
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
	ASSERT(m> 10, " genome seq too short");
	
	int i,j;//,c;
	float** M = 0;
	//float** C = hmm->fC;
	float** X = 0;
	float** Y = 0;
	
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
	
	
	char* seqa = a -1;
	char* seqb = b -1;
	float* prob = hmm->prob;
	
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
		//Y[i][j] = logsum(Y[i][j], prob[(seqa[i] << 9) | (seqb[j-1] << 6) | (5 << 3) | seqb[j]] );
	}
	
	
	M[0][m-1] = prob2scaledprob(0.0f);
	X[0][m-1] = prob2scaledprob(0.0f);
	Y[0][m-1] = prob2scaledprob(0.0f);
	//Y[0][m-1] = prob[   (EGAP << 12) | (seqb[m-2] << 8) | (EGAP << 4) | seqb[m-1]] + Y[0][m-2];
	
	
	M[0][m] = prob2scaledprob(0.0f);
	X[0][m] = prob2scaledprob(0.0f);
	Y[0][m] = prob2scaledprob(0.0f);
	//Y[0][m] = prob[   (EGAP << 12) | (seqb[m-1] << 8) | (EGAP << 4) | seqb[m]] + Y[0][m-1];
	
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
	int i,j,c,g;
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
	
	c = 0;
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
ERROR:
	
	return FAIL;
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


void init_thread_hmms(struct thread_data* thread_data_array,struct hmm* localhmm,int numthreads)
{
	int i,c;
	for(c = 0; c < numthreads;c++){
		for(i = 0; i < 65536;i++){
			thread_data_array[c].hmm->prob[i] = localhmm->prob[i];
			thread_data_array[c].hmm->tprob[i] = prob2scaledprob(0.0f);
		}
		for(i = 0; i < 25;i++){
			thread_data_array[c].hmm->random_genome_model[i] = localhmm->random_genome_model[i];
			thread_data_array[c].hmm->random_genome_model_t[i] = prob2scaledprob(0.0f);
			thread_data_array[c].hmm->random_read_model[i] = localhmm->random_read_model[i];
			thread_data_array[c].hmm->random_read_model_t[i] = prob2scaledprob(0.0f);
			
		}
	}
	
}

void entangle_hmms(struct thread_data* thread_data_array,struct hmm* localhmm,int numthreads)
{
	int i,j,c;
	
	//fprintf(stderr,"0:%d\n",thread_data_array[0].hmm->train_count);
	for(i = 0; i < 65536;i++){
		localhmm->tprob[i] = thread_data_array[0].hmm->tprob[i];
	}
	for(j = 0; j < 25;j++){
		localhmm->random_genome_model_t[j] = logsum(localhmm->random_genome_model_t[j] ,thread_data_array[0].hmm->random_genome_model_t[j]);
		localhmm->random_read_model_t[j] = logsum(localhmm->random_read_model_t[j] ,thread_data_array[0].hmm->random_read_model_t[j]);
	}
	for(c = 1; c < numthreads;c++){
		
		
		for(i = 0; i < 65536;i++){
			localhmm->tprob[i] = logsum(localhmm->tprob[i], thread_data_array[c].hmm->tprob[i]);
		}
		/*for(i = 0; i < 31;i++){
		 for(j = 0; j < 25;j++){
		 localhmm->genome_model_t[i][j] = logsum(localhmm->genome_model_t[i][j] ,thread_data_array[c].hmm->genome_model_t[i][j]);
		 }
		 }*/
		for(j = 0; j < 25;j++){
			localhmm->random_genome_model_t[j] = logsum(localhmm->random_genome_model_t[j] ,thread_data_array[c].hmm->random_genome_model_t[j]);
			localhmm->random_read_model_t[j] = logsum(localhmm->random_read_model_t[j] ,thread_data_array[c].hmm->random_read_model_t[j]);
		}
	}
}



struct hmm* init_hmm(int x, int y)
{
	struct hmm* hmm = NULL;
	int i,j,c;
	int r1,r2,g1,g2,key;
	float* sum_prob = 0;
	float sum;
	x = x + 2;
	y = y + 2;
	
	MMALLOC(hmm, sizeof(struct hmm));
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
	
	MMALLOC(hmm->tfM, sizeof(float**)*LIST_STORE_SIZE);
	MMALLOC(hmm->tfX, sizeof(float**)*LIST_STORE_SIZE);
	MMALLOC(hmm->tfY, sizeof(float**)*LIST_STORE_SIZE);
	
	
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
	
	
	
	hmm->x = x;
	hmm->y = y;
	
	for(i = 0; i < LIST_STORE_SIZE;i++){
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

void align_to_sam(struct genome_interval* g_int,struct sam_bam_entry* entry, char* aln,unsigned flag,float score)
{
	char nuc[] = "ACGTN";
	char cigarline[128];
	char mdline[128];
	char pline[500];
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
			pseq[i] =  "ACGTN"[entry->sequence[i]];//   nuc[(int)ri->seq[i] & 3 ];
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
		sprintf(pline,"%s\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%s",
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
			"*"//11 QUAL query QUALity (ASCII-33=Phred base quality)
			);
	}else{
		sprintf(pline,"%s\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%s",
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
			entry->qual//11 QUAL query QUALity (ASCII-33=Phred base quality)
			);
	}
	
	i = (int)strlen(pline);
	pline[i] = '\t';
	i++;
	sprintf(pline+i,"NM:i:%d",  mismatches);
	
	i = (int)strlen(pline);
	pline[i] = '\t';
	i++;
	sprintf(pline+i,"MD:Z:%s",  mdline);
	
	fprintf(stdout, "%s\n",pline);
}



void free_hmm(struct hmm* hmm)
{
	int i,j;
	int x;//,y;
	
	
	//y = hmm->y;
	
	if(hmm){
		x = hmm->x;
		if(hmm->tfM){
			for(i = 0; i < LIST_STORE_SIZE;i++){
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
			for(i = 0; i < LIST_STORE_SIZE;i++){
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
			for(i = 0; i < LIST_STORE_SIZE;i++){
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

int reverse(char* p,int len)
{
	int c, i, j;
	
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
ERROR:
	return FAIL;
}


int add_genome_sequences(struct sam_bam_file* sb_file,struct genome_sequences** gc, char* genome)
{
	faidx_t*  index = NULL;

	
	struct genome_interval* g_int = NULL;

	int i,j,len,c;
	
	int len_a,len_b,len_c;
	RUNP(index = get_faidx(genome));
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
			len_b = len;
//			sb_file->buffer[i]->g_len[j] = len;
			gc[i]->g_len[j] = len;
			//DPRINTF2("%s	(len:%d)",sb_file->buffer[i]->genomic_sequences[j],sb_file->buffer[i]->g_len[j]);
			
			/*if(i == 128621){
				for(c =0; c < len;c++){
					fprintf(stdout,"%c ",sb_file->buffer[i]->genomic_sequences[j][c]);
				}
				fprintf(stdout,"\n");

			}*/

			RUN(ACGT_to_0123(gc[i]->genomic_sequences[j] ,&len));
			
//			RUN(ACGT_to_0123(sb_file->buffer[i]->genomic_sequences[j] ,&len),"ACGT_to_0123 failed");
			
			len_c = len;
			/*f/or(c =0; c < len;c++){
				fprintf(stdout,"%d ",sb_file->buffer[i]->genomic_sequences[j][c]);
				
				
				
			}
			fprintf(stdout,"\n");*/
			
			//if(i == 128621){
			//fprintf(stdout,"%d %d %d %d %d\n",i,j, len_a,len_b,len_c);
			//}
			
			ASSERT(len_a == len_b,"%d %d %d %d %d\n",i,j, len_a,len_b,len_c);
			ASSERT(len_b == len_c,"%d %d %d %d %d\n",i,j, len_a,len_b,len_c);
			ASSERT(len_a == len_c,"%d %d %d %d %d\n",i,j, len_a,len_b,len_c);

			
			for(c =0; c < len;c++){
				ASSERT(gc[i]->genomic_sequences[j][c] < 5, "%d %d pos:%d  nuc:%d  ", i,j, c, gc[i]->genomic_sequences[j][c] ); 
			}

			
		}
	}
	
	/*for(i = 0; i < sb_file->num_read;i++){
		for(j = 0; j < sb_file->buffer[i]->num_hits;j++){
			len =sb_file->buffer[i]->g_len[j];
			for(c =0; c < len;c++){
				ASSERT(sb_file->buffer[i]->genomic_sequences[j][c] > 4, "");
			}
		}
	}*/
	free_genome_interval(g_int);
	free_faidx(index);
	return OK;
ERROR:
	if(g_int){
		free_genome_interval(g_int);
	}
	if(index){
		free_faidx(index);
	}
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
		gc[i]->g_len = NULL;
		MMALLOC(gc[i]->g_len,sizeof(int) * num_maxhits);
		MMALLOC(gc[i]->genomic_sequences,sizeof(char*) * num_maxhits);
		for (j = 0; j < num_maxhits; j++) {
			gc[i]->g_len[j] = 0;
			gc[i]->genomic_sequences[j] = NULL;
			
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
ERROR:
	return FAIL;
		
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
				MFREE(gc[i]->g_len);
				MFREE(gc[i]);
			}
				      
		}
	}
}





