
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
#include "rtr.h"
#include "pwrite.h"
#include "delve.h"

#include "delve_struct.h"

#include "shared_data.h"
#include "genome_rtree.h"
#include "genome_sequences.h"
#include "genome_priors.h"
#include "thread_data.h"
#include "hmm.h"
#include "io.h"

#define OPT_NTHREAD 1
#define OPT_OUT 2
#define OPT_SITER 3
#define OPT_GITER 4 
#define OPT_PSEUDOCOUNTS 5
#define OPT_GSPEUDO 6


#define STAT_ALL 0
#define STAT_IGNORE_UNALIGNED 1 

/* help message */
int print_help(char **argv);

/* main function running the main steps in Sequence*/
int run_delve(struct parameters* param);

/* estimate random model */
int run_estimate_random_models(struct shared_data* bsd);

/* estimate sequence model */
int run_estimate_sequence_model(struct shared_data* bsd);

/* score alignments and writed output to file */
int run_score_alignments(struct shared_data* bsd);

/* estimate genome model */
int run_estimate_genome_model(struct shared_data* bsd);

void* do_baum_welch_thread_random_model(void *threadarg);
void* do_baum_welch_thread(void *threadarg);
void* do_score_alignments_thread_hmm(void *threadarg);


int calculate_scores(struct hmm* hmm,struct sam_bam_entry* sam_entry, struct genome_sequences* gc,float temperature);
int alignment_stats(struct hmm* hmm,int num_hits, int mode);

/* Function to reverse base qualities */
int reverse(uint8_t* p,int len);

/* Function to write header - should really use the htslib library for this but couldn't figure out how.. */
int write_sam_header(struct sam_bam_file* sb_file,FILE* out);

int main (int argc, char *argv[]) 
{		
	struct parameters* param = NULL;
	int i,c;
	int tmp_val;
	tlog.echo_build_config();
     	MMALLOC(param, sizeof(struct parameters));
	param->infiles = NULL;
	param->genome = NULL;
	param->aln_infile = NULL;
	param->hmm_file = NULL;
	param->outdir = NULL;
	param->num_infiles = 0;
	param->num_threads = 4;
	param->num_maxhits = 10;
	param->pseudocounts = 10;
	param->genome_pseudocounts = 1;
	param->nogp = 0;
	param->siter = 5;
	param->giter = 5;
	param->devel = 0;
	
	while (1){	
		static struct option long_options[] ={
			{"t",required_argument,0,OPT_NTHREAD},
			{"siter",required_argument,0,OPT_SITER},
			{"giter",required_argument,0,OPT_GITER},			
			{"spseudo",required_argument,0,OPT_PSEUDOCOUNTS},
			{"gpseudo",required_argument,0,OPT_GSPEUDO},
			{"nogp",0,0,'n'},
			{"outdir",required_argument,0,OPT_OUT},			
			{"help",0,0,'h'},
			{"devel",0,0,'d'},
			{0, 0, 0, 0}
		};
		
		int option_index = 0;
		c = getopt_long_only (argc, argv,"nhd",long_options, &option_index);
		
		if (c == -1){
			break;
		}
		
		switch(c) {
		case OPT_PSEUDOCOUNTS:
			tmp_val = atoi(optarg);
			if(tmp_val == 0){
				ERROR_MSG("option spseudo cannot be zero!");
			}
			if(tmp_val >= 256){
				ERROR_MSG("option spseudo cannot be greater than 256!");
			}
			param->pseudocounts = (uint8_t) tmp_val;
			break;
		case OPT_GSPEUDO:
			tmp_val = atoi(optarg);
			if(tmp_val == 0){
				ERROR_MSG("option gpseudo cannot be zero!");
			}
			if(tmp_val >= 256){
				ERROR_MSG("option gpseudo cannot be greater than 256!");
			}
			param->genome_pseudocounts =  (uint8_t) tmp_val;
			break;
		case OPT_SITER:
			tmp_val = atoi(optarg);
			if(tmp_val == 0){
				ERROR_MSG("option siter cannot be zero!");
			}
			if(tmp_val >= 256){
				ERROR_MSG("option siter cannot be greater than 256!");
			}
			param->siter =  (uint8_t) tmp_val;
			break;
		case OPT_GITER:
			tmp_val = atoi(optarg);
			if(tmp_val == 0){
				ERROR_MSG("option giter cannot be zero!");
			}
			if(tmp_val >= 256){
				ERROR_MSG("option giter cannot be greater than 256!");
			}
			param->giter = (uint8_t) tmp_val;
			break;
		case OPT_NTHREAD:
			tmp_val = atoi(optarg);
			if(tmp_val == 0){
				ERROR_MSG("option thread cannot be zero!");
			}
			if(tmp_val >= 256){
				ERROR_MSG("option thread cannot be greater than 256!");
			}
			param->num_threads = (uint8_t) tmp_val;
			break;
		case OPT_OUT:
			param->outdir = optarg;
			break;
		case 'n':
			param->nogp = 1;
			break;
		case 'd':
			param->devel = 1;
			break;
		case 'h':
			RUN(print_help(argv));
			MFREE(param);
			exit(EXIT_SUCCESS);
			break;
		default:
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

	LOG_MSG("Starting run");
	
	if(!param->num_infiles){
		RUN(print_help(argv));
		ERROR_MSG("delve requires at least one input\n");

	}

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
		MFREE(param->infiles);
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

int print_help(char **argv)
{
	//const char description[] = "Aligns HMMs to new sequences.";
	const char usage[] = " <genome> <SAM/BAM file>";
	fprintf(stdout,"\nUsage: %s [-options] %s\n\n",basename(argv[0]) ,usage);
	
	fprintf(stdout,"Options:\n\n");

	//example of int option.. 
	//printf(opt_name,BUFFER_LEN,"%c%c%s <n>",'-','-',"name");
	
	fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--siter","Seq model training iterations." ,"[5]"  );
      	fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--giter","Genome model training iterations." ,"[5]"  );
      	fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--spseudo","Pseudocounts in Seq model." ,"[?]"  );
      	fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--gpseudo","Pseudocount in genome modes" ,"[?]"  );
      	fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--nogp","skip genome training." ,""  );
      	fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--t","number of threads." ,"[4]"  );
      	fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--outdir","prefix of output directory." ,""  );
	return OK;
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
        LOG_MSG("Quantifying read depth at putatitve mapping locations.");
	RUNP(bsd->sb_file = open_SAMBAMfile(bsd->param->aln_infile,bsd->buffer_size,bsd->num_maxhits , 0,0));
	RUNP(bsd->rtree = build_rtree(bsd->sb_file));
	
	/* assign max_seq_len */
	RUN(get_max_seq_len_from_sb_buffer(bsd->sb_file,&bsd->max_seq_len));
	
	RUN(close_SAMBAMfile(bsd->sb_file));
	//bsd->free(bsd);
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
	LOG_MSG("done");	
	
	RUN(clear_genome_sequences(bsd->gc,bsd->buffer_size, bsd->num_maxhits));
	RUN(close_SAMBAMfile(bsd->sb_file));
	
	/* Generating alignments  */
	LOG_MSG("Generating alignments...");
	RUN(set_output_file(bsd,"delve.sam"));
	RUNP(bsd->sb_file = open_SAMBAMfile(bsd->param->aln_infile,bsd->buffer_size,bsd->num_maxhits,0,0));
	RUN(write_sam_header(bsd->sb_file, bsd->pw->out_ptr ));	
	while(1){
		RUN(read_SAMBAM_chunk(bsd->sb_file,1,0));
		LOG_MSG("read:%d",bsd->sb_file->num_read);	
		if(!bsd->sb_file->num_read){
			break;
		}
	        RUN(add_genome_sequences(bsd));
		RUN(convert_buffer_ACGT_to_0123(bsd->sb_file));
		RUN(run_score_alignments(bsd));	       
		RUN(clear_genome_sequences(bsd->gc,bsd->buffer_size, bsd->num_maxhits));
        }
	
	RUN(close_SAMBAMfile(bsd->sb_file));
	LOG_MSG("Done.");
	
	if(param->nogp == 0){
		/* "gently" add in estimation of genome priors using sumuklated annealing..  */
		LOG_MSG("Estimate genome model");
		/* need to run on all data - this way is wrong... ??? */
		RUN(run_estimate_genome_model(bsd));
		
		LOG_MSG("done");
		
		/* Generating alignments  */
		LOG_MSG("Generating alignments...");
		
		RUN(set_output_file(bsd,"delve_gp.sam"));
		
		RUNP(bsd->sb_file = open_SAMBAMfile(bsd->param->aln_infile,bsd->buffer_size,bsd->num_maxhits,0,0));
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
	}
	//RUN(free_delve_region_data_from_tree(bsd->rtree));
	bsd->free(bsd);
	return OK;
ERROR:
	if(bsd){
		bsd->free(bsd);
	}
	return FAIL;
}

int calculate_scores(struct hmm* hmm, struct sam_bam_entry* sam_entry,struct genome_sequences* gc, float temperature)
{
	int i,len;
	
	for(i = 0; i < sam_entry->num_hits;i++){
		//genomic_sequence = gc->genomic_sequences[i];// data->sb_file->buffer[i]->genomic_sequences[c];
		len = gc->g_len[i];//    data->sb_file->buffer[i]->g_len[c];
		hmm = glocal_forward_log_Y(hmm,sam_entry->sequence,  gc->genomic_sequences[i],sam_entry->len, len ,i);				
		//hmm = random_model_calc(hmm,buffer[i]->sequence,  genomic_sequence,buffer[i]->len,  len, prob2scaledprob(1.0),0);
		RUN(random_model_genome_score(hmm,gc->genomic_sequences[i],len,prob2scaledprob(1.0)));
		hmm->alignment_scores[i] = hmm->score + gc->prior[i] * temperature;
		hmm->unaligned_scores[i] = hmm->unaligned_genome_score;
		//fprintf(stdout,"%d %f %f\n", i, gc->prior[i], hmm->alignment_scores[i]); 	   
	}
	RUN(random_model_read_score(hmm, sam_entry->sequence,sam_entry->len,prob2scaledprob(1.0)));
	//unaligned = hmm->unaligned_read_score;
	return OK;
ERROR:
	return FAIL;
}

/* takes raw scores, combines and normalized them. 
   there are two modes: 
   1)normalization between genome hots 
   2) normalization with all possibilities including random (unaligned...) 
 */

int alignment_stats(struct hmm* hmm,int num_hits, int mode)
{
	int i,j;
	float sum = 0.0f;

	/* Step 1: calculate raw probability of read not aligned  */
        for(i = 0; i < num_hits;i++){
		hmm->unaligned_read_score += hmm->unaligned_scores[i];//   genome_scores[c];// + prob2scaledprob(1.0 - scaledprob2prob(genome_scores[c]));
	}

	/* Step 2: calcualte raw probability of read hitting location "i" */
	for(i = 0; i < num_hits;i++){
		sum = prob2scaledprob(1.0f);
		for(j = 0; j < num_hits;j++){
			if(i == j){
				//read is aligned at this seed position
				sum += hmm->alignment_scores[i];//  scores[j];// + genome_scores[j];
			}else{
				sum += hmm->unaligned_scores[i];//    genome_scores[j];// + prob2scaledprob(1.0 - scaledprob2prob(genome_scores[j])) ;
			}
		}
		hmm->tmp_scores[i] = sum;
		
	}


	/* Ok - got the raw probabilities  */

	sum = hmm->unaligned_read_score;
	if(mode == STAT_IGNORE_UNALIGNED){ // i.e. contribution for unaligned possibility is ignored. 
		sum = prob2scaledprob(0.0f);
	}
	/* get sum */
	for(i = 0; i < num_hits;i++){
		sum = logsum(sum, hmm->tmp_scores[i]);
			     
	}
	/* norm */
	hmm->unaligned_read_score = hmm->unaligned_read_score - sum;
	for(i = 0; i < num_hits;i++){
		hmm->alignment_scores[i] = hmm->tmp_scores[i] - sum;
	}
	
	return OK;
}


int run_estimate_random_models(struct shared_data* bsd)
{
	struct thread_data** td = NULL;
	int i;
	int num_threads;
	int status;
	num_threads = bsd->param->num_threads;

	/* I think I should add pseudocounts  */

	RUN(add_pseudo_count(bsd->master_hmm,(float) bsd->param->pseudocounts));
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

	iterations = bsd->param->siter;
	num_threads = bsd->param->num_threads;
	/* initialize datastructs to pass bsd (shared) and thread_id's (private)  */
	RUNP(td = init_thread_data(bsd,num_threads));
	/* I think I should add pseudocounts  */
	for(iter = 0 ; iter< iterations;iter++){
		LOG_MSG("Iteration %d.",iter);
		RUN(add_pseudo_count(bsd->master_hmm,(float) bsd->param->pseudocounts));
	 	/* copy parameters from master hmm into copies used by threads...  */
		RUN(init_thread_hmms(bsd));

		/* set temperature */
		bsd->temperature = log10f((float) iter / (float) iterations  * 9.0f +1.0f);
		
		for(i = 0; i < num_threads;i++){
			bsd->thread_forward[i] = 0.0; // (i.e P = 1.0);
			//LOG_MSG("%d score %f.",i,bsd->thread_forward[i]);

		}
		/* kick off jobs  */
		for(i = 0; i < num_threads;i++){
			if((status = thr_pool_queue(bsd->pool,do_baum_welch_thread,td[i])) == -1) fprintf(stderr,"Adding job to queue failed.");	
		}
		/* wait for all jobs to finish */
		thr_pool_wait(bsd->pool);
		
		//for(i = 0; i < num_threads;i++){
			//	fprintf(stderr,"%d score %f\n",i,bsd->thread_forward[i]);
//		        LOG_MSG("%d score: %f\n", bsd->sb_file->buffer[i]->fscore);
		//}
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

int run_estimate_genome_model(struct shared_data* bsd)
{
	struct thread_data** td = NULL;
	int i,iter;
	int num_threads;
	int status;
	int iterations;
	
	iterations = bsd->param->giter;
	num_threads = bsd->param->num_threads;

	/* init_genome priors */

	RUN(add_genome_pseudocounts(bsd));
	RUN(re_estimate_genome_priors(bsd));
        /* initialize datastructs to pass bsd (shared) and thread_id's (private)  */
	RUNP(td = init_thread_data(bsd,num_threads));
	/* I think I should add pseudocounts  */
	for(iter = 0 ; iter <= iterations;iter++){
		LOG_MSG("Iteration %d.",iter);
		RUN(add_pseudo_count(bsd->master_hmm,(float) bsd->param->pseudocounts));
		RUN(add_genome_pseudocounts(bsd));
	 	/* copy parameters from master hmm into copies used by threads...  */
		RUN(init_thread_hmms(bsd));
		//RUN(copy_genome_priors_to_gc(bsd));

		/* set temperature */
		bsd->temperature = log10f((float) iter / (float) iterations  * 9.0f +1.0f);
		
		RUNP(bsd->sb_file = open_SAMBAMfile(bsd->param->aln_infile,bsd->buffer_size,bsd->num_maxhits,0,0));
		while(1){
			RUN(read_SAMBAM_chunk(bsd->sb_file,1,0));
			if(!bsd->sb_file->num_read){
				break;
			}
			LOG_MSG("GP est. iter %d read:%d",iter, bsd->sb_file->num_read);
			RUN(add_genome_sequences(bsd));
			RUN(convert_buffer_ACGT_to_0123(bsd->sb_file));
			RUN(copy_genome_priors_to_gc(bsd));

			/* kick off jobs  */
			for(i = 0; i < num_threads;i++){
				if((status = thr_pool_queue(bsd->pool,do_baum_welch_thread,td[i])) == -1) fprintf(stderr,"Adding job to queue failed.");	
			}

			/* wait for all jobs to finish */
			thr_pool_wait(bsd->pool);

			/* entangle HMMs - copy estimated counts from thread hmm copies back into the master HMM. */
			RUN(entangle_hmms(bsd));
			RUN(entangle_genome_priors(bsd));
			
			RUN(clear_genome_sequences(bsd->gc,bsd->buffer_size, bsd->num_maxhits));
		}
		RUN(close_SAMBAMfile(bsd->sb_file));
       
		
                /* re-estimate parameters.. */
		RUN(re_estimate(bsd->master_hmm));
		
		/* re-estimate random model more... */
		RUN(re_estimate_random(bsd->master_hmm));

		LOG_MSG("re-estimate genome priors...");
		/* re-estimate genome priors.. */
		RUN(re_estimate_genome_priors(bsd));
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
			RUN(random_model_genome_score(hmm,genomic_sequence,len,prob2scaledprob(1.0)));
			
		
		}
		RUN(random_model_read_score(hmm,buffer[i]->sequence,buffer[i]->len,prob2scaledprob(1.0)));
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
	float temperature = 1.0f;
	int thread_id = -1;
	int num_threads = 0;
	int num_sequences = 0;
	int i;
	int start;
	int stop;
	
	data = (struct thread_data *) threadarg;

	thread_id = data->thread_id;

	num_sequences = data->bsd->sb_file->num_read;
	num_threads = data->bsd->num_threads;
	
	hmm = data->bsd->thread_hmm[thread_id];
	buffer = data->bsd->sb_file->buffer;
	gc = data->bsd->gc;

	temperature = data->bsd->temperature;
	
	start = 0;
	stop = 0;
	RUN(get_start_stop_for_threads(num_threads, num_sequences,thread_id,&start,&stop));

	int len;
	//int strand = 0;
	char* genomic_sequence = NULL;
	int c;
	
	//LOG_MSG("thread %d : %d -%d.",thread_id,start,stop);
	
	for(i = start; i < stop;i++){
		if(buffer[i]->num_hits){
			//if(ri[i]->identity[0] >= 0.0f   ){
			//hit = 0;0
			//sum = prob2scaledprob(0.0);
			
			RUN(calculate_scores(hmm, buffer[i],gc[i],temperature));
			RUN(alignment_stats(hmm,buffer[i]->num_hits,STAT_IGNORE_UNALIGNED));
			/* add max prob to thread score for bookkeeping. */
			//fprintf(stdout,"score: best: %f\n",scores[max - sum);
			//fprintf(stdout,"thread %d: %f adding %f\n",thread_id,  data->bsd->thread_forward[thread_id],scores[max] - sum);
			for(c = 0; c < buffer[i]->num_hits;c++){
				data->bsd->thread_forward[thread_id] += hmm->alignment_scores[c];//   scores[c] - sum;
			}
			//hit = 0;
			//if(i < 5){
			
			///}
			//max =  prob2scaledprob(scaledprob2prob(max - sum ) -scaledprob2prob(max2 - sum ) );	
			for(c = 0; c < buffer[i]->num_hits;c++){
				//DPRINTF2("%s:%d-%d \n",g_int->chromosome,g_int->start,g_int->stop);
				genomic_sequence = gc[i]->genomic_sequences[c];//  data->sb_file->buffer[i]->genomic_sequences[c];
				len = gc[i]->g_len[c];//   data->sb_file->buffer[i]->g_len[c];
				
				DPRINTF2("len: %d and %d",buffer[i]->len,  len);
				hmm = glocal_backward_log_Y(hmm,buffer[i]->sequence,  genomic_sequence,buffer[i]->len,  len);
				/* weigth sequence by: 1) mapping quality (scores -sum) and 2) weigth of sequences at the loci.. (priors)   */
				hmm = get_prob_log_Y(hmm,buffer[i]->sequence,  genomic_sequence,buffer[i]->len,len,hmm->alignment_scores[c] + gc[i]->alignment_weigth[c]  ,c);// max ,c);//  hmm->score + hmm->score - sum + ri[i]->priors[hit]  );
				gc[i]->prior_e[c] = logsum(gc[i]->prior_e[c], hmm->alignment_scores[c]);
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
	char* genomic_sequence = NULL;
	char* aln = 0;
	float max;
	int thread_id = -1;
	int num_threads = 0;
	int num_sequences = 0;
	int i,j,c,k,f;
	int start,stop;
	int num_scores = 0;
	int flag;
       	int len;
		
	data = (struct thread_data *) threadarg;

	thread_id = data->thread_id;

	num_sequences = data->bsd->sb_file->num_read;
	num_threads = data->bsd->num_threads;
	
	hmm = data->bsd->thread_hmm[thread_id];
	pw = data->bsd->pw;
	buffer = data->bsd->sb_file->buffer;
	g_int = data->bsd->g_int_working[thread_id];
	gc = data->bsd->gc;

	start = 0;
	stop = 0;
	
	RUN(get_start_stop_for_threads(num_threads, num_sequences,thread_id,&start,&stop));

//	LOG_MSG("Looking at %d-%d.",start,stop);
	for(i = start;i < stop;i++){
	
		//hit = 0;
		num_scores = 0;
		if(!buffer[i]->num_hits){
		//	unaligned_to_sam(ri[i]);
		}else{

			RUN(calculate_scores(hmm, buffer[i],gc[i], 1.0)); /* temp is 1 - i.e. no effect... */
			RUN(alignment_stats(hmm,buffer[i]->num_hits,STAT_ALL));
		        
			max = prob2scaledprob(0.0);
			
			k = 1;
			
			for(c = 0; c < buffer[i]->num_hits;c++){
				if(hmm->alignment_scores[c] > max){
					max = hmm->alignment_scores[c];
					k = 1;
				}else if(hmm->alignment_scores[c] == max){
					k++;
				}
			}
			f = 0;
			if(k > 1){
				k = 0;
				f = 1;
			}else{
				k = 0;
			}
			
			if(hmm->unaligned_read_score  > max){
				//fprintf(stdout,"Unaligned... \n");
				//	unaligned_to_sam(ri[i]);
			}else{
				for(j = 0; j < buffer[i]->num_hits;j++){
					flag = 0;
					if(hmm->alignment_scores[j]  == max){
						if(k){
							if(k != f ){
								flag |= 0x100;
							}
							f++;
						}
						/* NEEDS REVIEW - NO IDEA WHAT THIS IS  */
					}else if(scaledprob2prob(hmm->alignment_scores[j]) > (1.0 / (double)buffer[i]->num_hits)){
						flag |= 0x100;
					}else{
						flag = -1;
					}
					if(flag != -1){
						RUN(get_chr_start_stop(data->bsd->sb_file->si,g_int,buffer[i]->start[j],buffer[i]->stop[j]));
						g_int->start -=ALIGNMENT_FLANKING_LENGTH;
						g_int->stop += ALIGNMENT_FLANKING_LENGTH;
					        
						genomic_sequence = gc[i]->genomic_sequences[j];//  data->sb_file->buffer[i]->genomic_sequences[j];
						len = gc[i]->g_len[j];//   data->sb_file->buffer[i]->g_len[j];
						
						aln = glocal_viterbi_log_Y(hmm,buffer[i]->sequence,  genomic_sequence,buffer[i]->len,  len);
												
						if(g_int->strand){	
							aln = reverse_path(aln);
							RUN(reverse_complement_sequence(buffer[i]->sequence , buffer[i]->len));//,"revcomp failed.");
							RUN(reverse(buffer[i]->base_qual, buffer[i]->len));//,"rev failed.");
							flag |= 16;
						}
						if(num_scores == 1){
							if( scaledprob2prob(hmm->alignment_scores[j]) >= buffer[i]->qual){
								hmm->alignment_scores[j] = prob2scaledprob(buffer[i]->qual);
							}
						}
						RUN(align_to_sam(pw, g_int, buffer[i], thread_id, aln, flag,scaledprob2prob(hmm->alignment_scores[j])));
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
