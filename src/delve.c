

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <string.h>
#include <inttypes.h>
#include <ctype.h>
#include <float.h>

#include "tldevel.h"
#include "rbtree.h"
#include "htsglue.h"
#include "rtr.h"
#include "pwrite.h"
#include "delve.h"

#include "delve_struct.h"
#include "dyn_prog.h"
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
#define OPT_ID 7

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
int generate_alignments(struct shared_data* bsd, char* suffix, int all);
int run_score_alignments(struct shared_data* bsd);

int determine_best_hit_in_case_of_equal_mapq(struct sam_bam_entry* entry, struct hmm* hmm, float max);


/* estimate genome model */
int run_estimate_genome_model(struct shared_data* bsd);

void* do_baum_welch_thread_random_model(void *threadarg);
void* do_baum_welch_thread(void *threadarg);
void* do_score_alignments_thread_hmm(void *threadarg);


int calculate_scores(struct hmm* hmm,struct sam_bam_entry* sam_entry, struct genome_sequences* gc,int mode, float temperature);
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
	param->id = NULL;
	param->num_infiles = 0;
	param->num_threads = 4;
	param->num_maxhits = 20;
	param->pseudocounts = 10;
	param->genome_pseudocounts = 10;
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
			{"outdir",required_argument,0,OPT_OUT},			
			{"id",required_argument,0,OPT_ID},					       
			{"nogp",0,0,'n'},
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
		case OPT_ID:
			param->id = optarg;
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
	RUNP(bsd->sb_file = open_SAMBAMfile(bsd->param->aln_infile,bsd->buffer_size,bsd->num_maxhits , 0,100));	
	RUNP(bsd->rtree = build_rtree(bsd->sb_file));
	
	/* assign max_seq_len */
	RUN(get_max_seq_len_from_sb_buffer(bsd->sb_file,&bsd->max_seq_len));
	
	RUN(close_SAMBAMfile(bsd->sb_file));
	//bsd->free(bsd);
	/* init all HMMs; now that I have the max_seq_len... */
        /* allocate HMMs... */
	RUN(init_shared_data_hmms(bsd));
	LOG_MSG("Iteration %d: AAAA: %f   real:%f.",-1 ,bsd->master_hmm->tprob[0] , scaledprob2prob(bsd->master_hmm->prob[0]));
	/* 2. build model 
	   Here I read in the first X (MAXNUMQUERY sequences and re-use them for model building 
	   etc. I.e. I open, read and only after all models are trained I close....
	*/
	/* open file again - read all alignments (param "1","0") */
	LOG_MSG("Opening alignment file...");
	RUNP(bsd->sb_file = open_SAMBAMfile(bsd->param->aln_infile,bsd->buffer_size,bsd->num_maxhits,0,100));	
	RUN(read_SAMBAM_chunk(bsd->sb_file,1,0));
	LOG_MSG("read: %d",bsd->sb_file->num_read);
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
	//RUN(run_estimate_random_models(bsd));
        
	LOG_MSG("done");
	/* Estmate Sequence Model  */
	LOG_MSG("Estimating Sequence model...");
	bsd->statmode = STAT_IGNORE_UNALIGNED;//  IGNORE_UNALIGNED;
	RUN(run_estimate_sequence_model(bsd));
	LOG_MSG("done");	

	/* Estmate Full Model  */
	LOG_MSG("Estimating Full model ...");
	bsd->statmode = STAT_ALL; 
	RUN(reset_prev_aln_priors(bsd));
	RUN(run_estimate_sequence_model(bsd));
	LOG_MSG("done");
	
	RUN(clear_genome_sequences(bsd->gc,bsd->buffer_size, bsd->num_maxhits));
	RUN(close_SAMBAMfile(bsd->sb_file));
	
	/* Generating alignments  */
	RUN(generate_alignments(bsd,"delve.sam",1));
	
	if(param->nogp == 0){
		/* "gently" add in estimation of genome priors using sumuklated annealing..  */
		LOG_MSG("Estimate genome model");
		/* need to run on all data - this way is wrong... ??? */
		RUN(run_estimate_genome_model(bsd));
		
		LOG_MSG("done");
		
	        /* Generating alignments  */
		RUN(generate_alignments(bsd,"delve_gp.sam",1));
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


/* This function calculates for hits 'i': */
/* - P(g_i, r | M) - alignment  */
/* - P(g_i | G)    - unaligned genome  */
/* - R(r   | R)    - unaligned read   */
/* That's it!  */

int calculate_scores(struct hmm* hmm, struct sam_bam_entry* sam_entry,struct genome_sequences* gc, int mode, float temperature)
{
	int i,j,len;
	float sum;
	/* STEP 1: calculate the probability of sequences given the various models  */
	for(i = 0; i < sam_entry->num_hits;i++){
		len = gc->g_len[i];

		/* - P(g_i, r | M) - alignment  */
		RUN(glocal_forward_log_Y(hmm,sam_entry->sequence,  gc->genomic_sequences[i],sam_entry->len, len ,i));
		hmm->alignment_scores[i] = hmm->score;//+ (gc->prior[i]-sum  + prob2scaledprob(0.7));// * temperature;

		/* - P(g_i | G)    - unaligned genome  */
		RUN(random_model_genome_score(hmm,gc->genomic_sequences[i],len,prob2scaledprob(0.0)));
		gc->genome_seq_score[i] = hmm->unaligned_genome_score;
		
	}
	

	/* - R(r   | R)    - unaligned read   */
	RUN(random_model_read_score(hmm, sam_entry->sequence,sam_entry->len,prob2scaledprob(0.0)));
	gc->read_score = hmm->unaligned_read_score;
	
	/* Step 2: mix priors...  */

	/*  prev_ priors do not work    */
	/* new idea - norm priors to sum to one  */
	sum = prob2scaledprob(0.0f); 
	for(i = 0; i < sam_entry->num_hits;i++){
		//gc->prior[i] = gc->prior[i] + gc->prev_aln_prior[i] ;// + prob2scaledprob(0.7f);
		sum = logsum(sum,gc->prior[i]);
	}
	
	for(i = 0; i < sam_entry->num_hits;i++){
		gc->prior[i] = gc->prior[i] -sum  + prob2scaledprob(0.5f);
	}
	
	/* what is the prior probability of aligning at i  */
	/* gc->prev_aln_prior * gc->genome_prior * p_aln   */
	/* what is the prior probability of not aligning at i  */
	/* 1 - (gc->prev_aln_prior * gc->genome_prior * p_aln)   */
	/*sum = prob2scaledprob(0.0f); 
	for(i = 0; i < sam_entry->num_hits;i++){
		gc->prior[i] = gc->prior[i] + gc->prev_aln_prior[i] ;// + prob2scaledprob(0.7f);
		sum = logsum(sum,gc->prior[i]);
	}
	u_p = gc->prev_unaln_prior +  prob2scaledprob(0.3f);
	sum = logsum(sum,u_p);
	*/
	//u_p = u_p - sum;TL is supported by a Fellowship from the Feilman Foundation.  and funds raised by the Sunsuper MACA Ride to Conquer Cancer. 
	

	/* Step 3: calculate the overall probability of alignment possibilities */
	
	// unaligned...

	hmm->unaligned_read_score = gc->read_score + prob2scaledprob(0.5f); //   + u_p;
	
	for(i = 0; i < sam_entry->num_hits;i++){
		hmm->unaligned_read_score += gc->genome_seq_score[i];// + prob2scaledprob(1.0f - scaledprob2prob(gc->prior[i])); 
		sum = prob2scaledprob(1.0f);// + prob2scaledprob(1.0f - scaledprob2prob(u_p)); 
		/* add genome prior and 0.7 prior of hitting the genome at all...  */
		for(j = 0; j < sam_entry->num_hits;j++){
			if(i == j){
				//read is aligned at this seed position
				sum += hmm->alignment_scores[j] + gc->prior[j];// scores[j];// + genome_scores[j];
			}else{
				sum += gc->genome_seq_score[j];// + prob2scaledprob(1.0f - scaledprob2prob(gc->prior[j]));
			}
		}
		hmm->tmp_scores[i] = sum;// + prob2scaledprob(0.7f);		
	}
	/* Ok - got the raw probabilities  */
	/* get sum */
	sum = prob2scaledprob(0.0f);
	for(i = 0; i < sam_entry->num_hits;i++){
		sum = logsum(sum, hmm->tmp_scores[i]);
	}
        
	if(mode == STAT_IGNORE_UNALIGNED){ // i.e. contribution for unaligned possibility is ignored. 
		for(i = 0; i < sam_entry->num_hits;i++){
			hmm->alignment_scores[i] = hmm->tmp_scores[i] - sum;
		}
		hmm->unaligned_read_score = prob2scaledprob(0.0);
		
	}else{
		sum = logsum(sum, hmm->unaligned_read_score);
		for(i = 0; i < sam_entry->num_hits;i++){
			hmm->alignment_scores[i] = hmm->tmp_scores[i] - sum;
		}
		hmm->unaligned_read_score = hmm->unaligned_read_score - sum;
	}
	/* record aln_priors for next round...  */
	for(i = 0; i < sam_entry->num_hits;i++){
		gc->prev_aln_prior[i] = hmm->alignment_scores[i];
	}
	gc->prev_unaln_prior = hmm->unaligned_read_score;
	
	if(!strcmp(sam_entry->name,"ORG8247_chr1_148556015_148556043_+_0")){
		for(i = 0; i < sam_entry->num_hits;i++){
			fprintf(stdout,"HIT:%d : %f\n",i,( hmm->alignment_scores[i]));
		}
		fprintf(stdout, "unaligned : %f\n", ( hmm->unaligned_read_score));
	}
	
	return OK;
ERROR:
	return FAIL;
}

/* takes raw scores, combines and normalized them. 
   there are two modes: 
   1) normalization between genome hits 
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
	//hmm->unaligned_read_score = hmm->unaligned_read_score + prob2scaledprob(0.3);
	
	
	/* Step 2: calcualte raw probability of read hitting location "i" */
	for(i = 0; i < num_hits;i++){
		sum = prob2scaledprob(1.0f);
		/* add genome prior and 0.7 prior of hitting the genome at all...  */
		for(j = 0; j < num_hits;j++){
			if(i == j){
				//read is aligned at this seed position
				sum += hmm->alignment_scores[j];// scores[j];// + genome_scores[j];
			}else{
				sum += hmm->unaligned_scores[j];//    genome_scores[j];// + prob2scaledprob(1.0 - scaledprob2prob(genome_scores[j])) ;
			}
		}
		hmm->tmp_scores[i] = sum;// + prob2scaledprob(0.7f);		
	}
	
	hmm->unaligned_read_score  = hmm->unaligned_read_score;// + prob2scaledprob(0.3);
	
	/* Ok - got the raw probabilities  */
	/* get sum */
	sum = prob2scaledprob(0.0f);
	for(i = 0; i < num_hits;i++){
		sum = logsum(sum, hmm->tmp_scores[i]);		     
	}
        
	if(mode == STAT_IGNORE_UNALIGNED){ // i.e. contribution for unaligned possibility is ignored. 
		for(i = 0; i < num_hits;i++){
			hmm->alignment_scores[i] = hmm->tmp_scores[i] - sum;
		}
		hmm->unaligned_read_score = prob2scaledprob(0.0);
	}else{
		sum = logsum(sum, hmm->unaligned_read_score);
		for(i = 0; i < num_hits;i++){
			hmm->alignment_scores[i] = hmm->tmp_scores[i] - sum;
		}
		hmm->unaligned_read_score = hmm->unaligned_read_score - sum;
	}	
	return OK;
}

int generate_alignments(struct shared_data* bsd, char* suffix, int all)
{
	char buffer[BUFFER_LEN];
       
	ASSERT(bsd != NULL, "No shared data!");

	if(suffix){
		snprintf(buffer,BUFFER_LEN,"%s",suffix);
	}else{
		snprintf(buffer,BUFFER_LEN,"%s","delve.sam");
	}
	
	LOG_MSG("Generating alignments...");		
	RUN(set_output_file(bsd,buffer));
		
	if(all){
		RUNP(bsd->sb_file = open_SAMBAMfile(bsd->param->aln_infile,bsd->buffer_size,bsd->num_maxhits,0,100));
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
	}else{
		
		/* this assumes a file is already open....  */
		RUN(write_sam_header(bsd->sb_file, bsd->pw->out_ptr ));
	
		RUN(run_score_alignments(bsd));
			
	}
	LOG_MSG("Done.");
	
	return OK;	
ERROR:
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

	RUN(add_pseudo_count_random(bsd->master_hmm,(float) bsd->param->pseudocounts));
	RUN(re_estimate_random(bsd->master_hmm)); /*  */

	/* copy parameters from master hmm into copies used by threads...  */
	RUN(init_thread_hmms(bsd));
	
	/* initialize datastructs to pass bsd (shared) and thread_id's (private)  */
	RUNP(td = init_thread_data(bsd,num_threads));
	
	/* kick off jobs */
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
	for(iter = 0 ; iter < iterations;iter++){
		LOG_MSG("Iteration %d: AAAA: %f   real:%f.",iter,bsd->master_hmm->tprob[0] , scaledprob2prob(bsd->master_hmm->prob[0]));
		RUN(add_pseudo_count(bsd->master_hmm,(float) bsd->param->pseudocounts));
	 	/* copy parameters from master hmm into copies used by threads...  */
		RUN(init_thread_hmms(bsd));

		/* set temperature */
		bsd->temperature = log10f((float) iter / (float) iterations  * 9.0f +1.0f);
		
		for(i = 0; i < num_threads;i++){
			bsd->thread_forward[i] = 0.0; // (i.e P = 1.0);
		}
		if(bsd->statmode == STAT_IGNORE_UNALIGNED){
			RUN(reset_prev_aln_priors(bsd));
		}
		
		/* kick off jobs  */
		for(i = 0; i < num_threads;i++){
			if((status = thr_pool_queue(bsd->pool,do_baum_welch_thread,td[i])) == -1) fprintf(stderr,"Adding job to queue failed.");	
		}
		/* wait for all jobs to finish */
		thr_pool_wait(bsd->pool);
		
		/* entangle HMMs - copy estimated counts from thread hmm copies back into the master HMM. */
		RUN(entangle_hmms(bsd));
		/* re-estimate parameters.. */
		RUN(re_estimate(bsd->master_hmm));
		/* re-estimate parameters.. */
		//RUN(re_estimate_random(bsd->master_hmm));
	
		if(bsd->param->devel){
			char buffer[BUFFER_LEN];
			if(bsd->statmode == STAT_ALL){
				snprintf(buffer,BUFFER_LEN,"aiter%d.delve.sam",iter+1);
			}else{
				snprintf(buffer,BUFFER_LEN,"siter%d.delve.sam",iter+1);
			}
			RUN(generate_alignments(bsd, buffer,0));
		}
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
	for(iter = 0; iter <= iterations;iter++){
		LOG_MSG("Iteration %d.",iter);
		RUN(add_pseudo_count(bsd->master_hmm,(float) bsd->param->pseudocounts));
		RUN(add_genome_pseudocounts(bsd));
	 	/* copy parameters from master hmm into copies used by threads...  */
		RUN(init_thread_hmms(bsd));
		//RUN(copy_genome_priors_to_gc(bsd));

		/* set temperature */
		bsd->temperature = log10f((float) iter / (float) iterations  * 9.0f +1.0f);
		
		RUNP(bsd->sb_file = open_SAMBAMfile(bsd->param->aln_infile,bsd->buffer_size,bsd->num_maxhits,0,100));
	
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
		if(bsd->param->devel){
			char buffer[BUFFER_LEN];
			snprintf(buffer,BUFFER_LEN,"giter%d.delve.sam",iter+1);
			RUN(generate_alignments(bsd, buffer,1));
		}
       
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

         /* copy parameters from master hmm into copies used by threads...  */
	RUN(init_thread_hmms(bsd));
	
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
	float max = 0;
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
	char* genomic_sequence = NULL;
	int len;
	for(i = start; i < stop;i++){
		max = prob2scaledprob(0.0);
       		for(hit = 0; hit < buffer[i]->num_hits;hit++){
		        genomic_sequence = gc[i]->genomic_sequences[hit];//  data->sb_file->buffer[i]->genomic_sequences[hit];
			len = gc[i]->g_len[hit];//  data->sb_file->buffer[i]->g_len[hit];
			if(gc[i]->alignment_weigth[hit] > max){
				max = gc[i]->alignment_weigth[hit];
			}
			RUN(random_model_genome_score(hmm,genomic_sequence,len, gc[i]->alignment_weigth[hit]));//  prob2scaledprob(1.0)));
			gc[i]->genome_seq_score[hit] = hmm->unaligned_genome_score;
		}
		RUN(random_model_read_score(hmm,buffer[i]->sequence,buffer[i]->len,max));// prob2scaledprob(1.0)));
		gc[i]->read_score = hmm->unaligned_read_score;
	}
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
	char* genomic_sequence = NULL;
	float temperature = 1.0f;
	int thread_id = -1;
	int num_threads = 0;
	int num_sequences = 0;
	int i;
	int c;
	int len;
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
	//LOG_MSG("thread %d : %d -%d.",thread_id,start,stop);
	
	for(i = start; i < stop;i++){
		if(buffer[i]->num_hits){
			//if(ri[i]->identity[0] >= 0.0f   ){
			//hit = 0;0
			//sum = prob2scaledprob(0.0);
			
			RUN(calculate_scores(hmm, buffer[i],gc[i],data->bsd->statmode,temperature));
			//RUN(alignment_stats(hmm,buffer[i]->num_hits,data->bsd->statmode));//   STAT_IGNORE_UNALIGNED));

		        
			//hit = 0;
			//if(i < 5){
			//	fprintf(stderr,"%d	%f	%f	%f	%f\n",c,sum,max, scaledprob2prob(max - sum ) ,scaledprob2prob(max2 - sum )  );
			//}
			//max =  prob2scaledprob(scaledprob2prob(max - sum ) -scaledprob2prob(max2 - sum ) );
			
			/* add max prob to thread score for bookkeeping. */
			//fprintf(stdout,"score: best: %f\n",scores[max - sum);
			//fprintf(stdout,"thread %d: %f adding %f\n",thread_id,  data->bsd->thread_forward[thread_id],scores[max] - sum);
		        for(c = 0; c < buffer[i]->num_hits;c++){
		        	data->bsd->thread_forward[thread_id] += hmm->alignment_scores[c];//   scores[c] - sum;
			}
			//max =  prob2scaledprob(scaledprob2prob(max - sum ) -scaledprob2prob(max2 - sum ) );	
			for(c = 0; c < buffer[i]->num_hits;c++){
				//DPRINTF2("%s:%d-%d \n",g_int->chromosome,g_int->start,g_int->stop);
				genomic_sequence = gc[i]->genomic_sequences[c];//  data->sb_file->buffer[i]->genomic_sequences[c];
				len = gc[i]->g_len[c];//   data->sb_file->buffer[i]->g_len[c];
				
				DPRINTF2("len: %d and %d",buffer[i]->len,  len);
				RUN(glocal_backward_log_Y(hmm,buffer[i]->sequence,  genomic_sequence,buffer[i]->len,  len));
				/* weigth sequence by: 1) mapping quality (scores -sum) and 2) weigth of sequences at the loci.. (priors)   */
				RUN(get_prob_log_Y(hmm,buffer[i]->sequence,  genomic_sequence,buffer[i]->len,len,hmm->alignment_scores[c] + gc[i]->alignment_weigth[c]  ,c));// max ,c);//  hmm->score + hmm->score - sum + ri[i]->priors[hit]  );
				//	RUN(random_model_genome_score(hmm,genomic_sequence,len,prob2scaledprob(1.0)));// -  scaledprob2prob (hmm->alignment_scores[c]))));//+ gc[i]->alignment_weigth[c]))));
				// unaligned genome model update... 
				gc[i]->prior_e[c] = logsum(gc[i]->prior_e[c], hmm->alignment_scores[c]);
				//if(i < 10){
				//	fprintf(stderr,"%s:%d	%f	%f\n",buffer[i]->name,c,hmm->alignment_scores[c],scaledprob2prob(hmm->alignment_scores[c]) );
				//}
			}
			//if(i < 10){
			//	fprintf(stderr,"%s:%s	%f	%f\n",buffer[i]->name,"NA",hmm->unaligned_read_score,scaledprob2prob(hmm->unaligned_read_score));
			//}
			//RUN(random_model_read_score(hmm, buffer[i]->sequence,buffer[i]->len,hmm->unaligned_read_score));
        		//unalign read model update... 
		}
	}
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

	LOG_MSG("Looking at %d-%d.",start,stop);
	for(i = start;i < stop;i++){
		if(!buffer[i]->num_hits){
			RUN(unaligned_to_sam(pw,buffer[i],thread_id));
	        }else{
			RUN(calculate_scores(hmm, buffer[i],gc[i],data->bsd->statmode,1.0)); /* temp is 1 - i.e. no effect... */
			//RUN(alignment_stats(hmm,buffer[i]->num_hits,STAT_ALL));
			max = prob2scaledprob(0.0);			
			k = 0;
			f = -1;
			for(c = 0; c < buffer[i]->num_hits;c++){				
				if(fabs(scaledprob2prob(hmm->alignment_scores[c]) - scaledprob2prob( max)) < 1e-4){
					k++;
				}else if(hmm->alignment_scores[c] > max){
					max = hmm->alignment_scores[c];
					k = 1;
					f = c;
				}
			}
			//uint32_t adler(const void* buf, size_t len)
			if(k > 2){
				f = determine_best_hit_in_case_of_equal_mapq(buffer[i], hmm, max);
			}
			/*if(!strncmp(buffer[i]->name,"ORG457658",9)){
				fprintf(stdout,"K:%d F:%d\n",k,f);
				for(c = 0; c < buffer[i]->num_hits;c++){
					fprintf(stdout,"%s %d %d %f  %f   \n",buffer[i]->name,  i,c, scaledprob2prob(hmm->alignment_scores[c]), FLT_EPSILON);
				}
				}*/
			//ASSERT(f != -1,"best hit cpould not be found...");
			if(hmm->unaligned_read_score  > max){
				RUN(unaligned_to_sam(pw,buffer[i],thread_id));
				//fprintf(stdout,"Unaligned - hits: %d \n", buffer[i]->num_hits);

				//RUN(pw->write(pw,thread_id,"%f\tunaligned\t%s \n",scaledprob2prob(hmm->unaligned_read_score),buffer[i]->name));
				
				//for(c = 0;c < buffer[i]->num_hits;c++){
				//	RUN(pw->write(pw,thread_id,"%f\t%d hits\t%s\n",scaledprob2prob(hmm->alignment_scores[c]), c,buffer[i]->name));
				//}
				//	unaligned_to_sam(ri[i]);
			}else{
				for(j = 0; j < buffer[i]->num_hits;j++){
					flag = 0x100;
					if(j == f){
						flag = 0;
					}
					RUN(get_chr_start_stop(data->bsd->sb_file->si,g_int,buffer[i]->start[j],buffer[i]->stop[j]));
					g_int->start -=ALIGNMENT_FLANKING_LENGTH;
					if( g_int->start < 0){
						g_int->start = 0;
					}
					g_int->stop += ALIGNMENT_FLANKING_LENGTH;
					        
					genomic_sequence = gc[i]->genomic_sequences[j];//  data->sb_file->buffer[i]->genomic_sequences[j];
					len = gc[i]->g_len[j];//   data->sb_file->buffer[i]->g_len[j];
						
					aln = glocal_viterbi_log_Y(hmm,buffer[i]->sequence,  genomic_sequence,buffer[i]->len,  len);
												
					if(g_int->strand){	
						aln = reverse_path(aln);
						RUN(reverse_complement_sequence(genomic_sequence , len));//,"revcomp failed.");
						RUN(reverse_complement_sequence(buffer[i]->sequence , buffer[i]->len));//,"revcomp failed.");
						RUN(reverse(buffer[i]->base_qual, buffer[i]->len));//,"rev failed.");
						flag |= 16;
					}					        
					RUN(align_to_sam(pw, g_int, buffer[i], thread_id, aln, flag,scaledprob2prob(hmm->alignment_scores[j])));
					if(g_int->strand){
						RUN(reverse_complement_sequence(genomic_sequence , len));
						RUN(reverse_complement_sequence(buffer[i]->sequence , buffer[i]->len));//,"revcomp failed.");
						RUN(reverse(buffer[i]->base_qual, buffer[i]->len));//,"rev failed.");
					}
					MFREE(aln);
				}
			}
		}
	}
	pw->flush(pw,thread_id);
	return NULL;
ERROR:
	return NULL;
}


/*

ORG: 
ORG465_chr1_2985727_2985755_+_1	0	chr1	2985728	23	28M	*	0	0	CGAAGGTGTCCAAACTGACAGTGCTGGG	*	XT:A:U	NM:i:1	X0:i:1	X1:i:1	XM:i:1	XO:i:0	XG:i:0	MD:Z:20A7	XA:Z:chr3,-169381438,28M,2;


ORG465_chr1_2985727_2985755_+_1	0	chr1	2985728	40	10M1I7M1I4M1I4M	*	0	0	CCCAGCACTGTCAGTTTGGACACCTTCG	*	NM:i:15	MD:Z:1G0A2G0T0G1C2A0A0C5A1G1T0


>hg19_dna range=chr1:2985728-2985758 5'pad=0 3'pad=0 strand=+ repeatMasking=none
CGAAGGTGTCCAAACTGACAGTGCTGGG
CCCAGCACTGTCAGTTTGGACACCTTCG
CGAAGGTGTCCAAACTGACAATGCTGGGGAG
*/


/* This funtion activates when there are multiple hist swith equal (near equal) scores   */
/* We need to report one hit - if we select one randomly the results will be inconsistent between runs  */
/* If we report the first hit we will introduce bias.  */
/* Therefore I calculate a has value of the read name and hit position and select the hit based on the value */
/* There *shoult be* reproducibility between runs but assignment of multi-mappers will be evenly distributed.... */
int determine_best_hit_in_case_of_equal_mapq(struct sam_bam_entry* entry, struct hmm* hmm, float max)
{
	char buffer[BUFFER_LEN];
	uint32_t hash;
	       
	int i;
	int k = 2;
	int loop_counter = 42;
	int best = -1;
	uint32_t max2 = 0;

	while(k > 1){
		max2 = 0;
		k = 0;
		/* loop through hits - calculate adler for hit pos, take max...  */
		for(i = 0; i < entry->num_hits;i++){
			if(fabs(scaledprob2prob(hmm->alignment_scores[i]) - scaledprob2prob( max)) < 1e-4){
				snprintf(buffer,BUFFER_LEN,"%d%s%" PRId64 "%" PRId64"\n",loop_counter+i,entry->name, entry->start[i],entry->stop[i]);
				hash =  adler(buffer,  strlen (buffer));
				if(hash > max2){
					max2 = hash;
					k = 1;
					best = i;
				}else if( hash == max2){ /* this detect a colision - extremely rare */
					k++;
				}
				//fprintf(stdout,"%d %ud %ud\n",i,hash,max2);
			}
		}
		loop_counter++;
	}
	//fprintf(stdout,"return best: %d\n", best);
	return best;
}
