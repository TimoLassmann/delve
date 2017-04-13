#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <ctype.h>
#include <string.h>
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "tldevel.h"
#include "htsglue.h"

#define MAX_LINE 10000

#define OPT_READLEN 1
#define OPT_AT_START_ERROR 2
#define OPT_AT_END_ERROR 3
#define OPT_INDEL 4
#define OPT_WINDOW_START 5
#define OPT_SERVER 6

#define RUNMODE_SIM 1
#define RUNMODE_EVAL 2

struct parameters {
	FILE* foutfile;
	char* outfile;
	char* bedfile;
	char* control_file;
	int infiles;
	int quiet_flag;
	int num_threads;
	unsigned int sim_length;
	unsigned int seed;
	float sim_start_error;
	float sim_end_error;
	float sim_mismatch_prob;
	float sim_GG_error_addition;
	float random_model_prior;
	int sim_numseq;
	int window_start;
	int server;
	int help;
	int mode;
};

int eval_sambam(struct parameters* param);
int simulate_reads_from_bed(faidx_t*  index ,struct parameters* param);

void usage();
int mutate_sequence(char* seq,struct parameters*  param, int* errors);
void convert_sequence(char* seq, int len);

int compare_genome_interval_mapping(struct genome_interval* a, struct genome_interval* b, int* res);
float* shuffle_arr_r(float* arr,int n,unsigned int* seed);

int main (int argc, char * argv[])
{
	int c;
	struct parameters* param = NULL;
	faidx_t*  index = NULL;
	char* fasta_file = NULL;
	char* mode_string = NULL;
	tlog.echo_build_config();

	MMALLOC(param, sizeof(struct parameters));
	param->mode = RUNMODE_SIM;
	param->help = 0;
	param->control_file = 0;
	param->bedfile = 0;
	param->outfile = 0;
	param->sim_length = 30;
	param->sim_start_error = 0.02f;
	param->sim_end_error = 0.03f;
	param->sim_mismatch_prob = 0.8f;
	param->sim_GG_error_addition = 0.0f;
	param->window_start = 0;
	param->sim_numseq = 1000000;
	param->server = 0;
	param->seed = (unsigned int) (time(NULL) * (42));
	while (1){
		static struct option long_options[] ={
			{"read-len",required_argument,0,OPT_READLEN},
			{"start-error",required_argument,0,OPT_AT_START_ERROR},
			{"stop-error",required_argument,0,OPT_AT_END_ERROR},
			{"mismatch-rate",required_argument,0,OPT_INDEL},
			{"sw",required_argument,0,OPT_WINDOW_START},
			{"server",0,0,OPT_SERVER},
			{"quiet",required_argument,0,'q'},
			{"help",0,0,'h'},
			{0, 0, 0, 0}
		};
		
		int option_index = 0;
		c = getopt_long_only (argc, argv,"o:qhn:",long_options, &option_index);
		
		if (c == -1){
			break;
		}
		
		switch(c) {
		case 0:
			break;
		case OPT_SERVER:
			param->server = 1;
			break;
		case OPT_READLEN:
			param->sim_length = atoi(optarg);
			break;
		case OPT_AT_START_ERROR:
			param->sim_start_error = atof(optarg);
			break;
		case OPT_AT_END_ERROR:
			param->sim_end_error = atof(optarg);
			break;
		case OPT_INDEL:
			param->sim_mismatch_prob = atof(optarg);
			break;
		case OPT_WINDOW_START:
			param->window_start = atoi(optarg);
			break;
		case 'n':
			if(optarg[0] == '-'){
				fprintf(stderr,"option -n requires an argument\n");
				exit(-1);
			}
			param->sim_numseq = atoi(optarg);
			break;
		case 'o':
			if(optarg[0] == '-'){
				fprintf(stderr,"option -o requires an argument\n");
				exit(-1);
			}
			param->outfile = optarg;
			break;
		case 'c':
			if(optarg[0] == '-'){
				fprintf(stderr,"option -c requires an argument\n");
				exit(-1);
			}
			param->control_file = optarg;
			break;
		case 'h':
			param->help = 1;
			break;
		case '?':
			exit(1);
			break;
		default:
			fprintf(stderr,"default\n\n\n\n");
			abort ();
		}
	}
	
	if(param->help || (argc - optind) < 1){
		usage();
		MFREE(param);
		return EXIT_SUCCESS;
	}
	/* 0. there is an outfile */
	if(param->outfile == NULL){
		ERROR_MSG("No outfile (-o) given.");
	}

	mode_string= argv[optind++];
	/* 1. mode  */

	param->mode = -1;
	if(strncmp(mode_string,"sim",3) == 0){
		param->mode = RUNMODE_SIM;
		LOG_MSG("Running in sim mode.");

		fasta_file = argv[optind++];
		/* 2. does fastafile exists and can it be opened...  */
		/* this will open the index and fail properly if not successfull */
		RUNP(index = get_faidx(fasta_file));
		LOG_MSG("Loaded \"%s\" faidx index.",fasta_file);
		
		param->bedfile = argv[optind++];
		/* 3. is bedfile there and readable */
		if(my_file_exists(param->bedfile) == 0){
			ERROR_MSG("File \"%s\" does not exist.",param->bedfile);
		}
		RUN(simulate_reads_from_bed(index,param));
	}
	if(strncmp(mode_string,"eval",4) == 0){
		param->mode = RUNMODE_EVAL;
		LOG_MSG("Running in eval mode.");

		param->bedfile = argv[optind++];
		if(my_file_exists(param->bedfile) == 0){
			ERROR_MSG("File \"%s\" does not exist.",param->bedfile);
		}
		RUN(eval_sambam(param));
	}	
	if(param->mode  == -1){
		ERROR_MSG("Mode \"%s\" not found", mode_string);
	}
	
	if(index){
		fai_destroy(index);
	}
	if(param){
		MFREE(param);
	}
	
	return EXIT_SUCCESS;
ERROR:
	if(index){
		fai_destroy(index);
	}
	if(param){
		MFREE(param);
	}
	return EXIT_FAILURE;
}

void usage()
{
	fprintf(stdout, "Usage:   bedsim [options] <tomedb> <in.bed> -o <suffix> \n\n");
	fprintf(stdout, "Options:\n");
	fprintf(stdout,"%15s%10s%7s%-30s\n","-read-len","INT" ,"", "Length of simulated reads [30].");
	fprintf(stdout,"%15s%10s%7s%-30s\n","-start-error","FLT" ,"", "Error rate at 5' end [0.02].");
	fprintf(stdout,"%15s%10s%7s%-30s\n","-stop-error","FLT" ,"", "Error rate at 3' end [0.03].");
	fprintf(stdout,"%15s%10s%7s%-30s\n","-mismatch-rate","FLT" ,"", "Fraction of simulated mismatches [0.8].");
	fprintf(stdout,"%15s%10s%7s%-30s\n","-sw","INT" ,"", "Simulate reads from +- INT start of feature [0].");
	fprintf(stdout,"%15s%10s%7s%-30s\n","-n","INT" ,"", "Number of reads to simulate [1000000].");
	fprintf(stdout,"%15s%10s%7s%-30s\n","-o","STR" ,"", "Output file name.");
	fprintf(stdout, "\n");	
}

int simulate_reads_from_bed(faidx_t*  index ,struct parameters* param)
{	
	struct genome_interval* g_int = NULL;
	struct seq_info* si = NULL;
	FILE* file = 0;
//	FILE* local_out_file = 0;
	FILE* outfile = 0;
	char line[MAX_LINE];
//	char seq_name[MAX_LINE];
	char* seq = NULL;
	int num_sequences = param->sim_numseq;
	int column = 0;
	int i,j,tmp;
	int start;
	int stop;
//	unsigned int chr_start;
	int strand = 0;
	int c = 0;
//	int hit = 0;
	int rank = 0;
	float y = 0.0f;
	float total = 0.0f;
	int errors;
	
	float* rank_arr = NULL;
	
	unsigned int sim_start;
	int seq_num = 1;
//	char* local_outfile;
//	char* local_control;
	RUNP(g_int = init_genome_interval(0,0,0));	
	RUNP(si =  make_si_info_from_fai(index));

	RUNP(file = fopen(param->bedfile , "r" ));
	RUNP(outfile = fopen(param->outfile , "w" ));
	while(fgets(line, MAX_LINE, file)){
		if(line[0] != '#'){
			total++;
		}
	}
	rewind(file);

	MMALLOC(rank_arr, sizeof(float) * total);
		
	for(i = 0; i < total;i++){
		rank_arr[i] = i+1;
	}
	RUNP(rank_arr = shuffle_arr_r(rank_arr ,total,&param->seed));
	y = 0.0;
	for(i = 0; i < total;i++){
		rank_arr[i] = (float)rank_arr[i] / (  pow((float) rank_arr[i], 1.2f + ((float) ((int)rank_arr[i] % (int)(total+1)) / total)));
		y+= rank_arr[i] ;
	}
	
	c = 0;
	for(i = 0; i < total;i++){
		rank_arr[i]  = (int)(rank_arr[i]  / y * (float)num_sequences);
		c += rank_arr[i] ;
	}
	
	while(c < num_sequences){
		rank_arr[(int) rand_r(&param->seed) % (int)total]++;
		c++;
	}
	
	for(c = 0;c < 1;c++){
		/*j = (int)strlen(param->outfile);
		local_outfile = malloc(sizeof(char)*(j+9));
		for(i = 0; i < j ;i++){
			local_outfile[i] = param->outfile[i];
		}
		local_outfile[j] = '.';
		local_outfile[j+1] = 'r';
		local_outfile[j+2] = 'e';
		local_outfile[j+3] = 'p';
		local_outfile[j+4] = (char)  c + 49;
		local_outfile[j+5] = '.';
		local_outfile[j+6] = 'f';
		local_outfile[j+7] = 'a';
		local_outfile[j+8] = 0;
		if (!(local_out_file = fopen(local_outfile , "w" ))){
			fprintf(stderr,"Cannot open bed file '%s'\n",local_outfile);
			exit(-1);
		}
		*/	
		rank = 0;
		
		while(fgets(line, MAX_LINE, file)){
			if(line[0] != '#'){
				column = 1; //<QNAME>
				tmp = 0;
				//hit = 0;
				for(j = 0;j < MAX_LINE;j++){
					tmp++;
					if(isspace((int)line[j])){
						break;
					}
				}
				
				for(j = 0;j < FIELD_BUFFER_LEN;j++){
					
					if(isspace((int)line[j])){
						g_int->chromosome[j] =0; 
						break;
					}
					g_int->chromosome[j] = line[j];
			        
				}
				for(i = 0; i < MAX_LINE;i++){
					if(line[i] == '\n'){
						break;
					}
					if(isspace((int)line[i])){
						column++;
						switch(column){
							case 2: // <FLAG>
								g_int->start = atoi(line + i + 1);
								break;
							case 3: // <RNAME>
								g_int->stop = atoi(line + i + 1);;
								break;
							case 6: // <RNAME>
								if(*(line+i+1) == '+'){
									g_int->strand = 0;
								}else{
									g_int->strand = 1;
								}
								break;
						}
					}
				}
				//g_int->stop = param->sim_length + g_int->start;
				for(j = 0; j < (int)(rank_arr[rank] );j++){
					start = g_int->start;
					stop = g_int->stop;
					if(!param->window_start){
						sim_start = (unsigned int) ( rand_r(&param->seed) % ((int)((g_int->stop-(param->sim_length+1))  - g_int->start)));
						sim_start += g_int->start;
						g_int->start = sim_start;
						g_int->stop = sim_start +param->sim_length;
				        }else{
						sim_start = (unsigned int) (rand_r(&param->seed) %((int)(2*param->window_start)));
						
						if(!strand){
							sim_start = sim_start - param->window_start;
							sim_start += g_int->start;
							g_int->start = sim_start;
							g_int->stop = sim_start +param->sim_length;
						}else{
							sim_start =  sim_start - param->window_start;
							sim_start += g_int->stop;
							g_int->start = sim_start;
							g_int->stop = sim_start +param->sim_length;
						}
					}

					DPRINTF3("lookign for: %s:%d-%d %c",g_int->chromosome,g_int->start,g_int->stop, "+-"[g_int->strand]);
					seq = get_sequence(index,g_int);
					errors = 0;
					RUN(mutate_sequence(seq,param ,&errors));

					fprintf(outfile,">ORG%d_%s_%d_%d_%c_%d\n%s\n",seq_num, g_int->chromosome,g_int->start,g_int->stop, "+-"[g_int->strand],errors,seq);
					seq_num++;
					MFREE(seq);
					g_int->start = start;
					g_int->stop = stop;
		 		}
				rank = rank + 1;
			}
		}
		rank = 0;	
		
	}
	fclose(file);
	fclose(outfile);
	free(rank_arr);
	free_sequence_info(si);
	free_genome_interval(g_int);
	return OK;
ERROR:
	return FAIL;

}


void convert_sequence(char* seq, int len)
{
	char alphabet[] = "ACGTN";
	int i;
	for(i = 0; i < len;i++){
		seq[i] = alphabet[(int)seq[i]];
	}
	seq[len] = 0;
}




int mutate_sequence(char* seq,struct parameters*  param,int* errors)
{
	int j,c;
	char tmp_seq[MAX_LINE];
	char alphabet[] = "ACGTN";
	float a,b, local_c;
	char subc;
	
	b  = ((param->sim_end_error - ((float)param->sim_length *  param->sim_start_error)) / (float)(param->sim_length -1)) * -1.0f;

	
	a = param->sim_start_error - b;
	c = 0;
	for(j = 0; j < param->sim_length;j++){
		local_c =0.0;
		if(j){
			if(seq[j-1] == 2 && seq[j] == 2){
				local_c += param->sim_GG_error_addition;
			}
		}
		if((float)rand_r(&param->seed)/ (float) RAND_MAX  <= (float)(j+1)* a +b + local_c){
			*errors = *errors + 1;
			if((float)rand_r(&param->seed)/(float)RAND_MAX  <=  param->sim_mismatch_prob){
				subc = alphabet[(int)rand_r(&param->seed)  % (4)];
				while(subc == seq[j]){
					subc = alphabet[(int)rand_r(&param->seed)  % (4)]; //subc = alphabet[random_int(4)];
				}
				tmp_seq[c] = subc;
				c++;
//				fprintf(stderr,"MISMATCH\n");
			}else{
				if(rand_r(&param->seed) <= 0.5f){
					tmp_seq[c] = alphabet[(int)rand_r(&param->seed)  % (4)];
					c++;
					tmp_seq[c] = seq[j];
					c++;
//						fprintf(stderr,"INSERT\n");
				}else{
//						fprintf(stderr,"DELETION\n");
				}
			}
		}else{
			tmp_seq[c] = seq[j];
			c++;
		}
	}
	for(j = 0; j < c;j++){
		seq[j] = tmp_seq[j];
	}
	seq[c] = 0;
	return OK;
}

float* shuffle_arr_r(float* arr,int n,unsigned int* seed)
{
	int i,j;
	float tmp;
	
	
	for (i = 0; i < n - 1; i++) {
		j = i +  (int) (rand_r(seed) % (int) (n-i));
		tmp = arr[j];
		arr[j] = arr[i];
		arr[i] = tmp;
	}
	return arr;
}




int eval_sambam(struct parameters* param)
{
	struct sam_bam_file* sb_file = NULL;
	struct sam_bam_entry* entry = NULL;
	struct genome_interval* g_int = NULL;
	struct genome_interval* g_org_int = NULL;
	struct results{
		int mapped[6];
		int mismapped[6];
		int error[6][6];
		int total_mapped;
		int unmapped;
	} res;
	int buffer_size = 1000000;
	int i = 0;
	int j = 0;
	int norm_qual = 0;
	int comparison = 0;
	int real_error = 0;
	char c = 0;


	/* Init res */
	res.unmapped = 0;
	res.total_mapped = 0;
	for(i = 0; i < 6;i++){
		res.mapped[i] = 0;
		res.mismapped[i] = 0;
		for(j = 0; j < 6;j++){
			res.error[i][j] = 0;
		}
	}
	
	RUNP(g_int = init_genome_interval(0,0,0));	

	RUNP(g_org_int = init_genome_interval(0,0,0));	

	
	RUNP(sb_file = open_SAMBAMfile(param->bedfile,buffer_size,10,0,0));
	echo_header(sb_file);
	while(1){
		RUN(read_SAMBAM_chunk(sb_file,1,0));
		//DPRINTF2("read:%d",sb_file->num_read );
		//RUN(read_sam_bam_chunk(infile,data->si, buffer,1,&numseq),"Read sam/bam chunk in thread %d failed", data->threadID);
		//if((status = read_sam_bam_chunk(infile,data->si, buffer,1,&numseq)) != kslOK)  exit(status);
		if(!sb_file->num_read){
			break;
		}
		res.total_mapped +=sb_file->num_read;
		res.unmapped += sb_file->total_entries_in_file - sb_file->num_read;
		for(i = 0; i < sb_file->num_read;i++){
			entry = sb_file->buffer[i];
	
			/*	 	int64_t* start;
			int64_t* stop;
			char* sequence;
			char* name;
			int qual;
			int len;
			int max_len;
			int num_hits;
			int max_num_hits;
			*/

			
			if(entry->num_hits){
				g_int->g_start = entry->start[0];
				g_int->g_stop = entry->stop[0];
				//fprintf(outfile,">ORG_%s_%d_%d_%c_%d\n%s\n",g_int->chromosome,g_int->start,g_int->stop, "+-"[g_int->strand],errors,seq);
				//ORG_chr22_51061767_51061797_+_1
				sscanf(entry->name,"%*[^_]_%[^_]_%d_%d_%c_%d",g_org_int->chromosome,&g_org_int->start,&g_org_int->stop, &c,&real_error);
				if(c == (int) '+'){
					g_org_int->strand = 0;
				}else{
					g_org_int->strand = 1;
				}
				RUN(internal_to_chr_start_stop_strand(sb_file->si,g_int));
//				RUN(get_chr_start_stop(sb_file->si  ,g_int,entry->start[0], entry->stop[0]));
//				fprintf(stdout,"%s\n%s\t%d\t%d\t%c\n", entry->name,g_int->chromosome,g_int->start,g_int->stop, "+-"[g_int->strand]);
//				fprintf(stdout,"%s\t%d\t%d\t%c\n", g_org_int->chromosome,g_org_int->start,g_org_int->stop, "+-"[g_org_int->strand]);
//				fprintf(stdout,"%d %f\n",entry->qual, (roundf((float) entry->qual /10.0)));
				norm_qual =  roundf((float) entry->qual /10.0);
				if(norm_qual > 5){
					norm_qual = 5;
				}
				if(real_error > 5){
					real_error = 5;
				}
					      
				RUN(compare_genome_interval_mapping(g_org_int,g_int,&comparison));
				if(comparison < 5){ /* distance criteria to declare a match */
					comparison = 1;
				}else{
					comparison = 0;
				}
				/* record stars */
				if(comparison){
					res.mapped[norm_qual]++;
					res.error[real_error][norm_qual]++;
				}else{
					res.mismapped[norm_qual]++;
				}
				
//			      fprintf(stdout,"%d %f RES:%d  error:%d \n",entry->qual, (roundf((float) entry->qual /10.0)), comparison ,real_error);
				
			}
		}
	}
	
	fprintf(stdout,"%d\n%d\n%f\n",res.total_mapped,res.unmapped,roundf ((float)( res.unmapped) / (float) (res.total_mapped+res.unmapped) *10000.0)/ 100.0);
	
	for(i = 0; i < 6;i++){
		fprintf(stdout,"%d\t",i*10);
	}
	fprintf(stdout,"mapping quality\n");
	for(i = 0; i < 6;i++){
		fprintf(stdout,"%d\t",res.mapped[i]);
	}
	fprintf(stdout,"correctly mapped\n");
	for(i = 0; i < 6;i++){
		fprintf(stdout,"%d\t",res.mismapped[i]);
	}
	fprintf(stdout,"in-correctly mapped\n");
	for(j = 0; j < 6;j++){
		for(i = 0; i < 6;i++){
			fprintf(stdout,"%d\t",res.error[j][i]);
		}
		fprintf(stdout,"%d error\n",j);
	}
	RUN(close_SAMBAMfile(sb_file));
	free_genome_interval(g_org_int);
	free_genome_interval(g_int);
	return OK;
ERROR:
	if(sb_file){
		RUN(close_SAMBAMfile(sb_file));
	}
	return FAIL;
}


int compare_genome_interval_mapping(struct genome_interval* a, struct genome_interval* b, int* res)
{
	if(a->strand != b->strand){
		*res = -1;
		return OK;
	}

	if(strcmp(a->chromosome, b->chromosome) != 0){
		*res = -1;
		return OK;
	}
	
	*res = MACRO_MAX(abs(a->start - b->start),abs( a->stop - b->stop)); 
	return OK;
}


