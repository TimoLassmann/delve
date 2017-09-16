#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <sys/stat.h>
#include <unistd.h>
#include <dirent.h>
#include <sys/wait.h>
//#include <string.h>
// #include <inttypes.h>
//#include <ctype.h>
//#include <float.h>


#include "tldevel.h"
#include "thr_pool.h"

#define OPT_WORKDIR 1
#define OPT_NUMTHREADS 2
#define OPT_GENOME 3
#define OPT_MAX_MEM 4
#define OPT_MAX_THR 5
#define OPT_OUT 6

struct parameters{
	char* work_dir;
	char* genome;
	char* out_file;
	int max_mem;
	int max_threads;
	int num_threads;
};

int run_pipeline(struct parameters* param);


int set_num_parallel(int num_in,int threads,int max_mem,int req_mem);
int set_threads_per_job(int num_p, int threads);

int print_help(char **argv);

char** get_files_names(char* directory, char* suffix, int* num_files);

int create_dir(char* name, int force);

int main (int argc, char *argv[])
{		
	struct stat st = {0};
	struct parameters* param = NULL;
	int c;
	tlog.echo_build_config();
	MMALLOC(param, sizeof(struct parameters));
	param->work_dir = NULL;
	param->num_threads = 4;
	param->max_mem = 8;
	param->max_threads = 4;
	param->out_file = NULL;
	while (1){	
		static struct option long_options[] ={
			{"threads",required_argument,0,OPT_NUMTHREADS},
			{"genome",required_argument,0,OPT_GENOME},
			{"maxm",required_argument,0,OPT_MAX_MEM},
			{"maxt",required_argument,0,OPT_MAX_THR},
			{"in",required_argument,0,OPT_WORKDIR},
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
		case OPT_WORKDIR:
			param->work_dir = optarg;
			break;
		case OPT_OUT:
			param->out_file = optarg;
			break;
		
		case OPT_GENOME:
			param->genome = optarg;
			break;
		case OPT_NUMTHREADS: 
			param->num_threads = atoi(optarg);
			break;
		case OPT_MAX_MEM: 
			param->max_mem = atoi(optarg);
			break;
		case OPT_MAX_THR: 
			param->max_threads = atoi(optarg);
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

	if(param->work_dir == NULL){
		ERROR_MSG("Direcory: %s does not exist.\nTry --help");
	}	
	if(param->out_file  == NULL){
		ERROR_MSG("No outfile! use -out XXX\nTry --help");
	}
	
	if(param->genome == NULL){
		ERROR_MSG("Direcory: %s does not exist.\nTry --help");
	}
	if(stat(param->work_dir, &st) == -1) {
		ERROR_MSG("Direcory: %s does not exist.\nTry --help");
		
	}

	if(param->num_threads > param->max_threads){
		LOG_MSG("You asked for more threads than max (%d). increase using -maxt XXX", param->max_threads);
		param->num_threads = param->max_threads;
	}
	
	LOG_MSG("Working on: %s   %d",param->work_dir,param->num_threads);

	RUN(run_pipeline(param));
	MFREE(param);
	return EXIT_SUCCESS;
ERROR:
	fprintf(stdout,"\n  Try run with  --help.\n\n");
	if(param){
		MFREE(param);
	}
	return EXIT_FAILURE;
}

int run_pipeline(struct parameters* param)
{
	struct thr_pool* pool = NULL;		
       	char** input_files = 0;
	FILE* fpr_out = NULL;
	int num_input_files = 0;
	int i;
	int t;
	int thread_break = 0;
	int threads_per_job = 1;

	/* Step 0 - find input files...  */
	
	RUNP(input_files = get_files_names(param->work_dir,".fa",&num_input_files));
	RUNP(fpr_out = fopen(param->out_file, "w"));
	
	/* bwa commands.. */
	/* assume bwa uses 4G  */
	thread_break =  set_num_parallel(num_input_files , param->num_threads,param->max_mem, 4);
	threads_per_job = set_threads_per_job(thread_break,param->num_threads);
	
	t = 0;
	for(i = 0 ; i < num_input_files;i++){
		fprintf(stdout,"%s\n",input_files[i]);
		fprintf(fpr_out,"bwa aln -n 3 -o 3  -t %d %s %s | ",threads_per_job, param->genome, input_files[i]);
                fprintf(fpr_out,"bwa samse %s  - %s -n 100   | ", param->genome, input_files[i]);
		fprintf(fpr_out,"samtools sort -@ %d -O bam -l 9 -T ~/tmp - >  %s.bwa.bam ",threads_per_job,input_files[i]);
		t++;
		if(t == thread_break){
			fprintf(fpr_out,"& \n wait;\n");
			fprintf(stdout,"jobs left: %d\n",(num_input_files-1) - i);
			thread_break =  set_num_parallel((num_input_files-1) -i , param->num_threads,param->max_mem, 4);
			threads_per_job = set_threads_per_job(thread_break,param->num_threads);	
			t = 0;
		}else{
			fprintf(fpr_out,"& \n");
		}
	}
	
	fprintf(fpr_out," wait;\n");
        /* delve commands.. */
	/* assume delve  uses 2G  */
	thread_break =  set_num_parallel(num_input_files , param->num_threads,param->max_mem, 2);
	threads_per_job = set_threads_per_job(thread_break,param->num_threads); 
	fprintf(stdout,"delve: %d jobs in parallel... \n", thread_break);
	
	int siter = 5;
	int giter = 5;
	
	t = 0;
			
	thread_break =  set_num_parallel(num_input_files , param->num_threads,param->max_mem, 2);
	threads_per_job = set_threads_per_job(thread_break,param->num_threads); 
	for(i = 0 ; i < num_input_files;i++){
		fprintf(fpr_out ,"delve  -siter %d -giter %d  -t %d  -d  %s %s.bwa.bam " ,siter,giter, threads_per_job, param->genome,input_files[i]);
		t++;
		if(t == thread_break){
			fprintf(fpr_out,"& \n wait;\n");
			fprintf(stdout,"jobs left: %d\n",(num_input_files-1) - i);
			thread_break =  set_num_parallel((num_input_files-1) -i , param->num_threads,param->max_mem, 4);
			threads_per_job = set_threads_per_job(thread_break,param->num_threads);
			
			t = 0;
		}else{
			fprintf(fpr_out,"& \n");
		}
	}
	
	fprintf(fpr_out," wait;\n");
	
	for(i = 0; i < num_input_files;i++){
		MFREE(input_files[i]);
	}
	MFREE(input_files);
	
	if(pool){
		thr_pool_destroy(pool);
	}
	return OK;
ERROR:
	if(pool){
		thr_pool_destroy(pool);
	}
	return FAIL;
}

int set_num_parallel(int num_in,int threads,int max_mem,int req_mem)
{
	int thread_break = 0;
	thread_break = threads;
	if( num_in < thread_break){
		thread_break =  num_in;
	}
	if(thread_break * req_mem > max_mem){
		thread_break = (int) floor((float) max_mem  / (float) req_mem);
	}
	return thread_break;
}

int set_threads_per_job(int num_p, int threads)
{
	int threads_per_job = 1;
	if(num_p  < threads){
		threads_per_job = (int)ceil((float) threads / (float)num_p);
	}
	return threads_per_job;
}


int print_help(char **argv)
{
	//const char description[] = "Aligns HMMs to new sequences.";
	const char usage[] = " -in <fasta dir> -out <script.sh> -maxm 500 -maxt 96 -t 48";
	fprintf(stdout,"\nUsage: %s [-options] %s\n\n",basename(argv[0]) ,usage);	
	fprintf(stdout,"Options:\n\n");
	
	fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--in","Target directory ." ,"[?]"  );
	fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--out","Shell script ." ,"[?]"  );
	fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--genome","which genome to align against." ,"[?]"  );
	fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--maxm","Max memory on system (in GB)." ,"[8]"  );
	fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--maxt","Number of cores on machine ." ,"[4]"  );
	fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--threads","Number of threads to use." ,"[4]"  );
	return OK;
}	


char** get_files_names(char* directory, char* suffix, int* num_files)
{
	char** list = NULL;
	
	ASSERT(directory !=NULL," directory is NULL");
	ASSERT(suffix != NULL, "suffix is NULL");
	ASSERT(strlen(suffix) >1, "suffic is too short: %s", suffix);
	
	DIR *dp;
	struct dirent *entry;
	struct stat statbuf;
	int len;

	*num_files = 0;
        // get current directory...
	//RUNP(getcwd(cwd, sizeof(cwd)));
	//open directory in question to count number of files...
	RUNP(dp = opendir(directory));
	
	if(chdir(directory)){
		ERROR_MSG("chdir failed.");
	}
	
	while((entry = readdir(dp)) != NULL) {
		lstat(entry->d_name,&statbuf);
		if(S_ISDIR(statbuf.st_mode)) {
			/* Found a directory, but ignore . and .. */
			if(strcmp(".",entry->d_name) == 0 ||
			   strcmp("..",entry->d_name) == 0)
				continue;
		}else{
			if(!strcmp(suffix, entry->d_name+ (strlen(entry->d_name) - strlen(suffix)))){
				MREALLOC(list, sizeof(char*) * ( *num_files + 1));
				//fprintf(stdout,"%p List\n", list);
				list[ *num_files] = NULL;
				len = strlen(entry->d_name)+strlen(directory)  +2;
				MMALLOC(list[ *num_files], sizeof(char) *len);
				snprintf(list[*num_files],len,"%s/%s",directory,entry->d_name);
				
				//strcpy(list[ *num_files], entry->d_name);
				//fprintf(stdout,"%p List\n", list);
				
				//for(j = 0; j <  *num_files+1;j++){
				//	fprintf(stdout,"%d\t%s\n",j,list[j]);
				//}
				*num_files = *num_files + 1;
			}
		}
	}
	if(closedir(dp)){
		ERROR_MSG("closedir failed.");
	}
	if(chdir("..")){
		ERROR_MSG("chdir failed.");
	}
	return list;
ERROR:
	return NULL;
}


int create_dir(char* name, int force)
{
	struct stat st = {0};
	
	if(stat(name, &st) == -1) {
		RUN(mkdir(name, 0700));
	}else{
		if(!force){
			ERROR_MSG("Directory %s already exists.",name);
		}
	}
	return OK;
ERROR:
	return FAIL;
}


