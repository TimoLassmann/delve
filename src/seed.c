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

#include "thr_pool.h"
#include "tldevel.h"

#include  <pthread.h>

#define OPT_INPUT 1
#define OPT_GENOME 2
#define OPT_OUTPUT 3
#define OPT_SEED_LEN 4
#define OPT_SEED_STEP 5
#define OPT_NUMQUERY 6
#define OPT_NTHREADS 7


#define MAX_LOOP_SEARCH_DEPTH 2

#define LIST_STORE_SIZE 10
#define FILE_END 0xFFFFFFFFFFFFFFFFu

#define MAX_LINE 512

int nuc_code[256];
unsigned int nuc_code5[256];
unsigned int rev_nuc_code5[16];

void init_nuc_code(void);
unsigned char* reverse_without_complement(unsigned char* p,int len);

struct chromosome{
        uint8_t* seq;
        char* name;
        int malloc_len;
        int seq_len;
};

struct genome{
        struct chromosome** chromosomes;
        int malloc_num;
        int num_chr;        
};

struct parameters{
        struct thr_pool* pool;        
        char* input;
        char* genome;
        char* output;
        int seed_len;
        int seed_step;
        int num_query;
        int num_threads;
};

int print_help(char **argv);
struct genome* load_genome(struct parameters* param);
void free_genome(struct genome* genome);




struct seq_info{
        unsigned long int* match_units;
        unsigned char* sn;
        unsigned char* seq;
        unsigned char* reverse_seq;
        unsigned char* qual;
        unsigned int id_num;
        int len;
        long int last_f_test;
        long int last_r_test;
};

struct qs_struct{
        struct seq_info** seq_info;
        unsigned char* string;	
        int* hash;
        int* hash_index;
        int** hash_t;
        int** hash_index_t;
        int alloc_num_query;
        int string_len;
        int q_num;
        int size;
        int min_len;
        int max_len;
};

struct seed_thread_data{
        struct qs_struct* qs;
        struct parameters* param;
        long int len;
        long int total_len;
        unsigned char* chromsome;
        int start;
        int end;
        int thread;
};


int seed_controller(struct parameters* param);

struct qs_struct* init_query_structure(struct parameters* param);
void free_query_structure(struct qs_struct * qs);

struct qs_struct* search_fasta_fastq(struct qs_struct* qs,struct parameters* param);
unsigned long int read_fasta_fastq(FILE *file, struct qs_struct* qs,struct parameters* param,unsigned long int position);
void reverse_complement(unsigned char* p,unsigned char* rc,int len);

struct qs_struct* make_hash(struct qs_struct* qs,int r);
struct qs_struct* search_hash(struct qs_struct* qs,unsigned char* t,int n, long int r,long int offset);
void print_sam(struct seq_info* si, unsigned char** t_names, long int* lengths,int size);

int seed_controller_thread(struct parameters* param);
struct qs_struct* search_fasta_fastq_thread(struct qs_struct* qs,struct parameters* param);
struct qs_struct* make_hash_thread(struct qs_struct* qs,int r,int step,int thread,int start,int stop);
void* search_hash_thread(void *threadarg);

struct qs_struct* make_B(struct qs_struct* qs,int n);
struct qs_struct* read_fasta_fastq2(FILE *file, struct qs_struct* qs,struct parameters* param) ;
unsigned int rev_bits(unsigned int x);
unsigned long int validate_bpm(unsigned char* t,unsigned char* p,int n,int m,long int offset);
unsigned int BinaryInsertionSortOne(unsigned long int *arr, unsigned int size, unsigned long int new);


int main (int argc, char *argv[]) 
{		
        struct parameters* param = NULL;
        int c;

        tlog.echo_build_config();

        MMALLOC(param, sizeof(struct parameters));
        param->input = NULL;
        param->genome = NULL;
        param->output = NULL;
        param->pool = NULL;
        param->seed_len = 12;
        param->seed_step = 8;
        param->num_query = 1000000;
        param->num_threads = 8;
        
        while (1){	
                static struct option long_options[] ={
                        {"in",required_argument,0,OPT_INPUT},
                        {"genome",required_argument,0,OPT_GENOME},
                        {"seed-len",required_argument,0,OPT_SEED_LEN},			
                        {"seed-step",required_argument,0,OPT_SEED_STEP},
                        {"nthread",required_argument,0,OPT_NTHREADS},
                        {"output",required_argument,0,OPT_OUTPUT},			       
                        {"help",0,0,'h'},
                        {0, 0, 0, 0}
                };
		
                int option_index = 0;
                c = getopt_long_only (argc, argv,"h",long_options, &option_index);
		
                if (c == -1){
                        break;
                }
		
                switch(c) {
                case OPT_INPUT:
                        param->input = optarg;
                        break;
                case OPT_OUTPUT:
                        param->output = optarg;
                        break;                                  
                case OPT_GENOME:
                        param->genome = optarg; 
                        break;
                case OPT_SEED_LEN:
                        param->seed_len = atoi(optarg);
                        break;
                case OPT_SEED_STEP:
                        param->seed_step = atoi(optarg);
                        break;
                case OPT_NTHREADS:
                        param->num_threads = atoi(optarg);
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
        	
        LOG_MSG("Starting run");

        if(!param->input){
                RUN(print_help(argv));
                ERROR_MSG("No input file! use --in <blah.fa>");
                
        }else{
                if(!my_file_exists(param->input)){
                        RUN(print_help(argv));
                        ERROR_MSG("The file <%s> does not exist.",param->input);               
                }           
        }
        if(!param->genome){
                RUN(print_help(argv));
                ERROR_MSG("No genome file! use -genome <blah.fa>\n");
        }else{
                if(!my_file_exists(param->genome)){
                        RUN(print_help(argv));
                        ERROR_MSG("The file <%s> does not exist.",param->genome);               
                }           
        }
        if(param->seed_len < 6){
                RUN(print_help(argv));
                ERROR_MSG("Seed len has to be at least 6.");   
        }
        if(param->seed_len < 12){
                RUN(print_help(argv));
                WARNING_MSG("With a seed len of %d the processing may take a while....",param->seed_len);     
        }        
        ASSERT(param->seed_step >= 1,"Step is too small.");
        init_nuc_code();
                
        /* If(!param->num_infiles){ */
        /*         RUN(print_help(argv)); */
        /*         ERROR_MSG("delve requires at least one input\n"); */
        /* } */
	
        /* RUN(run_delve(param)); */
	      RUN(seed_controller_thread(param));
        if(param){
                if(param->pool){
                        thr_pool_destroy(param->pool);
                }
                MFREE(param);
        }
        return EXIT_SUCCESS;
ERROR:
        fprintf(stdout,"\n  Try run with  --help.\n\n");
        if(param){
                if(param->pool){
                        thr_pool_destroy(param->pool);
                }
                MFREE(param);
        }
        return EXIT_FAILURE;
}

int print_help(char **argv)
{
        const char usage[] = " -in <fasta> -genome <genome>";
        fprintf(stdout,"\nUsage: %s [-options] %s\n\n",basename(argv[0]) ,usage);	
        fprintf(stdout,"Options:\n\n");
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--seed-len","Length of seeds." ,"[12]"  );
      	fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--seed-step","Distance between seeds." ,"[8]"  );
      	fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--nthread","Number of threads." ,"[8]"  );
      	fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--output","Output file name." ,"[?]"  );
        return OK;
}

int seed_controller_thread(struct parameters* param)
{
        struct qs_struct* qs = NULL;
        int t;
        int interval;
        FILE *fastafile = 0; 
        struct genome* genome = NULL;

        if((param->pool = thr_pool_create(param->num_threads, param->num_threads, 0, 0)) == NULL) ERROR_MSG("Creating pool thread failed.");


        LOG_MSG("Loading genome");
        RUNP(genome = load_genome(param));
        
        RUNP(qs = init_query_structure(param));
        
        RUNP(fastafile = fopen(param->input, "r" ));
        
        qs->size = param->num_query;
        while(qs->size ==  param->num_query ){	
                qs = read_fasta_fastq2(fastafile,qs,param);
                fprintf(stderr,"%d read\n",qs->size);
                //qs = make_B(qs,qs->size);
                interval =  (int)((double)qs->size /(double)param->num_threads);
		
                for(t = 0;t < param->num_threads - 1;t++) {
                        qs = make_hash_thread(qs,param->seed_len,param->seed_step,t, t*interval,t*interval+interval);

                }
                qs = make_hash_thread(qs,param->seed_len,param->seed_step,param->num_threads - 1, t*interval,qs->size);
		
                //qs = search_fasta_fastq_thread(qs,param);
		
                for(t = 0;t < param->num_threads ;t++){
                        free(qs->hash_t[t]);
                        free(qs->hash_index_t[t]);
                }

        }


        fclose(fastafile);

        
        
        free_query_structure(qs);
        
        free_genome(genome);
        return OK;
        

	
        /*if(!isatty(0)){
          fastafile = stdin;
          }else{*/
        if (!(fastafile = fopen(param->input, "r" ))){
                fprintf(stderr,"Cannot open query file '%s'\n", param->input);
                exit(-1);
        }
        //}
	
        qs->hash_t = malloc(sizeof(int*) * param->num_threads);
        qs->hash_index_t = malloc(sizeof(int*) * param->num_threads);
	
        qs->size = param->num_query;
        //while(fpos != FILE_END){
        //fpos = read_fasta_fastq(fastafile,qs,param,fpos);
        while(qs->size ==  param->num_query ){	
                qs = read_fasta_fastq2(fastafile,qs,param);
                fprintf(stderr,"%d read\n",qs->size);
                //qs = make_B(qs,qs->size);
                interval =  (int)((double)qs->size /(double)param->num_threads);
		
                for(t = 0;t < param->num_threads - 1;t++) {
                        qs = make_hash_thread(qs,param->seed_len,param->seed_step,t, t*interval,t*interval+interval);

                }
                qs = make_hash_thread(qs,param->seed_len,param->seed_step,param->num_threads - 1, t*interval,qs->size);
		
                qs = search_fasta_fastq_thread(qs,param);
		
                for(t = 0;t < param->num_threads ;t++){
                        free(qs->hash_t[t]);
                        free(qs->hash_index_t[t]);
                }
		
        }
        rewind(fastafile);
	
        fclose(fastafile);
        free(qs->hash_t);
        free(qs->hash_index_t);
        free_query_structure(qs);
        return OK;
ERROR:
        return FAIL;
}


struct qs_struct* init_query_structure(struct parameters* param)
{
        struct qs_struct* qs = NULL;
        int i,j;

        ASSERT(param != NULL, "No parameters.");
        ASSERT(param->num_query > 0,"Numquery too small"); 
        
        MMALLOC(qs,sizeof(struct qs_struct));	
        qs->q_num = 0;
        qs->min_len = 100000000;
        qs->max_len = -1;
        qs->string_len = 0;
        qs->seq_info =NULL;
        qs->alloc_num_query = param->num_query;
        qs->hash_t = NULL;
        qs->hash_index_t = NULL; 

        MMALLOC(qs->hash_t,sizeof(int*) * param->num_threads);
        MMALLOC(qs->hash_index_t,sizeof(int*) * param->num_threads);

        MMALLOC(qs->seq_info,sizeof(struct seq_info*) *param->num_query);
        //qs->string = malloc(sizeof(unsigned char)*param->num_query*255);
        for(i = 0; i < param->num_query ;i++){
                qs->seq_info[i] = NULL;
                
                MMALLOC(qs->seq_info[i],sizeof(struct seq_info));
                qs->seq_info[i]->sn = 0;
                qs->seq_info[i]->reverse_seq = 0;
                qs->seq_info[i]->seq = 0;
                qs->seq_info[i]->qual = 0;
                qs->seq_info[i]->last_f_test = -1000;
                qs->seq_info[i]->last_r_test = -1000;
                qs->seq_info[i]->match_units = NULL;
                MMALLOC(qs->seq_info[i]->match_units,sizeof(unsigned long int) *  LIST_STORE_SIZE);//param->reported_seed_hits ); 
                for(j = 0; j < LIST_STORE_SIZE;j++){
                        qs->seq_info[i]->match_units[j] = 0ul;
                }
        }
        return qs;
ERROR:
        free_query_structure(qs);
        return NULL;
}

void free_query_structure(struct qs_struct * qs)
{
        int i;
        if(qs){
                for(i = 0; i < qs->alloc_num_query ;i++){
                        if(qs->seq_info[i]->sn){
                                MFREE(qs->seq_info[i]->sn);//
                                MFREE(qs->seq_info[i]->seq);
                                MFREE(qs->seq_info[i]->reverse_seq);
                        }
                        if(qs->seq_info[i]->qual){
                                MFREE(qs->seq_info[i]->qual);
                        }
                        MFREE(qs->seq_info[i]->match_units); 
                        MFREE(qs->seq_info[i]);
		
                }

                MFREE(qs->hash_t);
                MFREE(qs->hash_index_t);

                MFREE(qs->seq_info);  
                MFREE(qs);
        }
}



struct qs_struct* make_B(struct qs_struct* qs,int n)
{
        int i,j;
        struct seq_info* si = 0;
        for(i = 0; i < n;i++){
                si = qs->seq_info[i];
                for(j = 0; j < 4;j++){
                        //	si->rB[j]= 0ul;
                        //	si->fB[j]= 0ul;
                }
			
                for(j = 0; j < si->len;j++){
                        //	si->fB[(int)(si->seq[j] &0x3u)] |= (1ul << (unsigned long int)(si->len-1-j));
                        //	si->rB[(int)(si->reverse_seq[j] &0x3u)] |= (1ul << (unsigned long int)(si->len-1-j));
                }
        }
        return qs;
}


struct qs_struct* make_hash_thread(struct qs_struct* qs,int r,int step,int thread,int start,int stop)
{
        unsigned int key = 0;
        int i,j,c;
        int elements = 0;
        unsigned char* tmp = 0;
        int m;
        //int n = qs->size;
        fprintf(stderr,"thread:%d\n",thread);
        qs->hash_index_t[thread] = malloc(sizeof(int) * ((1 << r * 2) +1));
	
        for(i = 0; i < (1 << r * 2)+1;i++){
                qs->hash_index_t[thread][i] = 0;
        }
	
        for(i = start; i < stop;i++){
                m = qs->seq_info[i]->len;
                if(m > 64){
                        m = 64;
                }
                tmp = qs->seq_info[i]->seq;//   sequences[i];
                for(j = 0;j + r  <= m;j+= step){
                        key = 0;
                        for(c = 0;c < r;c++){
                                key = (key << 2u) | ((int)tmp[c]   & 0x3u);			
                        }
                        qs->hash_index_t[thread][key+1]++;
                        tmp += step;
                }	
                //m = qs->seq_info[i]->len;
                tmp = qs->seq_info[i]->reverse_seq;//   sequences[i];
                j = r;
		
                while(j+step <= m){
                        j += step;
                }
		
                tmp += m - j;
		
		
                //tmp += (m % r) % step;
                for(j = 0;j + r  <= m;j+=step){
                        key = 0;
                        for(c = 0;c < r;c++){
                                key = (key << 2u) | ((int)tmp[c]& 0x3u);		
                        }
                        qs->hash_index_t[thread][key+1]++;
                        tmp += step;
                }	
        }
	
        for(i = 1; i < (1 << r * 2)+1;i++){
                qs->hash_index_t[thread][i] +=  qs->hash_index_t[thread][i-1];
        }
        elements = qs->hash_index_t[thread][(1 << r * 2) ];
	
        qs->hash_t[thread]  = malloc(sizeof(int) * (elements));
        for(i = 0; i < elements;i++){
                qs->hash_t[thread][i] = -1;
        }
        for(i = start; i < stop;i++){
                m = qs->seq_info[i]->len;
                if(m > 64){
                        m = 64;
                }
                tmp = qs->seq_info[i]->seq;//   sequences[i];
                for(j = 0;j + r  <= m;j+= step){
                        key = 0;
                        for(c = 0;c < r;c++){
                                key = (key << 2u) | ((int)tmp[c]   & 0x3u);
				
                        }
			
                        for(c = qs->hash_index_t[thread][key];c < qs->hash_index_t[thread][key+1u];c++){
                                if(qs->hash_t[thread] [c] == -1){
                                        qs->hash_t[thread] [c] = (i << 8u) | j;
                                        //qs->hash[c] = i;
                                        break;
                                }
                        }
                        tmp += step;
                }	
                //m = qs->seq_info[i]->len;
                tmp = qs->seq_info[i]->reverse_seq;//   sequences[i];
		
                //for(j = 0; j < m;j++){
                //	fprintf(stderr,"%c",tmp[j]+65);
                //}
                //fprintf(stderr,"\n");
		
                j = r;
		
                while(j+step <= m){
                        //fprintf(stderr,"%d	%d\n",j,m);
                        j += step;
                }
                //elements = m - j;
                //fprintf(stderr,"offset:%d	%d	%d\n",elements,j,m);
		
                tmp += m - j;
		
                //tmp += (m % r) % step;
		
                for(j = 0;j + r  <= m;j+= step){
                        key = 0;
                        //for(c = 0;c < elements;c++){
                        //	fprintf(stderr,"-");
                        //}
                        //for(c = 0;c < j;c++){
                        //	fprintf(stderr,"-");
                        //}
                        for(c = 0;c < r;c++){
                                key = (key << 2u) | ((int)tmp[c]& 0x3u);
                                //	fprintf(stderr,"%c",tmp[c]+65);
                        }
                        //fprintf(stderr,"\n");
                        for(c = qs->hash_index_t[thread][key];c < qs->hash_index_t[thread][key+1u];c++){
                                if(qs->hash_t[thread][c] == -1){
                                        qs->hash_t[thread][c] = (i << 8u) | j | 0x80u;
                                        break;
                                }
                        }
                        tmp += step;
                }	
                //fprintf(stderr,"\n");
        }
        //exit(0);
        return qs;
}


void* search_hash_thread(void *threadarg)
{
	
        struct seed_thread_data *data;
        data = (struct seed_thread_data *) threadarg;
        struct qs_struct* qs = data->qs;
        struct parameters* param = data->param;
	
	
        unsigned char* t = data->chromsome;
        long int n = data->len;
        long int r = (long int) param->seed_len;
        long int offset = data->total_len;
	
        int thread = data->thread;
	
        int j;//,g;
        int i;
        unsigned char* target = 0;
        struct seq_info* si = 0;
        unsigned int key = 0;
        unsigned int mask = 0;
        int t_len = 0;
        long int seed_offset = 0;
        int strand = 0;
        unsigned long int new = 0;
        long int add = 0;
	
        int* index = qs->hash_index_t[thread];
	
        for(i = 0; i <  r;i++){
                mask = (mask << 2u) | 3u;
        }
        key = 0u;
        for(i = 0; i <  r-1;i++){
                key = (key << 2u) | (t[i] & 0x3u);
        }
        for(i = r-1; i < n;i++){
                key = ((key << 2u) | (t[i]  & 0x3u)) & mask;
                if(key){
                        //g = 0;
                        for(j = index[key];j < index[key+1u];j++){
                                seed_offset = qs->hash_t[thread][j] & 0x7fu;
				
                                si = qs->seq_info[ qs->hash_t[thread][j] >> 8u];
                                //id = qs->hash[j] >> 8u;
                                strand = qs->hash_t[thread][j] & 0x80u;
                                add = ( (long int)i + ((long int) si->len - seed_offset)  -r + 10l); /// should be 10 ??? 
                                //if( add > 0){
                                target = t +add;  
                                //}else{
                                //	target = t; /// WRONG????? 
                                //}
				
                                //if( add + si->len < n){
                                t_len = si->len + 10;
                                //if(add - t_len < 0){
                                //	t_len = add;
                                //}
                                //}else{
                                //	t_len = n - add;
                                //}
                                if( add <= n && add > ((long int) (si->len + 20))){
                                        if(!strand){
                                                if(si->last_f_test <= i){
                                                        //fprintf(stderr,"ADD:%ld	LEN:%d\n",add,t_len);
                                                        new = validate_bpm(target ,si->seq , t_len , si->len, offset+ add);
							
                                                        //new = validate_bpm_pre(target ,si->fB , t_len , si->len, offset+ add);
							
                                                        BinaryInsertionSortOne(si->match_units,LIST_STORE_SIZE ,new);
                                                        si->last_f_test = i +  si->len  ;
                                                }
                                        }else{
                                                if(si->last_r_test <= i){
                                                        new = validate_bpm(target ,si->reverse_seq , t_len , si->len,offset + add);
                                                        //new = validate_bpm_pre(target ,si->rB , t_len , si->len, offset+ add);
                                                        new |= 1ul ;
                                                        BinaryInsertionSortOne(si->match_units,LIST_STORE_SIZE,new);
                                                        si->last_r_test = i +   si->len ;
                                                }
                                        }
                                }
                        }
                }
                //if(i % 100000 == 0 && thread == 0){
                //	fprintf(stderr,"\r%8.2f percent done",(double)i /(double)n * 100);
                //}
        }
	
        //fprintf(stderr,"\r%8.2f percent done\n",100.0);
        pthread_exit((void *) 0);
}	




struct qs_struct* search_fasta_fastq_thread(struct qs_struct* qs,struct parameters* param) 
{
        static int header = 1;
        struct seed_thread_data* thread_data_array = 0; 
        thread_data_array = malloc(sizeof(struct seed_thread_data)* param->num_threads);
        pthread_t threads[param->num_threads];
        pthread_attr_t attr;
        int t;
        int interval = 0;
        int rc;
	
        FILE *file;
        unsigned int string_pos = 0;
        int park_pos = -1;
        char line[MAX_LINE];
        int i,j,c,n;
        int seq_p = 0;
        int set = 0;
        long int len = 0;
        long int total_len = 0;
        unsigned char* chromsome =  malloc(sizeof(unsigned char)* 300000000);
	
        long int* lengths = 0; 
        lengths = malloc(sizeof(long int) * 500);
        unsigned char** chr_names = 0;
        chr_names = malloc(sizeof(unsigned char*) * 500);
	
	
        if (!(file = fopen(param->genome, "r" ))){
                fprintf(stderr,"Cannot open query file '%s'\n", param->genome);
                exit(-1);
        }
        c= -1;
	
        while(fgets(line, MAX_LINE, file)){
                //reading name .... 
                if((line[0] == '@' && !set)|| (line[0] == '>' && !set)){
                        //set sequence length of previous read
			
                        if(park_pos != -1){ // to avoid setting something before the first read was read
				
				
                                interval =  (int)((double)qs->size /(double)param->num_threads);
				
                                for(t = 0;t < param->num_threads - 1;t++) {
                                        thread_data_array[t].qs = qs;
                                        thread_data_array[t].param = param;
                                        thread_data_array[t].start = t*interval;
                                        thread_data_array[t].end = t*interval + interval;
                                        thread_data_array[t].thread = t;
                                        thread_data_array[t].chromsome = chromsome;
                                        thread_data_array[t].total_len = total_len;
                                        thread_data_array[t].len = len;
                                }
                                thread_data_array[param->num_threads - 1].qs = qs;
                                thread_data_array[param->num_threads - 1].param = param;
                                thread_data_array[param->num_threads - 1].start = t*interval;
                                thread_data_array[param->num_threads - 1].end = qs->size;
                                thread_data_array[param->num_threads - 1].thread = param->num_threads - 1;
                                thread_data_array[param->num_threads - 1].chromsome = chromsome;
                                thread_data_array[param->num_threads - 1].total_len = total_len;
                                thread_data_array[param->num_threads - 1].len = len;
				
				
				
				
                                rc = pthread_attr_init(&attr);
                                if(rc){
                                        fprintf(stderr,"ERROR; return code from pthread_attr_init() is %d\n", rc);
                                        exit(-1);
                                }
                                pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
				
                                for(t = 0;t < param->num_threads;t++) {
                                        rc = pthread_create(&threads[t], &attr, search_hash_thread, (void *) &thread_data_array[t]);
                                        if (rc) {
                                                fprintf(stderr,"ERROR; return code from pthread_create() is %d\n", rc);
                                                exit(-1);
                                        }
                                }
				
                                pthread_attr_destroy(&attr);
				
                                for (t = 0;t < param->num_threads;t++){
                                        rc = pthread_join(threads[t], NULL);
                                        if (rc){
                                                fprintf(stderr,"ERROR; return code from pthread_join()is %d\n", rc);
                                                exit(-1);
                                        }
                                }
				
				
				
				
                                //qs =  search_hash(qs,chromsome,len,param->r , total_len);
				
                                //exit(0);
				
                                for(i = 0; i< qs->size;i++){
                                        qs->seq_info[i]->last_f_test = -1000;
                                        qs->seq_info[i]->last_r_test = -1000;
                                }
                                total_len += len;
                                //fprintf(stderr,"	%d	%ud	%ud	(%f)\n",len,mers,total_len,(float)mers/(float)total_len*100.0f);
                        }
                        fprintf(stderr,"%s	%ld\n",line,total_len);
                        //fprintf(stderr,"%s",line);
			
                        n = strlen(line);
                        c++;
			
                        chr_names[c] = malloc(sizeof(unsigned char) * (n+1));
                        for(i = 1; i < n-1;i++){
                                chr_names[c][i-1] = line[i];
                        }
                        chr_names[c][n-2] = 0; 
                        lengths[c] = len;
			
                        string_pos = 0;			
                        len = 0;
                        seq_p = 1;
                        set = 1;			
                        park_pos++;
                        //get ready to read quality if present  
                }else if(line[0] == '+' && !set){
                        seq_p = 0;
                        set = 1;
                        //reading sequence or quality  
                }else{			
                        if(!set){
                                string_pos--;
                        }
			
                        if(seq_p){
                                for(i = 0;i < MAX_LINE;i++){
                                        if(iscntrl((int)line[i])){
                                                chromsome[string_pos] = 0;
                                                //qt->string[string_pos] = 0;
                                                string_pos++;
                                                break;
                                        }
					
                                        len++;
                                        chromsome[string_pos] = nuc_code5[(int)line[i]];
                                        string_pos++;
                                }
                        }
                        set = 0;
			
                }
        }
		
        interval =  (int)((double)qs->size /(double)param->num_threads);
	
        for(t = 0;t < param->num_threads - 1;t++) {
                thread_data_array[t].qs = qs;
                thread_data_array[t].param = param;
                thread_data_array[t].start = t*interval;
                thread_data_array[t].end = t*interval + interval;
                thread_data_array[t].thread = t;
                thread_data_array[t].chromsome = chromsome;
                thread_data_array[t].total_len = total_len;
                thread_data_array[t].len = len;
        }
        thread_data_array[param->num_threads - 1].qs = qs;
        thread_data_array[param->num_threads - 1].param = param;
        thread_data_array[param->num_threads - 1].start = t*interval;
        thread_data_array[param->num_threads - 1].end = qs->size;
        thread_data_array[param->num_threads - 1].thread = param->num_threads - 1;
        thread_data_array[param->num_threads - 1].chromsome = chromsome;
        thread_data_array[param->num_threads - 1].total_len = total_len;
        thread_data_array[param->num_threads - 1].len = len;
	
	
	
	
        rc = pthread_attr_init(&attr);
        if(rc){
                fprintf(stderr,"ERROR; return code from pthread_attr_init() is %d\n", rc);
                exit(-1);
        }
        pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	
        for(t = 0;t < param->num_threads;t++) {
                rc = pthread_create(&threads[t], &attr, search_hash_thread, (void *) &thread_data_array[t]);
                if (rc) {
                        fprintf(stderr,"ERROR; return code from pthread_create() is %d\n", rc);
                        exit(-1);
                }
        }
	
        pthread_attr_destroy(&attr);
	
        for (t = 0;t < param->num_threads;t++){
                rc = pthread_join(threads[t], NULL);
                if (rc){
                        fprintf(stderr,"ERROR; return code from pthread_join()is %d\n", rc);
                        exit(-1);
                }
        }
        c++;
	
        lengths[c] = len;
        chr_names[c] = 0;
        if(header){
                fprintf(stdout, "@HD\tVN:1.0\n");
		
                for(i = 0; i < c ;i++){
                        fprintf(stdout, "@SQ\tSN:%s\tLN:%u",chr_names[i], (unsigned int) lengths[i+1]);
                        
                        fprintf(stdout,"\n");
                        
                }
                fprintf(stdout, "@PG\tID:Delve\tVN:%s\tCL:%s\n",VERSION,"dseed something") ;
                header = 0;
        }
	
	
	
        for(i = 1; i <= c;i++){
                lengths[i] +=  lengths[i-1];
        }
	
	
        //total_len +=  len;
        for(i = 0; i < qs->size;i++){
                print_sam(qs->seq_info[i], chr_names, lengths ,c);
        }
	
        //exit(0);
        for(i = 0; i< qs->size;i++){
                for(j = 0; j < LIST_STORE_SIZE;j++){
                        qs->seq_info[i]->match_units[j] = 0ul;
                }
        }
        //string_pos = 0;
        qs->q_num = 0;
	
        for(i = 0; i < c;i++){
                free(chr_names[i]);
        }
	
        free(thread_data_array);
        free(chr_names);
        free(lengths);
        free(chromsome);
        fclose(file);
        return qs;
}





int seed_controller(struct parameters* param)
{
        struct qs_struct* qs = 0;
        unsigned long int fpos = 0u;
        //int c;
        FILE *fastafile = 0; 
        qs = init_query_structure(param);
	
	
        if (!(fastafile = fopen(param->input, "r" ))){
                fprintf(stderr,"Cannot open query file '%s'\n", param->input);
                exit(-1);
        }
        //c = 0;
        while(fpos != FILE_END){
                fpos = read_fasta_fastq(fastafile,qs,param,fpos);
                /*fprintf(stderr,"Chumk:%d\n",c);
                  for(i = 0; i < qs->size;i++){
                  fprintf(stderr,"%s	%s	",qs->seq_info[i]->sn,qs->seq_info[i]->qual);
                  for(j = 0;j < qs->seq_info[i]->len;j++){
                  fprintf(stderr,"%c",qs->seq_info[i]->seq[j] + 65);
                  }
                  fprintf(stderr,"\t");
                  for(j = 0;j < qs->seq_info[i]->len;j++){
                  fprintf(stderr,"%c",qs->seq_info[i]->reverse_seq[j] + 65);
                  }
                  fprintf(stderr,"\n");
		 
                  }*/
                //c++;
		
                qs = make_hash(qs,param->seed_len);
                qs = search_fasta_fastq(qs,param);
                free(qs->hash);
                free(qs->hash_index);
		
        }
        rewind(fastafile);
	
        fclose(fastafile);
        return 1;
}



struct qs_struct* read_fasta_fastq2(FILE *file, struct qs_struct* qs,struct parameters* param) 
{
        int park_pos = -1;
        char line[MAX_LINE];
        unsigned char rc[MAX_LINE];
        int i,j;
        int seq_p = 0;
        int set = 0;
        int len = 0;
        int fastq = 0;
        qs->size = 0;
        for(i = 0; i < param->num_query ;i++){
                //qs->seq_info[i] = malloc(sizeof(struct seq_info));
                if(qs->seq_info[i]->sn){
                        free(qs->seq_info[i]->sn);//
                        free(qs->seq_info[i]->seq);
                        free(qs->seq_info[i]->reverse_seq);
                }
                if(qs->seq_info[i]->qual){
                        //fprintf(stderr,"%s\n",qs->seq_info[i]->qual);
                        free(qs->seq_info[i]->qual);
                }

                qs->seq_info[i]->sn = 0;
                qs->seq_info[i]->seq = 0;
                qs->seq_info[i]->reverse_seq = 0;
                qs->seq_info[i]->qual = 0;
		
		
                qs->seq_info[i]->last_f_test = -1000;
                qs->seq_info[i]->last_r_test = -1000;
                //qs->seq_info[i]->match_units= malloc(sizeof(unsigned long int) *  LIST_STORE_SIZE);//param->reported_seed_hits ); 
                for(j = 0; j < LIST_STORE_SIZE;j++){
                        qs->seq_info[i]->match_units[j] = 0ul;
                }
        }
        while(fgets(line, MAX_LINE, file)){
                //reading name .... 
                //fprintf(stderr,"%s",line);
                //if(line[0] != '#'){ 
                if((line[0] == '@' && !set)|| (line[0] == '>' && !set)){
                        //set sequence length of previous read
                        //if(park_pos != -1){
                        //	fprintf(stderr,"%s	%s	%s\n",ri[park_pos]->name,ri[park_pos]->seq);
                        //}
				
                        //check if there is still space.... 
                        //if(param->num_query == size){
                        //	fseek (file , -  strlen(line) , SEEK_CUR);
				
                        //	return size;
                        //}
                        park_pos++;
				
                        qs->seq_info[park_pos]->last_f_test = -1000;
                        qs->seq_info[park_pos]->last_r_test = -1000;
                        qs->seq_info[park_pos]->sn = 0;
                        qs->seq_info[park_pos]->seq = 0;
                        qs->seq_info[park_pos]->reverse_seq = 0;
                        qs->seq_info[park_pos]->qual = 0;
                        len = 0;
                        seq_p = 1;
                        for(i = 1;i < MAX_LINE;i++){
                                len++;
                                if(isspace((int)line[i])){
                                        break;
                                }
					
                        }
                        qs->seq_info[park_pos]->sn = malloc(sizeof(unsigned char)* (len+1));
                        for(i = 1;i < MAX_LINE;i++){
					
                                if(isspace((int)line[i])){
                                        qs->seq_info[park_pos]->sn[i-1] = 0;
                                        break;
                                }
                                qs->seq_info[park_pos]->sn[i-1] = line[i];
                        }
                        //fprintf(stderr,"LEN:%d	%s\n",len,ri[park_pos]->name);
				
                        set = 1;
                        qs->size++;
                        //get ready to read quality if present  
                }else if(line[0] == '+' && !set){
                        seq_p = 0;
                        set = 1;
                        fastq = 1;
                        //fprintf(stderr,"fastq	%d\n",fastq);
                        //reading sequence or quality  
                }else{	
                        if(set){
                                if(seq_p){
                                        len = 0;
                                        for(i = 0;i < MAX_LINE;i++){
							
                                                if(isspace((int)line[i])){
                                                        break;
                                                }
                                                len++;
							
                                        }
                                        //fprintf(stderr,"SEQ LEN:%d	%s\n",len,line);
                                        qs->seq_info[park_pos]->seq = malloc(sizeof(unsigned char)* (len+1));
                                        for(i = 0;i < MAX_LINE;i++){
							
                                                if(isspace((int)line[i])){
                                                        qs->seq_info[park_pos]->seq[i] = 0;
                                                        break;
                                                }
                                                qs->seq_info[park_pos]->seq[i] = nuc_code5[(int)line[i]];
                                        }
                                        qs->seq_info[park_pos]->len = len;
						
						
                                        reverse_complement(qs->seq_info[park_pos]->seq ,rc,len);
                                        qs->seq_info[park_pos]->reverse_seq = malloc(sizeof(unsigned char)* (len+1));
                                        for(i = 0; i < len;i++){
                                                qs->seq_info[park_pos]->reverse_seq[i] = rc[i];
                                        }
                                        qs->seq_info[park_pos]->reverse_seq[len] = 0;
						
						
                                        if(param->num_query == qs->size && fastq == 0){
                                                return qs;
                                        }
                                }else{
                                        len = 0;
                                        for(i = 0;i < MAX_LINE;i++){
                                                len++;
                                                if(isspace((int)line[i])){
                                                        break;
                                                }
							
                                        }
                                        //fprintf(stderr,"QUAL LEN:%d\n",len);
                                        qs->seq_info[park_pos]->qual = malloc(sizeof(unsigned char)* (len+1));
                                        for(i = 0;i < MAX_LINE;i++){
							
                                                if(isspace((int)line[i])){
                                                        qs->seq_info[park_pos]->qual[i] = 0;
                                                        break;
                                                }
                                                qs->seq_info[park_pos]->qual[i] = line[i];
                                        }
						
                                        if(param->num_query == qs->size){
                                                return qs;
                                        }
                                }
                        }
                        set = 0;
                }
                //}
        }
        return qs;
}

unsigned long int read_fasta_fastq(FILE *file, struct qs_struct* qs,struct parameters* param,unsigned long int position) 
{
	
        unsigned int string_pos = 0;
        int park_pos = -1;
        char line[MAX_LINE];
        unsigned char rc[MAX_LINE];
        int i;//,j;
        int seq_p = 0;
        int set = 0;
        int len = 0;
        int size = 0;
        unsigned long int fpos = 0;
	
        fseek (file , position , SEEK_SET);
        //size = ftell (file);
        //fprintf(stderr,"Position = %d\n",size);
	
        while(fgets(line, MAX_LINE, file)){
                //reading name .... 
                //fprintf(stderr,"%s",line);
		
                if(line[0] != '#'){
		
                        if((line[0] == '@' && !set)|| (line[0] == '>' && !set)){
                                //set sequence length of previous read
			
                                if(park_pos != -1){ // to avoid setting something before the first read was read
                                        //len = (int)strlen((char*)qs->seq_info[park_pos]->seq);
                                        qs->seq_info[park_pos]->len = len;
				
                                        if(len < qs->min_len){
                                                qs->min_len = len;
                                        }
				
                                        if(len > qs->max_len){
                                                qs->max_len = len;
                                        }			
				
                                        reverse_complement(qs->seq_info[park_pos]->seq ,rc,len);
                                        qs->seq_info[park_pos]->reverse_seq = qs->string + string_pos;
                                        for(i = 0; i < len;i++){
                                                qs->string[string_pos] = rc[i];
                                                string_pos++;
                                        }
                                        qs->string[string_pos] = 0;
                                        string_pos++;
				
				
                                }
			
                                //check if there is still space.... 
                                if(string_pos + MAX_LINE >= (param->num_query *255) || size >= param->num_query){
                                        //fprintf(stderr,"\nLast_line = %s",line);
                                        //size = ftell (file);
                                        //fprintf(stderr,"Position = %d\n",size);
                                        fseek (file , -  strlen(line) , SEEK_CUR);
                                        fpos = ftell (file);
                                        qs->size = size;
                                        //fgets(line, MAX_LINE, file);
                                        //fprintf(stderr,"Last_line = %s",line);
                                        //fprintf(stderr,"Position = %d\n",fpos);
                                        return fpos;
                                        for(i= 0; i < park_pos;i++){
                                                fprintf(stderr,"%d	%s	%s	%s\n",qs->seq_info[i]->id_num,qs->seq_info[i]->sn,qs->seq_info[i]->seq,qs->seq_info[i]->qual);
                                        }
                                        size = 0;
                                        string_pos = 0;
                                        park_pos = -1;
                                }
                                len = 0;
                                seq_p = 1;
                                park_pos++;
                                qs->q_num++;
                                qs->seq_info[park_pos]->id_num = qs->q_num;
                                qs->seq_info[park_pos]->sn =  qs->string + string_pos;
                                qs->seq_info[park_pos]->qual = 0;
                                qs->seq_info[park_pos]->seq = 0;
                                qs->seq_info[park_pos]->last_f_test = -1000;
                                qs->seq_info[park_pos]->last_r_test = -1000;
                                for(i = 1;i < MAX_LINE;i++){
                                        if(iscntrl((int)line[i])){
                                                qs->string[string_pos] =0;
                                                string_pos++;
                                                break;
                                        }
                                        if(isspace((int)line[i])){
                                                qs->string[string_pos] = '_';
                                        }else{
                                                qs->string[string_pos] = line[i];
                                        }
                                        string_pos++;
                                }
                                set = 1;
                                size++;
                                //get ready to read quality if present  
                        }else if(line[0] == '+' && !set){
                                seq_p = 0;
                                set = 1;
                                //reading sequence or quality  
                        }else{	
                                if(set){
                                        if(seq_p){
                                                qs->seq_info[park_pos]->seq = qs->string + string_pos;
                                        }else{
                                                qs->seq_info[park_pos]->qual = qs->string + string_pos;
                                        }
                                }else{
                                        string_pos--;
                                }
			
                                for(i = 0;i < MAX_LINE;i++){
                                        if(iscntrl((int)line[i])){
                                                qs->string[string_pos] = 0;
                                                string_pos++;
                                                break;
                                        }
                                        if(seq_p){
                                                //fprintf(stderr,"Reading seq:%c	%d\n",line[i],nuc_code5[line[i]]);
                                                len++;
                                                qs->string[string_pos] = nuc_code5[(int)line[i]];
                                        }else{
                                                qs->string[string_pos] = line[i];
                                        }
                                        string_pos++;
                                }
                                set = 0;
                        }
                }
        }
        if(park_pos != -1){ // to avoid setting something before the first read was read
                //len = (int)strlen((char*)qs->seq_info[park_pos]->seq);
                qs->seq_info[park_pos]->len = len;
                //qs->seq_info[park_pos]->id_num = qs->q_num;
                if(len < qs->min_len){
                        qs->min_len = len;
                }
		
                if(len > qs->max_len){
                        qs->max_len = len;
                }		
                reverse_complement(qs->seq_info[park_pos]->seq ,rc,len);
                qs->seq_info[park_pos]->reverse_seq = qs->string + string_pos;
                for(i = 0; i < len;i++){
                        qs->string[string_pos] = rc[i];
                        string_pos++;
                }
                qs->string[string_pos] = 0;
                //string_pos++;
        }
        qs->size = size;
        return FILE_END;

}

struct qs_struct* search_fasta_fastq(struct qs_struct* qs,struct parameters* param) 
{
        FILE *file;
        unsigned int string_pos = 0;
        int park_pos = -1;
        char line[MAX_LINE];
        int i,j,c,n;
        int seq_p = 0;
        int set = 0;
        long int len = 0;
        long int total_len = 0;
        unsigned char* chromsome =  malloc(sizeof(unsigned char)* 300000000);
	
        long int* lengths = 0; 
        lengths = malloc(sizeof(unsigned int) * 500);
        unsigned char** chr_names = 0;
        chr_names = malloc(sizeof(unsigned char*) * 500);
	
	
        if (!(file = fopen(param->genome, "r" ))){
                fprintf(stderr,"Cannot open query file '%s'\n", param->genome);
                exit(-1);
        }
        c= -1;
	
        while(fgets(line, MAX_LINE, file)){
                //reading name .... 
                if((line[0] == '@' && !set)|| (line[0] == '>' && !set)){
                        //set sequence length of previous read
			
                        if(park_pos != -1){ // to avoid setting something before the first read was read
				
                                qs =  search_hash(qs,chromsome,len,param->seed_len , total_len);
                                for(i = 0; i< qs->size;i++){
                                        qs->seq_info[i]->last_f_test = -1000;
                                        qs->seq_info[i]->last_r_test = -1000;
                                }
                                total_len += len;
                                //fprintf(stderr,"	%d	%ud	%ud	(%f)\n",len,mers,total_len,(float)mers/(float)total_len*100.0f);
                        }
                        //fprintf(stderr,"%s	%ld",line,total_len);
			
                        n = strlen(line);
                        c++;
			
                        chr_names[c] = malloc(sizeof(unsigned char) * (n+1));
                        for(i = 1; i < n-1;i++){
                                chr_names[c][i-1] = line[i];
                        }
                        chr_names[c][n-1] = 0; 
                        lengths[c] = len;
			
                        string_pos = 0;			
                        len = 0;
                        seq_p = 1;
                        set = 1;			
                        park_pos++;
                        //get ready to read quality if present  
                }else if(line[0] == '+' && !set){
                        seq_p = 0;
                        set = 1;
                        //reading sequence or quality  
                }else{			
                        if(!set){
                                string_pos--;
                        }
			
                        if(seq_p){
                                for(i = 0;i < MAX_LINE;i++){
                                        if(iscntrl((int)line[i])){
                                                chromsome[string_pos] = 0;
                                                //qt->string[string_pos] = 0;
                                                string_pos++;
                                                break;
                                        }
					
                                        len++;
                                        chromsome[string_pos] = nuc_code5[(int)line[i]];
                                        string_pos++;
                                }
                        }
                        set = 0;
			
                }
        }
        c++;

        lengths[c] = len;
        chr_names[c] = 0;
	
        for(i = 1; i <= c;i++){
                lengths[i] +=  lengths[i-1];
        }
        qs =  search_hash(qs,chromsome,len,param->seed_len , total_len);

        //for(i = 0; i <= c;i++){
        ///	fprintf(stderr,"%s	%ld\n",chr_names[i],lengths[i]);
        //}
        //qs = search_chromosome_hash(qs,param,chromsome,len,total_len);
        //total_len +=  len;
	
        //fwrite(&param->upper_limit, sizeof(unsigned int),1,binfile);
        for(i = 0; i< qs->size;i++){
                print_sam(qs->seq_info[i], chr_names, lengths ,c);
        }
	
        for(i = 0; i< qs->size;i++){
                for(j = 0; j < LIST_STORE_SIZE;j++){
                        qs->seq_info[i]->match_units[j] = 0ul;
                }
        }
        //string_pos = 0;
        qs->q_num = 0;
	
        for(i = 0; i < c;i++){
                free(chr_names[i]);
        }
	
        free(chr_names);
        free(lengths);
        free(chromsome);
        fclose(file);
        return qs;
}

void print_sam(struct seq_info* si, unsigned char** t_names, long int* lengths,int size)
{
        char alphabet[] = "ACGTN";
        char* pline = malloc(sizeof (char)* MAX_LINE*10);
        char seq[MAX_LINE];
        unsigned char star[] = "*";
        int strand = 0;
        int no_qual = 0;
        int i,j;
        unsigned char* t_name =  0;
        long int pos = 0;
        int matches = 0;
	
        pline[0] = 0 ;
        strand = si->match_units[0] & 1ul;
        if(strand){
                strand = 16;
        }
	
        //fprintf(stderr,"before : %s\n",si->qual);
	
        if(!si->qual){
                si->qual = star;
                no_qual = 1;
        }
	
        //fprintf(stderr,"after : %s\n",si->qual);
	
        if(!si->match_units[0]){
		
		
		
                for(i = 0;i < si->len;i++){
                        seq[i] = alphabet[si->seq[i]];
                }
                seq[si->len] = 0;
		
                sprintf(pline,"%s\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%s",
                        si->sn ,//1 QNAME Query NAME of the read or the read pair 
                        4,//2 FLAG bitwise FLAG (pairing, strand, mate strand, etc.) 
                        "*",//3 RNAME Reference sequence NAME 
                        0,//4 POS 1-based leftmost POSition of clipped alignment 
                        0,//5 MAPQ MAPping Quality (Phred-scaled) 
                        "*",//6 CIGAR extended CIGAR string (operations: MIDNSHP) 
                        "*",//7 MRNM Mate Reference NaMe (‘=’ if same as RNAME) 
                        0,//8 MPOS 1-based leftmost Mate POSition 
                        0,//9 ISIZE inferred Insert SIZE 
                        seq,//10 SEQ query SEQuence on the same strand as the reference 
                        si->qual//11 QUAL query QUALity (ASCII-33=Phred base quality) 
                        );
                i = strlen(pline);
                pline[i] = 0;
                fprintf(stdout, "%s\n",pline);
        }else{
                pos =(long int) ((si->match_units[0] >> 1ul )& 0xFFFFFFFFul);
                matches = (int) (si->match_units[0] >> 56ul);
                //fprintf(stderr,"printing:%lu	%ld\n",si->match_units[0],pos);
                j = 0;
                for(i = 0; i <= size;i++){
                        if( pos < lengths[i]){
                                t_name = t_names[i-1];
                                pos = pos - lengths[i-1];
                                j = i-1 + 1;
                                break;
                        }
                }
                if(!j){
                        fprintf(stderr,"not found: %ld\n",pos);
                        for(i = 0; i <= size;i++){
                                fprintf(stderr,"%s	%ld	%s\n", t_names[i],lengths[i],si->sn);
                        }
                        exit(-1);
                }
                //fprintf(stderr,"NAME found : %d\n",j-1);
		
                if(!strand){
                        for(i = 0;i < si->len;i++){
                                seq[i] = alphabet[si->seq[i]];
                        }
                }else{
                        for(i = 0;i < si->len;i++){
                                seq[i] = alphabet[si->reverse_seq[i]];
                        }
                }
		
                seq[si->len] = 0;
                //fprintf(stdout, "LINe:\n\n");
                if(!strand){
                        sprintf(pline,"%s\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%s",
                                si->sn ,//1 QNAME Query NAME of the read or the read pair 
                                strand,//2 FLAG bitwise FLAG (pairing, strand, mate strand, etc.) 
                                t_name,//3 RNAME Reference sequence NAME 
                                (int)(pos+1l),//4 POS 1-based leftmost POSition of clipped alignment 
                                255,//5 MAPQ MAPping Quality (Phred-scaled) 
                                "*",//si->len,//6 CIGAR extended CIGAR string (operations: MIDNSHP) 
                                "*",//7 MRNM Mate Reference NaMe (‘=’ if same as RNAME) 
                                0,//8 MPOS 1-based leftmost Mate POSition 
                                0,//9 ISIZE inferred Insert SIZE 
                                seq,//10 SEQ query SEQuence on the same strand as the reference 
                                si->qual//11 QUAL query QUALity (ASCII-33=Phred base quality) 
                                );
                }else{
                        if(no_qual){
                                sprintf(pline,"%s\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%s",
                                        si->sn ,//1 QNAME Query NAME of the read or the read pair 
                                        strand,//2 FLAG bitwise FLAG (pairing, strand, mate strand, etc.) 
                                        t_name,//3 RNAME Reference sequence NAME 
                                        (int)(pos),//4 POS 1-based leftmost POSition of clipped alignment 
                                        255,//5 MAPQ MAPping Quality (Phred-scaled) 
                                        "*",//si->len,//6 CIGAR extended CIGAR string (operations: MIDNSHP) 
                                        "*",//7 MRNM Mate Reference NaMe (‘=’ if same as RNAME) 
                                        0,//8 MPOS 1-based leftmost Mate POSition 
                                        0,//9 ISIZE inferred Insert SIZE 
                                        seq,//10 SEQ query SEQuence on the same strand as the reference 
                                        si->qual//11 QUAL query QUALity (ASCII-33=Phred base quality) 
                                        );
                        }else{
                                sprintf(pline,"%s\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%s",
                                        si->sn ,//1 QNAME Query NAME of the read or the read pair 
                                        strand,//2 FLAG bitwise FLAG (pairing, strand, mate strand, etc.) 
                                        t_name,//3 RNAME Reference sequence NAME 
                                        (int)(pos),//4 POS 1-based leftmost POSition of clipped alignment 
                                        255,//5 MAPQ MAPping Quality (Phred-scaled) 
                                        "*",//si->len,//6 CIGAR extended CIGAR string (operations: MIDNSHP) 
                                        "*",//7 MRNM Mate Reference NaMe (‘=’ if same as RNAME) 
                                        0,//8 MPOS 1-based leftmost Mate POSition 
                                        0,//9 ISIZE inferred Insert SIZE 
                                        seq,//10 SEQ query SEQuence on the same strand as the reference 
                                        reverse_without_complement(si->qual, si->len)//11 QUAL query QUALity (ASCII-33=Phred base quality) 
                                        );
                        }
                }
		
                i = strlen(pline);
                pline[i] = '\t';
                i++;
                sprintf(pline+i,"NM:i:%d",  si->len - matches);
		
		
                if(si->match_units[1]){
			
                        i = strlen(pline);
                        pline[i] = '\t';
                        i++;
                        sprintf(pline+i,"XA:Z:");
                        j = 1;
                        while( j < LIST_STORE_SIZE){
                                if(! si->match_units[j] ){
                                        break;
                                }
                                strand = si->match_units[j] & 1ul;
                                pos = (unsigned int)  ((si->match_units[j] >> 1ul )& 0xFFFFFFFFul);
                                for(i = 0; i <= size;i++){
                                        //fprintf(stderr," looking at: %d	%d	%d	%s\n",i,pos,lengths[i],t_names[i]);
                                        if( pos < lengths[i]){
                                                t_name = t_names[i-1];
                                                pos = pos - lengths[i-1];
                                                break;
                                        }
                                }
				
                                i = strlen(pline);
                                if(!strand){
                                        sprintf(pline+i,"%s,+%d,%dM,%d;", t_name,(int)(pos+1l),si->len, si->len - (int)(si->match_units[j] >> 56ul));
                                }else{
                                        sprintf(pline+i,"%s,-%d,%dM,%d;", t_name,(int)(pos+1l),si->len, si->len - (int)(si->match_units[j] >> 56ul));
                                }
                                j++;
                        }
			
                }
		
                i = strlen(pline);
                pline[i] = 0;
                fprintf(stdout, "%s\n",pline);
        }
        if(no_qual){
                si->qual = 0;
        }
	
        free(pline);
}


void reverse_complement(unsigned char* p,unsigned char* rc,int len)
{
        int rev[5] = {3,2,1,0,4};
	
        int i,c;
        c = 0;
        for(i =len-1; i >= 0;i--){
                rc[c] = rev[(int)p[i]];
                c++;
        }
        rc[c] = 0;
}

struct qs_struct* make_hash(struct qs_struct* qs,int r)
{
        unsigned int key = 0;
        int i,j,c;
        int elements = 0;
        unsigned char* tmp = 0;
        int m;
        int n = qs->size;
	
        qs->hash_index = malloc(sizeof(int) * ((1 << r * 2) +1));
	
        for(i = 0; i < (1 <<r * 2)+1;i++){
                qs->hash_index[i] = 0;
        }
	
        for(i = 0; i < n;i++){
                m = qs->seq_info[i]->len;
                if(m > 31){
                        m = 31;
                }
                tmp = qs->seq_info[i]->seq;//   sequences[i];
                for(j = 0;j + r  <= m;j+= r){
                        key = 0;
                        for(c = 0;c < r;c++){
                                key = (key << 2u) | ((int)tmp[c]   & 0x3u);			
                        }
                        qs->hash_index[key+1]++;
                        tmp += r;
                }	
                //m = qs->seq_info[i]->len;
                tmp = qs->seq_info[i]->reverse_seq;//   sequences[i];
                for(j = 0;j + r  <= m;j+= r){
                        key = 0;
                        for(c = 0;c < r;c++){
                                key = (key << 2u) | ((int)tmp[c]& 0x3u);		
                        }
                        qs->hash_index[key+1]++;
                        tmp += r;
                }	
        }
	
        for(i = 1; i < (1 << r * 2)+1;i++){
                qs->hash_index[i] +=  qs->hash_index[i-1];
        }
        elements = qs->hash_index[(1 << r * 2) ];
	
        qs->hash = malloc(sizeof(int) * (elements));
        for(i = 0; i < elements;i++){
                qs->hash[i] = -1;
        }
        for(i = 0; i < n;i++){
                m = qs->seq_info[i]->len;
                if(m > 31){
                        m = 31;
                }
                tmp = qs->seq_info[i]->seq;//   sequences[i];
                for(j = 0;j + r  <= m;j+= r){
                        key = 0;
                        for(c = 0;c < r;c++){
                                key = (key << 2u) | ((int)tmp[c]   & 0x3u);			
                        }
                        for(c = qs->hash_index[key];c < qs->hash_index[key+1u];c++){
                                if(qs->hash[c] == -1){
                                        qs->hash[c] = (i << 8u) | j;
                                        //qs->hash[c] = i;
                                        break;
                                }
                        }
                        tmp += r;
                }	
		
		
		
                //m = qs->seq_info[i]->len;
                tmp = qs->seq_info[i]->reverse_seq;//   sequences[i];
		
		
		
                for(j = 0;j + r  <= m;j+= r){
                        key = 0;
                        for(c = 0;c < r;c++){
                                key = (key << 2u) | ((int)tmp[c]& 0x3u);
			
                        }
                        for(c = qs->hash_index[key];c < qs->hash_index[key+1u];c++){
                                if(qs->hash[c] == -1){
                                        qs->hash[c] = (i << 8u) | j | 0x80u;
                                        break;
                                }
                        }
                        tmp += r;
                }	
        }
        //exit(0);
        return qs;
}

struct qs_struct* search_hash(struct qs_struct* qs,unsigned char* t,int n, long int r,long int offset)
{
        int j;
        int i;
        unsigned char* target = 0;
        struct seq_info* si = 0;
        unsigned int key = 0;
        unsigned int mask = 0;
        int t_len = 0;
        long int seed_offset = 0;
        int strand = 0;
        unsigned long int new = 0;
        //static int test = 0;
        long int add = 0;

        for(i = 0; i <  r;i++){
                mask = (mask << 2u) | 3u;
        }
        //timer(0,0);
        key = 0u;
        for(i = 0; i <  r-1;i++){
                key = (key << 2u) | (t[i] & 0x3u);
        }
        for(i = r-1; i < n;i++){
                key = ((key << 2u) | (t[i] & 0x3u)) & mask;
                if(key){
                        for(j = qs->hash_index[key];j < qs->hash_index[key+1u];j++){
                                seed_offset = qs->hash[j] & 0x7fu;
				
                                si = qs->seq_info[ qs->hash[j] >> 8u];
                                //id = qs->hash[j] >> 8u;
                                strand = qs->hash[j] & 0x80u;
                                add = ( (long int)i + ((long int) si->len - seed_offset)  -r + 5l);
                                if( add > 0){
                                        target = t +add;  
                                }else{
                                        target = t;
                                        add = 0;
                                }
				
				
                                t_len = si->len + 10;
                                if( add <= n){
                                        if(!strand){
                                                if(si->last_f_test <= i){
                                                        new = validate_bpm(target ,si->seq , t_len , si->len, offset+ add);
                                                        BinaryInsertionSortOne(si->match_units,LIST_STORE_SIZE ,new);
                                                        si->last_f_test = i + si->len  ;
                                                }
                                        }else{
                                                if(si->last_r_test <= i){
                                                        new = validate_bpm(target ,si->reverse_seq , t_len , si->len,offset + add);
                                                        new |= 1ul ;
                                                        BinaryInsertionSortOne(si->match_units,LIST_STORE_SIZE,new);
                                                        si->last_r_test = i + si->len ;
                                                }
                                        }
                                }
                        }
                }
                //if(i % 100000 == 0){
                //	fprintf(stderr,"\r%8.2f percent done",(double)i /(double)n * 100);
                //}
        }

        //fprintf(stderr,"\r%8.2f percent done\n",100.0);
        //exit(0);
	
        //test++;
        //timer(qs->size,n);
        return qs;
}



unsigned int rev_bits(unsigned int x)
{
        x = ((x & 0x55555555) << 1) | ((x >>1) & 0x55555555);
        x = ((x & 0x33333333) << 2) | ((x >>2) & 0x33333333);
        x = ((x & 0x0F0F0F0F) << 4) | ((x >>4) & 0x0F0F0F0F);
        x = (x << 24) | ((x & 0xFF00) << 8) | ((x >>8) & 0xFF00) | (x >> 24);
        return x;	
}




unsigned long int validate_bpm(unsigned char* t,unsigned char* p,int n,int m,long int offset)
{
        unsigned long int B[4];
        
        if(m > 64){
                m = 64;
        }
        
        int diff = m;
        long int c;
        int i;
        int k = 256;
        register unsigned long int VP,VN,D0,HN,HP,X;
        long int MASK = 0;
        for(i = 0; i < 4;i++){
                B[i] = 0ul;
        }
	
        for(i = 0; i < m;i++){
                B[(int)(p[i] &0x3u)] |= (1ul << (unsigned long int)(m-1-i));
        }
        c = 0ul;
        VP = 0xFFFFFFFFFFFFFFFFul;
        VN = 0ul;
        m--;
        MASK = 1ul << m;
        for(i = 0; i < n;i++){
                X = (B[ (int)(*t  & 0x3u)] | VN);
                D0 = ((VP+(X&VP)) ^ VP) | X ;
                HN = VP & D0;
                HP = VN | ~(VP | D0);
                X = HP << 1ul;
                VN = X & D0;
                VP = (HN << 1ul) | ~(X | D0);
                diff += (HP & MASK) >> m;
                diff -= (HN & MASK) >> m;
                if(diff < k){
                        k = diff;
                        c = i;			
                }
                t--;
        }
        return (((unsigned long int) (m+1 -k)) << 56ul ) | ((unsigned long int)(( offset - c ) << 1ul));
}



unsigned int BinaryInsertionSortOne(unsigned long int *arr, unsigned int size, unsigned long int new){
        unsigned int X;
        unsigned int ret = 1;
        unsigned int l = 0;
        unsigned int r = size - 1;
        unsigned int next;
        if(new <= arr[r]){
                return ret;
        }
        r -= 1;
        ret += 1;
        if(new <= arr[r]){
                arr[size - 1] = new;
                return ret;
        }
        ret += 1;
        if(new >  arr[l]){
                memmove(arr + 1, arr, sizeof(unsigned long int)*(size - 1));
                arr[l] = new;
                return ret;
        }else if(new == arr[l] ){
                return ret;
        }
        
        next = l + ((r - l)>>1);
        for ( X = 0; X < MAX_LOOP_SEARCH_DEPTH; X++ ) {
                ret += 1;
                if(arr[next] >  new){
                        l = next;
                }else if(arr[next] <  new){
                        ret += 1;
                        r = next;
                }else{
                        return ret;
                        ret += 1;
                        memmove(arr + next + 1, arr + next, sizeof(unsigned long int)*(size - 1 - next));
                        arr[next] = new;
                        return ret;
                }
                next = l + ((r - l)>>1);
        }
	
        ret += 1;
        if(arr[next] <  new){
                memmove(arr + next + 1, arr + next, sizeof(unsigned long int)*(size - 1 - next));
                arr[next] = new;
                return ret;
        }else if(arr[next] >   new){
                memmove(arr + r + 1, arr + r, sizeof(unsigned long int)*(size - 1 - r));
                arr[r] = new;
                return ret;
        }
        return ret;
}


void init_nuc_code()
{
        int i;
        for(i = 0;i < 256;i++){
                nuc_code[i] = 0;
                nuc_code5[i] = 4;
                //	rev_nuc_code16[i] = 15;
        }
        nuc_code[65] = 0;//A Adenine
        nuc_code[67] = 1;//C	Cytosine
        nuc_code[71] = 2;//G	Guanine
        nuc_code[84] = 3;//T	Thymine
        nuc_code[85] = 3;//U	Uracil
	
		
        nuc_code[65+32] = 0;//A Adenine
        nuc_code[67+32] = 1;//C	Cytosine
        nuc_code[71+32] = 2;//G	Guanine
        nuc_code[84+32] = 3;//T	Thymine
        nuc_code[85+32] = 3;//U	Uracil
	
        nuc_code5[65] = 0;//A Adenine
        nuc_code5[67] = 1;//C	Cytosine
        nuc_code5[71] = 2;//G	Guanine
        nuc_code5[84] = 3;//T	Thymine
        nuc_code5[85] = 3;//U	Uracil
	
        rev_nuc_code5[0] = 3;//A Adenine
        rev_nuc_code5[1] = 2;//C	Cytosine
        rev_nuc_code5[2] = 1;//G	Guanine
        rev_nuc_code5[3] = 0;//T	Thymine
        rev_nuc_code5[4] = 4;//U	Uracil
	
        nuc_code5[65+32] = 0;//A Adenine
        nuc_code5[67+32] = 1;//C	Cytosine
        nuc_code5[71+32] = 2;//G	Guanine
        nuc_code5[84+32] = 3;//T	Thymine
        nuc_code5[85+32] = 3;//U	Uracil
}


unsigned char* reverse_without_complement(unsigned char* p,int len)
{
        unsigned char* tmp = malloc(sizeof(unsigned char)*MAX_LINE);
        int i,c;
        c = 0;
        for(i =len-1; i >= 0;i--){
                tmp[c] = (int)p[i];
                c++;
        }
        tmp[c] = 0;
        for(i= 0; i < len;i++){
                p[i] = tmp[i];
        }
        free(tmp);
        return p;
}


struct genome* load_genome(struct parameters* param)
{
        struct genome* genome = NULL;
        struct chromosome* chr = NULL;
        FILE* f_ptr = NULL;
        char line[MAX_LINE];
        int i;
        ASSERT(param != NULL, "No parameters.");


        MMALLOC(genome,sizeof(struct genome));
        genome->malloc_num = 64;
        genome->num_chr = -1;
        genome->chromosomes = NULL;
        MMALLOC(genome->chromosomes, sizeof(struct chromosome*) *genome->malloc_num );
        for(i = 0; i < genome->malloc_num;i++){
                genome->chromosomes[i] = NULL;
                chr = NULL;
                MMALLOC(chr,sizeof(struct chromosome));
                chr->seq = NULL;
                chr->name = NULL;
                chr->malloc_len = 65536;
                chr->seq_len = 0;
                MMALLOC(chr->seq, sizeof(uint8_t) * chr->malloc_len);
                MMALLOC(chr->name, sizeof(char) * 256);
                genome->chromosomes[i] = chr;
        }
        
        RUNP(f_ptr = fopen(param->genome, "r" ));
        
        chr = NULL;
        while(fgets(line, MAX_LINE, f_ptr)){
                //reading name .... 
                if(line[0] == '>'){
                        line[strlen(line)-1] = 0;                       
                        genome->num_chr++;
                        if(genome->num_chr){
                                LOG_MSG("read %d bases.",chr->seq_len);                                        
                        }
                        LOG_MSG("Reading in: %s",line);
                        if(genome->num_chr == genome->malloc_num){
                                genome->malloc_num = genome->malloc_num << 1;
                                MREALLOC(genome->chromosomes,sizeof(struct chromosome*) * genome->malloc_num);
                                for(i = genome->num_chr; i < genome->malloc_num;i++){
                                        genome->chromosomes[i] = NULL;
                                        chr = NULL;
                                        MMALLOC(chr,sizeof(struct chromosome));
                                        chr->seq = NULL;
                                        chr->name = NULL;
                                        chr->seq_len = 0;
                                        chr->malloc_len = 65536;
                                        MMALLOC(chr->seq, sizeof(uint8_t) * chr->malloc_len);
                                        MMALLOC(chr->name, sizeof(char) * 256);
                                        genome->chromosomes[i] = chr;
                                }
                                
                        }
                        chr = genome->chromosomes[genome->num_chr];
                        snprintf(chr->name,256,"%s",line);
                }else{
                        for(i = 0;i < MAX_LINE;i++){
                                if(iscntrl((int)line[i])){
                                        break;
                                }
                                chr->seq[chr->seq_len] = nuc_code5[(int)line[i]];
                                chr->seq_len++;
                                if(chr->seq_len == chr->malloc_len){
                                        chr->malloc_len = chr->malloc_len << 1;
                                        MREALLOC(chr->seq, sizeof(uint8_t) *chr->malloc_len);
                                }
                        }			
                }
        }
        LOG_MSG("read %d bases.",chr->seq_len);
        LOG_MSG("read in %d sequences.",genome->num_chr+1);

        fclose(f_ptr);
        return genome;
ERROR:
        if(f_ptr){
                fclose(f_ptr);
        }
        free_genome(genome);
        return NULL;
}

void free_genome(struct genome* genome)
{
        int i;
        struct chromosome* chr = NULL; 
        if(genome){
                for(i =0; i < genome->malloc_num;i++){
                        chr = genome->chromosomes[i];
                        MFREE(chr->seq);
                        MFREE(chr->name);
//                        MFREE(chr->
                        MFREE(chr);
                }
                MFREE(genome->chromosomes);
                MFREE(genome);
        }
}

