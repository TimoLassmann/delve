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



#define MAX_LOOP_SEARCH_DEPTH 10 /*  */

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

struct seq_hit{
        uint64_t pos : 48;
        uint8_t strand : 1;
        uint8_t error: 4;
        uint16_t chr: 11;
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
        //unsigned long int* match_units;
        struct seq_hit hits[LIST_STORE_SIZE];
        //int hit_insert_pos; 
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
        struct genome* genome;
        long int len;
        long int total_len;
        unsigned char* chromsome;
        int start;
        int end;
        int thread;
};


int seed_controller(struct parameters* param);

int run_search(struct parameters* param, struct qs_struct* qs, struct genome* genome);
void* do_run_search(void *threadarg);
struct qs_struct* init_query_structure(struct parameters* param);
void free_query_structure(struct qs_struct * qs);

//struct qs_struct* search_fasta_fastq(struct qs_struct* qs,struct parameters* param);

void reverse_complement(unsigned char* p,unsigned char* rc,int len);

struct qs_struct* make_hash(struct qs_struct* qs,int r);
struct qs_struct* search_hash(struct qs_struct* qs,unsigned char* t,int n, long int r,long int offset);
void print_sam(struct seq_info* si, struct genome* genome,FILE* f_ptr);//  unsigned char** t_names, long int* lengths,int size);
int insert_hit(struct seq_info* si, struct seq_hit* new);
int seed_controller_thread(struct parameters* param);
//struct qs_struct* search_fasta_fastq_thread(struct qs_struct* qs,struct parameters* param);
struct qs_struct* make_hash_thread(struct qs_struct* qs,int r,int step,int thread,int start,int stop);
void* search_hash_thread(void *threadarg);

struct qs_struct* make_B(struct qs_struct* qs,int n);
struct qs_struct* read_fasta_fastq2(FILE *file, struct qs_struct* qs,struct parameters* param) ;
unsigned int rev_bits(unsigned int x);
int validate_bpm(unsigned char* t,unsigned char* p,int n,int m,long int offset,struct seq_hit* hit);
//unsigned int BinaryInsertionSortOne(unsigned long int *arr, unsigned int size, unsigned long int new);


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
                if(qs->size){
                        LOG_MSG("%d read.",qs->size);
                        //qs = make_B(qs,qs->size);
                        interval =  (int)((double)qs->size /(double)param->num_threads);
		
                        for(t = 0;t < param->num_threads - 1;t++) {
                                qs = make_hash_thread(qs,param->seed_len,param->seed_step,t, t*interval,t*interval+interval);

                        }
                        qs = make_hash_thread(qs,param->seed_len,param->seed_step,param->num_threads - 1, t*interval,qs->size);
                        RUN(run_search(param, qs,genome));
                        //qs = search_fasta_fastq_thread(qs,param);
                
                        for(t = 0;t < param->num_threads ;t++){
                                free(qs->hash_t[t]);
                                free(qs->hash_index_t[t]);
                        }
                }
                
        }


        fclose(fastafile);

        
        
        free_query_structure(qs);
        
        free_genome(genome);
        return OK;
ERROR:
        return FAIL;
}


int run_search(struct parameters* param, struct qs_struct* qs, struct genome* genome)
{
        static int header = 1;
        FILE* f_ptr = NULL;
        struct seed_thread_data** thread_data_array = NULL;

        int i,j;
        int status;
        ASSERT(param != NULL,"No param.");
        ASSERT(qs != NULL,"No param.");
        ASSERT(genome != NULL,"No param.");
        
        MMALLOC(thread_data_array ,sizeof(struct seed_thread_data*)* param->num_threads);

        for(i = 0;i < param->num_threads;i++){
                thread_data_array[i] = NULL;
                MMALLOC(thread_data_array[i],sizeof(struct seed_thread_data));                
                thread_data_array[i]->param = param;
                thread_data_array[i]->qs = qs;
                thread_data_array[i]->genome = genome;
                thread_data_array[i]->thread = i;
                if((status = thr_pool_queue(param->pool,do_run_search,thread_data_array[i])) == -1) fprintf(stderr,"Adding job to queue failed.");	
        }
        thr_pool_wait(param->pool);

        
        if(header){
                RUNP(f_ptr = fopen(param->output, "w"));
                fprintf(f_ptr, "@HD\tVN:1.0\n");
		
                for(i = 0; i < genome->num_chr ;i++){
                        fprintf(f_ptr, "@SQ\tSN:%s\tLN:%u\n", genome->chromosomes[i]->name, (unsigned int) genome->chromosomes[i]->seq_len);
                        
                        //fprintf(stdout,"\n");
                        
                }
                fprintf(f_ptr, "@PG\tID:Delve\tVN:%s\tCL:%s\n",VERSION,"dseed something") ;
                header = 0;
        }else{
                RUNP(f_ptr = fopen(param->output, "a"));
        }
        
        
        //total_len +=  len;
        for(i = 0; i < qs->size;i++){
                print_sam(qs->seq_info[i], genome,f_ptr);//_names, lengths ,c);
        }
        fclose(f_ptr);
	
        //exit(0);
        for(i = 0; i< qs->size;i++){
                for(j = 0; j < LIST_STORE_SIZE;j++){
                        qs->seq_info[i]->hits[j].error = 0xF;
                }
        }
        //string_pos = 0;
        qs->q_num = 0;

        for (i = 0; i < param->num_threads; i++){
                MFREE(thread_data_array[i]);
        }
                     
        MFREE(thread_data_array); 
        
        return OK;
ERROR:
        if(f_ptr){
                fclose(f_ptr);
        }
        return FAIL;
}

void* do_run_search(void *threadarg)
{
        struct seed_thread_data* data = NULL;
        struct parameters* param =  NULL;
        struct genome* genome = NULL;
        struct qs_struct* qs = NULL;
        int thread_id; 

        ASSERT(threadarg != NULL,"No threadarg");
        data = (struct seed_thread_data*) threadarg;
        param = data->param;
        thread_id = data->thread;
        genome = data->genome;
        qs= data->qs;

        unsigned char* t = NULL;// = data->chromsome;
        long int n;//  data->len;
        long int r = (long int) param->seed_len;
        
        	
        int j;//,g;
        int i;
        uint8_t* target = 0;
        struct seq_info* si = 0;
        unsigned int key = 0;
        unsigned int mask = 0;
        int t_len = 0;
        long int seed_offset = 0;
        int strand = 0;
        long int add = 0;
	
        int* index = qs->hash_index_t[thread_id];
        int chr;

        struct seq_hit tmp_hit;
        
        for(chr = 0; chr <=  genome->num_chr;chr++){

                LOG_MSG("Thread %d is working on %s.", thread_id, genome->chromosomes[chr]->name);
                t = genome->chromosomes[chr]->seq;
                n = genome->chromosomes[chr]->seq_len;

                for(i = 0; i< qs->size;i++){
                        qs->seq_info[i]->last_f_test = -1000;
                        qs->seq_info[i]->last_r_test = -1000;
                }
                
                mask = 0;
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
                                        seed_offset = qs->hash_t[thread_id][j] & 0x7fu;
				
                                        si = qs->seq_info[ qs->hash_t[thread_id][j] >> 8u];

                                        strand = qs->hash_t[thread_id][j] & 0x80u;
                                        add = ( (long int)i + ((long int) si->len - seed_offset)  -r + 10l); /// should be 10 ??? 

                                        target = t +add;  

                                        t_len = si->len + 10;

                                        if( add <= n && add > ((long int) (si->len + 20))){
                                                if(!strand){
                                                        if(si->last_f_test <= i){
                                                                //fprintf(stderr,"ADD:%ld	LEN:%d\n",add,t_len);
                                                                RUN(validate_bpm(target ,si->seq , t_len , si->len, add, &tmp_hit));
                                                                tmp_hit.chr = chr;
                                                                tmp_hit.strand = 0;
                                                                RUN(insert_hit(si,&tmp_hit));
                                                                //BinaryInsertionSortOne(si->match_units,LIST_STORE_SIZE ,new);
                                                                si->last_f_test = i + si->len  ;
                                                        }
                                                }else{
                                                        if(si->last_r_test <= i){
                                                                RUN(validate_bpm(target ,si->reverse_seq , t_len , si->len, add,&tmp_hit));
                                                                tmp_hit.chr = chr;                                                                
                                                                tmp_hit.strand = 1;
                                                                RUN(insert_hit(si,&tmp_hit));
                                                                //BinaryInsertionSortOne(si->match_units,LIST_STORE_SIZE,new);
                                                                si->last_r_test = i + si->len ;
                                                        }
                                                }
                                        }
                                }
                        }
                }
                LOG_MSG("Thread %d is done.",thread_id);
        }
        return NULL;
ERROR:
        return NULL;
}


int insert_hit(struct seq_info* si, struct seq_hit* new)
{
        int p = 0;
        /* check if place to insert is free / hit there is  */
        ASSERT(si != NULL,"No sequence info.");

        p = LIST_STORE_SIZE-1;

        if(si->hits[p].error > new->error){
                si->hits[p].chr = new->chr;
                si->hits[p].strand = new->strand;
                si->hits[p].pos = new->pos;
                si->hits[p].error = new->error;
        }

        while( p != 0){
                if(si->hits[p-1].error > si->hits[p].error){
                        new->chr = si->hits[p-1].chr;
                        new->strand = si->hits[p-1].strand;
                        new->pos = si->hits[p-1].pos;
                        new->error = si->hits[p-1].error;

                        si->hits[p-1].chr = si->hits[p].chr;
                        si->hits[p-1].strand = si->hits[p].strand;
                        si->hits[p-1].pos = si->hits[p].pos;
                        si->hits[p-1].error = si->hits[p].error;
                        
                        si->hits[p].chr = new->chr;
                        si->hits[p].strand = new->strand;
                        si->hits[p].pos = new->pos;
                        si->hits[p].error = new->error;
                }else{
                        break;
                }
                p--;
        }
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
                

                for(j = 0; j < LIST_STORE_SIZE;j++){
                        qs->seq_info[i]->hits[j].chr = 0;
                        qs->seq_info[i]->hits[j].pos = 0;
                        qs->seq_info[i]->hits[j].strand = 0;
                        qs->seq_info[i]->hits[j].error = 0xF;
                        
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
                        MFREE(qs->seq_info[i]);
		
                }

                MFREE(qs->hash_t);
                MFREE(qs->hash_index_t);

                MFREE(qs->seq_info);  
                MFREE(qs);
        }
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
                        qs->seq_info[i]->hits[j].chr = 0;
                        qs->seq_info[i]->hits[j].pos = 0;
                        qs->seq_info[i]->hits[j].strand = 0;
                        qs->seq_info[i]->hits[j].error = 0xF;
                        
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


void print_sam(struct seq_info* si, struct genome* genome,FILE* f_ptr)// char** t_names, long int* lengths,int size)
{
        char alphabet[] = "ACGTN";
        char* pline = malloc(sizeof (char)* MAX_LINE*10);
        char seq[MAX_LINE];
        unsigned char star[] = "*";
        int strand = 0;
        int no_qual = 0;
        int i,j;
        char* t_name =  0;
        long int pos = 0;
        
        int errors = 0xf;
	
        pline[0] = 0 ;
        /*strand = si->match_units[0] & 1ul;
        if(strand){
                strand = 16;
                }*/
	
        //fprintf(stderr,"before : %s\n",si->qual);
	
        if(!si->qual){
                si->qual = star;
                no_qual = 1;
        }
	
        //fprintf(stderr,"after : %s\n",si->qual);
	
        if(si->hits[0].error == 0xF){
		
		
		
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
                fprintf(f_ptr, "%s\n",pline);
        }else{
                t_name = genome->chromosomes[si->hits[0].chr ]->name;
                pos = si->hits[0].pos;//   (long int) ((si->match_units[0] >> 11ul )& 0xFFFFFFFFul);
                errors = si->hits[0].error;
                strand = si->hits[0].strand;
                if(strand){
                        strand = 16;
                }
                //matches = (int) (si->match_units[0] >> 56ul);
                //fprintf(stderr,"printing:%lu	%ld\n",si->match_units[0],pos);
                		
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
                        sprintf(pline,"%s\t%d\t%s\t%d\t%d\t%dM\t%s\t%d\t%d\t%s\t%s",
                                si->sn ,//1 QNAME Query NAME of the read or the read pair 
                                strand,//2 FLAG bitwise FLAG (pairing, strand, mate strand, etc.) 
                                t_name,//3 RNAME Reference sequence NAME 
                                (int)(pos+1l),//4 POS 1-based leftmost POSition of clipped alignment 
                                255,//5 MAPQ MAPping Quality (Phred-scaled) 
                                si->len,//si->len,//6 CIGAR extended CIGAR string (operations: MIDNSHP) 
                                "*",//7 MRNM Mate Reference NaMe (‘=’ if same as RNAME) 
                                0,//8 MPOS 1-based leftmost Mate POSition 
                                0,//9 ISIZE inferred Insert SIZE 
                                seq,//10 SEQ query SEQuence on the same strand as the reference 
                                si->qual//11 QUAL query QUALity (ASCII-33=Phred base quality) 
                                );
                }else{
                        if(no_qual){
                                sprintf(pline,"%s\t%d\t%s\t%d\t%d\t%dM\t%s\t%d\t%d\t%s\t%s",
                                        si->sn ,//1 QNAME Query NAME of the read or the read pair 
                                        strand,//2 FLAG bitwise FLAG (pairing, strand, mate strand, etc.) 
                                        t_name,//3 RNAME Reference sequence NAME 
                                        (int)(pos),//4 POS 1-based leftmost POSition of clipped alignment 
                                        255,//5 MAPQ MAPping Quality (Phred-scaled) 
                                        si->len,//si->len,//6 CIGAR extended CIGAR string (operations: MIDNSHP) 
                                        "*",//7 MRNM Mate Reference NaMe (‘=’ if same as RNAME) 
                                        0,//8 MPOS 1-based leftmost Mate POSition 
                                        0,//9 ISIZE inferred Insert SIZE 
                                        seq,//10 SEQ query SEQuence on the same strand as the reference 
                                        si->qual//11 QUAL query QUALity (ASCII-33=Phred base quality) 
                                        );
                        }else{
                                sprintf(pline,"%s\t%d\t%s\t%d\t%d\t%dM\t%s\t%d\t%d\t%s\t%s",
                                        si->sn ,//1 QNAME Query NAME of the read or the read pair 
                                        strand,//2 FLAG bitwise FLAG (pairing, strand, mate strand, etc.) 
                                        t_name,//3 RNAME Reference sequence NAME 
                                        (int)(pos),//4 POS 1-based leftmost POSition of clipped alignment 
                                        255,//5 MAPQ MAPping Quality (Phred-scaled) 
                                        si->len,//6 CIGAR extended CIGAR string (operations: MIDNSHP) 
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
                sprintf(pline+i,"NM:i:%d",  errors);// si->len - matches);
		
		
                if(si->hits[1].error != 0xF){
			
                        i = strlen(pline);
                        pline[i] = '\t';
                        i++;
                        sprintf(pline+i,"XA:Z:");
                        j = 1;
                        while( j < LIST_STORE_SIZE){
                                if(si->hits[j].error == 0xF){
                                        break;
                                }
                                t_name = genome->chromosomes[si->hits[j].chr]->name;
                                pos = si->hits[j].pos;//   (long int) ((si->match_units[0] >> 11ul )& 0xFFFFFFFFul);
                                errors = si->hits[j].error;
                                strand = si->hits[j].strand;
                                
                             
                 
                                i = strlen(pline);
                                if(!strand){
                                        sprintf(pline+i,"%s,+%d,%dM,%d;", t_name,(int)(pos+1l),si->len, errors );
                                }else{
                                        sprintf(pline+i,"%s,-%d,%dM,%d;", t_name,(int)(pos+1l),si->len, errors);
                                }
                                j++;
                        }
			
                }
		
                i = strlen(pline);
                pline[i] = 0;
                fprintf(f_ptr, "%s\n",pline);
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

unsigned int rev_bits(unsigned int x)
{
        x = ((x & 0x55555555) << 1) | ((x >>1) & 0x55555555);
        x = ((x & 0x33333333) << 2) | ((x >>2) & 0x33333333);
        x = ((x & 0x0F0F0F0F) << 4) | ((x >>4) & 0x0F0F0F0F);
        x = (x << 24) | ((x & 0xFF00) << 8) | ((x >>8) & 0xFF00) | (x >> 24);
        return x;	
}

int validate_bpm(unsigned char* t,unsigned char* p,int n,int m,long int offset ,struct seq_hit* hit)
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
        hit->error = 0xf;
        if(k < 0xF){
                hit->error = k;
        }
        hit->pos = (uint64_t) (offset - c);
        return OK;
        // return (((unsigned long int) (m+1 -k)) << 56ul ) | ((unsigned long int)(( offset - c ) << 11ul));
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
                        for(i =0 ; i<strlen(line);i++){
                                if(isspace(line[i])){
                                        line[i] = 0;
                                                
                                }
                                        
                        }

                        genome->num_chr++;
                        if(genome->num_chr){
                                LOG_MSG("read %d bases.",chr->seq_len);                                        
                        }
                        LOG_MSG("Reading in: %s",line+1);
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
                        snprintf(chr->name,256,"%s",line+1);
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

