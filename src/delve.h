
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <pthread.h>


#if HAVE_CONFIG_H
#include <config.h>
#endif


#ifdef DEBUG
#define MAXNUMQUERY 1000
#else
#define MAXNUMQUERY 1000000
#endif

#define LIST_STORE_SIZE 10

#define NCODE  4
#define IGAP 13
#define EGAP 14
#define BORDER_REGION 2
#define PSEUDO_CONSTANT 1000
#define MAX_LOCAL_POSTERIORS 100000
#define MAX_ERROR_LIMIT_FOR_TRAINING 100
#define MAX_LENGTH_LIMIT_FOR_TRAINING 100

#define RUN_FULL_HMM 1
#define RUN_RANDOM_HMM 2
#define RUN_SCORE_HMM 3

#define ALIGNMENT_FLANKING_LENGTH 10
#define MAX_LINE 10000


#ifdef HUGE_VAL
#define SCALEINFTY HUGE_VAL
#endif





struct hmm{
	float** fM;
	float** fX;
	float** fY;
	
	float*** tfM;
	float*** tfX;
	float*** tfY;
	
	float** bM;
	float** bX;
	float** bY;
	
	unsigned short int** traceback;
	
	float* prob;
	float* tprob;
	float* print_prob;
	float genome_score;
	float random_score;
	float score;
	float* random_genome_model;
	float* random_read_model;
	float* random_genome_model_t;
	float* random_read_model_t;
	
	int x;
	int y;
	int num_hits;
};

/*struct thread_data{
	struct sam_bam_file* sb_file;
	faidx_t*  index;
	struct hmm* hmm;
	struct db* db;
	struct rescue_probabilities* rp;
	struct parameters* param;
	struct genome_probabilities* gp;
	struct genome_sequences** gc;
	struct genome_interval* g_int; 
	int* p_value_ranks;
	FILE* fout;
	int start;
	int end;
};
*/

struct genome_sequences{
	char** genomic_sequences;
	int* g_len;
};


extern struct hmm* init_hmm(int x, int y, int z);

void write_hmm_to_file(struct hmm* hmm,char* filename);
void read_hmm_from_file(struct hmm* hmm,char* filename);
char* transform_16to5_sequence(char* seq,int len);

//int run_delve(struct parameters* param);
//int run_pHMM(struct hmm* localhmm,struct sam_bam_file* sb_file,struct genome_sequences** gc, faidx_t*  index,int num_threads ,int size, int mode, FILE* fout);

//struct hmm* run_pHMM(struct hmm* localhmm,struct read_info** ri,struct db* db,int num_threads ,int size, int mode,FILE* fout);
//int init_thread_hmms(struct thread_data* thread_data_array,struct hmm* localhmm,int numthreads);
//int  entangle_hmms(struct thread_data* thread_data_array,struct hmm* localhmm,int numthreads);
void* do_baum_welch_thread(void *threadarg);
void* do_baum_welch_thread_random_model(void *threadarg);
void free_hmm(struct hmm* hmm);

char* glocal_viterbi_log_Y(struct hmm* hmm,char* a, char* b, int n,int m);

struct hmm* glocal_forward_log_Y(struct hmm* hmm, char* a,  char* b, int n,int m,int frame );
struct hmm* glocal_backward_log_Y(struct hmm* hmm, char* a, char* b, int n,int m);
struct hmm* get_prob_log_Y(struct hmm* hmm, char* a, char* b, int n,int m,float frac,int frame );
struct hmm* random_model_calc(struct hmm* hmm, char* rseq, char* gseq,int rlen, int glen,float frac,int mode);

struct read_info** run_score_alignments_hmm(struct hmm* hmm,struct read_info** ri,struct db* db, int num_threads,int size);
void* do_score_alignments_thread_hmm(void *threadarg);

//void print_samheader(TOMEDB* db,char* command_line);
void unaligned_to_sam(struct read_info* ri);

int align_to_sam(struct pwrite_main* pw,struct genome_interval* g_int,struct sam_bam_entry* entry,int id, char* aln,unsigned flag,float score);

//void align_to_sam(char* aln, TOMEDB* db,unsigned int pos, unsigned flag,float score, struct read_info* ri,unsigned char* posteriors);
char* reverse_path(char* org );
int reverse_complement_sequence(char* p,int len);
//int reverse(char* p,int len);
int re_estimate_random(struct hmm* hmm);

int re_estimate(struct hmm* hmm);
int add_pseudo_count(struct hmm* hmm, float total);

int ACGT_to_0123(char* seq,int* len);

