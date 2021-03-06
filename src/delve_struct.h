#ifndef DELVE_STRUCT_H

#define DELVE_STRUCT_H


#define ALIGNMENT_FLANKING_LENGTH 5LL

#define NCODE 4
#define IGAP 5
#define EGAP 6


#include "htslib/faidx.h"

struct thread_data{
        struct shared_data* bsd;
        int thread_id;
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
        float temperature;
        int statmode; 
        int iteration;
        int buffer_size;
        int num_threads;
        int num_maxhits;
        int max_seq_len; 
        //float pseudo_counts;
};

struct genome_sequences{
        char** genomic_sequences;
        float* alignment_weigth;
        float* genome_seq_score;
        float* prev_aln_prior;
        float* prior;
        float* prior_e;
        int* g_len;
        float read_score;
        float prev_unaln_prior;
};


struct delve_genome_region_data{
        float count;
        float alignment_weigth; 
        float prior;
        float prior_e;
};

struct parameters{
        char** infiles;
        char* genome;
        char* aln_infile;
        char* hmm_file;
        char* outdir;
        char* id;
	
        int num_infiles;
        int num_maxhits;

        uint8_t num_threads;
        uint8_t devel;
        uint8_t pseudocounts;  
        uint8_t genome_pseudocounts;
        uint8_t siter;
        uint8_t giter;
        uint8_t nogp;
	

};

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
        float* random_genome_model;
        float* random_read_model;
        float* random_genome_model_t;
        float* random_read_model_t;

        float* alignment_scores;
        float* unaligned_scores;
        float* tmp_scores;
	
        float unaligned_genome_score;
        float unaligned_read_score;
        float score;

	
        int x;
        int y;
        int num_hits;
};

#endif
