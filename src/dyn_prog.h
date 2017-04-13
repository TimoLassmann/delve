#ifndef DYN_PROG_H

#define DYN_PROG_H





struct hmm* glocal_forward_log_Y(struct hmm* hmm, char* a,  char* b, int n,int m,int frame);
struct hmm* glocal_backward_log_Y(struct hmm* hmm, char* a, char* b, int n,int m);

char* glocal_viterbi_log_Y(struct hmm* hmm,char* a, char* b, int n,int m);

struct hmm* get_prob_log_Y(struct hmm* hmm,char* a, char* b, int n,int m,float frac,int frame);

int random_model_genome_score(struct hmm* hmm, char* seq,int len,float frac);
int random_model_read_score(struct hmm* hmm, char* seq,int len,float frac);

#endif
