#ifndef HMM_H

#define HMM_H


struct hmm* init_hmm(int x, int y, int z);
void free_hmm(struct hmm* hmm);

int init_shared_data_hmms(struct shared_data* bsd);
int init_thread_hmms(struct shared_data* bsd);
int entangle_hmms(struct shared_data* bsd);


int add_pseudo_count(struct hmm* hmm, float total);
int re_estimate(struct hmm* hmm);
int add_pseudo_count_random(struct hmm* hmm, float total);
int re_estimate_random(struct hmm* hmm);

#endif

