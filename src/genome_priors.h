#ifndef GENOME_PRIOR_H

#define GENOME_PRIOR_H


int reset_prev_aln_priors( struct shared_data* bsd);

int add_genome_pseudocounts(struct shared_data* bsd);
int re_estimate_genome_priors(struct shared_data* bsd);

int set_uniform_genome_priors(struct shared_data* bsd);
int copy_genome_priors_to_gc(struct shared_data* bsd);
int entangle_genome_priors(struct shared_data* bsd);

#endif
