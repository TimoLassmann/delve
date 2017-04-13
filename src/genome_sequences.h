#ifndef genome_sequence_header

#define genome_sequence_header


/* Structures to cache genome seuqences.. */
struct genome_sequences** init_genome_sequences(int num, int num_maxhits);
int add_genome_sequences(struct shared_data* bsd);
int clear_genome_sequences(struct genome_sequences** gc, int num,int num_maxhits);
void free_genome_sequences(struct genome_sequences** gc, int num,int num_maxhits);

/* misc... */

int convert_buffer_ACGT_to_0123(struct sam_bam_file* sb_file);


#endif
