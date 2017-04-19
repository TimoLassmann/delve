#ifndef INPUT_OUTPUT_H

#define INPUT_OUTPUT_H

int set_output_file(struct shared_data* bsd, char* suffix);

int align_to_sam(struct pwrite_main* pw,struct genome_interval* g_int,struct sam_bam_entry* entry,int id, char* aln,unsigned flag,float score);
int write_sam_header(struct sam_bam_file* sb_file, FILE* out);
int reverse_complement_sequence(char* p,int len);
int reverse(uint8_t* p,int len);
char* reverse_path(char* org);

#endif

