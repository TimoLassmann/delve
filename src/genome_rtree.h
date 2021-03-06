#ifndef genome_rtree_header

#define genome_rtree_header


struct rtr_data* build_rtree(struct sam_bam_file* sb_file );
int set_sequence_weigth(struct shared_data* bsd);
int free_delve_region_data_from_tree(struct rtr_data* rtree);


#endif
