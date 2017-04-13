#ifndef shared_data_header

#define shared_data_header 


/* related to shared memory */
struct shared_data* init_shared_data(struct parameters* param, int buffer_size);
void free_shared_data(struct shared_data* bsd);


#endif
