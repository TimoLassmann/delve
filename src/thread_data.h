#ifndef THREAD_DATA_H

#define THREAD_DATA_H

/* related to thread data  */
struct thread_data** init_thread_data(struct shared_data* bsd,int num_threads);
void free_thread_data(struct thread_data** td,int num_threads);

int get_start_stop_for_threads(const int num_threads, const int num_seq, const int id, int *start,int *stop);

#endif
