

#include "tldevel.h"

#include "delve_struct.h"

struct thread_data** init_thread_data(struct shared_data* bsd,int num_threads)
{
	int i;
	struct thread_data** td = NULL;
		
	MMALLOC(td,sizeof(struct thread_data*) * num_threads);

	for(i = 0; i < num_threads;i++){
		td[i] = NULL;
		MMALLOC(td[i],sizeof(struct thread_data));
		td[i]->thread_id = i;
		td[i]->bsd = bsd;
 	}
	return td;
ERROR:
	return NULL;
}


void free_thread_data(struct thread_data** td, int num_threads)
{
	int i;
	if(td){
		for(i = 0; i < num_threads;i++){
			if(td[i]){
				MFREE(td[i]);
			}
		}
		MFREE(td);
	}
}



int get_start_stop_for_threads(const int num_threads, const int num_seq, const int id, int *start,int *stop)
{

	int interval = 0;

	if(num_threads > num_seq){
		if(id > num_seq){
			*start = 0;
			*stop = 0;
		}else{
			*start = id;
			*stop = id+1;
		}
		return OK;
	}
	interval = num_seq / num_threads;
	*start = id* interval;
	*stop =  id* interval + interval;
//	LOG_MSG("Trying to spread work: %d total_threads %d.",id, num_threads);
	if(id +1 == num_threads){
		*stop = num_seq;
	}	
	return OK;
}
