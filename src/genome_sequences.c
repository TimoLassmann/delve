

/* my libs  */
#include "tldevel.h"
#include "htsglue.h"

/* delve headers  */

#include "delve_struct.h"
#include "genome_sequences.h"



static int ACGT_to_0123(char* seq,int* len);

struct genome_sequences** init_genome_sequences(int num, int num_maxhits)
{
	struct genome_sequences** gc = NULL;
	int i,j;

	MMALLOC(gc,sizeof(struct genome_sequences*) * num);
	for(i = 0; i < num;i++){
		gc[i] = NULL;
		MMALLOC(gc[i], sizeof(struct genome_sequences));
		gc[i]->genomic_sequences = NULL;
		gc[i]->g_len = NULL;
		gc[i]->alignment_weigth = NULL;
		gc[i]->prior = NULL;
		gc[i]->prior_e = NULL;
		
		MMALLOC(gc[i]->g_len,sizeof(int) * num_maxhits);
		MMALLOC(gc[i]->genomic_sequences,sizeof(char*) * num_maxhits);
		MMALLOC(gc[i]->alignment_weigth,sizeof(float) * num_maxhits);
		MMALLOC(gc[i]->prior,sizeof(float) * num_maxhits);
		MMALLOC(gc[i]->prior_e,sizeof(float) * num_maxhits);
		for (j = 0; j < num_maxhits; j++) {
			gc[i]->g_len[j] = 0;
			gc[i]->genomic_sequences[j] = NULL;
			gc[i]->alignment_weigth[j] = 0.0f;
			gc[i]->prior[j] = 0.0f;
			gc[i]->prior_e[j] = 0.0f;
		}
	}
	return gc;
ERROR:
	free_genome_sequences(gc,num,num_maxhits);			    
	return NULL;
	 
}

int add_genome_sequences(struct shared_data* bsd)
{
//	faidx_t*  index = NULL;
	struct sam_bam_file* sb_file = NULL;
	struct genome_sequences** gc = NULL;
	faidx_t*  index = NULL;
	
	struct genome_interval* g_int = NULL;

	int i,j,len,c;

	
	sb_file = bsd->sb_file;
	gc = bsd->gc;
	index = bsd->index;
	RUNP(g_int = init_genome_interval(0,0,0));


	for(i = 0; i < sb_file->num_read;i++){
		for(j = 0; j < sb_file->buffer[i]->num_hits;j++){
			DPRINTF2("NEWHIT:");
			DPRINTF2("%s: %d ->%d", sb_file->buffer[i]->name,sb_file->buffer[i]->start[j],sb_file->buffer[i]->stop[j] );
			
			RUN(get_chr_start_stop(sb_file->si  ,g_int,sb_file->buffer[i]->start[j], sb_file->buffer[i]->stop[j]));
			g_int->start -=ALIGNMENT_FLANKING_LENGTH;
			g_int->stop += ALIGNMENT_FLANKING_LENGTH;
			DPRINTF2("%s:%d-%d ",g_int->chromosome,g_int->start,g_int->stop);

			gc[i]->genomic_sequences[j] = get_sequence(index,g_int);
			
			len = (int) strlen(gc[i]->genomic_sequences[j]);

			gc[i]->g_len[j] = len;
			//DPRINTF2("%s	(len:%d)",sb_file->buffer[i]->genomic_sequences[j],sb_file->buffer[i]->g_len[j]);

			RUN(ACGT_to_0123(gc[i]->genomic_sequences[j] ,&len));
			
			for(c =0; c < len;c++){
				ASSERT(gc[i]->genomic_sequences[j][c] < 5, "%d %d pos:%d  nuc:%d  ", i,j, c, gc[i]->genomic_sequences[j][c] ); 
			}		
		}
	}
	free_genome_interval(g_int);
	return OK;
ERROR:
	if(g_int){
		free_genome_interval(g_int);
	}
	return FAIL;
}


int clear_genome_sequences(struct genome_sequences** gc, int num,int num_maxhits)
{
	int i,j;
	for(i = 0;i< num;i++){
		if(gc[i]){
			for(j = 0; j < num_maxhits;j++){
				if(gc[i]->genomic_sequences[j]){
					MFREE(gc[i]->genomic_sequences[j]);
					gc[i]->genomic_sequences[j] = NULL;
				}
				gc[i]->g_len[j] = 0;
				       
			}
		}
				      
	}
	
	
	return OK;		
}


void free_genome_sequences(struct genome_sequences** gc, int num,int num_maxhits)
{
	int i,j;
	if(gc){
		for(i = 0;i< num;i++){
			if(gc[i]){
				for(j = 0; j < num_maxhits;j++){
					if(gc[i]->genomic_sequences[j]){
						MFREE(gc[i]->genomic_sequences[j]);
					}
				}
				MFREE(gc[i]->alignment_weigth);
				MFREE(gc[i]->prior);
				MFREE(gc[i]->prior_e);
				MFREE(gc[i]->genomic_sequences);
				MFREE(gc[i]->g_len);
				MFREE(gc[i]);
			}
				      
		}
		MFREE(gc);
	}
}

int convert_buffer_ACGT_to_0123(struct sam_bam_file* sb_file)

{
	int i;
	for(i = 0; i < sb_file->num_read;i++){
		RUN(ACGT_to_0123(sb_file->buffer[i]->sequence,&sb_file->buffer[i]->len));
	}
	return OK;
ERROR:
	return FAIL;
}




int ACGT_to_0123(char* seq,int* len)
{
	int i = 0;
	int c = 0;
	while(seq[i] !=0){
		switch(seq[i]){
			case 'A':
			case 'a':
				seq[c] = 0;
				c++;
				break;
			case 'C':
			case 'c':
				seq[c] = 1;
				c++;
				break;
			case 'G':
			case 'g':
				seq[c] = 2;
				c++;
				break;
			case 'T':
			case 't':
				seq[c] = 3;
				c++;
				break;
			default:
				seq[c] = NCODE;
				c++;
				//ERROR_MSG("Non ACGT letter in sequence:%s.",seq);
				break;
		}
		i++;
	}
	seq[c] = 0;
	*len = c;
	return OK;
}
