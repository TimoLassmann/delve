

#include <ctype.h>

#include "tldevel.h"
#include "htsglue.h"


#include "delve_struct.h"

#include "pwrite.h"


int set_output_file(struct shared_data* bsd, char* suffix)
{
	struct parameters* param = NULL;
	char buffer[BUFFER_LEN];
	char suf_buf[BUFFER_LEN];
	int64_t develcode;
	ASSERT(bsd != NULL,"No shared data.");
	ASSERT(suffix != NULL,"No suffix.");
	
	param = bsd->param;
	
	if(bsd->pw){
		bsd->pw->free(bsd->pw);
		bsd->pw = NULL;
	}
	/* add additional suffix to capture parameters while developing  */
	if(param->devel){
		develcode = 0;
		develcode |= ((int64_t) bsd->param->num_threads) << (int64_t)(8*4);
		develcode |= ((int64_t)bsd->param->siter) << (int64_t)(8*3);
		develcode |= ((int64_t)bsd->param->giter) << (int64_t)(8*2);
		develcode |= ((int64_t)bsd->param->pseudocounts) << (int64_t)(8*1);
		develcode |= ((int64_t)bsd->param->genome_pseudocounts) << (int64_t)(8*0);
	        
		snprintf(suf_buf,BUFFER_LEN,"DCODE%016jX.%s",develcode, suffix);
	}else{
		snprintf(suf_buf,BUFFER_LEN,"%s",suffix);
	}


	if(param->outdir){
		snprintf(buffer,BUFFER_LEN,"%s/%s.%s",param->outdir,basename(param->aln_infile),suf_buf);		
	}else{
		snprintf(buffer,BUFFER_LEN,"%s.%s",param->aln_infile,suf_buf);
	}
	
	LOG_MSG("Set output file to: %s",buffer);
	RUNP(bsd->pw = init_pwrite_main(buffer,param->num_threads,BUFFER_P_WRITE_SIZE));

	return OK;
ERROR:
	return FAIL; 
}

int align_to_sam(struct pwrite_main* pw,struct genome_interval* g_int,struct sam_bam_entry* entry,int id, char* aln,unsigned flag,float score)
{
	char nuc[] = "ACGTN";
	char cigarline[128];
	char mdline[128];
	char pseq[500];
	//unsigned char* qual = 0;
	
	//unsigned char noqual[] = "*";
	int mismatches = 0;
	int i,j,c,g,t;
	unsigned int add = 0;
	c = 0;
	int state = 0;
	mdline[0] = 0;
	c = aln[0];
	j = 0;
	t = 1;
        
	while((aln[c]  >> 4) == 5 ){
	//	fprintf(stderr,"%d	%d	len: %d\n",((aln[c]   >> 4)& 0xF ),(aln[c]  & 0xF),entry->len);
		c--;
	}
	for(i = c;i >= 1;i--){
	//	fprintf(stderr,"%d	%d	%d\n",((aln[i]   >> 4)& 0xF ),(aln[i]  & 0xF),t);
		if((aln[i]  >> 4) != 5 && (aln[i]  & 0xF) !=5){
			if((aln[i]  >> 4) != (aln[i]  & 0xF)){
				g = (int)strlen(mdline);
				snprintf(mdline+g,128,  "%d%c",j,nuc[(aln[i]  & 0xF)]);
				if((aln[i]  & 0xF) > 4 ){
					
					fprintf(stderr,"ERROR:%d	%d\n",aln[i] ,aln[i]  & 0xF) ;
				}
	//			fprintf(stderr,"%d	%c\n",j,nuc[(aln[i]  & 0xF)]);
				j = 0;
			}else{
				j++;
			}
			
			if(t == entry->len){
				i = -1;
				break;
			}
			t++;
			state = 0;
		}else if((aln[i]  >> 4) != 5 && (aln[i]  & 0xF) == 5){
			
			
			
			if(t == entry->len){
				i = -1;
				break;
			}
			t++;
			
	//		fprintf(stderr,"%c	\n",nuc[(aln[i]  & 0xF)]);
			//j = 0;
			state = 0;
		}else if((aln[i]  >> 4) == 5  && (aln[i]  & 0xF) != 5){
			if(state == 0){
				g = (int)strlen(mdline);
				snprintf(mdline+g,128,"%d^",j);
	//			fprintf(stderr,"%d	\n",j);
				state = 1;
			}
			g = (int)strlen(mdline);
			snprintf(mdline+g,128,"%c",nuc[(aln[i]  & 0xF )]);
			j = 0;
			//state = 0;
		}
	}
	
	if(j){
		g = (int)strlen(mdline);
		snprintf(mdline+g,128,"%d",j);
	//		fprintf(stderr,"%d	\n",j);
	}
	
	g = (int)strlen(mdline);
	if(!isdigit((int)  mdline[g-1])){
		snprintf(mdline+g,128,"%d",0);
	}
	
	//fprintf(stderr,"LINE:%s\n\n",mdline);
	
	
	i =  aln[0];
	t = 0;
	cigarline[0] = 0;
	while((aln[i]  >> 4) == 5){
		i--;
		add++;
	}
	
	while(i >= 1){
		j = 0;
		if((aln[i]  >> 4) != 5 && (aln[i]  & 0xF) != 5){
			while ( (aln[i]  >> 4) != 5 && (aln[i]  & 0xF) != 5 && i >= 1) {
				if((aln[i]  >> 4) != (aln[i]  & 0xF) ){
					mismatches++;
				}
				i--;
				j++;
				t++;
				
			}
			c = (int)strlen(cigarline);
			snprintf(cigarline+c,128,"%dM",j);
			if(t ==entry->len){
				break;
			}
			
		}else if((aln[i]  >> 4) != 5 && (aln[i]  & 0xF) == 5 ){
			while ((aln[i]  >> 4) != 5 && (aln[i]  & 0xF) == 5&& i >= 1) {
				i--;
				j++;
				t++;
				mismatches++;
			}
			c = (int)strlen(cigarline);
			snprintf(cigarline+c,128,"%dI",j);
			if(t ==entry->len){
				break;
			}
		}else if((aln[i]  >> 4) == 5 && (aln[i]  & 0xF) != 5 ){
			while ((aln[i]  >> 4) == 5 && (aln[i]  & 0xF) != 5&& i >= 1) {
				i--;
				j++;
				mismatches++;
			}
			c = (int)strlen(cigarline);
			snprintf(cigarline+c,128,"%dD",j);
		}
		//fprintf(stderr,"%d	%s	%d	%d	%d\n",i,cigarline,aln[i]>>4, aln[i] & 0xF,add);
	}
	//fprintf(stderr,"CIGAR ::: %s\n",cigarline);
	
	for(i = 0;i < entry->len;i++){
		if(entry->sequence[i] == 12){
			pseq[i] = 'N';
		}else{
			pseq[i] =  "ACGTN"[(int)entry->sequence[i]];//   nuc[(int)ri->seq[i] & 3 ];
		}
	}
	pseq[entry->len]  = 0;
	
	if(1.0f - score < 0.0001f){
		score = 0.0001f;
	}else{
		score = 1.0f - score;
	}
	
	//get_chr_name (seq_info,(unsigned int) pos),
	//get_chr_pos (seq_info,(unsigned int) pos),
	
	if(flag & 0x100){
		RUN(pw->write(pw,id,"%s\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%s\tNM:i:%d\tMD:Z:%s\n",
			      entry->name,//1 QNAME Query NAME of the read or the read pair
			      flag,//2 FLAG bitwise FLAG (pairing, strand, mate strand, etc.)
			      g_int->chromosome,
			      g_int->start+add +1,
			      //get_chr_name (db,(unsigned int) pos + add + 1), //3 RNAME Reference sequence NAME
			      //get_chr_pos (db,(unsigned int) pos+add +1) ,//4 POS 1-based leftmost POSition of clipped alignment
			      (int)((-10.0f * log10(score ))+0.5f),//5 MAPQ MAPping Quality (Phred-scaled)
			      cigarline,//6 CIGAR extended CIGAR string (operations: MIDNSHP)
			      "*",//7 MRNM Mate Reference NaMe (‘=’ if same as RNAME)
			      0,//8 MPOS 1-based leftmost Mate POSition
			      0,//9 ISIZE inferred Insert SIZE
			      "*" ,//10 SEQ query SEQuence on the same strand as the reference
			      "*",//11 QUAL query QUALity (ASCII-33=Phred base quality)
			      mismatches,
			      mdline
			    ));
	}else{
		RUN(pw->write(pw,id,"%s\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%s\tNM:i:%d\tMD:Z:%s\n",
			      entry->name,//1 QNAME Query NAME of the read or the read pair
			      flag,//2 FLAG bitwise FLAG (pairing, strand, mate strand, etc.)
			      g_int->chromosome,
			      g_int->start+add +1,
			      (int)((-10.0f * log10(score))+0.5f),//5 MAPQ MAPping Quality (Phred-scaled)
			      cigarline,//6 CIGAR extended CIGAR string (operations: MIDNSHP)
			      "*",//7 MRNM Mate Reference NaMe (‘=’ if same as RNAME)
			      0,//8 MPOS 1-based leftmost Mate POSition
			      0,//9 ISIZE inferred Insert SIZE
			      pseq ,//10 SEQ query SEQuence on the same strand as the reference
			      entry->base_qual,//11 QUAL query QUALity (ASCII-33=Phred base quality)
			      mismatches,
			      mdline
			    ));
	}

	
	return OK;
ERROR:
	return FAIL;
}





char* reverse_path(char* org )
{
	char* new = NULL;
	
	int i,c;
	int nuc[5] = {3,2,1,0,4};
	MMALLOC(new, sizeof(char)* (org[0]+1));
	
	new[0] = org[0];
	c = 1;
	for (i = org[0]; i > 0; i--){
		if(( org[i] >> 4) < 5){
			if((org[i] & 0xF) < 5){
				new[c] = (nuc[org[i] >> 4] << 4) |  nuc[org[i] & 0xF];
			}else{
				new[c] = (nuc[org[i] >> 4] << 4) |  (org[i] & 0xF);
			}
		}else{
			if((org[i] & 0xF) < 5){
				new[c] = (org[i] & 0xF0) |  nuc[org[i] & 0xF];
			}else{
				new[c] =org[i];
			}
		}
		//new[c] = org[i];
		c++;
	}
	MFREE(org);
	return new;
ERROR:
	return NULL;
}

int reverse(uint8_t* p,int len)
{
	int c, i, j;

	if(p[0] == '*'){
		return OK;
	}
	
	for (i = 0, j = len - 1; i < j; i++, j--)
	{
		c = p[i];
		p[i] = p[j];
		p[j] = c;
	}
	return OK;
}

int reverse_complement_sequence(char* p,int len)
{
	int c, i, j;
	
	for (i = 0, j = len - 1; i < j; i++, j--)
	{
		c = p[i];
		p[i] = p[j];
		p[j] = c;
	}
	for (i = 0; i <  len; i++){
		
		switch (p[i]) {
			case 0:
				
				p[i] = 3;
				break;
			case 1:
				p[i] = 2;
				break;
			case 2:
				p[i] = 1;
				break;
			case 3:
				p[i] = 0;
				break;
			default:
				break;
		}
		
	}
	return OK;
}

int write_sam_header(struct sam_bam_file* sb_file, FILE* out)
{
	int i;
	fprintf(out,"@HD\tVN:1.3\tSO:coordinate\n");
	for(i = 0; i < sb_file->header->n_targets;i++){
		fprintf(out,"@SQ\tSN:%s\tLN:%d\n", sb_file->header->target_name[i],(int)sb_file->header->target_len[i]);
	}
	return OK;
}

