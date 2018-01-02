
#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <unistd.h>
#include <dirent.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <ctype.h>
#include <string.h>
#include <inttypes.h>

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "tldevel.h"
#include "htsglue.h"

#define MAX_LINE 10000

#define OPT_READLEN 1
#define OPT_AT_START_ERROR 2
#define OPT_AT_END_ERROR 3
#define OPT_INDEL 4
#define OPT_WINDOW_START 5
#define OPT_SERVER 6
#define OPT_INDIR 7
#define OPT_OUT 8

#define RUNMODE_SIM 1
#define RUNMODE_EVAL 2

struct parameters {
        FILE* foutfile;
        char* outfile;
        char* bedfile;
        char* control_file;
        char* indir;
        int infiles;
        int quiet_flag;
        int num_threads;
        unsigned int sim_length;
        unsigned int seed;
        float sim_start_error;
        float sim_end_error;
        float sim_mismatch_prob;
        float sim_GG_error_addition;
        float random_model_prior;
        int sim_numseq;
        int window_start;
        int server;
        int help;
        int mode;
};

int eval_sambam(struct parameters* param);
int report_param(char* filename,FILE* fptr_out);

int simulate_reads_from_bed(faidx_t*  index ,struct parameters* param);

void usage();
char*  mutate_sequence(char* seq,struct parameters*  param, int* errors);
void convert_sequence(char* seq, int len);

int compare_genome_interval_mapping(struct genome_interval* a, struct genome_interval* b, int* res);
float* shuffle_arr_r(float* arr,int n);

char** get_files_names(char* directory, char* suffix, int* num_files);

static int str_cmp(const void * a, const void* b);



int main (int argc, char * argv[])
{

        struct stat st = {0};
        int i,c;
        struct parameters* param = NULL;
        faidx_t*  index = NULL;
        char* fasta_file = NULL;
        char* mode_string = NULL;
        //tlog.echo_build_config();

        MMALLOC(param, sizeof(struct parameters));
        param->mode = RUNMODE_SIM;
        param->help = 0;
        param->control_file = 0;
        param->bedfile = 0;
        param->outfile = 0;
        param->sim_length = 30;
        param->sim_start_error = 0.02f;
        param->sim_end_error = 0.03f;
        param->sim_mismatch_prob = 0.8f;
        param->sim_GG_error_addition = 0.0f;
        param->window_start = 0;
        param->sim_numseq = 1000000;
        param->server = 0;

        param->indir = NULL;
        while (1){
                static struct option long_options[] ={
                        {"read-len",required_argument,0,OPT_READLEN},
                        {"start-error",required_argument,0,OPT_AT_START_ERROR},
                        {"stop-error",required_argument,0,OPT_AT_END_ERROR},
                        {"mismatch-rate",required_argument,0,OPT_INDEL},
                        {"sw",required_argument,0,OPT_WINDOW_START},
                        {"server",0,0,OPT_SERVER},
                        {"in",required_argument,0,OPT_INDIR},
                        {"out",required_argument,0,OPT_OUT},
                        {"quiet",required_argument,0,'q'},
                        {"help",0,0,'h'},
                        {0, 0, 0, 0}
                };
		
                int option_index = 0;
                c = getopt_long_only (argc, argv,"qhn:",long_options, &option_index);
		
                if (c == -1){
                        break;
                }
		
                switch(c) {
                case 0:
                        break;
                case OPT_INDIR:
                        param->indir = optarg;
                        break;
                case OPT_SERVER:
                        param->server = 1;
                        break;
                case OPT_READLEN:
                        param->sim_length = atoi(optarg);
                        break;
                case OPT_AT_START_ERROR:
                        param->sim_start_error = atof(optarg);
                        break;
                case OPT_AT_END_ERROR:
                        param->sim_end_error = atof(optarg);
                        break;
                case OPT_INDEL:
                        param->sim_mismatch_prob = atof(optarg);
                        break;
                case OPT_WINDOW_START:
                        param->window_start = atoi(optarg);
                        break;
                case 'n':
                        if(optarg[0] == '-'){
                                fprintf(stderr,"option -n requires an argument\n");
                                exit(-1);
                        }
                        param->sim_numseq = atoi(optarg);
                        break;
                case OPT_OUT:
                        if(optarg[0] == '-'){
                                fprintf(stderr,"option -o requires an argument\n");
                                exit(-1);
                        }
                        param->outfile = optarg;
                        break;
                case 'c':
                        if(optarg[0] == '-'){
                                fprintf(stderr,"option -c requires an argument\n");
                                exit(-1);
                        }
                        param->control_file = optarg;
                        break;
                case 'h':
                        param->help = 1;
                        break;
                case '?':
                        exit(1);
                        break;
                default:
                        fprintf(stderr,"default\n\n\n\n");
                        abort ();
                }
        }
	
        if(param->help || (argc - optind) < 1){
                usage();
                MFREE(param);
                return EXIT_SUCCESS;
        }
        /* 0. there is an outfile */
        if(param->outfile == NULL){
                usage();
                ERROR_MSG("No outfile (-o) given.");
        }

        mode_string= argv[optind++];
        /* 1. mode  */

	

	
        param->mode = -1;
        if(strncmp(mode_string,"sim",3) == 0){
                param->mode = RUNMODE_SIM;
                LOG_MSG("Running in sim mode.");

                fasta_file = argv[optind++];
                /* 2. does fastafile exists and can it be opened...  */
                /* this will open the index and fail properly if not successfull */
                RUNP(index = get_faidx(fasta_file));
                LOG_MSG("Loaded \"%s\" faidx index.",fasta_file);
		
                param->bedfile = argv[optind++];
                /* 3. is bedfile there and readable */
                if(my_file_exists(param->bedfile) == 0){
                        usage();
                        ERROR_MSG("File \"%s\" does not exist.",param->bedfile);
                }
                RUN(simulate_reads_from_bed(index,param));
        }
        if(strncmp(mode_string,"eval",4) == 0){
                param->mode = RUNMODE_EVAL;
                LOG_MSG("Running in eval mode.");

                if(param->indir){
                        if(stat(param->indir, &st) == -1) {
                                ERROR_MSG("Directory %s does not  exists.", param->indir);
                        }
                        char** input_files = 0;
                        int num_input_files = 0;

                        /* Step 0 - find input files...  */
	
                        RUNP(input_files = get_files_names(param->indir,".bam",&num_input_files));
                        /* qsort array.  */
                        qsort(input_files,num_input_files,sizeof(char*), str_cmp); 
                        for(i = 0; i < num_input_files;i++){
                                LOG_MSG("Working on :%s",input_files[i]);
                                param->bedfile = input_files[i];
                                RUN(eval_sambam(param));
                        }

                        for(i = 0; i < num_input_files;i++){
                                MFREE(input_files[i]);
                        }
                        MFREE(input_files);
                        input_files = NULL;
                        num_input_files= 0;
                        RUNP(input_files = get_files_names(param->indir,".sam",&num_input_files));
                        for(i = 0; i < num_input_files;i++){
                                param->bedfile = input_files[i];
                                RUN(eval_sambam(param));
                        }

                        for(i = 0; i < num_input_files;i++){
                                MFREE(input_files[i]);
                        }
                        MFREE(input_files);
                }else{
                        param->bedfile = argv[optind++];
                        if(my_file_exists(param->bedfile) == 0){
                                ERROR_MSG("File \"%s\" does not exist.",param->bedfile);
                        }
                        RUN(eval_sambam(param));
                }
        }	
        if(param->mode  == -1){
                ERROR_MSG("Mode \"%s\" not found", mode_string);
        }
	
        if(index){
                fai_destroy(index);
        }
        if(param){
                MFREE(param);
        }
	
        return EXIT_SUCCESS;
ERROR:
        if(index){
                fai_destroy(index);
        }
        if(param){
                MFREE(param);
        }
        return EXIT_FAILURE;
}

void usage()
{
        fprintf(stdout, "Usage 1: bedsim sim  <tomedb> <in.bed> -o <suffix> \n\n");
        fprintf(stdout, "Options:\n");
        fprintf(stdout,"%15s%10s%7s%-30s\n","-read-len","INT" ,"", "Length of simulated reads [30].");
        fprintf(stdout,"%15s%10s%7s%-30s\n","-start-error","FLT" ,"", "Error rate at 5' end [0.02].");
        fprintf(stdout,"%15s%10s%7s%-30s\n","-stop-error","FLT" ,"", "Error rate at 3' end [0.03].");
        fprintf(stdout,"%15s%10s%7s%-30s\n","-mismatch-rate","FLT" ,"", "Fraction of simulated mismatches [0.8].");
        fprintf(stdout,"%15s%10s%7s%-30s\n","-n","INT" ,"", "Number of reads to simulate [1000000].");
        fprintf(stdout,"%15s%10s%7s%-30s\n","-o","STR" ,"", "Output file name.");
        fprintf(stdout, "\n");
	
        fprintf(stdout, "Usage 2: bedsim eval <XXX.bam>  \n\n");
        fprintf(stdout, "Options:\n");
        fprintf(stdout,"%15s%10s%7s%-30s\n","-in","STR" ,"", "Input directory.");
        fprintf(stdout,"%15s%10s%7s%-30s\n","-out","STR" ,"", "Output file..");
        fprintf(stdout, "\n");	

	
}

int simulate_reads_from_bed(faidx_t*  index ,struct parameters* param)
{	
        struct genome_interval* g_int = NULL;
        struct seq_info* si = NULL;
        FILE* file = 0;
//	FILE* local_out_file = 0;
        FILE* outfile = 0;
        char line[MAX_LINE];
//	char seq_name[MAX_LINE];
        char* seq = NULL;
        int num_sequences = param->sim_numseq;
        int column = 0;
        int i,j,tmp;
        int start;
        int stop;
//	unsigned int chr_start;
        int strand = 0;
        int c = 0;
//	int hit = 0;
        int rank = 0;
        float y = 0.0f;
        float total = 0.0f;
        int errors;
	
        float* rank_arr = NULL;
	
        unsigned int sim_start;
        int seq_num = 1;
//	char* local_outfile;
//	char* local_control;
        RUNP(g_int = init_genome_interval(0,0,0));	
        RUNP(si =  make_si_info_from_fai(index));

        RUNP(file = fopen(param->bedfile , "r" ));
        RUNP(outfile = fopen(param->outfile , "w" ));
        while(fgets(line, MAX_LINE, file)){
                if(line[0] != '#'){
                        total++;
                }
        }
        rewind(file);

        MMALLOC(rank_arr, sizeof(float) * total);
		
        for(i = 0; i < total;i++){
                rank_arr[i] = i+1;
        }
        RUNP(rank_arr = shuffle_arr_r(rank_arr ,total));
        y = 0.0;
        for(i = 0; i < total;i++){
                rank_arr[i] = (float)rank_arr[i] / (  pow((float) rank_arr[i], 1.2f + ((float) ((int)rank_arr[i] % (int)(total+1)) / total)));
                y+= rank_arr[i] ;
        }
	
        c = 0;
        for(i = 0; i < total;i++){
                rank_arr[i]  = (int)(rank_arr[i]  / y * (float)num_sequences);
//		fprintf(stdout,"%d\t%f\n",i,rank_arr[i]);
                c += rank_arr[i] ;
        }
	
        while(c < num_sequences){
                rank_arr[random_int_zero_to_x(total-1)]++;
                c++;
        }
	
        for(c = 0;c < 1;c++){
                rank = 0;
		
                while(fgets(line, MAX_LINE, file)){
                        if(line[0] != '#'){
                                column = 1; //<QNAME>
                                tmp = 0;
                                //hit = 0;
                                for(j = 0;j < MAX_LINE;j++){
                                        tmp++;
                                        if(isspace((int)line[j])){
                                                break;
                                        }
                                }
				
                                for(j = 0;j < FIELD_BUFFER_LEN;j++){
					
                                        if(isspace((int)line[j])){
                                                g_int->chromosome[j] =0; 
                                                break;
                                        }
                                        g_int->chromosome[j] = line[j];
			        
                                }
                                for(i = 0; i < MAX_LINE;i++){
                                        if(line[i] == '\n'){
                                                break;
                                        }
                                        if(isspace((int)line[i])){
                                                column++;
                                                switch(column){
                                                case 2: // <FLAG>
                                                        g_int->start = atoi(line + i + 1);
                                                        break;
                                                case 3: // <RNAME>
                                                        g_int->stop = atoi(line + i + 1);;
                                                        break;
                                                case 6: // <RNAME>
                                                        if(*(line+i+1) == '+'){
                                                                g_int->strand = 0;
                                                        }else{
                                                                g_int->strand = 1;
                                                        }
                                                        break;
                                                }
                                        }
                                }
                                //g_int->stop = param->sim_length + g_int->start;
                                for(j = 0; j < (int)(rank_arr[rank] );j++){
                                        start = g_int->start;
                                        stop = g_int->stop;
                                        if(!param->window_start){
						
                                                sim_start = random_int_zero_to_x(g_int->stop - g_int->start);
                                                sim_start += g_int->start;
                                                g_int->start = sim_start;
                                                g_int->stop = sim_start +param->sim_length;
                                        }else{
                                                sim_start = random_int_zero_to_x(2*param->window_start);
						
                                                if(!strand){
                                                        sim_start = sim_start - param->window_start;
                                                        sim_start += g_int->start;
                                                        g_int->start = sim_start;
                                                        g_int->stop = sim_start +param->sim_length;
                                                }else{
                                                        sim_start =  sim_start - param->window_start;
                                                        sim_start += g_int->stop;
                                                        g_int->start = sim_start;
                                                        g_int->stop = sim_start +param->sim_length;
                                                }
                                        }

                                        DPRINTF3("lookign for: %s:%d-%d %c",g_int->chromosome,g_int->start,g_int->stop, "+-"[g_int->strand]);
                                        seq = get_sequence(index,g_int);
                                        errors = 0;
                                        RUNP(seq = mutate_sequence(seq,param ,&errors));

                                        fprintf(outfile,">ORG%d_%s_%d_%d_%c_%d\n%s\n",seq_num, g_int->chromosome,g_int->start,g_int->stop, "+-"[g_int->strand],errors,seq);
                                        seq_num++;
                                        MFREE(seq);
                                        g_int->start = start;
                                        g_int->stop = stop;
                                }
                                rank = rank + 1;
                        }
                }
                rank = 0;	
		
        }
        fclose(file);
        fclose(outfile);
        free(rank_arr);
        free_sequence_info(si);
        free_genome_interval(g_int);
        return OK;
ERROR:
        return FAIL;

}


void convert_sequence(char* seq, int len)
{
        char alphabet[] = "ACGTN";
        int i;
        for(i = 0; i < len;i++){
                seq[i] = alphabet[(int)seq[i]];
        }
        seq[len] = 0;
}




char* mutate_sequence(char* seq,struct parameters*  param,int* errors)
{
        int j,c;
	
        char alphabet[] = "ACGTN";
        float a,b, local_c;
        char subc;
        char* tmp_seq = NULL;
        char* ptr;
        int local_len = 64;

        MMALLOC(tmp_seq,sizeof(char) * local_len);
	

	
        b  = ((param->sim_end_error - ((float)param->sim_length *  param->sim_start_error)) / (float)(param->sim_length -1)) * -1.0f;

	
        a = param->sim_start_error - b;
        c = 0;
        for(j = 0; j < param->sim_length;j++){
                local_c =0.0;
                if(j){
                        if(seq[j-1] == 2 && seq[j] == 2){
                                local_c += param->sim_GG_error_addition;
                        }
                }
		
                if(random_float_zero_to_x(1.0)  <= (float)(j+1)* a +b + local_c){
                        *errors = *errors + 1;
                        if(random_float_zero_to_x(1.0)  <=  param->sim_mismatch_prob){
                                subc = alphabet[(int) random_int_zero_to_x(3)];
                                while(subc == seq[j]){
                                        subc = alphabet[random_int_zero_to_x(3)]; //subc = alphabet[random_int(4)];
                                }
                                if(c == local_len){
                                        local_len = local_len << 1;
                                        MREALLOC(tmp_seq,sizeof(char) * local_len);
                                }
		
                                tmp_seq[c] = subc;
                                c++;
//				fprintf(stderr,"MISMATCH\n");
                        }else{
                                if(random_float_zero_to_x(1.0)  <= 0.5f){
                                        if(c == local_len){
                                                local_len = local_len << 1;
                                                MREALLOC(tmp_seq,sizeof(char) * local_len);
                                        }
                                        tmp_seq[c] = alphabet[random_int_zero_to_x(3)];
                                        c++;
                                        if(c == local_len){
                                                local_len = local_len << 1;
                                                MREALLOC(tmp_seq,sizeof(char) * local_len);
                                        }
                                        tmp_seq[c] = seq[j];
                                        c++;
//						fprintf(stderr,"INSERT\n");
                                }else{
//						fprintf(stderr,"DELETION\n");
                                }
                        }
                }else{
                        if(c == local_len){
                                local_len = local_len << 1;
                                MREALLOC(tmp_seq,sizeof(char) * local_len);
                        }
                        tmp_seq[c] = seq[j];
                        c++;
                }
        }

        if(c == local_len){
                local_len = local_len << 1;
                MREALLOC(tmp_seq,sizeof(char) * local_len);
        }
		
        tmp_seq[c] = 0;
        MFREE(seq);
        return tmp_seq;
ERROR:
        return NULL;
}

float* shuffle_arr_r(float* arr,int n)
{
        int i,j,c;
        float tmp;
	

        for(c = 0; c < n*3;c++){
		
                i = random_int_zero_to_x(n-1);
                j = random_int_zero_to_x(n-1);
                if(i != j){
                        tmp = arr[j];
                        arr[j] = arr[i];
                        arr[i] = tmp;

                }			
        }
        return arr;
}




int eval_sambam(struct parameters* param)
{
        FILE* fptr_out = NULL;
        struct sam_bam_file* sb_file = NULL;
        struct sam_bam_entry* entry = NULL;
        struct genome_interval* g_int = NULL;
        struct genome_interval* g_org_int = NULL;
        struct results{
                int mapped[6];
                int mismapped[6];
                int error[6][6];
                int total;
                int total_mapped;
                int total_unmapped;
        } res;
        int buffer_size = 1000000;
        int i = 0;
        int j = 0;
        int norm_qual = 0;
        int comparison = 0;
        int real_error = 0;
        int sum = 0;
        char c = 0;
        float tmp = 0.0f;

       	

	
        /* Init res */
        res.total = 0;
        res.total_mapped = 0;
        res.total_unmapped = 0;

        for(i = 0; i < 6;i++){
                res.mapped[i] = 0;
                res.mismapped[i] = 0;
                for(j = 0; j < 6;j++){
                        res.error[i][j] = 0;
                }
        }

        if(my_file_exists(param->outfile) == 0){
                RUNP(fptr_out = fopen(param->outfile, "w"));
                fprintf(fptr_out,"Filename");
                fprintf(fptr_out,"\tDPAR_threads");
                fprintf(fptr_out,"\tDPAR_siter");
                fprintf(fptr_out,"\tDPAR_giter");
                fprintf(fptr_out,"\tDPAR_spseudo");
                fprintf(fptr_out,"\tDPAR_gspeudo");
                for(i = 0; i < 6;i++){
                        fprintf(fptr_out,"\tQ%d_correct",i*10);
                }
                for(i = 0; i < 6;i++){
                        fprintf(fptr_out,"\tQ%d_wrong",i*10);
                }
                for(j = 0; j < 6;j++){
                        for(i = 0; i < 6;i++){
                                fprintf(fptr_out,"\tQ%d_E%d",i*10,res.error[j][i]);
                        }
                }
                fprintf(fptr_out,"\n");
		


//	ERROR_MSG("File \"%s\" does not exist.",param->outfile);
        }else{
                RUNP(fptr_out = fopen(param->outfile, "a"));
        }
        fprintf(fptr_out,"%s",param->bedfile);
	
        RUN(report_param(param->bedfile,fptr_out));
        RUNP(g_int = init_genome_interval(0,0,0));	

        RUNP(g_org_int = init_genome_interval(0,0,0));	

	
        RUNP(sb_file = open_SAMBAMfile(param->bedfile,buffer_size,10,0,0));
        //echo_header(sb_file);
        while(1){
                RUN(read_SAMBAM_chunk(sb_file,1,0));
                //DPRINTF2("read:%d",sb_file->num_read );
                //RUN(read_sam_bam_chunk(infile,data->si, buffer,1,&numseq),"Read sam/bam chunk in thread %d failed", data->threadID);
                //if((status = read_sam_bam_chunk(infile,data->si, buffer,1,&numseq)) != kslOK)  exit(status);
                if(!sb_file->num_read){
                        break;
                }
                res.total += sb_file->num_read;
                for(i = 0; i < sb_file->num_read;i++){
                        entry = sb_file->buffer[i];
			
                        /*	 	int64_t* start;
                              int64_t* stop;
                              char* sequence;
                              char* name;
                              int qual;
                              int len;
                              int max_len;
                              int num_hits;
                              int max_num_hits;
                        */

			
                        if(entry->num_hits){
                                res.total_mapped++;
                                g_int->g_start = entry->start[0];
                                g_int->g_stop = entry->stop[0];
                                //fprintf(outfile,">ORG_%s_%d_%d_%c_%d\n%s\n",g_int->chromosome,g_int->start,g_int->stop, "+-"[g_int->strand],errors,seq);
                                //ORG_chr22_51061767_51061797_+_1
                                sscanf(entry->name,"%*[^_]_%[^_]_%d_%d_%c_%d",g_org_int->chromosome,&g_org_int->start,&g_org_int->stop, &c,&real_error);
                                if(c == (int) '+'){
                                        g_org_int->strand = 0;
                                }else{
                                        g_org_int->strand = 1;
                                }
                                RUN(internal_to_chr_start_stop_strand(sb_file->si,g_int));
//				RUN(get_chr_start_stop(sb_file->si  ,g_int,entry->start[0], entry->stop[0]));
//				fprintf(stdout,"%s\n%s\t%d\t%d\t%c\n", entry->name,g_int->chromosome,g_int->start,g_int->stop, "+-"[g_int->strand]);
//				fprintf(stdout,"%s\t%d\t%d\t%c\n", g_org_int->chromosome,g_org_int->start,g_org_int->stop, "+-"[g_org_int->strand]);
//				fprintf(stdout,"%d %f\n",entry->qual, (roundf((float) entry->qual /10.0)));
				
                                norm_qual =   entry->qual /10;
                                if(norm_qual > 5){
                                        norm_qual = 5;
                                }
                                if(real_error > 5){
                                        real_error = 5;
                                }
					      
                                RUN(compare_genome_interval_mapping(g_org_int,g_int,&comparison));
                                if(comparison < 5){ /* distance criteria to declare a match */
                                        comparison = 1;
                                }else{
                                        //	fprintf(stdout,"%s\t%d\n",entry->name,entry->num_hits );
                                        comparison = 0;
                                }
                                /* record stars */
                                if(comparison){
                                        res.mapped[norm_qual]++;
                                        res.error[real_error][norm_qual]++;
                                }else{
					
                                        res.mismapped[norm_qual]++;
                                }
				
//			      fprintf(stdout,"%d %f RES:%d  error:%d \n",entry->qual, (roundf((float) entry->qual /10.0)), comparison ,real_error);
				
                        }
                }
        }

        sum = 0;
        for(i = 0; i < 6;i++){
                fprintf(stdout,"%d sum:%d\n",res.mismapped[i],c);
                sum += res.mismapped[i];
        }
        res.total_unmapped = res.total - res.total_mapped;
       
        fprintf(stdout,"%d\tTotal\n", res.total);
        fprintf(stdout,"%d\tMapped\n", res.total_mapped);
        fprintf(stdout,"%d\tUn-mapped\n", res.total_unmapped);

        fprintf(stdout,"%f\tPercentage mapped\n", roundf((float)res.total_mapped / (float) res.total * 100.0f));
        fprintf(stdout,"%f\tPercentage mis-mapped %d\n", roundf((float)sum  / (float) res.total * 10000.0f)/100.0f,sum);
	
        j = 0;
        sum = 0;

        for(i = 1; i < 6;i++){
                //fprintf(stdout,"%d sum:%d\n",res.mismapped[i],c);
                sum += res.mapped[i];
                j += res.mismapped[i];
        }
        fprintf(stdout,"%d\tMapped\n", sum);
	
        fprintf(stdout,"%f\tPercentage mis-mapped %d\n", roundf((float)j  / (float) sum  * 10000.0f)/100.0f,j);
	      
	
        for(i = 0; i < 6;i++){
                fprintf(stdout,"%d\t",i*10);
        }
        fprintf(stdout,"mapping quality\n");

        sum =0;
        for(i = 0; i < 6;i++){
                sum+= res.mapped[i];
                fprintf(stdout,"%d\t",res.mapped[i]);
                fprintf(fptr_out,"\t%d",res.mapped[i]);
        }
        fprintf(stdout,"correctly mapped\n");
        for(i = 0; i < 6;i++){
                sum += res.mismapped[i];
                fprintf(stdout,"%d\t",res.mismapped[i]);
                fprintf(fptr_out,"\t%d",res.mismapped[i]);
        }
        fprintf(stdout,"in-correctly mapped total:%d\n",sum);
        for(i = 0; i < 6;i++){
                if(res.mismapped[i] == 0){
                        fprintf(stdout,"100\t");
                }else if(res.mismapped[i] + res.mapped[i] == 0){
                        fprintf(stdout,"-1\t");
                }else{
                        tmp = (float) res.mismapped[i] / (float)(res.mapped[i] + res.mismapped[i]);
                        fprintf(stdout,"%d\t", (int)(-10.0f * log10(tmp)));
                }
        }
        fprintf(stdout,"Q estimation\n");

        for(j = 0; j < 6;j++){
                for(i = 0; i < 6;i++){
                        fprintf(stdout,"%d\t",res.error[j][i]);
                        fprintf(fptr_out,"\t%d",res.error[j][i]);
                }
                fprintf(stdout,"%d error\n",j);
        }
        RUN(close_SAMBAMfile(sb_file));

        fprintf(fptr_out,"\n");
        if(fptr_out){
                fclose(fptr_out);
        }
        free_genome_interval(g_org_int);
        free_genome_interval(g_int);
        return OK;
ERROR:
        if(sb_file){
                RUN(close_SAMBAMfile(sb_file));
        }
        return FAIL;
}

int report_param(char* filename, FILE* fptr_out)
{
        char* ptr = NULL;
        char* end = NULL;
	
        int64_t develcode = 0;
	
        ptr = strstr(filename,"DCODE");
        if(ptr){
                ptr +=5;
                develcode = strtoull(ptr,&end,16);
                if(develcode == 0 && end == ptr){
                        LOG_MSG("could not read a number");
                }else if(develcode == INT64_MAX){
                        LOG_MSG("too long");
                }else if(*end){
                        LOG_MSG("remaining: %s", end);
                }
		
                fprintf(stdout,"%s\n%016jX\n",filename,develcode);	
                fprintf(stdout,"%d\tNumber of threads\n",(int) (  develcode >>(int64_t)(8*4) & 0xFF));
                fprintf(stdout,"%d\tsiter\n",(int) (  develcode >>(int64_t)(8*3) & 0xFF));
                fprintf(stdout,"%d\tgiter\n",(int) (  develcode >>(int64_t)(8*2) & 0xFF));
                fprintf(stdout,"%d\tspseudo\n",(int) (  develcode >>(int64_t)(8*1) & 0xFF));
                fprintf(stdout,"%d\tgpseudo\n",(int) (  develcode >>(int64_t)(8*0) & 0xFF));

                fprintf(fptr_out,"\t%d",(int) (  develcode >>(int64_t)(8*4) & 0xFF));
                fprintf(fptr_out,"\t%d",(int) (  develcode >>(int64_t)(8*3) & 0xFF));
                fprintf(fptr_out,"\t%d",(int) (  develcode >>(int64_t)(8*2) & 0xFF));
                fprintf(fptr_out,"\t%d",(int) (  develcode >>(int64_t)(8*1) & 0xFF));
                fprintf(fptr_out,"\t%d",(int) (  develcode >>(int64_t)(8*0) & 0xFF));
		
		
        }else{
                fprintf(fptr_out,"\t%d",-1);
                fprintf(fptr_out,"\t%d",-1);
                fprintf(fptr_out,"\t%d",-1);
                fprintf(fptr_out,"\t%d",-1);
                fprintf(fptr_out,"\t%d",-1);
        }

        return OK;
}


int compare_genome_interval_mapping(struct genome_interval* a, struct genome_interval* b, int* res)
{
        if(a->strand != b->strand){
                *res = -1;
                return OK;
        }

        if(strcmp(a->chromosome, b->chromosome) != 0){
                *res = -1;
                return OK;
        }
	
        *res = MACRO_MAX(abs(a->start - b->start),abs( a->stop - b->stop)); 
        return OK;
}


char** get_files_names(char* directory, char* suffix, int* num_files)
{
        char** list = NULL;
        char cwd[BUFFER_LEN];
        ASSERT(directory !=NULL," directory is NULL");
        ASSERT(suffix != NULL, "suffix is NULL");
        ASSERT(strlen(suffix) >1, "suffic is too short: %s", suffix);
	
        DIR *dp;
        struct dirent *entry;
        struct stat statbuf;
        int len;

        int malloc_len = 16;

        MMALLOC(list,sizeof(char*) * malloc_len); 
	
        *num_files = 0;
	
        //get current directory...
        RUNP(getcwd(cwd, sizeof(cwd)));
        //open directory in question to count number of files...
        RUNP(dp = opendir(directory));
	
        if(chdir(directory)){
                ERROR_MSG("chdir failed.");
        }
	
        while((entry = readdir(dp)) != NULL) {
                lstat(entry->d_name,&statbuf);
                if(S_ISDIR(statbuf.st_mode)) {
                        /* Found a directory, but ignore . and .. */
                        if(strcmp(".",entry->d_name) == 0 ||
                           strcmp("..",entry->d_name) == 0)
                                continue;
                }else{
                        if(!strcmp(suffix, entry->d_name+ (strlen(entry->d_name) - strlen(suffix)))){
                                if(*num_files == malloc_len){
                                        malloc_len = malloc_len <<2;


                                        MREALLOC(list, sizeof(char*) * malloc_len);
                                }
                                //fprintf(stdout,"%p List\n", list);
                                list[ *num_files] = NULL;
                                len = strlen(entry->d_name)+strlen(directory)  +2;
                                MMALLOC(list[ *num_files], sizeof(char) *len);
                                snprintf(list[*num_files],len,"%s/%s",directory,entry->d_name);
                                *num_files = *num_files + 1;
                        }
                }
        }
        if(closedir(dp)){
                ERROR_MSG("closedir failed.");
        }
        if(chdir(cwd)){
                ERROR_MSG("chdir failed.");
        }
        return list;
ERROR:
        return NULL;
}



static int str_cmp(const void * a, const void* b)
{
        return strcmp((const char*) a, (const char*)b);
}
