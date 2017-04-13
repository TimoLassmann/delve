
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <pthread.h>


#if HAVE_CONFIG_H
#include <config.h>
#endif


#include "tldevel.h"
#include "thr_pool.h"
#include "dyn_prog.h"


#ifdef DEBUG
#define MAXNUMQUERY 1000
#else
#define MAXNUMQUERY 1000000
#endif

//#define LIST_STORE_SIZE 10

//#define NCODE  4
//#define BORDER_REGION 2
//#define PSEUDO_CONSTANT 1000
//#define MAX_LOCAL_POSTERIORS 100000
//#define MAX_ERROR_LIMIT_FOR_TRAINING 100
//#define MAX_LENGTH_LIMIT_FOR_TRAINING 100

//#define ALIGNMENT_FLANKING_LENGTH 10
//#define MAX_LINE 10000


#ifdef HUGE_VAL
#define SCALEINFTY HUGE_VAL
#endif




