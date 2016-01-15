
#ifndef CONFIGURATION
#define CONFIGURATION

#define UNIX

#ifdef UNIX

#include <fstream>
#include <cstring>
#include <iostream>
#include <cmath>
#include <cstdlib> 
#include <ctime>
#include <cstdio>
#include <vector>
#include <assert.h>
#include <sstream>

#else

#include <stdio.h>
#include <fstream>
#include <string>
#include <iostream>
#include <math.h>
#include <stdlib.h> 
#include <time.h>
#include <vector>
#include <assert.h>
#include <sstream>

#endif

using namespace std;



//#define OUTPUT_PAIRWISE_JSD


/************************************************************************/
/* GENERAL SETTINGS                                                     */
/************************************************************************/

extern bool COUNT_INTRON_RETENTION; //whether to take intron retention into account

extern double THRESHOLD_MIN_ASM_COVERAGE; //all groups must have mean expression larger than this threshold to make differential transcription
extern double false_discovery_rate;
extern double thresh_JSD; //minimum sqrt of max pairwise JSD for an ASM to be called significant, default 0.1
extern double thresh_stat_d; //minimum stat_difftrans for an ASM to be called significant, considered as the square of the mean difference-variance ratio, default as 3^2=9

//general configuration
extern int global_total_sample_num; //total number of samples
extern int global_total_group_num; //number of subtype groups
extern int* global_group_size; //number of samples in each group, group index from 1
extern int* global_group_base_index; //base index for samples in each group, i.e. sample index of the first sample in each group - 1, group index from 1, sample index from 1
extern string* global_group_name; //name of each group
extern int global_num_sample_mask; //number of samples to be omitted when read in ASMs, count from first one

//maximum permutation count
const long MAX_PERMUTATION_NUMBER = 1000; //MAX_NUMBER if not limited

/************************************************************************/
/* FIXED SETTINGS	                                                    */
/************************************************************************/

const double MAX_NUMBER = 2000000000;
const long DEFAULT_TESTASM_NUM = 20000;

const double THRESH_BAD_ESTIMATION = MAX_NUMBER; //an ASM will not be considered if more than half samples have estimation error ratio larger than this cutoff

//this one does not make sense and has been disabled
const double THRESHOLD_MIN_ASMPATH_COVERAGE = -1; //all groups in all time points must have mean expression larger than this threshold to make differential transcription


#endif


