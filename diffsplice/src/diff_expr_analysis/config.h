/*    
 *    config.h		
 *    DiffSplice
 *
 *    Copyright (C) 2012 University of Kentucky and
 *                       Yin Hu
 *
 *    Authors: Yin Hu
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

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

extern double THRESHOLD_MIN_EXPR_COVERAGE; //all groups must have mean expression larger than this threshold to make differential transcription
extern double false_discovery_rate;
extern double thresh_foldchange_up;
extern double thresh_foldchange_down;
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

const double MAX_NUMBER = 20000000000;
const long DEFAULT_TESTASM_NUM = 20000;

const double THRESH_BAD_ESTIMATION = MAX_NUMBER; //an ASM will not be considered if more than half samples have estimation error ratio larger than this cutoff


#endif


