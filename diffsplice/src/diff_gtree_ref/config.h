/*    
 *    genome_tree/config.h		
 *    DiffSplice
 *
 *    Copyright (C) 2013 University of Kentucky and
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

/************************************************************************/
/* GENERAL SETTINGS                                                     */
/************************************************************************/

#define UNIX

#ifdef UNIX

#include <fstream>
#include <cstring>
#include <iostream>
#include <cmath>
#include <cstdlib> 
#include <ctime>
#include <cstdio>
#include <sstream>
#include <vector>
#include <queue>
#include <assert.h>
#include <utility>
#include <algorithm>

#else

#include <stdio.h>
#include <fstream>
#include <string>
#include <iostream>
#include <math.h>
#include <stdlib.h> 
#include <time.h>
#include <sstream>
#include <vector>
#include <queue>
#include <assert.h>

#endif

using namespace std;

//IMPORTANT: BAM files are 0-based coordinates, so need to add 1 to all junction and exonic sorted files in order to match up with annotation (if provided)
const int READ_COORDINATE_OFFSET = 1; //the position that will be used throughout the program will be the position read in from fragment files + this offest

#define FILTER_JUNCTION
//#define NORMALIZE_COVERAGE

// for debug
//#define DEBUG_OUTPUT
extern bool OUTPUT_PER_SAMPLE_EXPR; //whether to give per sample expression in gTree output (GTree.txt), not recommended for large dataset; change to automatically disable if # of samples is > 20
const int debug_total_sample_num = 1; //may not be used, depending on the purpose, check the usage
const int debug_sample_index_offset = 0; //the index of the sample to be debugged - 1

//#define JUNCTIONONLY
//#define DO_MERGEALTERSITES
//#define COUTSTEPS

extern int SUPPORT_VECTOR_SIZE;
extern int SUPPORT_VECTOR_SIZE_WHOLE_DATA; //number of samples for the whole data, different with SUPPORT_VECTOR_SIZE only when debugging with a few samples
// extern double thresh_junction_filter_max_read_support;
// extern double thresh_junction_filter_mean_read_support;
// extern double thresh_junction_filter_num_samples_presence;
extern string junctionfilter_mode;
extern double thresh_junctionfilter_groupwise;
extern double thresh_num_present_samples;
extern double thresh_present_support;
extern string junction_annotation_as_white_list;

extern bool COUNT_MAJOR_PATHS_ONLY;
extern double coverageThreshold_exon; //throw-away := <= coverageThreshold * SUPPORT_VECTOR_SIZE
extern double coverageThreshold_intron; //throw-away := <= coverageThreshold * SUPPORT_VECTOR_SIZE

// effective index starts from 1 for these group info vectors, to keep consistent with test
extern int global_total_group_num; //number of subtype groups
extern vector<int> global_group_size; //number of samples in each group, group index from 1
extern vector<int> global_group_base_index; //base index for samples in each group, i.e. sample index of the first sample in each group - 1, group index from 1, sample index from 1
extern vector<string> global_group_name; //name of each group

/************************************************************************/
/* SETTINGS		                                                        */
/************************************************************************/

//for the input of splice junctions
const int min_junction_intron_length = 5; //consider junctions whose end_position - start_position - 1 < min_intron_length as small deletion and do not include them in splice graph

const long MIN_EXON_LENGTH = 5;
const long MIN_ALTER_SPLICE_SITE_LENGTH = 1; //minimum exonic length for alternative splice sites. need a maximum??
const long MAX_ALTER_SPLICE_SITE_LENGTH = 200; //maximum exonic length for alternative splice sites. longer one should 

const long MAX_NOCOVERAGE_LENGTH = 1; //maximal length for an exon to have no coverage
const double MAX_NOCOVERAGE_THRESH_LOW_COV_REGION = 5;
const long MAX_NOCOVERAGE_LENGTH_LOW_COV_REGION = 50;

//filter exon in gene if its coverage is less than 0.05 of mean gene coverage in more than 0.95 samples
const double coverageThreshold_GeneExon = 0; //percentage of mean gene coverage

const double MIN_COVERAGE_FILTER_FRAGMENT = 0; //minimum coverage for fragments when filtering standalone fragments.

const int MAXJUNCNUMINDEPPATH = 100; //maximum number of junctions when separating dependent paths. if #junctions > this number, do not enumerate paths
const int MAX_NUM_ENUMERATED_PATH = 50; //maximum number of dependent paths when enumerating, this is set in the purpose to replace the above one. try to enumerate first, stop and keep NOTHING if exceeding this limit
const int MAX_NUM_PATH_FOR_NO_DECOMPOSITION = 5; //4/7/2013, do not try to decompose if the number of paths is <= this value

const int COVERAGE_CHANGE_WINDOW = 50;
const double COVERAGE_CHANGE_THRESH_ALTER_SITE = 5.0;
const double COVERAGE_CHANGE_THRESH_EXON = 10.0;

const int thresh_path_expressed = 5; //for the selection of major paths


const long MAX_CHR_LENGTH = 1000000000; //maximal chromosome length


extern vector <double> DatasetNormalizationRatio;
extern string resultPath;
extern string chromosomeName;


//other running-time settings
extern bool setting_junc_anno_provided;


#endif

