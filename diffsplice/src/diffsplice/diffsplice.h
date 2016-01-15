/*    
 *    diffsplice.h		
 *    DiffSplice
 *
 *    Copyright (C) 2014 Yin Hu
 *
 *    Authors: Yin Hu
 *
 *    All rights reserved.
 */


#define  UNIX

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
#include <unistd.h>
#include <getopt.h>

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

#endif


using namespace std;



//bool SUMMARIZE = false;
//bool DELETE_OLD_RESULT = true;


vector <string> chrName;
vector <string> chrLength;


void inputChrName(string filename);
bool exist_frag(string fname_check);

void parse_inputfiles(string filename); //extract the data files from the given file.
void parse_settings(string filename); //extract the settings from the given file.
string itostr(int); //convert integer to string
void write_configfile(string target_path);


class sample
{
public:
	string name;
	int groupID;
	int sampleID_orig; // id as in the original order specified in datafile.cfg
	int sampleID_dfs; // id used in diffsplice process, note that this may be different from the one above if the original order does not put samples of the same group together
	string datafilename;

	sample(string, int, int, string);
	~sample();
};

class sampleGroup
{
public:
	string name;
	int groupID;
	vector <sample> samples;

	sampleGroup(string, int);
	~sampleGroup();
};


int totalSampleNum = 0;
vector <sampleGroup> allGroups;

// config options for gtree
string junctionfilter_mode = "none";
double thresh_junctionfilter_groupwise = 0;
double thresh_junctionfilter_num_present_samples = 0;
double thresh_junctionfilter_present_support = 0;
string junction_annotation_as_white_list = "yes";
bool COUNT_MAJOR_PATHS_ONLY = false;
double coverageThreshold_exon = 0; //throw-away := <= coverageThreshold * SUPPORT_VECTOR_SIZE
double coverageThreshold_intron = 1; //throw-away := <= coverageThreshold * SUPPORT_VECTOR_SIZE

// config options for differential tests
bool COUNT_INTRON_RETENTION = true;
double config_false_discovery_rate = 0.01;
double config_thresh_foldchange_up = 2.;
double config_thresh_foldchange_down = 0.5;
double config_thresh_sqrtJSD = 0.1;
double thresh_stat_d = 3;
double THRESHOLD_MIN_EXPR_COVERAGE = 10;//at least one group must have mean expression larger than this threshold to make differential expression;
double THRESHOLD_MIN_ASM_COVERAGE = 5; //all groups must have mean expression larger than this threshold to make differential transcription


int flag_verbose = 0;

bool PREPROCESS = false;
bool SEP_TREE = false;
bool SELECT_SIGNIFICANCE = false;

static struct option long_options[] =
{
		{"verbose", no_argument,		&flag_verbose, 1},
		{"help",	no_argument,		0, 'h'},

		{"mode",	required_argument,	0, 'm'},
		{"chrlist",	required_argument,	0, 'c'},
		{"thread",	required_argument,	0, 'p'},
		{"outdir",	required_argument,	0, 'o'},
		{"anno_splice",	required_argument,	0, 'a'},
		{"runname",	required_argument,	0, 'r'},
		{"settings", required_argument, 0, 's'},
		{"gaf", required_argument, 0, 'g'},

		{0, 0, 0, 0}
};

enum RUNMODE {run_nothing, run_default, run_full, run_rerun};
