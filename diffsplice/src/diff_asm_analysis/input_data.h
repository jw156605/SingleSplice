#include "config.h"
#include "asm_def.h"
#include "calc_function.h"


#ifndef INPUT_DATA
#define INPUT_DATA


void input_config(string filename);
void input_ASMs(string resultPath, string outputPath);
void input_annotation_gaf(string filename);
bool override_sample_groups(string filename);
void input_ASMs_regroup(string resultPath, string outputPath);


class group {
public:
	string group_name;
	vector<string> sample_name;
	vector<long> sample_idx; // sample index in original datafile sheet

	group() {};
	group(string name) : group_name(name) {};
};

extern vector<group> global_sample_groups;

#endif

