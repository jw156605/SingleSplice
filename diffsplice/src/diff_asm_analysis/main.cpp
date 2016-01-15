
#include "config.h"
#include "asm_def.h"
#include "input_data.h"
#include "calc_function.h"
#include "difftrans_test.h"
#include "write_output.h"

vector <ASM*> sortList_ASM; //fragment list to be sorted
vector<annotation_chr_genelist*> list_anno_chr_genelist;


//general configuration
int global_total_sample_num = 0; //total number of samples
int global_total_group_num = 0; //number of subtype groups
int* global_group_size; //number of samples in each group
int* global_group_base_index; //base index for samples in each group, i.e. sample index of the first sample in each group - 1
string* global_group_name; //name of each group
int global_num_sample_mask = 0; //number of samples to be omitted when read in ASMs, count from first one

bool COUNT_INTRON_RETENTION = false; //whether to take intron retention into account
double false_discovery_rate = 0.01;
double thresh_JSD = 0.1;
double thresh_stat_d = 9;
double THRESHOLD_MIN_ASM_COVERAGE = 0; //all groups must have mean expression larger than this threshold to make differential transcription



void initialization()
{
	mergeSort_initialization(sortList_ASM.size()-1);
		
	return;
}

void clearall()
{
	mergeSort_cleanup();

	for (unsigned long geneLoopCnt = 1; geneLoopCnt < sortList_ASM.size(); ++geneLoopCnt)
		delete sortList_ASM[geneLoopCnt];
	sortList_ASM.clear();
	
	return;
}



int main(int argc, char *argv[])
{
	string resultPath, resultFileSuffix, config_filename, outputPath, samplegroup_filename;

#ifdef UNIX
	if (argc < 4 || argc > 6)
	{
		cout << argv[0] << "\t<result_path>\t<config_file>\t<gene_annotation_file>\t[<output_dir>\t<sample_group_file>]" << endl;
		exit(1);
	}
	resultPath = argv[1];
	config_filename = argv[2];

	if (argc >= 5){
		outputPath = argv[4];
		if (argc == 6) {
			samplegroup_filename = argv[5];
		}

		string comd;
		comd = "mkdir -p " + outputPath + "/stat/asm"; system(comd.c_str());
		comd = "mkdir -p " + outputPath + "/stat/temp"; system(comd.c_str());
	}
	
#else
	resultPath = "result\\"; // "E:\\Work\\Cluster_BRCA\\diff_trans\\dfs_result\\";
	config_filename = "config_testtrans"; // "config_testtrans_tumor";
#endif

	if (outputPath.empty()){
		outputPath = resultPath;
	}

	input_config(config_filename);

	/* initialize random seed: */
	srand ( time(NULL) );

	sortList_ASM.reserve(DEFAULT_TESTASM_NUM);
	sortList_ASM.push_back(NULL);

	if (!samplegroup_filename.empty()){
		override_sample_groups(samplegroup_filename);
		global_num_sample_mask = 0;
		input_ASMs_regroup(resultPath, outputPath);
	}
	else{
		input_ASMs(resultPath, outputPath);
	}


#ifdef UNIX
	input_annotation_gaf(argv[3]);
#else
	input_annotation_gaf("gaf.txt");
#endif

	initialization();
		
	//test for all individuals
	resultFileSuffix = "multi_class";
	difftrans_test_multi_class(outputPath, resultFileSuffix);

	write_output(outputPath);

	clearall();

	return 0;
}

