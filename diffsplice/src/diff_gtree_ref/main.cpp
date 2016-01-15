#include "config.h"
#include "input_data.h"
#include "decomposition.h"
#include "estimation.h"
#include "output.h"
#include "splice_graph.h"

//global parameters - general settings from config file
int SUPPORT_VECTOR_SIZE = 0;
int SUPPORT_VECTOR_SIZE_WHOLE_DATA = 0;
// double thresh_junction_filter_max_read_support = 0;
// double thresh_junction_filter_mean_read_support = 0;
// double thresh_junction_filter_num_samples_presence = 0;
string junctionfilter_mode = "none";
double thresh_junctionfilter_groupwise = 0.0;
double thresh_num_present_samples = 0.0;
double thresh_present_support = 0.0;
string junction_annotation_as_white_list = "yes";

bool COUNT_MAJOR_PATHS_ONLY = false;
double coverageThreshold_exon = 0.0; //throw-away := <= coverageThreshold * SUPPORT_VECTOR_SIZE
double coverageThreshold_intron = 1.0; //throw-away := <= coverageThreshold * SUPPORT_VECTOR_SIZE

int global_total_group_num; //number of subtype groups
vector<int> global_group_size; //number of samples in each group, group index from 1
vector<int> global_group_base_index; //base index for samples in each group, i.e. sample index of the first sample in each group - 1, group index from 1, sample index from 1
vector<string> global_group_name; //name of each group

//other global parameters
vector <double> DatasetNormalizationRatio; //per sample expression will be multiplied by these numbers for normalization
string resultPath;
string chromosomeName;

bool setting_junc_anno_provided = false;
bool OUTPUT_PER_SAMPLE_EXPR = true;

//parameters for this cpp only
const long default_chr_length = 300000000; 
string inputPath_frag;


void initialization()
{

// 	//read statFile
// 	ifstream statfile;
// #ifdef UNIX
// 	filename = inputPath + "stat.txt";
// 	statfile.open(filename.c_str());
// #else
// 	statfile.open("tmp\\stat.txt");
// #endif
// 
// 	long tmp = -1, tmp_small = MAX_CHR_LENGTH, tmp_large = -1;
// 	while (statfile >> tmp)
// 	{
// 		if (tmp < tmp_small)
// 			tmp_small = tmp;
// 
// 		statfile >> tmp;
// 		if (tmp > tmp_large)
// 			tmp_large = tmp;
// 	}
// 	CHROMOSOME_START = tmp_small;
// 	CHROMOSOME_END = tmp_large;
// 	statfile.close();

	if (SUPPORT_VECTOR_SIZE > 20) {
		OUTPUT_PER_SAMPLE_EXPR = false;
	}
	else {
		OUTPUT_PER_SAMPLE_EXPR = true;
	}

	srand ( time(NULL) );
	
	//read dataset stat file
	DatasetNormalizationRatio.assign(SUPPORT_VECTOR_SIZE, 1.0);

	ifstream file_normalization_ratio;
	string filename = inputPath_frag + "/Normalization_Ratio.txt";
	file_normalization_ratio.open(filename.c_str());
	if (file_normalization_ratio.is_open())
	{
		for (int tmpCnt = 0; tmpCnt < SUPPORT_VECTOR_SIZE; ++tmpCnt)
		{
			file_normalization_ratio >> DatasetNormalizationRatio[tmpCnt];
		}
		file_normalization_ratio.close();
	}

	initialize_outputfiles();

	return;
}



void cleanAll()
{
	DatasetNormalizationRatio.clear(); //free_vector(DatasetNormalizationRatio);

	close_outputfiles();

	return;
}


void input_config(string filename)
{
	ifstream config_file;
	config_file.open(filename.c_str());

	string parameter, info;

	if (config_file.is_open() == true)
	{
		while (config_file >> parameter)
		{
			if (parameter.compare("SUPPORT_VECTOR_SIZE") == 0)
				config_file >> SUPPORT_VECTOR_SIZE;
			// 			else if (parameter.compare("thresh_junction_filter_max_read_support") == 0)
			// 				config_file >> thresh_junction_filter_max_read_support;
			// 			else if (parameter.compare("thresh_junction_filter_mean_read_support") == 0)
			// 				config_file >> thresh_junction_filter_mean_read_support;
			// 			else if (parameter.compare("thresh_junction_filter_num_samples_presence") == 0)
			// 				config_file >> thresh_junction_filter_num_samples_presence;
			else if (parameter.compare("junctionfilter_mode") == 0)
				config_file >> junctionfilter_mode;
			else if (parameter.compare("thresh_junctionfilter_groupwise") == 0)
				config_file >> thresh_junctionfilter_groupwise;
			else if (parameter.compare("thresh_junctionfilter_num_present_samples") == 0)
				config_file >> thresh_num_present_samples;
			else if (parameter.compare("thresh_junctionfilter_present_support") == 0)
				config_file >> thresh_present_support;
			else if (parameter.compare("junction_annotation_as_white_list") == 0)
				config_file >> junction_annotation_as_white_list;
			else if (parameter.compare("COUNT_MAJOR_PATHS_ONLY") == 0)
			{
				config_file >> info;
				if (info.compare("yes") == 0)
					COUNT_MAJOR_PATHS_ONLY = true;
				else
					COUNT_MAJOR_PATHS_ONLY = false;
			}
			else if (parameter.compare("coverageThreshold_exon") == 0)
				config_file >> coverageThreshold_exon;
			else if (parameter.compare("coverageThreshold_intron") == 0)
				config_file >> coverageThreshold_intron;

			else if (parameter.compare("global_total_group_num") == 0)
				config_file >> global_total_group_num;
			else if (parameter.compare("global_group_size") == 0)
			{
				global_group_size.assign((size_t)(1+global_total_group_num), 0);
				for (int tmpCnt = 1; tmpCnt <= global_total_group_num; ++tmpCnt)
					config_file >> global_group_size[tmpCnt];
			}
			else if (parameter.compare("global_group_base_index") == 0)
			{
				global_group_base_index.assign((size_t)(1+global_total_group_num), 0);
				for (int tmpCnt = 1; tmpCnt <= global_total_group_num; ++tmpCnt)
					config_file >> global_group_base_index[tmpCnt];
			}
			else if (parameter.compare("group_name") == 0)
			{
				global_group_name.assign((size_t)(1+global_total_group_num), "unspecified");
				for (int tmpCnt = 1; tmpCnt <= global_total_group_num; ++tmpCnt)
					config_file >> global_group_name[tmpCnt];
			}

			getline(config_file, info);
		}
	} 
	else
	{
		cout << "Error: fail to open the config file for gtree" << endl;
		exit(1);
	}

	config_file.close();

	if (SUPPORT_VECTOR_SIZE <= 0)
	{
		cout << "Error: incomplete config file for gtree" << endl;
		exit(1);
	}

	if (junctionfilter_mode != "group" && junctionfilter_mode != "nogroup" && junctionfilter_mode != "none") {
		junctionfilter_mode = "none";
		cout << "Warning: junction filter mode not recognized, no filter will be applied" << endl;
	}

#ifdef UNIX
	SUPPORT_VECTOR_SIZE_WHOLE_DATA = SUPPORT_VECTOR_SIZE;
#else
	SUPPORT_VECTOR_SIZE_WHOLE_DATA = SUPPORT_VECTOR_SIZE;
	SUPPORT_VECTOR_SIZE = debug_total_sample_num;
#endif

	return;
}



int main(int argc, char *argv[])
{
	string config_filename, inputPath_annotation, inputPath_exon_annotation, inputPath_junction_annotation;
	long chromosome_length = 0;

#ifdef UNIX
	if (argc == 7)
	{
		inputPath_frag = argv[1];
		inputPath_annotation = argv[2];
		if (inputPath_annotation.compare("no_anno") == 0) {
			inputPath_exon_annotation = inputPath_annotation;
			inputPath_junction_annotation = inputPath_annotation;
		}
		else {
			inputPath_exon_annotation = inputPath_annotation + "/annotation_exon/";
			inputPath_junction_annotation = inputPath_annotation + "/annotation_junction/";
		}
		resultPath = argv[3];
		chromosomeName = argv[4];
		config_filename = argv[5];
		chromosome_length = atol(argv[6]);

		if (chromosome_length <= 0)
			chromosome_length = default_chr_length;
		else
			chromosome_length += 100;

		cout << chromosomeName << "\t";
	}	
	else
	{
		cout << argv[0] << "\t<inputPath_frag>\t<inputPath_annotation>\t<resultPath>\t<chromosomeName>\t<config_file>\t<chromosome_length>\n";
		exit(1);
	}
#else
	inputPath_frag = "frag\\";
	inputPath_exon_annotation = "annotation_exon\\";
	inputPath_junction_annotation = "annotation_junction\\";
	resultPath = "result\\";
	chromosomeName = "chr11";
	config_filename = "config_gtree";
	chromosome_length = 135006516;

	cout << chromosomeName << "\t";
#endif	

	input_config(config_filename);
	initialization();	

	/************************************************************************/
	/* INPUT COVERAGE                                                       */
	/************************************************************************/
	RangeJunctionList *list_fragments = new RangeJunctionList;
	list_fragments->rangeLow = chromosome_length;
	list_fragments->rangeHigh = 0;
	if (!input_data_ref(list_fragments, inputPath_frag, inputPath_junction_annotation, inputPath_exon_annotation, chromosome_length)) {
		input_data(list_fragments, inputPath_frag, inputPath_junction_annotation, chromosome_length);
	}


	/************************************************************************/
	/* CONSTRUCTION                                                         */
	/************************************************************************/

	//build tree
	cout << "Decomposing the splice graph... " << flush;
	GTvertex *gTree = decomposition(list_fragments);


//	/************************************************************************/
//	/* ESTIMATION                                                           */
//	/************************************************************************/

	//recount
	cout << "Estimating the transcript abundance... " << flush;
	countGTree(gTree);

	cout << "Writing output files... " << flush;
	output(gTree);

	cout << "Finished" << flush;
	cleanAll();

	//	cin >> filename;
	return 0;
}


