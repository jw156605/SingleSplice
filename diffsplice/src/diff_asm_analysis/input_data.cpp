#include "input_data.h"

const long default_annotation_gene_num_perchr = 1000;
const long default_annotation_chr_num = 50;

vector<group> global_sample_groups;

int total_sample_num_initial = 0;

/************************************************************************/
/* INPUT DATA                                                           */
/************************************************************************/


void input_ASMs(string resultPath, string outputPath)
{
	string filename, cur_chrNM, cur_ASMid;
	int groupLoopCnt, sampleLoopCnt, pathLoopCnt, pathCnt, badEstimationCnt, sampleIndexCnt;
	bool flag_filter_ASMexpr, flag_filter_pathexpr;
	double groupSum_expr, groupSum_prop, grandSum_expr, grandSum_prop, minPathExpression, MSE_estimation, minASMExpression, tmp_entry, grand_SS, max_path_var;
	long ASMrangeLow, ASMrangeHigh;
	ASM *newASM;
	string info;

	ifstream testASM_file;
	filename = resultPath + "/stat/asm.txt";
	testASM_file.open(filename.c_str());

	ofstream filteredASM_file;
	filename = outputPath + "/stat/filteredASM.txt";
	filteredASM_file.open(filename.c_str());

	ofstream tmp_ASMid_file;
	filename = outputPath + "/stat/ASM_id.txt";
	tmp_ASMid_file.open(filename.c_str());

	while (testASM_file >> cur_ASMid)
	{		
		testASM_file >> cur_chrNM;
		testASM_file >> ASMrangeLow;
		testASM_file >> ASMrangeHigh;
		testASM_file >> pathCnt;

		newASM = new ASM(pathCnt);
		flag_filter_ASMexpr = false;
		minASMExpression = MAX_NUMBER;

		newASM->ASM_id = cur_ASMid;
		newASM->testID = sortList_ASM.size();
		newASM->chrNM = cur_chrNM;
		newASM->rangeLow = ASMrangeLow;
		newASM->rangeHigh = ASMrangeHigh;

		testASM_file >> newASM->ASMcategory;

		tmp_ASMid_file << cur_ASMid << "\t" << pathCnt << "\t" << newASM->ASMcategory << endl;

		//input ASM expression
		for (sampleLoopCnt = 1; sampleLoopCnt <= global_num_sample_mask*2; ++sampleLoopCnt)
			testASM_file >> tmp_entry; //omit first global_num_sample_mask samples
		sampleIndexCnt = 1;
		grandSum_expr = 0.0;
		for (groupLoopCnt = 1; groupLoopCnt <= global_total_group_num; ++groupLoopCnt)
		{
			groupSum_expr = 0.0;
			badEstimationCnt = 0;
			for (sampleLoopCnt = 1; sampleLoopCnt <= global_group_size[groupLoopCnt]; ++sampleLoopCnt, ++sampleIndexCnt)
			{
				testASM_file >> newASM->sample_expr[sampleIndexCnt];
				testASM_file >> MSE_estimation;
				newASM->error_ratio[sampleIndexCnt] = MSE_estimation;

				groupSum_expr += newASM->sample_expr[sampleIndexCnt];
				if (MSE_estimation > THRESH_BAD_ESTIMATION)
					++badEstimationCnt;
			}
			grandSum_expr += groupSum_expr;
			newASM->group_mean_expr[groupLoopCnt] = groupSum_expr / global_group_size[groupLoopCnt];
			if (newASM->group_mean_expr[groupLoopCnt] < minASMExpression)
				minASMExpression = newASM->group_mean_expr[groupLoopCnt];

			if (newASM->group_mean_expr[groupLoopCnt] < THRESHOLD_MIN_ASM_COVERAGE)// || badEstimationCnt > SAMPLE_CNT_PER_GROUP/2)
			{
				flag_filter_ASMexpr = true;
			}
		}	
		newASM->min_group_mean_expr = minASMExpression;
		newASM->grand_mean_expr = grandSum_expr / global_total_sample_num;
		getline(testASM_file, info);

		if (COUNT_INTRON_RETENTION == false && newASM->ASMcategory == 3)
		{
			//throw intron retention
			flag_filter_ASMexpr = true;
		}

		//filtering by mean expression
		if (flag_filter_ASMexpr == true)
		{
			filteredASM_file << newASM->chrNM << "\t" << newASM->rangeLow << "\t" << newASM->rangeHigh << "\t" << newASM->pathNum << "\t" << newASM->ASMcategory;
			for (groupLoopCnt = 1; groupLoopCnt <= global_total_group_num; ++groupLoopCnt)
				filteredASM_file << "\t" << newASM->group_mean_expr[groupLoopCnt];
			filteredASM_file << endl;
			for (pathLoopCnt = 1; pathLoopCnt <= pathCnt; ++pathLoopCnt)
			{
				getline(testASM_file, info);
				filteredASM_file << info << endl;
			}
			delete newASM;
		}
		else
		{
			//input ASM paths
			flag_filter_pathexpr = false;
			minPathExpression = MAX_NUMBER;
			max_path_var = 0;

			for (pathLoopCnt = 1; pathLoopCnt <= pathCnt; ++pathLoopCnt)
			{
				string cur_path_id;
				testASM_file >> cur_path_id;
				newASM->path_id.push_back(cur_path_id);

				testASM_file >> ASMrangeLow;
				testASM_file >> ASMrangeHigh;

				for (sampleLoopCnt = 1; sampleLoopCnt <= global_num_sample_mask*2; ++sampleLoopCnt)
					testASM_file >> tmp_entry; //omit first global_num_sample_mask samples

				sampleIndexCnt = 1;
				grandSum_expr = 0.0;
				grandSum_prop = 0.0;
				grand_SS = 0.0;
				for (groupLoopCnt = 1; groupLoopCnt <= global_total_group_num; ++groupLoopCnt)
				{
					groupSum_expr = 0.0;
					groupSum_prop = 0.0;
					for (sampleLoopCnt = 1; sampleLoopCnt <= global_group_size[groupLoopCnt]; ++sampleLoopCnt, ++sampleIndexCnt)
					{
						testASM_file >> newASM->sample_path_expr[sampleIndexCnt][pathLoopCnt];
						testASM_file >> newASM->sample_path_prop[sampleIndexCnt][pathLoopCnt];
						groupSum_expr += newASM->sample_path_expr[sampleIndexCnt][pathLoopCnt];
						groupSum_prop += newASM->sample_path_prop[sampleIndexCnt][pathLoopCnt];
						grand_SS += pow(newASM->sample_path_expr[sampleIndexCnt][pathLoopCnt], 2);
					}

					grandSum_expr += groupSum_expr;
					grandSum_prop += groupSum_prop;
					groupSum_expr /= global_group_size[groupLoopCnt]; //calculate the group mean expression of this path
					groupSum_prop /= global_group_size[groupLoopCnt]; //calculate the group mean proportion of this path
					if (groupSum_expr < minPathExpression)
						minPathExpression = groupSum_expr;

					newASM->group_mean_path_expr[groupLoopCnt][pathLoopCnt] = groupSum_expr;
					newASM->group_mean_path_prop[groupLoopCnt][pathLoopCnt] = groupSum_prop;

					if (groupSum_expr < THRESHOLD_MIN_ASMPATH_COVERAGE)
					{
						flag_filter_pathexpr = true;
					}
				}
				newASM->grand_mean_path_expr[pathLoopCnt] = grandSum_expr / global_total_sample_num;
				newASM->grand_mean_path_prop[pathLoopCnt] = grandSum_prop / global_total_sample_num;

				newASM->grand_var_path_expr[pathLoopCnt] = grand_SS / global_total_sample_num - pow(newASM->grand_mean_path_expr[pathLoopCnt], 2);
				if (newASM->grand_var_path_expr[pathLoopCnt] > max_path_var){
					max_path_var = newASM->grand_var_path_expr[pathLoopCnt];
					newASM->rep_path = pathLoopCnt;
				}

				getline(testASM_file, info);
			}
			newASM->min_group_mean_path_expr = minPathExpression;

			if (flag_filter_pathexpr == true)
			{
				filteredASM_file << newASM->chrNM << "\t" << newASM->rangeLow << "\t" << newASM->rangeHigh << "\t" << newASM->pathNum << "\t" << newASM->ASMcategory;
				for (groupLoopCnt = 1; groupLoopCnt <= global_total_group_num; ++groupLoopCnt)
					filteredASM_file << "\t" << newASM->group_mean_expr[groupLoopCnt] << endl;
				for (pathLoopCnt = 1; pathLoopCnt <= pathCnt; ++pathLoopCnt)
				{
					filteredASM_file << ASMrangeLow << "\t" << ASMrangeHigh;
					for (groupLoopCnt = 1; groupLoopCnt <= global_total_group_num; ++groupLoopCnt)
						filteredASM_file << "\t" << newASM->group_mean_path_expr[groupLoopCnt][pathLoopCnt] << "\t" << newASM->group_mean_path_prop[groupLoopCnt][pathLoopCnt];
					filteredASM_file << endl;
				}
				delete newASM;
			}
			else
			{
				if (sortList_ASM.size() >= sortList_ASM.capacity())
					sortList_ASM.reserve(sortList_ASM.size() + DEFAULT_TESTASM_NUM);
				sortList_ASM.push_back(newASM);
			}
		}
	}

	testASM_file.close();
	filteredASM_file.close();
	tmp_ASMid_file.close();

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
			if (parameter.compare("global_total_sample_num") == 0)
				config_file >> global_total_sample_num;
			else if (parameter.compare("global_total_group_num") == 0)
				config_file >> global_total_group_num;
			else if (parameter.compare("COUNT_INTRON_RETENTION") == 0)
			{
				config_file >> info;
				if (info.compare("true") == 0)
					COUNT_INTRON_RETENTION = true;
				else
					COUNT_INTRON_RETENTION = false;
			}
			else if (parameter.compare("num_samples_mask") == 0)
				config_file >> global_num_sample_mask;
			else if (parameter.compare("global_group_size") == 0)
			{
				global_group_size = new int [1+global_total_group_num];
				global_group_size[0] = 0;
				for (int tmpCnt = 1; tmpCnt <= global_total_group_num; ++tmpCnt)
					config_file >> global_group_size[tmpCnt];
			}
			else if (parameter.compare("global_group_base_index") == 0)
			{
				global_group_base_index = new int [1+global_total_group_num];
				global_group_base_index[0] = 0;
				for (int tmpCnt = 1; tmpCnt <= global_total_group_num; ++tmpCnt)
					config_file >> global_group_base_index[tmpCnt];
			}
			else if (parameter.compare("group_name") == 0)
			{
				global_group_name = new string [1+global_total_group_num];
				global_group_name[0] = "";
				for (int tmpCnt = 1; tmpCnt <= global_total_group_num; ++tmpCnt)
					config_file >> global_group_name[tmpCnt];
			}
			else if (parameter.compare("false_discovery_rate") == 0)
				config_file >> false_discovery_rate;
			else if (parameter.compare("thresh_JSD") == 0)
				config_file >> thresh_JSD;
			else if (parameter.compare("thresh_stat_d") == 0)
				config_file >> thresh_stat_d;
			else if (parameter.compare("THRESHOLD_MIN_ASM_COVERAGE") == 0)
				config_file >> THRESHOLD_MIN_ASM_COVERAGE;

			getline(config_file, info);
		}
	} 
	else
	{
		cout << "Error: fail to open the config file for differential transcription test" << endl;
		exit(1);
	}

	config_file.close();

	if (global_total_sample_num <= 0 || global_total_group_num <= 0)
	{
		cout << "Error: incomplete config file for differential transcription test" << endl;
		exit(1);
	}

	return;
}







void input_annotation_gaf(string filename)
{
	string cur_exonboundary, tmpStr;
	ifstream genelistfile;
	long tmpboundary; 
	annotation_gene *new_annogene;
	annotation_chr_genelist *new_anno_chr;
	unsigned long cur_chr_index = 0;

	genelistfile.open(filename.c_str());

	if (!genelistfile.is_open())
	{
		cout << "warning: fail to open gene annotation file" << endl;
		return;
	}

	getline(genelistfile, tmpStr);
	while (getline(genelistfile, tmpStr, '\t'))
	{
		new_annogene = new annotation_gene;

		getline(genelistfile, new_annogene->featureID, '\t');
		for (int tmpLoop = 1; tmpLoop <= 12; ++tmpLoop)
			getline(genelistfile, tmpStr, '\t');

		// 		getline(genelistfile, tmpStr, ':');
		// 		getline(genelistfile, cur_exonboundary, ':');
		// 		getline(genelistfile, tmpStr, '\t');
		// 		getline(genelistfile, new_annogene->geneNm, '\t');
		// 		getline(genelistfile, new_annogene->chr, ':');
		// 		getline(genelistfile, tmpStr, ':');
		// 		genelistfile >> new_annogene->strand;
		// 		getline(genelistfile, tmpStr);

		getline(genelistfile, new_annogene->chr, ':');
		getline(genelistfile, cur_exonboundary, ':');
		getline(genelistfile, new_annogene->strand, '\t');
		getline(genelistfile, new_annogene->geneNm, '\t');
		getline(genelistfile, tmpStr);

		while (cur_exonboundary.empty() == false)
		{
			size_t found = cur_exonboundary.find_first_of(" -,\n");
			if (found != string::npos)
			{
				tmpStr = cur_exonboundary.substr(0, found);
				cur_exonboundary.erase(0, found+1);
			}
			else
			{
				tmpStr = cur_exonboundary.substr(0);
				cur_exonboundary.clear();
			}

			tmpboundary = atol(tmpStr.c_str());
			tmpboundary = abs(tmpboundary);

// 			if (new_annogene->exonBoundary.size() >= new_annogene->exonBoundary.capacity())
// 				new_annogene->exonBoundary.reserve(new_annogene->exonBoundary.capacity() + 10);
// 			new_annogene->exonBoundary.push_back(tmpboundary);

			if (tmpboundary < new_annogene->range_start)
				new_annogene->range_start = tmpboundary;
			if (tmpboundary > new_annogene->range_end)
				new_annogene->range_end = tmpboundary;
		}

		for (cur_chr_index = 0; cur_chr_index < list_anno_chr_genelist.size(); ++cur_chr_index)
			if (new_annogene->chr.compare(list_anno_chr_genelist[cur_chr_index]->chrname) == 0)
				break;

		if (cur_chr_index >= list_anno_chr_genelist.size())
		{
			new_anno_chr = new annotation_chr_genelist;
			new_anno_chr->chrname = new_annogene->chr;

			if (list_anno_chr_genelist.size() >= list_anno_chr_genelist.capacity())
				list_anno_chr_genelist.reserve(list_anno_chr_genelist.capacity() + default_annotation_chr_num);
			list_anno_chr_genelist.push_back(new_anno_chr);

			if (new_anno_chr->genelist.size() >= new_anno_chr->genelist.capacity())
				new_anno_chr->genelist.reserve(new_anno_chr->genelist.capacity() + default_annotation_gene_num_perchr);
			new_anno_chr->genelist.push_back(new_annogene);
		} 
		else
		{
			if (list_anno_chr_genelist[cur_chr_index]->genelist.size() >= list_anno_chr_genelist[cur_chr_index]->genelist.capacity())
				list_anno_chr_genelist[cur_chr_index]->genelist.reserve(list_anno_chr_genelist[cur_chr_index]->genelist.capacity() + default_annotation_gene_num_perchr);
			list_anno_chr_genelist[cur_chr_index]->genelist.push_back(new_annogene);
		}
	}


	//Nov 12, 2012, sort gene list based on their start position

	for (cur_chr_index = 0; cur_chr_index < list_anno_chr_genelist.size(); ++cur_chr_index)
	{
		void **sortlist = new void * [list_anno_chr_genelist[cur_chr_index]->genelist.size() + 1];
		double *sortkey = new double [list_anno_chr_genelist[cur_chr_index]->genelist.size() + 1];

		for (unsigned long iLoop = 0; iLoop < list_anno_chr_genelist[cur_chr_index]->genelist.size(); ++iLoop)
		{
			sortkey[iLoop + 1] = list_anno_chr_genelist[cur_chr_index]->genelist[iLoop]->range_start;
			sortlist[iLoop + 1] = (void*) list_anno_chr_genelist[cur_chr_index]->genelist[iLoop];
		}
		mergeSort(sortlist, sortkey, list_anno_chr_genelist[cur_chr_index]->genelist.size());

		for (unsigned long iLoop = 0; iLoop < list_anno_chr_genelist[cur_chr_index]->genelist.size(); ++iLoop)
			list_anno_chr_genelist[cur_chr_index]->genelist[iLoop] = (annotation_gene*) sortlist[iLoop + 1];

		delete [] sortlist;
		delete [] sortkey;
	}

	return;
}

bool override_sample_groups(string filename){
	ifstream infile;
	infile.open(filename.c_str());
	if (!infile.is_open()){
		cout << "Error: a new sample group specification has been provided to override the existing groups, but the file " << filename << "cannot be opened." << endl;
		return false;
	}

	string group_name, sample_name, info;
	long sample_idx, sample_total_cnt = 0;

	while (infile >> group_name){
		infile >> sample_name;
		infile >> sample_idx;
		getline(infile, info);
		++sample_total_cnt;

		if (global_sample_groups.empty() || group_name.compare(global_sample_groups.back().group_name) != 0){
			global_sample_groups.push_back(group (group_name));
		}

		global_sample_groups.back().sample_name.push_back(sample_name);
		global_sample_groups.back().sample_idx.push_back(sample_idx);
	}

	// now update the configurations
	total_sample_num_initial = global_total_sample_num;
	global_total_sample_num = sample_total_cnt;
	global_total_group_num = global_sample_groups.size();
	global_num_sample_mask = 0;
	
	delete [] global_group_size;
	delete [] global_group_base_index;
	delete [] global_group_name;
	global_group_size = new int [1+global_total_group_num];
	global_group_size[0] = 0;
	global_group_base_index = new int [1+global_total_group_num];
	global_group_base_index[0] = 0; global_group_base_index[1] = 0;
	global_group_name = new string [1+global_total_group_num];
	global_group_name[0] = "";
	for (int tmpCnt = 0; tmpCnt < global_total_group_num; ++tmpCnt){
		global_group_size[tmpCnt+1] = global_sample_groups[tmpCnt].sample_idx.size();		
		global_group_name[tmpCnt+1] = global_sample_groups[tmpCnt].group_name;
	}
	for (int tmpCnt = 1; tmpCnt < global_total_group_num; ++tmpCnt){
		global_group_base_index[tmpCnt+1] = global_group_base_index[tmpCnt] + global_group_size[tmpCnt];
	}

	// output overridden configuration
	cout << "Sample groups have been overridden: " << global_total_sample_num << " samples with " << global_total_group_num << "groups; groups are name(size, base_index): ";
	for (int tmpCnt = 0; tmpCnt < global_total_group_num; ++tmpCnt){
		cout << global_group_name[tmpCnt+1] << "(" << global_group_size[tmpCnt+1] << "," << global_group_base_index[tmpCnt+1] << ") ";		
	}
	cout << endl;
		
	return true;
}

// input ASM for regrouped samples (sample groups are overridden)
void input_ASMs_regroup(string resultPath, string outputPath)
{
	string filename, cur_chrNM, cur_ASMid;
	int groupLoopCnt, sampleLoopCnt, pathLoopCnt, pathCnt, badEstimationCnt, sampleIndexCnt;
	bool flag_filter_ASMexpr, flag_filter_pathexpr;
	double groupSum_expr, groupSum_prop, grandSum_expr, grandSum_prop, minPathExpression, MSE_estimation, minASMExpression, tmp_entry, grand_SS, max_path_var;
	long ASMrangeLow, ASMrangeHigh;
	ASM *newASM;
	string info;

	double *tmp_store_1 = new double [total_sample_num_initial+1]; // use the initial count !!
	double *tmp_store_2 = new double [total_sample_num_initial+1];
	tmp_store_1[0] = -1;
	tmp_store_2[0] = -1;

	ifstream testASM_file;
	filename = resultPath + "/stat/asm.txt";
	testASM_file.open(filename.c_str());

	ofstream filteredASM_file;
	filename = outputPath + "/stat/filteredASM.txt";
	filteredASM_file.open(filename.c_str());

	ofstream tmp_ASMid_file;
	filename = outputPath + "/stat/ASM_id.txt";
	tmp_ASMid_file.open(filename.c_str());

	while (testASM_file >> cur_ASMid)
	{		
		testASM_file >> cur_chrNM;
		testASM_file >> ASMrangeLow;
		testASM_file >> ASMrangeHigh;
		testASM_file >> pathCnt;

		newASM = new ASM(pathCnt);
		flag_filter_ASMexpr = false;
		minASMExpression = MAX_NUMBER;

		newASM->ASM_id = cur_ASMid;
		newASM->testID = sortList_ASM.size();
		newASM->chrNM = cur_chrNM;
		newASM->rangeLow = ASMrangeLow;
		newASM->rangeHigh = ASMrangeHigh;

		testASM_file >> newASM->ASMcategory;

		tmp_ASMid_file << cur_ASMid << "\t" << pathCnt << "\t" << newASM->ASMcategory << endl;

		//input ASM expression
		for (sampleLoopCnt = 1; sampleLoopCnt <= total_sample_num_initial; ++sampleLoopCnt) {
			testASM_file >> tmp_store_1[sampleLoopCnt]; // first read the values to the temporary array
			testASM_file >> tmp_store_2[sampleLoopCnt]; // depreciated
		}

		sampleIndexCnt = 1;
		grandSum_expr = 0.0;
		for (groupLoopCnt = 1; groupLoopCnt <= global_total_group_num; ++groupLoopCnt)
		{
			groupSum_expr = 0.0;
			badEstimationCnt = 0;
			for (sampleLoopCnt = 1; sampleLoopCnt <= global_group_size[groupLoopCnt]; ++sampleLoopCnt, ++sampleIndexCnt)
			{
				// testASM_file >> newASM->sample_expr[sampleIndexCnt];
				// testASM_file >> MSE_estimation;
				newASM->sample_expr[sampleIndexCnt] = tmp_store_1[global_sample_groups[groupLoopCnt-1].sample_idx[sampleLoopCnt-1]];
				newASM->error_ratio[sampleIndexCnt] = tmp_store_2[global_sample_groups[groupLoopCnt-1].sample_idx[sampleLoopCnt-1]];

				groupSum_expr += newASM->sample_expr[sampleIndexCnt];
				if (MSE_estimation > THRESH_BAD_ESTIMATION)
					++badEstimationCnt;
			}
			grandSum_expr += groupSum_expr;
			newASM->group_mean_expr[groupLoopCnt] = groupSum_expr / global_group_size[groupLoopCnt];
			if (newASM->group_mean_expr[groupLoopCnt] < minASMExpression)
				minASMExpression = newASM->group_mean_expr[groupLoopCnt];

			if (newASM->group_mean_expr[groupLoopCnt] < THRESHOLD_MIN_ASM_COVERAGE)// || badEstimationCnt > SAMPLE_CNT_PER_GROUP/2)
			{
				flag_filter_ASMexpr = true;
			}
		}	
		newASM->min_group_mean_expr = minASMExpression;
		newASM->grand_mean_expr = grandSum_expr / global_total_sample_num;
		getline(testASM_file, info);

		if (COUNT_INTRON_RETENTION == false && newASM->ASMcategory == 3)
		{
			//throw intron retention
			flag_filter_ASMexpr = true;
		}

		//filtering by mean expression
		if (flag_filter_ASMexpr == true)
		{
			filteredASM_file << newASM->chrNM << "\t" << newASM->rangeLow << "\t" << newASM->rangeHigh << "\t" << newASM->pathNum << "\t" << newASM->ASMcategory;
			for (groupLoopCnt = 1; groupLoopCnt <= global_total_group_num; ++groupLoopCnt)
				filteredASM_file << "\t" << newASM->group_mean_expr[groupLoopCnt];
			filteredASM_file << endl;
			for (pathLoopCnt = 1; pathLoopCnt <= pathCnt; ++pathLoopCnt)
			{
				getline(testASM_file, info);
				filteredASM_file << info << endl;
			}
			delete newASM;
		}
		else
		{
			//input ASM paths
			flag_filter_pathexpr = false;
			minPathExpression = MAX_NUMBER;
			max_path_var = 0;

			for (pathLoopCnt = 1; pathLoopCnt <= pathCnt; ++pathLoopCnt)
			{
				string cur_path_id;
				testASM_file >> cur_path_id;
				newASM->path_id.push_back(cur_path_id);

				testASM_file >> ASMrangeLow;
				testASM_file >> ASMrangeHigh;

				for (sampleLoopCnt = 1; sampleLoopCnt <= total_sample_num_initial; ++sampleLoopCnt) {
					testASM_file >> tmp_store_1[sampleLoopCnt]; // first read the values to the temporary array
					testASM_file >> tmp_store_2[sampleLoopCnt]; // depreciated
				}

				sampleIndexCnt = 1;
				grandSum_expr = 0.0;
				grandSum_prop = 0.0;
				grand_SS = 0.0;
				for (groupLoopCnt = 1; groupLoopCnt <= global_total_group_num; ++groupLoopCnt)
				{
					groupSum_expr = 0.0;
					groupSum_prop = 0.0;
					for (sampleLoopCnt = 1; sampleLoopCnt <= global_group_size[groupLoopCnt]; ++sampleLoopCnt, ++sampleIndexCnt)
					{
						//testASM_file >> newASM->sample_path_expr[sampleIndexCnt][pathLoopCnt];
						//testASM_file >> newASM->sample_path_prop[sampleIndexCnt][pathLoopCnt];
						newASM->sample_path_expr[sampleIndexCnt][pathLoopCnt] = tmp_store_1[global_sample_groups[groupLoopCnt-1].sample_idx[sampleLoopCnt-1]];
						newASM->sample_path_prop[sampleIndexCnt][pathLoopCnt] = tmp_store_2[global_sample_groups[groupLoopCnt-1].sample_idx[sampleLoopCnt-1]];

						groupSum_expr += newASM->sample_path_expr[sampleIndexCnt][pathLoopCnt];
						groupSum_prop += newASM->sample_path_prop[sampleIndexCnt][pathLoopCnt];
						grand_SS += pow(newASM->sample_path_expr[sampleIndexCnt][pathLoopCnt], 2);
					}

					grandSum_expr += groupSum_expr;
					grandSum_prop += groupSum_prop;
					groupSum_expr /= global_group_size[groupLoopCnt]; //calculate the group mean expression of this path
					groupSum_prop /= global_group_size[groupLoopCnt]; //calculate the group mean proportion of this path
					if (groupSum_expr < minPathExpression)
						minPathExpression = groupSum_expr;

					newASM->group_mean_path_expr[groupLoopCnt][pathLoopCnt] = groupSum_expr;
					newASM->group_mean_path_prop[groupLoopCnt][pathLoopCnt] = groupSum_prop;

					if (groupSum_expr < THRESHOLD_MIN_ASMPATH_COVERAGE)
					{
						flag_filter_pathexpr = true;
					}
				}
				newASM->grand_mean_path_expr[pathLoopCnt] = grandSum_expr / global_total_sample_num;
				newASM->grand_mean_path_prop[pathLoopCnt] = grandSum_prop / global_total_sample_num;

				newASM->grand_var_path_expr[pathLoopCnt] = grand_SS / global_total_sample_num - pow(newASM->grand_mean_path_expr[pathLoopCnt], 2);
				if (newASM->grand_var_path_expr[pathLoopCnt] > max_path_var){
					max_path_var = newASM->grand_var_path_expr[pathLoopCnt];
					newASM->rep_path = pathLoopCnt;
				}

				getline(testASM_file, info);
			}
			newASM->min_group_mean_path_expr = minPathExpression;

			if (flag_filter_pathexpr == true)
			{
				filteredASM_file << newASM->chrNM << "\t" << newASM->rangeLow << "\t" << newASM->rangeHigh << "\t" << newASM->pathNum << "\t" << newASM->ASMcategory;
				for (groupLoopCnt = 1; groupLoopCnt <= global_total_group_num; ++groupLoopCnt)
					filteredASM_file << "\t" << newASM->group_mean_expr[groupLoopCnt] << endl;
				for (pathLoopCnt = 1; pathLoopCnt <= pathCnt; ++pathLoopCnt)
				{
					filteredASM_file << ASMrangeLow << "\t" << ASMrangeHigh;
					for (groupLoopCnt = 1; groupLoopCnt <= global_total_group_num; ++groupLoopCnt)
						filteredASM_file << "\t" << newASM->group_mean_path_expr[groupLoopCnt][pathLoopCnt] << "\t" << newASM->group_mean_path_prop[groupLoopCnt][pathLoopCnt];
					filteredASM_file << endl;
				}
				delete newASM;
			}
			else
			{
				if (sortList_ASM.size() >= sortList_ASM.capacity())
					sortList_ASM.reserve(sortList_ASM.size() + DEFAULT_TESTASM_NUM);
				sortList_ASM.push_back(newASM);
			}
		}
	}

	testASM_file.close();
	filteredASM_file.close();
	tmp_ASMid_file.close();

	delete [] tmp_store_1;
	delete [] tmp_store_2;

	return;
}