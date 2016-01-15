#include "input_data.h"

const long default_annotation_gene_num_perchr = 1000;
const long default_annotation_chr_num = 50;

/************************************************************************/
/* INPUT DATA                                                           */
/************************************************************************/


void input_genes(string resultPath)
{
	string filename, cur_chrNM;
	int groupLoopCnt, sampleLoopCnt, sampleIndexCnt;
	bool flag_filter;
	double groupSum_expr, grandSum_expr, tmp_entry;
	gene *newGene;
	string info;

	ifstream testGene_file;
	filename = resultPath + "stat/expression.txt";
	testGene_file.open(filename.c_str());

	ofstream filteredASM_file;
	filename = resultPath + "stat/filteredGene.txt";
	filteredASM_file.open(filename.c_str());


	while (testGene_file >> cur_chrNM)
	{		
		newGene = new gene;
		flag_filter = true;
		
		newGene->testID = sortList_ASM.size();
		newGene->chrNM = cur_chrNM;
		testGene_file >> newGene->rangeLow;
		testGene_file >> newGene->rangeHigh;

		//input gene expression
		for (sampleLoopCnt = 1; sampleLoopCnt <= global_num_sample_mask; ++sampleLoopCnt)
			testGene_file >> tmp_entry; //omit first global_num_sample_mask samples
		sampleIndexCnt = 1;
		grandSum_expr = 0.0;
		for (groupLoopCnt = 1; groupLoopCnt <= global_total_group_num; ++groupLoopCnt)
		{
			groupSum_expr = 0.0;
			for (sampleLoopCnt = 1; sampleLoopCnt <= global_group_size[groupLoopCnt]; ++sampleLoopCnt, ++sampleIndexCnt)
			{
				testGene_file >> newGene->sample_expr[sampleIndexCnt];
				
				groupSum_expr += newGene->sample_expr[sampleIndexCnt];
			}
			grandSum_expr += groupSum_expr;
			newGene->group_mean_expr[groupLoopCnt] = groupSum_expr / global_group_size[groupLoopCnt];
			
			if (newGene->group_mean_expr[groupLoopCnt] >= THRESHOLD_MIN_EXPR_COVERAGE)// || badEstimationCnt > SAMPLE_CNT_PER_GROUP/2)
			{
				flag_filter = false;
			}
		}	
		newGene->grand_mean_expr = grandSum_expr / global_total_sample_num;
		getline(testGene_file, info);

		
		//filtering by mean expression
		if (flag_filter == true)
		{
			filteredASM_file << newGene->chrNM << "\t" << newGene->rangeLow << "\t" << newGene->rangeHigh;
			for (groupLoopCnt = 1; groupLoopCnt <= global_total_group_num; ++groupLoopCnt)
				filteredASM_file << "\t" << newGene->group_mean_expr[groupLoopCnt];
			filteredASM_file << endl;
			delete newGene;
		}
		else
		{
			//input ASM paths
			if (sortList_ASM.size() >= sortList_ASM.capacity())
				sortList_ASM.reserve(sortList_ASM.size() + DEFAULT_TESTASM_NUM);
			sortList_ASM.push_back(newGene);
		}
	}

	testGene_file.close();
	filteredASM_file.close();

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
			else if (parameter.compare("thresh_foldchange_up") == 0)
				config_file >> thresh_foldchange_up;
			else if (parameter.compare("thresh_foldchange_down") == 0)
				config_file >> thresh_foldchange_down;
			else if (parameter.compare("thresh_stat_d") == 0)
				config_file >> thresh_stat_d;
			else if (parameter.compare("THRESHOLD_MIN_ASM_COVERAGE") == 0)
				config_file >> THRESHOLD_MIN_EXPR_COVERAGE;

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

