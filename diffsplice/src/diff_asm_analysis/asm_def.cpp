#include "asm_def.h"

pairdiff::pairdiff()
{
	stat_d = 0;
	stat_d_expected = 0;
	stat_difftrans = 0;
	jsd = 0;
	groupindex_1 = 0;
	groupindex_2 = 2;
	significant_stat = false;
	significant_both = false;
}

ASM::ASM(int pathnum)
{
	//basic information
	ASM_id = "ASM0";
	testID = 0;
	rangeLow = 0;
	rangeHigh = 0;
	pathNum = pathnum;
	ASMcategory = 0;
	in_gene = "na";
	
	//expression and proportion information
	long groupLoop, sampleLoop;

	sample_expr = new double [global_total_sample_num+1];
	error_ratio = new double [global_total_sample_num+1];
	group_mean_expr = new double [global_total_group_num+1];
	
	grand_mean_expr = 0;
	min_group_mean_expr = 0;
	min_group_mean_path_expr = 0;

	sample_path_expr = new double* [global_total_sample_num+1];
	sample_path_prop = new double* [global_total_sample_num+1];
	group_mean_path_expr = new double* [global_total_group_num+1];
	group_mean_path_prop = new double* [global_total_group_num+1];
	for (sampleLoop = 1; sampleLoop <= global_total_sample_num; ++sampleLoop)
	{
		sample_path_expr[sampleLoop] = new double [pathnum+1];
		sample_path_prop[sampleLoop] = new double [pathnum+1];
	}
	for (groupLoop = 1; groupLoop <= global_total_group_num; ++groupLoop)
	{
		group_mean_path_expr[groupLoop] = new double [pathnum+1];
		group_mean_path_prop[groupLoop] = new double [pathnum+1];
	}
	grand_mean_path_expr = new double [pathnum+1];
	grand_mean_path_prop = new double [pathnum+1];
	grand_var_path_expr = new double [pathnum+1];
	rep_path = 0;
		
	statistic_d = 0;
	statistic_s = 0;
	statistic_s_expr = 0;
	statistic_d_expected = 0;

	//whether significantly different or not
	significant_stat = false;
	significant_both = false;

	//practical difference
	max_pairwise_JSD_betwgroup = 0;
	max_grouppair_min_dist_ratio = 0;
	max_percent_large_dist = 0;

	//sample-sample clustering
//	sample_dist_mat_trans = new double* [global_total_sample_num+1];
//	for (sampleLoop = 1; sampleLoop <= global_total_sample_num; ++sampleLoop)
//		sample_dist_mat_trans[sampleLoop] = new double [global_total_sample_num+1];
}


ASM::~ASM()
{
	delete [] sample_expr;
	delete [] error_ratio;
	delete [] group_mean_expr;

	long groupLoop, sampleLoop;

	for (sampleLoop = 1; sampleLoop <= global_total_sample_num; ++sampleLoop)
	{
		delete [] sample_path_expr[sampleLoop];
		delete [] sample_path_prop[sampleLoop];
	}
	for (groupLoop = 1; groupLoop <= global_total_group_num; ++groupLoop)
	{
		delete [] group_mean_path_expr[groupLoop];
		delete [] group_mean_path_prop[groupLoop];
	}
	delete [] sample_path_expr;
	delete [] sample_path_prop;
	delete [] group_mean_path_expr;
	delete [] group_mean_path_prop;

	for (vector<pairdiff*>::iterator it_tmp = list_pairdiff.begin(); it_tmp != list_pairdiff.end(); ++it_tmp)
	{
		delete *it_tmp;
	}
	list_pairdiff.clear();


//	for (sampleLoop = 1; sampleLoop <= global_total_sample_num; ++sampleLoop)
//	{
//		delete [] sample_dist_mat_trans[sampleLoop];
//	}
//	delete [] sample_dist_mat_trans;
}












