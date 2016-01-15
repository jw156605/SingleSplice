#include "gene_def.h"

pairdiff::pairdiff()
{
	stat_d = 0;
	stat_d_expected = 0;
	stat_diffexpr = 0;
	foldchange = 0;
	groupindex_1 = 0;
	groupindex_2 = 2;
	significant_stat = false;
	significant_both = false;
}

gene::gene()
{
	//basic information
	gene_id = "gene0";
	testID = 0;
	rangeLow = 0;
	rangeHigh = 0;
	in_gene = "na";
	
	//expression and proportion information
	sample_expr = new double [global_total_sample_num+1];
	group_mean_expr = new double [global_total_group_num+1];
	grand_mean_expr = 0;

	statistic_d = 0;
	statistic_s = 0;
	statistic_s_expr = 0;
	statistic_d_expected = 0;

	//whether significantly different or not
	significant_stat = false;
	significant_both = false;

	max_pairwise_foldchange = 0;
}


gene::~gene()
{
	delete [] sample_expr;
	delete [] group_mean_expr;

	for (vector<pairdiff*>::iterator it_tmp = list_pairdiff.begin(); it_tmp != list_pairdiff.end(); ++it_tmp)
	{
		delete *it_tmp;
	}
	list_pairdiff.clear();
}












