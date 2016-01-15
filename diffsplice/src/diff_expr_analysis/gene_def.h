#include "config.h"

#ifndef ASM_DEF
#define ASM_DEF

//difference between a pair of groups
class pairdiff
{
public:
	double stat_d;
	double stat_d_expected;
	double stat_diffexpr;
	double foldchange;
	int groupindex_1;
	int groupindex_2;
	bool significant_stat; //statistically significant
	bool significant_both; //both statistically and practically significant

	pairdiff();
};

class gene
{
public:
	//basic information
	string gene_id;
	long testID;
	string chrNM;
	long rangeLow;
	long rangeHigh;
	string in_gene;

	vector <string> path_id;

	//expression information
	double* sample_expr; //previously individualExpression, absolute expression of the entire gene for each sample 
	double* group_mean_expr; //previously meanExpression, mean expression of each group
	double grand_mean_expr; //mean expression of all samples (of all groups)
		
	//statistics of all-together
	double statistic_d;
	double statistic_s;
	double statistic_s_expr;
	double statistic_d_expected;

	//whether significantly different or not
	bool significant_stat;
	bool significant_both;

	//statistics of group pair comparison
	vector <pairdiff*> list_pairdiff;	

	//maximum pairwise difference
	double max_pairwise_foldchange;
	string max_pairwise_foldchange_groupname;
	string info_large_pairwise_foldchange;
		
	gene();
	~gene();
}; 

extern vector <gene*> sortList_ASM; //ASM list to be sorted

class annotation_gene
{
public:
	string featureID;
	string geneNm;
	string chr;
	string strand;
	long range_start;
	long range_end;
	//vector<long> exonBoundary;

	annotation_gene() {range_start = long(MAX_NUMBER); range_end = 0;}
};

class annotation_chr_genelist
{
public:
	string chrname;
	vector<annotation_gene*> genelist;
};

extern vector <annotation_chr_genelist*> list_anno_chr_genelist;

#endif

