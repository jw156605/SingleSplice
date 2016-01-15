#include "config.h"

#ifndef ASM_DEF
#define ASM_DEF


const char alterSpliceCategory[][30] = {"combination", "exon_skipping", "mutual_exclusive", "intron_retention", "alter_splice_site", "alter_trans_start/end", "alter_trans_start/end", "diff_expression", "nested"};

//difference between a pair of groups
class pairdiff
{
public:
	double stat_d;
	double stat_d_expected;
	double stat_difftrans;
	double jsd;
	int groupindex_1;
	int groupindex_2;
	bool significant_stat; //statistically significant
	bool significant_both; //both statistically and practically significant

	pairdiff();
};

class ASM
{
public:
	//basic information
	string ASM_id;
	long testID;
	string chrNM;
	long rangeLow;
	long rangeHigh;
	int pathNum;
	int ASMcategory;
	string in_gene;

	vector <string> path_id;

	//expression and proportion information
	double* sample_expr; //previously individualExpression, absolute expression of the entire ASM for each sample 
	double* error_ratio; //relative ratio of the estimation error against the expression
	double* group_mean_expr; //previously meanExpression, mean expression of each group
	double grand_mean_expr; //mean expression of all samples (of all groups)
	double min_group_mean_expr; //previously minExpression, minimum of all group mean expression
	double min_group_mean_path_expr; //previously minPathExpression

	double** sample_path_expr; //previously individualPathExpression
	double** sample_path_prop; //previously individualPathProportion
	double** group_mean_path_expr; //previously meanPathExpression
	double** group_mean_path_prop; //previously meanPathProportion
	double* grand_mean_path_expr; //mean path expression of all groups
	double* grand_mean_path_prop; //mean path proportion of all groups

	// for select a representative path 10/7/2014
	double* grand_var_path_expr; //variance of path expression of all groups
	int rep_path;

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

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//practical difference, measured by non-probabilistic (non-distributional) statistics
	//1. group mean JSD, default
	double max_pairwise_JSD_betwgroup;
	string max_pairwise_JSD_betwgroup_groupname;
	string info_large_pairwise_JSD_betwgroup;
	//2. min distance from group members to other group's mean, used for group separation, difference between the two whole groups
	double max_grouppair_min_dist_ratio; //for a sample in group 1, calculate the ratio of the distance from this point to group mean of group 2 over the distance from this point to group mean of group 1, take the min of all point in a group for the group, as the gap; for two groups, calculate both min ratio, then take the min
	string max_grouppair_min_dist_ratio_groupname;
	string info_large_grouppair_min_dist_ratio;
	//3. percentage of group members with large distance to other group's members, i.e., the percentage of members in one group that have distances to all members of another group larger than a threshold, used when there may be multiple modes in a group such as in the TCGA analysis
	double max_percent_large_dist;
	string max_percent_large_dist_groupname;
	string info_large_percent_large_dist;


	//sample-sample clustering
//	double** sample_dist_mat_trans; //sample-sample differential transcription distance matrix

	ASM(int pathnum);
	~ASM();
}; 

extern vector <ASM*> sortList_ASM; //ASM list to be sorted

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

