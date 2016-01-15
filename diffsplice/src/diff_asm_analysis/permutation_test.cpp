#include "difftrans_test.h"
#include "permutation_test.h"

/************************************************************************/
/* Statistic for a single gene                                          */
/************************************************************************/

// double calculate_single_gene_r_two_class(double *meanProp_group1, double *meanProp_group2, int dimension)
// {
// 	return sqrt(calculate_JSD(meanProp_group1, meanProp_group2, dimension));	
// }

//calculate the variance among groups
double calculate_single_gene_r_multi_class(double **group_mean_prop, double *grand_mean_prop, int dimension, int* group_size, int group_num)
{
	if (group_num <= 1)
	{
		cout << "warning: less than 2 groups were specified" << endl;
		return 0;
	}

	double stat_r = 0.0;
	for (int groupLoop = 1; groupLoop <= group_num; ++groupLoop)
	{
		stat_r += group_size[groupLoop] * calculate_JSD(group_mean_prop[groupLoop], grand_mean_prop, dimension);
	}	
	stat_r /= group_num - 1;

	return stat_r;
}

// double calculate_single_gene_s_two_class(double** individualProp_group1, double* meanProp_group1, int n1, double** individualProp_group2, double* meanProp_group2, int n2, int dimension)
// {
// 	if (n1 + n2 <= 2)
// 		return 0;
// 
// 	double SumJSD_group1 = 0.0, SumJSD_group2 = 0.0;
// 	int sampleLoopCnt;
// 
// 	for (sampleLoopCnt = 1; sampleLoopCnt <= n1; ++sampleLoopCnt)
// 	{
// 		SumJSD_group1 += calculate_JSD(individualProp_group1[sampleLoopCnt], meanProp_group1, dimension);
// 	}
// 	for (sampleLoopCnt = 1; sampleLoopCnt <= n2; ++sampleLoopCnt)
// 	{
// 		SumJSD_group2 += calculate_JSD(individualProp_group2[sampleLoopCnt], meanProp_group2, dimension);
// 	}
// 
// 	return sqrt((1.0/n1 + 1.0/n2) * (SumJSD_group1 + SumJSD_group2) / (n1+n2-2));
// }

double calculate_single_gene_s_multi_class(double **individual_prop, double **group_mean_prop, int dimension, int* group_size, int* group_assign, int group_num, int sample_num)
{
	if (group_num <= 1)
	{
		cout << "warning: less than 2 groups were specified" << endl;
		return 0;
	}

	double stat_s = 0.0;
	int individualIndex = 1;
	if (group_assign)
	{
		//with altered group assignment
		for (int groupLoop = 1; groupLoop <= group_num; ++groupLoop)
		{
			for (int individualLoop = 1; individualLoop <= group_size[groupLoop]; ++individualLoop, ++individualIndex)
			{
				stat_s += calculate_JSD(individual_prop[group_assign[individualIndex]], group_mean_prop[groupLoop], dimension);
			}
		}
	}
	else
	{
		//no group assignment provided, use the default
		for (int groupLoop = 1; groupLoop <= group_num; ++groupLoop)
		{
			for (int individualLoop = 1; individualLoop <= group_size[groupLoop]; ++individualLoop, ++individualIndex)
			{
				stat_s += calculate_JSD(individual_prop[individualIndex], group_mean_prop[groupLoop], dimension);
			}
		}
	}
		
	stat_s /= (sample_num - group_num);

	return stat_s;
}

// double calculate_single_gene_d(double** individualProp_group1, double* meanProp_group1, int n1, double** individualProp_group2, double* meanProp_group2, int n2, int dimension, double s0, double value_s)
// {
// 	//if (value_s < 0)
// 		//value_s has not been calculated, calculate now
// 		value_s = calculate_single_gene_s_two_class(individualProp_group1, meanProp_group1, n1, individualProp_group2, meanProp_group2, n2, dimension);
// 
// 	if (value_s + s0 < 1e-15)
// 	{
// 		//cout << "dividing a number ~= 0 in calculate_single_gene_d" << endl;
// 		return 1e-15;
// 	}
// 	return calculate_single_gene_r_two_class(meanProp_group1, meanProp_group2, dimension) / (value_s + s0);
// }
double calculate_single_gene_d(double** individual_prop, double** group_mean_prop, double* grand_mean_prop, int dimension, int* group_size, int* group_assign, int group_num, int sample_num, double s0, double value_s)
{
	if (value_s < 0)
		//value_s has not been calculated, calculate now
		value_s = calculate_single_gene_s_multi_class(individual_prop, group_mean_prop, dimension, group_size, group_assign, group_num, sample_num);

	if (value_s + s0 < 1e-15)
	{
		//variance is too small, do not count the difference
		//cout << "dividing a number ~= 0 in calculate_single_gene_d" << endl;
		return 1e-15;
	}
	return calculate_single_gene_r_multi_class(group_mean_prop, grand_mean_prop, dimension, group_size, group_num) / (value_s + s0);
}


/************************************************************************/
/* Statistic for all genes                                              */
/************************************************************************/

long permutationTest(double s0, string resultPath, double* permuted_d[], string resultFileSuffix)
{
	if (global_total_group_num < 2)
	{
		//one group test makes no sense
		return 0;
	}

//	ofstream debugfile_shuffle_index("debug_shuffled_index.txt");

	double **group_mean_prop_permuted, *permuted_value_array;
	long permutation_cnt, asm_loop, group_loop, sample_loop;
	int cur_dimension, *group_assign_index;

	group_mean_prop_permuted = new double* [global_total_group_num+1];
	group_assign_index = new int [global_total_sample_num+1];
	permuted_value_array = new double [sortList_ASM.size()];
	for (sample_loop = 1; sample_loop <= global_total_sample_num; ++sample_loop) {group_assign_index[sample_loop] = sample_loop;}

	for (permutation_cnt = 1; permutation_cnt <= MAX_PERMUTATION_NUMBER; ++permutation_cnt)
	{
		//cout << "permutation " << permutation_cnt << endl;
		shuffle_index_array(group_assign_index, global_total_sample_num);
		
// 		for (sample_loop = 1; sample_loop <= global_total_sample_num; ++sample_loop) {debugfile_shuffle_index << group_assign_index[sample_loop] << "\t";}
// 		debugfile_shuffle_index << endl;

		//calculate the statistic d
		for (asm_loop = 1; asm_loop < sortList_ASM.size(); ++asm_loop)
		{
			cur_dimension = sortList_ASM[asm_loop]->pathNum;
			for (group_loop = 1; group_loop <= global_total_group_num; ++group_loop)
				group_mean_prop_permuted[group_loop] = new double [cur_dimension+1];

			calculate_group_mean_prop_permutated(sortList_ASM[asm_loop]->sample_path_prop, group_mean_prop_permuted, global_group_size, group_assign_index, global_total_group_num, cur_dimension);

//			permuted_value_array[asm_loop] = calculate_single_gene_d(sortList_ASM[asm_loop]->sample_path_prop, group_mean_prop_permuted, sortList_ASM[asm_loop]->grand_mean_path_prop, cur_dimension, global_group_size, group_assign_index,
//				global_total_group_num, global_total_sample_num, s0 + sortList_ASM[asm_loop]->statistic_s_expr, -1);
			permuted_value_array[asm_loop] = calculate_single_gene_d(sortList_ASM[asm_loop]->sample_path_prop, group_mean_prop_permuted, sortList_ASM[asm_loop]->grand_mean_path_prop, cur_dimension, global_group_size, group_assign_index,
				global_total_group_num, global_total_sample_num, s0 + sortList_ASM[asm_loop]->statistic_s_expr, sortList_ASM[asm_loop]->statistic_s);

			for (group_loop = 1; group_loop <= global_total_group_num; ++group_loop)
				delete [] group_mean_prop_permuted[group_loop];
		}

		//sort and write the permuted values
		quicksort(permuted_value_array, sortList_ASM.size()-1);
		for (asm_loop = 1; asm_loop < sortList_ASM.size(); ++asm_loop)
			permuted_d[asm_loop][permutation_cnt] = permuted_value_array[asm_loop];
	}

	delete [] group_mean_prop_permuted;
	delete [] group_assign_index;
	delete [] permuted_value_array;

//	debugfile_shuffle_index.close();

	allpermutedvalueCnt = MAX_PERMUTATION_NUMBER * (sortList_ASM.size()-1);
	return permutation_cnt-1;
}

double compute_s_expr(double expr)
{
	//given an expression level, return the penalty parameter using a logistic function
	if (expr <= 1)
		return 1;
	//return 2.0 * (1.0 - 1.0 / (1.0 + exp(-1 * log10(expr)))); //P(-t)
	return 2.0 * (1.0 - 1.0 / (1.0 + exp(-1 * expr / 3))); //P(-t)
}

double compute_s0(long sampleSize)
{
	//compute value of s0
	double s0 = 0.0;
	unsigned long asm_loop;

	//compute s^\alpha for \alpha from 0 to 1.0 in step size = 0.05
	const int ALPHANUM = 21;
	int alphaLoopCnt;
	double alpha_values[50];
	
	for (asm_loop = 1; asm_loop < sortList_ASM.size(); ++asm_loop)
	{
		sortKey_ASM[asm_loop] = sortList_ASM[asm_loop]->statistic_s;
	}
	mergeSort_ASM_sort(sortList_ASM.size()-1);
	for (alphaLoopCnt = 1; alphaLoopCnt <= ALPHANUM; ++alphaLoopCnt)
	{
		//alpha_values[alphaLoopCnt] = sortKey_testGene[long(ceil((sortList_testGene_Num - 1) * (0.05 * (alphaLoopCnt - 1))) + 1)];
		alpha_values[alphaLoopCnt] = calculate_percentile(sortKey_ASM, sortList_ASM.size()-1, 0.05 * (alphaLoopCnt - 1));
	}


	//compute the 100 quantiles of the si values: q1, ... q100
	double quantiles[101];
	int quantileLoopCnt;

	for (quantileLoopCnt = 1; quantileLoopCnt <= 100; ++quantileLoopCnt)
	{
		//quantiles[quantileLoopCnt] = sortKey_testGene[long(ceil((sortList_testGene_Num - 1) * (0.01 * quantileLoopCnt)) + 1)];
		quantiles[quantileLoopCnt] = calculate_percentile(sortKey_ASM, sortList_ASM.size()-1, 0.01 * quantileLoopCnt);
	}


	//for every alpha, compute the coefficient of variation of the vj values
	double vj_values[101], cv_values[ALPHANUM+1];
	double* di_alpha = new double [sortList_ASM.size()];
	long di_alpha_cnt;
	for (alphaLoopCnt = 1; alphaLoopCnt <= ALPHANUM; ++alphaLoopCnt)
	{
		asm_loop = 1;

		//compute vj
		for (quantileLoopCnt = 1; quantileLoopCnt <= 100; ++quantileLoopCnt)
		{
			di_alpha_cnt = 0;
			//compute d_i^alpha
			for (; asm_loop < sortList_ASM.size() && sortList_ASM[asm_loop]->statistic_s < quantiles[quantileLoopCnt]; ++asm_loop)
			{				
				di_alpha[++di_alpha_cnt] = calculate_single_gene_d(sortList_ASM[asm_loop]->sample_path_prop, sortList_ASM[asm_loop]->group_mean_path_prop, sortList_ASM[asm_loop]->grand_mean_path_prop, 
					sortList_ASM[asm_loop]->pathNum, global_group_size, NULL, global_total_group_num, global_total_sample_num, alpha_values[alphaLoopCnt] + sortList_ASM[asm_loop]->statistic_s_expr, sortList_ASM[asm_loop]->statistic_s);
			}

			if (di_alpha_cnt > 0)
			{
				vj_values[quantileLoopCnt] = calculate_mad(di_alpha, di_alpha_cnt) / 0.64;
			} 
			else
			{
				vj_values[quantileLoopCnt] = 0.0;
			}
		}

		//compute coefficient of variation of the vj values
		cv_values[alphaLoopCnt] = calculate_coefficient_variation(vj_values, 100);
	}
	delete [] di_alpha;

	//find the min of the c_v's
	double min_cv = MAX_NUMBER;
	int min_cv_index = 1; 
	for (alphaLoopCnt = 1; alphaLoopCnt <= ALPHANUM; ++alphaLoopCnt)
	{
		if (cv_values[alphaLoopCnt] < min_cv)
		{
			min_cv = cv_values[alphaLoopCnt];
			min_cv_index = alphaLoopCnt;
		}
	}
	s0 = alpha_values[min_cv_index];

	return s0;
}



//output permuted d values
void output_permuted_d(string resultPath, long permutationCnt, double* permuted_d[], string resultFileSuffix)
{
	string outputfilename;
	ofstream outputfile;
	long totalCnt = 0; //should be testGeneNum * permutationCnt

	//output all permuted d values
	outputfilename = resultPath + "/stat/asm/d_permutation_" + resultFileSuffix + ".txt";
	outputfile.open(outputfilename.c_str());
	for (unsigned long geneLoopCnt = 1; geneLoopCnt < sortList_ASM.size(); ++geneLoopCnt)
	{
		for (long loopCnt = 1; loopCnt <= permutationCnt; ++loopCnt)
		{
			++totalCnt;
			outputfile << permuted_d[geneLoopCnt][loopCnt] << endl; 
		}
	}
	outputfile.close();
	
	//output count of permuted values, genes, and permutations
	outputfilename = resultPath + "/stat/asm/d_permutation_cnt_" + resultFileSuffix + ".txt";
	outputfile.open(outputfilename.c_str());
	outputfile << totalCnt << "\t" << sortList_ASM.size()-1 << "\t" << permutationCnt << endl;
	outputfile.close();

	return;
}


//sampleSize is the number of samples in each group, which should be the sum of #tech-replicates of every individual; also, assume sampleSize is fixed for each group, for the convenience of permutation test
bool diff_transcription_analysis_multi_class(double *statistics_expected_d, double* statistics_permuted_d[], string resultPath, string resultFileSuffix, double &s0)
{
	long ASMloopCnt, permutationCnt_actual, cnt_small_s = 0;;
		
	s0 = 0;
	double s_max = 0.0;

	//compute s
	for (ASMloopCnt = 1; ASMloopCnt < sortList_ASM.size(); ++ASMloopCnt)
	{
		sortList_ASM[ASMloopCnt]->statistic_s = calculate_single_gene_s_multi_class(sortList_ASM[ASMloopCnt]->sample_path_prop, sortList_ASM[ASMloopCnt]->group_mean_path_prop, sortList_ASM[ASMloopCnt]->pathNum, 
			global_group_size, NULL, global_total_group_num, global_total_sample_num);
		if (sortList_ASM[ASMloopCnt]->statistic_s > s_max)
			s_max = sortList_ASM[ASMloopCnt]->statistic_s;
		if (sortList_ASM[ASMloopCnt]->statistic_s < 1e-3)
		{
			sortList_ASM[ASMloopCnt]->statistic_s = 1e-3; //MAX_NUMBER/2;
			++cnt_small_s;
		}
	}
	//add expression term using min path expression
	for (ASMloopCnt = 1; ASMloopCnt < sortList_ASM.size(); ++ASMloopCnt)
	{
		sortList_ASM[ASMloopCnt]->statistic_s_expr = compute_s_expr(sortList_ASM[ASMloopCnt]->min_group_mean_expr) * s_max * 1.3; //5/9/2013 change min_group_mean_path_expr to min_group_mean_expr, allow low expression paths because they may correspond to isoform switch
	}
	//compute s0
	s0 = compute_s0(global_total_sample_num);
	cout << "s0 = " << s0 << "\t";

	//compute di
	for (ASMloopCnt = 1; ASMloopCnt < sortList_ASM.size(); ++ASMloopCnt)
	{
		sortList_ASM[ASMloopCnt]->statistic_d = calculate_single_gene_d(sortList_ASM[ASMloopCnt]->sample_path_prop, sortList_ASM[ASMloopCnt]->group_mean_path_prop, sortList_ASM[ASMloopCnt]->grand_mean_path_prop,
			sortList_ASM[ASMloopCnt]->pathNum, global_group_size, NULL, global_total_group_num, global_total_sample_num, s0 + sortList_ASM[ASMloopCnt]->statistic_s_expr, sortList_ASM[ASMloopCnt]->statistic_s);
	}


	//permutation 
	permutationCnt_actual = permutationTest(s0, resultPath, statistics_permuted_d, resultFileSuffix);

	
	//output permutation result
	output_permuted_d(resultPath, permutationCnt_actual, statistics_permuted_d, resultFileSuffix);

	//calculate the expected d
	for (ASMloopCnt = 1; ASMloopCnt < sortList_ASM.size(); ++ASMloopCnt)
	{
		statistics_expected_d[ASMloopCnt] = calculate_mean(statistics_permuted_d[ASMloopCnt], permutationCnt_actual);
	}

	return true;
}

