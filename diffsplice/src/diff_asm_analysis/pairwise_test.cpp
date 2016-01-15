#include "difftrans_test.h"
#include "pairwise_test.h"


vector <test_ASM*> list_test_asm;

test_ASM::test_ASM()
{
	orig_ASM_id = 0;
	statistic = NULL;
}



/************************************************************************/
/* Statistic for a single gene                                          */
/************************************************************************/

double calculate_single_gene_r_two_class(double *meanProp_group1, double *meanProp_group2, int dimension)
{
	return sqrt(calculate_JSD(meanProp_group1, meanProp_group2, dimension));	
}

double calculate_single_gene_s_two_class(double** individualProp_group1, double* meanProp_group1, int n1, double** individualProp_group2, double* meanProp_group2, int n2, int dimension)
{
	if (n1 + n2 <= 2)
		return 0;

	double SumJSD_group1 = 0.0, SumJSD_group2 = 0.0;
	int sampleLoopCnt;

	for (sampleLoopCnt = 1; sampleLoopCnt <= n1; ++sampleLoopCnt)
	{
		SumJSD_group1 += calculate_JSD(individualProp_group1[sampleLoopCnt], meanProp_group1, dimension);
	}
	for (sampleLoopCnt = 1; sampleLoopCnt <= n2; ++sampleLoopCnt)
	{
		SumJSD_group2 += calculate_JSD(individualProp_group2[sampleLoopCnt], meanProp_group2, dimension);
	}

	return (1.0/n1 + 1.0/n2) * (SumJSD_group1 + SumJSD_group2) / (n1+n2-2);
}

double calculate_single_gene_d_two_class(double** individualProp_group1, double* meanProp_group1, int n1, double** individualProp_group2, double* meanProp_group2, int n2, int dimension, double s0, double value_s)
{
	if (value_s < 0)
		//value_s has not been calculated, calculate now
		value_s = calculate_single_gene_s_two_class(individualProp_group1, meanProp_group1, n1, individualProp_group2, meanProp_group2, n2, dimension);

	if (value_s + s0 < 1e-15)
	{
		//variance is too small, do not count the difference
		//cout << "dividing a number ~= 0 in calculate_single_gene_d" << endl;
		return 1e-15;
	}
	return calculate_single_gene_r_two_class(meanProp_group1, meanProp_group2, dimension) / (sqrt(value_s) + sqrt(s0));
}


/************************************************************************/
/* Statistic for all genes                                              */
/************************************************************************/


void calculate_test_group_mean_prop(double** sample_dist, double** group_mean_prop, int group_size_1, int group_size_2, int* group_assign, int dimension)
{
	int pathLoopCnt, sampleLoopCnt, groupLoopCnt;
	double sum;

	for (pathLoopCnt = 1; pathLoopCnt <= dimension; ++pathLoopCnt)
	{
		sum = 0.0;
		for (sampleLoopCnt = 1; sampleLoopCnt <= group_size_1; ++sampleLoopCnt)
		{
			sum += sample_dist[group_assign[sampleLoopCnt]][pathLoopCnt];
		}
		group_mean_prop[1][pathLoopCnt] = sum / group_size_1;

		sum = 0.0;
		for (sampleLoopCnt = group_size_1+1; sampleLoopCnt <= group_size_1+group_size_2; ++sampleLoopCnt)
		{
			sum += sample_dist[group_assign[sampleLoopCnt]][pathLoopCnt];
		}
		group_mean_prop[2][pathLoopCnt] = sum / group_size_2;
	}

	//normalize the proportions 
	for (groupLoopCnt = 1; groupLoopCnt <= 2; ++groupLoopCnt)
	{
		sum = 0.0;
		for (pathLoopCnt = 1; pathLoopCnt <= dimension; ++pathLoopCnt)
			sum += group_mean_prop[groupLoopCnt][pathLoopCnt];
		if (sum == 0.0)
			sum = 1.0;
		for (pathLoopCnt = 1; pathLoopCnt <= dimension; ++pathLoopCnt)
			group_mean_prop[groupLoopCnt][pathLoopCnt] = group_mean_prop[groupLoopCnt][pathLoopCnt] / sum;
	}	

	return;
}


long permutationTest_two_class(double s0, double* permuted_d[], int group_index_1, int group_index_2)
{
	if (global_total_group_num < 2)
	{
		//one group test makes no sense
		return 0;
	}

	//	ofstream debugfile_shuffle_index("debug_shuffled_index.txt");

	long permutation_cnt, asm_loop, group_loop, sample_loop;
	long group_size_1, group_size_2;
	int cur_dimension, *group_assign_index, cur_orig_id;
	const int test_group_num = 2;

	group_size_1 = global_group_size[group_index_1];
	group_size_2 = global_group_size[group_index_2];
	
	group_assign_index = new int [group_size_1 + group_size_2 + 1];
	for (sample_loop = 1; sample_loop <= group_size_1; ++sample_loop) {group_assign_index[sample_loop] = global_group_base_index[group_index_1] + sample_loop;}
	for (sample_loop = 1; sample_loop <= group_size_2; ++sample_loop) {group_assign_index[group_size_1 + sample_loop] = global_group_base_index[group_index_2] + sample_loop;}

	double **group_mean_prop_permuted, *permuted_value_array;

	group_mean_prop_permuted = new double* [test_group_num+1];
	permuted_value_array = new double [list_test_asm.size()];

	for (permutation_cnt = 1; permutation_cnt <= MAX_PERMUTATION_NUMBER; ++permutation_cnt)
	{
		shuffle_index_array(group_assign_index, group_size_1+group_size_2);

		//calculate the statistic d
		for (asm_loop = 1; asm_loop < list_test_asm.size(); ++asm_loop)
		{
			cur_orig_id = list_test_asm[asm_loop]->orig_ASM_id;
			cur_dimension = sortList_ASM[cur_orig_id]->pathNum;

			for (group_loop = 1; group_loop <= test_group_num; ++group_loop)
				group_mean_prop_permuted[group_loop] = new double [cur_dimension+1];

			calculate_test_group_mean_prop(sortList_ASM[cur_orig_id]->sample_path_prop, group_mean_prop_permuted, group_size_1, group_size_2, group_assign_index, cur_dimension);

			permuted_value_array[asm_loop] = calculate_single_gene_d_two_class(NULL, group_mean_prop_permuted[1], group_size_1, NULL, group_mean_prop_permuted[2], group_size_2, cur_dimension,
				s0 + sortList_ASM[cur_orig_id]->statistic_s_expr, sortList_ASM[cur_orig_id]->statistic_s);

			for (group_loop = 1; group_loop <= test_group_num; ++group_loop)
				delete [] group_mean_prop_permuted[group_loop];
		}

		//sort and write the permuted values
		quicksort(permuted_value_array, list_test_asm.size()-1);
		for (asm_loop = 1; asm_loop < list_test_asm.size(); ++asm_loop)
			permuted_d[asm_loop][permutation_cnt] = permuted_value_array[asm_loop];
	}

	delete [] group_mean_prop_permuted;
	delete [] group_assign_index;
	delete [] permuted_value_array;

	//	debugfile_shuffle_index.close();

	allpermutedvalueCnt = MAX_PERMUTATION_NUMBER * (list_test_asm.size()-1);
	return permutation_cnt-1;
}


/************************************************************************/
/* CALCULATE FDR                                                        */
/************************************************************************/

long count_falseGeneCnt_onePermutation_pairwise(double delta, long permutationIndex, double *statistics_expected_d, double *statistics_permuted_d[])
{
	long falseGeneCnt = 0, asm_loop;

	for (asm_loop = 1; asm_loop < list_test_asm.size(); ++asm_loop)
	{
		if (fabs(statistics_permuted_d[asm_loop][permutationIndex] - statistics_expected_d[asm_loop]) > delta)
		{
			//falsely called gene
			++falseGeneCnt;
		}
	}

	return falseGeneCnt;
}


void calculate_falseGeneCnt_pairwise(double delta, double &falseGeneCnt_median, double &falseGeneCnt_90percentile, double &falseGeneCnt_mean, double *falseGeneCntList, double *statistics_expected_d, double *statistics_permuted_d[])
{
	for (long permutationLoopCnt = 1; permutationLoopCnt <= MAX_PERMUTATION_NUMBER; ++permutationLoopCnt)
	{
		falseGeneCntList[permutationLoopCnt] = count_falseGeneCnt_onePermutation_pairwise(delta, permutationLoopCnt, statistics_expected_d, statistics_permuted_d);
	}

	quicksort(falseGeneCntList, MAX_PERMUTATION_NUMBER);

	falseGeneCnt_median = calculate_percentile(falseGeneCntList, MAX_PERMUTATION_NUMBER, 0.5);
	falseGeneCnt_90percentile = calculate_percentile(falseGeneCntList, MAX_PERMUTATION_NUMBER, 0.9);
	falseGeneCnt_mean = calculate_mean(falseGeneCntList, MAX_PERMUTATION_NUMBER);

	return;
}

//given a cutoff delta, count the number of significant calls
long count_significance_pairwise(double delta)
{
	long asm_loop, significanceCnt = 0;

	for (asm_loop = 1; asm_loop < list_test_asm.size(); ++asm_loop)
	{
		if (fabs(list_test_asm[asm_loop]->statistic->stat_d - list_test_asm[asm_loop]->statistic->stat_d_expected) > delta)
		{
			//significant 
			++significanceCnt;
		}
	}

	return significanceCnt;
}

long set_significance_pairwise(double delta, int group_index_1, int group_index_2)
{
	long asm_loop, significanceCnt = 0;

	for (asm_loop = 1; asm_loop < sortList_ASM.size(); ++asm_loop)
	{
		list_test_asm[asm_loop]->statistic->stat_difftrans = fabs(list_test_asm[asm_loop]->statistic->stat_d - list_test_asm[asm_loop]->statistic->stat_d_expected);
		if (list_test_asm[asm_loop]->statistic->stat_difftrans > delta)
		{
			//significant gene
			list_test_asm[asm_loop]->statistic->significant_stat = true;

			double cur_JSD = calculate_JSD(sortList_ASM[list_test_asm[asm_loop]->orig_ASM_id]->group_mean_path_prop[group_index_1], sortList_ASM[list_test_asm[asm_loop]->orig_ASM_id]->group_mean_path_prop[group_index_2], 
				sortList_ASM[list_test_asm[asm_loop]->orig_ASM_id]->pathNum);
			list_test_asm[asm_loop]->statistic->jsd = cur_JSD;

			if (sqrt(cur_JSD) >= thresh_JSD)// && list_test_asm[asm_loop]->statistic->stat_difftrans >= thresh_stat_d)
			{
				list_test_asm[asm_loop]->statistic->significant_both = true;
			}

			++significanceCnt;
		}
	}

	return significanceCnt;
}

void calculate_FDR_pairwise(double delta, long significanceCnt, double pi0, double &FDR_median, double &FDR_90percentile, double &FDR_mean, double *falseGeneCntList, double *statistics_expected_d, double *statistics_permuted_d[])
{
	double falseGeneCnt_median = 0.0, falseGeneCnt_90percentile = 0.0, falseGeneCnt_mean = 0.0;

	calculate_falseGeneCnt_pairwise(delta, falseGeneCnt_median, falseGeneCnt_90percentile, falseGeneCnt_mean, falseGeneCntList, statistics_expected_d, statistics_permuted_d);

	falseGeneCnt_median = falseGeneCnt_median * pi0;
	falseGeneCnt_90percentile = falseGeneCnt_90percentile * pi0;
	falseGeneCnt_mean = falseGeneCnt_mean * pi0;

	if (significanceCnt != 0)
	{
		FDR_median = falseGeneCnt_median / significanceCnt;
		FDR_90percentile = falseGeneCnt_90percentile / significanceCnt;
		FDR_mean = falseGeneCnt_mean / significanceCnt;
	}
	else
	{
		FDR_median = 0.0;
		FDR_90percentile = 0.0;
		FDR_mean = 0.0;
	}

	return;
}


//sampleSize is the number of samples in each group, which should be the sum of #tech-replicates of every individual; also, assume sampleSize is fixed for each group, for the convenience of permutation test
void diff_transcription_test_two_class(double *statistics_expected_d, double* statistics_permuted_d[], double* falseGeneCntList, int group_index_1, int group_index_2, double s0, double pi0)
{
	long ASMloopCnt, permutationCnt_actual;
	pairdiff *new_pairdiff;

	for (ASMloopCnt = 1; ASMloopCnt < list_test_asm.size(); ++ASMloopCnt)
	{
		new_pairdiff = new pairdiff;
		
		int cur_orig_id = list_test_asm[ASMloopCnt]->orig_ASM_id;
		new_pairdiff->stat_d = calculate_single_gene_d_two_class(NULL, sortList_ASM[cur_orig_id]->group_mean_path_prop[group_index_1], global_group_size[group_index_1],
			NULL, sortList_ASM[cur_orig_id]->group_mean_path_prop[group_index_2], global_group_size[group_index_2], sortList_ASM[cur_orig_id]->pathNum,
			s0 + sortList_ASM[cur_orig_id]->statistic_s_expr, sortList_ASM[cur_orig_id]->statistic_s);
		new_pairdiff->groupindex_1 = group_index_1;
		new_pairdiff->groupindex_2 = group_index_2;
		
		list_test_asm[ASMloopCnt]->statistic = new_pairdiff;
	}
	
	//permutation 
	permutationCnt_actual = permutationTest_two_class(s0, statistics_permuted_d, group_index_1, group_index_2);
		
	//calculate the expected d
	for (ASMloopCnt = 1; ASMloopCnt < list_test_asm.size(); ++ASMloopCnt)
	{
		statistics_expected_d[ASMloopCnt] = calculate_mean(statistics_permuted_d[ASMloopCnt], permutationCnt_actual);
	}

	/************************************************************************/
	/* Sort All Genes                                                       */
	/************************************************************************/
	double *sort_test_asm_key = new double[list_test_asm.size()];
	void **sort_test_asm_list = new void * [list_test_asm.size()];

	for (ASMloopCnt = 1; ASMloopCnt < list_test_asm.size(); ++ASMloopCnt)
	{
		sort_test_asm_key[ASMloopCnt] = list_test_asm[ASMloopCnt]->statistic->stat_d;
		sort_test_asm_list[ASMloopCnt] = (void*) list_test_asm[ASMloopCnt];
	}
	mergeSort(sort_test_asm_list, sort_test_asm_key, list_test_asm.size()-1);

	for (ASMloopCnt = 1; ASMloopCnt < list_test_asm.size(); ++ASMloopCnt)
	{
		list_test_asm[ASMloopCnt] = (test_ASM*) sort_test_asm_list[ASMloopCnt];
		list_test_asm[ASMloopCnt]->statistic->stat_d_expected = statistics_expected_d[ASMloopCnt];
	}


	/************************************************************************/
	/* Compute False Discovery Rate                                         */
	/************************************************************************/
	double delta, sig_cutoff = 0.0, FDR_median = 0.0, FDR_90percentile = 0.0, FDR_mean = 0.0;
	long significantCnt;

	for (delta = 0; delta <= 10; delta += 0.01)
	{
		significantCnt = count_significance_pairwise(delta);
		calculate_FDR_pairwise(delta, significantCnt, pi0, FDR_median, FDR_90percentile, FDR_mean, falseGeneCntList, statistics_expected_d, statistics_permuted_d);

		if (sig_cutoff < 1e-5 && FDR_median <= false_discovery_rate)
		{
			sig_cutoff = delta;
			break;
		}
		if (significantCnt <= 0)
			break;
	}
	if (sig_cutoff < 1e-5)
	{
		sig_cutoff = MAX_NUMBER;
	}
	set_significance_pairwise(sig_cutoff, group_index_1, group_index_2);


	for (ASMloopCnt = 1; ASMloopCnt < list_test_asm.size(); ++ASMloopCnt)
	{
		sortList_ASM[list_test_asm[ASMloopCnt]->orig_ASM_id]->list_pairdiff.push_back(list_test_asm[ASMloopCnt]->statistic);
		list_test_asm[ASMloopCnt]->statistic = NULL;
	}

	for (ASMloopCnt = 1; ASMloopCnt < list_test_asm.size(); ++ASMloopCnt)
	{
		sort_test_asm_key[ASMloopCnt] = list_test_asm[ASMloopCnt]->orig_ASM_id;
		sort_test_asm_list[ASMloopCnt] = (void*) list_test_asm[ASMloopCnt];
	}
	mergeSort(sort_test_asm_list, sort_test_asm_key, list_test_asm.size()-1);

	for (ASMloopCnt = 1; ASMloopCnt < list_test_asm.size(); ++ASMloopCnt)
	{
		list_test_asm[ASMloopCnt] = (test_ASM*) sort_test_asm_list[ASMloopCnt];
	}

	delete [] sort_test_asm_key;
	delete [] sort_test_asm_list;

	return;
}



//test the differences between pairs of groups
void diff_transcription_pairwise_test(double *statistics_expected_d, double* statistics_permuted_d[], double* falseGeneCntList, double s0, double pi0)
{
	test_ASM *new_asm;
	
	long asm_loop, num_comparison = 0, pair_loop;
	int groupCnt_1, groupCnt_2;

	list_test_asm.push_back(NULL);
	for (asm_loop = 1; asm_loop < sortList_ASM.size(); ++asm_loop)
	{
		new_asm = new test_ASM;
		new_asm->orig_ASM_id = asm_loop;
		list_test_asm.push_back(new_asm);
	}
	
	for (groupCnt_1 = 1; groupCnt_1 < global_total_group_num; ++groupCnt_1)
	{
		for (groupCnt_2 = groupCnt_1 + 1; groupCnt_2 <= global_total_group_num; ++groupCnt_2)
		{
			diff_transcription_test_two_class(statistics_expected_d, statistics_permuted_d, falseGeneCntList, groupCnt_1, groupCnt_2, s0, pi0);
			++num_comparison;
		}
	}

	for (asm_loop = 1; asm_loop < list_test_asm.size(); ++asm_loop)
	{
		delete list_test_asm[asm_loop];
	}
	list_test_asm.clear();


	double *sort_pairdiff_key = new double [num_comparison+1];
	void **sort_pairdiff_list = new void * [num_comparison+1];

	for (asm_loop = 1; asm_loop < sortList_ASM.size(); ++asm_loop)
	{
		for (pair_loop = 1; pair_loop <= num_comparison; ++pair_loop)
		{
			sort_pairdiff_key[pair_loop] = - sortList_ASM[asm_loop]->list_pairdiff[pair_loop-1]->stat_difftrans;
			sort_pairdiff_list[pair_loop] = (void*) sortList_ASM[asm_loop]->list_pairdiff[pair_loop-1];
		}
		mergeSort(sort_pairdiff_list, sort_pairdiff_key, num_comparison);

		for (pair_loop = 1; pair_loop <= num_comparison; ++pair_loop)
		{
			sortList_ASM[asm_loop]->list_pairdiff[pair_loop-1] = (pairdiff*) sort_pairdiff_list[pair_loop];
		}
	}

	delete [] sort_pairdiff_key;
	delete [] sort_pairdiff_list;

	return;
}

