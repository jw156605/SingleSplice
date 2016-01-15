#include "difftrans_test.h"

double q0 = 0.0; //25% point of the permuted d values
double q50 = 0.0; //75% point of the permuted d values
double pi0 = 0.0; //pi_0, estimate for the proportion of true null genes

long allpermutedvalueCnt; //count for all permuted values, which is geneCnt * permutationCnt

double* statistics_expected_d;
double** statistics_permuted_d; //all permuted value

double* falseGeneCntList;
double significance_cutoff = 0; //cutoff for significant differences


/************************************************************************/
/* CALCULATE STATISTICS                                                 */
/************************************************************************/

void calculate_q0_q50(string resultPath, string resultFileSuffix)
{
	//find the 0% point and 50%point
	long index_q0, index_q50, loopCnt;
	string filename;	
	ifstream permutationValue_file;
	string info;
	filename = resultPath + "stat/asm/d_permutation_sorted_" + resultFileSuffix + ".txt";
	permutationValue_file.open(filename.c_str());

	index_q0 = long(ceil((allpermutedvalueCnt - 1) * 0.0) + 1);
	index_q50 = long(ceil((allpermutedvalueCnt - 1) * 0.50) + 1);

	for (loopCnt = 1; loopCnt < index_q0; ++loopCnt)
	{
		getline(permutationValue_file, info); 
	}
	permutationValue_file >> q0;

	for (++loopCnt; loopCnt < index_q50; ++loopCnt)
	{
		getline(permutationValue_file, info); 
	}
	permutationValue_file >> q50;

	permutationValue_file.close();

	return;
}



void calculate_pi0()
{
	long geneLoopCnt, diCnt = 0;

	for (geneLoopCnt = 1; geneLoopCnt < sortList_ASM.size(); ++geneLoopCnt)
	{
		if (sortList_ASM[geneLoopCnt]->statistic_d > q0 && sortList_ASM[geneLoopCnt]->statistic_d < q50)
		{
			++diCnt;
		}
	}

	pi0 = double(diCnt) / (0.5 * (sortList_ASM.size()-1));
	pi0 = pi0 > 1? 1 : pi0;
	pi0 = pi0 < 0.5? 0.5 : pi0;

	return;
}


void calculate_statistics(string resultPath, string resultFileSuffix)
{
	//	getCounts(resultPath);

	calculate_q0_q50(resultPath, resultFileSuffix);
#ifndef UNIX
	q0 = -0.5422;
	q50 = 0.5422;
#endif
	calculate_pi0();

	//cout << resultFileSuffix << ":\tq0 = " << q0 << "\tq50 = " << q50 << "\tpi0 = " << pi0 << endl;

	return;
}


/************************************************************************/
/* CALCULATE FDR                                                        */
/************************************************************************/

long count_falseGeneCnt_onePermutation(double delta, long permutationIndex)
{
	long falseGeneCnt = 0, geneLoopCnt;

	for (geneLoopCnt = 1; geneLoopCnt < sortList_ASM.size(); ++geneLoopCnt)
	{
		if (fabs(statistics_permuted_d[geneLoopCnt][permutationIndex] - statistics_expected_d[geneLoopCnt]) > delta)
		{
			//falsely called gene
			++falseGeneCnt;
		}
	}

	return falseGeneCnt;
}

void calculate_falseGeneCnt(double delta, double &falseGeneCnt_median, double &falseGeneCnt_90percentile, double &falseGeneCnt_mean, string resultPath)
{
	for (long permutationLoopCnt = 1; permutationLoopCnt <= MAX_PERMUTATION_NUMBER; ++permutationLoopCnt)
	{
		falseGeneCntList[permutationLoopCnt] = count_falseGeneCnt_onePermutation(delta, permutationLoopCnt);
	}

	quicksort(falseGeneCntList, MAX_PERMUTATION_NUMBER);

	falseGeneCnt_median = calculate_percentile(falseGeneCntList, MAX_PERMUTATION_NUMBER, 0.5);
	falseGeneCnt_90percentile = calculate_percentile(falseGeneCntList, MAX_PERMUTATION_NUMBER, 0.9);
	falseGeneCnt_mean = calculate_mean(falseGeneCntList, MAX_PERMUTATION_NUMBER);

	//	cout << falseGeneCnt_mean << "\t" << falseGeneCnt_median << '\t' <<falseGeneCnt_90percentile << "\t" << falseGeneCntList[permutationCnt] << '\t';
	//	for (permutationLoopCnt = 1; permutationLoopCnt <= permutationCnt; ++permutationLoopCnt)
	//	{
	//		cout << falseGeneCntList[permutationLoopCnt] << "\t";
	//	}

	return;
}

//given a cutoff delta, count the number of significant calls
long count_significance(double delta)
{
	long geneLoopCnt, significanceCnt = 0;

	for (geneLoopCnt = 1; geneLoopCnt < sortList_ASM.size(); ++geneLoopCnt)
	{
		if (fabs(sortList_ASM[geneLoopCnt]->statistic_d - statistics_expected_d[geneLoopCnt]) > delta)
		{
			//significant gene
			++significanceCnt;
		}
	}

	return significanceCnt;
}

long set_significance(double delta)
{
	long geneLoopCnt, significanceCnt = 0;

	for (geneLoopCnt = 1; geneLoopCnt < sortList_ASM.size(); ++geneLoopCnt)
	{
		if (fabs(sortList_ASM[geneLoopCnt]->statistic_d - statistics_expected_d[geneLoopCnt]) > delta)
		{
			//significant gene
			sortList_ASM[geneLoopCnt]->significant_stat = true;
			++significanceCnt;
		}
	}

	return significanceCnt;
}

void calculate_FDR(double delta, long significanceCnt, double &FDR_median, double &FDR_90percentile, double &FDR_mean, string resultPath)
{
	double falseGeneCnt_median = 0.0, falseGeneCnt_90percentile = 0.0, falseGeneCnt_mean = 0.0;

	calculate_falseGeneCnt(delta, falseGeneCnt_median, falseGeneCnt_90percentile, falseGeneCnt_mean, resultPath);

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


void difftrans_test_multi_class(string resultPath, string resultFileSuffix)
{
	string filename, comd;
	double delta = 0.0, FDR_median = 0.0, FDR_90percentile = 0.0, FDR_mean = 0.0, s0 = 0;
	long significanceCnt = 0, geneLoopCnt, group;
	int groupLoopCnt, sampleLoopCnt, distLoopCnt, curDimension;


	statistics_permuted_d = new double* [sortList_ASM.size()];
	for (geneLoopCnt = 1; geneLoopCnt < sortList_ASM.size(); ++geneLoopCnt)
		statistics_permuted_d[geneLoopCnt] = new double [MAX_PERMUTATION_NUMBER + 1];
	statistics_expected_d = new double [sortList_ASM.size()];

	/************************************************************************/
	/* Differential analysis                                                */
	/************************************************************************/
		
	diff_transcription_analysis_multi_class(statistics_expected_d, statistics_permuted_d, resultPath, resultFileSuffix, s0);
	cout << resultFileSuffix << ":\t#asm = " << sortList_ASM.size()-1 << "\t#permutation = " << MAX_PERMUTATION_NUMBER << "\t#allpermutedvalue = " << allpermutedvalueCnt << endl;
	
#ifdef UNIX
	comd = "sort -n +0 -1 " + resultPath + "stat/asm/d_permutation_" + resultFileSuffix + ".txt > " + resultPath + "stat/asm/d_permutation_sorted_" + resultFileSuffix + ".txt";
	system(comd.c_str());
#endif

	calculate_statistics(resultPath, resultFileSuffix);


	/************************************************************************/
	/* Sort All Genes                                                       */
	/************************************************************************/
	for (geneLoopCnt = 1; geneLoopCnt < sortList_ASM.size(); ++geneLoopCnt)
		sortKey_ASM[geneLoopCnt] = sortList_ASM[geneLoopCnt]->statistic_d;
	mergeSort_ASM_sort(sortList_ASM.size()-1);

	for (geneLoopCnt = 1; geneLoopCnt < sortList_ASM.size(); ++geneLoopCnt)
		sortList_ASM[geneLoopCnt]->statistic_d_expected = statistics_expected_d[geneLoopCnt];

	/************************************************************************/
	/* Compute False Discovery Rate                                         */
	/************************************************************************/
	//stack = new long [permutationCnt+1];
	falseGeneCntList = new double [MAX_PERMUTATION_NUMBER+1];

	filename = resultPath + "FDR_transcription_" + resultFileSuffix + ".txt";
	ofstream FDRfile(filename.c_str());
	FDRfile << "delta is the significance cutoff of the differential transcription statistics (calculated as |stat_d_expected-stat_d|)" << endl; 
	FDRfile << "shifting delta will call different sets of significant ASMs and will result in different FDRs" << endl;
	FDRfile << "you may choose your desired FDR according to this list" << endl << endl;
	for (delta = 0; delta <= 10; delta += 0.01)
	{
		significanceCnt = count_significance(delta);
		calculate_FDR(delta, significanceCnt, FDR_median, FDR_90percentile, FDR_mean, resultPath);

		FDRfile << "delta = " << delta << ": " << significanceCnt << " ASMs picked, with FDR(median) = "
			<< FDR_median << " and FDR(mean) = " << FDR_mean << " and FDR(90percentile) = " << FDR_90percentile << endl;

		//cout << "delta = " << delta << ": " << significanceCnt << " ASMs picked, with FDR(median) = "
		//	<< FDR_median << " and FDR(mean) = " << FDR_mean << " and FDR(90percentile) = " << FDR_90percentile << endl;

		if (significance_cutoff < 1e-5 && FDR_median <= false_discovery_rate)
		{
			significance_cutoff = delta;
		}
		if (significanceCnt <= 0)
			break;
	}
	FDRfile.close();

	if (significance_cutoff < 1e-5)
	{
		significance_cutoff = delta;
	}
	set_significance(significance_cutoff);


	/************************************************************************/
	/* Pairwise differential tests                                          */
	/************************************************************************/
	diff_transcription_pairwise_test(statistics_expected_d, statistics_permuted_d, falseGeneCntList, s0, pi0);


	/************************************************************************/
	/* Clean up allocated arrays                                            */
	/************************************************************************/
	for (long geneLoopCnt = 1; geneLoopCnt < sortList_ASM.size(); ++geneLoopCnt)
		delete [] statistics_permuted_d[geneLoopCnt];
	delete [] falseGeneCntList;
	delete [] statistics_expected_d;
	delete [] statistics_permuted_d;

	return;
}


