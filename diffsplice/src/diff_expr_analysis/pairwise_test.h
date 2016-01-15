#include "config.h"
#include "calc_function.h"

#ifndef PAIRWISE_TEST
#define PAIRWISE_TEST


class test_ASM
{
public:
	long orig_ASM_id;
	pairdiff *statistic;

	test_ASM();
};


void diff_transcription_pairwise_test(double *statistics_expected_d, double* statistics_permuted_d[], double* falseGeneCntList, double s0, double pi0);


#endif
