#include "config.h"
#include "gene_def.h"

#ifndef CALC_FUNCTION
#define CALC_FUNCTION

extern double* sortKey_ASM;

double calculate_JSD(double *array_P, double *array_Q, long length);
void quicksort(double *sortArray, long length);
double calculate_mean(double *dataarray, long length);
void calculate_meanDistribution(double** distributions, double* meanDist, int distNum, int dimension);
void calculate_group_mean_expr_permutated(double* sample_expr, double* group_mean_expr, int* group_size, int* group_assign, int group_num);
double calculate_percentile(double *dataarray, long length, double percentile_0to1);
double calculate_SSE(double *dataarray, long length);
double calculate_mad(double *dataArray, long arrayLength);
double calculate_coefficient_variation(double *dataArray, long arrayLength);
void shuffle_index_array(int *index_array, int array_length);

void mergeSort_ASM_sort(long sortList_size);
void mergeSort_initialization(long sortList_ASM_Num);
void mergeSort_cleanup();

void mergeSort(void** sortlist, double* sortkey, long sortList_size);

#endif

