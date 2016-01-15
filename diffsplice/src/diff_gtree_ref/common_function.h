
#ifndef COMMON_FUNCTION
#define COMMON_FUNCTION

#include "config.h"

string itostr(long t); //convert a long integer to string
void mergeSort(void** sortlist, double* sortkey, long sortList_size); //sort a list in ascending order of sort key
void quicksort(double *sortArray, long length, long *quicksortStack);
double calculate_mean(double *dataarray, long length);

#endif

