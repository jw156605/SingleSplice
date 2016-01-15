#include "config.h"
#include "asm_def.h"
#include "calc_function.h"
#include "permutation_test.h"
#include "pairwise_test.h"


#ifndef DIFFTRANS_TEST
#define DIFFTRANS_TEST

//permutation count
extern long allpermutedvalueCnt; //count for all permuted values, which is geneCnt * permutationCnt
void difftrans_test_multi_class(string resultPath, string resultFileSuffix);

#endif


