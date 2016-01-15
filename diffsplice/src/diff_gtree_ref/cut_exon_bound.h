#include "config.h"

#ifndef CUT_EXON_BOUND
#define CUT_EXON_BOUND

#ifdef UNIX
bool refine_exon_bound(vector<double> &signal, int num_cut_pt, int exon_code, vector<long> &cut_pt, vector<double> &avg_coverage);
#endif

#endif


