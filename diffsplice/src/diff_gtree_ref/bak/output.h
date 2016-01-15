
#ifndef OUTPUT
#define OUTPUT

const int unknown = 0;
const int exon_skipping = 1;
const int mutual_exclusive = 2;
const int intron_retention = 3;
const int alter_splice_site = 4;
const int alter_start = 5;
const int alter_end = 6;
const int diff_expression = 7;

void output_ASMpath_gtf(GTvertex *targetVertex, RangeJunctionList *pathList); //output ASM paths in GTF format, used in constructGTree
void output(GTvertex *rootVertex);
void initialize_outputfiles();
void close_outputfiles();


#endif

