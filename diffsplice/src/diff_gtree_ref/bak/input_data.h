
#ifndef INPUT_DATA
#define INPUT_DATA

#include <unordered_map>

#include "config.h"
#include "splice_graph.h"

void input_data(RangeJunctionList *list_fragments, string dir_input_frag, string dir_junction_annotation, long chromosome_length);
bool input_data_ref(RangeJunctionList *list_fragments, string dir_input_frag, string dir_junction_annotation, string dir_exon_annotation, long chromosome_length);

extern ofstream outfile_junction_all; //bed track of all splice junctions
extern ofstream outfile_junction_filtered; //bed track of splice junctions after filtering
extern ofstream outfile_debuginfo;

#endif

