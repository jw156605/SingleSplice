#include "config.h"
#include "splice_graph.h"
#include "common_function.h"
#include "output.h"

long exon_skipping_cnt = 0;
long mutual_exclusive_cnt = 0;
long intron_retention_cnt = 0;

//For collecting d statistics
vector <GTvertex*> vertexListForStatistics;

//main output
ofstream outfile_gtree; //the tree structure
ofstream outfile_stats_expr; //for differential expression analysis
ofstream outfile_stat_asm; //for differential transcription analysis
ofstream outfile_gtf_splice_graph; //gtf track of splice graph
ofstream outfile_gtf_asm_path; //gtf track of asm paths (first level decomposition)
ofstream outfile_gtf_asm_path_decomposed; //gtf track of asm paths AFTER decomposition
ofstream outfile_junction_all; //bed track of all splice junctions
ofstream outfile_junction_filtered; //bed track of splice junctions after filtering
ofstream outfile_asm_composition; //output composition of ASMs

//new version output, for later R process
ofstream outfile_stats_expr_new; //for differential expression analysis
ofstream outfile_stat_asm_new_absexpr; //for differential transcription analysis
ofstream outfile_stat_asm_new_proportion; //for differential transcription analysis
ofstream outfile_stat_asmexpr_new; //asm expression

ofstream outfile_debuginfo; // debug information for the lifetime of the program

/*//determine the category of an alternative splicing module
int alterSpliceCate(GTvertex *curVertex)
{
	//categorize alternative splicing events

	//if not, and the ASM has been assigned a category, then use that
	if (curVertex->ASMcategory > 0)
	{
		return curVertex->ASMcategory;
	}

	if ((curVertex->childType == 2 || curVertex->childType == 3) && curVertex->childNum > 1)
	{
		//ASM
	}
	else
	{
		return -1;
	}

	//SECOND, categorize the alternative splicing
	alternative_path *pathA, *pathB;

	if (curVertex->childType == 2 && curVertex->major_alter_paths_num > 1)
	{
		pathA = curVertex->major_alter_paths;
		while (pathA != NULL)
		{
			pathB = pathA->next;
			while (pathB != NULL)
			{
				if (abs(pathA->path_start - pathB->path_start) < 2 && abs(pathA->path_end - pathB->path_end) < 2)
				{
					if (pathA->exonNum == 0 && pathA->junctionNum == 1 && pathB->junctionNum >= 2 && pathB->exonNum >= 1 || pathA->junctionNum >= 2 && pathA->exonNum >= 1 && pathB->exonNum == 0 && pathB->junctionNum == 1)
					{
						curVertex->ASMcategory = exon_skipping;						
						return exon_skipping;
					}
					if (pathA->exonNum >= 1 && pathA->junctionNum >= 2 && pathB->exonNum >= 1 && pathB->junctionNum >= 2)
					{
						curVertex->ASMcategory = mutual_exclusive;						
					}
					if (curVertex->major_alter_paths_num == 2 && (pathA->exonNum == 1 && pathA->junctionNum == 0 && pathB->exonNum == 0 && pathB->junctionNum == 1 || pathA->exonNum == 0 && pathA->junctionNum == 1 && pathB->exonNum == 1 && pathB->junctionNum == 0))
					{
						curVertex->ASMcategory = intron_retention;
						return intron_retention;
					}
				}
				else if (abs(pathA->path_start - pathB->path_start) < 2)
				{
					curVertex->ASMcategory = alter_end;	
				}
				else if (abs(pathA->path_end - pathB->path_end) < 2)
				{
					curVertex->ASMcategory = alter_start;	
				}

				pathB = pathB->next;
			}

			pathA = pathA->next;
		}
	}

	if (curVertex->ASMcategory == mutual_exclusive || curVertex->ASMcategory == alter_end || curVertex->ASMcategory == alter_start)
	{
		//let mutual exclusive have lower priority than exon_skipping, because in the loop above we check children pair by pair
		return curVertex->ASMcategory;
	}

	// 	if (curVertex->childType == 2 && curVertex->childNum > 1)
	// 	{
	// 		edgeA = curVertex->child;
	// 		while (edgeA != NULL)
	// 		{
	// 			vertexA = edgeA->linkedVertex;
	// 
	// 			edgeB = edgeA->next;
	// 			while (edgeB != NULL)
	// 			{
	// 				vertexB = edgeB->linkedVertex;
	// 
	// 				if (vertexA->rangeLow == vertexB->rangeLow && vertexA->rangeHigh == vertexB->rangeHigh)
	// 				{
	// 					if (vertexA->childNum == 0 && vertexA->junctionNum == 1 && vertexB->childNum == 3 && vertexB->junctionNum == 2 || vertexA->childNum == 3 && vertexA->junctionNum == 2 && vertexB->childNum == 0 && vertexB->junctionNum == 1)
	// 					{
	// 						output_vertex(curVertex, &exon_skipping_file);
	// 						exon_skipping_file << endl;
	// 						return exon_skipping;
	// 					}
	// 					if (vertexA->childNum == 3 && vertexA->junctionNum == 2 && vertexB->childNum == 3 && vertexB->junctionNum == 2)
	// 					{
	// 						//if (vertexA->junctionInRange->list->junc->end > vertexB->junctionInRange->list->next->junc->start || vertexB->junctionInRange->list->junc->end > vertexA->junctionInRange->list->next->junc->start)
	// 						{
	// 							output_vertex(curVertex, &mutual_exclusive_file);
	// 							mutual_exclusive_file << endl;
	// 							return mutual_exclusive;
	// 						}
	// 					}
	// 				}
	// 
	// 				edgeB = edgeB->next;
	// 			}
	// 
	// 			edgeA = edgeA->next;
	// 		}
	// 	}

	curVertex->ASMcategory = unknown;
	return unknown;
}*/


void extract_exons(alternative_path *asm_path, vector<pair <long, long> > &exonList)
{
	if (!asm_path){
		exonList.clear();
		return;
	}

	if (asm_path->pathVertex->childType == 0 || asm_path->pathVertex->childType == 3)
	{
		rangeJunction *curJunc = asm_path->pathVertex->junctionInRange->list;
		while (curJunc != NULL)
		{
			if (curJunc->junc->type == frag_exon || curJunc->junc->type == frag_retained_intron)
				exonList.push_back(make_pair(curJunc->junc->start, curJunc->junc->end));

			curJunc = curJunc->next;
		}
	}
	else if (asm_path->pathVertex->childType == 1 || asm_path->pathVertex->childType == 2)
	{
		GTedge *path_component = asm_path->pathVertex->child;
		while (path_component != NULL)
		{
			GTvertex *cur_component = path_component->linkedVertex;
			if (cur_component->exonNum > 0)
				exonList.push_back(make_pair(cur_component->rangeLow, cur_component->rangeHigh));

			path_component = path_component->next;
		}
	}
}

// check whether two paths are disjoint (no overlapping exons)
bool overlapping_path(alternative_path *path_A, alternative_path *path_B)
{
	vector<pair <long, long> > exonList_A, exonList_B;
	extract_exons(path_A, exonList_A);
	if (exonList_A.empty()){
		return false;
	}
	extract_exons(path_B, exonList_B);
	if (exonList_B.empty()){
		return false;
	}

	for (vector<pair <long, long> >::iterator it_A = exonList_A.begin(); it_A != exonList_A.end(); ++it_A){
		for (vector<pair <long, long> >::iterator it_B = exonList_B.begin(); it_B != exonList_B.end(); ++it_B){
			if (it_A->second >= it_B->first && it_A->first <= it_B->second){
				return true;
			}
		}
	}
	return false;
}

// create a union list of exons from two lists, return the union list and the composition (1 - only A, 2 - only B, 3 - both)
void union_exons(vector<pair <long, long> > &exonList_A, vector<pair <long, long> > &exonList_B, vector<pair <long, long> > &exonList_union, vector<int> &inList) {
	exonList_union.clear();
	inList.clear();

	vector<pair <long, long> >::iterator it_A = exonList_A.begin(), it_B = exonList_B.begin();
	while (it_A != exonList_A.end() && it_B != exonList_B.end()) {
		// no two exons will overlap in gtree
		if (it_A->first < it_B->first) {
			exonList_union.push_back(make_pair(it_A->first, it_A->second));
			inList.push_back(1);
			++it_A;
		}
		else if (it_A->first > it_B->first) {
			exonList_union.push_back(make_pair(it_B->first, it_B->second));
			inList.push_back(2);
			++it_B;
		}
		else if (it_A->first == it_B->first && it_A->second == it_B->second){
			exonList_union.push_back(make_pair(it_A->first, it_A->second));
			inList.push_back(3);
			++it_A;
			++it_B;
		}
		else {
			cout << "abnormal case in union_exons" << endl;
		}
	}

	while (it_A != exonList_A.end())
	{
		exonList_union.push_back(make_pair(it_A->first, it_A->second));
		inList.push_back(1);
		++it_A;
	}
	while (it_B != exonList_B.end())
	{
		exonList_union.push_back(make_pair(it_B->first, it_B->second));
		inList.push_back(2);
		++it_B;
	}
}

// compare two paths and return the ASM category code, only consider the situation where two paths have the same range (i.e. do not consider alternative transcription start/end)
int compare_two_paths(alternative_path *path_A, alternative_path *path_B)
{
	vector<pair <long, long> > exonList_A, exonList_B, exonList_union;
	extract_exons(path_A, exonList_A);
	extract_exons(path_B, exonList_B);

	long path_low = max(path_A->path_start, path_B->path_start);
	long path_high = min(path_A->path_end, path_B->path_end);
	long bound_low, bound_high; // boundary of an exon from other exons or path ends on both sides
	bool shared_low, shared_high; // boundary is from an shared exon, either an shared exon in these 2 paths or the path ends
	
	vector<int> inList;
	union_exons(exonList_A, exonList_B, exonList_union, inList);

	if (exonList_union.empty()) {
		return unknown;
	}

	// modify the union exons by merging adjacent exons (no junction in between) from the same list
	for (int i = 0; i+1 < exonList_union.size();) {
		if (inList[i] == inList[i+1] && exonList_union[i].second + 1 >= exonList_union[i+1].first) {
			exonList_union[i].second = exonList_union[i+1].second;
			exonList_union.erase(exonList_union.begin() + i+1);
			inList.erase(inList.begin() + i+1);
		}
		else {
			 ++i;
		}
	}

	int num_exon = exonList_union.size();
	bool has_shared_exon = false;
	for (int i = 0; i < num_exon; ++i){
		if (inList[i] == 3) {
			has_shared_exon = true;
		}
		else {
			bound_low = i > 0 ? exonList_union[i-1].second : path_low;
			bound_high = (i+1) < num_exon ? exonList_union[i+1].first : path_high;
			if (i == 0 || (i > 0 && inList[i-1] == 3)) {
				shared_low = true;
			}
			else {
				shared_low = false;
			}
			if (i == num_exon-1 || (i < num_exon-1 && inList[i+1] == 3)) {
				shared_high = true;
			}
			else {
				shared_high = false;
			}

			if (exonList_union[i].first > bound_low + 1 && exonList_union[i].second < bound_high - 1) {
				return exon_skipping;
			}
			else if (shared_low && shared_high) {
				if (exonList_union[i].first <= bound_low + 1 && exonList_union[i].second >= bound_high - 1) {
					return intron_retention;
				}
				else if ((exonList_union[i].first > bound_low + 1 && exonList_union[i].second >= bound_high - 1) 
					|| (exonList_union[i].first <= bound_low + 1 && exonList_union[i].second < bound_high - 1)) {
						return alter_splice_site;
				}
			}			
		}
	}

	if (has_shared_exon == false && path_A->exonNum >= 1 && path_A->junctionNum >= 2 && path_B->exonNum >= 1 && path_B->junctionNum >= 2) {
		return mutual_exclusive;
	}

	return unknown;
}

bool dominating_pair(const pair<alternative_path*, alternative_path*> &pair_A, const pair<alternative_path*, alternative_path*> &pair_B){
	if (!pair_A.first || !pair_A.second || !pair_B.first || !pair_B.second) {
		cout << "exception: null pointer in dominating pair comparison" << endl;
	}
	return ((pair_A.first->avg_support + pair_A.second->avg_support) >= (pair_B.first->avg_support + pair_B.second->avg_support));
}
// try to determine the category of an alternative splicing module using dominant paths
int asm_category_dominant_path(alternative_path *major_paths, int major_path_num)
{
	if (!major_paths || major_path_num < 2)
	{
		return -1;
	}

	vector<pair<alternative_path*, alternative_path*> > path_pairs;
	alternative_path *path_A, *path_B;
	for (path_A = major_paths; path_A; path_A = path_A->next){
		for (path_B = path_A->next; path_B; path_B = path_B->next){
			path_pairs.push_back(make_pair(path_A, path_B));
		}
	}

	//sort the path pairs from the most dominant pair to the least important
	sort(path_pairs.begin(), path_pairs.end(), dominating_pair);

	int ASMcategory = 0; // set unknown by default, since we already have at least 2 major paths

	for (vector<pair<alternative_path*, alternative_path*> >::iterator it = path_pairs.begin(); ASMcategory <= 0 && it != path_pairs.end(); ++it){
		path_A = it->first;
		path_B = it->second;

// 		if (overlapping_path(path_A, path_B)){
// 			continue; // only consider disjoint pairs
// 		}

		if (abs(path_A->path_start - path_B->path_start) < 2 && abs(path_A->path_end - path_B->path_end) < 2)
		{
			ASMcategory = compare_two_paths(path_A, path_B);

// 			if ((path_A->exonNum == 0 && path_A->junctionNum == 1 && path_B->junctionNum >= 2 && path_B->exonNum >= 1) || (path_A->junctionNum >= 2 && path_A->exonNum >= 1 && path_B->exonNum == 0 && path_B->junctionNum == 1))
// 			{
// 				ASMcategory = exon_skipping;
// 			}
// 			else if (path_A->exonNum >= 1 && path_A->junctionNum >= 2 && path_B->exonNum >= 1 && path_B->junctionNum >= 2)
// 			{
// 				ASMcategory = mutual_exclusive;
// 			}
// 			else if ((path_A->exonNum == 1 && path_A->junctionNum == 0 && path_B->exonNum == 0 && path_B->junctionNum == 1) || (path_A->exonNum == 0 && path_A->junctionNum == 1 && path_B->exonNum == 1 && path_B->junctionNum == 0))
// 			{
// 				ASMcategory = intron_retention; // need to enforce num of paths to be =2?
// 			}
// 			else if ((path_A->exonNum == 0 && path_A->junctionNum == 1 && path_B->junctionNum == 1 && path_B->exonNum == 1) || (path_A->junctionNum == 1 && path_A->exonNum == 1 && path_B->exonNum == 0 && path_B->junctionNum == 1))
// 			{
// 				ASMcategory = alter_splice_site;
// 			}
		}
		else if (abs(path_A->path_start - path_B->path_start) < 2)
		{
			ASMcategory = alter_end;
		}
		else if (abs(path_A->path_end - path_B->path_end) < 2)
		{
			ASMcategory = alter_start;
		}
	}

	return ASMcategory;
}

//determine the category of an alternative splicing module
int asm_category(GTvertex *curVertex)
{
	//cout << curVertex->vertex_id << endl;
	//categorize alternative splicing events, structure only, no expression
// 	if (curVertex->rangeLow == 75348718 && curVertex->rangeHigh == 75355797){
// 		cout << "a";
// 	}

	//if the ASM has been assigned a category, then use that
	if (curVertex->ASMcategory > 0)
	{
		return curVertex->ASMcategory;
	}

	if ((curVertex->childType == 2 || curVertex->childType == 3) && curVertex->childNum > 1)// && curVertex->major_alter_paths_num > 1)
	{
		//ASM
	}
	else
	{
		return -1;
	}

	//categorize the alternative splicing
	if (curVertex->major_alter_paths && curVertex->major_alter_paths_num >= 2)
	{
		curVertex->ASMcategory = asm_category_dominant_path(curVertex->major_alter_paths, curVertex->major_alter_paths_num);
		return curVertex->ASMcategory;
	}
	else {
		curVertex->ASMcategory = unknown;
		return unknown;
	}

//	alternative_path *pathA, *pathB;
//
//	if (curVertex->childType == 2 || curVertex->childType == 3)
//	{
//		pathA = curVertex->major_alter_paths;
//		while (pathA != NULL)
//		{
//			pathB = pathA->next;
//			while (pathB != NULL)
//			{
//				if (abs(pathA->path_start - pathB->path_start) < 2 && abs(pathA->path_end - pathB->path_end) < 2)
//				{
//					if (pathA->exonNum == 0 && pathA->junctionNum == 1 && pathB->junctionNum >= 2 && pathB->exonNum >= 1 || pathA->junctionNum >= 2 && pathA->exonNum >= 1 && pathB->exonNum == 0 && pathB->junctionNum == 1)
//					{
//						curVertex->ASMcategory = exon_skipping;
//						return exon_skipping;
//					}
//					if (curVertex->ASMcategory <= 0 && curVertex->major_alter_paths_num == 2 && pathA->exonNum >= 1 && pathA->junctionNum >= 2 && pathB->exonNum >= 1 && pathB->junctionNum >= 2)
//					{
//						curVertex->ASMcategory = mutual_exclusive;
//					}
//					if (curVertex->major_alter_paths_num == 2 && (pathA->exonNum == 1 && pathA->junctionNum == 0 && pathB->exonNum == 0 && pathB->junctionNum == 1 || pathA->exonNum == 0 && pathA->junctionNum == 1 && pathB->exonNum == 1 && pathB->junctionNum == 0))
//					{
//						curVertex->ASMcategory = intron_retention;
//						return intron_retention;
//					}
//					if (curVertex->ASMcategory <= 0 && (pathA->exonNum == 0 && pathA->junctionNum == 1 && pathB->junctionNum == 1 && pathB->exonNum == 1 || pathA->junctionNum == 1 && pathA->exonNum == 1 && pathB->exonNum == 0 && pathB->junctionNum == 1))
//					{
//						curVertex->ASMcategory = alter_splice_site;
//					}
//				}
//				else if (abs(pathA->path_start - pathB->path_start) < 2)
//				{
//					curVertex->ASMcategory = alter_end;
//				}
//				else if (abs(pathA->path_end - pathB->path_end) < 2)
//				{
//					curVertex->ASMcategory = alter_start;
//				}
//
//				pathB = pathB->next;
//			}
//
//			pathA = pathA->next;
//		}
//	}
//
//	if (curVertex->ASMcategory > 0)
//	{
//		return curVertex->ASMcategory;
//	}
//
//	curVertex->ASMcategory = unknown;
//	return unknown;
}


void output_alterJunc(GTvertex *alterSite, ofstream *outputfile)
{
	//output alternative splice sites
	GTedge *childedge;
	rangeJunction *curJunc;
	int tmp;

	childedge = alterSite->child;
	while (childedge != NULL)
	{
		(*outputfile) << endl;
		for (tmp = 1; tmp <= alterSite->level + 7; tmp++)
		{
			(*outputfile) << "  ";
		}
		(*outputfile) << "| ";

		for (tmp = 0; tmp < SUPPORT_VECTOR_SIZE; tmp++)
		{
			(*outputfile) << childedge->linkedVertex->proportion[tmp] << ";";
		}
		(*outputfile) << "=  ";		

		curJunc = childedge->linkedVertex->junctionInRange->list;
		while (curJunc != NULL)
		{
			(*outputfile) << "(";
			if (curJunc->junc->type == frag_exon || curJunc->junc->type == frag_retained_intron)
			{
				(*outputfile) << "exon: ";
			} 
			else
			{
				(*outputfile) << "junc: ";
			}
			(*outputfile) << curJunc->junc->start << ", " << curJunc->junc->end << "; ";
			for (tmp = 0; tmp < SUPPORT_VECTOR_SIZE; tmp++)
				(*outputfile) << (curJunc->junc->support)[tmp] << "/ ";
			(*outputfile) << ") ";

			curJunc = curJunc->next;
		}

		childedge = childedge->next;
	}

	return;
}

bool GTree_output(GTvertex *rootVertex, int sign, ofstream *outputfile)
{
	//output GTree in pre-order
	//return true if still in a minimum ASM, i.e., still no type 2 or 3 nodes
	GTvertex *curVertex;
	GTedge *curEdge;
	rangeJunction *curJunc;
	long i;
	int vertexCategory, tmp;
	bool isMinASM, tmpFlag;


	for (i = 1; i <= rootVertex->level; i++)
	{
		(*outputfile) << "  ";
	}

	if (rootVertex->level > 0)
	{
		if (sign == 1)
		{
			//children are independent regions
			(*outputfile) << "-";
		} 
		else if (sign == 2)
		{
			//children are independent paths
			(*outputfile) << "+";
		}
		else if (sign == 3)
		{
			//children are dependent paths
			(*outputfile) << "x";
		}
	}

	(*outputfile) << "[" << rootVertex->rangeLow << ", " << rootVertex->rangeHigh << "] " << rootVertex->vertex_id << " " << rootVertex->childNum << '/' << rootVertex->junctionNum << " ";


	vertexCategory = asm_category(rootVertex);

	if (vertexCategory == exon_skipping)
	{
		exon_skipping_cnt++;
		(*outputfile) << "exon_skipping ";
	}
	else if (vertexCategory == mutual_exclusive)
	{
		mutual_exclusive_cnt++;
		(*outputfile) << "mutual_exclusive ";
	}
	else if (vertexCategory == intron_retention)
	{
		intron_retention_cnt++;
		(*outputfile) << "retained_intron ";
	}
	else if (vertexCategory == alter_splice_site)
	{
		(*outputfile) << "alter_splice_site ";
	}
	else if (vertexCategory == alter_start || vertexCategory == alter_end)
	{
		(*outputfile) << "alter_start/end ";
	}


	if (rootVertex->child == NULL)
	{
		if (OUTPUT_PER_SAMPLE_EXPR) {
			(*outputfile) << "{";
			for (tmp = 0; tmp < SUPPORT_VECTOR_SIZE; tmp++)
			{
				(*outputfile) << rootVertex->support[tmp]  << "/" << rootVertex->proportion[tmp] << "; ";
			}
			(*outputfile) << "}";
		}
		else {
			(*outputfile) << "{" << rootVertex->support[SUPPORT_VECTOR_SIZE]  << "/" << rootVertex->proportion[SUPPORT_VECTOR_SIZE] << "}";
		}

		//print fragment list
		(*outputfile) << ": ";
		curJunc = rootVertex->junctionInRange->list;
		while (curJunc != NULL)
		{
			(*outputfile) << "(";
			if (curJunc->junc->type == frag_exon || curJunc->junc->type == frag_retained_intron)
			{
				(*outputfile) << "exon: ";
			} 
			else
			{
				(*outputfile) << "junc: ";
			}
			(*outputfile) << curJunc->junc->start << "-" << curJunc->junc->end;

			if (setting_junc_anno_provided && curJunc->junc->type == frag_junction && !curJunc->junc->in_junc_annotation)
			{
				(*outputfile) << "n";
			}
			(*outputfile) << "; ";
			
			if (OUTPUT_PER_SAMPLE_EXPR) {
				for (tmp = 0; tmp < SUPPORT_VECTOR_SIZE; tmp++)
					(*outputfile) << (curJunc->junc->support)[tmp] << "/ ";
			}
			else {
				(*outputfile) << (curJunc->junc->support)[SUPPORT_VECTOR_SIZE] << "";
			}
			(*outputfile) << ") ";

			curJunc = curJunc->next;
		}



		// 		if (rootVertex->anovaScore_support > 3.46)
		// 		{
		// 			(*outputfile) << "DIFF";
		// 		}
		(*outputfile) << endl;

		if (rootVertex->alterSpliceSite == NULL)
		{

		}
		else
		{
			//output alternative splice sites

			curEdge = rootVertex->alterSpliceSite;
			while (curEdge != NULL)
			{
				curVertex = curEdge->linkedVertex;
				for (i = 1; i <= curVertex->level; i++)
				{
					(*outputfile) << "  ";
				}
				(*outputfile) << "  alter: " << curVertex->childNum << ". ";
				output_alterJunc(curVertex, outputfile);
				(*outputfile) << endl;

				if (curVertex->ASMcategory > 0 && rootVertex->ASMcategory != diff_expression)
				{
					asm_category(curVertex);
				}

				curEdge = curEdge->next;
			}
		}


		// 		if (rootVertex->junctionNum > 1 && sign != 3)
		// 		{
		// 			curJunc = rootVertex->junctionInRange->list;
		// 			while (curJunc != NULL)
		// 			{
		// 				if (curJunc->junc->type == frag_exon)
		// 				{
		// 					abnormalfile << "1\t";
		// 				} 
		// 				else
		// 				{
		// 					abnormalfile << "0\t";
		// 				}
		// 				abnormalfile << curJunc->junc->start << "\t" << curJunc->junc->end << endl;
		// 				curJunc = curJunc->next;
		// 			}
		// 			cout << "error: rootVertex->junctionNum > 1 && sign != 3" << endl;
		// 			exit(1);
		// 			//abnormalfile << endl << endl << endl;
		// 		}

		return true;
	} 
	else
	{
		//extend children
		// 		(*outputfile) << "[";
		// 		for (tmp = 0; tmp < SUPPORT_VECTOR_SIZE; tmp++)
		// 		{
		// 			(*outputfile) << rootVertex->representative->support[tmp] << "/ ";
		// 		}
		// 		(*outputfile) << "]  ";

		if (OUTPUT_PER_SAMPLE_EXPR) {
			(*outputfile) << "{";
			for (tmp = 0; tmp < SUPPORT_VECTOR_SIZE; tmp++)
			{
				(*outputfile) << rootVertex->support[tmp] << "/" << rootVertex->proportion[tmp] << "; ";
			}
			(*outputfile) << "}";
		}
		else {
			(*outputfile) << "{"  << rootVertex->support[SUPPORT_VECTOR_SIZE] << "/" << rootVertex->proportion[SUPPORT_VECTOR_SIZE] << "}";
		}

		if (rootVertex->ASMcategory != diff_expression)
		{
			(*outputfile) << endl;
		}


		isMinASM = true;

		curEdge = rootVertex->child;
		while (curEdge != NULL)
		{
			curVertex = curEdge->linkedVertex;
			tmpFlag = GTree_output(curVertex, rootVertex->childType, outputfile);

			if (tmpFlag == false)
			{
				isMinASM = false;
			}

			curEdge = curEdge->next;
		}


		//		if (isMinASM == true)
		//		{
		if (rootVertex->childType == 2 || rootVertex->childType == 3)
		{
			return false;
		} 
		else
		{
			return true;
		}
		//		} 
		//		else
		//		{
		//			return false;
		//		}

		//output alternative splice sites
		if (rootVertex->alterSpliceSite != NULL)
		{
			curEdge = rootVertex->alterSpliceSite;
			while (curEdge != NULL)
			{
				curVertex = curEdge->linkedVertex;
				for (i = 1; i <= curVertex->level; i++)
				{
					(*outputfile) << "  ";
				}
				(*outputfile) << "  alter: " << curVertex->childNum << ". ";
				output_alterJunc(curVertex, outputfile);
				(*outputfile) << endl;

				if (curVertex->ASMcategory > 0 && rootVertex->ASMcategory != diff_expression)
				{
					asm_category(curVertex);
				}

				curEdge = curEdge->next;
			}
		}
	}
}


void output_ASMpath_gtf(GTvertex *targetVertex, RangeJunctionList *pathList)
{
	long path_cnt = 0;

	RangeJunctionList *curPath;
	rangeJunction *curJunc;
	char strand;

	curPath = pathList;
	while (curPath != NULL)
	{
		++path_cnt;
		curPath = curPath->nextList;
	}

	curPath = pathList;
	while (curPath != NULL)
	{
		if (curPath->transDirection == antisense)
			strand = '-';
		else
			strand = '+';

		if (targetVertex->prevSibling != NULL && abs(curPath->rangeLow - targetVertex->prevSibling->rangeHigh) < 2)
			outfile_gtf_asm_path << chromosomeName << "\t" << "ASM" << "\t" << "exon" << "\t" << targetVertex->prevSibling->rangeLow << "\t" << targetVertex->prevSibling->rangeHigh << "\t" << "." << "\t" << strand << "\t" << "." << "\t"
			<< "gene_id \"" << chromosomeName << "." << targetVertex->vertex_id << "\"; " << "transcript_id \"" << chromosomeName << "." << targetVertex->vertex_id << ".p" << path_cnt << "\";" << endl;

		curJunc = curPath->list;
		while (curJunc != NULL)
		{
			if (curJunc->junc->type == frag_exon || curJunc->junc->type == frag_retained_intron)
				outfile_gtf_asm_path << chromosomeName << "\t" << "ASM" << "\t" << "exon" << "\t" << curJunc->junc->start << "\t" << curJunc->junc->end << "\t" << "." << "\t" << strand << "\t" << "." << "\t"
				<< "gene_id \"" << chromosomeName << "." << targetVertex->vertex_id << "\"; " << "transcript_id \"" << chromosomeName << "." << targetVertex->vertex_id << ".p" << path_cnt << "\";" << endl;

			curJunc = curJunc->next;
		}

		if (targetVertex->nextSibling != NULL && abs(curPath->rangeHigh - targetVertex->nextSibling->rangeLow) < 2)
			outfile_gtf_asm_path << chromosomeName << "\t" << "ASM" << "\t" << "exon" << "\t" << targetVertex->nextSibling->rangeLow << "\t" << targetVertex->nextSibling->rangeHigh << "\t" << "." << "\t" << strand << "\t" << "." << "\t"
			<< "gene_id \"" << chromosomeName << "." << targetVertex->vertex_id << "\"; " << "transcript_id \"" << chromosomeName << "." << targetVertex->vertex_id << ".p" << path_cnt << "\";" << endl;


		--path_cnt;
		curPath = curPath->nextList;
	}

	return;
}

void output_asm_path(GTvertex *targetVertex)
{
// 	if (targetVertex->vertex_id.compare("chr14.asm71") == 0)
// 	{
// 		cout << "a";
// 	}
	
	alternative_path *asm_path;
	GTedge *path_component;
	GTvertex *cur_component;
	rangeJunction *curJunc;
	long path_cnt = 0, exon_cnt;
	char strand;
	string name_suffix1, name_suffix2;

	if (setting_junc_anno_provided && targetVertex->has_novel_junction_major)
		name_suffix1 = "n";
	else
		name_suffix1 = "";

	asm_path = targetVertex->major_alter_paths;
	while (asm_path != NULL)
	{
		++path_cnt;
		exon_cnt = 0;

		if (asm_path->transDirection == antisense)
			strand = '-';
		else
			strand = '+';

		if (setting_junc_anno_provided && asm_path->pathVertex->has_novel_junction)
			name_suffix2 = "n";
		else
			name_suffix2 = "";

		if (targetVertex->prevSibling != NULL && targetVertex->prevSibling->exonNum > 0 && abs(asm_path->path_start - targetVertex->prevSibling->rangeHigh) < 2)
			outfile_gtf_asm_path_decomposed << chromosomeName << "\t" << "ASM" << "\t" << "exon" << "\t" << targetVertex->prevSibling->rangeLow << "\t" << targetVertex->prevSibling->rangeHigh << "\t" << "." << "\t" << strand << "\t" << "." << "\t"
				<< "gene_id \"" << chromosomeName << "." << targetVertex->vertex_id + name_suffix1 << "\"; " << "transcript_id \"" << chromosomeName << "." << asm_path->pathVertex->vertex_id + name_suffix2 << "\";" << endl;

		if (asm_path->pathVertex->childType == 0 || asm_path->pathVertex->childType == 3)
		{
			curJunc = asm_path->pathVertex->junctionInRange->list;
			while (curJunc != NULL)
			{
				if (curJunc->junc->type == frag_exon || curJunc->junc->type == frag_retained_intron)
					outfile_gtf_asm_path_decomposed << chromosomeName << "\t" << "ASM" << "\t" << "exon" << "\t" << curJunc->junc->start << "\t" << curJunc->junc->end << "\t" << "." << "\t" << strand << "\t" << "." << "\t"
					<< "gene_id \"" << chromosomeName << "." << targetVertex->vertex_id + name_suffix1 << "\"; " << "transcript_id \"" << chromosomeName << "." << asm_path->pathVertex->vertex_id + name_suffix2 << "\";" << endl;

				curJunc = curJunc->next;
			}
		}
		else if (asm_path->pathVertex->childType == 1 || asm_path->pathVertex->childType == 2)
		{
			path_component = asm_path->pathVertex->child;
			while (path_component != NULL)
			{
				cur_component = path_component->linkedVertex;
				if (cur_component->exonNum > 0)
					outfile_gtf_asm_path_decomposed << chromosomeName << "\t" << "ASM" << "\t" << "exon" << "\t" << cur_component->rangeLow << "\t" << cur_component->rangeHigh << "\t" << "." << "\t" << strand << "\t" << "." << "\t"
					<< "gene_id \"" << chromosomeName << "." << targetVertex->vertex_id + name_suffix1 << "\"; " << "transcript_id \"" << chromosomeName << "." << asm_path->pathVertex->vertex_id + name_suffix2 << "\";" << endl;

				path_component = path_component->next;
			}
		}

		if (targetVertex->nextSibling != NULL && targetVertex->nextSibling->exonNum > 0 && abs(asm_path->path_end - targetVertex->nextSibling->rangeLow) < 2)
			outfile_gtf_asm_path_decomposed << chromosomeName << "\t" << "ASM" << "\t" << "exon" << "\t" << targetVertex->nextSibling->rangeLow << "\t" << targetVertex->nextSibling->rangeHigh << "\t" << "." << "\t" << strand << "\t" << "." << "\t"
				<< "gene_id \"" << chromosomeName << "." << targetVertex->vertex_id + name_suffix1 << "\"; " << "transcript_id \"" << chromosomeName << "." << asm_path->pathVertex->vertex_id + name_suffix2 << "\";" << endl;

		asm_path = asm_path->next;
	}

	return;
}



/************************************************************************/
/* DIFFERENTIAL EXPRESSION ANALYSIS                                     */
/************************************************************************/

void get_gene_expression_analysis_list(GTvertex *rootVertex, bool is_GTree_root)
{
	GTedge *curEdge;
	//int groupLoopCnt, individualLoopCnt, techRepLoopCnt, sampleSizeCnt, supportIndex;

	if (rootVertex->level < 1 && !is_GTree_root && rootVertex->exonNum >= 1)
		vertexListForStatistics.push_back(rootVertex);

	if (rootVertex->child == NULL)
	{
		//do nothing
	}
	else
	{
		//extend children
		curEdge = rootVertex->child;
		while (curEdge != NULL)
		{
			get_gene_expression_analysis_list(curEdge->linkedVertex, false);
			curEdge = curEdge->next;
		}
	}

	return;
}

void differential_analysis_expression_level(GTvertex *rootVertex)
{
	unsigned long geneLoopCnt;
	GTvertex *curVertex;
	int sampleLoopCnt;

	//get expression matrix
	vertexListForStatistics.clear();
	get_gene_expression_analysis_list(rootVertex, true);

	for (geneLoopCnt = 0; geneLoopCnt < vertexListForStatistics.size(); ++geneLoopCnt)
	{
		curVertex = vertexListForStatistics[geneLoopCnt];
		
		outfile_stats_expr << chromosomeName << "\t" << curVertex->rangeLow << "\t" << curVertex->rangeHigh;
		for (sampleLoopCnt = 0; sampleLoopCnt < SUPPORT_VECTOR_SIZE; ++sampleLoopCnt)
		{
			outfile_stats_expr << "\t" << curVertex->support[sampleLoopCnt];
		}
		for (sampleLoopCnt = 0; sampleLoopCnt < SUPPORT_VECTOR_SIZE; ++sampleLoopCnt)
		{
			outfile_stats_expr << "\t" << curVertex->proportion[sampleLoopCnt];
		}
		outfile_stats_expr << endl;

		outfile_stats_expr_new << chromosomeName << "." << curVertex->vertex_id;
		for (sampleLoopCnt = 0; sampleLoopCnt < SUPPORT_VECTOR_SIZE; ++sampleLoopCnt)
		{
			outfile_stats_expr_new << "\t" << curVertex->support[sampleLoopCnt];
		}
		outfile_stats_expr_new << endl;
	}

	vertexListForStatistics.clear();
	vector <GTvertex*> tmp;
	tmp.swap(vertexListForStatistics);

	return;
}


/************************************************************************/
/* DIFFERENTIAL TRANSCRIPTION ANALYSIS                                  */
/************************************************************************/


void output_one_asm(GTvertex *targetVertex)
{
	int sampleLoopCnt;
	alternative_path *curAlterPath;
	string ASM_id, name_suffix;
	
	if (setting_junc_anno_provided && targetVertex->has_novel_junction_major)
		name_suffix = "n";
	else
		name_suffix = "";
	
	ASM_id = chromosomeName + "." + targetVertex->vertex_id + name_suffix;

	//4/22/2013, only output high expression ASMs for ANOVA analysis, will move this step to later modules (diff analysis)
	bool output_for_stat_asm_new = false;
	if (calculate_mean(targetVertex->support, SUPPORT_VECTOR_SIZE) >= 15)
	{
		output_for_stat_asm_new = true;

		outfile_stat_asmexpr_new << ASM_id;
		for (sampleLoopCnt = 0; sampleLoopCnt < SUPPORT_VECTOR_SIZE; ++sampleLoopCnt)
		{
			outfile_stat_asmexpr_new << "\t" << targetVertex->support[sampleLoopCnt];
		}
		outfile_stat_asmexpr_new << endl;
	}

	//output vertex information
	outfile_stat_asm << ASM_id << "\t" << chromosomeName << "\t" << targetVertex->rangeLow << "\t" << targetVertex->rangeHigh << "\t" << targetVertex->major_alter_paths_num << "\t" << targetVertex->ASMcategory;
	for (sampleLoopCnt = 0; sampleLoopCnt < SUPPORT_VECTOR_SIZE; ++sampleLoopCnt)
	{
		outfile_stat_asm << "\t" << targetVertex->support[sampleLoopCnt] << "\t" << targetVertex->MSE_estimation[sampleLoopCnt];
	}
	outfile_stat_asm << endl;

	//output alternative paths
	curAlterPath = targetVertex->major_alter_paths;
	while (curAlterPath != NULL)
	{
		if (setting_junc_anno_provided && curAlterPath->pathVertex->has_novel_junction)
			name_suffix = "n";
		else
			name_suffix = "";

		outfile_stat_asm << chromosomeName << "." << curAlterPath->pathVertex->vertex_id + name_suffix << "\t";
		outfile_stat_asm << curAlterPath->whole_path_start << "\t" << curAlterPath->whole_path_end;
		for (sampleLoopCnt = 0; sampleLoopCnt < SUPPORT_VECTOR_SIZE; ++sampleLoopCnt)
		{
			outfile_stat_asm << "\t" << curAlterPath->support[sampleLoopCnt] << "\t" << curAlterPath->proportion[sampleLoopCnt];
		}
		outfile_stat_asm << endl;


		if (output_for_stat_asm_new)
		{
			outfile_stat_asm_new_absexpr << chromosomeName << "." << curAlterPath->pathVertex->vertex_id + name_suffix;
			outfile_stat_asm_new_proportion << chromosomeName << "." << curAlterPath->pathVertex->vertex_id + name_suffix;
			for (sampleLoopCnt = 0; sampleLoopCnt < SUPPORT_VECTOR_SIZE; ++sampleLoopCnt)
			{
				outfile_stat_asm_new_absexpr << "\t" << curAlterPath->support[sampleLoopCnt];
				outfile_stat_asm_new_proportion << "\t" << curAlterPath->proportion[sampleLoopCnt];
			}
			outfile_stat_asm_new_absexpr << endl;
			outfile_stat_asm_new_proportion << endl;
		}
		

		curAlterPath = curAlterPath->next;
	}

	output_asm_path(targetVertex);

	return;
}

void output_ASM_composition(GTvertex *targetVertex)
{
	//output the composition of an ASM, including its children and extending all its type-1 children
	double totalExpr = 0.0;
	long composCnt = 0;
	rangeJunction *curJunc;
	GTvertex *curVertex, *outputVertex;
	GTedge *curEdge, *outputEdge;

	//output basic information
	outfile_asm_composition << targetVertex->vertex_id << "\t" << chromosomeName << "\t" << targetVertex->rangeLow << "\t" << targetVertex->rangeHigh << "\t" << targetVertex->ASMcategory << "\t" << targetVertex->childNum << "\t";

	for (int sampleLoopCnt = 0; sampleLoopCnt < SUPPORT_VECTOR_SIZE; ++sampleLoopCnt)
		totalExpr += targetVertex->support[sampleLoopCnt];
	outfile_asm_composition << totalExpr / SUPPORT_VECTOR_SIZE << "\t";

	//count number of compositions & output compositions
	if (targetVertex->childType == 3)
	{
		curJunc = targetVertex->junctionInRange->list;
		while (curJunc != NULL)
		{
			if (curJunc->junc->type == frag_exon || curJunc->junc->type == frag_retained_intron || curJunc->junc->type == frag_junction)
				++composCnt;
			curJunc = curJunc->next;
		}

		outfile_asm_composition << composCnt << "\t";

		curJunc = targetVertex->junctionInRange->list;
		while (curJunc != NULL)
		{
			if (curJunc->junc->type == frag_exon || curJunc->junc->type == frag_retained_intron || curJunc->junc->type == frag_junction)
			{
				if (curJunc->junc->type == frag_junction)
					outfile_asm_composition << "1\t";
				else
					outfile_asm_composition << "0\t";
				if (curJunc->junc->start_real > 0)
					outfile_asm_composition << curJunc->junc->start_real << "\t";
				else 
					outfile_asm_composition << curJunc->junc->start << "\t";
				if (curJunc->junc->end_real > 0)
					outfile_asm_composition << curJunc->junc->end_real << "\t";
				else 
					outfile_asm_composition << curJunc->junc->end << "\t";
			}
			curJunc = curJunc->next;
		}
	}
	else if (targetVertex->childType == 2)
	{
		curEdge = targetVertex->child;
		while (curEdge != NULL)
		{
			curVertex = curEdge->linkedVertex;
			if (curVertex->estimated == false)
			{
				curJunc = curVertex->junctionInRange->list;
				while (curJunc != NULL)
				{
					if (curJunc->junc->type == frag_exon || curJunc->junc->type == frag_retained_intron || curJunc->junc->type == frag_junction)
						++composCnt;
					curJunc = curJunc->next;
				}
			} 
			else
			{
				if (curVertex->childType == 1)
				{
					composCnt += curVertex->childNum;
				}
				else
					++composCnt;
			}

			curEdge = curEdge->next;
		}

		outfile_asm_composition << composCnt << "\t";

		curEdge = targetVertex->child;
		while (curEdge != NULL)
		{
			curVertex = curEdge->linkedVertex;
			if (curVertex->estimated == false)
			{
				curJunc = curVertex->junctionInRange->list;
				while (curJunc != NULL)
				{
					if (curJunc->junc->type == frag_exon || curJunc->junc->type == frag_retained_intron || curJunc->junc->type == frag_junction)
					{
						if (curJunc->junc->type == frag_junction)
							outfile_asm_composition << "1\t";
						else
							outfile_asm_composition << "0\t";
						if (curJunc->junc->start_real > 0)
							outfile_asm_composition << curJunc->junc->start_real << "\t";
						else 
							outfile_asm_composition << curJunc->junc->start << "\t";
						if (curJunc->junc->end_real > 0)
							outfile_asm_composition << curJunc->junc->end_real << "\t";
						else 
							outfile_asm_composition << curJunc->junc->end << "\t";
					}
					curJunc = curJunc->next;
				}
			} 
			else
			{
				if (curVertex->childType == 1)
				{
					outputEdge = curVertex->child;
					while (outputEdge != NULL)
					{
						outputVertex = outputEdge->linkedVertex;
						if (outputVertex->childType == 0)
						{
							if (outputVertex->junctionInRange->list->junc->type == frag_junction)
								outfile_asm_composition << "1\t";
							else
								outfile_asm_composition << "0\t";
							if (outputVertex->junctionInRange->list->junc->start_real > 0)
								outfile_asm_composition << outputVertex->junctionInRange->list->junc->start_real << "\t";
							else
								outfile_asm_composition << outputVertex->junctionInRange->list->junc->start << "\t";
							if (outputVertex->junctionInRange->list->junc->end_real > 0)
								outfile_asm_composition << outputVertex->junctionInRange->list->junc->end_real << "\t";
							else
								outfile_asm_composition << outputVertex->junctionInRange->list->junc->end << "\t";
						}
						else
						{
							outfile_asm_composition << "2\t" << outputVertex->rangeLow << "\t" << outputVertex->rangeHigh << "\t";
						}

						outputEdge = outputEdge->next;
					}
				}
				else
					outfile_asm_composition << "2\t" << curVertex->rangeLow << "\t" << curVertex->rangeHigh << "\t";
			}

			curEdge = curEdge->next;
		}
	}

	outfile_asm_composition << endl;

	return;
}

// determine the strand of an exon in the gtf output
string determine_exon_strand(fragment *cur_exon)
{
	string strand;
	if (cur_exon->exon_strand_start == sense && cur_exon->exon_strand_end == sense)
		strand = "+";
	else if (cur_exon->exon_strand_start == antisense && cur_exon->exon_strand_end == antisense)
		strand = "-";
	else if (cur_exon->exon_strand_start == undetermined && cur_exon->exon_strand_end == antisense)
		strand = "-";
	else if (cur_exon->exon_strand_start == antisense && cur_exon->exon_strand_end == undetermined)
		strand = "-";
	else
		strand = "+";
	
	return strand;
}

void output_asm_analysis(GTvertex *rootVertex, bool is_gtree_root)
{
	//output all asms for differential transcription analysis
	GTvertex *curVertex;
	GTedge *curEdge;

	asm_category(rootVertex);

// 	if (rootVertex->level < 1)
// 	{
// 		++GENEcount;
// 	}

	if (rootVertex->level < 1 && !is_gtree_root && rootVertex->childType <= 1) // 10/13/2013 if the gene has a childtype>1, it should be enumerated for paths and these paths should be in the asm_path gtf, so no output for these genes
	{
		long segmentCnt = 1;

		if (rootVertex->child == NULL && rootVertex->exonNum > 0)
		{
// 			//rootVertex is a exon
// 			outfile_gtf_splice_graph << chromosomeName << "\t" << "ASM" << "\t" << "exon" << "\t" << rootVertex->rangeLow << "\t" << rootVertex->rangeHigh << "\t" << "." << "\t" << "+\t" << "." << "\t"
// 				<< "gene_id \"" << chromosomeName << "_" << rootVertex->vertex_id << "\"; " << "transcript_id \"" << chromosomeName << "_" << rootVertex->vertex_id << "_" << segmentCnt << "\";" << endl;

			// 10/13/2013 rootVertex is either an exon or a complex region with no path enumeration
			rangeJunction *curJunc = rootVertex->junctionInRange->list;
			while (curJunc != NULL)
			{
				string strand = "+";// determine_exon_strand(curJunc->junc);
				if (curJunc->junc->type == frag_exon || curJunc->junc->type == frag_retained_intron)
					outfile_gtf_splice_graph << chromosomeName << "\t" << "ASM" << "\t" << "exon" << "\t" << curJunc->junc->start << "\t" << curJunc->junc->end << "\t" << "." << "\t" << strand << "\t" << "." << "\t"
					<< "gene_id \"" << chromosomeName << "_" << rootVertex->vertex_id << "\"; " << "transcript_id \"" << chromosomeName << "_" << rootVertex->vertex_id << "_" << segmentCnt << "\";" << endl;

				curJunc = curJunc->next;
			}
		}
		else if (rootVertex->child != NULL)
		{
			curEdge = rootVertex->child;
			while (curEdge != NULL)
			{
				if (curEdge->linkedVertex->child == NULL)
				{
					// 10/13/2013 this child is either an exon or a complex region with no path enumeration
					rangeJunction *curJunc = curEdge->linkedVertex->junctionInRange->list;
					while (curJunc != NULL)
					{
						string strand = "+";// determine_exon_strand(curJunc->junc);
						if (curJunc->junc->type == frag_exon || curJunc->junc->type == frag_retained_intron)
							outfile_gtf_splice_graph << chromosomeName << "\t" << "ASM" << "\t" << "exon" << "\t" << curJunc->junc->start << "\t" << curJunc->junc->end << "\t" << "." << "\t" << strand << "\t" << "." << "\t"
							<< "gene_id \"" << chromosomeName << "_" << rootVertex->vertex_id << "\"; " << "transcript_id \"" << chromosomeName << "_" << rootVertex->vertex_id << "_" << segmentCnt << "\";" << endl;

						curJunc = curJunc->next;
					}
				} 
				else
				{
					++segmentCnt;
				}
				curEdge = curEdge->next;
			}
		}
	}

	if (rootVertex->child == NULL)
	{
		//alternative splice site
		if (rootVertex->alterSpliceSite != NULL)
		{
			curEdge = rootVertex->alterSpliceSite;
			while (curEdge != NULL)
			{
				curVertex = curEdge->linkedVertex;
				if (curVertex->ASMcategory > 0 && curVertex->major_alter_paths_num > 1)
				{
					output_one_asm(curVertex);
					output_ASM_composition(curVertex);
				}

				curEdge = curEdge->next;
			}
		}

		return;
	} 
	else
	{
		//extend children
		curEdge = rootVertex->child;
		while (curEdge != NULL)
		{
			output_asm_analysis(curEdge->linkedVertex, false);
			curEdge = curEdge->next;
		}

		if (rootVertex->childType == 2 || rootVertex->childType == 3)
		{
			if (rootVertex->major_alter_paths_num > 1)
				output_one_asm(rootVertex);
			output_ASM_composition(rootVertex);
		}


		if (rootVertex->alterSpliceSite != NULL)
		{
			curEdge = rootVertex->alterSpliceSite;
			while (curEdge != NULL)
			{
				curVertex = curEdge->linkedVertex;
				if (curVertex->ASMcategory > 0 && curVertex->major_alter_paths_num > 1)
				{
					output_one_asm(curVertex);
					output_ASM_composition(curVertex);
				}

				curEdge = curEdge->next;
			}
		}

		return;
	}
}

void differential_analysis_asm(GTvertex *rootVertex)
{
	output_asm_analysis(rootVertex, true);

	return;
}


/************************************************************************/
/* GENERAL PROCESSES				                                    */
/************************************************************************/

void initialize_outputfiles()
{	
	string filename;

#ifdef UNIX
	filename = resultPath + "detail/" + chromosomeName + "_GTree.txt";
	outfile_gtree.open(filename.c_str());

	filename = resultPath + "stat/expression.txt";
	outfile_stats_expr.open(filename.c_str(), fstream::app);
	outfile_stats_expr.precision(10);
	filename = resultPath + "stat/asm.txt";
	outfile_stat_asm.open(filename.c_str(), fstream::app);
	outfile_stat_asm.precision(10);

	filename = resultPath + "asm_path.gtf";
	outfile_gtf_asm_path.open(filename.c_str(), fstream::app);
	filename = resultPath + "splice_graph.gtf";
	outfile_gtf_splice_graph.open(filename.c_str(), fstream::app);
	filename = resultPath + "asm_path_decomposed.gtf";
	outfile_gtf_asm_path_decomposed.open(filename.c_str(), fstream::app);

	filename = resultPath + "stat/splice_all.bed";
	outfile_junction_all.open(filename.c_str(), fstream::app);
	filename = resultPath + "stat/splice_filtered.bed";
	outfile_junction_filtered.open(filename.c_str(), fstream::app);

	filename = resultPath + "detail/" + chromosomeName + "_asm.txt";
	outfile_asm_composition.open(filename.c_str());

	filename = resultPath + "stat/temp/expr_new.txt";
	outfile_stats_expr_new.open(filename.c_str(), fstream::app);
	outfile_stats_expr_new.precision(10);
	filename = resultPath + "stat/temp/asm_new_absexpr.txt";
	outfile_stat_asm_new_absexpr.open(filename.c_str(), fstream::app);
	outfile_stat_asm_new_absexpr.precision(10);
	filename = resultPath + "stat/temp/asm_new_proportion.txt";
	outfile_stat_asm_new_proportion.open(filename.c_str(), fstream::app);
	outfile_stat_asm_new_proportion.precision(10);
	filename = resultPath + "stat/temp/asm_expr_new.txt";
	outfile_stat_asmexpr_new.open(filename.c_str(), fstream::app);
	outfile_stat_asmexpr_new.precision(10);

	filename = resultPath + "debug_info.txt";
	outfile_debuginfo.open(filename.c_str(), fstream::app); outfile_debuginfo << chromosomeName << endl;
#else
	filename = resultPath + "GTree.txt";
	outfile_gtree.open(filename.c_str());

	filename = resultPath + "stat/expression.txt";
	outfile_stats_expr.open(filename.c_str());
	outfile_stats_expr.precision(10);
	filename = resultPath + "stat/asm.txt";
	outfile_stat_asm.open(filename.c_str());
	outfile_stat_asm.precision(10);


	filename = resultPath + "asm_path.gtf";
	outfile_gtf_asm_path.open(filename.c_str());
	outfile_gtf_asm_path << "browser full ASM\ntrack name=\"ASM\" description=\"ASM\" visibility=2 useScore=1\n";
	filename = resultPath + "non_asm.gtf";
	outfile_gtf_splice_graph.open(filename.c_str());
	outfile_gtf_splice_graph << "browser dense Gene\ntrack name=\"Gene structure\" description=\"Gene structure\" visibility=1 useScore=1\n";
	filename = resultPath + "asm_path_decomposed.gtf";
	outfile_gtf_asm_path_decomposed.open(filename.c_str());
	outfile_gtf_asm_path_decomposed << "browser full ASM_decomposed\ntrack name=\"ASM_decomposed\" description=\"ASM_decomposed\" visibility=2 useScore=1\n";

	filename = resultPath + "stat/expr_new.txt";
	outfile_stats_expr_new.open(filename.c_str());
	outfile_stats_expr_new.precision(10);
	filename = resultPath + "stat/asm_new_absexpr.txt";
	outfile_stat_asm_new_absexpr.open(filename.c_str());
	outfile_stat_asm_new_absexpr.precision(10);
	filename = resultPath + "stat/asm_new_proportion.txt";
	outfile_stat_asm_new_proportion.open(filename.c_str());
	outfile_stat_asm_new_proportion.precision(10);
	filename = resultPath + "stat/asm_expr_new.txt";
	outfile_stat_asmexpr_new.open(filename.c_str());
	outfile_stat_asmexpr_new.precision(10);

	filename = resultPath + "debug_info.txt";
	outfile_debuginfo.open(filename.c_str());
#endif

	return;
}

void close_outputfiles()
{
	outfile_gtree.close();
	outfile_stats_expr.close();
	outfile_stat_asm.close();
	outfile_gtf_asm_path.close();
	outfile_gtf_splice_graph.close();
	outfile_gtf_asm_path_decomposed.close();
	outfile_junction_all.close();
	outfile_junction_filtered.close();
	outfile_asm_composition.close();

	outfile_stats_expr_new.close();
	outfile_stat_asm_new_absexpr.close();
	outfile_stat_asm_new_proportion.close();
	outfile_stat_asmexpr_new.close();

	outfile_debuginfo << endl; outfile_debuginfo.close();

	return;
}

//deal with all sorts of outputs of the procedure
void output(GTvertex *rootVertex)
{

	/************************************************************************/
	/* DIFFERENTIAL ANALYSIS                                                */
	/************************************************************************/
	differential_analysis_expression_level(rootVertex); //output estimation results for differential gene expression analysis
	differential_analysis_asm(rootVertex); //output estimation results for differential transcription analysis


	/************************************************************************/
	/* OUTPUT                                                               */
	/************************************************************************/

	//output
	GTree_output(rootVertex, 1, &outfile_gtree);


	return;
}




bool ESGviewer_output(GTvertex *rootVertex, int sign, ofstream *outputfile)
{
	//output GTree in pre-order
	//return true if still in a minimum ASM, i.e., still no type 2 or 3 nodes
	GTvertex *curVertex;
	GTedge *curEdge;
	rangeJunction *curJunc;
	long i;
	int vertexCategory, tmp;
	bool isMinASM, tmpFlag;


	for (i = 1; i <= rootVertex->level; i++)
	{
		(*outputfile) << "  ";
	}

	if (rootVertex->level > 0)
	{
		if (sign == 1)
		{
			//children are independent regions
			(*outputfile) << "-";
		} 
		else if (sign == 2)
		{
			//children are independent paths
			(*outputfile) << "+";
		}
		else if (sign == 3)
		{
			//children are dependent paths
			(*outputfile) << "x";
		}
	}

	(*outputfile) << "[" << rootVertex->rangeLow << ", " << rootVertex->rangeHigh << "] " << rootVertex->vertex_id << " " << rootVertex->childNum << '/' << rootVertex->junctionNum << " ";


	vertexCategory = asm_category(rootVertex);

	if (vertexCategory == exon_skipping)
	{
		exon_skipping_cnt++;
		(*outputfile) << "exon_skipping ";
	}
	else if (vertexCategory == mutual_exclusive)
	{
		mutual_exclusive_cnt++;
		(*outputfile) << "mutual_exclusive ";
	}
	else if (vertexCategory == intron_retention)
	{
		intron_retention_cnt++;
		(*outputfile) << "retained_intron ";
	}
	else if (vertexCategory == diff_expression)
	{
		(*outputfile) << "diff_expression ";
	}


	if (rootVertex->child == NULL)
	{
		if (OUTPUT_PER_SAMPLE_EXPR) {
			(*outputfile) << "{";
			for (tmp = 0; tmp < SUPPORT_VECTOR_SIZE; tmp++)
			{
				(*outputfile) << rootVertex->support[tmp]  << "/" << rootVertex->proportion[tmp] << "; ";
			}
			(*outputfile) << "}";
		}

		//print fragment list
		(*outputfile) << ": ";
		curJunc = rootVertex->junctionInRange->list;
		while (curJunc != NULL)
		{
			(*outputfile) << "(";
			if (curJunc->junc->type == frag_exon || curJunc->junc->type == frag_retained_intron)
			{
				(*outputfile) << "exon: ";
			} 
			else
			{
				(*outputfile) << "junc: ";
			}
			(*outputfile) << curJunc->junc->start << ", " << curJunc->junc->end << "; ";

			if (setting_junc_anno_provided && curJunc->junc->type == frag_junction && !curJunc->junc->in_junc_annotation)
			{
				(*outputfile) << "n";
			}

			if (OUTPUT_PER_SAMPLE_EXPR) {
				for (tmp = 0; tmp < SUPPORT_VECTOR_SIZE; tmp++)
					(*outputfile) << (curJunc->junc->support)[tmp] << "/ ";
			}
			(*outputfile) << ") ";

			curJunc = curJunc->next;
		}



		// 		if (rootVertex->anovaScore_support > 3.46)
		// 		{
		// 			(*outputfile) << "DIFF";
		// 		}
		(*outputfile) << endl;

		if (rootVertex->alterSpliceSite == NULL)
		{

		}
		else
		{
			//output alternative splice sites

			curEdge = rootVertex->alterSpliceSite;
			while (curEdge != NULL)
			{
				curVertex = curEdge->linkedVertex;
				for (i = 1; i <= curVertex->level; i++)
				{
					(*outputfile) << "  ";
				}
				(*outputfile) << "  alter: " << curVertex->childNum << ". ";
				output_alterJunc(curVertex, outputfile);
				(*outputfile) << endl;

				if (curVertex->ASMcategory > 0 && rootVertex->ASMcategory != diff_expression)
				{
					asm_category(curVertex);
				}

				curEdge = curEdge->next;
			}
		}


		// 		if (rootVertex->junctionNum > 1 && sign != 3)
		// 		{
		// 			curJunc = rootVertex->junctionInRange->list;
		// 			while (curJunc != NULL)
		// 			{
		// 				if (curJunc->junc->type == frag_exon)
		// 				{
		// 					abnormalfile << "1\t";
		// 				} 
		// 				else
		// 				{
		// 					abnormalfile << "0\t";
		// 				}
		// 				abnormalfile << curJunc->junc->start << "\t" << curJunc->junc->end << endl;
		// 				curJunc = curJunc->next;
		// 			}
		// 			cout << "error: rootVertex->junctionNum > 1 && sign != 3" << endl;
		// 			exit(1);
		// 			//abnormalfile << endl << endl << endl;
		// 		}

		return true;
	} 
	else
	{
		//extend children
		// 		(*outputfile) << "[";
		// 		for (tmp = 0; tmp < SUPPORT_VECTOR_SIZE; tmp++)
		// 		{
		// 			(*outputfile) << rootVertex->representative->support[tmp] << "/ ";
		// 		}
		// 		(*outputfile) << "]  ";

		if (OUTPUT_PER_SAMPLE_EXPR) {
			(*outputfile) << "{";
			for (tmp = 0; tmp < SUPPORT_VECTOR_SIZE; tmp++)
			{
				(*outputfile) << rootVertex->support[tmp] << "/" << rootVertex->proportion[tmp] << "; ";
			}
			(*outputfile) << "}";
		}

		if (rootVertex->ASMcategory != diff_expression)
		{
			(*outputfile) << endl;
		}


		isMinASM = true;

		curEdge = rootVertex->child;
		while (curEdge != NULL)
		{
			curVertex = curEdge->linkedVertex;
			tmpFlag = GTree_output(curVertex, rootVertex->childType, outputfile);

			if (tmpFlag == false)
			{
				isMinASM = false;
			}

			curEdge = curEdge->next;
		}


		//		if (isMinASM == true)
		//		{
		if (rootVertex->childType == 2 || rootVertex->childType == 3)
		{
			return false;
		} 
		else
		{
			return true;
		}
		//		} 
		//		else
		//		{
		//			return false;
		//		}

		//output alternative splice sites
		if (rootVertex->alterSpliceSite != NULL)
		{
			curEdge = rootVertex->alterSpliceSite;
			while (curEdge != NULL)
			{
				curVertex = curEdge->linkedVertex;
				for (i = 1; i <= curVertex->level; i++)
				{
					(*outputfile) << "  ";
				}
				(*outputfile) << "  alter: " << curVertex->childNum << ". ";
				output_alterJunc(curVertex, outputfile);
				(*outputfile) << endl;

				if (curVertex->ASMcategory > 0 && rootVertex->ASMcategory != diff_expression)
				{
					asm_category(curVertex);
				}

				curEdge = curEdge->next;
			}
		}
	}
}
