
#include "splice_graph.h"




/************************************************************************/
/* GENOME TREE	                                                        */
/************************************************************************/


alter_junction::alter_junction()
{
	juncInfo = NULL;
	//	category = 0;
	//	proportion = 0.0;
	next = NULL;
}


fragment::fragment(int support_vec_size = SUPPORT_VECTOR_SIZE)
{
	frag_id = 0;
	transDirection = undetermined;
	
	vis_exonid = 0;
	splice_exonhead = NULL;
	splice_exontail = NULL;

	start = 0;
	end = 0;

	type = frag_junction;
	altersite = 0;
	in_exon_annotation = na_exon_ref;
	in_junc_annotation = false;
	bound_start = na_exon_bound;
	bound_end = na_exon_bound;

	exon_strand_start = undetermined;
	exon_strand_end = undetermined;

	support = new double [support_vec_size+1];
	for (int tmpCnt = 0; tmpCnt <= support_vec_size; ++tmpCnt)
		support[tmpCnt] = 0.0;

	alter = NULL;
	//	fivePalterTail = NULL;
	alterFragCnt = 0;
	coverage = NULL;

	start_real = 0;
	end_real = 0;
}

fragment* fragment::clone(int support_vec_size = SUPPORT_VECTOR_SIZE)
{
	fragment *newFrag = new fragment;

	newFrag->frag_name = frag_name;
	newFrag->frag_id = frag_id;
	newFrag->chromosome_start = chromosome_start;
	newFrag->chromosome_end = chromosome_end;
	newFrag->start = start;
	newFrag->end = end;
	newFrag->type = type;

	newFrag->transDirection = transDirection;
	newFrag->exon_strand_start = exon_strand_start;
	newFrag->exon_strand_end = exon_strand_end;


	for (int iLoop = 0; iLoop <= support_vec_size; iLoop++)
		newFrag->support[iLoop] = support[iLoop];

	newFrag->alter = NULL;
	//	newFrag->fivePalterTail = NULL;
	newFrag->coverage = NULL;
	newFrag->start_real = start_real;
	newFrag->end_real = end_real;

	return newFrag;
}

fragment::~fragment()
{
	alter_junction *delJunc;
	delJunc = alter;
	while (delJunc != NULL)
	{
		alter = delJunc->next;
		delete delJunc;
		delJunc = alter;
	}

	if (coverage != NULL)
		delete [] coverage;

	delete [] support;
}

spliceSite::spliceSite()
{
	position = 0;
	directionOut = true;
	strand = undetermined;
	is_ref = false;
}

rangeJunction::rangeJunction()
{
	junc = NULL;
	next = NULL;
}

RangeJunctionList::RangeJunctionList()
{
	rangeLow = 0;
	rangeHigh = 0;
	cnt_rangeJunction = 0;
	transDirection = undetermined;
	list = NULL;
	listtail = NULL;
	nextList = NULL;
}

RangeJunctionList* RangeJunctionList::clone()
{
	//clone a same RangeJunctionList
	RangeJunctionList *resultList;
	resultList = new RangeJunctionList;

	resultList->rangeLow = rangeLow;
	resultList->rangeHigh = rangeHigh;
	resultList->transDirection = transDirection;
	resultList->list = NULL;
	resultList->listtail = NULL;
	resultList->nextList = NULL;

	rangeJunction *curList, *newList;
	curList = list;
	while (curList != NULL)
	{
		newList = new rangeJunction;
		newList->junc = curList->junc;
		newList->next = NULL;

		if (resultList->list == NULL)
		{
			resultList->list = newList;
			resultList->listtail = newList;
		}
		else
		{
			resultList->listtail->next = newList;
			resultList->listtail = newList;
		}

		curList = curList->next;
	}

	return resultList;
}


void RangeJunctionList::insert_feature(rangeJunction *new_rangeJunction)
{
	//insert the new feature to tail of the feature list
	if (listtail == NULL)
	{
		assert(list == NULL);

		list = new_rangeJunction;
		listtail = new_rangeJunction;
	}
	else
	{
		assert(listtail->next == NULL);

		listtail->next = new_rangeJunction;
		listtail = new_rangeJunction;
	}

	++cnt_rangeJunction;

	return;
}

void RangeJunctionList::count_featurelist()
{
	//count number of features in the feature list
	cnt_rangeJunction = 0;
	rangeJunction *cur_feature = list;
	while (cur_feature != NULL)
	{
		++cnt_rangeJunction;
		cur_feature = cur_feature->next;
	}

	return;
}

RangeJunctionList::~RangeJunctionList()
{
	nextList = NULL;

	rangeJunction *delList;

	while (list != NULL)
	{
		delList = list;
		list = delList->next;
		delete delList;
	}
}

GTedge::GTedge()
{
	linkedVertex = NULL;
	next = NULL;
}

alternative_path::alternative_path()
{
	path_start = 0;
	path_end = 0;
	whole_path_start = 0;
	whole_path_end = 0;
	transDirection = undetermined;
	support = new double [SUPPORT_VECTOR_SIZE];
	proportion = new double [SUPPORT_VECTOR_SIZE];
	for (int i = 0; i < SUPPORT_VECTOR_SIZE; i++)
	{
		support[i] = 0.0;
		proportion[i] = 1.0;
	}
	junctionNum = 0;
	exonNum = 0;
	pathVertex = NULL;
	next = NULL;

	avg_support = 0.0;
}

alternative_path::~alternative_path()
{
	delete [] support;
	delete [] proportion;
}

GTvertex::GTvertex()
{
	level = 0;
	rangeLow = 0;
	rangeHigh = 0;
	child = NULL;
	childType = 0;
	childNum = 0;
	junctionInRange = NULL;
	junctionNum = 0;
	exonNum = 0;
	prevSibling = NULL;
	nextSibling = NULL;
	has_novel_junction = false;
	has_novel_junction_major = false;

	alterSpliceSite = NULL;

	support = new double [SUPPORT_VECTOR_SIZE];
	proportion = new double [SUPPORT_VECTOR_SIZE];
	MSE_estimation = new double [SUPPORT_VECTOR_SIZE];
	min_path_support = new double [SUPPORT_VECTOR_SIZE];
	obs_support = new double [SUPPORT_VECTOR_SIZE];
	for (int i = 0; i < SUPPORT_VECTOR_SIZE; i++)
	{
		support[i] = 0.0;
		proportion[i] = 1.0;
		MSE_estimation[i] = 0.0;
		min_path_support[i] = 0.0;
		obs_support[i] = 0.0;
	}
	
	estimated = false;
	path_extended = false;
	representative = NULL;
	estimate_exonNum = 0;

	ASMcategory = -1;

	major_alter_paths = NULL;
	major_alter_paths_num = 0;
}

GTvertex::~GTvertex()
{
	// 	if (junctionInRange != NULL)
	// 	{
	// 		delete junctionInRange;
	// 	}
	// 
	// 	delete representative;
	// 
	// 	//... delete children

	delete [] support;
	delete [] proportion;
	delete [] MSE_estimation;
	delete [] min_path_support;
	delete [] obs_support;
}

// GenomeTree::GenomeTree()
// {
// 	root = NULL;
// }
// 
// GenomeTree::~GenomeTree()
// {
// 	//
// }


//compare two fragments and return if fragment 1 has smaller genomic coordinate than fragment 2
bool comp_frag_ptr(const fragment *frag_1, const fragment *frag_2) {
	assert(frag_1 && frag_2);
	return frag_1->start < frag_2->start || (frag_1->start == frag_2->start && frag_1->end < frag_2->end);
}

// compare two splice sites, with the possibility of having two sites at the same position but with different directionOut
bool comp_ssite_withdup(const spliceSite &site1, const spliceSite &site2) {
	assert(site1.chromosome.compare(site2.chromosome) == 0);
	return site1.position < site2.position || (site1.position == site2.position && !site1.directionOut && site2.directionOut);
}
