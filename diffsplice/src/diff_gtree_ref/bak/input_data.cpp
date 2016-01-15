#include "input_data.h"
#include "common_function.h"
#include "splice_graph.h"
#include "cut_exon_bound.h"

vector<long> chromosome_expr; //array for counting expression on the chromosomes

vector <fragment*> junction_array;

unsigned long fragmentID_Cnt = 0;


//clear a vector and release its memory
template <class T>
void free_vector(T &t)
{
	T tmp;
	tmp.swap(t);
}

//add a candidate junction 
void input_junction_merged(string chromosome, string filename_junc_merged)
{
	fragment *new_junction;
	long occurrence, position_start, position_end, sample_cnt;
	int field_xs;
	double tmp_val;
	string info;

	ifstream juncfile(filename_junc_merged.c_str());
	if (!juncfile.is_open())
	{
		cout << "error: fail to open merged junction file" << endl;
		exit(1);
	}

	while (juncfile >> occurrence)
	{
		juncfile >> position_start;
		juncfile >> position_end;
		juncfile >> field_xs;

		position_start += READ_COORDINATE_OFFSET;
		position_end += READ_COORDINATE_OFFSET;

		if (position_end - position_start - 1 < min_junction_intron_length)
		{
			getline(juncfile, info); 
			continue; //filter small indels
		}

		if (junction_array.size() >= junction_array.capacity())
			junction_array.reserve(junction_array.capacity() + default_max_junction_num);

#ifdef UNIX
		new_junction = new fragment (SUPPORT_VECTOR_SIZE);
#else
		new_junction = new fragment (debug_total_sample_num);
#endif
		// 11/12/2014 rectify the splice junction strand
		if (!junction_array.empty() && position_start == junction_array.back()->start && position_end == junction_array.back()->end) {
#ifdef UNIX
			if (occurrence > junction_array.back()->support[SUPPORT_VECTOR_SIZE]) {
				junction_array.back()->transDirection = field_xs == 0 ? antisense : sense;
			}

			for (sample_cnt = 0; sample_cnt < SUPPORT_VECTOR_SIZE; ++sample_cnt) {
				juncfile >> tmp_val;
				junction_array.back()->support[sample_cnt] += tmp_val;
			}
			junction_array.back()->support[SUPPORT_VECTOR_SIZE] += occurrence;
#else
			if (occurrence > junction_array.back()->support[debug_total_sample_num]) {
				junction_array.back()->transDirection = field_xs == 0 ? antisense : sense;
			}
			for (sample_cnt = 0; sample_cnt < debug_total_sample_num; ++sample_cnt) {
				juncfile >> tmp_val;
				junction_array.back()->support[sample_cnt] += tmp_val;
			}
			junction_array.back()->support[debug_total_sample_num] += occurrence;
#endif

		}
		else {
			new_junction->type = frag_junction;
			new_junction->frag_id = ++fragmentID_Cnt;
			new_junction->chromosome_start = chromosome;
			new_junction->chromosome_end = chromosome;
			new_junction->start = position_start;
			new_junction->end = position_end;
			if (field_xs == 0)
				new_junction->transDirection = antisense;
			else
				new_junction->transDirection = sense;

#ifdef UNIX
			for (sample_cnt = 0; sample_cnt < SUPPORT_VECTOR_SIZE; ++sample_cnt)
				juncfile >> new_junction->support[sample_cnt];
			new_junction->support[SUPPORT_VECTOR_SIZE] = occurrence;
#else
			for (sample_cnt = 0; sample_cnt < debug_total_sample_num; ++sample_cnt)
				juncfile >> new_junction->support[sample_cnt];
			new_junction->support[debug_total_sample_num] = occurrence;
#endif

			junction_array.push_back(new_junction);
		}
		
		getline(juncfile, info);
	}

	juncfile.close();
	return;
}


// return true if a junction passes the group-based junction filter
bool pass_filter_group(const double *support)
{
	//vector <double> group_mean_support (global_total_group_num+1, 0.0);
	bool pass_filter = false;

	for (int groupLoopCnt = 1; groupLoopCnt <= global_total_group_num; ++groupLoopCnt)
	{
		assert(global_group_size[groupLoopCnt] > 0);

		double groupSum_expr = 0.0;
		for (int sampleLoopCnt = 0; sampleLoopCnt < global_group_size[groupLoopCnt]; ++sampleLoopCnt)
		{
			groupSum_expr += support[global_group_base_index[groupLoopCnt] + sampleLoopCnt];
		}
		groupSum_expr /= global_group_size[groupLoopCnt];

		if (groupSum_expr >= thresh_junctionfilter_groupwise)
		{
			pass_filter = true;
			break;
		}
	}

	return pass_filter;
}

// return true if a junction passes the no_group-based junction filter
bool pass_filter_no_group(const double *support)
{
	//vector <double> group_mean_support (global_total_group_num+1, 0.0);
	bool pass_filter = false;

	int num_present_sample = 0;
	for (long sampleLoop = 0; sampleLoop < SUPPORT_VECTOR_SIZE_WHOLE_DATA; ++sampleLoop)
	{
		if (support[sampleLoop] >= thresh_present_support)
		{
			++num_present_sample;
		}
	}

	if (num_present_sample >= thresh_num_present_samples){
		pass_filter = true;
	}
	return pass_filter;
}


//select junctions with the following criteria:
// 1. using annotation as a white list enabled and the junction in the given junction annotation or
// 2. pass the expression-based filter according to the mode and threshold specified by the user, specifically
// 2.1 mode 'group', at least one sample group should have an averaged expression >= threshold
// 2.2 mode 'nogroup', junction is present (with support >= the given threshold) in at least the given number of samples
void label_select_junctions(string filename_annotation)
{
	//vector <bool> junction_in_annotation (junction_array.size(), false); //consider as in annotation as default to tolerate the case when the annotation is not available

	unsigned long junctionLoop = 0, cnt_junction_in_anno = 0;

	bool use_white_list = true;
	if (junction_annotation_as_white_list.compare("no") == 0){
		use_white_list = false;
	}

	if (!filename_annotation.empty()) {
		ifstream infile(filename_annotation.c_str());

		if (infile.is_open())
		{
			setting_junc_anno_provided = true;

			long start_pos, end_pos;
			string info, chrname;

			while (infile >> chrname)
			{
				infile >> start_pos;
				infile >> end_pos;
				getline(infile, info);

				while (junctionLoop < junction_array.size() && (junction_array[junctionLoop]->start < start_pos || (junction_array[junctionLoop]->start == start_pos && junction_array[junctionLoop]->end < end_pos)))
				{
					//junction_in_annotation[junctionLoop] = false;
					++junctionLoop;
				}

				if (junctionLoop >= junction_array.size())
					break;

				if (junction_array[junctionLoop]->start == start_pos && junction_array[junctionLoop]->end == end_pos)
				{
					//find a match in junction annotation
					//junction_in_annotation[junctionLoop] = true;
					junction_array[junctionLoop]->in_junc_annotation = true;
					++junctionLoop;
					++cnt_junction_in_anno;
				}
			}

			infile.close();
		}
		else
		{
			cout << "warning: junction annotation file has been specified but cannot be opened." << endl;
		}
	}

#ifdef FILTER_JUNCTION
	vector <fragment*> junction_array_filtered;

	for (junctionLoop = 0; junctionLoop < junction_array.size(); ++junctionLoop)
	{
		bool pass_filter = true;

		if (junctionfilter_mode == "group"){
			pass_filter = pass_filter_group(junction_array[junctionLoop]->support);
		}
		else if (junctionfilter_mode == "nogroup"){
			pass_filter = pass_filter_no_group(junction_array[junctionLoop]->support);
		}
		else if (junctionfilter_mode == "none"){
			// do nothing
		}
		else {
			cout << "Error: unrecognized junction filter mode" << endl; // should never happen, because unrecognized mode will be automatically set as 'none'
			exit(1);
		}

		//if (num_present_sample >= thresh_num_present_samples || junction_in_annotation[junctionLoop])
		if (pass_filter || (use_white_list && junction_array[junctionLoop]->in_junc_annotation))
		{
			//select this junction
			junction_array_filtered.push_back(junction_array[junctionLoop]);
#ifndef UNIX
			//only keep the samples under debug
			double *tmp_support = new double [SUPPORT_VECTOR_SIZE+1];
			for (int sampleLoop = 0; sampleLoop < SUPPORT_VECTOR_SIZE; ++sampleLoop)
			{
				tmp_support[sampleLoop] = junction_array[junctionLoop]->support[debug_sample_index_offset+sampleLoop];
			}
			tmp_support[SUPPORT_VECTOR_SIZE] = junction_array[junctionLoop]->support[debug_total_sample_num];
			delete [] junction_array[junctionLoop]->support;
			junction_array[junctionLoop]->support = tmp_support;
#endif
		}
		else
		{
			//discard this junction
			delete junction_array[junctionLoop];
		}
	}

	//save the selected junctions
	junction_array.swap(junction_array_filtered);
	junction_array_filtered.clear(); 
	free_vector(junction_array_filtered);
#endif

	return;
}


void extract_exons(RangeJunctionList *cur_range, string chromosome)
{
	spliceSite *curSite;
	fragment *new_feature;
	rangeJunction *new_rangeJunction;
	unsigned long tmpLoop;

	if (junction_array.size() == 0)
	{
		new_feature = new fragment(SUPPORT_VECTOR_SIZE);
		new_feature->type = frag_exon;
		new_feature->altersite = 0;
		new_feature->frag_name = "Exon_1";
		new_feature->frag_id = ++fragmentID_Cnt;

		new_feature->start = cur_range->rangeLow;
		new_feature->end = cur_range->rangeHigh;
		new_feature->chromosome_start = chromosome;
		new_feature->chromosome_end = chromosome;

		//insert exonic fragments to list
		new_rangeJunction = new rangeJunction;
		new_rangeJunction->junc = new_feature;

		cur_range->insert_feature(new_rangeJunction);

		return;
	}

	spliceSite** splicesite_list;
	splicesite_list = new spliceSite* [junction_array.size() * 2 + 2];
	unsigned long splicesite_list_cnt = 0;

	//input splice sites
	for (tmpLoop = 0; tmpLoop < junction_array.size(); ++tmpLoop)
	{
		curSite = new spliceSite;		
		curSite->chromosome = junction_array[tmpLoop]->chromosome_start;
		curSite->position = junction_array[tmpLoop]->start;
		curSite->directionOut = true;
		curSite->strand = junction_array[tmpLoop]->transDirection;
		splicesite_list[++splicesite_list_cnt] = curSite;

		curSite = new spliceSite;		
		curSite->chromosome = junction_array[tmpLoop]->chromosome_end;
		curSite->position = junction_array[tmpLoop]->end;
		curSite->directionOut = false;
		curSite->strand = junction_array[tmpLoop]->transDirection;
		splicesite_list[++splicesite_list_cnt] = curSite;
	}

	///////////////////////////////
	//sort splice sites
	void** sortlist_splicesite = new void* [junction_array.size() * 2 + 2];
	double* sortkey_splicesite = new double [junction_array.size() * 2 + 2];

	for (tmpLoop = 1; tmpLoop <= splicesite_list_cnt; ++tmpLoop)
	{
		sortkey_splicesite[tmpLoop] = splicesite_list[tmpLoop]->position;
		sortlist_splicesite[tmpLoop] = (void*) splicesite_list[tmpLoop];
	}

	mergeSort(sortlist_splicesite, sortkey_splicesite, splicesite_list_cnt);

	for (tmpLoop = 1; tmpLoop <= splicesite_list_cnt; ++tmpLoop)
		splicesite_list[tmpLoop] = (spliceSite*) sortlist_splicesite[tmpLoop];

	delete [] sortlist_splicesite;
	delete [] sortkey_splicesite;

	///////////////////////////////
	//build exonic fragments from adjacent splice sites 
	//exonic if not out-in

	//first exon
	new_feature = new fragment(SUPPORT_VECTOR_SIZE);
	new_feature->type = frag_exon;
	new_feature->altersite = 0;
	new_feature->frag_name = "Exon_1";
	new_feature->frag_id = ++fragmentID_Cnt;

	curSite = splicesite_list[1];
	new_feature->chromosome_start = curSite->chromosome;
	new_feature->chromosome_end = curSite->chromosome;
	new_feature->start = cur_range->rangeLow;
	new_feature->end = curSite->position;
	new_feature->exon_strand_start = terminal;
	new_feature->exon_strand_end = curSite->strand;
	
	//insert exonic fragments to list
	new_rangeJunction = new rangeJunction;
	new_rangeJunction->junc = new_feature;

	cur_range->insert_feature(new_rangeJunction);


	for (tmpLoop = 1; tmpLoop < splicesite_list_cnt; ++tmpLoop)
	{
		if (splicesite_list[tmpLoop]->position == splicesite_list[tmpLoop+1]->position && splicesite_list[tmpLoop]->directionOut == splicesite_list[tmpLoop+1]->directionOut)
		{
			delete splicesite_list[tmpLoop];
		} 
		else
		{
			new_feature = new fragment(SUPPORT_VECTOR_SIZE);
			new_feature->type = frag_exon;

			if (splicesite_list[tmpLoop]->directionOut == false && splicesite_list[tmpLoop+1]->directionOut == true)
				new_feature->altersite = 0;
			else if (splicesite_list[tmpLoop]->directionOut == true && splicesite_list[tmpLoop+1]->directionOut == true)
				new_feature->altersite = 1;
			else if (splicesite_list[tmpLoop]->directionOut == false && splicesite_list[tmpLoop+1]->directionOut == false)
				new_feature->altersite = 2;
			else if (splicesite_list[tmpLoop]->directionOut == true && splicesite_list[tmpLoop+1]->directionOut == false)
			{
				new_feature->type = frag_retained_intron;
				new_feature->altersite = -1;
			}

			new_feature->frag_name = "Exon_" + itostr(tmpLoop/2 + 2);
			new_feature->frag_id = ++fragmentID_Cnt;

			curSite = splicesite_list[tmpLoop];
			new_feature->chromosome_start = curSite->chromosome;
			if (curSite->directionOut == true)
				new_feature->start = curSite->position + 1;
			else
				new_feature->start = curSite->position;
			new_feature->exon_strand_start = curSite->strand;
			delete curSite;

			curSite = splicesite_list[tmpLoop+1];
			new_feature->chromosome_end = curSite->chromosome;
			if (curSite->directionOut == true)
				new_feature->end = curSite->position;
			else
				new_feature->end = curSite->position - 1;
			new_feature->exon_strand_end = curSite->strand;

			if (new_feature->end >= new_feature->start)
			{
				//insert exonic fragments to list
				new_rangeJunction = new rangeJunction;
				new_rangeJunction->junc = new_feature;
				cur_range->insert_feature(new_rangeJunction);
			}
			else
			{
				delete new_feature;
			}
		}
	}

	//last exon
	new_feature = new fragment(SUPPORT_VECTOR_SIZE);
	new_feature->type = frag_exon;
	new_feature->altersite = 0;
	new_feature->frag_name = "Exon_" + itostr(splicesite_list_cnt/2 + 2);
	new_feature->frag_id = ++fragmentID_Cnt;

	curSite = splicesite_list[splicesite_list_cnt];
	new_feature->chromosome_start = curSite->chromosome;
	new_feature->chromosome_end = curSite->chromosome;
	new_feature->start = splicesite_list[splicesite_list_cnt]->position;
	new_feature->end = cur_range->rangeHigh;
	new_feature->exon_strand_start = curSite->strand;
	new_feature->exon_strand_end = terminal;
	delete curSite;
	
	//insert exonic fragments to list
	new_rangeJunction = new rangeJunction;
	new_rangeJunction->junc = new_feature;
	cur_range->insert_feature(new_rangeJunction);

	delete [] splicesite_list;

	return;
}


//compute relative change ratio at the curPos-th position of the region ([startIndex+curPos, endIndex]/[startIndex, startIndex+curPos)) defined by [startIndex, endIndex], curPos starts at 0
double windowChangeRatio(long startIndex, long endIndex, long curPos)
{
	double ratio = 1.0;
	long vectorLength = endIndex - startIndex + 1;

	if (curPos <= 0)
		return 2*COVERAGE_CHANGE_THRESH_EXON; // 2 * COVERAGE_CHANGE_THRESH;
	else if (curPos >= vectorLength - 1)
		return 0.0;

	double curSupp_left = 0, curSupp_right = 0, maxSupp_left, minSupp_left, maxSupp_right, minSupp_right; //support of the left window and the right window
	int windowSize_left = 0, windowSize_right = 0; //size of the left window and the right window 

	for (maxSupp_left = minSupp_left = chromosome_expr[startIndex + curPos - 1]; windowSize_left < COVERAGE_CHANGE_WINDOW && curPos - windowSize_left > 0; ++windowSize_left)
	{
		curSupp_left += chromosome_expr[startIndex + curPos - 1 - windowSize_left];
		if (chromosome_expr[startIndex + curPos - 1 - windowSize_left] > maxSupp_left)
			maxSupp_left = chromosome_expr[startIndex + curPos - 1 - windowSize_left];
		if (chromosome_expr[startIndex + curPos - 1 - windowSize_left] < minSupp_left)
			minSupp_left = chromosome_expr[startIndex + curPos - 1 - windowSize_left];
	}

	for (maxSupp_right = minSupp_right = chromosome_expr[startIndex + curPos]; windowSize_right < COVERAGE_CHANGE_WINDOW && curPos + windowSize_right < vectorLength; ++windowSize_right)
	{
		curSupp_right += chromosome_expr[startIndex + curPos + windowSize_right];
		if (chromosome_expr[startIndex + curPos + windowSize_right] > maxSupp_right)
			maxSupp_right = chromosome_expr[startIndex + curPos + windowSize_right];
		if (chromosome_expr[startIndex + curPos + windowSize_right] < minSupp_right)
			minSupp_right = chromosome_expr[startIndex + curPos + windowSize_right];
	}

	if (curSupp_left < 1e-10)
		ratio = 2*COVERAGE_CHANGE_THRESH_EXON; // 2 * COVERAGE_CHANGE_THRESH;
	else
		ratio = (curSupp_right / windowSize_right) / (curSupp_left / windowSize_left);

	if (minSupp_left > maxSupp_right && ratio < 1 && curSupp_left/windowSize_left >= coverageThreshold_exon)
		return ratio;
	else if (maxSupp_left < minSupp_right && ratio > 1 && curSupp_right/windowSize_right >= coverageThreshold_exon)
		return ratio;
	else
		return 1.0;
}


double cutAlterSite(long startIndex, long endIndex, int altersitetype, long &totalSupport, long &cutFragSupport, long &cutFragStart, long &cutFragEnd)
{
	double maxRatio = 0, supp_left, supp_right, curRatio;
	long size_left, size_right, curPos, vectorLength = endIndex - startIndex + 1;
	totalSupport = 0; cutFragSupport = 0; cutFragStart = startIndex; cutFragEnd = endIndex;

	if (altersitetype == 1)
	{
		for (size_left = 0, supp_left = 0.0; size_left < vectorLength - 1; ++size_left)
			supp_left += chromosome_expr[startIndex + size_left];

		size_right = 1;
		supp_right = chromosome_expr[endIndex];

		totalSupport = long(supp_left + supp_right);

		for (curPos = vectorLength - 1; curPos > 0; --curPos)
		{
			if (supp_left < 1e-10)
			{
				if (curPos < MAX_NOCOVERAGE_LENGTH)
				{
					//tolerate this gap
					return maxRatio>1? maxRatio : 0;
				} 
				else
				{
					//return this position
					cutFragStart = startIndex + curPos;
					cutFragSupport = supp_right;
					return COVERAGE_CHANGE_THRESH_ALTER_SITE + 1; //definitely accept
				}				
			}
			else
				curRatio = (supp_right / size_right) / (supp_left / size_left);

			if (curRatio > maxRatio)
			{
				maxRatio = curRatio;
				cutFragStart = startIndex + curPos;				
				cutFragSupport = long(supp_right);
			}

			--size_left;
			supp_left -= chromosome_expr[startIndex + curPos - 1];
			++size_right;
			supp_right += chromosome_expr[startIndex + curPos - 1];
		}

		return maxRatio>1? maxRatio : 0;
	}
	else if (altersitetype == 2)
	{
		size_left = 0;
		supp_left = 0.0;

		for (size_right = 0, supp_right = 0.0; size_right < vectorLength; ++size_right)
			supp_right += chromosome_expr[startIndex + size_right];

		totalSupport = long(supp_left + supp_right);

		for (curPos = 0; curPos < vectorLength; ++curPos)
		{
			if (supp_right < 1e-10)
			{
				if (vectorLength - curPos <= MAX_NOCOVERAGE_LENGTH)
				{
					//tolerate this gap
					return maxRatio>1? maxRatio : 0;
				} 
				else
				{
					//return this position
					cutFragEnd = startIndex + curPos - 1;
					cutFragSupport = supp_left;
					return COVERAGE_CHANGE_THRESH_ALTER_SITE + 1; //definitely accept
				}				
			}
			else if (size_left == 0)
				curRatio = 0;
			else
				curRatio = (supp_left / size_left) / (supp_right / size_right);

			if (curRatio > maxRatio)
			{
				maxRatio = curRatio;
				cutFragEnd = startIndex + curPos - 1;
				cutFragSupport = long(supp_left);
			}

			++size_left;
			supp_left += chromosome_expr[startIndex + curPos];
			--size_right;
			supp_right -= chromosome_expr[startIndex + curPos];
		}

		return maxRatio>1? maxRatio : 0;
	}
	else
	{
		return 0;
	}
}


void trim_exons(RangeJunctionList *cur_range, int index)
{
	//trim exonic regions based on expression

	long startPosition, tmp, tmpVectorLength, iLoop, curFragStart, curSupport, tmpExonAddLength, altersiteFragStart, altersiteFragEnd, altersiteFragSupport, max_nocoverage_length_in_use = MAX_NOCOVERAGE_LENGTH;
	double changeRatio, thresh_exon_cov_in_use = coverageThreshold_exon;
	rangeJunction *prev_feature, *cur_feature, *del_feature, *new_feature, *feature_addlist = NULL, *feature_addlist_tail = NULL;
	fragment *new_fragment;
	bool makingNewFrag, largeGap;
	startPosition = 0;

	prev_feature = NULL;
	cur_feature = cur_range->list;
	while (cur_feature != NULL)
	{
		//trim current feature

		// 		if (cur_feature->pos_start <= 26138276 && cur_feature->pos_end >= 26138276 || cur_feature->pos_start <= 26140391 && cur_feature->pos_end > 26140391)
		// 			cout << "a";

		// handle the case for reference-based reconstruction
		if (cur_feature->junc->type == frag_exon && cur_feature->junc->in_exon_annotation == annotated) {
			prev_feature = cur_feature;
			cur_feature = cur_feature->next;
			continue;
		}

		curSupport = 0;
		makingNewFrag = false;
		curFragStart = cur_feature->junc->start;
		tmpVectorLength = cur_feature->junc->end - cur_feature->junc->start + 1;
		feature_addlist = NULL;
		feature_addlist_tail = NULL;

		if (cur_feature->junc->type == frag_exon && cur_feature->junc->altersite > 0 && tmpVectorLength <= MAX_ALTER_SPLICE_SITE_LENGTH) // 10/13/2013 add a maximum length, longer exons should go to the exon cut. important in real data, e.g. chr1:113946277-114059752 in brca
		{
			changeRatio = cutAlterSite(cur_feature->junc->start, cur_feature->junc->end, cur_feature->junc->altersite, curSupport, altersiteFragSupport, altersiteFragStart, altersiteFragEnd);

			if (changeRatio < COVERAGE_CHANGE_THRESH_ALTER_SITE)
			{
				//alternative splice site
				tmpExonAddLength = tmpVectorLength;
				if (tmpExonAddLength >= MIN_ALTER_SPLICE_SITE_LENGTH && curSupport >= coverageThreshold_exon * tmpExonAddLength * SUPPORT_VECTOR_SIZE_WHOLE_DATA)
				{
					//get an end, make a new fragment
					new_fragment = new fragment(SUPPORT_VECTOR_SIZE);
					new_fragment->type = cur_feature->junc->type; //frag_exon;
					new_fragment->frag_name = "ExonAdd_" + itostr(fragmentID_Cnt + 1);
					new_fragment->frag_id = ++fragmentID_Cnt;
					new_fragment->altersite = cur_feature->junc->altersite;

					new_fragment->chromosome_start = cur_feature->junc->chromosome_start;
					new_fragment->chromosome_end = cur_feature->junc->chromosome_end;
					new_fragment->start = cur_feature->junc->start;
					new_fragment->end = cur_feature->junc->end;
					new_fragment->support[index] += curSupport;

					new_fragment->exon_strand_start = cur_feature->junc->exon_strand_start;
					new_fragment->exon_strand_end = cur_feature->junc->exon_strand_end;

					new_feature = new rangeJunction;
					new_feature->junc = new_fragment;

					if (feature_addlist == NULL)
					{
						feature_addlist = new_feature;
						feature_addlist_tail = new_feature;
					}
					else
					{
						feature_addlist_tail->next = new_feature;
						feature_addlist_tail = new_feature;
					}					 
				}
			}
			else
			{
				//alternative start/end
				tmpExonAddLength = altersiteFragEnd - altersiteFragStart + 1;
				if (tmpExonAddLength >= MIN_EXON_LENGTH && altersiteFragSupport >= coverageThreshold_exon * tmpExonAddLength * SUPPORT_VECTOR_SIZE_WHOLE_DATA)
				{
					//get an end, make a new fragment
					new_fragment = new fragment(SUPPORT_VECTOR_SIZE);
					new_fragment->type = cur_feature->junc->type; //frag_exon;
					new_fragment->frag_name = "ExonAdd_" + itostr(fragmentID_Cnt + 1);
					new_fragment->frag_id = ++fragmentID_Cnt;

					new_fragment->chromosome_start = cur_feature->junc->chromosome_start;
					new_fragment->chromosome_end = cur_feature->junc->chromosome_end;
					new_fragment->start = altersiteFragStart;
					new_fragment->end = altersiteFragEnd;
					new_fragment->support[index] += altersiteFragSupport;

					if (cur_feature->junc->altersite == 1)
					{
						new_fragment->exon_strand_start = terminal;
						new_fragment->exon_strand_end = cur_feature->junc->exon_strand_end;
					}
					else if (cur_feature->junc->altersite == 2)
					{
						new_fragment->exon_strand_start = cur_feature->junc->exon_strand_start;
						new_fragment->exon_strand_end = terminal;
					}

					new_feature = new rangeJunction;
					new_feature->junc = new_fragment;

					if (feature_addlist == NULL)
					{
						feature_addlist = new_feature;
						feature_addlist_tail = new_feature;
					}
					else
					{
						feature_addlist_tail->next = new_feature;
						feature_addlist_tail = new_feature;
					}
				}
			}
		}
		else if (cur_feature->junc->type == frag_exon)
		{
			double max_cov_this_exon = 0;
			for (iLoop = 0; iLoop < tmpVectorLength; ++iLoop)
			{
				if (chromosome_expr[cur_feature->junc->start + iLoop] / double(SUPPORT_VECTOR_SIZE_WHOLE_DATA) > max_cov_this_exon)
				{
					max_cov_this_exon = chromosome_expr[cur_feature->junc->start + iLoop] / double(SUPPORT_VECTOR_SIZE_WHOLE_DATA);
				}
			}
			if (max_cov_this_exon > MAX_NOCOVERAGE_THRESH_LOW_COV_REGION)
				max_nocoverage_length_in_use = MAX_NOCOVERAGE_LENGTH;
			else
				max_nocoverage_length_in_use = MAX_NOCOVERAGE_LENGTH_LOW_COV_REGION;
			
			thresh_exon_cov_in_use = max_cov_this_exon * 0.01;
			if (thresh_exon_cov_in_use < coverageThreshold_exon)
				thresh_exon_cov_in_use = coverageThreshold_exon;

			for (iLoop = 0; iLoop < tmpVectorLength; ++iLoop)
			{
				if (iLoop + max_nocoverage_length_in_use >= tmpVectorLength)
				{
					largeGap = true; //old - false; always cut short gaps 
				} 
				else
				{
					largeGap = true;
					for (tmp = 0; iLoop + tmp < tmpVectorLength && tmp < max_nocoverage_length_in_use; ++tmp)
					{
						if (chromosome_expr[cur_feature->junc->start + iLoop + tmp] > thresh_exon_cov_in_use * SUPPORT_VECTOR_SIZE_WHOLE_DATA)
						{
							largeGap = false;
							break;
						}
					}
				}

				changeRatio = windowChangeRatio(cur_feature->junc->start, cur_feature->junc->end, iLoop);

				//if ((chromosome_expr[cur_feature->junc->start + iLoop] <= coverageThreshold_exon * SUPPORT_VECTOR_SIZE && largeGap == true) || iLoop == tmpVectorLength - 1)
				//5/15/2013 add one condition, cut boundary at sharp changes also
				if (((chromosome_expr[cur_feature->junc->start + iLoop] <= thresh_exon_cov_in_use * SUPPORT_VECTOR_SIZE_WHOLE_DATA || changeRatio < 1/COVERAGE_CHANGE_THRESH_EXON) && largeGap == true) || iLoop == tmpVectorLength - 1)
				{
					if (makingNewFrag == true)
					{
						tmpExonAddLength = cur_feature->junc->start + iLoop - curFragStart + 1;
						if (tmpExonAddLength >= MIN_EXON_LENGTH && curSupport >= thresh_exon_cov_in_use * tmpExonAddLength * SUPPORT_VECTOR_SIZE_WHOLE_DATA)
						{
							//get an end, make a new fragment
							new_fragment = new fragment(SUPPORT_VECTOR_SIZE);
							new_fragment->type = cur_feature->junc->type;
							new_fragment->frag_name = "ExonAdd_" + itostr(fragmentID_Cnt + 1);
							new_fragment->frag_id = ++fragmentID_Cnt;

							new_fragment->chromosome_start = cur_feature->junc->chromosome_start;
							new_fragment->chromosome_end = cur_feature->junc->chromosome_end;
							new_fragment->start = curFragStart;
							new_fragment->end = cur_feature->junc->start + iLoop;
							new_fragment->support[index] += curSupport;

							if (new_fragment->start == cur_feature->junc->start)
								new_fragment->exon_strand_start = cur_feature->junc->exon_strand_start;
							else
								new_fragment->exon_strand_start = terminal;
							if (new_fragment->end == cur_feature->junc->end)
								new_fragment->exon_strand_end = cur_feature->junc->exon_strand_end;
							else
								new_fragment->exon_strand_end = terminal;

							new_feature = new rangeJunction;
							new_feature->junc = new_fragment;

							if (feature_addlist == NULL)
							{
								feature_addlist = new_feature;
								feature_addlist_tail = new_feature;
							}
							else
							{
								feature_addlist_tail->next = new_feature;
								feature_addlist_tail = new_feature;
							}
						}

						curSupport = 0;
						makingNewFrag = false;								
					}
					else
						curSupport += chromosome_expr[cur_feature->junc->start + iLoop];
				}
				else
				{
					//if (makingNewFrag == false && chromosome_expr[cur_feature->junc->start + iLoop] > coverageThreshold_exon)// && changeRatio > COVERAGE_CHANGE_THRESH)
					//5/15/2013 add the changeRatio condition
					if (makingNewFrag == false && chromosome_expr[cur_feature->junc->start + iLoop] > thresh_exon_cov_in_use)// * SUPPORT_VECTOR_SIZE_WHOLE_DATA && changeRatio > COVERAGE_CHANGE_THRESH_EXON)
						//if (makingNewFrag == false)
					{
						curFragStart = cur_feature->junc->start + iLoop;
						makingNewFrag = true;
						curSupport = 0;
					}
					curSupport += chromosome_expr[cur_feature->junc->start + iLoop];
				}
			}
		}
		else if (cur_feature->junc->type == frag_retained_intron)
		{
			for (iLoop = 0; iLoop < tmpVectorLength; ++iLoop)
			{
				if (iLoop + MAX_NOCOVERAGE_LENGTH >= tmpVectorLength)
				{
					largeGap = true; //old - false; always cut short gaps 
				} 
				else
				{
					largeGap = true;
					for (tmp = 0; iLoop + tmp < tmpVectorLength && tmp < MAX_NOCOVERAGE_LENGTH; ++tmp)
					{
						if (chromosome_expr[cur_feature->junc->start + iLoop + tmp] > coverageThreshold_exon * SUPPORT_VECTOR_SIZE_WHOLE_DATA)
						{
							largeGap = false;
							break;
						}
					}
				}

				changeRatio = windowChangeRatio(cur_feature->junc->start, cur_feature->junc->end, iLoop);

				if (((changeRatio < 1/COVERAGE_CHANGE_THRESH_EXON || chromosome_expr[cur_feature->junc->start + iLoop] <= coverageThreshold_intron * SUPPORT_VECTOR_SIZE_WHOLE_DATA) && largeGap == true) || iLoop == tmpVectorLength - 1)
				{
					if (makingNewFrag == true)
					{
						tmpExonAddLength = cur_feature->junc->start + iLoop - curFragStart + 1;
						if (tmpExonAddLength >= MIN_EXON_LENGTH && curSupport >= coverageThreshold_intron * tmpExonAddLength * SUPPORT_VECTOR_SIZE_WHOLE_DATA) //&& tmpExonAddLength < tmpVectorLength - 1
						{
							//get an end, make a new fragment
							new_fragment = new fragment(SUPPORT_VECTOR_SIZE);
							new_fragment->type = cur_feature->junc->type;
							new_fragment->frag_name = "ExonAdd_" + itostr(fragmentID_Cnt + 1);
							new_fragment->frag_id = ++fragmentID_Cnt;

							new_fragment->chromosome_start = cur_feature->junc->chromosome_start;
							new_fragment->chromosome_end = cur_feature->junc->chromosome_end;
							new_fragment->start = curFragStart;
							new_fragment->end = cur_feature->junc->start + iLoop;
							new_fragment->support[index] += curSupport;

							if (new_fragment->start == cur_feature->junc->start)
								new_fragment->exon_strand_start = cur_feature->junc->exon_strand_start;
							else
								new_fragment->exon_strand_start = terminal;
							if (new_fragment->end == cur_feature->junc->end)
								new_fragment->exon_strand_end = cur_feature->junc->exon_strand_end;
							else
								new_fragment->exon_strand_end = terminal;

							new_feature = new rangeJunction;
							new_feature->junc = new_fragment;

							if (feature_addlist == NULL)
							{
								feature_addlist = new_feature;
								feature_addlist_tail = new_feature;
							}
							else
							{
								feature_addlist_tail->next = new_feature;
								feature_addlist_tail = new_feature;
							}
						}

						curSupport = 0;
						makingNewFrag = false;								
					}
				}
				else
				{
					if (makingNewFrag == false && changeRatio > COVERAGE_CHANGE_THRESH_EXON && chromosome_expr[cur_feature->junc->start + iLoop] >= coverageThreshold_intron * SUPPORT_VECTOR_SIZE_WHOLE_DATA)
						//if (makingNewFrag == false)
					{
						curFragStart = cur_feature->junc->start + iLoop;
						makingNewFrag = true;
						curSupport = 0;
					}
					curSupport += chromosome_expr[cur_feature->junc->start + iLoop];
				}
			}
		}

		//insert added features into list
		if (feature_addlist != NULL)
		{
			//delete cur_feature
			//add new features
			if (prev_feature == NULL)
			{
				cur_range->list = feature_addlist;
				feature_addlist_tail->next = cur_feature->next;
				prev_feature = feature_addlist_tail;
			}
			else
			{
				prev_feature->next = feature_addlist;
				feature_addlist_tail->next = cur_feature->next;
				prev_feature = feature_addlist_tail;
			}

			if (cur_feature == cur_range->listtail)
			{
				cur_range->listtail = feature_addlist_tail;
			}
		}
		else
		{
			//delete cur_feature
			//prev_feature is unchanged
			if (prev_feature == NULL)
			{
				cur_range->list = cur_feature->next;
				if (cur_feature == cur_range->listtail)
					cur_range->listtail = cur_range->list;
			}
			else
			{
				prev_feature->next = cur_feature->next;
				if (cur_feature == cur_range->listtail)
					cur_range->listtail = prev_feature;
			}

		}

		del_feature = cur_feature;
		cur_feature = cur_feature->next;
		delete del_feature;
	}

	return;
}

#ifdef UNIX
// this function has two jobs: 1. to refine the boundary of a transcription start/end exon; 2. to separate mixed exon of two genes
void cut_exon_bound(RangeJunctionList *list_exon)
{
	rangeJunction *cur_feature;
	int num_cut_pt, exon_code;
	vector <double> exon_coverage, avg_coverage;
	vector <long> cut_pt;

#ifdef DEBUG_OUTPUT
	ofstream debug_out_cut_bound("debug_output/debug_cut_boundary.txt", fstream::app);
#endif
	
	cur_feature = list_exon->list;
	while (cur_feature != NULL)
	{
		exon_code = -1;
		if (cur_feature->junc->type == frag_exon && cur_feature->junc->altersite == 0 && (cur_feature->junc->end - cur_feature->junc->start + 1 <= 50000)) // 10/13/2013 add a maximum length, longer exons will need too much memory (O(n^2))
		{
			if (cur_feature->junc->exon_strand_start == terminal && cur_feature->junc->exon_strand_end != terminal)
			{	// transcription start exon of a gene
				exon_code = 1;
				num_cut_pt = 1;
			}
			else if (cur_feature->junc->exon_strand_start != terminal && cur_feature->junc->exon_strand_end == terminal)
			{	// transcription end exon of a gene
				exon_code = 2;
				num_cut_pt = 1;
			}
			else if ((cur_feature->junc->exon_strand_start == sense && cur_feature->junc->exon_strand_end == antisense)
				|| (cur_feature->junc->exon_strand_start == antisense && cur_feature->junc->exon_strand_end == sense))
			{	// mixed exon of two genes
				exon_code = 3;
				num_cut_pt = 2;
			}
		}

		if (exon_code > 0)
		{
			// collect coverage on this exon
			for (long pos = cur_feature->junc->start; pos <= cur_feature->junc->end; ++pos)
				exon_coverage.push_back(double(chromosome_expr[pos]));

			if (refine_exon_bound(exon_coverage, num_cut_pt, exon_code, cut_pt, avg_coverage))
			{
				if (exon_code == 1)
				{
#ifdef DEBUG_OUTPUT
					debug_out_cut_bound << cur_feature->junc->chromosome_start << "\t" << cur_feature->junc->start << "\t" << cur_feature->junc->end << "\t" << exon_code << "\t" << cur_feature->junc->start - 1 + cut_pt[0] << endl;
#endif
					cur_feature->junc->start = cur_feature->junc->start - 1 + cut_pt[0];
					cur_feature->junc->support[SUPPORT_VECTOR_SIZE_WHOLE_DATA] = avg_coverage[1];
				}
				else if (exon_code == 2)
				{
#ifdef DEBUG_OUTPUT
					debug_out_cut_bound << cur_feature->junc->chromosome_start << "\t" << cur_feature->junc->start << "\t" << cur_feature->junc->end << "\t" << exon_code << "\t" << cur_feature->junc->start - 1 + cut_pt[0] << endl;
#endif
					cur_feature->junc->end = cur_feature->junc->start - 1 + cut_pt[0];
					cur_feature->junc->support[SUPPORT_VECTOR_SIZE_WHOLE_DATA] = avg_coverage[0];
				}
				else if (exon_code == 3)
				{		
#ifdef DEBUG_OUTPUT
					debug_out_cut_bound << cur_feature->junc->chromosome_start << "\t" << cur_feature->junc->start << "\t" << cur_feature->junc->end << "\t" << exon_code << "\t" << cur_feature->junc->start - 1 + cut_pt[0] << "\t" << cur_feature->junc->start - 1 + cut_pt[1] << endl;
#endif
			
					fragment *new_fragment = new fragment(SUPPORT_VECTOR_SIZE);
					new_fragment->type = cur_feature->junc->type;
					new_fragment->frag_name = "ExonAdd_" + itostr(fragmentID_Cnt + 1);
					new_fragment->frag_id = ++fragmentID_Cnt;

					new_fragment->chromosome_start = cur_feature->junc->chromosome_start;
					new_fragment->chromosome_end = cur_feature->junc->chromosome_end;
					new_fragment->start = cur_feature->junc->start - 1 + cut_pt[1];
					new_fragment->end = cur_feature->junc->end;
					new_fragment->support[SUPPORT_VECTOR_SIZE_WHOLE_DATA] = avg_coverage[2];

					new_fragment->exon_strand_start = terminal;
					new_fragment->exon_strand_end = cur_feature->junc->exon_strand_end;
					
					cur_feature->junc->end = cur_feature->junc->start - 1 + cut_pt[0];
					cur_feature->junc->support[SUPPORT_VECTOR_SIZE_WHOLE_DATA] = avg_coverage[0];
					cur_feature->junc->exon_strand_end = terminal;

					rangeJunction *new_feature = new rangeJunction;
					new_feature->junc = new_fragment;
					new_feature->next = cur_feature->next;
					cur_feature->next = new_feature;
					cur_feature = new_feature;
				}
			}

			exon_coverage.clear();
		}

		cur_feature = cur_feature->next;
	}

#ifdef DEBUG_OUTPUT
	debug_out_cut_bound.close();
#endif

	return;
}
#endif 

//collect average read coverage on every exon
void collect_exon_expression_sequential(RangeJunctionList *cur_range, int index, string infilename_frag)
{
	ifstream infile_frag;
	infile_frag.open(infilename_frag.c_str());

	if (!infile_frag.is_open())
	{
		cout << "warning: fail to open fragment file " << infilename_frag << endl;
		return;
	}

	long position_start, position_end, occurrence, iLoop;
	rangeJunction *cur_feature = cur_range->list, *tmp_feature;
	string info;

	while (infile_frag >> occurrence)
	{
		infile_frag >> position_start;
		infile_frag >> position_end;
		getline(infile_frag, info);

		position_start += READ_COORDINATE_OFFSET;
		position_end += READ_COORDINATE_OFFSET;

		while (cur_feature != NULL && position_start > cur_feature->junc->end)
		{
			cur_feature = cur_feature->next;
		}

		if (cur_feature == NULL)
			break;

		tmp_feature = cur_feature;

		while (tmp_feature != NULL && position_start <= tmp_feature->junc->end && position_end >= tmp_feature->junc->start)
		{

			if (position_start >= tmp_feature->junc->start)
			{
				if (position_end <= tmp_feature->junc->end)
				{
					((tmp_feature->junc->support)[index]) += occurrence * (position_end - position_start + 1);
				}
				else if (position_end > tmp_feature->junc->end)
				{
					((tmp_feature->junc->support)[index]) += occurrence * (tmp_feature->junc->end - position_start + 1);
				}
			} 
			else
			{
				if (position_end <= tmp_feature->junc->end)
				{
					((tmp_feature->junc->support)[index]) += occurrence * (position_end - tmp_feature->junc->start + 1);
				}
				else if (position_end > tmp_feature->junc->end)
				{
					((tmp_feature->junc->support)[index]) += occurrence * (tmp_feature->junc->end - tmp_feature->junc->start + 1);
				}
			}
			tmp_feature = tmp_feature->next;
		}
	}

	infile_frag.close();

	cur_feature = cur_range->list;
	while (cur_feature != NULL)
	{
		assert(cur_feature->junc->type == frag_exon || cur_feature->junc->type == frag_retained_intron);
		cur_feature->junc->support[index] = double(cur_feature->junc->support[index]) / (cur_feature->junc->end - cur_feature->junc->start + 1);

		cur_feature = cur_feature->next;
	}

	return;
}


//read the input file, collect the expression on the chromosome
void collect_chromosome_expr_frag(string chromosome, string dir_frag, long &chr_start, long &chr_end)
{	
	string filename, info;

	//read in the chromosome start and end position
	filename = dir_frag + "/stat.txt";
	ifstream infile_stat(filename.c_str());
	if (!infile_stat.is_open())
	{
		cout << "error: fail to open file " << filename << endl;
		exit(1);
	}
	long sample_index, sample_startpos, sample_endpos;
	while (infile_stat >> sample_index)
	{
		infile_stat >> sample_startpos;
		infile_stat >> sample_endpos;
		getline(infile_stat, info);

		sample_startpos += READ_COORDINATE_OFFSET;
		sample_endpos += READ_COORDINATE_OFFSET;

		if (sample_startpos < chr_start)
			chr_start = sample_startpos;
		if (sample_endpos > chr_end)
			chr_end = sample_endpos;
	}

	infile_stat.close();

	//read in merged read expression and record on expression array
#ifdef UNIX
	filename = dir_frag + "/merged_exonic_sort.txt";
#else
	filename = dir_frag + "/" + itostr(debug_sample_index_offset+1) + "_exonic_sort.txt";
	//filename = dir_frag + "/merged_exonic_sort.txt";
#endif
	ifstream infile_frag(filename.c_str());

	if (infile_frag.is_open() == false)
	{
		cout << "error: fail to open input file " << filename << endl;
		exit(1);
	}

	long occurrence, startpos, endpos, tmploop;
	while (infile_frag >> occurrence)
	{
		infile_frag >> startpos;
		infile_frag >> endpos;
		getline(infile_frag, info);

		startpos += READ_COORDINATE_OFFSET;
		endpos += READ_COORDINATE_OFFSET;

		for (tmploop = startpos; tmploop <= endpos; ++tmploop)
			chromosome_expr[tmploop] += occurrence;
	}

	infile_frag.close();

// #ifdef UNIX
// #ifdef DEBUG_OUTPUT
// 	filename = "debug_output/chr_coverage_" + chromosome;
// 	ofstream outfile_chrcov(filename.c_str());
// 	for (tmploop = 1; tmploop <= chr_end; ++tmploop)
// 		outfile_chrcov << tmploop << "\t" << chromosome_expr[tmploop] << endl;
// 	outfile_chrcov.close();
// #endif
// #endif

// #ifndef UNIX
// 	filename = "tmp/tmp.txt";
// 	ifstream infile_chrcov(filename.c_str());
// 	long position;
// 	while (infile_chrcov >> position)
// 	{
// 		position += READ_COORDINATE_OFFSET;
// 		infile_chrcov >> chromosome_expr[position];
// 	}
// 	infile_chrcov.close();
// #endif

	return;
}


void drop_incomplete_junctions(RangeJunctionList *list_fragments)
{
	//filter standalone fragments (exons and junctions)
	if (list_fragments->cnt_rangeJunction <= 1)
	{
		//cout << "Warning: less than 2 fragments." << endl;
		return;
	}

#ifdef DEBUG_OUTPUT
	ofstream debug_out_deleted_junc("debug_output/debug_deleted_junctions.txt");
#endif

	rangeJunction *curJunc, *prevJunc, *tmpJunc;
	bool tobedeleted;

	void** sortlist_feature = new void* [list_fragments->cnt_rangeJunction + 2];
	double* sortkey_feature = new double [list_fragments->cnt_rangeJunction + 2];

	//first, both ends of a junction must appear. filter false ones
	//filter junctions with no start exon
		
	rangeJunction *tmp_feature = list_fragments->list;
	long num_features = 0, tmpLoop;
	while (tmp_feature != NULL)
	{
		sortlist_feature[++num_features] = (void*) tmp_feature;
		tmp_feature = tmp_feature->next;
	}
	for (tmpLoop = 1; tmpLoop <= num_features; ++tmpLoop)
		sortkey_feature[tmpLoop] = ((rangeJunction*)sortlist_feature[tmpLoop])->junc->end;
	mergeSort(sortlist_feature, sortkey_feature, num_features);

	for (tmpLoop = num_features; tmpLoop > 1; --tmpLoop)
		((rangeJunction*)sortlist_feature[tmpLoop])->next = (rangeJunction*) sortlist_feature[tmpLoop-1];
	((rangeJunction*)sortlist_feature[1])->next = NULL;
	list_fragments->list = ((rangeJunction*)sortlist_feature[num_features]);
	list_fragments->listtail = ((rangeJunction*)sortlist_feature[1]);
		

	curJunc = list_fragments->list;
	prevJunc = NULL;
	while (curJunc != NULL)
	{
		if (curJunc->junc->type == frag_junction)
		{
			tobedeleted = true;
			tmpJunc = curJunc->next;
			while (tmpJunc != NULL && tmpJunc->junc->end >= curJunc->junc->start)
			{
				if ((tmpJunc->junc->type == frag_exon || tmpJunc->junc->type == frag_retained_intron) && tmpJunc->junc->end == curJunc->junc->start)
				{
					tobedeleted = false;
					curJunc->junc->splice_exonhead = tmpJunc->junc;
					break;
				}
				tmpJunc = tmpJunc->next;
			}
			if (tobedeleted == true)
			{
#ifdef DEBUG_OUTPUT
				debug_out_deleted_junc << curJunc->junc->start << "\t" << curJunc->junc->end << endl;
#endif
				if (prevJunc != NULL)
				{
					prevJunc->next = curJunc->next;
					delete curJunc;
					--list_fragments->cnt_rangeJunction;
					curJunc = prevJunc->next;
				}
				else
				{
					list_fragments->list = curJunc->next;
					delete curJunc;
					--list_fragments->cnt_rangeJunction;
					curJunc = list_fragments->list;
				}
			}
			else
			{
				prevJunc = curJunc;
				curJunc = curJunc->next;
			}
		}
		else
		{
			prevJunc = curJunc;
			curJunc = curJunc->next;
		}
	}


	//filter junctions with no end exon

	tmp_feature = list_fragments->list;
	num_features = 0;
	while (tmp_feature != NULL)
	{
		sortlist_feature[++num_features] = (void*) tmp_feature;
		tmp_feature = tmp_feature->next;
	}
	for (tmpLoop = 1; tmpLoop <= num_features; ++tmpLoop)
		sortkey_feature[tmpLoop] = ((rangeJunction*)sortlist_feature[tmpLoop])->junc->start;
	mergeSort(sortlist_feature, sortkey_feature, num_features);

	for (tmpLoop = 1; tmpLoop < num_features; ++tmpLoop)
		((rangeJunction*)sortlist_feature[tmpLoop])->next = (rangeJunction*) sortlist_feature[tmpLoop+1];
	((rangeJunction*)sortlist_feature[num_features])->next = NULL;
	list_fragments->list = ((rangeJunction*)sortlist_feature[1]);
	list_fragments->listtail = ((rangeJunction*)sortlist_feature[num_features]);

	curJunc = list_fragments->list;
	prevJunc = NULL;
	while (curJunc != NULL)
	{
		if (curJunc->junc->type == frag_junction)
		{
			tobedeleted = true;
			tmpJunc = curJunc->next;
			while (tmpJunc != NULL && tmpJunc->junc->start <= curJunc->junc->end)
			{
				if ((tmpJunc->junc->type == frag_exon || tmpJunc->junc->type == frag_retained_intron) && tmpJunc->junc->start == curJunc->junc->end)
				{
					tobedeleted = false;
					curJunc->junc->splice_exontail = tmpJunc->junc;
					break;
				}
				tmpJunc = tmpJunc->next;
			}
			if (tobedeleted == true)
			{
#ifdef DEBUG_OUTPUT
				debug_out_deleted_junc << curJunc->junc->start << "\t" << curJunc->junc->end << endl;
#endif
				if (prevJunc != NULL)
				{
					prevJunc->next = curJunc->next;
					delete curJunc;
					--list_fragments->cnt_rangeJunction;
					curJunc = prevJunc->next;
				}
				else
				{
					list_fragments->list = curJunc->next;
					delete curJunc;
					--list_fragments->cnt_rangeJunction;
					curJunc = list_fragments->list;
				}
			}
			else
			{
				prevJunc = curJunc;
				curJunc = curJunc->next;
			}
		}
		else
		{
			prevJunc = curJunc;
			curJunc = curJunc->next;
		}
	}

	delete [] sortlist_feature;
	delete [] sortkey_feature;

#ifdef DEBUG_OUTPUT
	debug_out_deleted_junc.close();
#endif

	return;
}


void output_junction_bed_format(ofstream *out_file)
{
	fragment *curFrag;

	for (long iLoop = 0; iLoop < junction_array.size(); ++iLoop)
	{
		if (junction_array[iLoop]->type == frag_junction)
		{
			curFrag = junction_array[iLoop];
			(*out_file) << curFrag->chromosome_start << "\t" << curFrag->start - 19 << "\t" << curFrag->end + 20 << "\talljunction\t1000\t";
			if (curFrag->transDirection == undetermined || curFrag->transDirection == sense)
				(*out_file) << "+\t";
			else
				(*out_file) << "-\t";
			(*out_file) << curFrag->start - 19 << "\t" << curFrag->end + 20 << "\t0\t2\t" << "20,20,\t0," << curFrag->end - curFrag->start + 19 << endl;
		}
	}

	return;
}

//filter intron with relatively low expression, relative to its nearby exons
void filter_intron_with_relative_expression(RangeJunctionList *cur_range, int index)
{
	//notice that the support should be divided by the exon length before this is called! have corrected here but need to pay attention if using this in other package like assembly
	rangeJunction *prev_feature, *cur_feature, *del_feature;
	double cov_nearby_exon_max, cov_last_exon = -1;
	int nearby_cnt;

	prev_feature = NULL;
	cur_feature = cur_range->list;

	while (cur_feature != NULL)
	{
		if (cur_feature->junc->type == frag_exon && (cur_feature->junc->in_exon_annotation == annotated || cur_feature->junc->in_exon_annotation == na_exon_ref))
		{
			//do nothing
			cov_last_exon = cur_feature->junc->support[index];
			prev_feature = cur_feature;
			cur_feature = cur_feature->next;
			continue;
		}

		cov_nearby_exon_max = 0.0;
		nearby_cnt = 0;
		if (cov_last_exon > 0) //(prev_feature != NULL && abs(cur_feature->junc->start - prev_feature->junc->end) < 2)
		{
			cov_nearby_exon_max = cov_last_exon; // prev_feature->junc->support[index];
			++nearby_cnt;
		}
		rangeJunction *tmp_feature = cur_feature->next;
		while (tmp_feature && (tmp_feature->junc->type != frag_exon || (tmp_feature->junc->in_exon_annotation != annotated && tmp_feature->junc->in_exon_annotation != na_exon_ref)))
		{
			tmp_feature = tmp_feature->next;
		}
		if (tmp_feature)
		{
			cov_nearby_exon_max = tmp_feature->junc->support[index] > cov_nearby_exon_max ? tmp_feature->junc->support[index] : cov_nearby_exon_max;
			++nearby_cnt;
		}

		if (nearby_cnt < 1 || cur_feature->junc->support[index] < cov_nearby_exon_max/2.)
		{
			//filter this retained intron
			if (prev_feature == NULL)
			{
				cur_range->list = cur_feature->next;
				if (cur_feature == cur_range->listtail)
					cur_range->listtail = cur_range->list;
			}
			else
			{
				prev_feature->next = cur_feature->next;
				if (cur_feature == cur_range->listtail)
					cur_range->listtail = prev_feature;
			}

			del_feature = cur_feature;
			cur_feature = cur_feature->next;
			delete del_feature;
		}
		else
		{
			prev_feature = cur_feature;
			cur_feature = cur_feature->next;
		}
	}
	
	return;	
}


// label the exons for visualization: exons will be ordered linearly; exons with smaller coordinates will have smaller id; introns will be considered into counting, which means that two exons will have consecutive id's only when no intron separates them
// assume all fragments in this list are exons, and they are sorted by coordinate
void label_exon_visID(RangeJunctionList *list_exon)
{
	rangeJunction *cur_exon = list_exon->list;
	unsigned long id = 0;
	long last_bound = -1;

	while (cur_exon)
	{
		if (abs(last_bound - cur_exon->junc->start) <= 1)
		{
			cur_exon->junc->vis_exonid = id;
		}
		else
		{
			cur_exon->junc->vis_exonid = ++id;
		}

		++id;
		last_bound = cur_exon->junc->end;

		cur_exon = cur_exon->next;
	}

	return;
}


//read in reads, collect junctions and hits, record exons
void input_data(RangeJunctionList *list_fragments, string dir_input_frag, string dir_junction_annotation, long chromosome_length)
{
	cout << " ab initio reconstruction: " << flush;
	unsigned long tmpLoop, num_features;
	rangeJunction *new_rangeJunction;
	int sampleLoop;
	string filename;

	chromosome_expr.assign(chromosome_length, 0);

	//////////////////////////////////////
	//get all junctions
	filename = dir_input_frag + "/merged_junction_sort.txt";
	input_junction_merged(chromosomeName, filename);

	//////////////////////////////////////
	//get expression on the chromosome
	collect_chromosome_expr_frag(chromosomeName, dir_input_frag, list_fragments->rangeLow, list_fragments->rangeHigh);
	cout << " expr " << flush;

	///////////////////////////////
	//sort splice junctions 
	void** sortlist_junction = new void* [junction_array.size() + 2];
	double* sortkey_junction = new double [junction_array.size() + 2];

	for (tmpLoop = 0; tmpLoop < junction_array.size(); ++tmpLoop)
	{
		sortkey_junction[tmpLoop+1] = junction_array[tmpLoop]->end;
		sortlist_junction[tmpLoop+1] = (void*) junction_array[tmpLoop];
	}
	mergeSort(sortlist_junction, sortkey_junction, junction_array.size());
	
	for (tmpLoop = 0; tmpLoop < junction_array.size(); ++tmpLoop)
	{
		sortkey_junction[tmpLoop+1] = ((fragment*)sortlist_junction[tmpLoop+1])->start;
	}
	mergeSort(sortlist_junction, sortkey_junction, junction_array.size());

	for (tmpLoop = 0; tmpLoop < junction_array.size(); ++tmpLoop)
		junction_array[tmpLoop] = (fragment*) sortlist_junction[tmpLoop+1];

	delete [] sortlist_junction;
	delete [] sortkey_junction;

	cout << "# of raw junctions = " << junction_array.size() << " " << flush;
	output_junction_bed_format(&outfile_junction_all);

	//////////////////////////////////////
	//select junction according to annotation and presence
	if (dir_junction_annotation.compare("no_anno") == 0){
		filename = "";
	}
	else{
		filename = dir_junction_annotation + chromosomeName + ".txt";
	}
	label_select_junctions(filename);
	cout << "after filtering = " << junction_array.size() << " " << flush; 
	output_junction_bed_format(&outfile_junction_filtered);

	//////////////////////////////////////
	//extract exon
	extract_exons(list_fragments, chromosomeName);

	//////////////////////////////////////
	//trim exon boundaries
	trim_exons(list_fragments, SUPPORT_VECTOR_SIZE);	

#ifdef UNIX
	cut_exon_bound(list_fragments);
#endif	

	//////////////////////////////////////
	//5/15/2013 filter introns
	rangeJunction *cur_feature = list_fragments->list;
	while (cur_feature != NULL)
	{
		assert(cur_feature->junc->type == frag_exon || cur_feature->junc->type == frag_retained_intron);
		cur_feature->junc->support[SUPPORT_VECTOR_SIZE] = double(cur_feature->junc->support[SUPPORT_VECTOR_SIZE]) / (cur_feature->junc->end - cur_feature->junc->start + 1);

		cur_feature = cur_feature->next;
	}
	filter_intron_with_relative_expression(list_fragments, SUPPORT_VECTOR_SIZE);			


	//////////////////////////////////////
	//collect coverage of every exon
#ifdef UNIX
	for (sampleLoop = 0; sampleLoop < SUPPORT_VECTOR_SIZE; ++sampleLoop)
	{
		filename = dir_input_frag + "/" + itostr(sampleLoop+1) + "_exonic_sort.txt";
		collect_exon_expression_sequential(list_fragments, sampleLoop, filename);
	}
#else
	for (sampleLoop = 0; sampleLoop < SUPPORT_VECTOR_SIZE; ++sampleLoop)
	{
		filename = dir_input_frag + "/" + itostr(sampleLoop+debug_sample_index_offset+1) + "_exonic_sort.txt";
		collect_exon_expression_sequential(list_fragments, sampleLoop, filename);
	}
#endif

	//////////////////////////////////////
	//label visualization id for all exons
	label_exon_visID(list_fragments);

	//////////////////////////////////////
	//add splice junctions to the feature list
	for (tmpLoop = 0; tmpLoop < junction_array.size(); ++tmpLoop)
	{
		new_rangeJunction = new rangeJunction;
		new_rangeJunction->junc = junction_array[tmpLoop];
		list_fragments->insert_feature(new_rangeJunction);
		junction_array[tmpLoop] = NULL;
	}
	junction_array.clear();

	//////////////////////////////////////
	//sort features in the ranges
	list_fragments->count_featurelist(); //need to re-count because features are not inserted one by one when trimming exons

	drop_incomplete_junctions(list_fragments);
	list_fragments->count_featurelist(); //should be redundant

	if (list_fragments->cnt_rangeJunction <= 0)
	{
		cout << "Warning: no expression input found, will exit now" << endl;
		exit(1);
	}

	void** sortlist_feature = new void* [list_fragments->cnt_rangeJunction + 2];
	double* sortkey_feature = new double [list_fragments->cnt_rangeJunction + 2];

	rangeJunction *tmp_feature = list_fragments->list;
	num_features = 0;
	while (tmp_feature != NULL)
	{
		sortlist_feature[++num_features] = (void*) tmp_feature;
		tmp_feature = tmp_feature->next;
	}

	assert(num_features == list_fragments->cnt_rangeJunction);

	for (tmpLoop = 1; tmpLoop <= num_features; ++tmpLoop)
		sortkey_feature[tmpLoop] = ((rangeJunction*)sortlist_feature[tmpLoop])->junc->end;
	mergeSort(sortlist_feature, sortkey_feature, num_features);

	for (tmpLoop = 1; tmpLoop <= num_features; ++tmpLoop)
		sortkey_feature[tmpLoop] = ((rangeJunction*)sortlist_feature[tmpLoop])->junc->start;
	mergeSort(sortlist_feature, sortkey_feature, num_features);

	for (tmpLoop = 1; tmpLoop < num_features; ++tmpLoop)
		((rangeJunction*)sortlist_feature[tmpLoop])->next = (rangeJunction*) sortlist_feature[tmpLoop+1];
	((rangeJunction*)sortlist_feature[num_features])->next = NULL;
	list_fragments->list = ((rangeJunction*)sortlist_feature[1]);
	list_fragments->listtail = ((rangeJunction*)sortlist_feature[num_features]);

	delete [] sortlist_feature;
	delete [] sortkey_feature;

	free_vector(junction_array);
	free_vector(chromosome_expr);

	cout << "total features = " << list_fragments->cnt_rangeJunction << " " << flush;

	return;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// NOW BEGIN THE FUNCTIONS THAT DEAL WITH REFERENCE-ASSISTED RECONSTRUCTION 03/09/2015

// read annotated exons from the annotation file specified by filename, and return the ref exons
// return true if successfully read the file, return false otherwise
bool input_ref_exons(string filename, vector<fragment *> *ref_exons)
{
	assert(ref_exons); // ref_exons should be initialized in the caller
	(*ref_exons).empty(); // clear existing content (redundant operation, just to make sure)

	ifstream infile(filename.c_str());
	if (!infile.is_open()) {
		return false; // could not open the annotation file
	}

	// file format: start_pos end_pos strand
	long pos_start, pos_end;
	string chr, strand, info;
	while (infile >> chr) {
		infile >> pos_start;
		infile >> pos_end;
		infile >> strand;
		getline(infile, info);

		fragment *new_exon = new fragment(0); // for now, only need to store the overall expression, use support[0] to retrieve the value
		new_exon->chromosome_start = chr;
		new_exon->chromosome_end = chr;
		new_exon->start = pos_start;
		new_exon->end = pos_end;
		new_exon->transDirection = strand.compare("-") == 0 ? antisense : sense;
		new_exon->type = frag_exon;

		(*ref_exons).emplace_back(new_exon);
	}
	infile.close();
	return true;
}

// collect the expression of reference exons; the expression is the averaged read coverage of an exon, averaging over all bases of the exon and over all samples
void collect_expr_ref_exons(vector<fragment *> *ref_exons) {
	assert(ref_exons);

	for (auto it_frag = ref_exons->begin(); it_frag != ref_exons->end(); ++it_frag) {
		long expr = 0;
		for (long pos = (*it_frag)->start; pos <= (*it_frag)->end; ++pos) {
			expr += chromosome_expr[pos];
		}
		(*it_frag)->support[0] = double(expr) / double(((*it_frag)->end - (*it_frag)->start + 1) * SUPPORT_VECTOR_SIZE_WHOLE_DATA);
	}
}

// for each junction determine how its splice site matches to the annotated exons
// note: borrow fragment's fields 'splice_exonhead' and 'splice_exontail' to store the associated ref exon, because the ref exon's expression will be necessary in junction filtering; these fields will be put back to NULL after the junction filtering
void match_junctions_to_ref_exons(const vector<fragment *> &ref_exons) {
	// here will just use the brute-force O(n^2) approach, because both lists are sorted and most junctions are short (meaning that we can either find the match or terminate simply by examining the next exon or next few)
	auto it_exon = ref_exons.cbegin();
	auto it_junc = junction_array.begin();
	while (it_junc != junction_array.end()) {
		// move anchor
		while (it_exon != ref_exons.cend() && (*it_exon)->end < (*it_junc)->start) {
			++it_exon;
		}
		if (it_exon == ref_exons.cend()) {
			break;
		}

		// lower splice site
		auto it_exon_tmp = it_exon;
		(*it_junc)->bound_start = out_exon;
		while (it_exon_tmp != ref_exons.cend() && (*it_exon_tmp)->start <= (*it_junc)->start) {
			if ((*it_junc)->start == (*it_exon_tmp)->end) { // ((*it_junc)->start == (*it_exon_tmp)->start || (*it_junc)->start == (*it_exon_tmp)->end) { // 3/15/2015 need the lower splice site to be on the higher boundary of the ref exon
				(*it_junc)->bound_start = boundary;
				(*it_junc)->splice_exonhead = *it_exon_tmp;
				break;
			}
			else if ((*it_junc)->start > (*it_exon_tmp)->start && (*it_junc)->start < (*it_exon_tmp)->end) { // not >= here because junction out at the first base is unlikely
				(*it_junc)->bound_start = in_exon;
				(*it_junc)->splice_exonhead = *it_exon_tmp;
				// break;
			}
			++it_exon_tmp;
		}

		// higher splice site
		(*it_junc)->bound_end = out_exon;
		while (it_exon_tmp != ref_exons.cend() && (*it_exon_tmp)->start <= (*it_junc)->end) {
			if ((*it_junc)->end == (*it_exon_tmp)->start) { // ((*it_junc)->end == (*it_exon_tmp)->start || (*it_junc)->end == (*it_exon_tmp)->end) {
				(*it_junc)->bound_end = boundary;
				(*it_junc)->splice_exontail = *it_exon_tmp;
				break;
			}
			else if ((*it_junc)->end > (*it_exon_tmp)->start && (*it_junc)->end < (*it_exon_tmp)->end) { // not <= here because junction in at the last base is unlikely
				(*it_junc)->bound_end = in_exon;
				(*it_junc)->splice_exontail = *it_exon_tmp;
				// break;
			}
			++it_exon_tmp;
		}

		++it_junc;
	}

	// label the rest junctions as out of exon
	while (it_junc != junction_array.end()) {
		(*it_junc)->bound_start = out_exon;
		(*it_junc)->bound_end = out_exon;
		++it_junc;
	}
}

// determine whether a splice site meets the filtering criteria (see below, label_select_junctions_ref_exons)
bool pass_filter_splicesite_ref(const exon_bound_type type, const double ref_exon_expr, const double *junc_supp) {
	bool pass_filter = true;

	if (type == boundary || type == in_exon) {
		double junc_supp_avg = 0.;
		for (long sampleLoop = 0; sampleLoop < SUPPORT_VECTOR_SIZE_WHOLE_DATA; ++sampleLoop)
		{
			junc_supp_avg += junc_supp[sampleLoop];
		}
		junc_supp_avg /= SUPPORT_VECTOR_SIZE_WHOLE_DATA;

		if (type == boundary) {
			pass_filter = junc_supp_avg >= ref_exon_expr * 0.05;
		}
		else {
			pass_filter = junc_supp_avg >= ref_exon_expr * 0.1;
		}
	}
	else { // type == out_exon || type == unknown
		if (junctionfilter_mode == "group"){
			pass_filter = pass_filter_group(junc_supp);
		}
		else if (junctionfilter_mode == "nogroup"){
			pass_filter = pass_filter_no_group(junc_supp);
		}
		else if (junctionfilter_mode == "none"){
			// do nothing
		}
		else {
			cout << "Error: unrecognized junction filter mode" << endl; // should never happen, because unrecognized mode will be automatically set as 'none'
			exit(1);
		}
	}

	return pass_filter;
}

//with reference exons available, select junctions with the following criteria:
// 1. using annotation as a white list enabled and the junction in the given junction annotation or
// 2. pass the location-and-expression-based filter according to how a junction matches to annotated exons, specifically
// 2.1 if the splice site matches a known boundary and the junction expression is >= 5% of the exon expression, keep it
// 2.2 if the splice site falls inside a known exon and the junction expression is >= 10% of the exon expression, keep it
// 2.3 if the splice site falls outside any known exon, then the junction must fulfill the expression filter defined by label_select_junctions
// note: both splice sites of a junction will be tested and need to pass the filter in order to keep the junction
void label_select_junctions_ref_exons(string filename_annotation, const vector<fragment *> &ref_exons)
{
	//vector <bool> junction_in_annotation (junction_array.size(), false); //consider as in annotation as default to tolerate the case when the annotation is not available

	unsigned long junctionLoop = 0, cnt_junction_in_anno = 0;

	bool use_white_list = true;
	if (junction_annotation_as_white_list.compare("no") == 0){
		use_white_list = false;
	}

	if (!filename_annotation.empty()) {
		ifstream infile(filename_annotation.c_str());

		if (infile.is_open())
		{
			setting_junc_anno_provided = true;

			long start_pos, end_pos;
			string info, chrname;

			while (infile >> chrname)
			{
				infile >> start_pos;
				infile >> end_pos;
				getline(infile, info);

				while (junctionLoop < junction_array.size() && (junction_array[junctionLoop]->start < start_pos || (junction_array[junctionLoop]->start == start_pos && junction_array[junctionLoop]->end < end_pos)))
				{
					//junction_in_annotation[junctionLoop] = false;
					++junctionLoop;
				}

				if (junctionLoop >= junction_array.size())
					break;

				if (junction_array[junctionLoop]->start == start_pos && junction_array[junctionLoop]->end == end_pos)
				{
					//find a match in junction annotation
					//junction_in_annotation[junctionLoop] = true;
					junction_array[junctionLoop]->in_junc_annotation = true;
					++junctionLoop;
					++cnt_junction_in_anno;
				}
			}

			infile.close();
		}
		else
		{
			cout << "warning: junction annotation file has been specified but cannot be opened." << endl;
		}
	}

	match_junctions_to_ref_exons(ref_exons);

#ifdef FILTER_JUNCTION
	vector <fragment*> junction_array_filtered;

	for (junctionLoop = 0; junctionLoop < junction_array.size(); ++junctionLoop)
	{
		double ref_exon_expr_start = junction_array[junctionLoop]->splice_exonhead ? junction_array[junctionLoop]->splice_exonhead->support[0] : 0.;
		double ref_exon_expr_end = junction_array[junctionLoop]->splice_exontail ? junction_array[junctionLoop]->splice_exontail->support[0] : 0.;
		junction_array[junctionLoop]->splice_exonhead = NULL; // these two fields were 'borrowed' to temporarily store the associated ref exons, need to clear them
		junction_array[junctionLoop]->splice_exontail = NULL;

		bool pass_filter = pass_filter_splicesite_ref(junction_array[junctionLoop]->bound_start, ref_exon_expr_start, junction_array[junctionLoop]->support) &&
				pass_filter_splicesite_ref(junction_array[junctionLoop]->bound_end, ref_exon_expr_end, junction_array[junctionLoop]->support);

		//if (num_present_sample >= thresh_num_present_samples || junction_in_annotation[junctionLoop])
		if (pass_filter || (use_white_list && junction_array[junctionLoop]->in_junc_annotation))
		{
			//select this junction
			junction_array_filtered.push_back(junction_array[junctionLoop]);
#ifndef UNIX
			//only keep the samples under debug
			double *tmp_support = new double [SUPPORT_VECTOR_SIZE_WHOLE_DATA+1];
			for (int sampleLoop = 0; sampleLoop < SUPPORT_VECTOR_SIZE_WHOLE_DATA; ++sampleLoop)
			{
				tmp_support[sampleLoop] = junction_array[junctionLoop]->support[debug_sample_index_offset+sampleLoop];
			}
			tmp_support[SUPPORT_VECTOR_SIZE_WHOLE_DATA] = junction_array[junctionLoop]->support[debug_total_sample_num];
			delete [] junction_array[junctionLoop]->support;
			junction_array[junctionLoop]->support = tmp_support;
#endif
		}
		else
		{
			//discard this junction
			delete junction_array[junctionLoop];
		}
	}

	//save the selected junctions
	junction_array.swap(junction_array_filtered);
	junction_array_filtered.clear();
	free_vector(junction_array_filtered);
#endif

	return;
}

// extract exon boundaries from filtered junctions and annotated exons
void extract_exons_ref(string chromosome, const vector<fragment *> &ref_exons, RangeJunctionList *cur_range) 
{
	fragment *new_feature;
	rangeJunction *new_rangeJunction;
	if (junction_array.empty() && ref_exons.empty())
	{
		new_feature = new fragment(SUPPORT_VECTOR_SIZE);
		new_feature->type = frag_exon;
		new_feature->in_exon_annotation = novel;
		new_feature->altersite = 0;
		new_feature->frag_name = "Exon_1";
		new_feature->frag_id = ++fragmentID_Cnt;

		new_feature->start = cur_range->rangeLow;
		new_feature->end = cur_range->rangeHigh;
		new_feature->chromosome_start = chromosome;
		new_feature->chromosome_end = chromosome;

		//insert exonic fragments to list
		new_rangeJunction = new rangeJunction;
		new_rangeJunction->junc = new_feature;

		cur_range->insert_feature(new_rangeJunction);

		return;
	}

	// collect splice sites from ref exon boundaries 
	vector<spliceSite> splicesite_vec;
	unordered_map<long, splicesite_type> ht_refsite; // hash table for ref splice sites: location -> directionOut type, will output debug info if see same location sites with different directions
	for (auto it_exon = ref_exons.cbegin(); it_exon != ref_exons.cend(); ++it_exon) {
		// add in the start position
		if (ht_refsite.find((*it_exon)->start) == ht_refsite.end()) {
			ht_refsite[(*it_exon)->start] = splice_in;

			spliceSite new_site (chromosome, (*it_exon)->start, false, (*it_exon)->transDirection, true);
			splicesite_vec.push_back(new_site);
		}
		else if (ht_refsite[(*it_exon)->start] == splice_out) {
			ht_refsite[(*it_exon)->start] = splice_both;
			outfile_debuginfo << "in&out exon boundary at " << (*it_exon)->start << endl;

			spliceSite new_site (chromosome, (*it_exon)->start, false, (*it_exon)->transDirection, true);
			splicesite_vec.push_back(new_site);
		}
		// else: ht_refsite[(*it_exon)->start] == splice_in || ht_refsite[(*it_exon)->start] == splice_both, do nothing

		// add in the end position
		if (ht_refsite.find((*it_exon)->end) == ht_refsite.end()) {
			ht_refsite[(*it_exon)->end] = splice_out;

			spliceSite new_site (chromosome, (*it_exon)->end, true, (*it_exon)->transDirection, true);
			splicesite_vec.push_back(new_site);
		}
		else if (ht_refsite[(*it_exon)->end] == splice_in) {
			ht_refsite[(*it_exon)->end] = splice_both;
			outfile_debuginfo << "in&out exon boundary at " << (*it_exon)->end << endl;

			spliceSite new_site (chromosome, (*it_exon)->end, true, (*it_exon)->transDirection, true);
			splicesite_vec.push_back(new_site);
		}
		// else: ht_refsite[(*it_exon)->end] == splice_out || ht_refsite[(*it_exon)->end] == splice_both, do nothing
	}

	//input splice sites
	for (auto it_junc = junction_array.cbegin(); it_junc != junction_array.cend(); ++it_junc)
	{
		// add in the start position
		if (ht_refsite.find((*it_junc)->start) == ht_refsite.end()) {
			ht_refsite[(*it_junc)->start] = splice_out;

			spliceSite new_site (chromosome, (*it_junc)->start, true, (*it_junc)->transDirection, false); // (*it_junc)->in_junc_annotation); decided to only trust the ref exons
			splicesite_vec.push_back(new_site);
		}
		else if (ht_refsite[(*it_junc)->start] == splice_in) {
			ht_refsite[(*it_junc)->start] = splice_both;
			outfile_debuginfo << "splice out on exon start at " << (*it_junc)->start << endl;

			spliceSite new_site (chromosome, (*it_junc)->start, true, (*it_junc)->transDirection, false); // (*it_junc)->in_junc_annotation);
			splicesite_vec.push_back(new_site);
		}
		// else: ht_refsite[(*it_junc)->start] == splice_out || ht_refsite[(*it_junc)->start] == splice_both, do nothing

		// add in the end position
		if (ht_refsite.find((*it_junc)->end) == ht_refsite.end()) {
			ht_refsite[(*it_junc)->end] = splice_in;

			spliceSite new_site (chromosome, (*it_junc)->end, false, (*it_junc)->transDirection, false); // (*it_junc)->in_junc_annotation);
			splicesite_vec.push_back(new_site);
		}
		else if (ht_refsite[(*it_junc)->end] == splice_out) {
			ht_refsite[(*it_junc)->end] = splice_both;
			outfile_debuginfo << "splice in on exon end at " << (*it_junc)->end << endl;

			spliceSite new_site (chromosome, (*it_junc)->end, false, (*it_junc)->transDirection, false); // (*it_junc)->in_junc_annotation);
			splicesite_vec.push_back(new_site);
		}
		// else: ht_refsite[(*it_junc)->end] == splice_in || ht_refsite[(*it_junc)->end] == splice_both, do nothing
	}

	assert(splicesite_vec.size() >= 2);
	sort(splicesite_vec.begin(), splicesite_vec.end(), comp_ssite_withdup);

	///////////////////////////////
	//build exonic fragments from adjacent splice sites 
	//exonic if not out-in

	long exon_count = 0;
	if (splicesite_vec.front().position > cur_range->rangeLow) { // coverage before the first annotated exon
		new_feature = new fragment(SUPPORT_VECTOR_SIZE);
		new_feature->type = frag_exon;
		new_feature->in_exon_annotation = splicesite_vec.front().is_ref ? extension : novel; // here simplify the conditions a lot: if the first site is not from a ref exon (hence a junction, no matter annotated or not), make this a novel exon (novel as not in the exon annotation); otherwise, the first site is a ref exon's start, so make this an extension
		new_feature->altersite = 0;
		new_feature->frag_name = "Exon_" + itostr(++exon_count);
		new_feature->frag_id = ++fragmentID_Cnt;

		new_feature->chromosome_start = splicesite_vec.front().chromosome;
		new_feature->chromosome_end = splicesite_vec.front().chromosome;
		new_feature->start = cur_range->rangeLow;
		new_feature->end = splicesite_vec.front().position;
		new_feature->exon_strand_start = terminal;
		new_feature->exon_strand_end = splicesite_vec.front().strand;

		//insert exonic fragments to list
		new_rangeJunction = new rangeJunction;
		new_rangeJunction->junc = new_feature;

		cur_range->insert_feature(new_rangeJunction);
	}

	for (vector<spliceSite>::iterator ssite_curr = splicesite_vec.begin(), ssite_next = splicesite_vec.begin() + 1; ssite_next != splicesite_vec.end(); ++ssite_curr, ++ssite_next)
	{
		assert(ssite_curr->position != ssite_next->position || ssite_curr->directionOut != ssite_next->directionOut);

		new_feature = new fragment(SUPPORT_VECTOR_SIZE);
		new_feature->type = frag_exon;

		if (ssite_curr->directionOut == false && ssite_next->directionOut == true)
			new_feature->altersite = 0;
		else if (ssite_curr->directionOut == true && ssite_next->directionOut == true)
			new_feature->altersite = 1;
		else if (ssite_curr->directionOut == false && ssite_next->directionOut == false)
			new_feature->altersite = 2;
		else if (ssite_curr->directionOut == true && ssite_next->directionOut == false)
		{
			new_feature->type = frag_retained_intron; // simplified here: can actually be intron or intergenic region. can be distinguished by adding a state to splicesite which indicates whether the site has a junction hit. but decided to leave as it is for now, to avoid additional work in trim_exons and other functions
			new_feature->altersite = -1;
		}

		if (new_feature->type == frag_exon) {
			if (ssite_curr->is_ref && ssite_next->is_ref)
				new_feature->in_exon_annotation = annotated;
			else if (ssite_curr->is_ref || ssite_next->is_ref)
				new_feature->in_exon_annotation = extension;
			else
				new_feature->in_exon_annotation = novel;
		}
		
		new_feature->frag_name = "Exon_" + itostr(++exon_count);
		new_feature->frag_id = ++fragmentID_Cnt;

		new_feature->chromosome_start = ssite_curr->chromosome;
		if (ssite_curr->directionOut)
			new_feature->start = ssite_curr->position + 1;
		else
			new_feature->start = ssite_curr->position;
		new_feature->exon_strand_start = ssite_curr->strand;

		new_feature->chromosome_end = ssite_next->chromosome;
		if (ssite_next->directionOut)
			new_feature->end = ssite_next->position;
		else
			new_feature->end = ssite_next->position - 1;
		new_feature->exon_strand_end = ssite_next->strand;

		if (new_feature->end >= new_feature->start)
		{
			//insert exonic fragments to list
			new_rangeJunction = new rangeJunction;
			new_rangeJunction->junc = new_feature;
			cur_range->insert_feature(new_rangeJunction);
		}
		else
		{
			delete new_feature;
		}
	}

	//last exon
	if (splicesite_vec.back().position < cur_range->rangeHigh) { // coverage before the first annotated exon
		new_feature = new fragment(SUPPORT_VECTOR_SIZE);
		new_feature->type = frag_exon;
		new_feature->in_exon_annotation = splicesite_vec.back().is_ref ? extension : novel;
		new_feature->altersite = 0;
		new_feature->frag_name = "Exon_" + itostr(++exon_count);
		new_feature->frag_id = ++fragmentID_Cnt;

		new_feature->chromosome_start = splicesite_vec.back().chromosome;
		new_feature->chromosome_end = splicesite_vec.back().chromosome;
		new_feature->start = splicesite_vec.back().position;
		new_feature->end = cur_range->rangeHigh;
		new_feature->exon_strand_start = splicesite_vec.back().strand;
		new_feature->exon_strand_end = terminal;

		//insert exonic fragments to list
		new_rangeJunction = new rangeJunction;
		new_rangeJunction->junc = new_feature;
		cur_range->insert_feature(new_rangeJunction);
	}

	return;
}


//read in reads, collect junctions and hits, record exons, use annotated exon with necessary filtering & extension
bool input_data_ref(RangeJunctionList *list_fragments, string dir_input_frag, string dir_junction_annotation, string dir_exon_annotation, long chromosome_length)
{
	unsigned long tmpLoop, num_features;
	rangeJunction *new_rangeJunction;
	int sampleLoop;
	string filename;

	//////////////////////////////////////
	//try to get annotated exons
	if (dir_exon_annotation.compare("no_anno") == 0){
		return false; //no exon annotation provided, should go to reconstruction module
	}
	else{
		filename = dir_exon_annotation + chromosomeName + ".ref";
	}
	vector<fragment *> ref_exons;
	if (!input_ref_exons(filename, &ref_exons)) {
		return false; // fail to read in annotated exons, will go to reconstruction module
	}
	sort(ref_exons.begin(), ref_exons.end(), comp_frag_ptr);
	cout << " ref-assisted reconstruction(" << ref_exons.size() << " ref exons received): " << flush;


	chromosome_expr.assign(chromosome_length, 0);

	//////////////////////////////////////
	//get all junctions
	filename = dir_input_frag + "/merged_junction_sort.txt";
	input_junction_merged(chromosomeName, filename);

	//////////////////////////////////////
	//get expression on the chromosome
	collect_chromosome_expr_frag(chromosomeName, dir_input_frag, list_fragments->rangeLow, list_fragments->rangeHigh);
	cout << " expr " << flush;

	//////////////////////////////////////
	//get expression of the annotated exons
	collect_expr_ref_exons(&ref_exons);


	///////////////////////////////
	//sort splice junctions
	void** sortlist_junction = new void* [junction_array.size() + 2];
	double* sortkey_junction = new double [junction_array.size() + 2];

	for (tmpLoop = 0; tmpLoop < junction_array.size(); ++tmpLoop)
	{
		sortkey_junction[tmpLoop+1] = junction_array[tmpLoop]->end;
		sortlist_junction[tmpLoop+1] = (void*) junction_array[tmpLoop];
	}
	mergeSort(sortlist_junction, sortkey_junction, junction_array.size());

	for (tmpLoop = 0; tmpLoop < junction_array.size(); ++tmpLoop)
	{
		sortkey_junction[tmpLoop+1] = ((fragment*)sortlist_junction[tmpLoop+1])->start;
	}
	mergeSort(sortlist_junction, sortkey_junction, junction_array.size());

	for (tmpLoop = 0; tmpLoop < junction_array.size(); ++tmpLoop)
		junction_array[tmpLoop] = (fragment*) sortlist_junction[tmpLoop+1];

	delete [] sortlist_junction;
	delete [] sortkey_junction;

	cout << "# of raw junctions = " << junction_array.size() << " " << flush;
	output_junction_bed_format(&outfile_junction_all);

	//////////////////////////////////////
	//select junction according to annotation and presence
	if (dir_junction_annotation.compare("no_anno") == 0){
		filename = "";
	}
	else{
		filename = dir_junction_annotation + chromosomeName + ".ref";
	}
	label_select_junctions_ref_exons(filename, ref_exons); 
	cout << "after filtering = " << junction_array.size() << " " << flush;
	output_junction_bed_format(&outfile_junction_filtered);

	//////////////////////////////////////
	//extract exon
	extract_exons_ref(chromosomeName, ref_exons, list_fragments);

	//////////////////////////////////////
	//trim exon boundaries // till here
	trim_exons(list_fragments, SUPPORT_VECTOR_SIZE);

// #ifdef UNIX
// 	cut_exon_bound(list_fragments);
// #endif

	//////////////////////////////////////
	//5/15/2013 filter introns
	rangeJunction *cur_feature = list_fragments->list;
	while (cur_feature != NULL)
	{
		assert(cur_feature->junc->type == frag_exon || cur_feature->junc->type == frag_retained_intron);
		cur_feature->junc->support[SUPPORT_VECTOR_SIZE] = double(cur_feature->junc->support[SUPPORT_VECTOR_SIZE]) / (cur_feature->junc->end - cur_feature->junc->start + 1);

		cur_feature = cur_feature->next;
	}
	filter_intron_with_relative_expression(list_fragments, SUPPORT_VECTOR_SIZE);


	//////////////////////////////////////
	//collect coverage of every exon
#ifdef UNIX
	for (sampleLoop = 0; sampleLoop < SUPPORT_VECTOR_SIZE; ++sampleLoop)
	{
		filename = dir_input_frag + "/" + itostr(sampleLoop+1) + "_exonic_sort.txt";
		collect_exon_expression_sequential(list_fragments, sampleLoop, filename);
	}
#else
	for (sampleLoop = 0; sampleLoop < SUPPORT_VECTOR_SIZE; ++sampleLoop)
	{
		filename = dir_input_frag + "/" + itostr(sampleLoop+debug_sample_index_offset+1) + "_exonic_sort.txt";
		collect_exon_expression_sequential(list_fragments, sampleLoop, filename);
	}
#endif

	//////////////////////////////////////
	//label visualization id for all exons
	label_exon_visID(list_fragments);

	//////////////////////////////////////
	//add splice junctions to the feature list
	for (tmpLoop = 0; tmpLoop < junction_array.size(); ++tmpLoop)
	{
		new_rangeJunction = new rangeJunction;
		new_rangeJunction->junc = junction_array[tmpLoop];
		list_fragments->insert_feature(new_rangeJunction);
		junction_array[tmpLoop] = NULL;
	}
	junction_array.clear();

	//////////////////////////////////////
	//sort features in the ranges
	list_fragments->count_featurelist(); //need to re-count because features are not inserted one by one when trimming exons

	drop_incomplete_junctions(list_fragments);
	list_fragments->count_featurelist(); //should be redundant

	if (list_fragments->cnt_rangeJunction <= 0)
	{
		cout << "Warning: no expression input found, will exit now" << endl;
		exit(1);
	}

	void** sortlist_feature = new void* [list_fragments->cnt_rangeJunction + 2];
	double* sortkey_feature = new double [list_fragments->cnt_rangeJunction + 2];

	rangeJunction *tmp_feature = list_fragments->list;
	num_features = 0;
	while (tmp_feature != NULL)
	{
		sortlist_feature[++num_features] = (void*) tmp_feature;
		tmp_feature = tmp_feature->next;
	}

	assert(num_features == list_fragments->cnt_rangeJunction);

	for (tmpLoop = 1; tmpLoop <= num_features; ++tmpLoop)
		sortkey_feature[tmpLoop] = ((rangeJunction*)sortlist_feature[tmpLoop])->junc->end;
	mergeSort(sortlist_feature, sortkey_feature, num_features);

	for (tmpLoop = 1; tmpLoop <= num_features; ++tmpLoop)
		sortkey_feature[tmpLoop] = ((rangeJunction*)sortlist_feature[tmpLoop])->junc->start;
	mergeSort(sortlist_feature, sortkey_feature, num_features);

	for (tmpLoop = 1; tmpLoop < num_features; ++tmpLoop)
		((rangeJunction*)sortlist_feature[tmpLoop])->next = (rangeJunction*) sortlist_feature[tmpLoop+1];
	((rangeJunction*)sortlist_feature[num_features])->next = NULL;
	list_fragments->list = ((rangeJunction*)sortlist_feature[1]);
	list_fragments->listtail = ((rangeJunction*)sortlist_feature[num_features]);

	delete [] sortlist_feature;
	delete [] sortkey_feature;

	free_vector(junction_array);
	free_vector(chromosome_expr);

	for (auto it_frag = ref_exons.begin(); it_frag != ref_exons.end(); ++it_frag) {
		delete *it_frag;
	}

	cout << "total features = " << list_fragments->cnt_rangeJunction << " " << flush;

	return true;
}
