#include "write_output.h"

extern double significance_cutoff;



string ntostr(double t)
{
	ostringstream oss;
	oss << t;
	return oss.str();
}


void assign_gene_for_asm()
{
	unsigned chrLoop, geneLoop;
	annotation_gene *curGene;
	ASM *cur_asm;

	for (unsigned long asm_loop = 1; asm_loop < sortList_ASM.size(); ++asm_loop)
	{
		cur_asm = sortList_ASM[asm_loop];

		for (chrLoop = 0; chrLoop < list_anno_chr_genelist.size(); ++chrLoop)
			if (list_anno_chr_genelist[chrLoop]->chrname.compare(cur_asm->chrNM) == 0)
				break;
		if (chrLoop >= list_anno_chr_genelist.size())
		{
			continue;
		}

		for (geneLoop = 0; geneLoop < list_anno_chr_genelist[chrLoop]->genelist.size(); ++geneLoop)
		{
			curGene = list_anno_chr_genelist[chrLoop]->genelist[geneLoop];

			if (cur_asm->rangeHigh < curGene->range_start)
				break;

			// 		if (curGene->geneNm.compare("MAST2|23139") == 0)
			// 			cout << "1";

			if (cur_asm->rangeLow < curGene->range_end && cur_asm->rangeHigh > curGene->range_start)
			{
				// ASM overlaps with the gene
				long overlap_pos_end = cur_asm->rangeHigh < curGene->range_end ? cur_asm->rangeHigh : curGene->range_end;
				long overlap_pos_start = cur_asm->rangeLow > curGene->range_start ? cur_asm->rangeLow : curGene->range_start;
				long overlap_length = overlap_pos_end - overlap_pos_start + 1;

				if (overlap_length > (cur_asm->rangeHigh-cur_asm->rangeLow+1) * 0.5 || overlap_length > (curGene->range_end-curGene->range_start+1) * 0.5)
				{
					//in this gene
					if (cur_asm->in_gene.compare("na") == 0)
						cur_asm->in_gene = curGene->geneNm;
					else
						cur_asm->in_gene += "," + curGene->geneNm;
				}
			}
		}
	}

	return;
}

void calculate_max_pairwise_JSD_betwgroup()
{
	for (unsigned long asm_loop = 1; asm_loop < sortList_ASM.size(); ++asm_loop)
	{
		double cur_max_JSD = 0;
		string cur_max_JSD_groupname = "", cur_large_JSD_info = "";

		for (int group1 = 1; group1 < global_total_group_num; ++group1)
		{
			for (int group2 = group1+1; group2 <= global_total_group_num; ++group2)
			{
				double cur_JSD = calculate_JSD(sortList_ASM[asm_loop]->group_mean_path_prop[group1], sortList_ASM[asm_loop]->group_mean_path_prop[group2], sortList_ASM[asm_loop]->pathNum);
				if (cur_JSD > cur_max_JSD)
				{
					cur_max_JSD = cur_JSD;
					cur_max_JSD_groupname = global_group_name[group1] + "-" + global_group_name[group2];
				}

				cur_JSD = sqrt(cur_JSD); 
				if (cur_JSD > thresh_JSD)
				{
					cur_large_JSD_info += global_group_name[group1] + "-" + global_group_name[group2] + "," + ntostr(cur_JSD) + "; ";
				}
			}
		}

		sortList_ASM[asm_loop]->max_pairwise_JSD_betwgroup = cur_max_JSD;
		sortList_ASM[asm_loop]->max_pairwise_JSD_betwgroup_groupname = cur_max_JSD_groupname;
		sortList_ASM[asm_loop]->info_large_pairwise_JSD_betwgroup = cur_large_JSD_info;
	}

	return;
}


void calculate_minmax_dist_to_groupmean_ratio()
{
	for (unsigned long asm_loop = 1; asm_loop < sortList_ASM.size(); ++asm_loop)
	{
		double cur_max_mindistratio = 0;
		string cur_max_mindistratio_groupname = "", cur_large_mindistratio_info = "";

		for (int group1 = 1; group1 < global_total_group_num; ++group1)
		{
			for (int group2 = group1+1; group2 <= global_total_group_num; ++group2)
			{
				double jsd1, jsd2, cur_ratio, cur_minratio = MAX_NUMBER;
				for (int sample_loop = global_group_base_index[group1]+1; sample_loop <= global_group_base_index[group1]+global_group_size[group1]; ++sample_loop)
				{
					jsd1 = calculate_JSD(sortList_ASM[asm_loop]->group_mean_path_prop[group1], sortList_ASM[asm_loop]->sample_path_prop[sample_loop], sortList_ASM[asm_loop]->pathNum);
					jsd2 = calculate_JSD(sortList_ASM[asm_loop]->group_mean_path_prop[group2], sortList_ASM[asm_loop]->sample_path_prop[sample_loop], sortList_ASM[asm_loop]->pathNum);
					cur_ratio = sqrt(jsd2) / sqrt(jsd1);
					if (cur_ratio < cur_minratio)
						cur_minratio = cur_ratio;
				}
				for (int sample_loop = global_group_base_index[group2]+1; sample_loop <= global_group_base_index[group2]+global_group_size[group2]; ++sample_loop)
				{
					jsd1 = calculate_JSD(sortList_ASM[asm_loop]->group_mean_path_prop[group1], sortList_ASM[asm_loop]->sample_path_prop[sample_loop], sortList_ASM[asm_loop]->pathNum);
					jsd2 = calculate_JSD(sortList_ASM[asm_loop]->group_mean_path_prop[group2], sortList_ASM[asm_loop]->sample_path_prop[sample_loop], sortList_ASM[asm_loop]->pathNum);
					cur_ratio = sqrt(jsd1) / sqrt(jsd2);
					if (cur_ratio < cur_minratio)
						cur_minratio = cur_ratio;
				}

				if (cur_minratio < MAX_NUMBER && cur_minratio > cur_max_mindistratio)
				{
					cur_max_mindistratio = cur_minratio;
					cur_max_mindistratio_groupname = global_group_name[group1] + "-" + global_group_name[group2];
				}

				if (cur_minratio < MAX_NUMBER && cur_minratio > 0.5) //should be >1 to make sense actually, but small values here suggest no perfect group separation
				{
					cur_large_mindistratio_info += global_group_name[group1] + "-" + global_group_name[group2] + "," + ntostr(cur_minratio) + "; ";
				}
			}
		}

		sortList_ASM[asm_loop]->max_grouppair_min_dist_ratio = cur_max_mindistratio;
		sortList_ASM[asm_loop]->max_grouppair_min_dist_ratio_groupname = cur_max_mindistratio_groupname;
		sortList_ASM[asm_loop]->info_large_grouppair_min_dist_ratio = cur_large_mindistratio_info;
	}

	return;
}



void calculate_max_percent_large_dist()
{
	const double thresh_large_distance = 0.3; //sqrt JSD
	const double thresh_large_distance_ratio = 0.90; //label large distance if this sample has distance large than *thresh_large_distance* with at least *ratio* of the samples in the other group

	for (unsigned long asm_loop = 1; asm_loop < sortList_ASM.size(); ++asm_loop)
	{
		double cur_max_percent = 0;
		string cur_max_percent_groupname = "", cur_large_percent_info = "";

		for (int group1 = 1; group1 <= global_total_group_num; ++group1)
		{
			for (int group2 = 1; group2 <= global_total_group_num; ++group2)
			{
				if (group1 == group2)
					continue;

				int cnt_large_dist = 0;
				for (int sample_loop1 = global_group_base_index[group1]+1; sample_loop1 <= global_group_base_index[group1]+global_group_size[group1]; ++sample_loop1)
				{
					int tmp_cnt_large_dist = 0;
					
					for (int sample_loop2 = global_group_base_index[group2]+1; sample_loop2 <= global_group_base_index[group2]+global_group_size[group2]; ++sample_loop2)
					{
						double jsd = calculate_JSD(sortList_ASM[asm_loop]->sample_path_prop[sample_loop1], sortList_ASM[asm_loop]->sample_path_prop[sample_loop2], sortList_ASM[asm_loop]->pathNum);
						if (sqrt(jsd) >= thresh_large_distance)
						{
							++tmp_cnt_large_dist;
						}
					}

					if (double(tmp_cnt_large_dist) / global_group_size[group2] >= thresh_large_distance_ratio)
					{
						++cnt_large_dist;
					}
				}

				double cur_percent = double(cnt_large_dist) / global_group_size[group1];
				if (cur_percent > cur_max_percent)
				{
					cur_max_percent = cur_percent;
					cur_max_percent_groupname = global_group_name[group1] + "-" + global_group_name[group2];
				}

				if (cur_percent > 0.3)
				{
					cur_large_percent_info += global_group_name[group1] + "-" + global_group_name[group2] + "," + ntostr(cur_percent) + "; ";
				}
			}
		}

		sortList_ASM[asm_loop]->max_percent_large_dist = cur_max_percent;
		sortList_ASM[asm_loop]->max_percent_large_dist_groupname = cur_max_percent_groupname;
		sortList_ASM[asm_loop]->info_large_percent_large_dist = cur_large_percent_info;
	}

	return;
}



void output_difftrans_table(string resultPath)
{
	string filename;
	unsigned long asm_loop;
	//bool sameDirection, sameDirectionVSexpected;

	/************************************************************************/
	/* Sort All Genes                                                       */
	/************************************************************************/
	for (asm_loop = 1; asm_loop < sortList_ASM.size(); ++asm_loop)
	{
		sortKey_ASM[asm_loop] = sortList_ASM[asm_loop]->statistic_d * (-1.0);
	}
	mergeSort_ASM_sort(sortList_ASM.size()-1);


	/************************************************************************/
	/* Diagnostics                                                          */
	/************************************************************************/
	//	sprintf(filename, "%sstat/compare_d_asm.txt", resultPath);
	//	ofstream d_statistics_outputfile(filename);
	//	ASM *curASM;
	// 	for (geneLoopCnt = 1; geneLoopCnt <= sortList_ASM_Num; ++geneLoopCnt)
	// 	{
	// 		curASM = sortList_ASM[geneLoopCnt];
	// 		sameDirection = true;
	// 		sameDirectionVSexpected = true;
	// 
	// 		d_statistics_outputfile << curASM->ASM_id << "\t" << curASM->chrNM << "\t" << curASM->rangeLow << "\t" <<curASM->rangeHigh << "\t" << curASM->ASMcategory << "\t" << curASM->pathNum << "\t"
	// 			<< curASM->statistic_d << '\t' << curASM->statistic_d_expected << "\t" << fabs(curASM->statistic_d - curASM->statistic_d_expected) << "\t" << curASM->statistic_s << "\t"
	// 			<< sqrt(calculate_JSD(curASM->meanPathProportion[1], curASM->meanPathProportion[2], curASM->pathNum)) << "\t";
	// 
	// 		for (individualLoopCnt = 1; individualLoopCnt <= num_blocks; ++individualLoopCnt)
	// 		{
	// 			d_statistics_outputfile << curASM->statistic_d_individual[individualLoopCnt] << "\t" << curASM->statistic_d_expected_individual[individualLoopCnt] << "\t"
	// 				<< fabs(curASM->statistic_d_individual[individualLoopCnt] - curASM->statistic_d_expected_individual[individualLoopCnt]) << "\t" << curASM->statistic_s_individual[individualLoopCnt] << "\t";
	// 	
	// 			if (curASM->statistic_d_individual[individualLoopCnt] * curASM->statistic_d < 0)
	// 			{
	// 				sameDirection = false;
	// 			}
	// 			if ((curASM->statistic_d_individual[individualLoopCnt] - curASM->statistic_d_expected_individual[individualLoopCnt]) * (curASM->statistic_d - curASM->statistic_d_expected) < 0)
	// 			{
	// 				sameDirectionVSexpected = false;
	// 			}
	// 		}
	// 
	// 		d_statistics_outputfile << sameDirection << "\t" << sameDirectionVSexpected << "\t" << curASM->meanExpression[1] << "\t" << curASM->meanExpression[2] << endl;
	// 	}
	// 	d_statistics_outputfile.close();
	filename = resultPath + "differential_transcription.txt";
	ofstream d_statistics_outputfile(filename.c_str());
	
	long **pair_diff_cnt = new long* [global_total_group_num+1];
	for (int group_loop = 1; group_loop <= global_total_group_num; ++group_loop)
	{
		pair_diff_cnt[group_loop] = new long [global_total_group_num+1];
		for (int group_loop_2 = 1; group_loop_2 <= global_total_group_num; ++group_loop_2)
			pair_diff_cnt[group_loop][group_loop_2] = 0;
	}

	ASM *curASM;
	long significantCnt = 0;
	double stat_difftrans, jsd;
	d_statistics_outputfile << "location\tasm_id\tgene\tnum of paths\tcategory\tstat_diff_trans(|stat_d_expected-stat_d|)\tmax pairwise sqrtJSD\tmax pairwise sqrtJSD group names\tgrand_mean_cov\tgroup_mean_cov\tsignificant\tlarge pairwise sqrtJSD groups\tsignificant changes between two groups" <<endl;
	for (asm_loop = 1; asm_loop < sortList_ASM.size(); ++asm_loop)
	{
		curASM = sortList_ASM[asm_loop];
		stat_difftrans = fabs(curASM->statistic_d - curASM->statistic_d_expected);
		jsd = sqrt(curASM->max_pairwise_JSD_betwgroup);

		d_statistics_outputfile << curASM->chrNM << " " << curASM->rangeLow << " " << curASM->rangeHigh << "\t" << curASM->ASM_id << "\t" << curASM->in_gene << "\t" << curASM->pathNum << '\t' << alterSpliceCategory[curASM->ASMcategory] << "\t"
			<< stat_difftrans << "\t" << jsd << "\t" << curASM->max_pairwise_JSD_betwgroup_groupname << "\t" << curASM->grand_mean_expr << "\t";
		for (int group_loop = 1; group_loop <= global_total_group_num; ++group_loop)
			d_statistics_outputfile << curASM->group_mean_expr[group_loop] << ", ";
		d_statistics_outputfile << "\t";
		if (stat_difftrans >= significance_cutoff && jsd >= thresh_JSD)// && stat_difftrans >= thresh_stat_d)
		{
			curASM->significant_both = true;
			d_statistics_outputfile << "yes\t";
			++significantCnt;
		}
		else
			d_statistics_outputfile << "no\t";
		d_statistics_outputfile << curASM->info_large_pairwise_JSD_betwgroup;

		d_statistics_outputfile << "\t";
		for (unsigned long pair_loop = 0; pair_loop < curASM->list_pairdiff.size(); ++pair_loop)
		{
			if (curASM->list_pairdiff[pair_loop]->significant_both)
			{
				d_statistics_outputfile << global_group_name[curASM->list_pairdiff[pair_loop]->groupindex_1] << "-" << global_group_name[curASM->list_pairdiff[pair_loop]->groupindex_2] << " ("
					<< curASM->list_pairdiff[pair_loop]->stat_difftrans << "," << sqrt(curASM->list_pairdiff[pair_loop]->jsd) << "); ";

				pair_diff_cnt[curASM->list_pairdiff[pair_loop]->groupindex_1][curASM->list_pairdiff[pair_loop]->groupindex_2] += 1;
			}
		}


#ifdef OUTPUT_PAIRWISE_JSD
		for (int samplecnt = 1; samplecnt <= SAMPLE_CNT_PER_GROUP; ++samplecnt)
		{
			if (curASM->individualExpression[1][samplecnt] > THRESHOLD_MIN_ASM_COVERAGE*2 && curASM->individualExpression[2][samplecnt] > THRESHOLD_MIN_ASM_COVERAGE*2
				&& curASM->error_ratio[1][samplecnt] < 0.1 && curASM->error_ratio[2][samplecnt] < 0.1)
				d_statistics_outputfile << "\t" << sqrt(calculate_JSD(curASM->individualPathProportion[1][samplecnt], curASM->individualPathProportion[2][samplecnt], curASM->pathNum));
			else
				d_statistics_outputfile << "\t0";
		}
#endif

		d_statistics_outputfile << endl;
	}
	d_statistics_outputfile.close();


	filename = resultPath + "differential_transcription_pairwise_table.txt";
	ofstream diff_trans_pair_table(filename.c_str());
	for (int group_loop = 1; group_loop <= global_total_group_num; ++group_loop)
	{
		for (int group_loop_2 = 1; group_loop_2 <= global_total_group_num; ++group_loop_2)
		{
			diff_trans_pair_table << pair_diff_cnt[group_loop][group_loop_2] << "\t";
		}
		diff_trans_pair_table << endl;
	}
	diff_trans_pair_table.close();
	for (int group_loop = 1; group_loop <= global_total_group_num; ++group_loop)
		delete [] pair_diff_cnt[group_loop];
	delete [] pair_diff_cnt;

	//cout << "under FDR(median)<=" << false_discovery_rate << ", requiring stat_difftrans >= " << thresh_stat_d << " and sqrtJSD >= " << thresh_JSD << ", select " << significantCnt << " ASMs with significant change on transcription" << endl;
	cout << "under FDR(median)<=" << false_discovery_rate << ", requiring sqrtJSD >= " << thresh_JSD << ", select " << significantCnt << " ASMs with significant change on transcription" << endl;
}


// void output_samplesample_jsd_matrix(string resultPath)
// {
// 	long output_asm_count = 0;
// 	string filename;
// 	ofstream outfile_matrix, outfile_summary, outfile_;
// 	
// 	 = resultPath + "/sample_jsd_matrix.txt";
// 	(filename.c_str());
// 
// 	for (unsigned long asm_loop = 1; asm_loop < sortList_ASM.size(); ++asm_loop)
// 	{
// 		if (sortList_ASM[asm_loop]->significant == false)
// 			continue;
// 
// 		++output_asm_count;
// 
// 		for (int sample1 = 1; sample1 <= global_total_sample_num; ++sample1)
// 		{
// 			for (int sample2 = 1; sample2 <= global_total_sample_num; ++sample2)
// 			{
// 				 outfile_matrix << sqrt(calculate_JSD(sortList_ASM[asm_loop]->sample_path_prop[sample1], sortList_ASM[asm_loop]->sample_path_prop[sample2], sortList_ASM[asm_loop]->pathNum)) << "\t";
// 			}
// 			outfile_matrix << endl;
// 		}
// 	}
// 
// 	return;
// }

// 
// void output_for_clustering(string resultPath)
// {
// 	long output_asm_count = 0;
// 	ASM *curASM;
// 	string filename;
// 	ofstream outfile_asm, outfile_summary;
// 
// 	filename = resultPath + "/asm_for_clustering.txt";
// 	outfile_asm.open(filename.c_str());
// 
// 	for (unsigned long asm_loop = 1; asm_loop < sortList_ASM.size(); ++asm_loop)
// 	{
// 		if (sortList_ASM[asm_loop]->significant == false)
// 			continue;
// 
// 		++output_asm_count;
// 		curASM = sortList_ASM[asm_loop];
// 		double stat_difftrans = fabs(curASM->statistic_d - curASM->statistic_d_expected);
// 		double jsd = sqrt(curASM->max_pairwise_JSD_betwgroup);
// 
// 		outfile_asm << curASM->ASM_id << "\t"  << curASM->chrNM << ":" << curASM->rangeLow << "-" << curASM->rangeHigh << "\t" << curASM->in_gene << "\t" << curASM->pathNum << '\t' << alterSpliceCategory[curASM->ASMcategory] << "\t"
// 			<< stat_difftrans << "\t" << jsd << "\t" << curASM->max_grouppair_min_dist_ratio << "\t" << curASM->max_percent_large_dist << '\t' << curASM->grand_mean_expr << endl;
// 		for (int path_loop = 1; path_loop <= curASM->pathNum; ++path_loop)
// 		{
// 			outfile_asm << curASM->ASM_id << "_" << path_loop;
// 			for (long sampleLoopCnt = 1; sampleLoopCnt <= global_total_sample_num; ++sampleLoopCnt)
// 			{
// 				outfile_asm << "\t" << curASM->sample_path_prop[sampleLoopCnt][path_loop];
// 			}
// 			outfile_asm << endl;
// 		}
// 	}
// 
// 	return;
// }



void output_for_clustering(string resultPath)
{
	long output_asm_count = 0;
	ASM *curASM;
	string filename;
	ofstream outfile_asm, outfile_summary;

	filename = resultPath + "/asm_for_clustering.txt";
	outfile_asm.open(filename.c_str());
	filename = resultPath + "/asm_info.txt";
	outfile_summary.open(filename.c_str());

	for (unsigned long asm_loop = 1; asm_loop < sortList_ASM.size(); ++asm_loop)
	{
		curASM = sortList_ASM[asm_loop];

		if (curASM->significant_stat == false || sqrt(curASM->max_pairwise_JSD_betwgroup) < 0.1)// || curASM->max_percent_large_dist < 0.05)
			continue;

		++output_asm_count;
		double stat_difftrans = fabs(curASM->statistic_d - curASM->statistic_d_expected);
		double jsd = sqrt(curASM->max_pairwise_JSD_betwgroup);

		outfile_summary << asm_loop << "\t" << curASM->ASM_id << "\t"  << curASM->chrNM << ":" << curASM->rangeLow << "-" << curASM->rangeHigh << "\t" << curASM->in_gene << "\t" << curASM->pathNum << '\t' << alterSpliceCategory[curASM->ASMcategory] << "\t"
			<< stat_difftrans << "\t" << jsd << "\t" << curASM->max_grouppair_min_dist_ratio << "\t" << curASM->max_percent_large_dist << '\t' << curASM->grand_mean_expr << endl;
		for (int path_loop = 1; path_loop <= curASM->pathNum; ++path_loop)
		{
			outfile_asm << curASM->sample_path_prop[1][path_loop];
			for (long sampleLoopCnt = 2; sampleLoopCnt <= global_total_sample_num; ++sampleLoopCnt)
			{
				outfile_asm << "\t" << curASM->sample_path_prop[sampleLoopCnt][path_loop];
			}
			outfile_asm << endl;
		}
	}

	return;
}

void output_for_visualization(string resultPath)
{
	ASM *curASM;
	string filename;
	ofstream outfile_vis;

	filename = resultPath + "/asm_path_expr_for_vis.txt";
	outfile_vis.open(filename.c_str());

	for (unsigned long asm_loop = 1; asm_loop < sortList_ASM.size(); ++asm_loop)
	{
		curASM = sortList_ASM[asm_loop];

		outfile_vis << curASM->ASM_id;
		for (long sampleLoopCnt = 1; sampleLoopCnt <= global_total_sample_num; ++sampleLoopCnt)
		{
			outfile_vis << "\t" << curASM->sample_expr[sampleLoopCnt];
		}
		outfile_vis << endl;
				
		for (int path_loop = curASM->pathNum; path_loop >= 1; --path_loop)
		{
			//reverse the order because the gtree output for asm paths is from the last path to the first one
			outfile_vis << curASM->path_id[path_loop-1];
			for (long sampleLoopCnt = 1; sampleLoopCnt <= global_total_sample_num; ++sampleLoopCnt)
			{
				outfile_vis << "\t" << curASM->sample_path_expr[sampleLoopCnt][path_loop];
			}
			outfile_vis << endl;
		}
	}

	return;
}


void output_rep_path(string resultPath)
{
	ASM *curASM;
	string filename;
	ofstream outfile_vis;

	filename = resultPath + "/asm_representative_path.txt";
	outfile_vis.open(filename.c_str());

	outfile_vis << "asmID\tpathID\tgene\tsig";
	for (long sampleLoopCnt = 1; sampleLoopCnt <= global_total_sample_num; ++sampleLoopCnt)
		outfile_vis << "\ts" << sampleLoopCnt;
	outfile_vis << endl;

	for (unsigned long asm_loop = 1; asm_loop < sortList_ASM.size(); ++asm_loop)
	{
		curASM = sortList_ASM[asm_loop];

		outfile_vis << curASM->ASM_id << "\t" << curASM->path_id[curASM->rep_path-1] << "\t" << curASM->in_gene << "\t" << curASM->significant_both;
		for (long sampleLoopCnt = 1; sampleLoopCnt <= global_total_sample_num; ++sampleLoopCnt)
		{
			outfile_vis << "\t" << curASM->sample_path_expr[sampleLoopCnt][curASM->rep_path];
		}
		outfile_vis << endl;
	}

	return;
}


void write_output(string resultPath)
{
	calculate_max_pairwise_JSD_betwgroup();
	calculate_minmax_dist_to_groupmean_ratio();
	calculate_max_percent_large_dist();
	assign_gene_for_asm();

	output_difftrans_table(resultPath);
	output_for_clustering(resultPath + "/stat/temp/");
	output_for_visualization(resultPath + "/stat/temp/");
	output_rep_path(resultPath);

	return;
}
