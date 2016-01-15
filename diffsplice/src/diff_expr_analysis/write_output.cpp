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
	gene *cur_asm;

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

void calculate_max_pairwise_foldchange_betwgroup()
{
	for (unsigned long asm_loop = 1; asm_loop < sortList_ASM.size(); ++asm_loop)
	{
		double cur_max_fc = 0;
		string cur_max_fc_groupname = "", cur_large_fc_info = "";

		for (int group1 = 1; group1 <= global_total_group_num; ++group1)
		{
			for (int group2 = 1; group2 <= global_total_group_num; ++group2)
			{
				// do all pairs, including repeated pairs like 1-2 and 2-1, to find the maximum fold change
				if (group1 == group2)
					continue;

				double cur_fc = sortList_ASM[asm_loop]->group_mean_expr[group2] / sortList_ASM[asm_loop]->group_mean_expr[group1];
				if (cur_fc > cur_max_fc)
				{
					cur_max_fc = cur_fc;
					cur_max_fc_groupname = global_group_name[group1] + "->" + global_group_name[group2];
				}

				if (cur_fc > thresh_foldchange_up)
				{
					cur_large_fc_info += global_group_name[group1] + "->" + global_group_name[group2] + "," + ntostr(cur_fc) + "; ";
				}
			}
		}

		sortList_ASM[asm_loop]->max_pairwise_foldchange = cur_max_fc;
		sortList_ASM[asm_loop]->max_pairwise_foldchange_groupname = cur_max_fc_groupname;
		sortList_ASM[asm_loop]->info_large_pairwise_foldchange = cur_large_fc_info;
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
		sortKey_ASM[asm_loop] = sortList_ASM[asm_loop]->statistic_d;
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
	filename = resultPath + "differential_expression.txt";
	ofstream d_statistics_outputfile(filename.c_str());

	long **pair_diff_cnt_up = new long* [global_total_group_num+1];
	for (int group_loop = 1; group_loop <= global_total_group_num; ++group_loop)
	{
		pair_diff_cnt_up[group_loop] = new long [global_total_group_num+1];
		for (int group_loop_2 = 1; group_loop_2 <= global_total_group_num; ++group_loop_2)
		{
			pair_diff_cnt_up[group_loop][group_loop_2] = 0;
		}
	}

	gene *curASM;
	long significantCnt = 0;
	double stat_diffexpr;
	d_statistics_outputfile << "location\tgene_id\tgene\tstat_diff_expr(|stat_d_expected-stat_d|)\tmax pairwise fold change\tmax pairwise fold change group names\tgrand_mean_cov\tgroup_mean_cov\tsignificant\tlarge pairwise fold change groups\tsignificant changes between two groups" <<endl;
	for (asm_loop = 1; asm_loop < sortList_ASM.size(); ++asm_loop)
	{
		curASM = sortList_ASM[asm_loop];
		stat_diffexpr = fabs(curASM->statistic_d - curASM->statistic_d_expected);
		
		d_statistics_outputfile << curASM->chrNM << " " << curASM->rangeLow << " " << curASM->rangeHigh << "\t" << curASM->gene_id << "\t" << curASM->in_gene << "\t" 
			<< stat_diffexpr << "\t" << curASM->max_pairwise_foldchange << "\t" << curASM->max_pairwise_foldchange_groupname << "\t" << curASM->grand_mean_expr << "\t";
		for (int group_loop = 1; group_loop <= global_total_group_num; ++group_loop)
			d_statistics_outputfile << curASM->group_mean_expr[group_loop] << ", ";
		d_statistics_outputfile << "\t";
		if (stat_diffexpr >= significance_cutoff && curASM->max_pairwise_foldchange >= thresh_foldchange_up && stat_diffexpr >= thresh_stat_d)
		{
			curASM->significant_both = true;
			d_statistics_outputfile << "yes\t";
			++significantCnt;
		}
		else
			d_statistics_outputfile << "no\t";
		d_statistics_outputfile << curASM->info_large_pairwise_foldchange;

		d_statistics_outputfile << "\t";
		for (unsigned long pair_loop = 0; pair_loop < curASM->list_pairdiff.size(); ++pair_loop)
		{
			if (curASM->list_pairdiff[pair_loop]->significant_both)
			{
				d_statistics_outputfile << global_group_name[curASM->list_pairdiff[pair_loop]->groupindex_1] << "-" << global_group_name[curASM->list_pairdiff[pair_loop]->groupindex_2] << " ("
					<< curASM->list_pairdiff[pair_loop]->stat_diffexpr << "," << curASM->list_pairdiff[pair_loop]->foldchange << "); ";

				++(pair_diff_cnt_up[curASM->list_pairdiff[pair_loop]->groupindex_1][curASM->list_pairdiff[pair_loop]->groupindex_2]);
			}
		}

		d_statistics_outputfile << endl;
	}
	d_statistics_outputfile.close();

	filename = resultPath + "differential_expression_pairwise_table_up.txt"; ofstream diff_trans_pair_table_up(filename.c_str());
	for (int group_loop = 1; group_loop <= global_total_group_num; ++group_loop)
	{
		for (int group_loop_2 = 1; group_loop_2 <= global_total_group_num; ++group_loop_2)
		{
			diff_trans_pair_table_up << pair_diff_cnt_up[group_loop][group_loop_2] << "\t";
		}
		diff_trans_pair_table_up << endl;
	}
	diff_trans_pair_table_up.close();

	for (int group_loop = 1; group_loop <= global_total_group_num; ++group_loop)
	{
		delete [] pair_diff_cnt_up[group_loop];
	}
	delete [] pair_diff_cnt_up;

	cout << "under FDR(median)<=" << false_discovery_rate << ", requiring stat_difftrans >= " << thresh_stat_d << " and fold change >= " << thresh_foldchange_up << ", select " << significantCnt << " loci with significant change on expression" << endl;
}


void output_for_visualization(string resultPath)
{
	gene *curASM;
	string filename;
	ofstream outfile_vis;

	filename = resultPath + "/gene_expr_for_vis.txt";
	outfile_vis.open(filename.c_str());

	for (unsigned long asm_loop = 1; asm_loop < sortList_ASM.size(); ++asm_loop)
	{
		curASM = sortList_ASM[asm_loop];
		if (curASM->in_gene.compare("na") == 0)
			continue;
		
		outfile_vis << curASM->gene_id;
		for (long sampleLoopCnt = 1; sampleLoopCnt <= global_total_sample_num; ++sampleLoopCnt)
		{
			outfile_vis << "\t" << curASM->sample_expr[sampleLoopCnt];
		}
		outfile_vis << endl;
	}

	return;
}

void write_output(string resultPath)
{
	calculate_max_pairwise_foldchange_betwgroup();
	assign_gene_for_asm();

	output_difftrans_table(resultPath);
	output_for_visualization(resultPath + "/stat/temp/");

	return;
}
