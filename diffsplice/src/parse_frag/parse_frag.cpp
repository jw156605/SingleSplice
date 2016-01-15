#include "parse_frag.h"

int record::rec_compare(const record &ref_record, bool is_exon)
{
	if (this->start_pos < ref_record.start_pos || this->start_pos == ref_record.start_pos && this->end_pos < ref_record.end_pos)
	{
		return -1;
	}
	else if (this->start_pos == ref_record.start_pos && this->end_pos == ref_record.end_pos)
	{
		if (is_exon)
			return 0;
		else
		{
			if (this->xs < ref_record.xs)
				return -1;
			else if (this->xs == ref_record.xs)
				return 0;
			else
				return 1;
		}
	}
	else
	{
		return 1;
	}
}


string itostr(int t)
{
	ostringstream oss;
	oss << t;
	return oss.str();
}


void parse(string inputfilename, string outputfilename_prefix, int sample_index)
{
	ifstream inputfile;
	inputfile.open(inputfilename.c_str());

	ofstream exonFile, junctionFile;
	exonFile.open((outputfilename_prefix + "_exonic.txt").c_str());
	junctionFile.open((outputfilename_prefix + "_junction.txt").c_str());

	string info, cigar, xs;
	long startPoint, endPoint, tmp, strLoop;
	long totalReadNum = 0;
	long stat_MinStart = MAX, stat_MaxEnd = 0;
	
	while (inputfile >> startPoint)
	{
		inputfile >> cigar;
		inputfile >> xs;
		getline(inputfile, info);

		++totalReadNum;
		if (cigar.compare("*") != 0)
		{
			tmp = 0;

			if (xs.compare("-") == 0)
				xs = "0";
			else
				xs = "1";

			if (startPoint < stat_MinStart)
				stat_MinStart = startPoint;

			for (strLoop = 0; strLoop < cigar.size(); ++strLoop)
			{
				if (cigar[strLoop] == 'M')
				{
					endPoint = startPoint + tmp;
					exonFile << startPoint << '\t' << endPoint - 1 << endl;
					startPoint = endPoint;
					tmp = 0;
				} 
				else if (cigar[strLoop] == 'N')
				{
					endPoint = startPoint + tmp;
					junctionFile << startPoint - 1 << '\t' << endPoint << "\t" << xs << endl;
					startPoint = endPoint;
					tmp = 0;
				}
				else if (cigar[strLoop] == 'I')
				{
					tmp = 0;
				}
				else if ((cigar[strLoop] == 's') || (cigar[strLoop] == 'S'))
				{
					tmp = 0;
				}
				else if (cigar[strLoop] == 'D')
				{
					endPoint = startPoint + tmp;
					startPoint = endPoint;
					tmp = 0;
				}
				else if (cigar[strLoop] >= '0' && cigar[strLoop] <= '9')
				{
					tmp = tmp * 10 + cigar[strLoop] - 48;
				}
				else
				{
					cout << "unrecognized CIGAR: " << cigar[strLoop] << " at sample " << sample_index << " (" << startPoint << ", " << cigar << ")" << endl;
					tmp = 0;
				}
			}

			if (endPoint - 1 > stat_MaxEnd)
				stat_MaxEnd = endPoint - 1;
		}
	}
	
	statFile << sample_index << "\t" << stat_MinStart << '\t' << stat_MaxEnd << "\t" << totalReadNum << endl;

	return;
}



//merge sorted files, add up occurrences
void merge_sorted_files(int num_sample, string filename_prefix, string filename_suffix, bool is_exon)
{
	if (num_sample > 1000)
	{
		cout << "warning: will open more than 1000 files at the same time." << endl;
	}

	vector <ifstream*> vec_infile;
	ofstream outfile;
	string filename, info;
	long sample_loop;

	vector <record*> vec_record;
	vector <bool> vec_file_ongoing (num_sample, true);
	record cur_min_record; //record for current minimum record
	vector <int> vec_min_index; //vector storing indices of files that have the minimum record
	vector <long> vec_support_persample;

	for (sample_loop = 0; sample_loop < num_sample; ++sample_loop)
	{
		filename = filename_prefix + "/" + itostr(sample_loop+1) + filename_suffix;
		ifstream *new_infile = new ifstream;
		new_infile->open(filename.c_str());
		if (!new_infile->is_open())
		{
			cout << "error: fail to open input file " << filename << endl;
			return;
		}
		else
			vec_infile.push_back(new_infile);

		record *new_record = new record;
		vec_record.push_back(new_record);
	}

	filename = filename_prefix + "/merged" + filename_suffix;
	outfile.open(filename.c_str());
	if (!outfile.is_open())
		cout << "warning: fail to open output file " << filename << endl;


	for (sample_loop = 0; sample_loop < num_sample; ++sample_loop)
	{
		*vec_infile[sample_loop] >> vec_record[sample_loop]->occurrence;
		*vec_infile[sample_loop] >> vec_record[sample_loop]->start_pos;
		*vec_infile[sample_loop] >> vec_record[sample_loop]->end_pos;
		if (!is_exon)
			*vec_infile[sample_loop] >> vec_record[sample_loop]->xs;
		getline(*vec_infile[sample_loop], info);
	}

	int num_finished_sample = 0, compare_result;
	while (num_finished_sample < num_sample)
	{
		//find the smallest record
		cur_min_record.start_pos = MAX;
		cur_min_record.end_pos = MAX;
		if (!is_exon)
			vec_support_persample.assign(num_sample, 0);

		for (sample_loop = 0; sample_loop < num_sample; ++sample_loop)
		{
			if (!vec_file_ongoing[sample_loop])
				continue;

			compare_result = vec_record[sample_loop]->rec_compare(cur_min_record, is_exon);
			if (compare_result < 0)
			{
				//this record is < current minimum record, update current minimum record
				cur_min_record = *(vec_record[sample_loop]);
				if (!is_exon)
				{
					vec_support_persample.assign(num_sample, 0);
					vec_support_persample[sample_loop] = vec_record[sample_loop]->occurrence;
				}
				vec_min_index.clear();
				vec_min_index.push_back(sample_loop);
			}
			else if (compare_result == 0)
			{
				//update occurrence of the current minimum record, and add this file
				cur_min_record.occurrence += vec_record[sample_loop]->occurrence;
				if (!is_exon)
					vec_support_persample[sample_loop] = vec_record[sample_loop]->occurrence;
				vec_min_index.push_back(sample_loop);
			}
		}

		//output the smallest record
		if (is_exon)
		{
			outfile << cur_min_record.occurrence << "\t" << cur_min_record.start_pos << "\t" << cur_min_record.end_pos << endl;

			for (sample_loop = 0; sample_loop < vec_min_index.size(); ++sample_loop)
			{
				if (*(vec_infile[vec_min_index[sample_loop]]) >> vec_record[vec_min_index[sample_loop]]->occurrence
					&& *(vec_infile[vec_min_index[sample_loop]]) >> vec_record[vec_min_index[sample_loop]]->start_pos
					&& *(vec_infile[vec_min_index[sample_loop]]) >> vec_record[vec_min_index[sample_loop]]->end_pos
					&& getline(*(vec_infile[vec_min_index[sample_loop]]), info))
				{
					//file is not ended
				}
				else
				{
					vec_file_ongoing[vec_min_index[sample_loop]] = false;
					++num_finished_sample;
				}
			}
		} 
		else
		{
			outfile << cur_min_record.occurrence << "\t" << cur_min_record.start_pos << "\t" << cur_min_record.end_pos << "\t" << cur_min_record.xs;
			for (sample_loop = 0; sample_loop < num_sample; ++sample_loop)
				outfile << "\t" << vec_support_persample[sample_loop];
			outfile << endl;

			for (sample_loop = 0; sample_loop < vec_min_index.size(); ++sample_loop)
			{
				if (*(vec_infile[vec_min_index[sample_loop]]) >> vec_record[vec_min_index[sample_loop]]->occurrence
					&& *(vec_infile[vec_min_index[sample_loop]]) >> vec_record[vec_min_index[sample_loop]]->start_pos
					&& *(vec_infile[vec_min_index[sample_loop]]) >> vec_record[vec_min_index[sample_loop]]->end_pos
					&& *(vec_infile[vec_min_index[sample_loop]]) >> vec_record[vec_min_index[sample_loop]]->xs 
					&& getline(*(vec_infile[vec_min_index[sample_loop]]), info))
				{
					//file is not ended
				}
				else
				{
					vec_file_ongoing[vec_min_index[sample_loop]] = false;
					++num_finished_sample;
				}
			}
		}
	}


	for (sample_loop = 0; sample_loop < num_sample; ++sample_loop)
	{
		vec_infile[sample_loop]->close();
		delete vec_infile[sample_loop];
		delete vec_record[sample_loop];
	}
	outfile.close();

	return;
}



void calculate_normalization_ratio(string dir_input_read, string dir_output_frag, int num_sample)
{
	vector <double> vec_sample_size(num_sample, 0);
	double max_sample_size = 0;
	string filename, info;
	ifstream infile_stat;

	for (int sample_cnt = 1; sample_cnt <= num_sample; ++sample_cnt)
	{
		filename = dir_input_read + "/" + itostr(sample_cnt) + "/DatasetStat.txt";
		infile_stat.open(filename.c_str());
		if (!infile_stat.is_open())
		{
			cout << "error: fail to open file " << filename << endl;
			exit(1);
		}

		infile_stat >> vec_sample_size[sample_cnt-1];
		getline(infile_stat, info);

		if (vec_sample_size[sample_cnt-1] > max_sample_size)
			max_sample_size = vec_sample_size[sample_cnt-1];

		infile_stat.close();
	}

	filename = dir_output_frag + "/Normalization_Ratio.txt";
	ofstream outfile_normalize(filename.c_str());
	if (!outfile_normalize.is_open())
	{
		cout << "error: fail to open file " << filename << endl;
		exit(1);
	}

	for (int sample_cnt = 0; sample_cnt < num_sample; ++sample_cnt)
	{
		outfile_normalize << max_sample_size / vec_sample_size[sample_cnt] << endl;
	}

	outfile_normalize.close();


	return;
}

//process all samples for a chromosome
int main(int argc, char* argv[])
{
	string path_input_read, chromosome_name, path_output_frag;
	int num_sample;

#ifdef UNIX
	if (argc != 5)
	{
		cout << argv[0] << "\t<readfile_path>\t<chr>\t<output_folder>\t<num_of_samples>" << endl;
		return 1;
	}
	path_input_read = argv[1];
	chromosome_name = argv[2];
	path_output_frag = argv[3];
	num_sample = atoi(argv[4]);
#else
	path_input_read = "";
	chromosome_name = "chr16";
	path_output_frag = "frag";
	num_sample = 3;
#endif
	
	string infilename, outfilename_prefix, comd;
	infilename = path_output_frag + "/stat.txt";
	statFile.open(infilename.c_str());

	for (int sample_cnt = 1; sample_cnt <= num_sample; ++sample_cnt)
	{
		infilename = path_input_read + "/" + itostr(sample_cnt) + "/" + chromosome_name + ".txt";
		outfilename_prefix = path_output_frag + "/" + itostr(sample_cnt);
		parse(infilename, outfilename_prefix, sample_cnt);

		comd = "sort -n " + outfilename_prefix + "_junction.txt | uniq -c > " + outfilename_prefix + "_junction_sort.txt";
		system(comd.c_str());
		comd = "sort -n " + outfilename_prefix + "_exonic.txt | uniq -c > " + outfilename_prefix + "_exonic_sort.txt";
		system(comd.c_str());

		comd = "rm -f " + outfilename_prefix + "_junction.txt";
		system(comd.c_str());
		comd = "rm -f " + outfilename_prefix + "_exonic.txt";
		system(comd.c_str());
	}

	statFile.close();

	merge_sorted_files(num_sample, path_output_frag, "_exonic_sort.txt", true);
	merge_sorted_files(num_sample, path_output_frag, "_junction_sort.txt", false);


	calculate_normalization_ratio(path_input_read, path_output_frag, num_sample);

	return 0;
}
