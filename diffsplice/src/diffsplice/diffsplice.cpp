/*    
 *    diffsplice.cpp		
 *    DiffSplice
 *
 *    Copyright (C) 2012 University of Kentucky and
 *                       Yin Hu
 *
 *    Authors: Yin Hu
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


#include "diffsplice.h"

sample::sample(string samplename, int groupid, int sampleid, string datafilename)
{
	this->name = samplename;
	this->groupID = groupid;
	this->sampleID_orig = sampleid;
	this->sampleID_dfs = 0;
	this->datafilename = datafilename;
}

sample::~sample() {}

sampleGroup::sampleGroup(string groupname, int ID)
{
	name = groupname;
	groupID = ID;
}

sampleGroup::~sampleGroup() {}

// decide the steps to run
void set_run_steps(RUNMODE runmode)
{
	if (runmode == run_full || runmode == run_default)
	{
		// run through all steps
		PREPROCESS = true;
	}
	else if (runmode == run_rerun)
	{
		// do not need to preprocess again
		PREPROCESS = false;
	}
	SEP_TREE = true;
	SELECT_SIGNIFICANCE = true;

	return;
}


// print the usage of this program
void print_usage()
{
	cout << "Usage: diffsplice [-o output_directory] [-m full/rerun] [-s settings]"
		 <<	" [-a annotated_splice_junctions] [-g gaf_annotation] [-r runname] <data_manifest>" << endl;
}

int main(int argc, char* argv[])
{
	time_t t;
	string target_path = "./", anno_splice_path = "no_anno", runName = "dfs_run", gaf_file = "no_file", settings_file, data_file; // control strings
	string cur_time, comd, file_path, result_path, bin_path = "./", filename, chrlist; // string variables
	int opt_case = 0;
	RUNMODE run_mode = run_default;

	time(&t); cur_time = ctime(&t); cur_time.erase(cur_time.size() - 1);
	cout << "[" << cur_time << "]\tDiffSplice 1.0.0" << endl;

	// try to determine the path of the binaries
	string tmpstring = argv[0];
	unsigned found = tmpstring.rfind("/");
	if (found != string::npos)
		bin_path = tmpstring.substr(0, found+1);

	while (1)
	{
		/* getopt_long stores the option index here. */
		int option_index = 0;
		opt_case = getopt_long (argc, argv, "hm:p:o:a:r:s:g:c:", long_options, &option_index);

		/* Detect the end of the options. */
		if (opt_case == -1)
			break;

		switch (opt_case)
		{
		case 0:
			/* If this option set a flag, do nothing else now. */
			if (long_options[option_index].flag != 0)
				break;
			printf ("option %s", long_options[option_index].name);
			if (optarg)
				printf (" with arg %s", optarg);
			printf ("\n");
			break;

		case 'm':
			printf ("Run mode = %s\n", optarg);
			if (strcmp(optarg, "full") == 0)
				run_mode = run_full;
			else if (strcmp(optarg, "rerun") == 0)
				run_mode = run_rerun;
			else
			{
				printf ("error: unrecognized running mode %s\n", optarg);
				run_mode = run_nothing;
			}
			break;

		case 'o':
			printf ("Results will be placed in folder %s\n", optarg);
			target_path = optarg;
			if (target_path[target_path.length()-1] != '/')
				target_path += "/";
			break;

		case 'a':
			printf ("Annotated splice junctions provided as in `%s'\n", optarg);
			anno_splice_path = optarg;
			if (anno_splice_path[anno_splice_path.length()-1] != '/')
				anno_splice_path += "/";
			break;

		case 'g':
			printf ("GAF annotation provided as `%s'\n", optarg);
			gaf_file = optarg;
			break;

		case 'r':
			printf ("Alias of this run is `%s'\n", optarg);
			runName = optarg;
			break;

		case 'c':
			printf ("Only work on chromosomes listed in `%s'\n", optarg);
			chrlist = optarg;
			break;

		case 's':
			printf ("Parameter settings provided as in `%s'\n", optarg);
			settings_file = optarg;
			break;

		case 'p':
			// do not support multi-thread yet
			break;

		case 'h':
			print_usage();
			exit(1);

		case '?':
			/* getopt_long already printed an error message. */
			break;

		default:
			print_usage();
			abort ();
		}
	}
	if (optind == argc - 1)
	{
		printf ("Take data manifest as %s\n", argv[optind]);
		data_file = argv[optind];
	}
	else if (optind == argc)
	{
		printf("No data manifest received. Please provide a list of input files.\n");
		print_usage();
		exit(1);
	}
	else
	{
		printf("Unrecognized input parameters. Please check the command line.\n");
		print_usage();
		exit(1);
	}

	if (run_mode == run_nothing)
		exit(1);

	if (settings_file.empty())
	{
		printf("No settings file received. Will use the default settings.\n");
	}
	else
	{
		parse_settings(settings_file);
	}

	parse_inputfiles(data_file);
	set_run_steps(run_mode);

	// create data folder that hosts intermediate files and internal config files
	comd = "mkdir -p " + target_path + "data"; system(comd.c_str());
	write_configfile(target_path + "data/");

	time(&t); cur_time = ctime(&t); cur_time.erase(cur_time.size() - 1);
	cout << "[" << cur_time << "]\tData set: " << runName << ", "<< allGroups.size() << " groups with a total of " << totalSampleNum << " samples." << endl;

	// pre-process bam files
	if (PREPROCESS)
	{
		time(&t); cur_time = ctime(&t); cur_time.erase(cur_time.size() - 1);
		cout << "[" << cur_time << "]\tBegin to parse BAM files" << endl;

		file_path = target_path + "data/bam_parsed";
		comd = "rm -f -r " + file_path; system(comd.c_str());
		comd = "mkdir -p " + file_path; system(comd.c_str());

		// write the bam file list
		filename = target_path + "data/bamfile_list.txt";
		ofstream outfile(filename.c_str());
		if (!outfile.is_open())
		{
			cout << "Error: cannot write to folder " << target_path << "data/." << endl;
			exit(1);
		}
		for (vector<sampleGroup>::iterator it_group = allGroups.begin(); it_group != allGroups.end(); ++it_group)
			for (vector<sample>::iterator it_sample = it_group->samples.begin(); it_sample != it_group->samples.end(); ++it_sample)
				outfile << it_sample->sampleID_dfs << "\t" << it_sample->datafilename << endl;
		outfile.close();

		comd = bin_path + "parse_bam " + filename + " " + file_path + "/"; system(comd.c_str());
	}

	// get list of chromosomes
	if (chrlist.empty())
		filename = target_path + "data/bam_parsed/ChromosomeName.txt";
	else
		filename = chrlist;
	inputChrName(filename);

	// collect expression fragment files
	if (PREPROCESS)
	{
		time(&t); cur_time = ctime(&t); cur_time.erase(cur_time.size() - 1);
		cout << "[" << cur_time << "]\tBegin to prepare input files for DiffSplice. This step may take long." << endl;

		file_path = target_path + "data/frag";
		comd = "rm -f -r " + file_path; system(comd.c_str());
		comd = "mkdir -p " + file_path; system(comd.c_str());

		for (unsigned int i = 0; i < chrName.size(); ++i)
		{
			time(&t); cur_time = ctime(&t); cur_time.erase(cur_time.size() - 1);
			cout << "[" << cur_time << "]\tPreparing " << chrName[i] << endl;

			comd = "mkdir -p " + file_path + "/" + chrName[i]; system(comd.c_str());
			comd = bin_path + "parse_frag " + target_path + "data/bam_parsed/" + " " + chrName[i] + " " + file_path+"/"+chrName[i]+"/ " + itostr(totalSampleNum);
			system(comd.c_str());
		}
	}

	if (SEP_TREE == true)
	{
		time(&t); cur_time = ctime(&t); cur_time.erase(cur_time.size() - 1);
		cout << "[" << cur_time << "]\tBegin to run DiffSplice" << endl;

		file_path = target_path + "data/frag/";
		result_path = target_path + "result/";
		comd = "rm -f -r " + result_path; system(comd.c_str());
		comd = "mkdir -p " + result_path; system(comd.c_str());
		comd = "mkdir -p " + result_path + "stat"; system(comd.c_str());
		comd = "mkdir -p " + result_path + "detail"; system(comd.c_str());
		comd = "mkdir -p " + result_path + "gene"; system(comd.c_str());
		comd = "mkdir -p " + result_path + "stat/expr"; system(comd.c_str());
		comd = "mkdir -p " + result_path + "stat/asm"; system(comd.c_str());
		comd = "mkdir -p " + result_path + "stat/temp"; system(comd.c_str());
		comd = "mkdir -p " + result_path + "asm"; system(comd.c_str());


		ofstream asm_path_file, non_asm_file, splice_all_file, splice_filter_file;
		filename = result_path + "asm_path.gtf";
		asm_path_file.open(filename.c_str());
		asm_path_file << "browser full " << runName << "_ASM_path\ntrack name=\"" << runName << "_ASM_path\" description=\"" << runName << "_ASM_path\" visibility=2 useScore=1 itemRgb=on" << endl;
		asm_path_file.close();

		filename = result_path + "splice_graph.gtf";
		non_asm_file.open(filename.c_str());
		non_asm_file << "browser dense " << runName << "_splice_graph\ntrack name=\"" << runName << "_splice_graph\" description=\"" << runName << "_splice_graph\" visibility=1 useScore=1 itemRgb=on" << endl;
		non_asm_file.close();

		filename = result_path + "asm_path_decomposed.gtf";
		asm_path_file.open(filename.c_str());
		asm_path_file << "browser full " << runName << "_ASM_decomposed\ntrack name=\"" << runName << "_ASM_decomposed\" description=\"" << runName << "_ASM_decomposed\" visibility=2 useScore=1 itemRgb=on" << endl;
		asm_path_file.close();

		filename = result_path + "stat/splice_all.bed";
		splice_all_file.open(filename.c_str());
		splice_all_file << "browser squish " << runName << "_junction_all\ntrack name=\"" << runName << "_junction_all\" description=\"" << runName << "_junction_all\" visibility=2 useScore=1 itemRgb=on" << endl;
		splice_all_file.close();

		filename = result_path + "stat/splice_filtered.bed";
		splice_filter_file.open(filename.c_str());
		splice_filter_file << "browser pack " << runName << "_junction_filtered\ntrack name=\"" << runName << "_junction_filtered\" description=\"" << runName << "_junction_filtered\" visibility=2 useScore=1 itemRgb=on" << endl;
		splice_filter_file.close();


		for (unsigned int i = 0; i < chrName.size(); ++i)
		{
			time(&t); cur_time = ctime(&t); cur_time.erase(cur_time.size() - 1); cout << "[" << cur_time << "]\t" << flush;
			comd = bin_path + "diff_gtree " + file_path + chrName[i] + "/ " + anno_splice_path + " " + result_path + " " + chrName[i] + " " + target_path + "data/config_gtree" + " " + chrLength[i];
			system(comd.c_str()); cout << endl;
		}
	}

	if (SELECT_SIGNIFICANCE == true)
	{
		result_path = target_path + "result/";

		time(&t); cur_time = ctime(&t); cur_time.erase(cur_time.size() - 1);
		cout << "[" << cur_time << "]\tBegin the differential expression analysis on genomic loci" << endl;

		comd = bin_path +  "diff_expr_analysis " + result_path + " " + target_path + "data/config_testexpr" + " " + gaf_file;
		system(comd.c_str());

		time(&t); cur_time = ctime(&t); cur_time.erase(cur_time.size() - 1);
		cout << "[" << cur_time << "]\tBegin the differential transcription analysis on ASMs" << endl;

		comd = bin_path +  "diff_asm_analysis " + result_path + " " + target_path + "data/config_testtrans" + " " + gaf_file;
		system(comd.c_str());
	}

	time(&t); cur_time = ctime(&t); cur_time.erase(cur_time.size() - 1);
	cout << "[" << cur_time << "]\tFinish DiffSplice run" << endl;

	return 0;
}



void inputChrName(string filename)
{
	ifstream infile;
	infile.open(filename.c_str());

	if (!infile.is_open())
	{
		cout << "error: fail to open chromosome name file " << filename << endl;
		exit(1);
	}

	string in_chrName, in_chrLength, info;

	while (infile >> in_chrName)
	{
		infile >> in_chrLength;
		getline(infile, info);

		if (in_chrName.empty() || in_chrName[0] == '#' || in_chrName.compare("*") == 0)
			continue;

		chrName.push_back(in_chrName);
		chrLength.push_back(in_chrLength);
	}

	return;
}

bool exist_frag(string fname_check)
{
	ifstream f_check(fname_check.c_str());
	if (f_check.good()) {
		f_check.close();
		return true;
	} else {
		f_check.close();
		return false;
	}
}

void parse_inputfiles(string filename)
{
	ifstream infile;
	infile.open(filename.c_str());

	string info, datafilename;

	if (infile.is_open())
	{
		while (infile >> info)
		{
			if (info[0] == '#' || info[0] == '\0')
			{
				// comment line
				getline(infile, info);
			}
			else
			{
				// data line

				// find the group
				unsigned long groupIndex = 0;
				for (; groupIndex < allGroups.size() && allGroups[groupIndex].name.compare(info) != 0; ++groupIndex) {}

				if (groupIndex >= allGroups.size())
				{
					// a new sample group
					sampleGroup newGroup(info, groupIndex);
					allGroups.push_back(newGroup);
				}

				// find the sample
				infile >> info;
				unsigned long sampleIndex = 0;
				for (; sampleIndex < allGroups[groupIndex].samples.size() && allGroups[groupIndex].samples[sampleIndex].name.compare(info) != 0; ++sampleIndex) {}

				if (sampleIndex < allGroups[groupIndex].samples.size())
				{
					cout << "Warning: possible duplicated sample " << info << " in group " << allGroups[groupIndex].name << ", will include this anyway." << endl;
				}

				// add this sample anyway
				infile >> datafilename;
				sample newSample(info, groupIndex, ++totalSampleNum, datafilename);
				allGroups[groupIndex].samples.push_back(newSample);

				getline(infile, info);
			}
		}
	}
	else
	{
		cout << "fail to open the data file" << endl;
		exit(1);
	}

	if (allGroups.empty())
	{
		cout << "fail to import datafile.cfg" << endl;
		exit(1);
	}

	// assign internal ID for every sample
	long internal_ID = 1;
	for (vector<sampleGroup>::iterator it_group = allGroups.begin(); it_group != allGroups.end(); ++it_group)
		for (vector<sample>::iterator it_sample = it_group->samples.begin(); it_sample != it_group->samples.end(); ++it_sample)
				it_sample->sampleID_dfs = internal_ID++;

	return;
}


void parse_settings(string filename)
{
	ifstream infile;
	infile.open(filename.c_str());

	string info;

	if (infile.is_open())
	{
		while (infile >> info)
		{
			if (info[0] == '#' || info[0] == '\0')
			{
				// comment line
				getline(infile, info);
			}
			else
			{
				//data line
				if (info.compare("ignore_minor_alternative_splicing_variants") == 0)
				{
					infile >> info;
					if (info.compare("yes") == 0)
						COUNT_MAJOR_PATHS_ONLY = true;
					else
						COUNT_MAJOR_PATHS_ONLY = false;
				}
				else if (info.compare("junctionfilter_mode") == 0)
				{
					infile >> junctionfilter_mode;
					if (junctionfilter_mode.compare("group") != 0 && junctionfilter_mode.compare("nogroup") != 0 && junctionfilter_mode.compare("none") != 0){
						cout << "Value of junctionfilter_mode not recognized in " << filename << ", will use no junction filter by default." << endl;
						junctionfilter_mode = "none";
					}
				}
				else if (info.compare("thresh_junctionfilter_groupwise") == 0)
				{
					infile >> thresh_junctionfilter_groupwise;
				}
				else if (info.compare("thresh_junctionfilter_num_present_samples") == 0)
				{
					infile >> thresh_junctionfilter_num_present_samples;
				}
				else if (info.compare("thresh_junctionfilter_present_support") == 0)
				{
					infile >> thresh_junctionfilter_present_support;
				}
				else if (info.compare("junction_annotation_as_white_list") == 0)
				{
					infile >> junction_annotation_as_white_list;
					if (junction_annotation_as_white_list.compare("yes") != 0 && junction_annotation_as_white_list.compare("no") != 0){
						cout << "Value of junction_annotation_as_white_list not recognized in " << filename << ", will use as a white list by default." << endl;
						junction_annotation_as_white_list = "yes";
					}
				}
				else if (info.compare("thresh_average_read_coverage_exon") == 0)
				{
					infile >> coverageThreshold_exon;
				}
				else if (info.compare("thresh_average_read_coverage_intron") == 0)
				{
					infile >> coverageThreshold_intron;
				}
				else if (info.compare("false_discovery_rate") == 0)
				{
					infile >> config_false_discovery_rate;
				}
				else if (info.compare("thresh_foldchange_up") == 0)
				{
					infile >> config_thresh_foldchange_up;
				}
				else if (info.compare("thresh_foldchange_down") == 0)
				{
					infile >> config_thresh_foldchange_down;
				}
				else if (info.compare("thresh_sqrtJSD") == 0)
				{
					infile >> config_thresh_sqrtJSD;
				}
				else
				{
					cout << "Warning: unrecognized option in settings.cfg - " << info << endl;
				}
				
				getline(infile, info);
			}
		}
	}
	else
	{
		cout << "Fail to open the setting file. Will try with default settings." << endl;
	}

	return;
}


string itostr(int t)
{
	ostringstream oss;
	oss << t;
	return oss.str();
}


void write_configfile_gtree(string target_path)
{
	ofstream conf_file;
	string filename;

	filename = target_path + "config_gtree";
	conf_file.open(filename.c_str());

	if (conf_file.is_open())
	{
		conf_file << "SUPPORT_VECTOR_SIZE\t" << totalSampleNum << endl; //total number of samples, for all groups
		conf_file << "junctionfilter_mode\t" << junctionfilter_mode << endl;
		conf_file << "thresh_junctionfilter_groupwise\t" << thresh_junctionfilter_groupwise << endl;
		conf_file << "thresh_junctionfilter_num_present_samples\t" << thresh_junctionfilter_num_present_samples << endl;
		conf_file << "thresh_junctionfilter_present_support\t" << thresh_junctionfilter_present_support << endl;
		conf_file << "junction_annotation_as_white_list\t" << junction_annotation_as_white_list << endl;

		conf_file << "COUNT_MAJOR_PATHS_ONLY\t";
		if (COUNT_MAJOR_PATHS_ONLY == true)
			conf_file << "yes" << endl;
		else
			conf_file << "no" << endl;
		conf_file << "coverageThreshold_exon\t" << coverageThreshold_exon << endl;
		conf_file << "coverageThreshold_intron\t" << coverageThreshold_intron << endl;

		conf_file << "global_total_group_num " << allGroups.size() << endl;
		conf_file << "global_group_size";
		for (vector<sampleGroup>::iterator it_group = allGroups.begin(); it_group != allGroups.end(); ++it_group)
			conf_file << " " << it_group->samples.size();
		conf_file << endl;

		conf_file << "global_group_base_index";
		for (vector<sampleGroup>::iterator it_group = allGroups.begin(); it_group != allGroups.end(); ++it_group)
			conf_file << " " << it_group->samples[0].sampleID_dfs - 1;
		conf_file << endl;

		conf_file << "group_name";
		for (vector<sampleGroup>::iterator it_group = allGroups.begin(); it_group != allGroups.end(); ++it_group)
			conf_file << " " << it_group->name;
		conf_file << endl;
	}
	else
	{
		cout << "Error: cannot write file " << filename << ". Please check permission on result folder." << endl;
		exit(1);
	}

	conf_file.close();

	return;
}

void write_configfile_testexpr(string target_path)
{
	ofstream conf_file;
	string filename;

	filename = target_path + "config_testexpr";
	conf_file.open(filename.c_str());

	if (conf_file.is_open())
	{
		conf_file << "global_total_sample_num " << totalSampleNum << endl;
		conf_file << "global_total_group_num " << allGroups.size() << endl;
		conf_file << "num_samples_mask 0" << endl;
		conf_file << "COUNT_INTRON_RETENTION " << (COUNT_INTRON_RETENTION ? "true" : "false") << endl;

		conf_file << "global_group_size";
		for (vector<sampleGroup>::iterator it_group = allGroups.begin(); it_group != allGroups.end(); ++it_group)
			conf_file << " " << it_group->samples.size();
		conf_file << endl;

		conf_file << "global_group_base_index";
		for (vector<sampleGroup>::iterator it_group = allGroups.begin(); it_group != allGroups.end(); ++it_group)
			conf_file << " " << it_group->samples[0].sampleID_dfs - 1;
		conf_file << endl;

		conf_file << "group_name";
		for (vector<sampleGroup>::iterator it_group = allGroups.begin(); it_group != allGroups.end(); ++it_group)
			conf_file << " " << it_group->name;
		conf_file << endl;

		conf_file << "false_discovery_rate " << config_false_discovery_rate << endl;
		conf_file << "thresh_foldchange_up " << config_thresh_foldchange_up << endl;
		conf_file << "thresh_foldchange_down " << config_thresh_foldchange_down << endl;
		conf_file << "thresh_stat_d " << thresh_stat_d << endl;
		conf_file << "THRESHOLD_MIN_ASM_COVERAGE " << THRESHOLD_MIN_EXPR_COVERAGE << endl; //at least one group must have mean expression larger than this threshold to make differential expression;
	}
	else
	{
		cout << "Error: cannot write file " << filename << ". Please check permission on result folder." << endl;
		exit(1);
	}

	conf_file.close();

	return;
}


void write_configfile_testtrans(string target_path)
{
	ofstream conf_file;
	string filename;

	filename = target_path + "config_testtrans";
	conf_file.open(filename.c_str());

	if (conf_file.is_open())
	{
		conf_file << "global_total_sample_num " << totalSampleNum << endl;
		conf_file << "global_total_group_num " << allGroups.size() << endl;
		conf_file << "num_samples_mask 0" << endl;
		conf_file << "COUNT_INTRON_RETENTION " << (COUNT_INTRON_RETENTION ? "true" : "false") << endl;

		conf_file << "global_group_size";
		for (vector<sampleGroup>::iterator it_group = allGroups.begin(); it_group != allGroups.end(); ++it_group)
			conf_file << " " << it_group->samples.size();
		conf_file << endl;

		conf_file << "global_group_base_index";
		for (vector<sampleGroup>::iterator it_group = allGroups.begin(); it_group != allGroups.end(); ++it_group)
			conf_file << " " << it_group->samples[0].sampleID_dfs - 1;
		conf_file << endl;

		conf_file << "group_name";
		for (vector<sampleGroup>::iterator it_group = allGroups.begin(); it_group != allGroups.end(); ++it_group)
			conf_file << " " << it_group->name;
		conf_file << endl;

		conf_file << "false_discovery_rate " << config_false_discovery_rate << endl;
		conf_file << "thresh_JSD " << config_thresh_sqrtJSD << endl;
		conf_file << "thresh_stat_d " << thresh_stat_d << endl;
		conf_file << "THRESHOLD_MIN_ASM_COVERAGE " << THRESHOLD_MIN_ASM_COVERAGE << endl; //all groups must have mean expression larger than this threshold to make differential transcription
	}
	else
	{
		cout << "Error: cannot write file " << filename << ". Please check permission on result folder." << endl;
		exit(1);
	}

	conf_file.close();

	return;
}

void write_configfile(string target_path)
{
	write_configfile_gtree(target_path);
	write_configfile_testexpr(target_path);
	write_configfile_testtrans(target_path);

	return;
}

