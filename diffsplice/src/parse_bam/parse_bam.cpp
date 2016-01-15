#include <stdio.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include "sam.h"


using namespace std;

class chr_entry
{
public:
	string chr_name;
	uint32_t chr_length;
};


vector <chr_entry*> targetChromosome; //chromosomes in this file, cleaned before processing a bam file
vector <chr_entry*> allChromosome; //chromosomes in all files

vector <unsigned long> readAlignIndex;
unsigned long totalReadNum;
vector <unsigned long> totalReadNum_persample;
vector <unsigned long> totalAlignNum_persample;
vector <string> fileindex_persample;

vector <ofstream*> outputfile;

string dirPrefix;



string itostr(int t)
{
	ostringstream oss;
	oss << t;
	return oss.str();
}

string get_cigar_str(uint32_t *cigar_array, uint32_t n_cigar)
{
	string cigar_str;

	for (unsigned int i = 0; i < n_cigar; ++i) 
	{
		int cigar_oper = cigar_array[i] & BAM_CIGAR_MASK; // operation
		int cigar_length = cigar_array[i] >> BAM_CIGAR_SHIFT; // length
		
		cigar_str += itostr(cigar_length);

		switch(cigar_oper)
		{
			case BAM_CMATCH: cigar_str += "M"; break;
			case BAM_CINS: cigar_str += "I"; break;
			case BAM_CDEL: cigar_str += "D"; break;
			case BAM_CREF_SKIP: cigar_str += "N"; break;
			case BAM_CSOFT_CLIP: cigar_str += "S"; break;
			case BAM_CHARD_CLIP: cigar_str += "R"; break;
			case BAM_CPAD: cigar_str += "P"; break;
			default: cout << "warning: uncataloged CIGAR " << cigar_oper; break;
		}
	}

	return cigar_str;
}

static int fetch_func(const bam1_t *b, void *data)
{
	++totalReadNum; //add 1 to the total number of reads

	samfile_t *fp = (samfile_t*)data;
	
	const bam1_core_t *c = &b->core;
	if (c->flag&BAM_FUNMAP) return 0; /* skip unmapped reads */

	//output a single read
	string chromosome, outputfilename;
	int tmpLoop;
	chr_entry *new_chr_entry;

	chromosome = fp->header->target_name[c->tid];

	for (tmpLoop = targetChromosome.size() - 1; tmpLoop >= 0 ; --tmpLoop)
	{
		if (chromosome.compare(targetChromosome[tmpLoop]->chr_name) == 0)
			break;
	}

	if (tmpLoop < 0)
	{
		//chromosome not found
		outputfilename = dirPrefix + chromosome + ".txt";

		outputfile.push_back(new ofstream);		
		tmpLoop = outputfile.size() - 1;
		outputfile[tmpLoop]->open(outputfilename.c_str());

		new_chr_entry = new chr_entry;
		new_chr_entry->chr_name = chromosome;
		new_chr_entry->chr_length = fp->header->target_len[c->tid];
		targetChromosome.push_back(new_chr_entry);		
		readAlignIndex.push_back(0);
	}

	readAlignIndex[tmpLoop] += 1;
	
	//output read information
	*outputfile[tmpLoop] << c->pos << "\t" << get_cigar_str(bam1_cigar(b), c->n_cigar);
	
	//output XS field
	uint8_t *xs = bam_aux_get(b, "XS");
	if (xs)
	{
		*outputfile[tmpLoop] << "\t" << bam_aux2A(xs) << endl;
	} 
	else
	{
		*outputfile[tmpLoop] << "\t+" << endl;
	}
	
	return 0;
}


void outputStats()
{
	unsigned int tmpLoop;
	string outputfilename;
	ofstream outputChrname, outputDataStat;
	long totalAlignNum = 0;

	outputfilename = dirPrefix + "ChromosomeName.txt";
	outputChrname.open(outputfilename.c_str());
	outputfilename = dirPrefix + "DatasetStat.txt";
	outputDataStat.open(outputfilename.c_str(), fstream::app);

	for (tmpLoop = 0; tmpLoop < targetChromosome.size(); ++tmpLoop)
	{
		outputChrname << targetChromosome[tmpLoop]->chr_name << '\t' << readAlignIndex[tmpLoop] << endl; 
		totalAlignNum += readAlignIndex[tmpLoop];

		(*outputfile[tmpLoop]).close();
		delete outputfile[tmpLoop];
	}

	outputChrname << "#" << totalReadNum << '\t' << totalAlignNum << endl;

	outputDataStat << totalReadNum << '\t' << totalAlignNum << endl;

	outputChrname.close();
	outputDataStat.close();

	totalReadNum_persample.push_back(totalReadNum);
	totalAlignNum_persample.push_back(totalAlignNum);

	return;
}

void process_one_file(string filename)
{
	//initialization
	totalReadNum = 0;

	//process the given bam file
	samfile_t *fp;
	if ((fp = samopen(filename.c_str(), "rb", 0)) == 0) {
		fprintf(stderr, "parse_bam: Fail to open BAM file %s\n", filename.c_str());
		return;
	}
	bam1_t *b = bam_init1();
	while (samread(fp, b) >= 0) fetch_func(b, fp);
	bam_destroy1(b);
	samclose(fp);

	//output read stats
	outputStats();

	//add chromosomes not already cataloged
	unsigned int iLoop, jLoop;
	for (iLoop = 0; iLoop < targetChromosome.size(); ++iLoop)
	{
		for (jLoop = 0; jLoop < allChromosome.size(); ++jLoop)
		{
			if (targetChromosome[iLoop]->chr_name.compare(allChromosome[jLoop]->chr_name) == 0)
			{
				break;
			}
		}

		if (jLoop >= allChromosome.size())
		{
			//new chromosome
			allChromosome.push_back(targetChromosome[iLoop]);
			targetChromosome[iLoop] = NULL;
		}
		else
		{
			if (targetChromosome[iLoop]->chr_length > allChromosome[jLoop]->chr_length)
			{
				allChromosome[jLoop]->chr_length = targetChromosome[iLoop]->chr_length;
			}
		}
	}

	//clear vectors
	for (iLoop = 0; iLoop < targetChromosome.size(); ++iLoop)
	{
		if (!targetChromosome[iLoop])
			delete targetChromosome[iLoop];
	}
	targetChromosome.clear();
	readAlignIndex.clear();
	outputfile.clear();

	return;
}

int main(int argc, char *argv[])
{
	if (argc != 3) {
		fprintf(stderr, "Usage: parse_bam <file_bamfilenames> <target_path>\n");
		return 1;
	}

	ifstream infile(argv[1]);
	if (infile.is_open() == false)
	{
		cout << "Error: Fail to open the file containing the bam file list." << endl;
		return 1;
	}

	string fileindex, filename, info, parse_result_dir, cmd;

	parse_result_dir = argv[2];
	if (parse_result_dir[parse_result_dir.size()-1] != '/')
		parse_result_dir += "/";

	while (infile >> fileindex)
	{
		infile >> filename;
		getline(infile, info);

		fileindex_persample.push_back(fileindex);

		cmd = "mkdir -p " + parse_result_dir + fileindex;
		system(cmd.c_str());

		dirPrefix = parse_result_dir + fileindex + "/";
		process_one_file(filename);
	}

	infile.close();
	

	//output all chromosome names in the file list
	string outputfilename = parse_result_dir + "ChromosomeName.txt";
	ofstream outputChrname;
	outputChrname.open(outputfilename.c_str(), fstream::app);
	
	if (outputChrname.is_open() == false)
	{
		cout << "Error: Fail to write the chromosome names to file." << endl;
		return 1;
	}

	for (unsigned int tmpLoop = 0; tmpLoop < allChromosome.size(); ++tmpLoop)
	{
		outputChrname << allChromosome[tmpLoop]->chr_name << "\t" << allChromosome[tmpLoop]->chr_length << endl; 
		delete allChromosome[tmpLoop];
	}
	outputChrname.close();


	//output the number of reads of each sample
	outputfilename = parse_result_dir + "DataStat.txt";
	ofstream outputDataStatAllSample;
	outputDataStatAllSample.open(outputfilename.c_str(), fstream::app);

	if (outputDataStatAllSample.is_open() == false)
	{
		cout << "Error: Fail to write the sample stats to file." << endl;
		return 1;
	}

	for (unsigned int tmpLoop = 0; tmpLoop < totalReadNum_persample.size(); ++tmpLoop)
	{
		outputDataStatAllSample << fileindex_persample[tmpLoop] << "\t" << totalReadNum_persample[tmpLoop] << "\t" << totalAlignNum_persample[tmpLoop] << endl; 
	}
	outputDataStatAllSample.close();


	allChromosome.clear();
	totalReadNum_persample.clear();
	totalAlignNum_persample.clear();

	return 0;
}
