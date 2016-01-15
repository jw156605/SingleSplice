

/************************************************************************/
/* Define the classes for splice graph and ASM                          */
/************************************************************************/

#ifndef SPLICE_GRAPH
#define SPLICE_GRAPH

#include "config.h"

/************************************************************************/
/* CLASS DECLARATION	                                                */
/************************************************************************/

class GTvertex;
class fragment;

class JuncGraphVertex;

/************************************************************************/
/* TYPE DEFINITION	                                                    */
/************************************************************************/

enum fragment_type {frag_exon, frag_intron, frag_retained_intron, frag_junction, frag_insertion, frag_deletion, frag_SNP};
enum trans_direction {undetermined, sense, antisense, terminal, mixed}; // terminal is only useful for strand of exon boundary, it stands for the start/end of a gene, a boundary not defined by a junction. this is useful when trying to cut boundary for transcription start/end exons
enum JuncGraphVertexType {virStart, virEnd, normal};

enum exon_bound_type {na_exon_bound, boundary, in_exon, out_exon}; // location of a splice site as compared to the annotated exons (only when reference exons are provided and input_data_ref is called)
enum splicesite_type {splice_out, splice_in, splice_both}; // directionOut of splice site
enum exon_ref_type {na_exon_ref, annotated, extension, novel}; // for frag_exon, label the reference category for the exon

/************************************************************************/
/* CLASS DEFINITION - FRAGMENT                                          */
/************************************************************************/
class alter_junction
{
public:
	fragment *juncInfo;
	//	int category; //0 for the merged target junction, -1 for 5', 1 for 3', 2 for exon
	//	double proportion;
	alter_junction *next;

	alter_junction();
};

class fragment
{
public:
	string frag_name;
	unsigned long frag_id; //unique ID throughout the pipeline
	string chromosome_start;
	string chromosome_end;
	trans_direction transDirection;
	long start;
	long end;
	fragment_type type; //fragment type
	int altersite; //if the fragment is an exon, then this field indicates whether this fragment is an alternative splice site: 0 for exon, 1 for donor site, 2 for acceptor site
	exon_ref_type in_exon_annotation;
	bool in_junc_annotation;
	exon_bound_type bound_start; // if this fragment is a junction, then the ref exon-matching type of the starting coordinate
	exon_bound_type bound_end;
	
	// for ESG visualization
	unsigned long vis_exonid; // exon id for visualization: 1. this id is linear, smaller id means smaller coordinate; 2. adjacent id's mean exons with no gap (intron), such as retained intron or alternative splice sites
	fragment* splice_exonhead; // head exon of a splice
	fragment* splice_exontail; // tail exon of a splice

	trans_direction exon_strand_start;
	trans_direction exon_strand_end;

	double *support; //support from the fragment file

	alter_junction* alter;
	long alterFragCnt; //fragment count of alternative splice sites, note it's the number of fragments
	//	alter_junction* fivePalterTail; //tail pointer of 5' alternative junction
	long *coverage;
	long start_real; //real start position in case the position is changed due to alternative splice sites
	long end_real; //real end position

	fragment(int support_vec_size);
	fragment* clone(int support_vec_size); //make a clone, used for exon, so alter_junction not included
	~fragment();
};

class spliceSite
{
public:
	string chromosome;
	long position;
	bool directionOut; //true for out-going junction, false for in-coming junction
	trans_direction strand;
	bool is_ref;

	spliceSite();
	spliceSite(string chr, long pos, bool dirOut, trans_direction str, bool in_ref) : chromosome(chr), position(pos), directionOut(dirOut), strand(str), is_ref(in_ref) {}
};


/************************************************************************/
/* CLASS DEFINITION - GenomeTree                                        */
/************************************************************************/
class rangeJunction
{
	//fragment within a range
public:
	fragment *junc;
	rangeJunction *next;

	rangeJunction();
};

class RangeJunctionList
{
	//fragment list within a range
public:
	long rangeLow;
	long rangeHigh;
	unsigned long cnt_rangeJunction;
	trans_direction transDirection;
	rangeJunction *list;
	rangeJunction *listtail;
	RangeJunctionList *nextList;

	RangeJunctionList();
	RangeJunctionList* clone();
	void insert_feature(rangeJunction *new_rangeJunction);
	void count_featurelist();
	~RangeJunctionList();
};


class GTedge
{
	//edge of Genome Tree
public:
	GTvertex *linkedVertex;
	GTedge *next;

	GTedge();
};

class alternative_path
{
public:
	long path_start; //site that starts diverging
	long path_end; //site that ends diverging
	long whole_path_start; //the whole path including the starting exon
	long whole_path_end; //the whole path including the ending exon
	trans_direction transDirection;
	double *support; //support at this vertex
	double *proportion; //proportion at this vertex
	long junctionNum;
	long exonNum;
	GTvertex *pathVertex;
	alternative_path *next;

	double avg_support;

	alternative_path();
	~alternative_path();
};

class GTvertex
{
	//vertex of Genome Tree
public:
	string vertex_id;
	long level;
	long rangeLow;
	long rangeHigh;
	GTedge *child;
	int childType; //0 for no child (i.e. leaf node), 1 for independent region, 2 for independent path, 3 for dependent path
	long childNum;
	RangeJunctionList *junctionInRange; //junctions within the vertex's range
	long junctionNum;
	long exonNum;
	GTvertex *prevSibling; //the sibling on the left (5') side
	GTvertex *nextSibling; //the sibling on the right (3') side
	bool has_novel_junction;
	bool has_novel_junction_major;

	GTedge *alterSpliceSite; //alternative splice sites involved in this vertex

	//for difference analysis
	double *support; //support at this vertex
	double *proportion; //proportion at this vertex (for genes, this array stores mean 25%-75% coverage)
	double *MSE_estimation; //mean squared error of estimated transcript abundance
	double *min_path_support;
	double *obs_support; //observed support directly calculated from raw support (with no estimation)

	bool estimated; //if true, the vertex has been estimated and will be treated as a whole
	bool path_extended; //if true, this vertex has been thrown into separate path processes for the extension of 
	fragment *representative;
	long estimate_exonNum; //number of blocks under estimation, equivalent to number of exons previously

	int ASMcategory;

	alternative_path* major_alter_paths;
	long major_alter_paths_num;

	GTvertex();
	~GTvertex();
};

// class GenomeTree
// {
// public:
// 	GTvertex* root;
// 
// 	GenomeTree();
// 	~GenomeTree();
// };


bool comp_frag_ptr(const fragment *, const fragment *);
bool comp_ssite_withdup(const spliceSite &, const spliceSite &);

const long default_max_junction_num = 30000;


#endif

