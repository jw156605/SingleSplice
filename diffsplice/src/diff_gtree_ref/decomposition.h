
#ifndef DECOMPOSITION
#define DECOMPOSITION

#include "config.h"
#include "splice_graph.h"

/************************************************************************/
/* CLASS DEFINITION - JUNCTION GRAPH		                            */
/************************************************************************/

class JuncGraphEdge
{
	//edge of fragment graph
public:
	JuncGraphVertex *linkedVertex;
	JuncGraphEdge *next;

	JuncGraphEdge();
};

class JuncGraphVertex
{
	//vertex of fragment graph
public:
	rangeJunction* corresJunc;
	JuncGraphEdge* edges;
	JuncGraphVertex* next;
	bool traversed;
	bool hasInEdge; //if hasInEdge is false, consider it as a start vertex when enumerating transcripts
	bool hasOutEdge;

	JuncGraphVertexType vertexType;

	JuncGraphVertex();
	~JuncGraphVertex();
};

class JuncGraphPath
{
	//path of fragment graph
public:
	RangeJunctionList *pathJuncList;
	JuncGraphVertex *arrivedVertex;

	JuncGraphPath();
	~JuncGraphPath();
};

class JuncGraph
{
	//fragment graph
public:
	JuncGraphVertex *vertices;

	JuncGraph();
	~JuncGraph();
};


//class for filtering
// class deletedSites
// {
// public:
// 	long sites;
// 	deletedSites *next;
// 
// 	deletedSites();
// };

GTvertex* decomposition(RangeJunctionList *original_list);
RangeJunctionList* separateDepPath(RangeJunctionList* origList, bool &separable, int &numPath, int run_code); //run_code: 0 = no output, 1 = gtf generation, 2 = decomposition

#endif


