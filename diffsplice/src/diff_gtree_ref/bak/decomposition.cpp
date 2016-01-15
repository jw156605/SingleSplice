#include "decomposition.h"
#include "common_function.h"
#include "output.h"

//junction graph path queue for computing all paths
queue <JuncGraphPath*> juncPathQueue;

//stack for DFS on the g-tree
vector <GTvertex*> stack_GTvertex;

//graph for the search of paths
JuncGraph *junctionGraph;

long GENEcount = 0;
long ASMcount = 0;

ofstream outfile_not_enumerated;
//ofstream outfile_debug;
ofstream outfile_gene;

//RangeJunctionList Queue
//for separating dependent paths
void juncPathQueue_initialization()
{
	while (juncPathQueue.empty() == false)
	{
		juncPathQueue.pop();
	}

	return;
}

bool juncPathQueue_enqueue(JuncGraphPath *x)
{
	juncPathQueue.push(x);

	return true;
}

JuncGraphPath* juncPathQueue_dequeue()
{
	if (juncPathQueue.empty() == true)
	{
		return NULL;
	} 
	else
	{
		JuncGraphPath *next_ele = juncPathQueue.front();
		juncPathQueue.pop();
		return next_ele;
	}
}


/************************************************************************/
/* JUNCTION GRAPH                                                       */
/************************************************************************/


JuncGraphEdge::JuncGraphEdge()
{
	linkedVertex = NULL;
	next = NULL;
}

JuncGraphVertex::JuncGraphVertex()
{
	corresJunc = NULL;
	edges = NULL;
	next = NULL;
	traversed = false;
	hasInEdge = false;
	hasOutEdge = false;
	vertexType = normal;
}

JuncGraphVertex::~JuncGraphVertex()
{
	if (corresJunc != NULL)
	{
		delete corresJunc;
	}

	JuncGraphEdge *delEdge;

	while (edges != NULL)
	{
		delEdge = edges;
		edges = delEdge->next;
		delete delEdge;
	}
}

JuncGraphPath::JuncGraphPath()
{
	pathJuncList = NULL;
	arrivedVertex = NULL;
}

JuncGraphPath::~JuncGraphPath()
{
	arrivedVertex = NULL;
	if (pathJuncList != NULL)
		delete pathJuncList;
}

JuncGraph::JuncGraph()
{
	vertices = NULL;
}

JuncGraph::~JuncGraph()
{
	JuncGraphVertex *delVertex;

	while (vertices != NULL)
	{
		delVertex = vertices;
		vertices = delVertex->next;
		delete delVertex;
	}
}


//void output_vertex_debug(const GTvertex *out_vertex)
//{
//	outfile_debug << out_vertex->vertex_id << "\t" << chromosomeName << ":" << out_vertex->rangeLow << '-' << out_vertex->rangeHigh << endl;
//	rangeJunction *countJunc = out_vertex->junctionInRange->list;
//	while (countJunc != NULL)
//	{
//		if (countJunc->junc->type == frag_junction)
//		{
//			outfile_debug << "0\t";
//		}
//		else
//		{
//			outfile_debug << "1\t";
//		}
//		outfile_debug << countJunc->junc->start << "\t" << countJunc->junc->end << endl;
//
//		countJunc = countJunc->next;
//	}
//	outfile_debug << endl;
//}

//fragment graph
//assume junctions have been sorted based on position

//separate genes 
RangeJunctionList* separateGene(RangeJunctionList* origList, bool &separable)
{
	//separate genes
	//return a set of lists, each of which corresponds to a gene
	//IMPORTANT: RangeJunctionList is head-inserting, for the convenience of stack

	if (origList == NULL)
	{
		separable = false;
		return NULL;
	}

	RangeJunctionList *resultList, *curList;
	resultList = NULL;
	curList = NULL;

	rangeJunction *curJunc, *curListTail;
	curJunc = NULL;
	curListTail = NULL;

	long endBoard = -2;
	int numIndep = 0; //number of independent regions

	curJunc = origList->list;
	while (curJunc != NULL)
	{
		if (curJunc->junc->start > endBoard+1)
		{
			//get an independent region
			numIndep++;

			//omit the first list, it is empty
			if (curList == NULL)
			{
				//do nothing
			} 
			else
			{
				curList->rangeHigh = endBoard;
			}

			curList = new RangeJunctionList;
			curList->nextList = resultList;

			// 			if (resultList == NULL)
			// 			{
			// 				//first list
			// 				curList->rangeLow = origList->rangeLow;
			// 			} 
			// 			else
			// 			{
			// 				curList->rangeLow = endBoard;
			// 			}
			curList->rangeLow = curJunc->junc->start;

			resultList = curList;

			curListTail = NULL;
		}

		//process current fragment
		if (curListTail == NULL)
		{
			curList->list = curJunc;
			curListTail = curJunc;
		}
		else
		{
			curListTail->next = curJunc;
			curListTail = curJunc;
		}

		if (curJunc->junc->end > endBoard)
		{
			endBoard = curJunc->junc->end;
		}

		curJunc = curJunc->next;
		curListTail->next = NULL;
	}

	if (curList == NULL)
	{
		cout << "curList == NULL" << endl;
		exit(1);
	} 
	else
	{
		//curList->rangeHigh = origList->rangeHigh;
		curList->rangeHigh = endBoard;
	}

	if (numIndep == 1)
		separable = false; //cannot be separated into multiple regions
	else if (numIndep > 1)
		separable = true; //can be separated
	else
	{
		cout << "numGene == 0" << endl;
		exit(1);
	}

	origList->list = NULL;
	delete origList;

	return resultList;
}




//separate independent paths

bool compatibleJunctions(rangeJunction* junctionA, rangeJunction* junctionB)
{
	//check the compatibility of two junctions
	if (junctionA->junc->end <= junctionB->junc->start && junctionB->junc->start - junctionA->junc->end <= 1)
	{
		return true;
	}
	else if (junctionA->junc->start >= junctionB->junc->end && junctionA->junc->start - junctionB->junc->end <= 1)
	{
		return true;
	}
	else
		return false;
}

void constructJuncGraph_undirected(RangeJunctionList* origList, bool virtualSE)
{
	//construct fragment graph based on given junctions
	//add virtual start and virtual end if virtualSE is true 

	JuncGraphVertex *curVertex, *tailVertex, *edgeVertex;
	JuncGraphEdge *newEdge;
	rangeJunction *curJunc;

	//build vertices

	curJunc = origList->list;
	while (curJunc != NULL)
	{
		curVertex = new JuncGraphVertex;
		curVertex->corresJunc = curJunc;
		curJunc = curJunc->next;
		curVertex->corresJunc->next = NULL;

		if (junctionGraph->vertices == NULL)
		{
			junctionGraph->vertices = curVertex;
			tailVertex = curVertex;
		} 
		else
		{
			tailVertex->next = curVertex;
			tailVertex = curVertex;
		}
	}

	if (virtualSE == true)
	{
		curVertex = new JuncGraphVertex;
		curVertex->vertexType = virStart;
		curVertex->next = junctionGraph->vertices;
		junctionGraph->vertices = curVertex;

		curVertex = new JuncGraphVertex;
		curVertex->vertexType = virEnd;
		tailVertex->next = curVertex;
		tailVertex = curVertex;
	}

	//build edges

	curVertex = junctionGraph->vertices;
	while (curVertex != NULL)
	{
		if (curVertex->vertexType == normal)
		{
			edgeVertex = curVertex->next;
			while (edgeVertex != NULL)
			{
				if (edgeVertex->vertexType == normal)
				{
					if (compatibleJunctions((curVertex->corresJunc), (edgeVertex->corresJunc)) == true)
					{
						//create an edge for curVertex
						newEdge = new JuncGraphEdge;
						newEdge->linkedVertex = edgeVertex;

						newEdge->next = curVertex->edges;
						curVertex->edges = newEdge;

						//create an edge for edgeVertex
						newEdge = new JuncGraphEdge;
						newEdge->linkedVertex = curVertex;

						newEdge->next = edgeVertex->edges;
						edgeVertex->edges = newEdge;

						edgeVertex->hasInEdge = true;
						curVertex->hasOutEdge = true;
					} 
					else
					{
						//break; 
					}
				}

				edgeVertex = edgeVertex->next;
			}
		}

		curVertex = curVertex->next;
	}

	if (virtualSE == true)
	{
		curVertex = junctionGraph->vertices;
		edgeVertex = curVertex->next;
		while (edgeVertex != NULL)
		{
			if (edgeVertex->vertexType == normal && edgeVertex->hasInEdge == false)
			{
				//create an edge for curVertex
				newEdge = new JuncGraphEdge;
				newEdge->linkedVertex = edgeVertex;

				newEdge->next = curVertex->edges;
				curVertex->edges = newEdge;

				//create an edge for edgeVertex
				newEdge = new JuncGraphEdge;
				newEdge->linkedVertex = curVertex;

				newEdge->next = edgeVertex->edges;
				edgeVertex->edges = newEdge;

			}

			edgeVertex = edgeVertex->next;
		}

		curVertex = tailVertex;
		edgeVertex = junctionGraph->vertices->next;
		while (edgeVertex != NULL)
		{
			if (edgeVertex->vertexType == normal && edgeVertex->hasOutEdge == false)
			{
				//create an edge for curVertex
				newEdge = new JuncGraphEdge;
				newEdge->linkedVertex = edgeVertex;

				newEdge->next = curVertex->edges;
				curVertex->edges = newEdge;

				//create an edge for edgeVertex
				newEdge = new JuncGraphEdge;
				newEdge->linkedVertex = curVertex;

				newEdge->next = edgeVertex->edges;
				edgeVertex->edges = newEdge;

			}

			edgeVertex = edgeVertex->next;
		}
	}


	return;
}

void constructJuncGraph_directed(RangeJunctionList* origList)
{
	//construct fragment graph based on given junctions

	JuncGraphVertex *curVertex, *tailVertex, *edgeVertex;
	JuncGraphEdge *newEdge;
	rangeJunction *curJunc;

	//build vertices

	curJunc = origList->list;
	while (curJunc != NULL)
	{
		curVertex = new JuncGraphVertex;
		curVertex->corresJunc = curJunc;
		curJunc = curJunc->next;
		curVertex->corresJunc->next = NULL;

		if (junctionGraph->vertices == NULL)
		{
			junctionGraph->vertices = curVertex;
			tailVertex = curVertex;
		} 
		else
		{
			tailVertex->next = curVertex;
			tailVertex = curVertex;
		}
	}

	//build edges

	curVertex = junctionGraph->vertices;
	while (curVertex != NULL)
	{
		edgeVertex = curVertex->next;
		while (edgeVertex != NULL)
		{
			if (compatibleJunctions((curVertex->corresJunc), (edgeVertex->corresJunc)) == true)
			{
				//create an edge for curVertex
				newEdge = new JuncGraphEdge;
				newEdge->linkedVertex = edgeVertex;

				newEdge->next = curVertex->edges;
				curVertex->edges = newEdge;

				// 				//create an edge for edgeVertex
				// 				newEdge = new JuncGraphEdge;
				// 				newEdge->linkedVertex = curVertex;
				// 
				// 				newEdge->next = edgeVertex->edges;
				// 				edgeVertex->edges = newEdge;

				edgeVertex->hasInEdge = true;
				curVertex->hasOutEdge = true;
			} 
			else
			{
				//break; 
			}
			edgeVertex = edgeVertex->next;
		}

		curVertex = curVertex->next;
	}


	return;
}


void DFS_visit(JuncGraphVertex *u, vector <rangeJunction*> &list_visited)
{
	//visit a vertex u during DFS
	u->traversed = true;

	JuncGraphVertex *v;
	JuncGraphEdge *curEdge;
	curEdge = u->edges;
	while (curEdge != NULL)
	{
		v = curEdge->linkedVertex;
		if (v->traversed == false)
		{
			DFS_visit(v, list_visited);
		}
		curEdge = curEdge->next;
	}

	//add u into sortList_Junction
	if (list_visited.size() >= list_visited.capacity())
		list_visited.reserve(list_visited.size()+1000);
	list_visited.push_back(u->corresJunc);
	u->corresJunc = NULL;

	return;
}

RangeJunctionList* search_conn_comp(int &numCC, long rangeLow, long rangeHigh)
{
	//search all connected component within junctionGraph
	JuncGraphVertex *curVertex;
	RangeJunctionList *resultList, *curList;
	long tmpLoop;
	vector <rangeJunction*> list_visited;
	void** sortlist_feature;
	double* sortkey_feature;

	numCC = 0; //number of connected component

	resultList = NULL;
	curVertex = junctionGraph->vertices;
	while (curVertex != NULL)
	{
		if (curVertex->traversed == false)
		{
			DFS_visit(curVertex, list_visited);

			numCC++;

			curList = new RangeJunctionList;

			//get a connected component
			sortlist_feature = new void* [list_visited.size() + 2];
			sortkey_feature = new double [list_visited.size() + 2];
			for (tmpLoop = 0; tmpLoop < list_visited.size(); ++tmpLoop)
			{
				sortlist_feature[tmpLoop+1] = (void*) list_visited[tmpLoop];
				sortkey_feature[tmpLoop+1] = list_visited[tmpLoop]->junc->end;
			}
			mergeSort(sortlist_feature, sortkey_feature, list_visited.size());
			curList->rangeHigh = ((rangeJunction*)sortlist_feature[list_visited.size()])->junc->end;

			for (tmpLoop = 1; tmpLoop <= list_visited.size(); ++tmpLoop)
				sortkey_feature[tmpLoop] = ((rangeJunction*)sortlist_feature[tmpLoop])->junc->start;
			mergeSort(sortlist_feature, sortkey_feature, list_visited.size());
			curList->rangeLow = ((rangeJunction*)sortlist_feature[1])->junc->start;


			for (tmpLoop = 1; tmpLoop < list_visited.size(); ++tmpLoop)
				((rangeJunction*)sortlist_feature[tmpLoop])->next = (rangeJunction*) sortlist_feature[tmpLoop+1];
			((rangeJunction*)sortlist_feature[list_visited.size()])->next = NULL;
			curList->list = ((rangeJunction*)sortlist_feature[1]);
			curList->listtail = ((rangeJunction*)sortlist_feature[list_visited.size()]);

			delete [] sortlist_feature;
			delete [] sortkey_feature;

			curList->nextList = resultList;
			resultList = curList;
			
			list_visited.clear();
		}

		curVertex = curVertex->next;
	}

	return resultList;
}

//separate independent regions
RangeJunctionList* separateIndepRegion(RangeJunctionList* origList, bool &separable)
{
	//separate independent regions
	//return a set of lists, each of which corresponds to an independent region
	//IMPORTANT: RangeJunctionList is head-inserting, for the convenience of stack

	if (origList == NULL)
	{
		separable = false;
		return NULL;
	}

	junctionGraph = new JuncGraph;
	constructJuncGraph_undirected(origList, true);

	RangeJunctionList *resultList, *curList;
	resultList = NULL;
	curList = NULL;

	rangeJunction *curJunc, *curListTail;
	curJunc = NULL;
	curListTail = NULL;

	JuncGraphVertex *curVertex;
	JuncGraphEdge *curEdge;

	long endBoard = 0, geneEnd = 0;
	int numIndep = 0; //number of independent regions

	//handling alternative start
	curEdge = junctionGraph->vertices->edges;
	while (curEdge != NULL)
	{
		if (curEdge->linkedVertex->corresJunc->junc->start > endBoard)
		{
			endBoard = curEdge->linkedVertex->corresJunc->junc->start;
		}
		curEdge = curEdge->next;
	}

	curVertex = junctionGraph->vertices->next;
	while (curVertex != NULL)
	{
		if (curVertex->vertexType == normal)
		{
			curJunc = curVertex->corresJunc;
			curVertex->corresJunc = NULL;

			if (curJunc->junc->start >= endBoard)
			{
				//omit the first list, it is empty
				if (curList != NULL)
				{
					if (endBoard < MAX_CHR_LENGTH)
						curList->rangeHigh = endBoard;
					else
						curList->rangeHigh = geneEnd;
				}

				curList = NULL;
			}

			if (curList == NULL)
			{
				//get an independent region
				numIndep++;

				curList = new RangeJunctionList;
				curList->nextList = resultList;
				curList->rangeLow = curJunc->junc->start;
				resultList = curList;

				curListTail = NULL;
			}

			//process current fragment
			if (curListTail == NULL)
			{
				curList->list = curJunc;
				curListTail = curJunc;
			}
			else
			{
				curListTail->next = curJunc;
				curListTail = curJunc;
			}
			curListTail->next = NULL;


			if (curVertex->edges->linkedVertex->vertexType == virEnd)
			{
				endBoard = MAX_CHR_LENGTH;
				if (curJunc->junc->end > geneEnd)
				{
					geneEnd = curJunc->junc->end;
				}
			}
			else if (curJunc->junc->end > endBoard)
			{
				endBoard = curJunc->junc->end;
			}
		}

		curVertex = curVertex->next;
	}

	if (curList == NULL)
	{
		cout << "curList == NULL" << endl;
		exit(1);
	} 
	else
	{
		//curList->rangeHigh = origList->rangeHigh;
		if (endBoard < MAX_CHR_LENGTH)
			curList->rangeHigh = endBoard;
		else
			curList->rangeHigh = geneEnd;
	}


	delete junctionGraph;

	if (numIndep == 1)
		separable = false; //cannot be separated into multiple regions
	else if (numIndep > 1)
		separable = true; //can be separated
	else
	{
		cout << "numIndep == 0" << endl;
		exit(1);
	}

	origList->list = NULL;
	delete origList;

	return resultList;
}


RangeJunctionList* separateIndepPath(RangeJunctionList* origList, bool &separable)
{
	//separate independent paths 
	//return a set of lists, each of which corresponds to an independent path
	//return decomposable
	if (origList == NULL)
	{
		separable = false;
		return NULL;
	}

	junctionGraph = new JuncGraph;
	constructJuncGraph_undirected(origList, false);
	origList->list = NULL;

	RangeJunctionList *resultList;
	int numCC = 0;
	resultList = search_conn_comp(numCC, origList->rangeLow, origList->rangeHigh);

	delete junctionGraph;

	if (numCC == 1)
		separable = false; //cannot be separated into multiple paths
	else if (numCC > 1)
		separable = true; //can be separated
	else
	{
		cout << "numCC == 0";
		exit(1);
	}

	origList->list = NULL;
	delete origList;

	return resultList;
}

void output_region_not_enumerated(RangeJunctionList *output_list, int run_code)
{
	if (run_code <= 0)
		return;
	
	string reason;
	if (run_code == 1)
		reason = "gtf_generation";
	else
		reason = "decomposition";

	outfile_not_enumerated << reason << "\t" << output_list->rangeLow << '\t' << output_list->rangeHigh << endl;
	rangeJunction *countJunc = output_list->list;
	while (countJunc != NULL)
	{
		if (countJunc->junc->type == frag_junction)
		{
			outfile_not_enumerated << "j\t";
		} 
		else
		{
			outfile_not_enumerated << "e\t";
		}
		outfile_not_enumerated << countJunc->junc->start << "\t" << countJunc->junc->end;
		if (OUTPUT_PER_SAMPLE_EXPR) {
			for (int tmp = 0; tmp < SUPPORT_VECTOR_SIZE; tmp++)
			{
				outfile_not_enumerated << "\t" << countJunc->junc->support[tmp];
			}
		}
		outfile_not_enumerated << endl;

		countJunc = countJunc->next;
	}
	outfile_not_enumerated << endl;

	return;
}

// output gene information for clustering (LUSC project)
void output_gene_composition(RangeJunctionList *output_list, string geneID)
{
	rangeJunction *thisJunc = output_list->list;
	while (thisJunc != NULL)
	{
		outfile_gene << geneID << "\t" << output_list->rangeLow << '\t' << output_list->rangeHigh << "\t";

		if (thisJunc->junc->type == frag_junction)
		{
			outfile_gene << "j\t";
		}
		else
		{
			outfile_gene << "e\t";
		}
		outfile_gene << thisJunc->junc->start << "\t" << thisJunc->junc->end;
		for (int tmp = 0; tmp < SUPPORT_VECTOR_SIZE; tmp++)
		{
			outfile_gene << "\t" << thisJunc->junc->support[tmp];
		}
		outfile_gene << endl;

		thisJunc = thisJunc->next;
	}

	return;
}


//separate dependent paths and return a set of lists, each of which corresponds to a dependent path; if not (fully) enumerated, return NULL
RangeJunctionList* separateDepPath(RangeJunctionList* origList, bool &separable, int &numPath, int run_code)
{	
	if (origList == NULL)
	{
		separable = false;
		return origList;
	}

	separable = true;
	numPath = 0;

	rangeJunction *countJunc;
	long junctionCntinList = 0;
	countJunc = origList->list;
	while (countJunc != NULL)
	{
		if (countJunc->junc->type == frag_junction)
		{
			junctionCntinList++;
		}
		countJunc = countJunc->next;
	}
	if (junctionCntinList > MAXJUNCNUMINDEPPATH)
	{
		output_region_not_enumerated(origList, run_code);

		separable = false;
		return origList;
	}

	RangeJunctionList *backupList;
	backupList = origList->clone();
	junctionGraph = new JuncGraph;
	constructJuncGraph_directed(backupList);
	backupList->list = NULL;

	RangeJunctionList *resultList, *newList;
	rangeJunction *newJunc, *curJunc;
	long startPosition, endPosition;
	bool pathExtended;

	JuncGraphVertex *curVertex;
	JuncGraphEdge *curEdge;
	JuncGraphPath *curPath, *newPath;
	//find start vertex and end vertex

	startPosition = origList->rangeLow;
	endPosition = origList->rangeHigh;

	resultList = NULL;
	juncPathQueue_initialization();


	//build initial queue
	curVertex = junctionGraph->vertices;
	while (curVertex != NULL)
	{
		//		if (curVertex->corresJunc->junc->start == startPosition)
		if (curVertex->hasInEdge == false)
		{
			newList = new RangeJunctionList;
			newList->rangeLow = curVertex->corresJunc->junc->start;
			newList->rangeHigh = curVertex->corresJunc->junc->end;
			newList->transDirection = curVertex->corresJunc->junc->transDirection;

			newJunc = new rangeJunction;
			newJunc->junc = curVertex->corresJunc->junc;
			newList->list = newJunc;
			newList->listtail = newJunc;

			newPath = new JuncGraphPath;
			newPath->arrivedVertex = curVertex;
			newPath->pathJuncList = newList;

			juncPathQueue_enqueue(newPath);
		}
		curVertex = curVertex->next;
	}

	//main part
	curPath = juncPathQueue_dequeue();
	while (curPath != NULL)
	{
		if (juncPathQueue.size() > 10 * MAX_NUM_ENUMERATED_PATH || numPath > MAX_NUM_ENUMERATED_PATH)
		{
			separable = false;
			while (curPath = juncPathQueue_dequeue())
			{
				delete curPath;
			}
			output_region_not_enumerated(origList, run_code);
			
			break;
		}
		
		//		if (curPath->arrivedVertex->corresJunc->junc->end == endPosition)
		pathExtended = false;

		if (curPath->arrivedVertex->edges != NULL)
		{
			curEdge = curPath->arrivedVertex->edges;
			while (curEdge != NULL)
			{
				if (curEdge->linkedVertex->corresJunc->junc->end >= curPath->arrivedVertex->corresJunc->junc->end)
				{
					if (curPath->pathJuncList->transDirection == undetermined || curEdge->linkedVertex->corresJunc->junc->transDirection == undetermined
						|| curEdge->linkedVertex->corresJunc->junc->transDirection == curPath->pathJuncList->transDirection)
					{
						pathExtended = true;
						newList = curPath->pathJuncList->clone();

						newJunc = new rangeJunction;
						newJunc->junc = curEdge->linkedVertex->corresJunc->junc;
						//newJunc->next = newList->list;
						//newList->list = newJunc;
						newList->listtail->next = newJunc;
						newList->listtail = newJunc;

						if (newList->transDirection == undetermined && (newJunc->junc->transDirection == sense || newJunc->junc->transDirection == antisense))
						{
							newList->transDirection = newJunc->junc->transDirection;
						}

						if (newJunc->junc->start < newList->rangeLow)
						{
							newList->rangeLow = newJunc->junc->start;
						}
						if (newJunc->junc->end > newList->rangeHigh)
						{
							newList->rangeHigh = newJunc->junc->end;
						}

						newPath = new JuncGraphPath;
						newPath->arrivedVertex = curEdge->linkedVertex;
						newPath->pathJuncList = newList;

						juncPathQueue_enqueue(newPath);
					}
				}

				curEdge = curEdge->next;
			}			
		}

		if (pathExtended == false)
		{
			//current path arrives destination
			numPath++;

			curPath->pathJuncList->nextList = resultList;
			resultList = curPath->pathJuncList;		
			
			curPath->pathJuncList = NULL;
		}

		delete curPath;

		curPath = juncPathQueue_dequeue();
	}


	delete junctionGraph;

	if (separable == true)
	{
		if (numPath == 1)
			separable = false; //cannot be separated into multiple paths
		else if (numPath > 1)
			separable = true; //can be separated
		else
		{
			cout << origList->rangeLow << " - " << origList->rangeHigh << "  numPath == 0" << endl;
			exit(1);
		}
	}	

	delete backupList;

	if (separable == true)
	{
		return resultList;
	} 
	else
	{
		numPath = 1;
		delete resultList;
		return origList;
	}
}


/************************************************************************/
/* Vertex Stack for DFS                                                 */
/************************************************************************/
void stack_initial()
{
	stack_GTvertex.clear();
	stack_GTvertex.reserve(default_max_junction_num);

	return;
}

bool stack_empty()
{
	return stack_GTvertex.empty();
}

void stack_push(GTvertex *x)
{
	if (stack_GTvertex.size() >= stack_GTvertex.capacity())
		stack_GTvertex.reserve(stack_GTvertex.size() + default_max_junction_num);

	stack_GTvertex.push_back(x);

	return;
}

GTvertex* stack_pop()
{
	if (stack_empty() == true)
	{
		return NULL;
	} 
	else
	{
		GTvertex *vertex_pop = stack_GTvertex[stack_GTvertex.size()-1];
		stack_GTvertex.pop_back();
		return vertex_pop;
	}
}


bool calculateGeneMeanCoverage(RangeJunctionList *fragmentList, double *meanCoverageList)
{
	rangeJunction *curJunc;

	//count number of exons in the gene
	long exonCnt = 0, exonCntLoop;
	int sampleLoop;

	curJunc = fragmentList->list;
	while (curJunc != NULL)
	{
		if (curJunc->junc->type == frag_exon) // 5/15/2013 do not count retained intron || curJunc->junc->type == frag_retained_intron)
			++exonCnt;
		curJunc = curJunc->next;
	}

	if (exonCnt < 5)
		return false; //no need to do filtering

	//get expression array
	double **coverageArray = new double* [SUPPORT_VECTOR_SIZE];
	for (sampleLoop = 0; sampleLoop < SUPPORT_VECTOR_SIZE; ++sampleLoop)
	{
		coverageArray[sampleLoop] = new double [exonCnt + 1];

		curJunc = fragmentList->list;
		exonCntLoop = 0;
		while (curJunc != NULL)
		{
			if (curJunc->junc->type == frag_exon) // || curJunc->junc->type == frag_retained_intron)
				coverageArray[sampleLoop][++exonCntLoop] = curJunc->junc->support[sampleLoop];
			curJunc = curJunc->next;
		}
	}
	long *quicksortStack = new long [2*exonCnt + 10];

	//calculate mean coverage in each sample
	//5/15/2013 use upper 3 quartiles, not just 25-75 quantile, i.e., keep highest coverage also
	long index_q25, index_q75, index_q100;
	double total_coverage;
	index_q25 = long(floor(0.25 * (exonCnt - 1) + 1));
	index_q75 = long(floor(0.75 * (exonCnt - 1) + 1));
	index_q100 = exonCnt;

	for (sampleLoop = 0; sampleLoop < SUPPORT_VECTOR_SIZE; ++sampleLoop)
	{
		quicksort(coverageArray[sampleLoop], exonCnt, quicksortStack);

		total_coverage = 0.0;
		for (exonCntLoop = index_q25 + 1; exonCntLoop <= index_q100; ++exonCntLoop)
			total_coverage += coverageArray[sampleLoop][exonCntLoop];

		meanCoverageList[sampleLoop] = total_coverage / (index_q100 - index_q25);
	}

	//clean up
	for (sampleLoop = 0; sampleLoop < SUPPORT_VECTOR_SIZE; ++sampleLoop)
	{
		delete [] coverageArray[sampleLoop];
	}
	delete [] coverageArray;
	delete [] quicksortStack;

	return true;
}




void filterGeneLowCoverageExon(GTvertex *targetVertex)
{
	double *meanCoverageList = new double [SUPPORT_VECTOR_SIZE];
	rangeJunction *curJunc, *prevJunc;
	long iLoop;
	int thresh_lowCovSampleCnt = int(ceil(SUPPORT_VECTOR_SIZE * (1 - coverageThreshold_GeneExon))), lowCovSampleCnt, sampleLoop;

	//find mean coverage of 25% to 75% percentile expression in this gene
	if (calculateGeneMeanCoverage(targetVertex->junctionInRange, meanCoverageList) == false)
	{
		delete [] meanCoverageList;
		return;
	}

	for (sampleLoop = 0; sampleLoop < SUPPORT_VECTOR_SIZE; ++sampleLoop)
	{
		targetVertex->proportion[sampleLoop] = meanCoverageList[sampleLoop];
	}

	//do not run this for now
	delete [] meanCoverageList;

	return;
}


//construct the tree hierarchy with the given feature list and decompose , return the root
GTvertex* constructGTree(RangeJunctionList* origList)
{
	//construct GTree
	//input: sorted range fragment list
	bool separable_region, separable_path, separable_depPath;
	RangeJunctionList *curList, *curPathList, *gtfPathList, *delList;
	long asm_path_cnt;
	GTvertex *newVertex;
	GTedge *newEdge, *curEdge;
	int num_enumerated_path;

	stack_initial();

	GTvertex *curVertex, *lastVertex, *rootVertex;
	curVertex = new GTvertex;
	curVertex->junctionInRange = origList;
	curVertex->level = 0;
	curVertex->rangeLow = origList->rangeLow; // CHROMOSOME_START;
	curVertex->rangeHigh = origList->rangeHigh; // CHROMOSOME_END;
	curVertex->child = NULL;
	curVertex->childType = 1;

	rootVertex = curVertex;

	//separate genes
	separable_region = true;
	curList = NULL;

	curList = separateGene(curVertex->junctionInRange, separable_region);

	if (separable_region == true)
	{
		while (curList != NULL)
		{
			//(1) if two genes are overlapping, e.g., positive strand and negative strand genes on the same locus
			//(2) if there are some standalone stuff, like false region, microRNA, etc.   they are not paths, but may encode some important information
			//therefore, separate them first.
			curPathList = curList;
			curList = curList->nextList;
			curPathList->nextList = NULL;
			separable_path = true;

			curPathList = separateIndepPath(curPathList, separable_path);

			if (separable_path == true)
			{
				while (curPathList != NULL)
				{
					newVertex = new GTvertex;
					newVertex->level = curVertex->level; //genes, same level as the chromosome, 0
					newVertex->junctionInRange = curPathList;
					newVertex->rangeLow = curPathList->rangeLow;
					newVertex->rangeHigh = curPathList->rangeHigh;
					newVertex->child = NULL;
					newVertex->childType = 1;

					curPathList = curPathList->nextList;
					newVertex->junctionInRange->nextList = NULL;

					newEdge = new GTedge;
					newEdge->linkedVertex = newVertex;
					newEdge->next = curVertex->child;
					curVertex->child = newEdge;

					filterGeneLowCoverageExon(newVertex);
					stack_push(newVertex);
				}
			}
			else
			{
				newVertex = new GTvertex;
				newVertex->level = curVertex->level; //genes, same level as the chromosome, 0
				newVertex->junctionInRange = curPathList;
				newVertex->rangeLow = curPathList->rangeLow;
				newVertex->rangeHigh = curPathList->rangeHigh;
				newVertex->child = NULL;
				newVertex->childType = 1;

				newEdge = new GTedge;
				newEdge->linkedVertex = newVertex;
				newEdge->next = curVertex->child;
				curVertex->child = newEdge;

				filterGeneLowCoverageExon(newVertex);
				stack_push(newVertex);
			}			
		}
	} 
	else
	{
		curVertex->junctionInRange = curList;
		curVertex->childType = 1;

		//(1) if two genes are overlapping, e.g., positive strand and negative strand genes on the same locus
		//(2) if there are some standalone stuff, like false region, microRNA, etc.   they are not paths, but may encode some important information
		//therefore, separate them first.
		separable_path = true;
		curList = NULL;

		curList = separateIndepPath(curVertex->junctionInRange, separable_path);

		if (separable_path == true)
		{
			while (curList != NULL)
			{
				newVertex = new GTvertex;
				newVertex->level = curVertex->level; //genes, same level as the chromosome, 0
				newVertex->junctionInRange = curList;
				newVertex->rangeLow = curList->rangeLow;
				newVertex->rangeHigh = curList->rangeHigh;
				newVertex->child = NULL;
				newVertex->childType = 1;

				curList = curList->nextList;
				newVertex->junctionInRange->nextList = NULL;

				newEdge = new GTedge;
				newEdge->linkedVertex = newVertex;
				newEdge->next = curVertex->child;
				curVertex->child = newEdge;

				filterGeneLowCoverageExon(newVertex);
				stack_push(newVertex);
			}
		}
		else
		{
			curVertex->junctionInRange = curList;
			filterGeneLowCoverageExon(curVertex);
			stack_push(curVertex);
		}
	}


	//separate ASMs
	curVertex = stack_pop();
	while (curVertex != NULL)
	{
// 		if (curVertex->rangeLow == 24563643)// && curVertex->rangeHigh == 28211806)
// 		{
//			output_vertex_debug(curVertex);
// 			cout << "a";
// 		}

		separable_region = true;
		separable_path = true;
		separable_depPath = true;

		lastVertex = NULL;
		curList = NULL;

		if (curVertex->level < 1 && curVertex->vertex_id.empty())
		{
			curVertex->vertex_id = "gene" + itostr(++GENEcount);

			// output gene composition information for clustering (LUSC project)
			output_gene_composition(curVertex->junctionInRange, curVertex->vertex_id);
		}

		if (curVertex->childType == 1)
		{
			//independent regions
			//			cout << curVertex->rangeLow << "\t" << curVertex->rangeHigh << endl;
			//			if (curVertex->rangeLow == 50569593)
			//				cout << "catch --";
			curList = separateIndepRegion(curVertex->junctionInRange, separable_region);
			//			cout << "done" << endl;

			if (separable_region == true)
			{
				while (curList != NULL)
				{
					newVertex = new GTvertex;
					newVertex->level = curVertex->level + 1;
					newVertex->junctionInRange = curList;
					newVertex->rangeLow = curList->rangeLow;
					newVertex->rangeHigh = curList->rangeHigh;
					newVertex->childType = 2;

					curList = curList->nextList;
					newVertex->junctionInRange->nextList = NULL;

					newEdge = new GTedge;
					newEdge->linkedVertex = newVertex;
					newEdge->next = curVertex->child;
					curVertex->child = newEdge;

					if (lastVertex != NULL)
					{
						lastVertex->prevSibling = newVertex;
						newVertex->nextSibling = lastVertex;
					}
					lastVertex = newVertex;

					stack_push(newVertex);
				}
			} 
			else
			{
				curVertex->junctionInRange = curList;
				curVertex->childType = 2;

				stack_push(curVertex);
			}
		} 
		else if (curVertex->childType == 2)
		{
			//independent paths
			if (curVertex->vertex_id.empty() && curVertex->junctionInRange->list->next != NULL)
				curVertex->vertex_id = "asm" + itostr(++ASMcount);

// 			if (curVertex->vertex_id.compare("asm254") == 0 || curVertex->vertex_id.compare("asm541") ==0 || curVertex->vertex_id.compare("asm821") == 0)
// 			{
// 				output_vertex_debug(curVertex);
// 			}
			

			/////////////////////////////////////////////////////////////////////////////
			//output gtf tracks for level 1 ASMs 
			//the 2nd condition check whether there are multiple fragments, i.e., an ASM
			if (curVertex->level <= 1 && curVertex->junctionInRange->list->next != NULL)
			{		
				gtfPathList = separateDepPath(curVertex->junctionInRange, separable_depPath, num_enumerated_path, 1);
				if (separable_depPath == true)
				{
					output_ASMpath_gtf(curVertex, gtfPathList);
					while (gtfPathList != NULL)
					{
						delList = gtfPathList;
						gtfPathList = gtfPathList->nextList;
						delete delList;
					}
				}
				else
				{
					output_ASMpath_gtf(curVertex, gtfPathList);
					curVertex->junctionInRange = gtfPathList;
				}
			}
			/////////////////////////////////////////////////////////////////////////////

			//4/7/2013 try to separate dependent path first, if the number of paths is small, take them
			curList = separateDepPath(curVertex->junctionInRange, separable_depPath, num_enumerated_path, 0);
			curVertex->path_extended = true;
			if (separable_depPath == true && num_enumerated_path <= MAX_NUM_PATH_FOR_NO_DECOMPOSITION)
			{
				//reach a leaf
				curVertex->childType = 3;
				
				while (curList != NULL)
				{
					newVertex = new GTvertex;
					newVertex->level = curVertex->level + 1;
					newVertex->junctionInRange = curList;
					newVertex->rangeLow = curList->rangeLow;
					newVertex->rangeHigh = curList->rangeHigh;
					newVertex->childType = 0;

					curList = curList->nextList;
					newVertex->junctionInRange->nextList = NULL;

					newEdge = new GTedge;
					newEdge->linkedVertex = newVertex;
					newEdge->next = curVertex->child;
					curVertex->child = newEdge;
				}

				//assign id to alternative paths
				asm_path_cnt = 0;
				curEdge = curVertex->child;
				while (curEdge != NULL)
				{
					curEdge->linkedVertex->vertex_id = curVertex->vertex_id + ".p" + itostr(++asm_path_cnt);
					curEdge = curEdge->next;
				}
			} 
			else
			{
				if (separable_depPath == true)
				{
					while (curList != NULL)
					{
						delList = curList;
						curList = curList->nextList;
						delete delList;
					}
				}

				curList = separateIndepPath(curVertex->junctionInRange, separable_path);

				if (separable_path == true)
				{
					while (curList != NULL)
					{
						newVertex = new GTvertex;
						newVertex->level = curVertex->level + 1;
						newVertex->junctionInRange = curList;
						newVertex->rangeLow = curList->rangeLow;
						newVertex->rangeHigh = curList->rangeHigh;
						newVertex->childType = 1;

						curList = curList->nextList;
						newVertex->junctionInRange->nextList = NULL;

						newEdge = new GTedge;
						newEdge->linkedVertex = newVertex;
						newEdge->next = curVertex->child;
						curVertex->child = newEdge;

						// 					if (lastVertex != NULL)
						// 					{
						// 						lastVertex->prevSibling = newVertex;
						// 						newVertex->nextSibling = lastVertex;
						// 					}
						// 					lastVertex = newVertex;

						stack_push(newVertex);
					}

					//assign id to alternative paths
					asm_path_cnt = 0;
					curEdge = curVertex->child;
					while (curEdge != NULL)
					{
						curEdge->linkedVertex->vertex_id = curVertex->vertex_id + ".p" + itostr(++asm_path_cnt);
						curEdge = curEdge->next;
					}
				} 
				else
				{
					//reach a leaf
					curVertex->junctionInRange = curList;
					curList = separateDepPath(curList, separable_depPath, num_enumerated_path, 2);

					if (separable_depPath == true)
					{
						curVertex->childType = 3;

						while (curList != NULL)
						{
							newVertex = new GTvertex;
							newVertex->level = curVertex->level + 1;
							newVertex->junctionInRange = curList;
							newVertex->rangeLow = curList->rangeLow;
							newVertex->rangeHigh = curList->rangeHigh;
							newVertex->childType = 0;

							curList = curList->nextList;
							newVertex->junctionInRange->nextList = NULL;

							newEdge = new GTedge;
							newEdge->linkedVertex = newVertex;
							newEdge->next = curVertex->child;
							curVertex->child = newEdge;

							// 						if (lastVertex != NULL)
							// 						{
							// 							lastVertex->prevSibling = newVertex;
							// 							newVertex->nextSibling = lastVertex;
							// 						}
							// 						lastVertex = newVertex;
						}

						//assign id to alternative paths
						asm_path_cnt = 0;
						curEdge = curVertex->child;
						while (curEdge != NULL)
						{
							curEdge->linkedVertex->vertex_id = curVertex->vertex_id + ".p" + itostr(++asm_path_cnt);
							curEdge = curEdge->next;
						}
					} 
					else
					{
						curVertex->junctionInRange = curList;
						curVertex->childType = 0;
					}
				}
			}			
		}

		curVertex = stack_pop();
	}

	return rootVertex;
}









//construct the ASM hierarchy from the given list of features (exons and junctions), return the root of the tree
GTvertex* decomposition(RangeJunctionList *original_list)
{
	//build tree
	//cout << "total exon and junction = " << sortList_Junction_Num << "... " << flush;

	string filename = resultPath + "detail/" + chromosomeName + "_not_enumerated.txt";
	outfile_not_enumerated.open(filename.c_str());
//	outfile_debug.open("debug.txt");
	filename = resultPath + "gene/" + chromosomeName + "_gene.txt";
	outfile_gene.open(filename.c_str());

	GTvertex *gTree_root = constructGTree(original_list);
	outfile_not_enumerated.close();
//	outfile_debug.close();
	outfile_gene.close();

	return gTree_root;
}

