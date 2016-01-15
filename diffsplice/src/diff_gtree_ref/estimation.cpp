#include "config.h"
#include "splice_graph.h"
#include "decomposition.h"
#include "output.h"
#include "estimation.h"

bool countGTree(GTvertex *rootVertex);


double logFactorial(long n)
{
	//calculate the log of the factorial of n
	double logfact = 0.0;

	for (long iLoop = 2; iLoop <= n; iLoop++)
	{
		logfact += log(double(iLoop));
	}

	return logfact;
}

long count_estimate_juncNum(GTvertex *targetVertex)
{
	long juncCnt = 0;
	rangeJunction *curJunc;
	GTedge *curEdge;
	GTvertex *curVertex;

	if (targetVertex->childType == 3)
	{
		curJunc = targetVertex->junctionInRange->list;
		while (curJunc != NULL)
		{
			if (curJunc->junc->type == frag_junction)
			{
				++juncCnt;
			}
			curJunc = curJunc->next;
		}
	}
	else if (targetVertex->childType == 2)
	{
		curEdge = targetVertex->child;
		while (curEdge != NULL)
		{
			curVertex = curEdge->linkedVertex;
			if (curVertex->childType == 0 && curVertex->exonNum == 0 && curVertex->junctionNum == 1)
			{
				++juncCnt;
			} 

			curEdge = curEdge->next;
		}
	}

	return juncCnt;
}

double abundance_estimation(GTvertex *rootVertex, int index_tissue, double &MSE)
{
	//return the log likelihood
	//8/16/2011 MSE measure the mean squared error of the estimated coverage versus the observed coverage in the alternatively spliced region (i.e., excluding exons adjacent to the ASM)
	//9/12/2011 New estimation method

	if (rootVertex->child == NULL || rootVertex->childType == 1 || rootVertex->childType == 0)
	{
		return 0.0;
	}
	
	//transcript abundance estimation
	long exonCnt = 0, junctionCnt = 0, fragmentCnt_est = 0, transcriptCnt = 0;
	GTvertex *curVertex, *siblingVertex;
	GTedge *curEdge;
	rangeJunction *curJunc;
	long iLoop, jLoop, kLoop;
	bool hasPrevExon = false, hasNextExon = false;
	const int MaxPrevFragNum = 3, MaxNextFragNum = 3; // up to 3 exons ahead / after
	int prevFragCnt = 0, nextFragCnt = 0, readlength = 50, fakeJunctionLength = 50;

	curEdge = rootVertex->child;
	while (curEdge != NULL)
	{
		if (curEdge->linkedVertex->childType == 0 && curEdge->linkedVertex->junctionNum > MAXJUNCNUMINDEPPATH)
		{
			return 0.0;
		}
		curEdge = curEdge->next;
	}

	if (rootVertex->estimate_exonNum > 0)
	{
		exonCnt = rootVertex->estimate_exonNum;
	} 
	else
	{
		exonCnt = rootVertex->exonNum;
	}
	junctionCnt = count_estimate_juncNum(rootVertex);
	fragmentCnt_est = exonCnt + junctionCnt;
	transcriptCnt = rootVertex->childNum;

	bool **A = new bool * [transcriptCnt]; //matrix A, transcript-exon matrix
	fragment **fragArray = new fragment * [fragmentCnt_est + MaxPrevFragNum + MaxNextFragNum]; //exon array. exonArray[0] is the exon prev to rootVertex, and exonArray[1+exonCnt] is the exon after the rootVertex.
	long *transFragCnt = new long [transcriptCnt]; //frag count for every transcript
	double *transLength = new double [transcriptCnt]; //transcript length
	for (iLoop = 0; iLoop < transcriptCnt; ++iLoop)
	{
		A[iLoop] = new bool [fragmentCnt_est + MaxPrevFragNum + MaxNextFragNum];
		transFragCnt[iLoop] = 0;
		transLength[iLoop] = 0.0;
	}

	double *errorRatio_transExpr = new double [transcriptCnt];
	long *errorRatio_fragCnt = new long [transcriptCnt];
	for (iLoop = 0; iLoop < transcriptCnt; ++iLoop)
	{
		errorRatio_transExpr[iLoop] = 0.0;
		errorRatio_fragCnt[iLoop] = 0;
	}


	//build exonArray
	siblingVertex = rootVertex;
	while (prevFragCnt < MaxPrevFragNum && siblingVertex->prevSibling != NULL && siblingVertex->prevSibling->childType == 0 && siblingVertex->prevSibling->prevSibling != NULL)
	{
		//there is an exon right before the siblingVertex; 6/24/2013, add the last condition to exclude the first or the last (see below) exon of the gene
		if (siblingVertex->prevSibling->childNum != 0 || siblingVertex->prevSibling->junctionInRange == NULL)// || rootVertex->prevSibling->junctionNum > 0)
		{
			cout << "Abnormal case in abundance estimation: Unrecognized prev exon." << endl;
			exit(1);
		} 
		else if (siblingVertex->prevSibling->junctionInRange->list->junc->type == frag_exon || siblingVertex->prevSibling->junctionInRange->list->junc->type == frag_retained_intron
			)//|| siblingVertex->prevSibling->junctionInRange->list->junc->type == frag_junction)
		{
			fragArray[prevFragCnt++] = siblingVertex->prevSibling->junctionInRange->list->junc;
			hasPrevExon = true;
		}

		siblingVertex = siblingVertex->prevSibling;
	}

	if (rootVertex->childType == 3)
	{
		if (hasPrevExon == true)
			iLoop = prevFragCnt;
		else
			iLoop = 0;
		//build fragArray from range junction list
		curJunc = rootVertex->junctionInRange->list;
		while (curJunc != NULL)
		{
			if (curJunc->junc->type == frag_exon || curJunc->junc->type == frag_retained_intron || curJunc->junc->type == frag_junction)
			{
				fragArray[iLoop++] = curJunc->junc;
			}
			curJunc = curJunc->next;
		}
	}
	else if (rootVertex->childType == 2)
	{
		if (hasPrevExon == true)
			iLoop = prevFragCnt;
		else
			iLoop = 0;
		//build fragArray from all children
		curEdge = rootVertex->child;
		while (curEdge != NULL)
		{
			curVertex = curEdge->linkedVertex;
			if (curVertex->estimated == false)
			{
				curJunc = curVertex->junctionInRange->list;
				while (curJunc != NULL)
				{
					if (curJunc->junc->type == frag_exon || curJunc->junc->type == frag_retained_intron || curJunc->junc->type == frag_junction)
					{
						fragArray[iLoop++] = curJunc->junc;
					}
					curJunc = curJunc->next;
				}
			} 
			else
			{
				fragArray[iLoop++] = curVertex->representative;
			}

			curEdge = curEdge->next;
		}
	}


	siblingVertex = rootVertex;
	while (nextFragCnt < MaxNextFragNum && siblingVertex->nextSibling != NULL && siblingVertex->nextSibling->childType == 0 && siblingVertex->nextSibling->nextSibling != NULL)
	{
		//there is an exon right after the siblingVertex
		if (siblingVertex->nextSibling->childNum != 0 || siblingVertex->nextSibling->junctionInRange == NULL)// || rootVertex->prevSibling->junctionNum > 0)
		{
			cout << "Abnormal case in abundance estimation: Unrecognized next exon." << endl;
			exit(1);
		} 
		else if (siblingVertex->nextSibling->junctionInRange->list->junc->type == frag_exon || siblingVertex->nextSibling->junctionInRange->list->junc->type == frag_retained_intron
			)//|| siblingVertex->nextSibling->junctionInRange->list->junc->type == frag_junction)
		{
			fragArray[prevFragCnt + fragmentCnt_est + nextFragCnt++] = siblingVertex->nextSibling->junctionInRange->list->junc;
			hasNextExon = true;
		}

		siblingVertex = siblingVertex->nextSibling;
	}


	//build transcript-exon matrix
	curEdge = rootVertex->child;
	iLoop = 0;
	while (curEdge != NULL)
	{
		curVertex = curEdge->linkedVertex;

		for (jLoop = 0; jLoop < prevFragCnt + fragmentCnt_est + nextFragCnt; jLoop++)
		{
			A[iLoop][jLoop] = false;
		}
		if (hasPrevExon == true)
		{
			if (abs(curVertex->rangeLow - rootVertex->rangeLow) <= 1)
			{
				for (jLoop = 0; jLoop < prevFragCnt; ++jLoop)
				{
					A[iLoop][jLoop] = true;
					transFragCnt[iLoop] += 1;
					if (fragArray[jLoop]->type == frag_exon || fragArray[jLoop]->type == frag_retained_intron)
						transLength[iLoop] += abs(fragArray[jLoop]->end - fragArray[jLoop]->start + 1);
					if (fragArray[jLoop]->type == frag_junction)
						transLength[iLoop] += 0;//fakeJunctionLength;
				}
			}
		} 
		if (hasNextExon == true)
		{
			if (abs(rootVertex->rangeHigh - curVertex->rangeHigh) <= 1)
			{
				for (jLoop = 0; jLoop < nextFragCnt; ++jLoop)
				{
					A[iLoop][prevFragCnt + fragmentCnt_est + jLoop] = true;
					transFragCnt[iLoop] += 1;
					if (fragArray[prevFragCnt + fragmentCnt_est + jLoop]->type == frag_exon || fragArray[prevFragCnt + fragmentCnt_est + jLoop]->type == frag_retained_intron)
						transLength[iLoop] += abs(fragArray[prevFragCnt + fragmentCnt_est + jLoop]->end - fragArray[prevFragCnt + fragmentCnt_est + jLoop]->start + 1);
					if (fragArray[prevFragCnt + fragmentCnt_est + jLoop]->type == frag_junction)
						transLength[iLoop] += 0;//fakeJunctionLength;
				}
			}
		} 

		if (curVertex->estimated == false)
		{
			curJunc = curVertex->junctionInRange->list;
			while (curJunc != NULL)
			{
				if (curJunc->junc->type == frag_exon || curJunc->junc->type == frag_retained_intron || curJunc->junc->type == frag_junction)
				{
					for (jLoop = prevFragCnt; jLoop < prevFragCnt + fragmentCnt_est; jLoop++)
					{
						if (fragArray[jLoop]->frag_id == curJunc->junc->frag_id)
						{
							A[iLoop][jLoop] = true;
							transFragCnt[iLoop] += 1;
							if (fragArray[jLoop]->type == frag_exon || fragArray[jLoop]->type == frag_retained_intron)
								transLength[iLoop] += abs(fragArray[jLoop]->end - fragArray[jLoop]->start + 1);
							if (fragArray[jLoop]->type == frag_junction)
								transLength[iLoop] += 0;//fakeJunctionLength;

							break;
						} 
					}

					errorRatio_transExpr[iLoop] += curJunc->junc->support[index_tissue];
					errorRatio_fragCnt[iLoop] += 1;
				}
				curJunc = curJunc->next;
			}
		} 
		else
		{
			for (jLoop = prevFragCnt; jLoop < prevFragCnt + fragmentCnt_est; jLoop++)
			{
				if (fragArray[jLoop]->frag_id == curVertex->representative->frag_id)
				{
					A[iLoop][jLoop] = true;
					transFragCnt[iLoop] += 1;
					transLength[iLoop] += abs(fragArray[jLoop]->end - fragArray[jLoop]->start + 1);

					break;
				} 
			}

			errorRatio_transExpr[iLoop] += curVertex->representative->support[index_tissue];
			errorRatio_fragCnt[iLoop] += 1;
		}

		iLoop++;
		curEdge = curEdge->next;
	}



	//estimation
	double *P = new double [fragmentCnt_est + MaxPrevFragNum + MaxNextFragNum]; // probability that reads fall into each exon
	double *Y = new double [fragmentCnt_est + MaxPrevFragNum + MaxNextFragNum]; // observed exon expression
	double *Q = new double [transcriptCnt]; //proportions of alternative transcripts
	double *Qold = new double [transcriptCnt]; //new Q vector
	double *Qtmp = new double [transcriptCnt];
	double *Qdelta = new double [transcriptCnt]; //change on Q
	double *Z = new double [transcriptCnt]; //expectation of read count, will be derived from C
	double *C = new double [transcriptCnt]; //expectation of read coverage
	double lambda; 


	//fill in Y
	for (iLoop = 0; iLoop < prevFragCnt + fragmentCnt_est + nextFragCnt; iLoop++)
	{
		if ((fragArray[iLoop]->support)[index_tissue] < 0.0)
			(fragArray[iLoop]->support)[index_tissue] = 0.0;
		Y[iLoop] = (fragArray[iLoop]->support)[index_tissue];
	}

	//initiate Q, Z, transExonCnt
	for (iLoop = 0; iLoop < transcriptCnt; iLoop++)
	{
		Q[iLoop] = 1.0 / transcriptCnt;
		Qtmp[iLoop] = 0.0;
		Qold[iLoop] = 0.0;
		Qdelta[iLoop] = 0.0;

		Z[iLoop] = 0.0;
		C[iLoop] = 0.0;
	}


	//calculate P
	long lengthSum = 0;
	for (iLoop = 0; iLoop < prevFragCnt + fragmentCnt_est + nextFragCnt; iLoop++)
	{
		if (fragArray[iLoop]->type == frag_exon || fragArray[iLoop]->type == frag_retained_intron)
		{
			lengthSum += abs(fragArray[iLoop]->end - fragArray[iLoop]->start + 1);
		} 
	}
	if (hasPrevExon == false) {
		lengthSum += fakeJunctionLength;
	}
	for (iLoop = 0; iLoop < prevFragCnt + fragmentCnt_est + nextFragCnt; iLoop++)
	{
		if (fragArray[iLoop]->type == frag_exon || fragArray[iLoop]->type == frag_retained_intron)
		{
			P[iLoop] = double(abs(fragArray[iLoop]->end - fragArray[iLoop]->start + 1)) / lengthSum;
		}
		else if (fragArray[iLoop]->type == frag_junction)
		{
			P[iLoop] = double(fakeJunctionLength) / lengthSum;
		}
	}



	//EM
	double sum_count, sum_denominator, sum_numerator, tmpvalue;
	bool stoppable = false, tmpFlag;
	long MAX_ITER = 100, iterCnt = 0;

	//12/18/2012 limit the max # of transcripts in an estimation
	if (transcriptCnt > 30)
		iterCnt = MAX_ITER;

	while (iterCnt < MAX_ITER && stoppable == false)
	{
		for (iLoop = 0; iLoop < transcriptCnt; iLoop++)
		{
			Qold[iLoop] = Q[iLoop];
		}

		//E-step
		for (iLoop = 0; iLoop < transcriptCnt; iLoop++)
		{
			sum_count = 0.0;

			for (jLoop = 0; jLoop < prevFragCnt + fragmentCnt_est + nextFragCnt; jLoop++)
			{
				if (A[iLoop][jLoop] == true)
				{
					sum_denominator = 0.0;
					for (kLoop = 0; kLoop < transcriptCnt; kLoop++)
					{
						if (A[kLoop][jLoop] == true)
						{
							sum_denominator += Q[kLoop]; // P[jLoop] * Q[kLoop];
						}
					}

					if (sum_denominator > 0.0)
					{
						tmpvalue = Q[iLoop] * Y[jLoop] / sum_denominator; // P[jLoop] * Q[iLoop] * Y[jLoop] / sum_denominator;
						sum_count += tmpvalue;
					}
				}
			}

			C[iLoop] = sum_count / transFragCnt[iLoop];
			Z[iLoop] = C[iLoop] * transLength[iLoop] / readlength;
		}


		//M-step
		//calculate lambda
		lambda = 0.0;
		for (iLoop = 0; iLoop < transcriptCnt; iLoop++)
		{
			lambda += Z[iLoop];
		}

		//calculate Q
		//get Q iteratively until stable
		if (transcriptCnt > 1)
		{
			tmpFlag = false;
			while (tmpFlag == false)
			{
				for (iLoop = 0; iLoop < transcriptCnt; iLoop++)
				{
					Qtmp[iLoop] = Q[iLoop];
				}
				for (iLoop = 0; iLoop < transcriptCnt; iLoop++)
				{
					sum_denominator = 0.0;
					for (jLoop = 0; jLoop < prevFragCnt + fragmentCnt_est + nextFragCnt; jLoop++)
					{
						if (A[iLoop][jLoop] == true)
						{
							sum_denominator += P[jLoop];
						}
					}
					sum_denominator *= (lambda - Z[iLoop]);

					sum_numerator = 0.0;
					for (jLoop = 0; jLoop < transcriptCnt; jLoop++)
					{
						if (jLoop == iLoop)
						{
							continue;
						}

						sum_count = 0.0;
						for (kLoop = 0; kLoop < prevFragCnt + fragmentCnt_est + nextFragCnt; kLoop++)
						{
							if (A[jLoop][kLoop] == true)
							{
								sum_count += P[kLoop];
							}
						}
						sum_count *= Qtmp[jLoop];

						sum_numerator += sum_count;
					}
					sum_numerator *= Z[iLoop];

					if (sum_denominator > 0.0)
					{
						Qtmp[iLoop] = sum_numerator / sum_denominator;
					} 
					else
					{
						if (lambda == Z[iLoop])
						{
							Qtmp[iLoop] = 1.0;
						}
						else
						{
							Qtmp[iLoop] = 0.0;
						}
					}
				}

				sum_count = 0.0;
				for (iLoop = 0; iLoop < transcriptCnt; iLoop++)
				{
					sum_count += Qtmp[iLoop];
				}
				for (iLoop = 0; iLoop < transcriptCnt; iLoop++)
				{
					//calculate change on Q
					if (sum_count > 0.0)
					{
						Qtmp[iLoop] = Qtmp[iLoop] / sum_count;
					} 
					else
					{
						Qtmp[iLoop] = 1.0 / transcriptCnt;
					}
					Qdelta[iLoop] = fabs(Qtmp[iLoop] - Q[iLoop]);
					Q[iLoop] = Qtmp[iLoop];
				}

				tmpFlag = true;
				for (iLoop = 0; iLoop < transcriptCnt; iLoop++)
				{
					if (Q[iLoop] > 0.0 && fabs(Qdelta[iLoop] / Q[iLoop]) > 0.001)
					{
						tmpFlag = false;
						break;
					}
				}
			}

			for (iLoop = 0; iLoop < transcriptCnt; iLoop++)
			{
				//calculate change on Q
				Qdelta[iLoop] = fabs(Qold[iLoop] - Q[iLoop]);
			}

			iterCnt++;
			tmpFlag = true;
			for (iLoop = 0; iLoop < transcriptCnt; iLoop++)
			{
				if (Q[iLoop] > 0.0 && fabs(Qdelta[iLoop] / Q[iLoop]) > 0.001)
				{
					tmpFlag = false;
					break;
				}
			}
			if (iterCnt >= 3 && tmpFlag == true)
			{
				stoppable = true;
			}	
		} 
		else
		{
			stoppable = true;
		}

	}

	//one more E-step
	for (iLoop = 0; iLoop < transcriptCnt; iLoop++)
	{
		sum_count = 0.0;

		for (jLoop = 0; jLoop < prevFragCnt + fragmentCnt_est + nextFragCnt; jLoop++)
		{
			if (A[iLoop][jLoop] == true)
			{
				sum_denominator = 0.0;
				for (kLoop = 0; kLoop < transcriptCnt; kLoop++)
				{
					if (A[kLoop][jLoop] == true)
					{
						sum_denominator += Q[kLoop]; // P[jLoop] * Q[kLoop];
					}
				}

				if (sum_denominator > 0.0)
				{
					sum_count += Q[iLoop] * Y[jLoop] / sum_denominator; // P[jLoop] * Q[iLoop] * Y[jLoop] / sum_denominator;
				}
			}
		}

		C[iLoop] = sum_count / transFragCnt[iLoop];
		Z[iLoop] = C[iLoop] * transLength[iLoop] / readlength;
	}


	for (iLoop = 0; iLoop < transcriptCnt; iLoop++)
	{


		if (Z[iLoop] < 0)
		{
			cout << "neg_Z\t";
		}
	}

#ifdef DEBUG_OUTPUT
	ofstream debug_out_estimation("debug_output/debug_estimation.txt", fstream::app);
#endif

	//copy results & calculate MSE
	MSE = 0.0;
	double min_path_support = MAX_CHR_LENGTH;
	curEdge = rootVertex->child;
	iLoop = 0;
	while (curEdge != NULL)
	{
		curVertex = curEdge->linkedVertex;

		//!!!Do not use support for the whole transcript fragment here! Use only the portion belonging to the vertex!
		//		curVertex->support[index_tissue] = Z[iLoop];
		if (transFragCnt[iLoop] > 0)
		{
			curVertex->support[index_tissue] = C[iLoop]; //Z[iLoop] / transFragCnt[iLoop]; //evenly divided, for coverage only

			if (C[iLoop] < min_path_support)
				min_path_support = C[iLoop];
		} 
		else
		{
			curVertex->support[index_tissue] = 0.0;
		}

		if (errorRatio_fragCnt[iLoop] > 0 && curVertex->support[index_tissue] > 10 && abs(errorRatio_transExpr[iLoop]/errorRatio_fragCnt[iLoop] - curVertex->support[index_tissue]) / curVertex->support[index_tissue] > MSE)
		{
			curVertex->obs_support[index_tissue] = errorRatio_transExpr[iLoop]/errorRatio_fragCnt[iLoop];
			MSE = abs(errorRatio_transExpr[iLoop]/errorRatio_fragCnt[iLoop] - curVertex->support[index_tissue]) / curVertex->support[index_tissue];
		}

		curVertex->proportion[index_tissue] = Q[iLoop];
#ifdef DEBUG_OUTPUT
		if (rootVertex->rangeLow == 57227143 && rootVertex->rangeHigh == 57242546)
			debug_out_estimation << curVertex->proportion[index_tissue] << "\t";
#endif

		if (Q[iLoop] < 0)
		{
			cout << "neg_Q\t";
		}

		iLoop++;
		curEdge = curEdge->next;
	}
#ifdef DEBUG_OUTPUT
	if (rootVertex->rangeLow == 57227143 && rootVertex->rangeHigh == 57242546)
		debug_out_estimation << endl;
#endif

	rootVertex->min_path_support[index_tissue] = min_path_support;

	//calculate the log likelihood
	double loglikelihood = 0.0, lambda_i;

	// 	lambda = 0.0;
	// 	for (iLoop = 0; iLoop < transcriptCnt; iLoop++)
	// 	{
	// 		lambda += Z[iLoop];
	// 	}
	// 
	// 	sum_numerator = 0.0;
	// 	for (jLoop = 0; jLoop < transcriptCnt; jLoop++)
	// 	{
	// 		sum_count = 0.0;
	// 		for (kLoop = 0; kLoop < prevFragCnt + fragmentCnt_est + nextFragCnt; kLoop++)
	// 		{
	// 			if (A[jLoop][kLoop] == true)
	// 			{
	// 				sum_count += P[kLoop];
	// 			}
	// 		}
	// 		sum_count *= Q[jLoop];
	// 
	// 		sum_numerator += sum_count;
	// 	}
	// 
	// 	for (iLoop = 0; iLoop < transcriptCnt; iLoop++)
	// 	{
	// 		sum_count = 0.0;
	// 		for (kLoop = 0; kLoop < prevFragCnt + fragmentCnt_est + nextFragCnt; kLoop++)
	// 		{
	// 			if (A[iLoop][kLoop] == true)
	// 			{
	// 				sum_count += P[kLoop];
	// 			}
	// 		}
	// 		sum_count *= Q[iLoop];
	// 		if (sum_numerator > 0.0)
	// 		{
	// 			lambda_i = sum_count / sum_numerator;
	// 		} 
	// 		else
	// 		{
	// 			lambda_i = 0.0;
	// 		}
	// 
	// 		sum_count = 0.0;
	// 		for (kLoop = 0; kLoop < prevFragCnt + fragmentCnt_est + nextFragCnt; kLoop++)
	// 		{
	// 			if (A[iLoop][kLoop] == true)
	// 			{
	// 				sum_count += P[kLoop];
	// 			}
	// 		}
	// 
	// 		loglikelihood += -lambda_i + Z[iLoop] * log(lambda_i);
	// 
	// 		for (kLoop = 0; kLoop < prevFragCnt + fragmentCnt_est + nextFragCnt; kLoop++)
	// 		{
	// 			if (A[iLoop][kLoop] == true)
	// 			{
	// 				if (long(Z[iLoop]/transFragCnt[iLoop] >= 0))
	// 				{
	// 					// 					if (ceil(Z[iLoop]/transExonCnt[iLoop]) > 1000)
	// 					// 					{
	// 					// 						cout << "Warning: calculating factorial of a number larger than 10" << endl;
	// 					// 					}
	// 					//loglikelihood += - log(double(factorial(long(ceil(Z[iLoop]/transExonCnt[iLoop])))));
	// 					loglikelihood += - logFactorial(long(ceil(Z[iLoop]/transFragCnt[iLoop])));
	// 				}
	// 				else
	// 				{
	// 					cout << "ERROR: estimation, factorial" << endl;
	// 				}
	// 				loglikelihood += Z[iLoop]/transFragCnt[iLoop] * log(P[kLoop] / sum_count);
	// 			}
	// 		}
	// 	}

#ifdef DEBUG_OUTPUT
	debug_out_estimation.close();
#endif

	//clean up
	delete [] fragArray;
	delete [] transFragCnt;
	delete [] errorRatio_fragCnt;
	delete [] errorRatio_transExpr;
	for (iLoop = 0; iLoop < transcriptCnt; iLoop++)
	{
		delete [] A[iLoop];
	}
	delete [] A;
	delete [] P;
	delete [] Y;
	delete [] Z;
	delete [] Q;
	delete [] Qold;
	delete [] Qtmp;
	delete [] Qdelta;
	delete [] C;
	delete [] transLength;

	return loglikelihood;
}


void vertexForAlterSpliceSite(fragment *mergedJunction, GTvertex *fatherVertex)
{
	//build a representative GTvertex for the merged alternative splice site junctions
	//add the vertex to the fatherVertex
	if (mergedJunction->alter == NULL)
	{
		return;
	}

	GTvertex *newVertex, *resultVertex;
	newVertex = new GTvertex;

	alter_junction *cur_alter_junction;
	GTedge *newEdge;
	RangeJunctionList *newRangeJunctionList;
	rangeJunction *newRangeJunction;
	bool separable = false;
	long iLoop, jLoop;
	fragment **alterFragList;

	alterFragList = new fragment * [mergedJunction->alterFragCnt + 1];

	//build the junction list for the vertex
	newRangeJunctionList = new RangeJunctionList;
	newRangeJunctionList->rangeLow = mergedJunction->start;
	newRangeJunctionList->rangeHigh = mergedJunction->end;
	iLoop = 1;
	cur_alter_junction = mergedJunction->alter;
	while (cur_alter_junction != NULL)
	{
		if (cur_alter_junction->juncInfo->alter != NULL)
		{
			cout << "Warning: nested alternative splice site." << endl;
		}
		alterFragList[iLoop] = cur_alter_junction->juncInfo->clone(SUPPORT_VECTOR_SIZE);

		//reload the real start and end!
		if (alterFragList[iLoop]->start_real != 0)
		{
			alterFragList[iLoop]->start = alterFragList[iLoop]->start_real;
		}
		if (alterFragList[iLoop]->end_real != 0)
		{
			alterFragList[iLoop]->end = alterFragList[iLoop]->end_real;
		}

		iLoop++;
		cur_alter_junction = cur_alter_junction->next;
	}

	//selection sort
	fragment *curFragment, *tmpFragment;
	long minIndex;
	for (iLoop = 1; iLoop < mergedJunction->alterFragCnt; iLoop++)
	{
		minIndex = iLoop;
		for (jLoop = iLoop + 1; jLoop <= mergedJunction->alterFragCnt; jLoop++)
		{
			curFragment = alterFragList[jLoop];
			if (curFragment->start < alterFragList[minIndex]->start || curFragment->start == alterFragList[minIndex]->start && curFragment->end < alterFragList[minIndex]->end)
			{
				//curFragment is less than alterFragList[minIndex], replace
				minIndex = jLoop;
			}
		}
		//exchange selection
		if (minIndex != iLoop)
		{
			tmpFragment = alterFragList[iLoop];
			alterFragList[iLoop] = alterFragList[minIndex];
			alterFragList[minIndex] = tmpFragment;
		}
	}

	for (iLoop = mergedJunction->alterFragCnt; iLoop >= 1 ; iLoop--)
	{
		newRangeJunction = new rangeJunction;
		newRangeJunction->junc = alterFragList[iLoop];
		newRangeJunction->next = newRangeJunctionList->list;
		newRangeJunctionList->list = newRangeJunction;
	}

	delete [] alterFragList;


	//fill the vertex
	newVertex->level = fatherVertex->level; //same level as fatherVertex
	newVertex->rangeLow = mergedJunction->start;
	newVertex->rangeHigh = mergedJunction->end;
	if (fatherVertex->childType == 0)
	{
		//the merged junction is isolated
		newVertex->level = fatherVertex->level; //same level as fatherVertex
		newVertex->prevSibling = fatherVertex->prevSibling; //same siblings as the father
		newVertex->nextSibling = fatherVertex->nextSibling; //same siblings as the father
	} 
	else
	{
		//the merged junction is in a path (dependent path)
		newVertex->level = fatherVertex->level + 1; //same level as fatherVertex
		newVertex->prevSibling = NULL;
		newVertex->nextSibling = NULL;
	}
	newVertex->junctionInRange = newRangeJunctionList;
	resultVertex = newVertex;

	//construct the paths
	int tmp_numPath;
	newRangeJunctionList = separateDepPath(newRangeJunctionList, separable, tmp_numPath, 0);

	if (separable == true)
	{
		resultVertex->childType = 3;

		while (newRangeJunctionList != NULL)
		{
			newVertex = new GTvertex;
			newVertex->level = resultVertex->level + 1;
			newVertex->junctionInRange = newRangeJunctionList;
			newVertex->rangeLow = newRangeJunctionList->rangeLow;
			newVertex->rangeHigh = newRangeJunctionList->rangeHigh;
			newVertex->childType = 0;

			newRangeJunctionList = newRangeJunctionList->nextList;
			newVertex->junctionInRange->nextList = NULL;

			newEdge = new GTedge;
			newEdge->linkedVertex = newVertex;
			newEdge->next = resultVertex->child;
			resultVertex->child = newEdge;
		}
	} 
	else
	{
		resultVertex->junctionInRange = newRangeJunctionList;
		resultVertex->childType = 0;
	}

	//count the vertex
	countGTree(resultVertex);

	resultVertex->ASMcategory = alter_splice_site;

	newEdge = new GTedge;
	newEdge->linkedVertex = resultVertex;
	newEdge->next = fatherVertex->alterSpliceSite;
	fatherVertex->alterSpliceSite = newEdge;

	return;
}




fragment* makeRepresentative(GTvertex *rootVertex)
{
	//make representative for rootVertex
	//the representative will consist of a junction (if needed), an exon (length = all exonic region, support = all sum), and a junction (if needed)
	//currently only use an exon to represent
	GTvertex *curVertex;
	GTedge *curEdge;
	rangeJunction *curJunc;
	long exonicLength = 0;
	double *exonicSupportSum = new double [SUPPORT_VECTOR_SIZE];
	int exonicChildNum = 0;
	int iLoop;

	for (iLoop = 0; iLoop < SUPPORT_VECTOR_SIZE; iLoop++)
	{
		exonicSupportSum[iLoop] = 0.0;
	}

	if (rootVertex->childType == 3 || rootVertex->childType == 2)
	{
		curEdge = rootVertex->child;
		while (curEdge != NULL)
		{
			curVertex = curEdge->linkedVertex;

			for (iLoop = 0; iLoop < SUPPORT_VECTOR_SIZE; iLoop++)
			{
				exonicSupportSum[iLoop] += curVertex->support[iLoop];
			}

			curEdge = curEdge->next;
		}
	}
	else if (rootVertex->childType == 1)
	{
		curEdge = rootVertex->child;
		while (curEdge != NULL)
		{
			curVertex = curEdge->linkedVertex;

			if (curVertex->estimated == false)
			{
				curJunc = curVertex->junctionInRange->list;
				while (curJunc != NULL)
				{
					if (curJunc->junc->type == frag_exon || curJunc->junc->type == frag_retained_intron || curJunc->junc->type == frag_junction)
					{
						for (iLoop = 0; iLoop < SUPPORT_VECTOR_SIZE; iLoop++)
						{
							exonicSupportSum[iLoop] += curJunc->junc->support[iLoop];
						}
						exonicChildNum++;
					}
					curJunc = curJunc->next;
				}
			} 
			else
			{
				for (iLoop = 0; iLoop < SUPPORT_VECTOR_SIZE; iLoop++)
				{
					exonicSupportSum[iLoop] += curVertex->representative->support[iLoop];
				}
				exonicChildNum++;
			}

			if (curVertex->childNum != 0 || curVertex->exonNum != 0)
				exonicLength += curVertex->rangeHigh - curVertex->rangeLow + 1;

			curEdge = curEdge->next;
		}

		for (iLoop = 0; iLoop < SUPPORT_VECTOR_SIZE; iLoop++)
		{
			if (exonicChildNum > 0)
			{
				exonicSupportSum[iLoop] = exonicSupportSum[iLoop] / exonicChildNum;
			} 
			else
			{
				exonicSupportSum[iLoop] = 0.0;
			}
		}
	}
	else if (rootVertex->childType == 0)
	{
		curJunc = rootVertex->junctionInRange->list;
		while (curJunc != NULL)
		{
			if (curJunc->junc->type == frag_exon || curJunc->junc->type == frag_retained_intron || curJunc->junc->type == frag_junction)
			{
				for (iLoop = 0; iLoop < SUPPORT_VECTOR_SIZE; ++iLoop)
				{
					if (curJunc->junc->support[iLoop] > exonicSupportSum[iLoop])
					{
						exonicSupportSum[iLoop] = curJunc->junc->support[iLoop];
					}
				}
			}
			curJunc = curJunc->next;
		}
	}

	fragment *newFragment;
	newFragment = new fragment(SUPPORT_VECTOR_SIZE);
	newFragment->frag_id = ++fragmentID_Cnt;
	newFragment->start = 0;
	newFragment->end = 200;
	newFragment->type = frag_exon;
	for (iLoop = 0; iLoop < SUPPORT_VECTOR_SIZE; iLoop++)
	{
		newFragment->support[iLoop] = exonicSupportSum[iLoop];
		rootVertex->support[iLoop] = exonicSupportSum[iLoop];
	}

	delete [] exonicSupportSum;

	return newFragment;
}


bool preCountGTree(GTvertex *rootVertex)
{
	//collect basic counts on GTree for further statistics
	//collect child count and fragment count for every vertex
	//return true for leaves
	GTvertex *curVertex;
	GTedge *curEdge;
	rangeJunction *curJunc;
	double *fragSupportSum = new double [SUPPORT_VECTOR_SIZE], *exonSupportSum = new double [SUPPORT_VECTOR_SIZE]; 
	int juncCnt, exonCnt, fragCnt, tmp, iLoop;

	if (rootVertex->child == NULL)
	{
		//leaf
		juncCnt = 0;
		exonCnt = 0;
		fragCnt = 0;
		for (iLoop = 0; iLoop < SUPPORT_VECTOR_SIZE; iLoop++)
		{
			fragSupportSum[iLoop] = 0;
			exonSupportSum[iLoop] = 0;
		}

		curJunc = rootVertex->junctionInRange->list;
		while (curJunc != NULL)
		{
			fragCnt++;
			if (curJunc->junc->type == frag_junction)
			{
				juncCnt++;
			}
			else if (curJunc->junc->type == frag_exon || curJunc->junc->type == frag_retained_intron)
			{
				exonCnt++;
				for (iLoop = 0; iLoop < SUPPORT_VECTOR_SIZE; iLoop++)
				{
					exonSupportSum[iLoop] += curJunc->junc->support[iLoop];
				}
			}

			for (iLoop = 0; iLoop < SUPPORT_VECTOR_SIZE; iLoop++)
			{
				fragSupportSum[iLoop] += curJunc->junc->support[iLoop];
			}

			curJunc = curJunc->next;
		}
		rootVertex->junctionNum = juncCnt;
		rootVertex->exonNum = exonCnt;

		rootVertex->estimated = false;
		rootVertex->estimate_exonNum = 0; //leaf has nothing to estimate

		// 		if (juncCnt <= 1)
		// 		{
		// 			//one junction or one exon
		// 			for (iLoop = 0; iLoop < SUPPORT_VECTOR_SIZE; iLoop++)
		// 			{
		// 				rootVertex->support[iLoop] = rootVertex->junctionInRange->list->junc->support[iLoop];
		// 			}
		// 		}
		// 		else
		// 		{
		// 			//more junctions
		// 			//just calculate average for now, will change to transcript proportion later
		// 			// 			for (tmp = 0; tmp < SUPPORT_VECTOR_SIZE; tmp++)
		// 			// 			{
		// 			// 				rootVertex->support[tmp] = fragSupportSum[tmp] / fragCnt;
		// 			// 			}
		// 		}
		for (iLoop = 0; iLoop < SUPPORT_VECTOR_SIZE; iLoop++)
		{
			//count the exonic expression
			if (exonCnt > 0)
			{
				rootVertex->support[iLoop] = double(exonSupportSum[iLoop]) / exonCnt;
			} 
			else
			{
				rootVertex->support[iLoop] = 0.0;
			}
		}

		//rootVertex->anovaScore_support = Anova_test(rootVertex->support);

		delete [] fragSupportSum; delete [] exonSupportSum;
		return true;
	} 
	else
	{
		//extend children
		rootVertex->estimate_exonNum = 0;

		curEdge = rootVertex->child;
		while (curEdge != NULL)
		{
			(rootVertex->childNum)++;
			curVertex = curEdge->linkedVertex;
			if (preCountGTree(curVertex) == false)
			{
			}
			if (rootVertex->childType != 3)
			{
				rootVertex->junctionNum += curVertex->junctionNum;
				rootVertex->exonNum += curVertex->exonNum;

				if (curVertex->estimated == false)
				{
					rootVertex->estimate_exonNum += curVertex->exonNum;
				} 
				else
				{
					rootVertex->estimate_exonNum += 1; //representative exon
				}
			} 

			curEdge = curEdge->next;
		}

		if (rootVertex->childType == 3)
		{
			juncCnt = 0;
			exonCnt = 0;
			fragCnt = 0;

			curJunc = rootVertex->junctionInRange->list;
			while (curJunc != NULL)
			{
				fragCnt++;
				if (curJunc->junc->type == frag_junction)
				{
					juncCnt++;
				}
				else if (curJunc->junc->type == frag_exon || curJunc->junc->type == frag_retained_intron)
				{
					exonCnt++;
				}
				curJunc = curJunc->next;
			}
			rootVertex->junctionNum = juncCnt;
			rootVertex->exonNum = exonCnt;

			rootVertex->estimate_exonNum = exonCnt;
		}

		delete [] fragSupportSum; delete [] exonSupportSum;
		return false;
	}
}

long select_major_paths(GTvertex *rootVertex)
{
	long major_path_num = 0;
	double path_proportion_sum = 0.0, path_expression_sum = 0.0, path_proportion_max = 0, path_expression_max;
	GTedge *curEdge;
	GTvertex *curVertex;
	alternative_path *newAlterPath;
	int iLoop, expressedCnt;

	curEdge = rootVertex->child;
	while (curEdge != NULL)
	{
		curVertex = curEdge->linkedVertex;
		curVertex->has_novel_junction_major = curVertex->has_novel_junction;

		if (COUNT_MAJOR_PATHS_ONLY == true)
		{
			path_proportion_sum = 0.0; path_expression_sum = 0.0; path_proportion_max = 0.0; path_expression_max = 0.0; expressedCnt = 0;
			for (iLoop = 0; iLoop < SUPPORT_VECTOR_SIZE; iLoop++)
			{
				path_proportion_sum += curVertex->proportion[iLoop];
				path_expression_sum += curVertex->support[iLoop];
				if (curVertex->proportion[iLoop] >= path_proportion_max)
					path_proportion_max = curVertex->proportion[iLoop];
				if (curVertex->support[iLoop] >= path_expression_max)
					path_expression_max = curVertex->support[iLoop];
				if (curVertex->support[iLoop] >= thresh_path_expressed)
					++expressedCnt;
			}

			if (path_proportion_max >= 0.01 || (path_proportion_sum / SUPPORT_VECTOR_SIZE >= 0.005 && path_expression_sum / SUPPORT_VECTOR_SIZE >= thresh_path_expressed) || expressedCnt >= SUPPORT_VECTOR_SIZE*0.01)
			{
				//major path
				newAlterPath = new alternative_path;
				newAlterPath->path_start = curVertex->rangeLow;
				newAlterPath->path_end   = curVertex->rangeHigh;
				if (curVertex->prevSibling != NULL)
					newAlterPath->whole_path_start = curVertex->prevSibling->rangeLow;
				else
					newAlterPath->whole_path_start = curVertex->rangeLow;
				if (curVertex->nextSibling != NULL)
					newAlterPath->whole_path_end = curVertex->nextSibling->rangeHigh;
				else
					newAlterPath->whole_path_end = curVertex->rangeHigh;
				double sum = 0.0;
				for (iLoop = 0; iLoop < SUPPORT_VECTOR_SIZE; iLoop++)
				{
					newAlterPath->support[iLoop] = curVertex->support[iLoop];
					newAlterPath->proportion[iLoop] = curVertex->proportion[iLoop];
					sum += newAlterPath->support[iLoop];
				}
				newAlterPath->junctionNum = curVertex->junctionNum;
				newAlterPath->exonNum = curVertex->exonNum;
				newAlterPath->pathVertex = curVertex;
				newAlterPath->avg_support = sum / SUPPORT_VECTOR_SIZE;
				newAlterPath->transDirection = curVertex->transDirection;

				++major_path_num;
				newAlterPath->next = rootVertex->major_alter_paths;
				rootVertex->major_alter_paths = newAlterPath;

				if (curVertex->has_novel_junction)
				{
					rootVertex->has_novel_junction_major = true;
				}
				
			}
		} 
		else
		{
			//all paths
			newAlterPath = new alternative_path;
			newAlterPath->path_start = curVertex->rangeLow;
			newAlterPath->path_end   = curVertex->rangeHigh;
			if (curVertex->prevSibling != NULL)
				newAlterPath->whole_path_start = curVertex->prevSibling->rangeLow;
			else
				newAlterPath->whole_path_start = curVertex->rangeLow;
			if (curVertex->nextSibling != NULL)
				newAlterPath->whole_path_end = curVertex->nextSibling->rangeHigh;
			else
				newAlterPath->whole_path_end = curVertex->rangeHigh;
			double sum = 0.0;
			for (iLoop = 0; iLoop < SUPPORT_VECTOR_SIZE; iLoop++)
			{
				newAlterPath->support[iLoop] = curVertex->support[iLoop];
				newAlterPath->proportion[iLoop] = curVertex->proportion[iLoop];
				sum += newAlterPath->support[iLoop];
			}
			newAlterPath->junctionNum = curVertex->junctionNum;
			newAlterPath->exonNum = curVertex->exonNum;
			newAlterPath->pathVertex = curVertex;
			newAlterPath->avg_support = sum / SUPPORT_VECTOR_SIZE;
			newAlterPath->transDirection = curVertex->transDirection;

			++major_path_num;
			newAlterPath->next = rootVertex->major_alter_paths;
			rootVertex->major_alter_paths = newAlterPath;
		}

		curEdge = curEdge->next;
	}

	if (!COUNT_MAJOR_PATHS_ONLY)
		rootVertex->has_novel_junction_major = rootVertex->has_novel_junction;

	rootVertex->major_alter_paths_num = major_path_num;
	return major_path_num;
}


bool countGTree(GTvertex *rootVertex)
{
	//collect basic counts on GTree for further statistics
	//collect child count and fragment count for every vertex
	//return true for leaves
	GTvertex *curVertex;
	GTedge *curEdge;
	alternative_path *curAlterPath;
	rangeJunction *curJunc;
	double *fragSupportSum = new double [SUPPORT_VECTOR_SIZE], *exonSupportSum = new double [SUPPORT_VECTOR_SIZE], totalSupport;
	int juncCnt, exonCnt, fragCnt, tmp, iLoop, jLoop;
	bool allChildrenAreLeaves;
	double *Parray_support, *Parray_proportion, *Qarray_support, *Qarray_proportion, **proportion_matrix, MSE_estimation;

	//for debug
//	time_t t; string cur_time;
//	time(&t); cur_time = ctime(&t); cur_time.erase(cur_time.size() - 1);
//	cout << "<" << cur_time << "> " << rootVertex->vertex_id << endl;
//	cout << rootVertex->vertex_id << endl;


	if (rootVertex->child == NULL)
	{
		//leaf
		juncCnt = 0;
		exonCnt = 0;
		fragCnt = 0;
		for (iLoop = 0; iLoop < SUPPORT_VECTOR_SIZE; iLoop++)
		{
			fragSupportSum[iLoop] = 0;
			exonSupportSum[iLoop] = 0;
		}

		curJunc = rootVertex->junctionInRange->list;
		while (curJunc != NULL)
		{
			fragCnt++;
			if (curJunc->junc->type == frag_junction)
			{
				juncCnt++;
			}
			else if (curJunc->junc->type == frag_exon || curJunc->junc->type == frag_retained_intron)
			{
				exonCnt++;
				for (iLoop = 0; iLoop < SUPPORT_VECTOR_SIZE; iLoop++)
				{
					exonSupportSum[iLoop] += curJunc->junc->support[iLoop];
				}
			}

			for (iLoop = 0; iLoop < SUPPORT_VECTOR_SIZE; iLoop++)
			{
				fragSupportSum[iLoop] += curJunc->junc->support[iLoop];
			}

			curJunc = curJunc->next;
		}
		rootVertex->junctionNum = juncCnt;
		rootVertex->exonNum = exonCnt;

		//4/11/2013 //6/8/2013 add path_extended to prevent the case that the alternative paths of a vertex are falsely labeled as estimated, which results from trying to separate dependent paths before trying to separate independent paths
		if (fragCnt == 1 || rootVertex->path_extended == false)
		{
			//single exon or single junction, label as not estimated, 
			rootVertex->estimated = false;
			rootVertex->estimate_exonNum = 0; //leaf has nothing to estimate

			for (iLoop = 0; iLoop < SUPPORT_VECTOR_SIZE; iLoop++)
			{
				//count the exonic expression
				if (exonCnt > 0)
				{
					rootVertex->support[iLoop] = double(exonSupportSum[iLoop]) / exonCnt;
				} 
				else
				{
					rootVertex->support[iLoop] = 0.0;
				}
			}
		}
		else if (fragCnt > 1)
		{
			//no child but multiple fragments => region whose path not enumerated, label as estimated
			rootVertex->estimated = true;
			rootVertex->representative = makeRepresentative(rootVertex);
		}

		// 		if (juncCnt <= 1)
		// 		{
		// 			//one junction or one exon
		// 			for (iLoop = 0; iLoop < SUPPORT_VECTOR_SIZE; iLoop++)
		// 			{
		// 				rootVertex->support[iLoop] = rootVertex->junctionInRange->list->junc->support[iLoop];
		// 			}
		// 		}
		// 		else
		// 		{
		// 			//more junctions
		// 			//just calculate average for now, will change to transcript proportion later
		// 			// 			for (tmp = 0; tmp < SUPPORT_VECTOR_SIZE; tmp++)
		// 			// 			{
		// 			// 				rootVertex->support[tmp] = fragSupportSum[tmp] / fragCnt;
		// 			// 			}
		// 		}

		//rootVertex->anovaScore_support = Anova_test(rootVertex->support);


		//alternative splice sites
		if (fragCnt == 1 && juncCnt == 1 && rootVertex->junctionInRange->list->junc->alter != NULL)
		{
			vertexForAlterSpliceSite(rootVertex->junctionInRange->list->junc, rootVertex);
		}

// 		Parray_support = new double [SUPPORT_VECTOR_SIZE / 2 + 1];
// 		Qarray_support = new double [SUPPORT_VECTOR_SIZE / 2 + 1];
// 
// 		for (iLoop = 1; iLoop <= SUPPORT_VECTOR_SIZE / 2; iLoop++)
// 		{
// 			Parray_support[iLoop] = rootVertex->support[iLoop - 1];
// 			Qarray_support[iLoop] = rootVertex->support[SUPPORT_VECTOR_SIZE / 2 + iLoop - 1];
// 		}
// 
// 		delete [] Parray_support;
// 		delete [] Qarray_support;	

		delete [] fragSupportSum; delete [] exonSupportSum;

		return true;
	} 
	else
	{
		//extend children
		allChildrenAreLeaves = true;
		rootVertex->estimate_exonNum = 0;
		rootVertex->childNum = 0;
		rootVertex->junctionNum = 0;
		rootVertex->exonNum = 0;

		curEdge = rootVertex->child;
		while (curEdge != NULL)
		{
			(rootVertex->childNum)++;
			curVertex = curEdge->linkedVertex;
			if (countGTree(curVertex) == false)
			{
				allChildrenAreLeaves = false;
			}
			if (rootVertex->childType != 3)
			{
				rootVertex->junctionNum += curVertex->junctionNum;
				rootVertex->exonNum += curVertex->exonNum;

				if (curVertex->estimated == false)
				{
					rootVertex->estimate_exonNum += curVertex->exonNum;
				} 
				else
				{
					rootVertex->estimate_exonNum += 1; //representative exon
				}
			} 

			if (setting_junc_anno_provided && !curVertex->estimated && curVertex->junctionInRange)
			{
				curJunc = curVertex->junctionInRange->list;
				while (curJunc != NULL)
				{
					if (curJunc->junc->type == frag_junction && !curJunc->junc->in_junc_annotation)
					{
						rootVertex->has_novel_junction = true;
						curVertex->has_novel_junction = true;
					}
					curJunc = curJunc->next;
				}
			}

			curEdge = curEdge->next;
		}

		if (rootVertex->childType == 3)
		{
			juncCnt = 0;
			exonCnt = 0;
			fragCnt = 0;

			curJunc = rootVertex->junctionInRange->list;
			while (curJunc != NULL)
			{
				fragCnt++;
				if (curJunc->junc->type == frag_junction)
				{
					//alternative splice sites
					if (curJunc->junc->alter != NULL)
					{
						vertexForAlterSpliceSite(curJunc->junc, rootVertex);
					}

					juncCnt++;

// 					if (setting_junc_anno_provided && !curJunc->junc->in_junc_annotation)
// 						rootVertex->has_novel_junction = true;
				}
				else if (curJunc->junc->type == frag_exon || curJunc->junc->type == frag_retained_intron)
				{
					exonCnt++;
				}
				curJunc = curJunc->next;
			}
			rootVertex->junctionNum = juncCnt;
			rootVertex->exonNum = exonCnt;

			rootVertex->estimate_exonNum = exonCnt;
		}

		//		if (allChildrenAreLeaves == true && (rootVertex->childType == 2 || rootVertex->childType == 3))
		if (rootVertex->childType == 2 || rootVertex->childType == 3)
		{
			for (iLoop = 0; iLoop < SUPPORT_VECTOR_SIZE; iLoop++)
			{
				abundance_estimation(rootVertex, iLoop, MSE_estimation);
				rootVertex->MSE_estimation[iLoop] = MSE_estimation;
			}

			select_major_paths(rootVertex);

			// 			if (SUPPORT_VECTOR_SIZE > 1 && rootVertex->major_alter_paths_num > 1)
			// 			{
			// 				Parray_support = new double [rootVertex->major_alter_paths_num + 1];
			// 				Parray_proportion = new double [rootVertex->major_alter_paths_num + 1];
			// 				Qarray_support = new double [rootVertex->major_alter_paths_num + 1];
			// 				Qarray_proportion = new double [rootVertex->major_alter_paths_num + 1];
			// 
			// 				proportion_matrix = new double * [rootVertex->major_alter_paths_num + 1];
			// 				for (iLoop = 0; iLoop <= rootVertex->major_alter_paths_num; iLoop++)
			// 				{
			// 					proportion_matrix[iLoop] = new double[SUPPORT_VECTOR_SIZE + 1];
			// 				}
			// 
			// 				for (iLoop = 0; iLoop < rootVertex->major_alter_paths_num + 1; iLoop++)
			// 				{
			// 					Parray_support[iLoop] = 0.0;
			// 					Parray_proportion[iLoop] = 0.0;
			// 					Qarray_support[iLoop] = 0.0;
			// 					Qarray_proportion[iLoop] = 0.0;
			// 
			// 					for (jLoop = 0; jLoop <= SUPPORT_VECTOR_SIZE; jLoop++)
			// 					{
			// 						proportion_matrix[iLoop][jLoop] = 0.0;
			// 					}
			// 				}
			// 
			// 				tmp = 1;
			// 				curAlterPath = rootVertex->major_alter_paths;
			// 				while (curAlterPath != NULL)
			// 				{
			// 					for (iLoop = 0; iLoop < SUPPORT_VECTOR_SIZE/2; iLoop++)
			// 					{
			// 						Parray_support[tmp] += curAlterPath->support[iLoop];
			// 						Parray_proportion[tmp] += curAlterPath->proportion[iLoop];
			// 						Qarray_support[tmp] += curAlterPath->support[SUPPORT_VECTOR_SIZE / 2 + iLoop];
			// 						Qarray_proportion[tmp] += curAlterPath->proportion[SUPPORT_VECTOR_SIZE / 2 + iLoop];
			// 					}
			// 
			// 					for (iLoop = 1; iLoop <= SUPPORT_VECTOR_SIZE; iLoop++)
			// 					{
			// 						proportion_matrix[tmp][iLoop] = curAlterPath->proportion[iLoop - 1];
			// 					}
			// 
			// 					tmp++;
			// 					curAlterPath = curAlterPath->next;
			// 				}
			// 
			// 				totalSupport = 0.0;
			// 				for (iLoop = 1; iLoop <= rootVertex->major_alter_paths_num; iLoop++)
			// 				{
			// 					Parray_support[iLoop] = Parray_support[iLoop] / (SUPPORT_VECTOR_SIZE / 2);
			// 					Parray_proportion[iLoop] = Parray_proportion[iLoop] / (SUPPORT_VECTOR_SIZE / 2);
			// 					Qarray_support[iLoop] = Qarray_support[iLoop] / (SUPPORT_VECTOR_SIZE / 2);
			// 					Qarray_proportion[iLoop] = Qarray_proportion[iLoop] / (SUPPORT_VECTOR_SIZE / 2);
			// 
			// 					totalSupport += Parray_support[iLoop] + Qarray_support[iLoop];
			// 				}
			// 
			// // 				outputToPool(proportion_matrix, SUPPORT_VECTOR_SIZE, rootVertex->major_alter_paths_num);
			// // 
			// // 				rootVertex->JSDsqrt_mean = JSD_test(Parray_proportion, Qarray_proportion, rootVertex->major_alter_paths_num);
			// // 				rootVertex->JSDsqrt_combination = JSDsqrt_matrix_combination(proportion_matrix, SUPPORT_VECTOR_SIZE, rootVertex->major_alter_paths_num);
			// // 
			// // 				rootVertex->JSD_Pvalue_mean = JSD_significance_level(rootVertex->JSDsqrt_mean, rootVertex->major_alter_paths_num, totalSupport);
			// // 				rootVertex->JSD_Pvalue_combination = JSD_significance_level(pow(rootVertex->JSDsqrt_combination, 2), rootVertex->major_alter_paths_num, totalSupport);
			// // 
			// // 				rootVertex->divergence_stat = D_statistics(proportion_matrix, SUPPORT_VECTOR_SIZE, rootVertex->major_alter_paths_num);
			// 
			// 				if (SUPPORT_VECTOR_SIZE > 2)
			// 				{
			// // 					rootVertex->JSD_withinGroupRank_mean = withinMatrixRankTest_fullspace_JSD(proportion_matrix, SUPPORT_VECTOR_SIZE, rootVertex->major_alter_paths_num, rootVertex->JSDsqrt_mean, totalSupport);
			// // 					rootVertex->divergence_stat_rank = withinMatrixRankTest_fullspace_Dstatistics(proportion_matrix, SUPPORT_VECTOR_SIZE, rootVertex->major_alter_paths_num, rootVertex->divergence_stat);
			// 				} 
			// 				else
			// 				{
			// 					rootVertex->JSD_withinGroupRank_mean = 1;
			// 				}
			// 
			// 
			// /*				rootVertex->JSDsqrt_mean = sqrt(rootVertex->JSDsqrt_mean);*/
			// 
			// 				//				rootVertex->spearmanCorr_support = Spearman_correlation(Parray_support, Qarray_support, rootVertex->childNum);
			// 				//				rootVertex->spearmanCorr_proportion = Spearman_correlation(Parray_proportion, Qarray_proportion, rootVertex->childNum);
			// // 				rootVertex->pearsonCorr_support = Pearson_correlation(Parray_support, Qarray_support, rootVertex->major_alter_paths_num);
			// // 				rootVertex->pearsonCorr_proportion = Pearson_correlation(Parray_proportion, Qarray_proportion, rootVertex->major_alter_paths_num);
			// 
			// 				delete [] Parray_support;
			// 				delete [] Parray_proportion;
			// 				delete [] Qarray_support;
			// 				delete [] Qarray_proportion;
			// 
			// 				for (iLoop = 0; iLoop <= rootVertex->major_alter_paths_num; iLoop++)
			// 				{
			// 					delete [] proportion_matrix[iLoop];
			// 				}
			// 				delete [] proportion_matrix;
			// 			}
		}
		else if (rootVertex->childType == 1)
		{
			// 			Parray_support = new double [rootVertex->childNum + 1];
			// 			Qarray_support = new double [rootVertex->childNum + 1];
			// 
			// 			for (iLoop = 0; iLoop < rootVertex->childNum + 1; iLoop++)
			// 			{
			// 				Parray_support[iLoop] = 0.0;
			// 				Qarray_support[iLoop] = 0.0;
			// 			}
			// 
			// 			tmp = 1;
			// 			curEdge = rootVertex->child;
			// 			while (curEdge != NULL)
			// 			{
			// 				curVertex = curEdge->linkedVertex;
			// 
			// 				for (iLoop = 0; iLoop < SUPPORT_VECTOR_SIZE/2; iLoop++)
			// 				{
			// 					Parray_support[tmp] += curVertex->support[iLoop];
			// 					Qarray_support[tmp] += curVertex->support[SUPPORT_VECTOR_SIZE / 2 + iLoop];
			// 				}
			// 
			// 				tmp++;
			// 				curEdge = curEdge->next;
			// 			}
			// 
			// 			for (iLoop = 1; iLoop < rootVertex->childNum + 1; iLoop++)
			// 			{
			// 				Parray_support[iLoop] = Parray_support[iLoop] / (SUPPORT_VECTOR_SIZE / 2);
			// 				Qarray_support[iLoop] = Qarray_support[iLoop] / (SUPPORT_VECTOR_SIZE / 2);
			// 			}
			// 
			// //			rootVertex->spearmanCorr_support = Spearman_correlation(Parray_support, Qarray_support, rootVertex->childNum);
			// 			rootVertex->pearsonCorr_support = Pearson_correlation(Parray_support, Qarray_support, rootVertex->childNum);
			// 
			// 			delete [] Parray_support;
			// 			delete [] Qarray_support;
		}

		rootVertex->estimated = true;
		rootVertex->representative = makeRepresentative(rootVertex);


		// 		Parray_support = new double [SUPPORT_VECTOR_SIZE / 2 + 1];
		// 		Qarray_support = new double [SUPPORT_VECTOR_SIZE / 2 + 1];
		// 
		// 		for (iLoop = 1; iLoop <= SUPPORT_VECTOR_SIZE / 2; iLoop++)
		// 		{
		// 			Parray_support[iLoop] = rootVertex->support[iLoop - 1];
		// 			Qarray_support[iLoop] = rootVertex->support[SUPPORT_VECTOR_SIZE / 2 + iLoop - 1];
		// 		}
		// 
		// 		delete [] Parray_support;
		// 		delete [] Qarray_support;	


		delete [] fragSupportSum; delete [] exonSupportSum;

		return false;
	}
}
