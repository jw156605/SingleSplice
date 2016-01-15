
#include "cut_exon_bound.h"
#include "wavelet2s.h"

const double MAT_INF = 2E20;

void center_vec(double *vec, long N)
{
	// calculate the vector mean
	double mean_vec = 0;
	for (long iLoop = 1; iLoop <= N; ++iLoop)
	{
		mean_vec += vec[iLoop];
	}
	mean_vec /= N;

	for (long iLoop = 1; iLoop <= N; ++iLoop)
	{
		vec[iLoop] -= mean_vec;
	}

	return;
}

// cumulative sum
void cumsum(double *vec, double *cumsum, long N)
{
	cumsum[0] = 0;
	for (long iLoop = 1; iLoop <= N; ++iLoop)
	{
		cumsum[iLoop] = cumsum[iLoop-1] + vec[iLoop];
	}
	return;
}

double calc_min(double *vec, long N, long &min_index)
{
	double min_value = vec[1];
	min_index = 1;
	for (long iLoop = 2; iLoop <= N; ++iLoop)
	{
		if (vec[iLoop] < min_value)
		{
			min_value = vec[iLoop];
			min_index = iLoop;
		}
	}
	return min_value;
}

// cut exons based on the wavelet decomposition
// input:	y		- wavelet signal
//			k		- maximum number of change points 
//			d		- minimum gap between two change points
//			num_cut - specified number of cut points
// output:	cut_pt		- the list of change points
void cut_exon(vector <double> y_vec, int K, int d, int num_cut, vector <long> &cut_pt)
{
	if (K <= 0)
		return;
	if (d <= 0)
		d = 0;

	long N = y_vec.size(), iLoop, jLoop;
	double *y = new double [N+1];
	for (iLoop = 1; iLoop <= N; ++iLoop)
		y[iLoop] = y_vec[iLoop-1];

	// center signal
	center_vec(y, N);

	// increment max_cut_pt
	++K;

	// t_est computation using dynamic programming
	double *y_squared = new double [N+1];
	for (iLoop = 1; iLoop <= N; ++iLoop)
		y_squared[iLoop] = y[iLoop] * y[iLoop];
	double **matD = new double* [N+1];
	for (iLoop = 1; iLoop <= N; ++iLoop)
		matD[iLoop] = new double [N+1];

	for (iLoop = 1; iLoop <= N; ++iLoop)
		for (jLoop = 1; jLoop <= N; ++jLoop)
			matD[iLoop][jLoop] = MAT_INF;

	for (iLoop = 1; iLoop <= N-d; ++iLoop)
	{
		double *vi = new double [N-iLoop+2];
		double *csum = new double [N-iLoop+2];

		// calculate the cumulative sum
		for (jLoop = 1; jLoop <= N-iLoop+1; ++jLoop)
			vi[jLoop] = y_squared[jLoop+iLoop-1];
		cumsum(vi, csum, N-iLoop+1);

		for (jLoop = 1; jLoop <= N-iLoop+1; ++jLoop)
		{
			if (csum[jLoop] > 1E-20 && jLoop < MAT_INF)
			{
				vi[jLoop] = jLoop * log( csum[jLoop]/jLoop );
			} 
			else
			{
				vi[jLoop] = MAT_INF;
			}
		}

		for (jLoop = 1; jLoop <= N-iLoop-d+2; ++jLoop)
		{
			matD[iLoop][iLoop+d-2+jLoop] = vi[d-1+jLoop];
		}

		delete [] vi;
		delete [] csum;
	}

	double **I = new double* [K+1];
	for (iLoop = 1; iLoop <= K; ++iLoop)
		I[iLoop] = new double [N+1];
	for (iLoop = 1; iLoop <= K; ++iLoop)
		for (jLoop = 1; jLoop <= N; ++jLoop)
			I[iLoop][jLoop] = 0;

	for (iLoop = 1; iLoop <= N; ++iLoop)
		I[1][iLoop] = matD[1][iLoop];

	long **t = new long* [K+1];
	for (iLoop = 1; iLoop <= K; ++iLoop)
		t[iLoop] = new long [N+1];
	for (iLoop = 1; iLoop <= K; ++iLoop)
		for (jLoop = 1; jLoop <= N; ++jLoop)
			t[iLoop][jLoop] = 0;

	if (K > 2)
	{
		for (long k = 2; k <= K-1; ++k)
		{
			for (iLoop = k; iLoop <= N; ++iLoop)
			{
				double *tmp_vec = new double [iLoop];
				for (jLoop = 1; jLoop <= iLoop-1; ++jLoop)
				{
					tmp_vec[jLoop] = I[k-1][jLoop] + matD[jLoop+1][iLoop];
				}
				I[k][iLoop] = calc_min(tmp_vec, iLoop-1, t[k-1][iLoop]);

				delete [] tmp_vec;
			}
		}
	}

	long **t_est = new long* [K+1];
	for (iLoop = 1; iLoop <= K; ++iLoop)
		t_est[iLoop] = new long [K+1];
	for (iLoop = 1; iLoop <= K; ++iLoop)
		for (jLoop = 1; jLoop <= K; ++jLoop)
			t_est[iLoop][jLoop] = 0;
	for (iLoop = 1; iLoop <= K; ++iLoop)
		t_est[iLoop][iLoop] = N;

	double *tmp_vec = new double [N];
	for (jLoop = 1; jLoop <= N-1; ++jLoop) { tmp_vec[jLoop] = I[K-1][jLoop] + matD[jLoop+1][N];}
	I[K][N] = calc_min(tmp_vec, N-1, t[K-1][N]);
	delete [] tmp_vec;

	for (iLoop = 2; iLoop <= K; ++iLoop)
	{
		for (jLoop = iLoop-1; jLoop >= 1; --jLoop)
		{
			long col = t_est[iLoop][jLoop+1];
			if (col > 0)
				t_est[iLoop][jLoop] = t[jLoop][col];
		}
	}

	for (iLoop = 1; iLoop <= num_cut+1; ++iLoop) // output one more, which is the total length, for convenience in the calling function
	{
		cut_pt.push_back(t_est[num_cut+1][iLoop]);
	}

// 	for (iLoop = 1; iLoop <= K; ++iLoop)
// 	{
// 		for (jLoop = 1; jLoop <= K; ++jLoop)
// 		{
// 			cout << t_est[iLoop][jLoop] << " ";
// 		}
// 		cout << endl;
// 	}


	delete [] y;
	delete [] y_squared;
	for (iLoop = 1; iLoop <= N; ++iLoop)
	{
		delete [] matD[iLoop];		
	}
	for (iLoop = 1; iLoop <= K; ++iLoop)
	{
		delete [] I[iLoop];
		delete [] t[iLoop];
		delete [] t_est[iLoop];
	}
	delete [] matD;
	delete [] I;
	delete [] t;
	delete [] t_est;

	return;
}

void find_exon_cut_point(vector<double> &signal, int num_cut_pt, vector <long> &cut_pt) 
{
	// 	cout << "********J- LEVEL DISCRETE WAVELET TRANSFORM IMPLEMENTATION*********" << endl; // prints
	// 	cout << "This program accepts signal from the user in a file format " << endl;
	// 	cout << "and performs Discrete Wavelet Transform with specified   " << endl;
	// 	cout << "wavelet. "                                               << endl;
	// 	cout << "                                                             " << endl;
	// 	cout << " The Following Wavelets are in the Database:                 " << endl;
	// 	cout << " haar, db1, db2, db3, db4, db5, db6, db7, db8, db9, db10,  "   << endl;
	// 	cout << " db11, db12, db13, db14, db15.                               "  << endl;
	// 	cout << " bior1.1, bio1.3, bior1.5, bior2.2, bior2.4,bior2.6,bior2.8, " << endl;
	// 	cout << " bior3.1, bior3.3, bior3.5, bior3.7, bior3.9, bior4.4,"        << endl;
	// 	cout << " bior5.5, bior6.8."                                            << endl;
	// 	cout << " coif1, coif2, coif3, coif4, coif5."                           << endl;
	// 	cout << "Please Enter the Wavelet Name at the Prompt( No quotes)     :" << endl;

	string nm = "db3"; // nm will store the name of Wavelet Family
	int J = 1;

	vector<double> dwt_output, flag;

	// perform J-Level DWT
	vector<int> length;

	dwt_sym(signal, J, nm, dwt_output,flag,length);
	
	//Perform J-Level IDWT
	vector<double> output;

	for (unsigned int i = 0; i < length[0]; ++i)
	{
		dwt_output[i] = 0;
	}
	idwt_sym(dwt_output, flag,nm,output,length);

	cut_exon(output, 5, 100, num_cut_pt, cut_pt);

	return;
}

// exon_code: 1. transcription start exon; 2. transcription end exon; 3. mixed genes
bool refine_exon_bound(vector<double> &signal, int num_cut_pt, int exon_code, vector<long> &cut_pt, vector<double> &avg_coverage)
{
	bool is_cutting = false;
	cut_pt.clear(); avg_coverage.clear();

	vector <long> range_bound;
	range_bound.push_back(0);

	find_exon_cut_point(signal, num_cut_pt, cut_pt);

	for (int range_index = 0; range_index <= num_cut_pt; ++range_index)
	{
		avg_coverage.push_back(0.0);
		range_bound.push_back(cut_pt[range_index]-1);
		for (long pos = range_bound[range_index]; pos <= range_bound[range_index+1]; ++pos)
		{
			avg_coverage[range_index] += signal[pos];
		}
		avg_coverage[range_index] /= range_bound[range_index+1] - range_bound[range_index] + 1;
	}

	if (exon_code == 1)
	{
		if (avg_coverage[1] / avg_coverage[0] > 10)
			is_cutting = true;
	}
	else if (exon_code == 2)
	{
		if (avg_coverage[0] / avg_coverage[1] > 10)
			is_cutting = true;
	}
	else if (exon_code == 3)
	{
		if ((avg_coverage[0] / avg_coverage[1] > 2 && avg_coverage[2] / avg_coverage[1] > 2) 
			&& (avg_coverage[0] / avg_coverage[1] > 5 || avg_coverage[2] / avg_coverage[1] > 5))
			is_cutting = true;
	}

	return is_cutting;
}