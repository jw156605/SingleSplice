#include <boost/random/bernoulli_distribution.hpp>
#include <boost/random/chi_squared_distribution.hpp>
#include <boost/random/gamma_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <iostream>
#include <fstream>
#include <iterator>
#include <algorithm>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <vector>
using namespace std;

void read_expression_levels(vector<double> & A, vector<double> & B, string expr_file);
vector<double> ratio_variance_from_noise(const vector<double> & A,
					 const vector<double> & B,
					 double a0, double a1,
					 double b0, double b1, int n);
void sample_gamma_with_dropouts(double p, double k, double theta, int n, 
				vector<double> & A, boost::random::mt19937 & generator);
vector<string> split(const string& s, const string& delim, const bool keep_empty=false);
void read_expression_levels(vector<double> & A, vector<double> & B, string expr_file);
void get_ratio(const vector<double> & A, const vector<double> & B, vector<double> &S);
double mean(const vector<double> & A);
double var(const vector<double> & A);

int main(int argc, char* argv[])
{
  string expr_file(argv[1]);
  double a0 = atof(argv[2]);
  double a1 = atof(argv[3]);
  double b0 = atof(argv[4]);
  double b1 = atof(argv[5]);
  int n = atoi(argv[6]);
  vector<double> A, B;
  read_expression_levels(A,B,expr_file);

  vector<double> result = ratio_variance_from_noise(A,B,a0,a1,b0,b1,n);
  for (int i = 0; i < result.size()-1; ++i)
    {
      cout << result[i] << " ";
    }
  cout << result[result.size()-1] << endl;
}

vector<string> split(const string& s, const string& delim, const bool keep_empty) {
    vector<string> result;
    if (delim.empty()) {
        result.push_back(s);
        return result;
    }
    string::const_iterator substart = s.begin(), subend;
    while (true) {
        subend = search(substart, s.end(), delim.begin(), delim.end());
        string temp(substart, subend);
        if (keep_empty || !temp.empty()) {
            result.push_back(temp);
        }
        if (subend == s.end()) {  
	  break;
        }
        substart = subend + delim.size();
    }
	result.pop_back();
    return result;
}

void read_expression_levels(vector<double> & A, vector<double> & B, string expr_file)
{
  ifstream fin;
  fin.open(expr_file.c_str());
  if (fin.fail())
    {
      cout << "Expression file failed to open" << endl;
      return;
    }
  string line;
  getline(fin, line);
  vector<string> levels = split(line,",");
  for (int i = 0; i < levels.size(); ++i)
    {
      A.push_back(atof(levels[i].c_str()));
    }
  getline(fin, line);
  levels = split(line,",");
  for (int i = 0; i < levels.size(); ++i)
    {
      B.push_back(atof(levels[i].c_str()));
    }
}

void get_ratio(const vector<double> & A, const vector<double> & B, vector<double> &S)
{
  for (int i = 0; i < A.size(); ++i)
    {
      if (A[i]+B[i] == 0)
	{
	  S[i] = 0.5;
	}
      else
	{
	  S[i] = A[i]/(A[i]+B[i]);
	}
    }
}

double mean(const vector<double> & A)
{
  double mu = 0;
  for (int i = 0; i < A.size(); ++i)
    {
      mu += A[i];
    }
  return mu / (double)A.size();
}

double var(const vector<double> & A)
{
  double sigma_sq = 0;
  double mu = mean(A);
  for (int i = 0; i < A.size(); ++i)
    {
      sigma_sq += (A[i]-mu)*(A[i]-mu);
    }
  return (sigma_sq / (double)(A.size()-1));
}

void sample_gamma_with_dropouts(double p, double k, double theta, int n, vector<double> & A, boost::random::mt19937 & generator)
{
  boost::random::bernoulli_distribution<>  dropout(p);
  boost::random::gamma_distribution<> noise_model(k, theta);
  for (int i = 0; i < n; ++i)
    {
      if (dropout(generator)) //Sample from Bernoulli(p) to determine if this is a dropout event
	{
	  A[i] = 0;
	}
      else                    //Not a dropout--sample from Gamma(k,theta)
	{
	  A[i] = noise_model(generator);
	}
    }
}

vector<double> ratio_variance_from_noise(const vector<double> & A,
					 const vector<double> & B,
					 double a0, double a1,
					 double b0, double b1, int n)
{
  int r = A.size();
  double mu1 = mean(A);
  double mu2 = mean(B);
  double p1 = 1/(1+exp(-(b0+b1*log(mu1+1))));
  double p2 = 1/(1+exp(-(b0+b1*log(mu2+1))));
  vector<double> S(r); //S stores ratio of A and B 
  get_ratio(A,B,S);
  double thresh = var(S);

  boost::random::mt19937 generator(time(0)); //Mersenne twister 19937 seeded with clock
  boost::random::chi_squared_distribution<> rchisq(r-1); //chi-squared with r-1 degrees of freedom

  vector<double> null_variances(n);
  for (int i = 0; i < n; ++i)
      {
	double cv2 = rchisq(generator)/(r-1)*(a1/mu1 + a0);
	//double cv2 = a1/mu1 + a0;
	double var1 = cv2*(mu1*mu1);
	double k1 = (mu1*mu1)/(var1*(1-p1)-p1*(mu1*mu1));
	double theta1 = (var1*(1-p1)-p1*(mu1*mu1))/(mu1*(1-p1));
	cv2 = rchisq(generator)/(r-1)*(a1/mu2 + a0);
	//cv2 = a1/mu2 + a0;
	double var2 = cv2*(mu2*mu2);
	double k2 = (mu2*mu2)/(var2*(1-p2)-p2*(mu2*mu2));
	double theta2 = (var2*(1-p2)-p2*(mu2*mu2))/(mu2*(1-p2));
	vector<double> A_sim(r), B_sim(r), S_sim(r);
	sample_gamma_with_dropouts(p1,k1,theta1,r,A_sim,generator);
	sample_gamma_with_dropouts(p2,k2,theta2,r,B_sim,generator);
	get_ratio(A_sim,B_sim,S_sim);
	null_variances[i] = var(S_sim);
      }
  int num_greater = 0;
  for (int i = 0; i < n; ++i)
    {
      if (null_variances[i] > thresh)
	++num_greater;
    }
  double p_val = num_greater/ (double)n;
  double exp_var = mean(null_variances);
  vector<double> result;
  result.push_back(mu1);
  result.push_back(mu2);
  result.push_back(exp_var);
  result.push_back(thresh);
  result.push_back(p_val);
  return (result);
}
