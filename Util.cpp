#pragma once
#include "Util.h"

//#include <algorithm>
#include <sstream>

using std::istringstream;
using std::make_pair;


vector<string> Util::StringTokenize(string line)
{
    	istringstream   strin(line);
    	vector<string>  result;
    	string          token;

    	while (strin >> token)
        	result.push_back(token);

    	return result;
}

vector<string> Util::StringSplit(string line, char separator)
{
    	vector<string>  result;
    	line += separator;

    	int p = 0;
    	for (int i = 0; i < line.length(); i ++)
        	if (line[i] == separator)
        	{
            		if (i - p > 0) result.push_back( line.substr(p, i-p) );
            		p = i + 1;
        	}

    	return result;
}



bool Util::FindNodeInVec(vector<int>::iterator begin, vector<int>::iterator end, int iNode)
{
    	while (begin!=end)
    	{
        	if (*begin==iNode)
           		break;
        	else
            		begin++;
    	}
 
    	if (begin!=end)
        	return true;
   	else
        	return false;
}




int Util::MultiSample(vector<double> prob)
{
	if (prob.size() == 0) 
	{
		fprintf(stderr, "Prob can not be empty.");
		exit(0);
	}
	srand( time(0) );
    	double sum = accumulate( prob.begin(), prob.end(), 0.0 );
	if (sum == 0.0)
	{
		fprintf(stderr, "Sum can not be zero.");
		exit(0);
	}
    	for (int i = 0; i < prob.size(); i++)
    	{
		prob[i] /= sum;
    	}

    	double P = rand() / (RAND_MAX + 0.0);
    	double p = 0.0;

    	for (int i = 0; i < prob.size(); i++)
    	{
		p += prob[i];
		if (P <= p)
			return i;
    	}

}


int Util::MultiSample(double* prob, int len)
{
	if (len <= 0)
	{
		fprintf(stderr, "Prob can not be empty.");
		exit(0);
	}
	srand( time(0) );
    	double sum = 0.0;
	for (int i = 0; i < len; i++)
		sum += prob[i];
    for (int i = 0; i < len; i++)
    {
		prob[i] /= sum;
    }

    double P = rand() / (RAND_MAX + 0.0);
    double p = 0.0;

    for (int i = 0; i < len; i++)
    {
		p += prob[i];
		if (P <= p)
			return i;
    }

}


pair<int*, double*>  Util::alias_setup(double* probs, int K)
{
	
	double* q;
	int* J;

	
	q = (double*)calloc(K, sizeof(double));
	J = (int*)calloc(K, sizeof(int));
	
	vector<int> smaller;
	vector<int> larger;
	for (int k = 0; k < K; k++)
	{
		q[k] = K * probs[k];
		if (q[k] < 1.0)
			smaller.push_back(k);
		else
			larger.push_back(k);
	}

	int small, large;
	while (smaller.size() > 0 && larger.size() > 0)
	{
		small = *(smaller.end() - 1);
		large = *(larger.end() - 1);
		smaller.pop_back();
		larger.pop_back();

		J[small] = large;
		q[large] -= (1.0 - q[small]);

		if (q[large] < 1.0)
			smaller.push_back(large);
		else
			larger.push_back(large);
	}

	return make_pair(J, q);

}


int Util::alias_draw( int* J, double* q, int K )
{
	srand( time(0) );
	int kk = rand() * (K-1) / RAND_MAX;
	if (rand()/(RAND_MAX + 0.0) < q[kk])
		return kk;
	else
		return J[kk];
}
