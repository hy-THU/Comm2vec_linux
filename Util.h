#pragma once

#include <string>
#include <vector>
#include <map>
#include <utility>
#include <stdio.h>
#include <numeric>
#include <time.h>
#include <algorithm>

using std::string;
using std::vector;
using std::map;
using std::pair;


class Util
{
public:
    	static vector<string> StringTokenize(string line);
    	static vector<string> StringSplit(string line, char separator);
    	static bool FindNodeInVec(vector<int>::iterator begin, vector<int>::iterator end, int iNode);
    	template <typename T> static void Norm_2D(T** arr, int nrow, int ncol, bool sym);
    	template <typename T> static void Norm(T* arr, int len, bool sym);
    	static int  MultiSample(vector<double> prob);
    	static int  MultiSample(double* prob, int len);
    	template <typename T> static T** init_Array2D(int nrow, int ncol, bool sym);
    	template <typename T> static T** init_Array2D(int nrow, int ncol, T val);
    	template <typename T> static void format_Array2D(T** arr, int nrow, int ncol);
    	static pair<int*, double*>  alias_setup(double* prob, int K);
    	static int alias_draw( int* J, double* q, int K );
		template <typename T> static bool HaveIntersection(vector<T> a, vector<T> b);
};

template <typename T>  
T** Util::init_Array2D(int nrow, int ncol, bool sym)  
{
    srand( time(0) );
    int size = sizeof(T);  
    int point_size = sizeof(T*);
     
    T ** arr = (T **) malloc(point_size * nrow + size * nrow * ncol);  
    if (arr != NULL)  
    {     
        memset(arr, 0, point_size * nrow + size * nrow * ncol);  
        T *head = (T*)((long)arr + point_size * nrow);  
       	for (int i = 0; i < nrow; i++)
		{
			arr[i] = (T*)((long)head + i * ncol * size);
			for (int j = 0; j < ncol; j++)
				arr[i][j] = (T)((rand() / (RAND_MAX + 0.0) - 0.5) / ncol);
		}
    }
    else
    {
		fprintf(stderr, "Memory allocation failed.");
		exit(0);
    }
    //Util::Norm_2D <T> (arr, nrow, ncol, sym);
    return (T**)arr;  
}


template <typename T>  
T** Util::init_Array2D(int nrow, int ncol, T val)  
{  
    int size = sizeof(T);  
    int point_size = sizeof(T*);  
     
    T ** arr = (T **) malloc(point_size * nrow + size * nrow * ncol);  
    if (arr != NULL)  
    {     
        memset(arr, 0, point_size * nrow + size * nrow * ncol);
		T *head = (T*)((long)arr + point_size * nrow);
		if (val == 0)
			while (nrow--)  
				arr[nrow] = (T*)((long)head + nrow * ncol * size);
		else
		{
			for (int i = 0; i < nrow; i++)
			{
				arr[i] = (T*)((long)head + i * ncol * size);
				for (int j = 0; j < ncol; j++)
					arr[i][j] = val;
			}
		}
    }
    else
    {
		fprintf(stderr, "Memory allocation failed.");
		exit(0);
    }

    return (T**)arr;  
}

template <typename T>  
void Util::format_Array2D(T** arr, int nrow, int ncol)
{
    int size = sizeof(T);
    int point_size = sizeof(T*);
    T* head = (T*)((long)arr + point_size * nrow);
    memset(head, 0, size * nrow * ncol);
}

template<typename T>
void Util::Norm_2D(T** arr, int nrow, int ncol, bool sym)
{
    for (int i = 0; i < nrow; i ++)
    {
		Util::Norm<T>(arr[i], ncol, sym);
    }
}

template<typename T>
void Util::Norm(T* arr, int len, bool sym)
{
    T sum;
    T adj = 1.0/len/2;
    memset(&sum, 0, sizeof(T));
    for (int i = 0; i < len; i ++)
		sum += arr[i];
    if (sym)
		for (int i = 0; i < len; i ++)
		{
			arr[i] /= sum;
			arr[i] -= adj;
		}
    else
		for (int i = 0; i < len; i ++)
			arr[i] /= sum;
}

template<typename T>
bool Util::HaveIntersection(vector<T> a, vector<T> b)
{
	typename vector<T>::iterator ap = a.begin();
	typename vector<T>::iterator bp = b.begin();
	while (ap != a.end() && bp != b.end())
	{
		if (*ap < *bp)
			ap ++;
		else if ( *ap > *bp )
			bp ++;
		else
			return true;
	}
	return false;
}