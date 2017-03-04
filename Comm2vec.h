#pragma once
#pragma comment(lib,"x86_pthread/pthreadVC2.lib")

#include "Util.h"
#include "DataSet.h"
#include "Config.h"
#include "Walks.h"


#include <vector>
#include <map>
#include <algorithm>
#include <cmath>
#include <time.h>
#include <numeric>
#include <pthread.h>



using std::vector;
using std::map;

#define MAX_EXP 6
#define EXP_TABLE_SIZE 1000


class Comm2Vec
{
public:
    DataSet* data;
    Config*  conf;
	
	NodeIDType** walks;
	int* J;
	double* q;
	double* node_cdf;
	int node_count_actual;

	typedef struct
	{
		pthread_mutex_t mutex;
		double** vector;
	}Comm_Vec;
	typedef struct
	{
		pthread_mutex_t mutex;
		double** vector;
	}Node_Vec;
	typedef struct
	{
		pthread_mutex_t mutex;
		double** vector;
	}Node_Comm_Vec;
	typedef struct
	{
		pthread_mutex_t mutex;
		double** vector;
	}Syn1neg_Node;
	Comm_Vec comm_vec;
	Node_Vec node_vec;
	Node_Comm_Vec node_comm_vec;
	Syn1neg_Node syn1neg_node;

	//double** comm_vec;
	//double** node_vec;
	//double** node_comm_vec;
	//double** syn1neg_node;

	double* exptable;

	void Train();
	static void* TrainThread(void* id);
	static void* TrainThreadSecond(void* id);
	static void* TrainThreadThird(void* id);
	static void* TrainThreadForth(void* id);
	static void* TrainThreadFifth(void* id);
	void InitNet();
	void UpdateNodeCommVec();
	static void UpdateNodeCommVec(Comm2Vec* self, vector<CommIDType> &node_need_update_comm);
	void MakeExptable();
	void Output();
	void GenJq();
	int NodeSampling();
	void NormalizeL2();

	void Do();
	
    template <typename T> T** init_Array2D(int row, int col);
	template <typename T> T** init_Array2D(int row, int col, T val);



    Comm2Vec(Config* conf){ this->conf = conf; };
	~Comm2Vec()
	{
		delete(data);
		delete(conf);
	};
};
