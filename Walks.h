#pragma once

#include "Util.h"
#include "DataSet.h"
#include "Config.h"


#include <vector>
#include <map>
#include <algorithm>
#include <cmath>
#include <time.h>
#include <numeric>
#include <pthread.h>


//#pragma comment(lib,"x86_pthread/pthreadVC2.lib")

using std::vector;
using std::map;


class Walks
{
public:
    //DataSet* data;
    //Config* conf;



	//vector<NodeIDType> OneWalk();
	//vector< vector<NodeIDType> > SimuWalks();

	//NodeIDType** OneWalk(NodeIDType startnode);
	static NodeIDType** SimuWalks(DataSet* data, Config* conf);
	static void* SimuWalksThread(void* args);


	static int node_count_actual;
    //template <typename T> T** init_Array2D(int row, int col);
	//template <typename T> T** init_Array2D(int row, int col, T val);

	//vector < vector<int> >   walks;


	//Walks(){ node_count = 0; };
	//~Walks(){};

};
