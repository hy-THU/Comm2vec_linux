#pragma once

#include "Util.h"
#include "Config.h"

#include <string>
#include <vector>
#include <map>
#include <numeric>

using namespace std;

typedef int										NodeIDType;
typedef int                                     CommIDType;
typedef double                                  CommWtType;
//typedef int										ActIDType;
typedef pair< NodeIDType, NodeIDType >			EdgeType;
typedef vector<NodeIDType> 		            	NbrType;
//typedef double*                                 AttrType;
//typedef vector< pair<NodeIDType, ActIDType> >   ActType;
//typedef map<CommIDType, CommWtType>             CommType;



class DataNode
{
public:

    NodeIDType            ID;

    int                   n_comm;
    //CommType              comm;
    vector<CommIDType>    comm_list;
    vector<CommWtType>    weight;

    int                   degree;
    NbrType               nbr;
    vector<double>        nbr_wt;
    //vector<double>        probs;
    //map<int, NodeIDType>  nbr_dict;
};

struct CommNode
{
    vector<NodeIDType>    node_list;
    vector<CommWtType>    weight;
};


class DataSet
{
public:

    int num_node;
    int num_edge;
    int num_comm;
    //int num_attr;
    //int num_act;


    vector<DataNode*>	       node;
    CommNode**                 comm;
    //vector<vector<CommIDType>> comm;
    vector<EdgeType>         edge;
	map< EdgeType, vector<double> >	edge_probs;


    void LoadData(const Config* conf);
    void LoadTopoData(const char* topofile, const int directed, const int weighted, const int walkpattern);
    void LoadCommData(const char* commfile);
	void AddProbs(const int walkpattern, const double p, const double q);
    vector<double> GetEdgeProbs(const EdgeType curedge, const double p, const double q);

    //void LoadActData(const char* actfile);
    //void DumpData(const char* output_file);
    //void ReturnData(const char* output_file);

    DataSet() { }
    ~DataSet()
    {
        for (int i = 0; i < node.size(); i++)
            delete( node[i] );
    }

};
