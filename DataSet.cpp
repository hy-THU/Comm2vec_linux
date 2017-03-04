#pragma once
#include "DataSet.h"
#include <algorithm>
#include <cstdio>


#define MAX_BUF_SIZE 65535
#define STEP_CMTY_COUNT 100000

void DataSet::LoadData(const Config* conf)
{
	time_t start_time, end_time;
	printf("Loading data ...\n");
	start_time = time(NULL);
	LoadCommData(conf -> commfile.c_str());
	LoadTopoData(conf -> topofile.c_str(), conf -> directed, conf -> weighted, conf -> walkpattern);
	AddProbs(conf -> walkpattern, conf -> p, conf -> q);
	end_time = time(NULL);
	printf("Data loaded.  Time costed is %d seconds\n", (end_time - start_time));
}

void DataSet::LoadCommData(const char* commfile)
{
	comm = (CommNode**)calloc(STEP_CMTY_COUNT, sizeof(CommNode*));
	int comm_size = STEP_CMTY_COUNT;
	CommIDType comm_id;
	CommWtType wt;
	char buf[MAX_BUF_SIZE];
	char* eof;

	vector<string> tokens;
	vector<string> comm_wt;

	FILE *fin = fopen(commfile, "r");
	NodeIDType nodeID = 0;
	num_comm = 0;
	int flag = 0;
	while ( eof = fgets(buf, MAX_BUF_SIZE, fin) )
	{

		if (strlen(buf) == 1) continue;//empty line
		if ( flag )
		{
			tokens = Util::StringTokenize(buf);
			if ( tokens.size() != 2 )
			{
				fprintf( stderr, "Wrong data format!\n");
				exit(0);
			}
			num_node = atoi(tokens[0].c_str());
			num_comm = atoi(tokens[1].c_str());
			flag = 0;
		}

		DataNode* temp_node = new DataNode();
		temp_node -> ID = nodeID;
		temp_node -> degree = 0;

		// Parse detail information
		tokens = Util::StringTokenize(buf);

		for( int i = 0 ; i < tokens.size() ; i++ )
		{
			comm_wt = Util::StringSplit(tokens[i], ':');
			if ( comm_wt.size() != 2 )
				continue;
			comm_id = atoi(comm_wt[0].c_str());
			wt = atof(comm_wt[1].c_str());
			temp_node -> comm_list.push_back( comm_id );
			temp_node -> weight.push_back( wt );
			if(comm_size <= comm_id)
			{
				CommNode** temp_comm = (CommNode**)realloc(comm, (comm_size + STEP_CMTY_COUNT) * sizeof(CommNode*));
				if(temp_comm == NULL)
				{
					fprintf(stderr, "Memory allocation failed.");
					exit(0);
				}
				comm = temp_comm;
				memset(comm + comm_size, 0, STEP_CMTY_COUNT * sizeof(CommNode*));
				comm_size += STEP_CMTY_COUNT;
			}
			if(comm[comm_id] == NULL)
				comm[comm_id] = new CommNode();
			comm[comm_id] -> node_list.push_back(nodeID);
			comm[comm_id] -> weight.push_back(wt);
			if (num_comm < comm_id)
				num_comm = comm_id;
		}
		temp_node -> n_comm = temp_node -> comm_list.size();
		node.push_back( temp_node );
		nodeID ++;
	}
	num_comm++;
	if ( num_node != nodeID )
		num_node = nodeID;

	fclose(fin);
	printf("#Nodes:%d\t#Communities:%d\n", num_node, num_comm);
}

void DataSet::LoadTopoData(const char* topofile, const int directed, const int weighted, const int walkpattern)
{
	char          buf[MAX_BUF_SIZE];
	char*         eof;
	NodeIDType    a,b; //endpoints of an edge
	num_edge = 0;


	vector<string>  tokens;

	FILE        *fin = fopen(topofile, "r");

	while ( eof = fgets(buf, MAX_BUF_SIZE, fin) )
	{

		if (strlen(buf) == 1) continue;//empty line

		// Parse detail information
        	tokens = Util::StringTokenize(buf);
		if (tokens[0] == "#") continue;

		if ( tokens.size() != 2 && tokens.size() != 3 )
		{
			fprintf( stderr, "Wrong topology data format!\n");
			exit( 0 );
		}

		a = atoi( tokens[0].c_str() );
		b = atoi( tokens[1].c_str() );
		if (a < 0 || a >= num_node || b < 0 || b >= num_node) continue; //out of range

		//ignore repeated edge
		NbrType::iterator iter = find( node[a]->nbr.begin(), node[a]->nbr.end(), b );
		if (iter != node[a] -> nbr.end())
			continue;
		if (directed)
		{
			if (weighted)
			{
				node[a]->nbr.push_back( b );
				node[a]->nbr_wt.push_back( atof(tokens[2].c_str()) );
			}
			else
			{
				node[a]->nbr.push_back( b );
				node[a]->nbr_wt.push_back( 1.0 );
			}
			//node[a] ->nbr_dict.insert( make_pair(node[a] -> degree, b) );
			node[a] -> degree ++;
			//edge.push_back( EdgeType(a, b) );
		}
		else
		{
			if (weighted)
			{
				//node[a]->nbr.insert( make_pair(b, atof(tokens[2].c_str())) );
				node[a]->nbr.push_back( b );
				node[a]->nbr_wt.push_back( atof(tokens[2].c_str()) );
				node[b]->nbr.push_back( a );
				node[b]->nbr_wt.push_back( atof(tokens[2].c_str()) );
				//node[b]->nbr.insert( make_pair(a, atof(tokens[2].c_str())) );
			}
			else
			{
				node[a]->nbr.push_back( b );
				node[a]->nbr_wt.push_back( 1.0 );
				node[b]->nbr.push_back( a );
				node[b]->nbr_wt.push_back( 1.0 );
			}
			//node[a] ->nbr_dict.insert( make_pair(node[a] -> degree, b) );
			node[a] -> degree ++;
			//node[b] ->nbr_dict.insert( make_pair(node[b] -> degree, a) );
			node[b] -> degree ++;

			//edge.push_back( EdgeType (a, b) );
			//edge.push_back( EdgeType (b, a) );
		}
		if (walkpattern)
		{
			edge.push_back( EdgeType (a, b) );
			if (!directed)
				edge.push_back( EdgeType (b, a) );
		}
		num_edge ++;

	}

	fclose(fin);

	/*
	num_edge = 0;
	for (int i = 0; i < num_node; i++)
		num_edge += node[i] -> degree;
	if (directed == 0)
		num_edge /= 2;*/

}

void DataSet::AddProbs(const int walkpattern, const double p, const double q)
{
	/*
	double sum;
	for (int i = 0; i < num_node; i++ )
	{
		sum = accumulate( node[i]->nbr_wt.begin(), node[i]->nbr_wt.end(), 0.0 );
		for (int j = 0; j < node[i]->degree; j++)
			node[i] -> probs.push_back( node[i]->nbr_wt[j] / sum );
	}*/

	if (walkpattern)
	{
		for (vector<EdgeType>:: iterator iter = edge.begin(); iter != edge.end(); iter++)
		{
			edge_probs.insert( make_pair(*iter, GetEdgeProbs(*iter, p, q)) );
		}
	}
}

vector<double> DataSet::GetEdgeProbs(const EdgeType curedge, const double p, const double q)
{
	vector<double> probs;
	NodeIDType a = curedge.first;
	NodeIDType b = curedge.second;
	NodeIDType target;
	for (int i = 0; i < node[b] -> nbr.size(); i ++)
	{
		target = node[b] -> nbr[i];
		if (target == a)
			probs.push_back(node[b] -> nbr_wt[i]/p);
		else if ( Util::FindNodeInVec(node[target]-> nbr.begin(), node[target] -> nbr.end(), a) || Util::FindNodeInVec(node[target]-> nbr.begin(), node[target] -> nbr.end(), a) )
			probs.push_back(node[b] -> nbr_wt[i]);
		else
			probs.push_back(node[b] -> nbr_wt[i]/q);
	}
	return probs;
}

