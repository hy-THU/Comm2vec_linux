#pragma once

#include "Walks.h"

int Walks::node_count_actual = 100;
NodeIDType** Walks::SimuWalks(DataSet* data, Config* conf)
{
	//printf("%d\n", Walks::node_count_actual);
	time_t start_time, end_time;
	printf("Walking ...\n");
	start_time = time(NULL);
	NodeIDType** walks = Util::init_Array2D<NodeIDType> ( conf -> numwalks * data -> num_node, conf -> walklength, 0);
	
	node_count_actual = 0;
	struct argsformat
	{
		int id;
		DataSet* data;
		Config* conf;
		NodeIDType** walks;
	};
	argsformat* args = (argsformat*)calloc(conf -> walkthreads, sizeof(argsformat));
	pthread_t *pt = (pthread_t *)malloc( conf -> walkthreads * sizeof(pthread_t));
	if (args == NULL || pt == NULL)
	{
		fprintf(stderr, "Memory allocation failed.");
		exit(0);
	}
	for (int i = 0; i < conf->walkthreads; i++)
	{
		args[i].id = i;
		args[i].data = data;
		args[i].conf = conf;
		args[i].walks = walks;
		pthread_create(&pt[i], NULL, Walks::SimuWalksThread, (void *)&args[i]); 
	}
	for (int i = 0; i < conf->walkthreads; i++) pthread_join(pt[i], NULL);
	
	end_time = time(NULL);
	printf("\nWalk finished.  Time costed is %d seconds\n", (end_time - start_time));
	return walks;
	
}

void* Walks::SimuWalksThread(void* args)
{
	NodeIDType curnodeID, prevnodeID, nextnodeID;
	DataNode* curnode;
	struct argsformat
	{
		int id;
		DataSet* data;
		Config* conf;
		NodeIDType** walks;
	};
	int id = ((argsformat*)args) -> id;
	DataSet* data = ((argsformat*)args) -> data;
	Config* conf = ((argsformat*)args) -> conf;
	NodeIDType** walks = ((argsformat*)args) -> walks;

	int start = (double)id / conf -> walkthreads * data -> num_node;
	int end = ((double)id + 1.0) / conf -> walkthreads * data -> num_node;

	for (NodeIDType i = start; i < end; i++)
	{
		if (i - start > 100)
		{
			Walks::node_count_actual += i - start;
			start = i;
			printf("\r	Progress: %.2f%%", (Walks::node_count_actual + 1) / double(data -> num_node) * 100);
			fflush(stdout);
		}
		for (int j = 0; j < conf -> numwalks; j++)
		{
			walks[j * data -> num_node + i][0] = i;
			curnodeID = i;
			if(conf -> walkpattern)
			{
				curnode = data -> node[i];
				prevnodeID = i;
				curnodeID = curnode -> nbr[ Util::MultiSample ( curnode -> nbr_wt ) ];
				walks[j * data -> num_node + i][1] = curnodeID;
				for (int k = 2; k < conf -> walklength; k++)
				{
					curnode = data -> node[curnodeID];
					nextnodeID = curnode -> nbr[ Util::MultiSample ( data -> edge_probs[EdgeType(prevnodeID, curnodeID)] ) ];
					walks[j * data -> num_node + i][k] = nextnodeID;
					prevnodeID = curnodeID;
					curnodeID = nextnodeID;
				}
			}
			else
				for (int k = 1; k < conf -> walklength; k++)
				{
					curnode = data -> node[curnodeID];
					curnodeID = curnode -> nbr[ Util::MultiSample ( curnode -> nbr_wt ) ];
					walks[j * data -> num_node + i][k] = curnodeID;
				}
		}
	}
	return args;
}
