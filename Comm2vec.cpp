#pragma once
#include "Comm2vec.h"
#include <math.h>

void Comm2Vec::GenJq()
{
	double* node_freq = (double*)calloc(data->num_node, sizeof(double));
	for (int i = 0; i < data -> num_node * conf -> numwalks; i++)
		for (int j = 0; j < conf -> walklength; j++)
		{
			node_freq[ walks[i][j] ] += 1.0;
		}
	Util::Norm<double>(node_freq, data->num_node, false);
	//pair<int*, double*> temp = Util::alias_setup(node_freq, data -> num_node);
	//J = temp.first;
	//q = temp.second;
	node_cdf = (double*)calloc(data->num_node, sizeof(double));
	for (int i = 0; i < data -> num_node; i ++)
		node_freq[i] = pow(node_freq[i], 0.75);
		//node_freq[i] = data -> node[i] -> degree;
	Util::Norm<double>(node_freq, data->num_node, false);
	node_cdf[0] = node_freq[0];
	for (int i =1; i < data ->num_node; i++)
	{
		node_cdf[i] = node_cdf[i-1] + node_freq[i];
	}
}

int Comm2Vec::NodeSampling()
{
	srand( time(0) );
	double P = rand() / (RAND_MAX + 0.0);
	for (int i =0; i < data ->num_node; i++)
		if (P <= node_cdf[i])
			return i;
}

void Comm2Vec::Do()
{
	data = new DataSet();
	data -> LoadData( conf );
    
	walks = Walks::SimuWalks(data, conf);
	
	printf("Preparing to train ...\n");
	comm_vec.vector = Util::init_Array2D<double>( data -> num_comm, conf -> commdim, true);
	node_vec.vector = Util::init_Array2D<double>( data -> num_node, conf -> nodedim, true);
	node_comm_vec.vector = Util::init_Array2D<double>( data -> num_node, conf -> commdim, 0.0 );
	syn1neg_node.vector = Util::init_Array2D<double>( data -> num_node, conf -> nodedim + conf -> commdim, 0.0 );
	//UpdateNodeCommVec();
	MakeExptable();
	GenJq();
	Train();	
	//NormalizeL2();
	Output();
	/*FILE *test = fopen("test", "w");
	for (int i = 0; i < data -> num_node * conf -> numwalks; i ++)
	{
		for (int j = 0; j < conf -> walklength; j++)
		{
			for (int k = 0; k < data -> node[walks[i][j]] ->comm_list.size(); k ++)
			{
				fprintf(test,"%d\t", data -> node[walks[i][j]] ->comm_list[k] );
			}
			fprintf(test,"\n");
		}
		fprintf(test,"\n");
	}
	fclose(test);*/
}

void Comm2Vec::NormalizeL2()
{
	double sum;
	for (int i = 0; i < data -> num_node; i ++)
	{
		sum = 0.0;
		for (int k = 0; k < conf -> nodedim; k ++)
		{
			sum += pow( node_vec.vector[i][k], 2);
		}
		sum = pow(sum, 0.5);
		for (int k = 0; k < conf -> nodedim; k ++)
		{
			node_vec.vector[i][k] /= sum;
		}
	}
	for (int i = 0; i < data -> num_comm; i ++)
	{
		sum = 0.0;
		for (int k = 0; k < conf -> commdim; k ++)
		{
			sum += pow( comm_vec.vector[i][k], 2);
		}
		sum = pow(sum, 0.5);
		for (int k = 0; k < conf -> commdim; k ++)
		{
			comm_vec.vector[i][k] /= sum;
		}
	}
}

void Comm2Vec::Train()
{
	time_t start_time, end_time;
	printf("Training ...\n");
	start_time = time(NULL);
	node_count_actual = 0;
	int train_method = 3;
	struct argsformat
	{
		int id;
		Comm2Vec* self;
	};
	argsformat* args = (argsformat*)calloc(conf -> trainthreads, sizeof(argsformat));
	pthread_t *pt = (pthread_t *)malloc( conf -> trainthreads * sizeof(pthread_t));
	if (args == NULL || pt == NULL)
	{
		fprintf(stderr, "Memory allocation failed.");
		exit(0);
	}

	//pthread_mutex_init(&comm_vec.mutex, NULL);
	//pthread_mutex_init(&node_vec.mutex, NULL);
	//pthread_mutex_init(&node_comm_vec.mutex, NULL);
	//pthread_mutex_init(&syn1neg_node.mutex, NULL);

	for (int i = 0; i < conf->trainthreads; i++)
	{
		args[i].id = i;
		args[i].self = this;
		if (train_method == 1)
			pthread_create(&pt[i], NULL, Comm2Vec::TrainThread, (void *)&args[i]);
		else if (train_method == 2)
			pthread_create(&pt[i], NULL, Comm2Vec::TrainThreadSecond, (void *)&args[i]);
		else if (train_method == 3)
			pthread_create(&pt[i], NULL, Comm2Vec::TrainThreadThird, (void *)&args[i]);
		else if (train_method == 4)
			pthread_create(&pt[i], NULL, Comm2Vec::TrainThreadForth, (void *)&args[i]);
		else
			pthread_create(&pt[i], NULL, Comm2Vec::TrainThreadForth, (void *)&args[i]);
	}
	for (int i = 0; i < conf->trainthreads; i++) pthread_join(pt[i], NULL);
	free(args);
	//pthread_mutex_destroy(&comm_vec.mutex);
	//pthread_mutex_destroy(&node_vec.mutex);
	//pthread_mutex_destroy(&node_comm_vec.mutex);
	//pthread_mutex_destroy(&syn1neg_node.mutex);
	end_time = time(NULL);
	printf("\nTrain finished.  Time costed is %d seconds\n", (end_time - start_time));
}

void* Comm2Vec::TrainThread(void* args)
{
	srand( time(0) );
	NodeIDType target;
	int label;
	vector<NodeIDType> node_context;
	CommIDType comm_id;
	NodeIDType node_id;

	struct argsformat
	{
		int id;
		Comm2Vec* self;
	};
	int id = ((argsformat*)args) -> id;
	Comm2Vec* self = ((argsformat*)args) -> self;
	double alpha = self -> conf -> alpha;
	double start_alpha = alpha;
	
	double* neu1_node = (double*)calloc(self -> conf -> nodedim, sizeof(double));
	double* neu1e_node = (double*)calloc(self -> conf -> nodedim, sizeof(double));
	double* neu1e_comm = (double*)calloc(self -> conf -> commdim, sizeof(double));

	if (neu1_node == NULL || neu1e_node == NULL || neu1e_comm == NULL)
	{
		fprintf(stderr, "Memory allocation failed.");
		exit(0);
	}
	
	double f, g;

	int start = (double)id / self -> conf -> trainthreads * self -> data -> num_node * self -> conf -> numwalks;
	int end = ((double)id + 1.0) / self -> conf-> trainthreads * self -> data->num_node * self -> conf->numwalks;
	for (int milestone = start; milestone < end; milestone++ )
	{
		if (milestone - start > 100)
		{
			self -> node_count_actual += milestone - start;
			start = milestone;
			alpha = start_alpha * (1 - double(self -> node_count_actual)/(self -> data -> num_node * self -> conf -> numwalks + 1));
			if (alpha < start_alpha * 0.01) alpha = start_alpha * 0.01;
			printf("\r	Alpha: %f    Progress: %.2f%%", alpha, (self -> node_count_actual + 1) / double(self -> data -> num_node * self -> conf -> numwalks) * 100);
			fflush(stdout);
		}
		
		for(int walk_position = 0; walk_position < self -> conf -> walklength; walk_position++)
		{

			memset(neu1_node, 0, sizeof(double) * self -> conf -> nodedim);
			memset(neu1e_node, 0, sizeof(double) * self -> conf -> nodedim);
			node_context.clear();
			int b = rand() % self -> conf -> window;
			for (int a = b; a < self -> conf->window *2 + 1 - b; a++) if (a != self -> conf -> window)
			{
				int c = walk_position - self -> conf->window + a;
				if ( c < 0 || c >= self -> conf -> walklength ) continue;
				for (int k = 0; k < self -> conf -> nodedim; k++)
				{
					neu1_node[k] += self -> node_vec.vector[ self -> walks[milestone][c] ][k];
					neu1_node[k] += self -> node_comm_vec.vector[ self -> walks[milestone][c] ][k];
				}
				node_context.push_back(self -> walks[milestone][c]);
			}

			for (int d = 0; d < self -> conf -> negative + 1; d ++)
			{
				if (d == 0)
				{
					target = self -> walks[milestone][walk_position];
					label = 1;
				}
				else
				{
					do
					{
						//target = rand() / (RAND_MAX + 0.0) * (self -> data -> num_node - 1);
						target = self -> NodeSampling();
					}
					while( Util::HaveIntersection<CommIDType>( self -> data -> node[target] ->comm_list, self -> data -> node[ self -> walks[milestone][walk_position] ] ->comm_list) );
					label = 0;
				}
				memset(&f, 0, sizeof(double));
				//pthread_mutex_lock(&self -> syn1neg_node.mutex);
				for (int k = 0; k < self -> conf -> nodedim; k++)
				{
					f += neu1_node[k] * self -> syn1neg_node.vector[target][k];
				}
				if (f > MAX_EXP) 
					g = (label - 1) * alpha;
				else if (f < -MAX_EXP) 
					g = (label - 0) * alpha;
				else 
					g = (label - self -> exptable[(int)((f + MAX_EXP) * (EXP_TABLE_SIZE / MAX_EXP / 2))]) * alpha;
				
				for ( int k = 0; k < self -> conf -> nodedim; k++)
				{
					neu1e_node[k] += g * self -> syn1neg_node.vector[target][k];
					self -> syn1neg_node.vector[target][k] += g * neu1_node[k];
				}
				//pthread_mutex_unlock(&self -> syn1neg_node.mutex);
			}

			for (vector<NodeIDType>::iterator it = node_context.begin(); it != node_context.end(); ++it)
			{
				
				//pthread_mutex_lock(&self -> node_vec.mutex);
				for (int k = 0; k < self -> conf -> nodedim; k++)
				{
					self -> node_vec.vector[ *it ][k] += neu1e_node[k];
				}
				//pthread_mutex_unlock(&self -> node_vec.mutex);
				//update comm vector
				//pthread_mutex_lock(&self -> comm_vec.mutex);
				for (int i = 0; i < self -> data -> node[*it] -> comm_list.size(); ++i)
				{
					comm_id = self -> data -> node[*it] -> comm_list[i];
					for (int k = 0; k < self -> conf -> commdim; k++)
					{
						neu1e_comm[k] = neu1e_node[k] * self -> data -> node[*it] -> weight[i];
						self -> comm_vec.vector[ comm_id ][k] += neu1e_comm[k];
						
					}
					for(int j = 0; j < self -> data -> comm[comm_id] -> node_list.size(); j++)
					{
						node_id = self -> data -> comm[comm_id] -> node_list[j];
						for (int k = 0; k < self -> conf -> commdim; k++)
						{
							self -> node_comm_vec.vector[node_id][k] += neu1e_comm[k] * self -> data -> comm[comm_id] -> weight[j];
						}
					}
					
				}
				//pthread_mutex_unlock(&self -> comm_vec.mutex);
			}
			//pthread_mutex_lock(&self -> comm_vec.mutex);
			//pthread_mutex_lock(&self -> node_comm_vec.mutex);
			//Comm2Vec::UpdateNodeCommVec(self, node_need_update_comm);
			//pthread_mutex_unlock(&self -> node_comm_vec.mutex);
			//pthread_mutex_unlock(&self -> comm_vec.mutex);
			
		}
	}
}

void Comm2Vec::UpdateNodeCommVec()
{
	Util::format_Array2D<double>(node_comm_vec.vector, data -> num_node, conf -> commdim);
	for (int i = 0; i < data -> num_node; i++)
	{
		for (int k = 0; k < conf -> commdim; k++)
		{
			for (int it = 0; it < data -> node[i] -> comm_list.size(); ++it)
			{
				node_comm_vec.vector[i][k] += comm_vec.vector[ data -> node[i] -> comm_list[it] ][k] * data -> node[i] -> weight[it];
			}
		}
	}
}

void Comm2Vec::UpdateNodeCommVec(Comm2Vec* self, vector<CommIDType> &node_need_update_comm)
{
	for (vector<CommIDType>::iterator i = node_need_update_comm.begin(); i != node_need_update_comm.end(); ++i)
	{
		for (int it = 0; it < self -> data -> node[*i] -> comm_list.size(); it++)
		{
			memset(self -> node_comm_vec.vector[*i], 0, sizeof(double) * self -> conf -> commdim);
			for (int k = 0; k < self -> conf -> commdim; k++)
			{
				self -> node_comm_vec.vector[*i][k] += self -> comm_vec.vector[ self -> data -> node[*i] -> comm_list[it] ][k] * self -> data -> node[*i] -> weight[it];
			}
		}
	}
}

void Comm2Vec::MakeExptable()
{
	exptable = (double*)calloc(EXP_TABLE_SIZE + 1, sizeof(double));
	if (exptable == NULL)
	{
		fprintf(stderr, "Memory allocation failed.");
		exit(0);
	}
	for (int i = 0; i < EXP_TABLE_SIZE + 1; i++)
	{
		exptable[i] = exp(((double)i / EXP_TABLE_SIZE * 2 - 1) * MAX_EXP);
		exptable[i] = exptable[i] / (exptable[i] + 1);
	}
}

void Comm2Vec::Output()
{
	printf("Outputing the model into the disk ...\n");
	FILE *fout_nodeemb = fopen(conf -> nodeembfile.c_str(), "w");
	FILE *fout_commemb = fopen(conf -> commembfile.c_str(), "w");
	
	//fprintf(fout_nodeemb, "#Nodes      %d\n", data -> num_node);
	for (int i = 0; i < data -> num_node; i++)
	{
		for (int k = 0; k < conf -> nodedim; k++)
		{
			fprintf(fout_nodeemb, "%f\t", node_vec.vector[i][k]);
		}
		fprintf(fout_nodeemb, "\n");
	}
	fclose(fout_nodeemb);

	//fprintf(fout_commemb, "#Communities      %d\n", data -> num_comm);
	for (int i = 0; i < data -> num_comm; i++)
	{
		for (int k = 0; k < conf -> commdim; k++)
		{
			fprintf(fout_commemb, "%f\t", comm_vec.vector[i][k]);
		}
		fprintf(fout_commemb, "\n");
	}

	fclose(fout_commemb);
	printf("Output finished.\n");
}

void* Comm2Vec::TrainThreadSecond(void* args)
{
	srand( time(0) );
	NodeIDType target;
	int label;
	vector<NodeIDType> node_context;

	CommIDType comm_id;
	NodeIDType node_id;

	struct argsformat
	{
		int id;
		Comm2Vec* self;
	};
	int id = ((argsformat*)args) -> id;
	Comm2Vec* self = ((argsformat*)args) -> self;
	double alpha = self -> conf -> alpha;
	double start_alpha = alpha;


	double* neu1_node = (double*)calloc(self -> conf -> nodedim, sizeof(double));
	double* neu1e_node = (double*)calloc(self -> conf -> nodedim, sizeof(double));
	double* neu1e_comm = (double*)calloc(self -> conf -> commdim, sizeof(double));

	if (neu1_node == NULL || neu1e_node == NULL || neu1e_comm == NULL)
	{
		fprintf(stderr, "Memory allocation failed.");
		exit(0);
	}
	
	double f, g;

	int start = (double)id / self -> conf -> trainthreads * self -> data -> num_node * self -> conf -> numwalks;
	int end = ((double)id + 1.0) / self -> conf-> trainthreads * self -> data->num_node * self -> conf->numwalks;
	for (int milestone = start; milestone < end; milestone++ )
	{
		if (milestone - start > 100)
		{
			self -> node_count_actual += milestone - start;
			start = milestone;
			alpha = start_alpha * (1 - double(self -> node_count_actual)/(self -> data -> num_node * self -> conf -> numwalks + 1));
			if (alpha < start_alpha * 0.01) alpha = start_alpha * 0.01;
			printf("\r	Alpha: %f    Progress: %.2f%%", alpha, (self -> node_count_actual + 1) / double(self -> data -> num_node * self -> conf -> numwalks) * 100);
		}
		
		for(int walk_position = 0; walk_position < self -> conf -> walklength; walk_position++)
		{
			memset(neu1_node, 0, sizeof(double) * self -> conf -> nodedim);
			memset(neu1e_node, 0, sizeof(double) * self -> conf -> nodedim);
			node_context.clear();
			int b = rand() % self -> conf -> window;
			for (int a = b; a < self -> conf->window *2 + 1 - b; a++) if (a != self -> conf -> window)
			{
				int c = walk_position - self -> conf->window + a;
				if ( c < 0 || c >= self -> conf -> walklength ) continue;
				node_id = self -> walks[milestone][c];
				for (int k = 0; k < self -> conf -> nodedim; k++)
				{
					neu1_node[k] += self -> node_vec.vector[ node_id ][k];
					for (int i = 0; i < self -> data -> node[ node_id ] -> comm_list.size(); i++)
						neu1_node[k] += self -> comm_vec.vector[ self -> data -> node[ node_id ] -> comm_list[i] ][k] * self -> data -> node[ node_id ] -> weight[i];
				}
				node_context.push_back(node_id);
			}

			for (int d = 0; d < self -> conf -> negative + 1; d ++)
			{
				if (d == 0)
				{
					target = self -> walks[milestone][walk_position];
					label = 1;
				}
				else
				{
					do
					{
						//target = rand() / (RAND_MAX + 0.0) * (self -> data -> num_node - 1);
						target = self -> NodeSampling();
					}
					while( Util::HaveIntersection<CommIDType>( self -> data -> node[target] ->comm_list, self -> data -> node[ self -> walks[milestone][walk_position] ] ->comm_list) );
						//while(target == self -> walks[milestone][walk_position]);
					label = 0;
				}
				memset(&f, 0, sizeof(double));
				for (int k = 0; k < self -> conf -> nodedim; k++)
				{
					f += neu1_node[k] * self -> syn1neg_node.vector[target][k];
				}
				if (f > MAX_EXP) 
					g = (label - 1) * alpha;
				else if (f < -MAX_EXP) 
					g = (label - 0) * alpha;
				else 
					g = (label - self -> exptable[(int)((f + MAX_EXP) * (EXP_TABLE_SIZE / MAX_EXP / 2))]) * alpha;
				
				for ( int k = 0; k < self -> conf -> nodedim; k++)
				{
					neu1e_node[k] += g * self -> syn1neg_node.vector[target][k];
					self -> syn1neg_node.vector[target][k] += g * neu1_node[k];
				}
			}

			for (vector<NodeIDType>::iterator it = node_context.begin(); it != node_context.end(); ++it)
			{
				
				for (int k = 0; k < self -> conf -> nodedim; k++)
				{
					self -> node_vec.vector[ *it ][k] += neu1e_node[k];
				}
				//update comm vector
				for (int i = 0; i < self -> data -> node[*it] -> comm_list.size(); ++i)
				{
					comm_id = self -> data -> node[*it] -> comm_list[i];
					for (int k = 0; k < self -> conf -> commdim; k++)
					{
						neu1e_comm[k] = neu1e_node[k] * self -> data -> node[*it] -> weight[i];
						self -> comm_vec.vector[ comm_id ][k] += neu1e_comm[k];
						
					}
					
				}
			}
			
		}
	}
}

void* Comm2Vec::TrainThreadThird(void* args)
{
	srand( time(0) );
	NodeIDType target;
	int label;
	vector<NodeIDType> node_context;

	CommIDType comm_id;
	NodeIDType node_id;

	struct argsformat
	{
		int id;
		Comm2Vec* self;
	};
	int id = ((argsformat*)args) -> id;
	Comm2Vec* self = ((argsformat*)args) -> self;
	double alpha = self -> conf -> alpha;
	double start_alpha = alpha;


	double* neu1_node = (double*)calloc(self -> conf -> nodedim, sizeof(double));
	double* neu1e_node = (double*)calloc(self -> conf -> nodedim, sizeof(double));
	double* neu1e_comm = (double*)calloc(self -> conf -> commdim, sizeof(double));

	if (neu1_node == NULL || neu1e_node == NULL || neu1e_comm == NULL)
	{
		fprintf(stderr, "Memory allocation failed.");
		exit(0);
	}
	
	double f, g;

	int start = (double)id / self -> conf -> trainthreads * self -> data -> num_node * self -> conf -> numwalks;
	int end = ((double)id + 1.0) / self -> conf-> trainthreads * self -> data->num_node * self -> conf->numwalks;
	for (int milestone = start; milestone < end; milestone++ )
	{
		if (milestone - start > 100)
		{
			self -> node_count_actual += milestone - start;
			start = milestone;
			alpha = start_alpha * (1 - double(self -> node_count_actual)/(self -> data -> num_node * self -> conf -> numwalks + 1));
			if (alpha < start_alpha * 0.01) alpha = start_alpha * 0.01;
			printf("\r	Alpha: %f    Progress: %.2f%%", alpha, (self -> node_count_actual + 1) / double(self -> data -> num_node * self -> conf -> numwalks) * 100);
		}
		
		for(int walk_position = 0; walk_position < self -> conf -> walklength; walk_position++)
		{
			memset(neu1_node, 0, sizeof(double) * self -> conf -> nodedim);
			memset(neu1e_node, 0, sizeof(double) * self -> conf -> nodedim);
			node_context.clear();
			int b = rand() % self -> conf -> window;
			for (int a = b - b; a < self -> conf->window *2 + 1 - b + b; a++) if (a != self -> conf -> window)
			{
				int c = walk_position - self -> conf->window + a;
				if ( c < 0 || c >= self -> conf -> walklength ) continue;
				node_id = self -> walks[milestone][c];
				for (int k = 0; k < self -> conf -> nodedim; k++)
				{
					neu1_node[k] += self -> node_vec.vector[ node_id ][k];
					//for (int i = 0; i < self -> data -> node[ node_id ] -> comm_list.size(); i++)
						//neu1_node[k] += self -> comm_vec.vector[ self -> data -> node[ node_id ] -> comm_list[i] ][k] * self -> data -> node[ node_id ] -> weight[i];
				}
				node_context.push_back(node_id);
			}
			//average the context nodes
			/*for (int k = 0; k < self -> conf -> nodedim; k++ )
				neu1_node[k] /= node_context.size();*/
			//add cmty context
			for (int i = 0; i < self -> data -> node[ self -> walks[milestone][walk_position] ] ->comm_list.size(); i++)
				for (int k = 0; k < self -> conf -> nodedim; k++ )
					//neu1_node[k] += self -> comm_vec.vector[ self -> data -> node[walk_position] ->comm_list[i] ][k] * self -> data -> node[walk_position] ->weight[i];
					neu1_node[k] += self -> comm_vec.vector[ self -> data -> node[ self -> walks[milestone][walk_position] ] ->comm_list[i] ][k];
			//int num = node_context.size() + self -> data -> node[walk_position] ->comm_list.size();
			/*for (int k = 0; k < self -> conf -> nodedim; k++)
			{
				neu1_node[k] /= num;
			}*/
			for (int d = 0; d < self -> conf -> negative + 1; d ++)
			{
				if (d == 0)
				{
					target = self -> walks[milestone][walk_position];
					label = 1;
				}
				else
				{
					do
					{
						//target = rand() / (RAND_MAX + 0.0) * (self -> data -> num_node - 1);
						//target = NodeIDType(target + 0.5);
						target = self -> NodeSampling();
					}
					while( Util::HaveIntersection<CommIDType>( self -> data -> node[target] ->comm_list, self -> data -> node[ self -> walks[milestone][walk_position] ] ->comm_list) );
					//while(target == self -> walks[milestone][walk_position]);
					label = 0;
				}
				//memset(&f, 0, sizeof(double));
				f = 0.0;
				for (int k = 0; k < self -> conf -> nodedim; k++)
				{
					f += neu1_node[k] * self -> syn1neg_node.vector[target][k];
				}
				if (f > MAX_EXP) 
					g = (label - 1) * alpha;
				else if (f < -MAX_EXP) 
					g = (label - 0) * alpha;
				else 
					g = (label - self -> exptable[(int)((f + MAX_EXP) * (EXP_TABLE_SIZE / MAX_EXP / 2))]) * alpha;
				
				for ( int k = 0; k < self -> conf -> nodedim; k++)
				{
					neu1e_node[k] += g * self -> syn1neg_node.vector[target][k];
					self -> syn1neg_node.vector[target][k] += g * neu1_node[k];										
				}
			}

			//update comm vector
			for (int i = 0; i < self -> data -> node[self -> walks[milestone][walk_position]] -> comm_list.size(); i++)
			{
				comm_id = self -> data -> node[self -> walks[milestone][walk_position]] -> comm_list[i];
				for (int k = 0; k < self -> conf -> commdim; k++)
				{
					//neu1e_comm[k] = neu1e_node[k] * self -> data -> node[walk_position] -> weight[i];
					self -> comm_vec.vector[ comm_id ][k] += neu1e_node[k];
					
				}
				
			}
			/*for (int k = 0; k < self -> conf -> nodedim; k++ )
				neu1e_node[k] /= node_context.size();*/
			//update node vector
			for (vector<NodeIDType>::iterator it = node_context.begin(); it != node_context.end(); ++it)
			{
				
				for (int k = 0; k < self -> conf -> nodedim; k++)
				{
					self -> node_vec.vector[ *it ][k] += neu1e_node[k];
				}
				
			}
			
		}
	}
}



//concatenate the node and the community
void* Comm2Vec::TrainThreadForth(void* args)
{
	srand( time(0) );
	NodeIDType target;
	int label;
	vector<NodeIDType> node_context;

	CommIDType comm_id;
	NodeIDType node_id;

	struct argsformat
	{
		int id;
		Comm2Vec* self;
	};
	int id = ((argsformat*)args) -> id;
	Comm2Vec* self = ((argsformat*)args) -> self;
	double alpha = self -> conf -> alpha;
	double start_alpha = alpha;
	int nodedim = self -> conf -> nodedim;
	int commdim = self -> conf -> commdim;


	double* neu1_node = (double*)calloc(nodedim + commdim, sizeof(double));
	double* neu1e_node = (double*)calloc(nodedim + commdim, sizeof(double));
	double* neu1e_comm = (double*)calloc(commdim, sizeof(double));

	if (neu1_node == NULL || neu1e_node == NULL || neu1e_comm == NULL)
	{
		fprintf(stderr, "Memory allocation failed.");
		exit(0);
	}
	
	double f, g;

	int start = (double)id / self -> conf -> trainthreads * self -> data -> num_node * self -> conf -> numwalks;
	int end = ((double)id + 1.0) / self -> conf-> trainthreads * self -> data->num_node * self -> conf->numwalks;
	for (int milestone = start; milestone < end; milestone++ )
	{
		if (milestone - start > 100)
		{
			self -> node_count_actual += milestone - start;
			start = milestone;
			alpha = start_alpha * (1 - double(self -> node_count_actual)/(self -> data -> num_node * self -> conf -> numwalks + 1));
			if (alpha < start_alpha * 0.01) alpha = start_alpha * 0.01;
			printf("\r	Alpha: %f    Progress: %.2f%%", alpha, (self -> node_count_actual + 1) / double(self -> data -> num_node * self -> conf -> numwalks) * 100);
		}
		
		for(int walk_position = 0; walk_position < self -> conf -> walklength; walk_position++)
		{
			memset(neu1_node, 0, sizeof(double) * (nodedim + commdim));
			memset(neu1e_node, 0, sizeof(double) * (nodedim + commdim));
			node_context.clear();
			int b = rand() % self -> conf -> window + 1;
			for (int a = b + b; a < self -> conf->window *2 + 1 - b + b; a++) if (a != self -> conf -> window)
			{
				int c = walk_position - self -> conf->window + a;
				if ( c < 0 || c >= self -> conf -> walklength ) continue;
				node_id = self -> walks[milestone][c];
				for (int k = 0; k < nodedim; k++)
				{
					neu1_node[k] += self -> node_vec.vector[ node_id ][k];
					//for (int i = 0; i < self -> data -> node[ node_id ] -> comm_list.size(); i++)
						//neu1_node[k] += self -> comm_vec.vector[ self -> data -> node[ node_id ] -> comm_list[i] ][k] * self -> data -> node[ node_id ] -> weight[i];
				}
				node_context.push_back(node_id);
			}
			//average the context nodes
			//for (int k = 0; k < self -> conf -> nodedim; k++ )
				//neu1_node[k] /= node_context.size();
			//add cmty context
			for (int i = 0; i < self -> data -> node[self -> walks[milestone][walk_position]] ->comm_list.size(); i++)
				for (int k = 0; k < commdim; k++ )
					//neu1_node[k] += self -> comm_vec.vector[ self -> data -> node[walk_position] ->comm_list[i] ][k] * self -> data -> node[walk_position] ->weight[i];
					neu1_node[k + nodedim] += self -> comm_vec.vector[ self -> data -> node[self -> walks[milestone][walk_position]] ->comm_list[i] ][k];
			//int num = node_context.size() + self -> data -> node[walk_position] ->comm_list.size();
			/*for (int k = 0; k < self -> conf -> nodedim; k++)
			{
				neu1_node[k] /= num;
			}*/
			for (int d = 0; d < self -> conf -> negative + 1; d ++)
			{
				if (d == 0)
				{
					target = self -> walks[milestone][walk_position];
					label = 1;
				}
				else
				{
					do
					{
						//target = rand() / (RAND_MAX + 0.0) * (self -> data -> num_node - 1);
						//target = NodeIDType(target + 0.5);
						target = self -> NodeSampling();
						
					}
					while( Util::HaveIntersection<CommIDType>( self -> data -> node[target] ->comm_list, self -> data -> node[ self -> walks[milestone][walk_position] ] ->comm_list) );
					//while(target == self -> walks[milestone][walk_position]);
					label = 0;
				}
				//memset(&f, 0, sizeof(double));
				f = 0.0;
				for (int k = 0; k < (nodedim + commdim); k++)
				{
					f += neu1_node[k] * self -> syn1neg_node.vector[target][k];
				}
				if (f > MAX_EXP) 
					g = (label - 1) * alpha;
				else if (f < -MAX_EXP) 
					g = (label - 0) * alpha;
				else 
					g = (label - self -> exptable[(int)((f + MAX_EXP) * (EXP_TABLE_SIZE / MAX_EXP / 2))]) * alpha;
				
				for ( int k = 0; k < (nodedim + commdim); k++)
				{
					neu1e_node[k] += g * self -> syn1neg_node.vector[target][k];
					self -> syn1neg_node.vector[target][k] += g * neu1_node[k];										
				}
			}

			//update comm vector
			for (int i = 0; i < self -> data -> node[self -> walks[milestone][walk_position]] -> comm_list.size(); i++)
			{
				comm_id = self -> data -> node[self -> walks[milestone][walk_position]] -> comm_list[i];
				for (int k = 0; k < commdim; k++)
				{
					//neu1e_comm[k] = neu1e_node[k] * self -> data -> node[walk_position] -> weight[i];
					self -> comm_vec.vector[ comm_id ][k] += neu1e_node[k + nodedim];
					
				}
				
			}
			//for (int k = 0; k < self -> conf -> nodedim; k++ )
				//neu1e_node[k] /= node_context.size();
			//update node vector
			for (vector<NodeIDType>::iterator it = node_context.begin(); it != node_context.end(); ++it)
			{
				
				for (int k = 0; k < nodedim; k++)
				{
					self -> node_vec.vector[ *it ][k] += neu1e_node[k];
				}
				
			}
			
		}
	}
}

void* Comm2Vec::TrainThreadFifth(void* args)
{
	srand( time(0) );
	NodeIDType target;
	int label;
	vector<NodeIDType> node_context;

	CommIDType comm_id;
	NodeIDType node_id;

	struct argsformat
	{
		int id;
		Comm2Vec* self;
	};
	int id = ((argsformat*)args) -> id;
	Comm2Vec* self = ((argsformat*)args) -> self;
	double alpha = self -> conf -> alpha;
	double start_alpha = alpha;


	double* neu1_node = (double*)calloc(self -> conf -> nodedim, sizeof(double));
	double* neu1e_node = (double*)calloc(self -> conf -> nodedim, sizeof(double));
	double* neu1e_comm = (double*)calloc(self -> conf -> commdim, sizeof(double));

	if (neu1_node == NULL || neu1e_node == NULL || neu1e_comm == NULL)
	{
		fprintf(stderr, "Memory allocation failed.");
		exit(0);
	}
	
	double f, g;

	int start = (double)id / self -> conf -> trainthreads * self -> data -> num_node * self -> conf -> numwalks;
	int end = ((double)id + 1.0) / self -> conf-> trainthreads * self -> data->num_node * self -> conf->numwalks;
	for (int milestone = start; milestone < end; milestone++ )
	{
		if (milestone - start > 100)
		{
			self -> node_count_actual += milestone - start;
			start = milestone;
			alpha = start_alpha * (1 - double(self -> node_count_actual)/(self -> data -> num_node * self -> conf -> numwalks + 1));
			if (alpha < start_alpha * 0.01) alpha = start_alpha * 0.01;
			printf("\r	Alpha: %f    Progress: %.2f%%", alpha, (self -> node_count_actual + 1) / double(self -> data -> num_node * self -> conf -> numwalks) * 100);
		}
		
		for(int walk_position = 0; walk_position < self -> conf -> walklength; walk_position++)
		{
			
			node_context.clear();
			int b = rand() % self -> conf -> window;
			for (int a = b - b; a < self -> conf->window *2 + 1 - b + b; a++) if (a != self -> conf -> window)
			{
				memset(neu1e_node, 0, sizeof(double) * self -> conf -> nodedim);
				int c = walk_position - self -> conf->window + a;
				if ( c < 0 || c >= self -> conf -> walklength ) continue;
				node_id = self -> walks[milestone][c];
				for (int k = 0; k < self -> conf -> nodedim; k++)
				{
					neu1_node[k] += self -> node_vec.vector[ node_id ][k];
					
				}
				node_context.push_back(node_id);
				
				for (int d = 0; d < self -> conf -> negative + 1; d ++)
				{
					if (d == 0)
					{
						target = self -> walks[milestone][walk_position];
						label = 1;
					}
					else
					{
						do
						{
							//target = rand() / (RAND_MAX + 0.0) * (self -> data -> num_node - 1);
							//target = NodeIDType(target + 0.5);
							target = self -> NodeSampling();
						}
						while( Util::HaveIntersection<CommIDType>( self -> data -> node[target] ->comm_list, self -> data -> node[ self -> walks[milestone][walk_position] ] ->comm_list) );
						//while(target == self -> walks[milestone][walk_position]);
						label = 0;
					}
					//memset(&f, 0, sizeof(double));
					f = 0.0;
					for (int k = 0; k < self -> conf -> nodedim; k++)
					{
						f += neu1_node[k] * self -> syn1neg_node.vector[target][k];
					}
					if (f > MAX_EXP) 
						g = (label - 1) * alpha;
					else if (f < -MAX_EXP) 
						g = (label - 0) * alpha;
					else 
						g = (label - self -> exptable[(int)((f + MAX_EXP) * (EXP_TABLE_SIZE / MAX_EXP / 2))]) * alpha;
				
					for ( int k = 0; k < self -> conf -> nodedim; k++)
					{
						neu1e_node[k] += g * self -> syn1neg_node.vector[target][k];
						self -> syn1neg_node.vector[target][k] += g * self -> node_vec.vector[ node_id ][k];										
					}
				}
				
				for (int k = 0; k < self -> conf -> nodedim; k++)
				{
					self -> node_vec.vector[ node_id ][k] += neu1e_node[k];
				}
				
			}
			
			
			//update comm vector
			/*for (int i = 0; i < self -> data -> node[walk_position] -> comm_list.size(); i++)
			{
				comm_id = self -> data -> node[self -> walks[milestone][walk_position]] -> comm_list[i];
				for (int d = 0; d < self -> conf -> negative + 1; d ++)
				{
					if (d == 0)
					{
						target = self -> walks[milestone][walk_position];
						label = 1;
					}
					else
					{
						do
						{
							//target = rand() / (RAND_MAX + 0.0) * (self -> data -> num_node - 1);
							//target = NodeIDType(target + 0.5);
							target = self -> NodeSampling();
						}
						while( Util::HaveIntersection<CommIDType>( self -> data -> node[target] ->comm_list, self -> data -> node[ self -> walks[milestone][walk_position] ] ->comm_list) );
						//while(target == self -> walks[milestone][walk_position]);
						label = 0;
					}
					//memset(&f, 0, sizeof(double));
					f = 0.0;
					for (int k = 0; k < self -> conf -> nodedim; k++)
					{
						f += self -> comm_vec.vector[ comm_id ][k] * self -> syn1neg_node.vector[target][k];
					}
					if (f > MAX_EXP) 
						g = (label - 1) * alpha;
					else if (f < -MAX_EXP) 
						g = (label - 0) * alpha;
					else 
						g = (label - self -> exptable[(int)((f + MAX_EXP) * (EXP_TABLE_SIZE / MAX_EXP / 2))]) * alpha;
				
					for ( int k = 0; k < self -> conf -> nodedim; k++)
					{
						neu1e_node[k] += g * self -> syn1neg_node.vector[target][k];
						self -> syn1neg_node.vector[target][k] += g * self -> comm_vec.vector[ comm_id ][k];										
					}
				}

				
				for (int k = 0; k < self -> conf -> commdim; k++)
				{
					//neu1e_comm[k] = neu1e_node[k] * self -> data -> node[walk_position] -> weight[i];
					self -> comm_vec.vector[ comm_id ][k] += neu1e_node[k];
					
				}
				
			}*/
					
			
		}
	}
}
