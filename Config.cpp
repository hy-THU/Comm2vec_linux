#include "Config.h"

void Config::SetDefault()
{

	topofile = "topo.txt";
	commfile = "cmty.txt";
	nodeembfile = "node.embedding";
	commembfile = "comm.embedding";

	nodedim = 20;
	commdim = 20;
	numwalks = 1;
	walklength = 100;
	window = 10;
	sample  = 100;
	hs  = 3;
	negative  = 5;
	trainthreads  = 2;
	walkthreads = 1;
	alpha = 0.025;
	binary = 0;
	directed = 0;
	weighted = 0;
	walkpattern = 1;
	p = 1.0;
	q = 1.0;

}

bool Config::LoadConfig(int argc, char* argv[])
{

	if (argc == 1) return 0;

	int i = 1;

	while (i < argc)
	{

        if (strcmp(argv[i], "-topofile") == 0)
        {
            this->topofile = argv[++i]; ++ i;
        }
        else if (strcmp(argv[i], "-commfile") == 0)
        {
            this->commfile = argv[++i]; ++i;
        }
        else if (strcmp(argv[i], "-nodeembfile") == 0)
        {
            this->nodeembfile = argv[++i]; ++i;
        }
		else if (strcmp(argv[i], "-commembfile") == 0)
        {
            this->commembfile = argv[++i]; ++i;
        }
		else if (strcmp(argv[i], "-nodedim") == 0)
        {
            this->nodedim = atoi(argv[++i]); ++i;
        }
		else if (strcmp(argv[i], "-commdim") == 0)
        {
            this->commdim = atoi(argv[++i]); ++i;
        }
		else if (strcmp(argv[i], "-numwalks") == 0)
        {
			this->numwalks = atoi(argv[++i]); ++i;
        }
		else if (strcmp(argv[i], "-walklength") == 0)
        {
			this->walklength = atoi(argv[++i]); ++i;
        }
		else if (strcmp(argv[i], "-window") == 0)
        {
			this->window = atoi(argv[++i]); ++i;
        }
		else if (strcmp(argv[i], "-sample") == 0)
        {
            this->sample = atof(argv[++i]); ++i;
        }
		else if (strcmp(argv[i], "-hs") == 0)
        {
            this->hs = atoi(argv[++i]); ++i;
        }
		else if (strcmp(argv[i], "-negative") == 0)
        {
            this->negative = atoi(argv[++i]); ++i;
        }
		else if (strcmp(argv[i], "-trainthreads") == 0)
        {
            this->trainthreads = atoi(argv[++i]); ++i;
        }
		else if (strcmp(argv[i], "-walkthreads") == 0)
        {
            this->walkthreads = atoi(argv[++i]); ++i;
        }
		else if (strcmp(argv[i], "-alpha") == 0)
        {
            this->alpha = atof(argv[++i]); ++i;
        }
		else if (strcmp(argv[i], "-binary") == 0)
        {
            this->binary = atoi(argv[++i]); ++i;
        }
		else if (strcmp(argv[i], "-directed") == 0)
        {
            this->directed = atoi(argv[++i]); ++i;
        }
		else if (strcmp(argv[i], "-weighted") == 0)
        {
			this->weighted = atoi(argv[++i]); ++i;
        }
		else if (strcmp(argv[i], "-walkpattern") == 0)
        {
			this->walkpattern = atoi(argv[++i]); ++i;
        }
        else if (strcmp(argv[i], "-p") == 0)
        {
			this->p = atof(argv[++i]); ++i;
        }
		else if (strcmp(argv[i], "-q") == 0)
        {
			this->q = atof(argv[++i]); ++i;
        }
		else ++ i;
	}
    return 1;
}

void Config::ShowUsage()
{
    printf("Comm2vec v0.1                                             \n");
    printf("     by Yu Han, Tsinghua University                       \n");
    printf("                                                          \n");
    printf("Usage: ./C2V [options]                                    \n");
    printf(" Options:                                                 \n");
    printf("   -topofile     string : topology file (default: 'topo')                  \n");
    printf("   -commfile     string : community file (default: 'comm')                \n");
	printf("   -nodeembfile  string : node embedding file (default: 'node.embedding')          \n");
	printf("   -commembfile  string : community embedding file (default: 'comm.embedding')     \n");
	printf("   -nodedim      int    : number of dimensions for nodes (default: 100) \n");
	printf("   -commdim      int    : number of dimensions for communities (default: 100)      \n");
	printf("   -numwalks     int    : number of walks (default: 1)    \n");
	printf("   -walklength   int    : walk length (default: 100)    \n");
	printf("   -window       int    : max skip length between nodes (default: 5)    \n");
	printf("   -sample       double : threshold for occurrence of words. (default: 0)             \n");
	printf("   -hs           int    : whether to use Hierarchical Softmax (default: 1)            \n");
    printf("   -negative     int    : number of negative examples (default: 0)                     \n");
	printf("   -trainthreads int    : number of threads to use for trainning (default: 1)                      \n");
	printf("   -walkthreads  int    : number of threads to use for walking (default: 1)                      \n");
	printf("   -alpha        double : starting learning rate (default: 0.01)                  \n");
	printf("   -binary       int    : whether to save the resulting vectors in binary mode (default: 0) \n");
	printf("   -directed     int    : whether the network is directed (default: 0) \n");
	printf("   -weighted     int    : whether the edges are weighted (default: 0) \n");
	printf("   -walkpattern  int    : walk in the way of Deepwalk(0) or Node2vec(1) (default: 1) \n");
	printf("   -p            double : p value for Node2vec walkpattern (default: 1.0) \n");
	printf("   -q            double : q value for Node2vec walkpattern (default: 1.0) \n");

    	//exit( 0 );
}
