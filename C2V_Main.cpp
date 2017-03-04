#pragma once
#include "Config.h"
#include "DataSet.h"
#include "Comm2vec.h"


#include <stdio.h>


//int Walks::node_count_actual = 0;

int main(int argc, char* argv[])
{
    
    // Load Configuartion
    Config* conf = new Config();
    if (! conf->LoadConfig(argc, argv))
    {
        conf->ShowUsage();
        exit( 0 );
    }
    
    Comm2Vec* model = new Comm2Vec(conf);
    model -> Do();
    delete(model);
	       
}
