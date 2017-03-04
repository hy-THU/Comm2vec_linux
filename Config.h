#pragma once

#include <string>
#include <stdio.h>
#include <iostream>
#include <cstdlib>
#include <cstring>

using std::string;



class Config
{
public:

	string topofile;
    string commfile;
	string nodeembfile;
	string commembfile;

	int nodedim;
	int commdim;
	int numwalks;
	int walklength;
	int window;
	double sample;
	int hs;
	int negative;
	int walkthreads;
	int trainthreads;
	double alpha;
	int binary;
	int directed;
	int weighted;
	int walkpattern;
	double p;
	double q;




    Config(){ SetDefault(); }
    void SetDefault();

    // false => parameter wrong
    bool LoadConfig(int argc, char* argv[]);
    static void ShowUsage();
};
