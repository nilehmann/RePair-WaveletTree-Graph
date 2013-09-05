/**
 * Test running times
 * */


#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <string>
#include <algorithm>

#include <BitSequenceBuilder.h>
#include <Mapper.h>
#include "GraphSequenceTest.h"


#define TOMICRO(C) C*1000000

using namespace std;




int main (int argc, char **argv)
{ 	
	srand(time(NULL));
	if( argc != 7) {
	    cout << "usage: " << argv[0] << " <graphname> <rrr_sample_rate> <repair_sample_rate> <max_len_rule> <occ_sample_rate> <alfa>"  << endl;
	    return 0;
	}
	string graphname(argv[1]);
	uint repair_sample_rate = transform(argv[3]);
	uint max_len_rule = transform(argv[4]);
	uint occ_sample_rate = transform(argv[5]);
	double alfa;
	sscanf(argv[6],"%lf",&alfa);


	//Loggin stream
	string logname;
	logname.append("log/");	
	logname.append(argv[1]);logname.append("_");
	logname.append(argv[2]);logname.append("_");
	logname.append(argv[3]);logname.append("_");
	logname.append(argv[4]);logname.append("_");
	logname.append(argv[5]);logname.append("_");
	logname.append(argv[6]);logname.append(".log");

	ofstream logstream(logname.c_str(), ofstream::out);
    logstream.setf(ios::fixed,ios::floatfield);
	logstream.precision(2);



	printf(" [MSG] Loading Graph\n");
	string repair("Crawls/"+graphname+"/"+graphname+".repair"+argv[2]);
	ifstream repairstream (repair.c_str(), ifstream::in);
	if(!repairstream.good()){
		cout << "File not found: " << repair << endl;
		return 0;
	}

	logstream << graphname << " RRR sample rate " << argv[2] << endl << endl;	

	GraphSequenceTest *graph = GraphSequenceTest::load(repairstream);
	repairstream.close();
	graph->rpwvt_conf(repair_sample_rate, max_len_rule, occ_sample_rate, alfa);

	size_t size = graph->getSize();
	printf(" [SIZ] Total Size [bytes]: %u\n", size);
	printf(" [SIZ] Total Size [MB]: %.2f\n", size/1024.0/1024);
	printf(" [SIZ] Bits/Edge: %.2f\n", size*8.0/graph->edgesCount());

	logstream << "[SIZE]" << endl;
	logstream << "Total Size [bytes]: " << size << endl;  
	logstream << "Total Size [MB]: " << size/1024.0/1024 << endl;
	logstream << "Bits/Edge: "<< size*8.0/graph->edgesCount() << endl << endl;


	printf(" [TEST] Testing direct neighbors\n");	
	int repetitions = 10;
	double totaltime;
	vector<double> times;

	logstream << "[Direct Neighbors]" << endl;
	logstream << "Node\tTime [μs/edge]" << endl;
	totaltime = 0;
	for(uint r = 0; r < repetitions; ++r ){
		uint node = rand() % graph->nodesCount();

		double time = graph->testDirectNeighbors(node, logstream);
		times.push_back(time);
		totaltime += time;

		logstream << node << "\t"<< TOMICRO(time) << endl;

	}
	sort(times.begin(), times.end());

	logstream << "Mean [µs/edge]: " << TOMICRO(totaltime / times.size());
	logstream << ", Median [µs/edge]:" << TOMICRO(times[times.size()/2]);
	logstream << endl << endl;

	printf("Mean [µs/edge]: %.2f, Median [µs/edge]:%.2f\n", TOMICRO(totaltime / times.size()), TOMICRO(times[times.size()/2]));


	printf(" [TEST] Testing reverse neighbors\n");

	logstream << "[Reverse Neighbors]" << endl;
	logstream << "Node\tTime [μs/edge]" << endl;
	totaltime = 0;
	times.clear();
	for(uint r = 0; r < repetitions; ++r ){
		uint node = rand() % graph->nodesCount();

		double time = graph->testDirectNeighbors(node, logstream);
		times.push_back(time);
		totaltime += time;

		logstream << node << "\t"<< TOMICRO(time) << endl;

	}
	sort(times.begin(), times.end());

	logstream << "Mean [µs/edge]: " << TOMICRO(totaltime / times.size());
	logstream << ", Median [µs/edge]:" << TOMICRO(times[times.size()/2]);
	logstream << endl;


	printf("Mean [µs/edge]: %.2f, Median [µs/edge]:%.2f\n", TOMICRO(totaltime / times.size()), TOMICRO(times[times.size()/2]));
}
