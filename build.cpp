/**
 * Builds and save a Graph represented by
 * a RePair WaveletTree. It loads a graph from
 * Crawls/<graphname>/<graphname>.plain and use
 * a fixed sampling rate for the RRR bitsequence
 * on the wavelettree.
 * */


#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <string>

#include <BitSequenceBuilder.h>
#include <Mapper.h>
#include "RePairWaveletTreeBuilder.h"
#include "GraphSequence.h"


using namespace std;





int main (int argc, char **argv)
{ 
	if(argc != 3) {
	    cout << "usage: " << argv[0] << " <graphname> <rrr_sample_rate>"  << endl;
	    return 0;
	}

	string graphname(argv[1]);
	uint sample_rate = transform(string(argv[2]));


	printf(" [MSG] Building adjacency list\n");
	
	string inputfile("Crawls/"+graphname+"/"+graphname+".plain");
	ifstream ifs (inputfile.c_str(), ifstream::in);

	uint nodes = loadValue<uint>(ifs); //First Int is the number of nodes

	uint **adjacencyList = new uint*[nodes];

	for(uint i = 0; i < nodes; ++i){
		uint d = loadValue<uint>(ifs);  
		adjacencyList[i] = new uint[d+1];
		adjacencyList[i][0] = d;            //Outdegree before the neighbors
		for(uint j = 1; j <= d; ++j){
			uint v = loadValue<uint>(ifs);  //The neighbors
			adjacencyList[i][j] = v;
		}
	}
	ifs.close();


	printf(" [MSG] Building graph...\n");

	BitSequenceBuilder *bsb = new BitSequenceBuilderRRR(sample_rate);
	Mapper *mp = new MapperNone();
	SequenceBuilder *sb = new RePairWaveletTreeBuilder(bsb, mp, true);
	GraphSequence *graph = new GraphSequence(adjacencyList, nodes, sb, true);


	size_t size = graph->getSize();
	printf(" [SIZ] Total Size [bytes]: %u\n", size);
	printf(" [SIZ] Total Size [MB]: %.2f\n", size/1024.0/1024);
	printf(" [SIZ] Bits/Edge: %.2f\n", size*8.0/graph->edgesCount());

	string outputfile("Crawls/"+graphname+"/"+graphname+".repair"+argv[2]);
	ofstream ofs (outputfile.c_str(), ofstream::out);
	graph->save(ofs);

}
