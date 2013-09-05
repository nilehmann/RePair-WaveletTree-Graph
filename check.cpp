/**
 * Check the implementation of the graph
 * */


#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <string>

#include <BitSequenceBuilder.h>
#include <Mapper.h>
#include "GraphSequence.h"


using namespace std;
int main (int argc, char **argv)
{ 	
	if(argc != 3) {
	    cout << "usage: " << argv[0] << " <graphname> <rrr_sample_rate>"  << endl;
	    return 0;
	}
	string graphname(argv[1]);



	printf(" [MSG] Loading Graph\n");
	string repair("Crawls/"+graphname+"/"+graphname+".repair"+argv[2]);
	ifstream repairstream (repair.c_str(), ifstream::in);
	GraphSequence *graph = GraphSequence::load(repairstream);
	repairstream.close();



	printf(" [MSG] Building adjacency list\n");
	string plain("Crawls/"+graphname+"/"+graphname+".plain");
	ifstream plainstream (plain.c_str(), ifstream::in);

	uint nodes = loadValue<uint>(plainstream); //First Int is the number of nodes
	
	uint **adjacencyList = new uint*[nodes];

	for(uint i = 0; i < nodes; ++i){
		uint d = loadValue<uint>(plainstream);  //For each node the number of neighbors
		adjacencyList[i] = new uint[d+1];
		adjacencyList[i][0] = d;
		for(uint j = 1; j <= d; ++j){
			uint v = loadValue<uint>(plainstream);  //The neighbors
			adjacencyList[i][j] = v;
		}
	}
	plainstream.close();

	size_t size = graph->getSize();
	printf(" [SIZ] Total Size [bytes]: %u\n", size);
	printf(" [SIZ] Total Size [MB]: %.2f\n", size/1024.0/1024);
	printf(" [SIZ] Bits/Edge: %.2f\n", size*8.0/graph->edgesCount());




	printf(" [CHK] Cheking direct neighbors\n");
	
	for(uint node = 0; node < nodes; node += 50  ){
		if(node %10000 == 0){
			printf(".");
			fflush(stdout);
		}
		vector<uint> *neighbors = graph->getNeighbors(node);
		for(uint i = 0; i < neighbors->size(); ++i){
			uint neighbor = neighbors->at(i);
			if(neighbor != adjacencyList[node][i+1]){
				printf("Error at node %u, neighbors[%u]=%u\n", node, i, neighbor);
			}
		}

	}
	printf("\n");



	printf(" [CHK] Checking reverse neighbors\n");

	for(uint node = 0; node < nodes; node += 50){
		if(node %10000 == 0){
			printf(".");
			fflush(stdout);
		}
		vector<uint> *neighbors = graph->getReverseNeighbors(node);
		for(uint i = 0; i < neighbors->size(); ++i){
			uint neighbor = neighbors->at(i);
			bool b = false;
			for(uint j = 0; j < adjacencyList[neighbor][0]; ++j)
				if(node == adjacencyList[neighbor][j+1]){
					b=true;
					break;
				}
			if(!b)
				printf("Error at node %u, neighbors[%u]=%u\n", node, i, neighbor);
		}
	}
	printf("\n");
	
}
