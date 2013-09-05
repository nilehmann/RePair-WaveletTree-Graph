/* WaveletTreeNoptrs.h
 * Copyright (C) 2013, Nicolás Lehmann.
 *
 * Graph implementation using a sequence representation
 * for the concatenation of the adjacency lists
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */

#include <ctime>

#include "GraphSequenceTest.h"


namespace cds_static
{

	GraphSequenceTest::GraphSequenceTest(){
	}

	GraphSequenceTest::~GraphSequenceTest(){
	}
	

	double GraphSequenceTest::testDirectNeighbors(uint n, ofstream &fp) const{
		clock_t begin = clock();

		size_t pos = bitseq->select1(n+1);
		size_t next = bitseq->select1(n+2);
		size_t len = (next == -1?bitseq->getLength() : next) - pos;


		int repetitions = 1024/len > 1 ? 1024/len : 1;



		for(int r = 0; r < repetitions; ++r){
			for(int i = 0; i < len; ++i)
				adjacencyList->access(pos + i);
		}


		clock_t end = clock();
		double elapsed_secs = ((double)(end - begin))/ CLOCKS_PER_SEC/repetitions;		

		if(elapsed_secs < 0.00001)
			fprintf(stderr,"To fast repetitions: %d, len: %d, secs: %.10f\n", repetitions, len,elapsed_secs);
		return elapsed_secs/len;
	}
	double GraphSequenceTest::testReverseNeighbors(uint n, ofstream &fp) const{
		clock_t begin = clock();

		size_t len = adjacencyList->rank(n, adjacencyList->getLength());


		int repetitions = 1024/len > 1 ? 1024/len : 1;


		for(int r = 0; r < repetitions; ++r){
			for(int i = 0; i < len; ++i)
				bitseq->rank1(adjacencyList->select(n, i+1)) - 1;
		}

		clock_t end = clock();
		double elapsed_secs = ((double)(end - begin))/ CLOCKS_PER_SEC/repetitions;	
		
		if(elapsed_secs < 0.00001)
			fprintf(stderr,"To fast repetitions: %d, len: %d, secs: %.10f\n", repetitions, len,elapsed_secs);

		return elapsed_secs/len;
	}

	GraphSequenceTest * GraphSequenceTest::load(ifstream & fp) {
		GraphSequenceTest *ret = new GraphSequenceTest();
		ret->nodes = loadValue<size_t>(fp);
		ret->edges = loadValue<long long>(fp);
		ret->falseEdges = loadValue<long long>(fp);
		ret->adjacencyList = RePairWaveletTree::load(fp);

		ret->bitseq = BitSequence::load(fp);


		return ret;
	}

	void GraphSequenceTest::rpwvt_conf(uint repair_sample, uint max_len, uint occ_sample, double alpha){
		RePairWaveletTree *wvt = (RePairWaveletTree *) this->adjacencyList;
		wvt->re_sample_repair(repair_sample);
		wvt->re_set_max_rule_len(max_len);
		if(occ_sample == 0)
			wvt->turn_occ_array();
		else{
			wvt->turn_occ_bitmap();
			wvt->re_set_occ_sample(occ_sample);
		}	

		wvt->re_set_alpha_factor(alpha);
	}
};
