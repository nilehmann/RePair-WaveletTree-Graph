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
#include "GraphSequence.h"


namespace cds_static
{
	GraphSequence::GraphSequence(uint **_adjacencyList, size_t nodes, SequenceBuilder *sb, bool debug){
		this->nodes = nodes;

		this->edges = 0;
		this->falseEdges = 0;
		for(uint n = 0; n < nodes; ++n){
			this->edges += _adjacencyList[n][0];
			if(_adjacencyList[n][0] == 0){
				this->falseEdges++;
			}
		}

		uint *seq = new uint[this->edges + this->falseEdges];

		BitString bitstring(this->edges + this->falseEdges);
		uint k = 0;
		for(uint n = 0; n < nodes; ++n){
			bitstring.setBit(k,1);
			uint neighbors = _adjacencyList[n][0];
			for(uint i = 0; i < neighbors; ++i){
				if(i > 0)
					bitstring.setBit(k,0);
				seq[k++] = _adjacencyList[n][i+1];
			}
			if(neighbors == 0)
				seq[k++] = this->nodes+1;
			delete [] _adjacencyList[n];
		}
		delete [] _adjacencyList;

		if(debug)
			printf(" [GRPH] Building Bitsequence\n");
		this->bitseq = new BitSequenceRRR(bitstring);
		if(debug)
			printf(" [GRPH] Building Sequence\n");
		this->adjacencyList = sb->build(seq, this->edges + this->falseEdges);

	}

	GraphSequence::GraphSequence(){
		nodes = 0;
		edges = falseEdges = 0;
		adjacencyList = NULL;
		bitseq = NULL;
	}

	GraphSequence::~GraphSequence(){
		delete adjacencyList;
		delete bitseq;
	}
	

	vector<uint> *GraphSequence::getNeighbors(uint n) const{
		size_t pos = bitseq->select1(n+1);
		size_t next = bitseq->select1(n+2);
		size_t len = (next == -1?bitseq->getLength() : next) - pos;

		vector<uint> *neighbors = new vector<uint>(len);


		for(int i = 0; i < len; ++i){
			(*neighbors)[i] = adjacencyList->access(pos + i);
			if(neighbors->at(i) == this->nodes +1)
				neighbors->clear();
		}

		return neighbors;
	}
	vector<uint> *GraphSequence::getReverseNeighbors(uint n) const{
		size_t len = adjacencyList->rank(n, adjacencyList->getLength());

		vector<uint> *neighbors = new vector<uint>(len);

		for(int i = 0; i < len; ++i){
			(*neighbors)[i] = bitseq->rank1(adjacencyList->select(n, i+1)) - 1;
		}

		return neighbors;
	}

	size_t GraphSequence::getSize() const{
		size_t ptrs = sizeof(Sequence)+ sizeof(BitSequence);
		size_t data = sizeof(long long)*2+sizeof(size_t);
		return adjacencyList->getSize()+bitseq->getSize()+ptrs+data;
	}
	

	void GraphSequence::save(ofstream & fp) const
	{
		saveValue<size_t>(fp, nodes);
		saveValue<long long>(fp, edges);
		saveValue<long long>(fp, falseEdges);
		adjacencyList->save(fp);
		bitseq->save(fp);
	}

	GraphSequence * GraphSequence::load(ifstream & fp) {
		GraphSequence *ret = new GraphSequence();
		ret->nodes = loadValue<size_t>(fp);
		ret->edges = loadValue<long long>(fp);
		ret->falseEdges = loadValue<long long>(fp);
		ret->adjacencyList = RePairWaveletTree::load(fp);

		ret->bitseq = BitSequence::load(fp);


		return ret;
	}
			
};
