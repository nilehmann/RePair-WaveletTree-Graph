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


#ifndef GRAPH_SEQUENCE_H
#define GRAPH_SEQUENCE_H

#include <libcdsBasics.h>

#include <Mapper.h>
#include <SequenceBuilder.h>
#include <BitSequenceBuilder.h>
#include "RePairWaveletTree.h"

namespace cds_static
{
	class GraphSequence
	{
		public:
			/**
			 * Builds a graph from a adjacency list represented by
			 * <it>nodes</it> sequences. Each sequence has to start with the number
			 * of adjacent nodes and then the corresponding neighbors.
			 * If a node has zero neighbors his respective sequence must be
			 * one containing only a zero
			 * @param adjacencyList	Adjacency list of the graph.
			 * @param nodes			Number of nodes
			 * @param sb			Builder for sequence that represent the adjacency sequence
			 * @param debug			True if you want to print debug messages
			 * */
			GraphSequence(uint **adjacencyList, size_t nodes, SequenceBuilder *sb, bool debug = false);
			virtual ~GraphSequence();
			
			vector<uint> *getNeighbors(uint n) const;
			vector<uint> *getReverseNeighbors(uint n) const;

			size_t neighborsCount(uint n) const;
			size_t reverseNeighborsCount(uint n) const;
		
			size_t nodesCount() const{return nodes;}
			size_t edgesCount() const{return edges;}

			virtual void save(ofstream & fp) const;
			static GraphSequence * load(ifstream & fp);

			
			virtual size_t getSize() const;

		protected:
			Sequence *adjacencyList;
			BitSequence *bitseq;
			size_t nodes;
			long long edges;
			long long falseEdges;

			GraphSequence();
			
	};
};
#endif
