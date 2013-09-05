/* WaveletTreeNoptrs.h
 * Copyright (C) 2013, Nicolás Lehmann.
 *
 * Class for testing query times for graph 
 * based on repair waveletree 
 *
 */


#ifndef GRAPH_SEQUENCE_TEST_H
#define GRAPH_SEQUENCE_TEST_H


#include <libcdsBasics.h>
#include <Mapper.h>
#include <SequenceBuilder.h>
#include <BitSequenceBuilder.h>
#include "GraphSequence.h"
#include "RePairWaveletTree.h"

namespace cds_static
{
	class GraphSequenceTest : public GraphSequence
	{
		public:
			virtual ~GraphSequenceTest();
			
			/**
			 * Retrieves all direct neighbors of node <var>n</var>
			 * and return the average tieme required per edge.
			 * @param n Node identifier
			 * @param fp Stream to write a log
			 * */
			double testDirectNeighbors(uint n, ofstream &fp) const;
			/**
			 * Retrieves all reverse neighbors of node <var>n</var>
			 * and return the average tieme required per edge.
			 * @param n Node identifier
			 * @param fp Stream to write a log
			 * */
			double testReverseNeighbors(uint n, ofstream &fp) const;

			static GraphSequenceTest * load(ifstream & fp);

			/*On demand configuration*/
			void rpwvt_conf(uint repair_sample, uint max_len, uint occ_sample, double alpha);

		protected:
			GraphSequenceTest();
	};
};
#endif
