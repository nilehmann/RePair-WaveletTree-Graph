/* WaveletTreeNoptrs.h
 * Copyright (C) 2013, Nicolás Lehmann 
 *
 * WaveletTree implementation that compress the bitmaps using RePair.
 * The implementation test at each level if the repair bitmap use
 * less space than another bitmap implementation.
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

#ifndef _REPAIR_WVTREE_H
#define _REPAIR_WVTREE_H

#include <iostream>
#include <cassert>
#include <libcdsBasics.h>
#include <BitSequence.h>
#include <BitSequenceBuilder.h>
#include <Sequence.h>
#include <Mapper.h>
#include <bpe.h>

using namespace std;

namespace cds_static
{

	class RePairWaveletTree : public Sequence
	{
		public:

			/** 
			 * Builds a Wavelet Tree for the string
			 * pointed by <var>symbols</var> assuming its length
			 * equals <var>n</var>.
			 * @param symbols		Sequence of symbols.
			 * @param n 			Length of the sequence.
			 * @param bmb			Builder for the bitmaps not compressed with RePair.
			 * @param deleteSymbols	True if you want to delete de sequence after construction.
			 * */
			RePairWaveletTree(uint * symbols, size_t n, BitSequenceBuilder * bmb, Mapper * am, bool deleteSymbols = false);


			/** Destroys the Wavelet Tree */
			virtual ~RePairWaveletTree();

			virtual size_t rank(uint symbol, size_t pos) const;
			virtual size_t select(uint symbol, size_t j) const;
			virtual uint access(size_t pos) const;
			virtual uint access(size_t pos, size_t &r) const;
			virtual size_t getSize() const;

			/* find the q-th smallest element in T[l..r] */
			virtual uint quantile(size_t left,size_t right,uint q);

			/* find the q-th smallest element in T[l..r] and return the freq */
			pair<uint32_t,size_t> quantile_freq(size_t left,size_t right,uint q);

			virtual size_t count(uint symbol) const;

			virtual void save(ofstream & fp) const;
			static RePairWaveletTree * load(ifstream & fp);


			/*On demand configuration methods*/
			void re_sample_repair(uint sample_rate);
			void re_set_max_rule_len(uint len);
			/*Must be called after re_sample_repair and re_set_max_rule*/
			void re_set_alpha_factor(double alpha);
			void re_set_occ_sample(uint sample_rate);
			void turn_occ_array();
			void turn_occ_bitmap();


		protected:
			RePairWaveletTree();

			Mapper * am;
			BitSequence **bitseq;

			BPE **repair;
			BitSequence **nonrepair;

			bool arrayocc;
			BitSequence *occ;
			uint *OCC;

			BitSequence *buildBitmaps(uint *bm, uint n, BitSequenceBuilder *bmb, uint level) const;

			/** Length of the string. */
			size_t n;

			/** Height of the Wavelet Tree. */
			uint height, max_v;

			/** Obtains the maximum value from the string
			 * symbols of length n */
			uint max_value(uint * symbols, size_t n);

			/** How many bits are needed to represent val */
			uint bits(uint val);

			/** Returns true if val has its ind-th bit set
			 * to one. */
			bool is_set(uint val, uint ind) const;

			/** Recursive function for building the bitmaps of the Wavelet Tree. */
			void build_level(uint **bm, uint *symbols, uint level, uint length, uint offset);

	};
};
#endif
