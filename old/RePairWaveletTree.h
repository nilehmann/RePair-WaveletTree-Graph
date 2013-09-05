/* WaveletTreeNoptrs.h
 * Copyright (C) 2013, Nicolás Lehmann 
 *
 * WaveletTree implementation that compress the bitmaps using RePair.
 * The implementation test at each level if the repair bitmap use
 * less space than another bitmap implementation.
 * 
 * Old implementation that use an int array for the occurrences
 * of the symbols instead of a Bitmap
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

	class RePairWaveletTree : public WaveletTreeNoptrs
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
			RePairWaveletTree(uint * symbols, size_t n, BitSequenceBuilder * bmb, Mapper *am, bool deleteSymbols = false);


			/** Destroys the Wavelet Tree */
			virtual ~RePairWaveletTree();

			virtual void save(ofstream & fp) const;
			static RePairWaveletTree * load(ifstream & fp);

			virtual size_t getSize() const;

		
		protected:
			RePairWaveletTree();
			BitSequence *chooseBitmap(uint *bm, uint n, BitSequenceBuilder *bmb) const;
	};
};
#endif
