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

#include "RePairWaveletTree.h"

namespace cds_static
{
	RePairWaveletTree::RePairWaveletTree(uint * symbols, size_t n, BitSequenceBuilder * bmb, Mapper * am, bool deleteSymbols)  {
		this->n = n;
		this->length = n;
		this->am = am;
		am->use();
		for(uint i = 0; i < n; i++)
			symbols[i] = am->map(symbols[i]);
		max_v = max_value(symbols, n);
		height = bits(max_v);

		OCC = new uint[max_v + 2];
		for (uint i = 0; i <= max_v + 1; i++)
			OCC[i] = 0;
		for (uint i = 0; i < n; i++)
			OCC[symbols[i] + 1]++;

		uint to_add = 0;
		for (uint i = 1; i <= max_v + 1; i++)
			if (OCC[i] == 0) to_add++;

		uint * new_symb = new uint[n + to_add];
		for (uint i = 0; i < n; i++)
			new_symb[i] = symbols[i];

		if (deleteSymbols) {
			delete [] symbols;
			symbols = 0;
		}

		to_add = 0;
		for (uint i = 1; i <= max_v + 1; i++)
		if (OCC[i] == 0) {
			OCC[i]++;
			new_symb[n + to_add] = i - 1;
			to_add++;
		}

		uint new_n = n + to_add;
		for(uint i = 1;i <= max_v + 1; i++)
			OCC[i] += OCC[i - 1];
		this->n = new_n;

		uint **_bm = new uint*[height];
		for(uint i = 0; i < height; i++) {
			_bm[i] = new uint[new_n / W + 1];
			for(uint j = 0; j < new_n / W + 1; j++)
				_bm[i][j] = 0;
		}
		build_level(_bm, new_symb, 0, new_n, 0);
		bitstring = new BitSequence*[height];
		for (uint i = 0; i < height; i++) {
			bitstring[i] = chooseBitmap(_bm[i], new_n, bmb);
			delete [] _bm[i];
		}
		delete [] _bm;

		if (!deleteSymbols)
			for(uint i = 0; i < n; i++)
				symbols[i] = am->unmap(symbols[i]);

	}

	RePairWaveletTree::RePairWaveletTree() {
		bitstring = NULL;
		OCC = NULL;
		am = NULL;
	}

	BitSequence * RePairWaveletTree::chooseBitmap(uint *bm, uint n, BitSequenceBuilder *bmb) const{
		ulong *seq = new ulong[n];
    	for(int i=0;i<n;i++)
        	seq[i] = bitget(bm, i);


		BPE *repair = new BPE(seq,n, false,-1,1024,1);
		int mp = repair->getMaxProf();
		int len = repair->getLargoMaxReal();
		repair->re_set_SOME_rules2(mp, 4);
//		repair->re_set_max_rule_len(len/2);

		BitSequence *bitseq = bmb->build(bm, n);

		size_t repairSize = repair->getSize();
		size_t bitseqSize = bitseq->getSize();
		
		printf("%d %d", repairSize, bitseqSize);

		if(repairSize < bitseqSize){
			printf(" REPAIR\n");
			delete bitseq;
			return repair;
		}
		else{
			printf(" RRR\n");
			delete repair;
			return bitseq;
		}
	}

	RePairWaveletTree::~RePairWaveletTree() {
	}

	size_t RePairWaveletTree::getSize() const {
		size_t ptrs = height * sizeof(BitSequence*);
		printf("OCC size %u\n", sizeof(uint) * (max_v + 2));
		size_t bytesBitstrings = 0;
		for(uint i = 0; i < height; i++)
			bytesBitstrings += bitstring[i]->getSize();
		return bytesBitstrings + sizeof(uint) * (max_v + 2) + ptrs;
	}

	void RePairWaveletTree::save(ofstream & fp) const
	{
		saveValue<size_t>(fp,n);
		saveValue(fp, max_v);
		saveValue(fp, height);
		am->save(fp);
		for (uint i = 0; i < height; i++)
			bitstring[i]->save(fp);
		saveValue<uint>(fp, OCC, max_v + 2);
	}

	RePairWaveletTree * RePairWaveletTree::load(ifstream & fp) {
		RePairWaveletTree * ret = new RePairWaveletTree();
		ret->n = loadValue<size_t>(fp);
		ret->length = ret->n;
		ret->max_v = loadValue<uint>(fp);
		ret->height = loadValue<uint>(fp);
		ret->am = Mapper::load(fp);
		if (ret->am == NULL) {
			delete ret;
			return NULL;
		}
		ret->am->use();
		ret->bitstring = new BitSequence*[ret->height];
		for(uint i = 0; i < ret->height; i++)
			ret->bitstring[i] = NULL;
		for(uint i = 0; i < ret->height; i++) {
			uint r = loadValue<uint>(fp);
			size_t pos = fp.tellg();
			fp.seekg(pos-sizeof(uint));
			if(r == REPAIR_HDR)
				ret->bitstring[i] = BPE::load(fp);
			else
				ret->bitstring[i] = BitSequence::load(fp);

			if (ret->bitstring[i] == NULL) {
				delete ret;
				return NULL;
			}
		}
		ret->OCC = loadValue<uint>(fp, ret->max_v + 2);
		return ret;
	}


};
