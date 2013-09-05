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

#include "RePairWaveletTree.h"

namespace cds_static
{


	RePairWaveletTree::RePairWaveletTree(uint * symbols, size_t n, BitSequenceBuilder * bmb, Mapper * am, bool deleteSymbols) : Sequence(n) {
		this->n = n;
		this->am = am;
		am->use();
		for(uint i = 0; i < n; i++)
			symbols[i] = am->map(symbols[i]);
		max_v = max_value(symbols, n);
		height = bits(max_v);

		this->OCC = new uint[max_v + 1];
		for (uint i = 0; i <= max_v; i++)
			OCC[i] = 0;
		//Count occurrence of each symbol
		for (uint i = 0; i < n; i++)
			OCC[symbols[i]]++;

		//Count the number of symbols not mentioned
		uint to_add = 0;
		for (uint symb = 0; symb <= max_v; symb++)
			if (OCC[symb] == 0) to_add++;

		//Copy the symbols to new sequence
		uint * new_symb = new uint[n + to_add];
		for (uint i = 0; i < n; i++)
			new_symb[i] = symbols[i];

		if (deleteSymbols) {
			delete [] symbols;
			symbols = 0;
		}

		//Add the symbols not mentioned
		to_add = 0;
		for (uint symb = 0; symb <= max_v; symb++)
		if (OCC[symb] == 0) {
			OCC[symb]++;
			new_symb[n + to_add] = symb;
			to_add++;
		}

		uint new_n = n + to_add;
		this->n = new_n;


		//Build occ bitmap
		BitString oc(new_n);
		uint pos = 0;
		for(uint symb = 0; symb <= max_v; ++symb){
			pos += OCC[symb];
			if(symb == 0) pos--;
			oc.setBit(pos, 1);
			OCC[symb] += symb > 0? OCC[symb-1] : 0;
		}
		occ = new BitSequenceRRR(oc,32);
		this->arrayocc = false;

		uint **_bm = new uint*[height];
		for(uint l = 0; l < height; l++) {
			_bm[l] = new uint[new_n / W + 1];
			for(uint i = 0; i < new_n / W + 1; i++)
				_bm[l][i] = 0;
		}
		build_level(_bm, new_symb, 0, new_n, 0);
		bitseq = new BitSequence*[height];
		repair = new BPE*[height];
		nonrepair = new BitSequence*[height];
		for (uint l = 0; l < height; l++) {
			bitseq[l] = buildBitmaps(_bm[l], new_n, bmb, l);
			delete [] _bm[l];
		}
		delete [] _bm;
		fprintf(stderr,"4\n");

		if (!deleteSymbols)
			for(uint i = 0; i < n; i++)
				symbols[i] = am->unmap(symbols[i]);

	}

	BitSequence * RePairWaveletTree::buildBitmaps(uint *bm, uint n, BitSequenceBuilder *bmb, uint level) const{
		ulong *seq = new ulong[n];
    	for(int i=0;i<n;i++)
        	seq[i] = bitget(bm, i);


		repair[level] = new BPE(seq,n, false,-1,32,5);
		int len = repair[level]->getLargoMaxReal();
		repair[level]->re_set_max_rule_len(len);

		nonrepair[level] = bmb->build(bm, n);

		size_t repairSize = repair[level]->getSize();
		size_t nonrepairSize = nonrepair[level]->getSize();
		
		printf("%u %u", repairSize, nonrepairSize);
		if(repairSize < nonrepairSize){
			printf(" REPAIR\n");
			return repair[level];
		}
		else{
			printf(" RRR\n");
			return nonrepair[level];
		}
	}

	RePairWaveletTree::RePairWaveletTree():Sequence(0) {
		bitseq = NULL;
		occ = NULL;
		am = NULL;
	}

	RePairWaveletTree::~RePairWaveletTree() {
		for (uint i = 0; i < height; i++){
			delete bitseq[i];
			delete repair[i];
			delete repair[i];
		}
		delete [] bitseq;
		delete [] repair;
		delete [] nonrepair;
		delete occ;
		delete [] OCC;
		if (am)
			am->unuse();
	}

	

	inline uint get_start(uint symbol, uint mask) {
		return symbol & mask;
	}

	inline uint get_end(uint symbol, uint mask) {
		return get_start(symbol, mask) + ~mask + 1;
	}

	bool RePairWaveletTree::is_set(uint val, uint ind) const
	{
		assert (ind < height);
		return (val & (1 << (height - ind - 1))) != 0;
	}

	uint RePairWaveletTree::access(size_t pos) const
	{
		uint ret = 0;
		
		size_t start = 0;
		for (uint level = 0; level < height; level++) {
			size_t onesToPos, onesBefore = 0;
			if (start > 0)
				onesBefore = bitseq[level]->rank1(start - 1);
			
			if (bitseq[level]->access(pos, onesToPos)) {
				ret |= (1 << (height - level - 1));
				pos = onesToPos - onesBefore - 1;
				if(arrayocc)
					start = ret == 0 ? 0 : OCC[ret-1];
				else
					start = ret == 0 ? 0 : occ->select1(ret) + 1;
				pos += start;
			} else {
				pos = onesToPos - 1 + onesBefore;
			}
		}

		return am->unmap(ret);
	}

	uint RePairWaveletTree::access(size_t pos, size_t &r) const
	{
		uint ret = 0;

		size_t start = 0;
		for (uint level = 0; level < height; level++) {
			size_t onesToPos, onesBefore=0;
			if (start > 0) 
				onesBefore = bitseq[level]->rank1(start-1);
			
			if(bitseq[level]->access(pos, onesToPos)) {
				ret |= (1 << (height - level - 1));
				r = onesToPos - onesBefore;
				if(arrayocc)
					start = ret == 0 ? 0 : OCC[ret-1];
				else
					start = ret == 0 ? 0 : occ->select1(ret) + 1;
				pos = r - 1 + start;
			}
			else {
				pos = onesToPos - 1 + onesBefore;
				r = pos + 1 - start;
			}
		}

		return am->unmap(ret);
	}

	size_t RePairWaveletTree::rank(uint symbol, size_t pos) const
	{
		symbol = am->map(symbol);

		size_t start = 0;
		size_t count = 0;
		
		for(uint level = 0; level < height; level++) {
			
			uint masked = (symbol >> (height - level - 1)) << (height - level - 1);
			
			size_t before = 0;
			if (start > 0)
				before = bitseq[level]->rank1(start - 1);
			
			if (is_set(symbol, level)) {
				count = bitseq[level]->rank1(pos) - before;
				if(arrayocc)
					start = masked == 0 ? 0 : OCC[masked-1];
				else
					start = masked == 0 ? 0 : occ->select1(masked) + 1;
				pos = count + start - 1;
			} else {
				count = pos - start + before - bitseq[level]->rank1(pos) + 1;
				masked += (1 << (height - level - 1)); 
				pos = count + start - 1;
			}

			if (count == 0) return 0;
		}
		return count;
	}

	size_t RePairWaveletTree::select(uint symbol, size_t j) const
	{
		symbol = am->map(symbol);
		
		uint mask = (1 << height) - 2;

		size_t pos = j;

		for (int level = height - 1; level >= 0; level--) {
			
			size_t start = get_start(symbol, mask);

			if(arrayocc)
				start = start == 0? 0: OCC[start-1];
			else
				start = start == 0 ? 0 : occ->select1(start) + 1;



			uint ones_start = 0;
			if (start > 0)
				ones_start = bitseq[level]->rank1(start - 1);


			if (is_set(symbol,level)) {
				pos = bitseq[level]->select1(ones_start + pos) - start + 1;
			} else {
				pos = bitseq[level]->select0(start - ones_start + pos) - start + 1;
			}


			mask <<= 1;
		}

		return pos - 1;



	}

	size_t RePairWaveletTree::getSize() const
	{
		size_t ptrs = sizeof(RePairWaveletTree) + height * sizeof(Sequence*);
		size_t bytesBitstrings = 0;
		for(uint i = 0; i < height; i++)
			bytesBitstrings += bitseq[i]->getSize();
		printf("Size occ: %d, OCC: %d\n", occ->getSize(), sizeof(uint)*(max_v+1));
		return bytesBitstrings + (arrayocc ? sizeof(uint)*(max_v+1) : occ->getSize());
	}

	void RePairWaveletTree::build_level(uint **bm, uint *symbols, uint level, uint length, uint offset) {
		if (level == height) {
			delete [] symbols;
			return;
		}

		uint cleft = 0;
		for (size_t i = 0; i < length; i++)
			if (!is_set(symbols[i],level))
				cleft++;

		uint cright = length - cleft;

		uint *left = new uint[cleft];
		uint *right = new uint[cright];
		cleft = cright = 0;
		for (size_t i = 0; i < length; i++) {
			if (!is_set(symbols[i], level)) {
				left[cleft++] = symbols[i];
				bitclean(bm[level], offset + i);
			} else {
				right[cright++] = symbols[i];
				bitset(bm[level], offset + i);
			}
		}

		delete [] symbols;
		symbols = NULL;

		build_level(bm, left, level + 1, cleft, offset);
		left = NULL;			 // Gets deleted in recursion.
		build_level(bm, right, level + 1, cright, offset + cleft);
		right = NULL;			 // Gets deleted in recursion.
	}

	uint RePairWaveletTree::max_value(uint *symbols, size_t n) {
		uint max_v = 0;
		for (size_t i = 0; i < n; i++)
			max_v = max(symbols[i], max_v);
		return max_v;
	}



	uint RePairWaveletTree::bits(uint val) {
		uint ret = 0;
		while (val!=0) {
			ret++;
			val >>= 1;
		}
		return ret;
	}

	size_t RePairWaveletTree::count(uint symbol) const
	{
		uint mapped = am->map(symbol);
		size_t a = mapped == 0? -1 : occ->select1(mapped);
		size_t b = occ->select1(mapped+1);
		return b-a;
	}

	uint RePairWaveletTree::quantile(size_t left,size_t right,uint q) {
		pair<uint,size_t> res = quantile_freq(left,right,q);
		return res.first;
	}

	pair<uint32_t,size_t> RePairWaveletTree::quantile_freq(size_t left,size_t right,uint q) {
		/* decrease q as the smallest element q=1 is
		 * found by searching for 0 */
		q--;

		assert( right >= left );
		assert( (right-left+1) >= q );
		assert( right < length );

		uint sym = 0;
		uint freq = 0;
		uint level = 0;
		size_t start = 0, end = n-1;
		size_t before;
		BitSequence* bs;

		while(level<height) {
			bs = bitseq[level];

			/* calc start of level bound */
			if(start == 0) before = 0;
			else before = bs->rank1(start-1);

			/* number of 1s before T[l..r] */
			size_t rank_before_left = bs->rank1(start+left-1);
			/* number of 1s before T[r] */
			size_t rank_before_right = bs->rank1(start+right);
			/* number of 1s in T[l..r] */
			size_t num_ones = rank_before_right - rank_before_left;
			/* number of 0s in T[l..r] */
			size_t num_zeros = (right-left+1) - num_ones;

			/* if there are more than q 0s we go right. left otherwise */
			if(q >= num_zeros) { /* go right */
				freq = num_ones; /* calc freq */
				/* set bit to 1 in sym */
				sym = 1 << (height - level - 1); //set(sym,level);
				/* number of 1s before T[l..r] within the current node */
				left = rank_before_left - before;
				/* number of 1s in T[l..r] */
				right = rank_before_right - before - 1;
				q = q - num_zeros;
				/* calc starting pos of right childnode */
				start = end - (bs->rank1(end)-before) + 1;
			}					 /* go left q = q // sym == sym */
			else {
				freq = num_zeros;/* calc freq */
				/* number of zeros before T[l..r] within the current node */
				left = left - (rank_before_left - before);
				/* number of zeros in T[l..r] + left bound */
				right = right - (rank_before_right - before);
				/* calc end pos of left childnode */
				end = end - (bs->rank1(end) - before);
			}
			level++;
		}

		/* unmap symbol */
		return pair<uint,size_t>(am->unmap(sym),static_cast<uint>(freq));
	}

	void RePairWaveletTree::save(ofstream & fp) const
	{
		saveValue<size_t>(fp,n);
		saveValue(fp, max_v);
		saveValue(fp, height);
		am->save(fp);
		for (uint i = 0; i < height; i++){
			repair[i]->save(fp);
			nonrepair[i]->save(fp);
		}
		occ->save(fp);
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
		ret->bitseq = new BitSequence*[ret->height];
		ret->repair = new BPE*[ret->height];
		ret->nonrepair = new BitSequence*[ret->height];
		for(uint i = 0; i < ret->height; i++){
			ret->bitseq[i] = NULL;			
			ret->repair[i] = NULL;		
			ret->nonrepair[i] = NULL;		
		}
		for(uint i = 0; i < ret->height; i++) {
			ret->repair[i] = BPE::load(fp);
			ret->nonrepair[i] = BitSequence::load(fp);

			size_t repairSize = ret->repair[i]->getSize();
			size_t nonrepairSize = ret->nonrepair[i]->getSize();

			if(repairSize < nonrepairSize)
				ret->bitseq[i] = ret->repair[i];
			else
				ret->bitseq[i] = ret->nonrepair[i];
		}
		ret->occ = BitSequence::load(fp);
		if (ret->occ == NULL) {
			delete ret;
			return NULL;
		}
		ret->OCC = loadValue<uint>(fp, ret->max_v + 2);
		ret->arrayocc = true;
		return ret;
	}

	void RePairWaveletTree::re_sample_repair(uint sample_rate) {
		for(uint l = 0; l < height; ++l)
			repair[l]->build_GN_samples(sample_rate);
	}
	
	void RePairWaveletTree::re_set_max_rule_len(uint len) {
		for(uint l = 0; l < height; ++l){
			uint maxlen = repair[l]->getLargoMaxReal();
			maxlen = maxlen / (1 << (len < 32 ? len : 31) ); 
			repair[l]->re_set_max_rule_len(maxlen);
		}
	}


	
	void RePairWaveletTree::re_set_alpha_factor(double alpha){
		size_t totalbitseq = 0, totalrepair = 0,totalnonrepair = 0;
		size_t repaircount = 0, nonrepaircount = 0;
		for(uint i = 0; i < height; i++) {
			size_t repairSize = repair[i]->getSize();
			size_t nonrepairSize = nonrepair[i]->getSize();

			if(repairSize < alpha*nonrepairSize){
				bitseq[i] = repair[i];
				totalbitseq += repairSize;
				totalrepair += repairSize;
				repaircount++;
			}
			else{
				bitseq[i] = nonrepair[i];
				totalbitseq += nonrepairSize;
				totalnonrepair += nonrepairSize;
				nonrepaircount++;
			}
		}
		printf("Total size [MB] = %.2f", totalbitseq/1024.0/1024);
		printf(", Total repair [MB] = %.2f(%u)", totalrepair/1024.0/1024, repaircount);
		printf(", Total nonrepair [MB] = %.2f(%u)\n", totalnonrepair/1024.0/1024, nonrepaircount);
	}
	void RePairWaveletTree::re_set_occ_sample(uint sample_rate){				
		if(this->occ)
			delete this->occ;
		BitString oc(this->n);
		for(uint symb = 0; symb <= max_v; ++symb)
			oc.setBit(OCC[symb]-1, 1);
		occ = new BitSequenceRRR(oc,sample_rate);
	}
	void RePairWaveletTree::turn_occ_array(){
		arrayocc = true;
	}
	void RePairWaveletTree::turn_occ_bitmap(){
		arrayocc = false;
	}

};
