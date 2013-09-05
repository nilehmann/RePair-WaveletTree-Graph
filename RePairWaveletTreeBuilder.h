#include <libcdsBasics.h>

#include <Mapper.h>
#include <SequenceBuilder.h>
#include <BitSequenceBuilder.h>
#include "RePairWaveletTree.h"

#ifndef REPAIR_WAVELET_TREE_BUILDER_H
#define REPAIR_WAVELET_TREE_BUILDER_H

namespace cds_static
{
	class RePairWaveletTreeBuilder : public SequenceBuilder
	{
		public:
			RePairWaveletTreeBuilder(BitSequenceBuilder * bsb, Mapper *am, bool deleteSymbols = false);
			virtual ~RePairWaveletTreeBuilder();
			virtual Sequence * build(uint * seq, size_t len);

			/**
			 * Not implemented
			 * */
			virtual Sequence * build(const Array & seq);

		protected:
			BitSequenceBuilder * bsb;
			Mapper *am;
			bool deleteSymbols;
	};
};
#endif
