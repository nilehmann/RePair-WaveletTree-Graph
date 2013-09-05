#include "RePairWaveletTreeBuilder.h"

namespace cds_static
{

	RePairWaveletTreeBuilder::RePairWaveletTreeBuilder(BitSequenceBuilder * bsb, Mapper * am,bool deleteSymbols) {
		this->bsb = bsb;
		this->deleteSymbols = deleteSymbols;
		this->am = am;
		bsb->use();
		am->use();
	}

	RePairWaveletTreeBuilder::~RePairWaveletTreeBuilder() {
		bsb->unuse();
		am->unuse();
	}

	Sequence * RePairWaveletTreeBuilder::build(uint * sequence, size_t len) {
		return new RePairWaveletTree(sequence, len, bsb, am, deleteSymbols);
	}

	Sequence * RePairWaveletTreeBuilder::build(const Array & seq) {
		return NULL;
	}
};
