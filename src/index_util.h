#ifndef SEQAN_HEADER_INDEX_UTILITY_H
#define SEQAN_HEADER_INDEX_UTILITY_H
#include <seqan/index.h>
#include "base.h"
#include "shape_extend.h"
#include "index_extend.h"

using namespace seqan;
typedef HIndex<25> LIndex;

bool create_index(StringSet<String<Dna5> > & seq, LIndex & index, unsigned & threads, bool efficient);


//inline uint64_t getXDir(LIndex const & index, uint64_t const & xval, uint64_t const & yval)

#endif