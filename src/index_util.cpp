#include "index_util.h"

using namespace seqan;
bool create_index(StringSet<String<Dna5> > & seq, LIndex & index, unsigned & threads, bool efficient)
{
    return _createQGramIndexDirSA_parallel(seq, index.xstr, index.ysa, index.shape, index.emptyDir, threads, efficient);
}
