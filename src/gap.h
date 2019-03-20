#ifndef LINEAR_HEADER_GAP_H
#define LINEAR_HEADER_GAP_H
#include <seqan/sequence.h>

using namespace seqan;

/**
 * Re-map gaps in cords.
 * Gaps at the beginning and the end of the read are also included.
 */
int mapGaps(StringSet<String<Dna5> > & seqs, 
            String<Dna5> & read, 
            String<Dna5> & comstr,
            String<uint64_t> & cords, 
            String<uint64_t> & g_hs,
            String<uint64_t> & g_anchor,
            String<uint64_t> & clips,
            StringSet<String<short> > & f1,
            StringSet<String<short> >& f2,
            int const thd_gap, 
            int const thd_tileSize
           );

#endif