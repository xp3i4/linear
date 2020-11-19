#ifndef LINEAR_HEADER_GAP_H
#define LINEAR_HEADER_GAP_H
#include <seqan/sequence.h>
#include "cluster_util.h"
#include "pmpfinder.h"
#include "gap_util.h"

using namespace seqan;
//ClipRecords : operator struct cord



/**
 * Re-map gaps in cords.
 * Gaps at the beginning and the end of the read are also included.
 */
int mapGaps(StringSet<String<Dna5> > & seqs, 
            String<Dna5> & read, 
            String<Dna5> & comstr,
            String<uint64_t> & cords_str, 
            String<uint64_t> & cords_end, 
            String<uint64_t> & clips,
            String<UPair> & apx_gaps,
            StringSet<FeaturesDynamic> & f1,
            StringSet<FeaturesDynamic>& f2,
            GapParms & gap_parms);

int print_clips_gvf_(StringSet<String<uint64_t> > & clips, 
              StringSet<CharString> & readsId, 
              StringSet<CharString> & genomesId,
              std::ofstream & of);

int c_stream_(String<Dna5> & seq,String<uint64_t> & g_hs, 
              uint64_t sq_str, uint64_t sq_end, int step, int shape_len, uint64_t type);
int stickMainChain(String<uint64_t> & chain1, String<uint64_t> & chain2, uint64_t(*getX1)(uint64_t), uint64_t(*getY1)(uint64_t), uint64_t(*getX2)(uint64_t), uint64_t(*getY2)(uint64_t));

#endif