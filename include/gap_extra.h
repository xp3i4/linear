#ifndef LINEAR_HEADER_GAP_EXTRA_H
#define LINEAR_HEADER_GAP_EXTRA_H
#include "gap_util.h"

int try_dup(String<Dna5> & seq,
            String<Dna5> & read,
            String<Dna5> & comstr,
            StringSet<FeaturesDynamic > & f1,
            StringSet<FeaturesDynamic > & f2,
            String<uint64_t> & tiles_left, //result
            String<uint64_t> & tiles_rght, //result
            uint64_t gap_str,
            uint64_t gap_end,
            float thd_err_rate,
            uint64_t thd_tile_size,
            GapParms & gap_parms)
            //float band_ratio,
            //int direction)
            ;
#endif