#ifndef LINEAR_HEADER_GAP_H
#define LINEAR_HEADER_GAP_H
#include <seqan/sequence.h>
#include "pmpfinder.h"

using namespace seqan;
//ClipRecords : operator struct cord

extern uint64_t EmptyClipConst;

uint64_t getClipStr(String<uint64_t> &, int); //seq1
uint64_t getClipEnd(String<uint64_t> &, int);
int getClipsLen(String<uint64_t> &);
void insertClipStr(String<uint64_t> &, uint64_t);
void insertClipEnd(String<uint64_t> &, uint64_t);
bool isClipEmpty(uint64_t);
int getGapAnchorsChainScore(uint64_t const &, uint64_t const &, ChainScoreParms &);
int getGapBlocksChainScore2(uint64_t const & cord11, uint64_t const & cord12, uint64_t const & cord21, uint64_t const & cord22, uint64_t const & read_len, ChainScoreParms & chn_score_parms);

struct GapParms
{
    ChainScoreMetric chn_score1; ///createTilesFromAnchors2_ ::getGapAnchorsChainScore
    ChainScoreMetric chn_score2; ///chainTiles::getGapBlocksChainScore

    GapParms(float thd_error_rate);
};

/**
 * Re-map gaps in cords.
 * Gaps at the beginning and the end of the read are also included.
 */
int mapGaps(StringSet<String<Dna5> > & seqs, 
            String<Dna5> & read, 
            String<Dna5> & comstr,
            String<uint64_t> & cords_str, 
            String<uint64_t> & cords_end, 
            String<uint64_t> & g_hs,
            String<uint64_t> & g_anchor,
            String<uint64_t> & clips,
            String<UPair> & apx_gaps,
            StringSet<FeaturesDynamic> & f1,
            StringSet<FeaturesDynamic>& f2,
            int64_t thd_gap, 
            int64_t thd_tileSize,
            float thd_err_rate
           );

int print_clips_gvf_(StringSet<String<uint64_t> > & clips, 
              StringSet<CharString> & readsId, 
              StringSet<CharString> & genomesId,
              std::ofstream & of);

#endif