#ifndef LINEAR_HEADER_GAP_H
#define LINEAR_HEADER_GAP_H
#include <seqan/sequence.h>
#include "cluster_util.h"
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
int getExtendClipScore(uint64_t const & anchor1, uint64_t const & anchor2, ChainScoreParms & chn_score_parms);

struct GapParms
{
    ChainScoreMetric chn_score1; ///createTilesFromAnchors2_ ::getGapAnchorsChainScore
    ChainScoreMetric chn_score2; ///chainTiles::getGapBlocksChainScore
    ChainScoreMetric chn_ext_clip_metric1;
    int direction; 

    //global parms
    int thd_tile_size;
    CharString read_id;

    //extendClipRange() 
    int thd_ecr_shape_len;
    int thd_ecr_reject_da;

    //reform_tiles_ 
    int thd_rfts_cord_gaps;

    //m_g_anchor2_
    uint thd_accept_score;

    //mapExtend
    int64_t f_me_map_extend;
    int64_t thd_me_reject_gap;

    //createTilesFromChains_
    uint thd_ctfc_accept_score;

    //g_mapHs_setAnchors_
    int f_gmsa_direction;
    float thd_gmsa_d_anchor_rate;

    //chainTiles
    uint64_t thd_cts_major_limit;

    //createTilesFromAnchors2_
    int64_t thd_ctfas2_connect_danchor;
    int64_t thd_ctfas2_connect_dy_dx;


    GapParms(float thd_error_rate);
};

uint64_t extendClip(String<Dna5> & seq1, String<Dna5> & seq2, String<uint64_t> & tiles, uint64_t ext_str, uint64_t ext_end, int direction, GapParms & gap_parms);

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
            int64_t thd_gap, 
            int64_t thd_tileSize,
            float thd_err_rate,
            GapParms & gap_parms
           );

int print_clips_gvf_(StringSet<String<uint64_t> > & clips, 
              StringSet<CharString> & readsId, 
              StringSet<CharString> & genomesId,
              std::ofstream & of);

#endif