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
    //global
    float thd_err;
    ChainScoreMetric chn_score1; ///createTilesFromAnchors2_ ::getGapAnchorsChainScore
    ChainScoreMetric chn_score2; ///chainTiles::getGapBlocksChainScore
    ChainScoreMetric chn_ext_clip_metric1;
    int direction; 
    uint64_t ref_len;
    uint64_t read_len;
    int int_precision;
    int thd_tile_size;
    CharString read_id;

    //extendClipRange() 
    int thd_ecr_shape_len;
    int thd_ecr_reject_da;

    //reform_tiles_ 
    int thd_rfts_cord_gaps;

    //m_g_anchor2_
    uint thd_accept_score;

    //mapExtend_
    int64_t f_me_map_extend;
    int64_t thd_me_reject_gap;

    //createTilesFromChains_
    uint thd_ctfcs_accept_score;
    uint thd_ctfcs_pattern_in_window;

    //g_mapHs_setAnchors_
    int f_gmsa_direction;
    float thd_gmsa_d_anchor_rate;

    //chainTiles
    uint64_t thd_cts_major_limit;

    //createTilesFromAnchors2_
    int64_t thd_ctfas2_connect_danchor;
    int64_t thd_ctfas2_connect_dy_dx;

    //extendsInterval
    int f_eis_raw_clip;     //whether to clip
    int f_eis_raw_clip_ins; //clip two parts of the extension as ins rather than dup
    int thd_eis_shape_len;
    int thd_eis_step1;
    int thd_eis_step2;

    //dropChainGapX
    int thd_dcgx_window_size;
    int thd_dcgx_Xdrop_peak;
    int thd_dcgx_Xdrop_sum;

    //mapTilesFromAnchors
    int thd_mtfas_pattern_in_window;
    int thd_mtfas_overlap_size;
    int thd_mtfas_gap_size;
    int thd_mtfas_overlap_tile;
    int thd_mtfas_swap_tile;
    int thd_mtfas_anchor_density;
    int thd_mtfas_min_segment;

    //trimTiles
    int thd_tts_overlap_size;
    int thd_tts_gap_size;

    //stickMainChain
    int64_t thd_smcn_danchor;

    //getExtendsIntervalChainsOverlaps
    uint64_t thd_dcomx_err_dx;
    uint64_t thd_dcomx_err_dy;

    //extendsIntervalClipOverlaps_
    uint64_t thd_eicos_clip_dxy;
    int thd_eicos_window_size;
    int thd_eicos_f_as_ins;
  
    //mapAlongChain 
    int thd_etfas_shape_len;
    int thd_etfas_step1;
    int thd_etfas_step2;

    //clipChain
    int thd_ccps_clip_min;
    int thd_ccps_clip_init;
    int thd_ccps_clip1_upper;
    int thd_ccps_clip2_lower;
    int thd_ccps_window_size;
    GapParms(float thd_error_rate);
    void clipChainParms(int shape_len, int step1, int step2, float thd_err_rate);
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

int c_stream_(String<Dna5> & seq,String<uint64_t> & g_hs, 
              uint64_t sq_str, uint64_t sq_end, int step, int shape_len, uint64_t type);
int stickMainChain(String<uint64_t> & chain1, String<uint64_t> & chain2, uint64_t(*getX1)(uint64_t), uint64_t(*getY1)(uint64_t), uint64_t(*getX2)(uint64_t), uint64_t(*getY2)(uint64_t));

#endif