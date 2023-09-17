#ifndef LINEAR_HEADER_GAP_UTIL_H
#define LINEAR_HEADER_GAP_UTIL_H
#include <utility>
#include "cluster_util.h"
#include "pmpfinder.h"

//NOTE::clip & map direction:: towards left < 0, right > 0, both 0
extern uint64_t EmptyClipConst;
int const g_sv_inv = 1;     
int const g_sv_ins = 2;     
int const g_sv_del = 4;     
int const g_sv_trs = 8;     //translocation
int const g_sv_dup = 16;    //duplication
int const g_sv_l = 32;      //clip towards left
int const g_sv_r = 64;      //clip towards right
int const g_sv_gap = 128;   //gap unmapped:inversion duplication translocation.
int const g_sv_shrink = 256; // -->|  |<-- clip norml to the gap  (gap shrink).
int const g_sv_extend = 512; // |<--  -->| clip the gap to normal (gap extend)
int const g_clip_semi_l = 1024; //left
int const g_clip_semi_r = 2048;

int const g_clip_left = -1;
int const g_clip_closed = 0;
int const g_clip_rght = 1;
int const g_map_left = -1;
int const g_map_closed = 0;
int const g_map_rght = 1;

uint64_t getClipStr(String<uint64_t> &, int); //seq1
uint64_t getClipEnd(String<uint64_t> &, int);
void insertClipStr(String<uint64_t> &, uint64_t);
void insertClipEnd(String<uint64_t> &, uint64_t);
bool isClipEmpty(uint64_t);
int getClipsLen(String<uint64_t> &);
int getGapAnchorsChainScore(uint64_t const &, uint64_t const &, ChainScoreParms &);
int getGapBlocksChainScore2(uint64_t const & cord11, uint64_t const & cord12, uint64_t const & cord21, uint64_t const & cord22, uint64_t const & read_len, ChainScoreParms & chn_score_parms);
int getExtendClipScore(uint64_t const & anchor1, uint64_t const & anchor2, ChainScoreParms & chn_score_parms);
//uint64_t extendClip(String<Dna5> & seq1, String<Dna5> & seq2, String<uint64_t> & tiles, uint64_t ext_str, uint64_t ext_end, int direction, GapParms & gap_parms);
void _updateCordsStrEndValue(String<uint64_t> & cords_str,
                             String<uint64_t> & cords_end,
                             unsigned i,
                             uint64_t cord1,
                             uint64_t cord2,
                             uint64_t thd_tile_size);
int64_t _getMaxGapsyOverlap(String<UPair> & gapsy, uint64_t gap_str, uint64_t gap_end);
void remove_tile_sgn (uint64_t & val);
uint64_t get_tile_strand (uint64_t val);
uint64_t is_tile_end(uint64_t val);
void set_tile_end (uint64_t & val);
void set_tile_start (uint64_t & val);
void remove_tile_sgn (uint64_t & val);
void copy_tile_sgn (uint64_t tile1, uint64_t & tile2);
void remove_tile_sgn_start(uint64_t &val);
void remove_tile_sgn_end(uint64_t & val);
bool is_tile_start(uint64_t val);
bool is_tile_body(uint64_t val);
uint64_t shift_tile(uint64_t const & val, int64_t x, int64_t y);
uint64_t get_tile_x (uint64_t val);
uint64_t get_tile_y (uint64_t val);
uint64_t get_tile_id(uint64_t val);
uint64_t create_tile (uint64_t id, uint64_t cordx, uint64_t cordy, uint64_t strand);
void set_tile_strand(uint64_t & val);

int getGapAnchorsChainScore2(uint64_t const & anchor1, uint64_t const & anchor2, ChainScoreParms &
    chn_score_parms);
int getGapBlocksChainScore3(uint64_t const & cord11, uint64_t const & cord12, uint64_t const & cord21,
 uint64_t const & cord22, uint64_t const & read_len, ChainScoreParms & chn_score_parms);
void g_print_tile (uint64_t tile, CharString str);
void g_print_tiles_(String<uint64_t> & tiles, CharString str = "print_tiles");
int insert_tiles2Cords_(String<uint64_t> & cords, unsigned & pos, String<uint64_t> & tiles,
                        int direction, int thd_max_segs_num);
int insert_tiles2Cords_(String<uint64_t> & cords_str, 
                        String<uint64_t> & cords_end,
                        unsigned & pos,
                        String<uint64_t> & tiles_str,
                        String<uint64_t> & tiles_end,
                        int direction,
                        int thd_cord_size,
                        int thd_max_segs_num);

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
    int f_nn1_anchor_sv; //flag of dl::nn1

    //extendClipRange() 
    int thd_ecr_shape_len;
    int thd_ecr_reject_da;

    //reform_tiles_ 
    int thd_rfts_cord_gaps;
    int f_rfts_clip;

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

    //mapGap_
    int64_t thd_mg1_danc_indel;
    int64_t thd_max_extend2;
    int64_t f_dup;

    //mapGaps 
    int64_t thd_gap_len_min;

    GapParms(float thd_error_rate);
    void clipChainParms(int shape_len, float thd_err_rate);
    void printParms(std::string);
};

/*----------  Map and extend components  ----------*/

int mapExtend(StringSet<String<Dna5> > & seqs, 
              String<Dna5> & read, String<Dna5> & comstr,
              StringSet<FeaturesDynamic > & f1, 
              StringSet<FeaturesDynamic > & f2,
              String<uint64_t> & tiles_str, 
              String<uint64_t> & tiles_end, 
              uint64_t gap_str, 
              uint64_t gap_end, 
              int direction,
              GapParms & gap_parms);

int mapExtends(StringSet<String<Dna5> > & seqs, 
               String<Dna5> & read, 
               String<Dna5> & comstr,
               StringSet<FeaturesDynamic > & f1, 
               StringSet<FeaturesDynamic > & f2,
               String<uint64_t> & tiles_str1, 
               String<uint64_t> & tiles_end1,  
               String<uint64_t> & tiles_str2, 
               String<uint64_t> & tiles_end2, 
               uint64_t gap_str1, uint64_t gap_end1, 
               uint64_t gap_str2, uint64_t gap_end2,
               int64_t thd_dxy_min,
               GapParms & gap_parms);
int mapInterval(String<Dna5> & seq1, //genome
                 String<Dna5> & seq2, //read
                 String<Dna5> & comstr,
                 String<uint64_t> & tiles,    //results
                 StringSet<FeaturesDynamic > & f1,  
                 StringSet<FeaturesDynamic > & f2,
                 uint64_t gap_str,
                 uint64_t gap_end, 
                 int64_t anchor_lower,
                 int64_t anchor_upper,
                 int direction,
                 GapParms & gap_parms) // extern parm
;
int mapGeneric(StringSet<String<Dna5> > & seqs, 
               String<Dna5> & read, 
               String<Dna5> & comstr,
               StringSet<FeaturesDynamic > & f1, 
               StringSet<FeaturesDynamic > & f2,
               String<uint64_t> & tiles_str1, 
               String<uint64_t> & tiles_end1,  
               uint64_t gap_str, 
               uint64_t gap_end, 
               GapParms & gap_parms);



#endif