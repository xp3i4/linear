#include <utility> 
#include <seqan/align.h>
#include "base.h"
#include "shape_extend.h"
#include "cords.h"
#include "gap_util.h"

/*=============================================
=               Global Variables             =
=============================================*/

int editDist(String<Dna5> & seq1, String<Dna5> seq2, uint64_t str1, uint64_t str2)
{
    //seqan::Score<int16_t, seqan::Simple> scoreAffine(2, -2, -1, -4);
    //typedef Align<String<Dna5>, ArrayGaps> TAlign;
    //TAlign align;
    //resize(rows(align), 2);
    //assignSource(row(align, 0), infix(seq1, str1, str1 + 96));
    //assignSource(row(align, 1), infix(seq2, str2, str2 + 96));
    String<Dna5> tmp1 = infix(seq1, str1, str1 + 96);
    String<Dna5> tmp2 = infix(seq2, str2, str2 + 96);
//    return globalAlignmentScore (tmp1, tmp2, MyersBitVector());
    double t1 = sysTime();
    unsigned score =  0; //globalAlignmentScore (tmp1, tmp2, scoreAffine);
    return score; 
}

GapParms::GapParms(float err_rate) : //estimated err_rate
    //global
    thd_err(err_rate),
    chn_score1(1, 50, &getGapAnchorsChainScore),
    chn_score2(1, 0, &getGapBlocksChainScore2),
    chn_ext_clip_metric1(1, 0, &getExtendClipScore),
    direction(0),
    thd_tile_size(96),
    int_precision(10000), //convert float to int to avoid float division, precision : 10^-4

    thd_ecr_shape_len(3),
    thd_ecr_reject_da(20),
    f_rfts_clip(1),
    f_me_map_extend(0),
    thd_me_reject_gap(200), //~192 = 96 x 2
    thd_accept_score(32),
    thd_ctfcs_accept_score(32),
    thd_ctfcs_pattern_in_window(1),
    thd_gmsa_d_anchor_rate(0.1), //max(ins_rate, del_rate)
    f_gmsa_direction(0),
    thd_cts_major_limit(1),
    thd_ctfas2_connect_danchor(50),
    thd_ctfas2_connect_dy_dx(150),
    f_eis_raw_clip(1),
    f_eis_raw_clip_ins(1),
    thd_eis_shape_len(9),
    thd_eis_step1(5),
    thd_eis_step2(1),
    thd_dcgx_window_size(5),   
    thd_dcgx_Xdrop_peak(125), //kmer size, step of kmer, err_rate related.
    thd_dcgx_Xdrop_sum(60*thd_dcgx_window_size), //kmer size, step, err_rate related.

    thd_tts_overlap_size (thd_tile_size * 0.85),
    thd_tts_gap_size (100),

    thd_smcn_danchor(12), //kmer size and step of kmer, err_rate related.
    thd_dcomx_err_dx(25),
    thd_dcomx_err_dy(25),

    thd_eicos_clip_dxy(30), //given k=5,step1=3,step2=1, 99% distance < 30 (test dataset GIAB)
    thd_eicos_window_size(8), //[weakly] @thd_etfas_shape_len, step1, step2 related.
    thd_eicos_f_as_ins(true),

    thd_etfas_shape_len(5),// gap_parms.thd_ecr_shape_len;
    thd_etfas_step1(3),
    thd_etfas_step2(1),

    //mapGap_
    thd_mg1_danc_indel(80),
    thd_max_extend2(5000)

{ 
    clipChainParms(5, 0.1);
}
//todo::fill in the parms 
void GapParms::clipChainParms(int shape_len, float thd_err_rate)
{
    thd_ccps_window_size = 5;
    //if (thd_err_rate > 0.20)
    //{

    //}
    //else if (thd_err_rate >=0.15 && thd_err_rate < 0.17)
    //{

    //}
    //else //regarded as err = 0.1
    //{
        thd_ccps_clip_min = std::min(thd_err_rate, float(0.1)) * int_precision;
        thd_ccps_clip_init = thd_err_rate * int_precision; 
        thd_ccps_clip1_upper = 8 * int_precision;
        thd_ccps_clip2_lower = 12 * int_precision;
    //}
}
void GapParms::printParms(std::string header)
{
    dout << header 
        << thd_err
        << thd_gap_len_min
        << "\n";
}
/*
 * NOTE! the following parameters are correlated.
 * Change them carefully 
 */
int const g_thd_anchor = 6;
float const g_thd_anchor_density = 0.03;
float const g_thd_error_percent = 0.2;
//----------------------------------------
///c_ functions to clip breakpoints by counting kmers
const unsigned c_shape_len = 8; //anchor shape
const unsigned c_shape_len2 = 4; //base-level clipping shape
const unsigned c_shape_len3 = 4; //base-level clipping gap shape

/**
 * Operations of coordinates of ClipRecords 
 * Based on struct Tile 
 * end[1]|..cord: end = 0 none empty; else empty 
 */
uint64_t EmptyClipConst = (~0) & ((1ULL << 50) - 1);

uint64_t getClipStr(String<uint64_t> & clips, int i) 
{
    if (i << 1 < length(clips)){
        return clips[i << 1];
    }
    else{
        return 0;
    }
}

uint64_t getClipEnd(String<uint64_t> & clips, int i)
{
    return i << 1 < length(clips) - 1 ? clips[(i << 1) + 1] : 0;
}

int getClipsLen(String<uint64_t> & clips)
{
    return length(clips) >> 1;
}

void insertClipStr(String<uint64_t> & clips, uint64_t clip)
{
    if (length(clips) >> 1 << 1 != length(clips))
    {
        appendValue (clips, EmptyClipConst);
        back(clips) &= ~(1023 << 50);
    }
    appendValue(clips, clip);
}

void insertClipEnd(String<uint64_t> & clips, uint64_t clip)
{
    if (length(clips) >> 1 << 1 == length(clips))
    {
        appendValue (clips, EmptyClipConst);
        set_cord_id (back(clips), get_cord_id(clip));
    }
    appendValue (clips, clip);
}

bool isClipEmpty(uint64_t clip) //NOTE: clip is the tile structure.
{
    return (clip & ((1ULL << 50) - 1)) == EmptyClipConst;
}

int isClipTowardsLeft (int clip_direction)
{
    return clip_direction <= 0; 
}

int isClipTowardsRight (int clip_direction)
{
    return clip_direction >= 0;
}

/*
 * String of SV types
 */
int insertSVType(String<int> & sv_types, int type)
{
    appendValue (sv_types, type);
}
int getSVType (String<int> & sv_types, int i)
{
    return sv_types[i];
}

int print_clips_gvf_(StringSet<String<uint64_t> > & clips, 
                     StringSet<CharString> & readsId, 
                     StringSet<CharString> & genomesId,
                     std::ofstream & of)
                     //std::string outputPrefix)
{
    //std::string file_path = outputPrefix + ".gvf";
    //of.open(toCString(file_path));
    of << "##gvf-version 1.10\n";
    std::string source = ".";
    std::string type = ".";
    for (unsigned i = 0; i < length(clips); i++)
    {
        for (unsigned j = 0; j < getClipsLen(clips[i]); j++)
        {
            uint64_t clip_str = getClipStr(clips[i], j);
            uint64_t clip_end = getClipEnd(clips[i], j);
            uint64_t cord_str1 = get_cord_x(clip_str);
            uint64_t cord_str2 = get_cord_y(clip_str);
            uint64_t cord_end1 = get_cord_x(clip_end);
            uint64_t cord_end2 = get_cord_y(clip_end);
            char strand = !(get_cord_strand(clip_str))?'+':'-';
            CharString genomeId = genomesId[get_cord_id(clips[i][j])];
            of  << genomeId << "\t" 
                << source << "\t" 
                << type << "\t"; 
            if (!isClipEmpty(clip_str))
            {
                of << cord_str1 << "\t";   
            }
            else 
            {
                of << ".\t";
            }
            if (!isClipEmpty(clip_end))
            {
                of << cord_end1 << "\t";   
            }
            else 
            {
                of << ".\t";
            }
            of << "readStrand=" << strand << ";";
            of << "readId=" << readsId[i] << ";";
            if (isClipEmpty(clip_str))
            {
                of << "read_clip_str=.;";   
            }
            else 
            {
                of << "read_clip_str=" << cord_str2 <<";";
            }

            if (isClipEmpty(clip_end))
            {
                of << "read_clip_end=.;";   
            }
            else 
            {
                of << "read_clip_end=" << cord_end2 << ";";
            }
            of << "\n";
        }
    }
    return 0;
}

//=======  End of interface function  =======*//

/**
 * Struct Tile  : Cord
 * tile_sign[2]|strand[1]|tileEnd[1](cordEnd)|x[40]|y[20]
 * tile_sign:=1 start, 2 end, 0 body;
 * 0-61 bits same as the Format::Cord
*/
struct TileBase
{
    uint64_t const xBitLen = 40;
    uint64_t const yBitLen = 20;
    uint64_t const strandBit = 61;
    uint64_t const xmask = (1ULL << xBitLen) - 1;
    uint64_t const ymask = (1ULL << yBitLen) - 1;
    uint64_t const sgnBit_str = 1ULL << 62;
    uint64_t const sgnBit_end = 2ULL << 62;
}_defaultTileBase;

struct Tile
{
    uint64_t getX (uint64_t val, 
                   uint64_t const & bit = _defaultTileBase.yBitLen,  
                   uint64_t const & mask = _defaultTileBase.xmask)
    {
        return (val >> bit) & mask;
    }
    uint64_t getY (uint64_t val, 
                   uint64_t const & mask = _defaultTileBase.ymask)
    {
        return val & mask;
    }
    uint64_t makeValue(uint64_t x, uint64_t y, uint64_t const & bit = _defaultTileBase.yBitLen)
    {
        return (x << bit) + y;
    }
    uint16_t getStrand (uint64_t val,
                         uint64_t bit = _defaultTileBase.strandBit)
    {
        return (val >> bit) & 1;
    }
    void setStrand(uint64_t &val, uint64_t bit = _defaultTileBase.strandBit)
    {
        val |= (1ULL << bit);
    }
    void setTileEnd(uint64_t & val, 
        uint64_t bit = _defaultTileBase.sgnBit_end)
    {
        val |= bit;
    }
    void setTileStart(uint64_t & val, uint64_t bit = _defaultTileBase.sgnBit_str)
    {
        val |= bit;
    }
    uint64_t isTileEnd(uint64_t & val, uint64_t bit = _defaultTileBase.sgnBit_end)
    {
        return val & bit;
    }
    bool isTileStart(uint64_t & val, uint64_t bit = _defaultTileBase.sgnBit_str)
    {
        return val & bit;
    }
    bool isTileBody (uint64_t & val)
    {
        return !isTileStart(val) && !isTileEnd(val);
    }
    void removeTileSgnStart(uint64_t & val, uint64_t bit = ~_defaultTileBase.sgnBit_str)
    {
        val &= bit;
    }
    void removeTileSgnEnd(uint64_t & val, uint64_t bit = ~_defaultTileBase.sgnBit_end)
    {
        val &= bit;
    }
    void removeTileSgn(uint64_t & val,
    uint64_t bit = ~(_defaultTileBase.sgnBit_str | _defaultTileBase.sgnBit_end))
    {
        val &= bit;
    }
    void copyTileSgn(uint64_t tile1, uint64_t & tile2, 
        uint64_t bit = _defaultTileBase.sgnBit_str | _defaultTileBase.sgnBit_end)
    {
        tile2 = (tile1 & bit) | (tile2 & (~bit));
    }
}_defaultTile;
/*
 * shortcut to get max overlap(y cord) between the given gap and those in the gaps list 
 * @gapsy: list of gaps(y cord)
 * @str_y, @end_y: start and end of the gap
 */
int64_t _getMaxGapsyOverlap(String<UPair> & gapsy, uint64_t gap_str, uint64_t gap_end)
{
    int64_t overlap = 0;
    int64_t gap_stry = get_cord_y(gap_str);
    int64_t gap_endy = get_cord_y(gap_end);
    for (unsigned i = 0 ; i < length(gapsy); i++)
    {
        int64_t ystr = gapsy[i].first;
        int64_t yend = gapsy[i].second;
        if (gap_stry >= ystr && gap_stry <= yend)
        {
            return std::min(gap_endy, yend) - gap_stry;
        }
        else if (gap_endy >= ystr && gap_endy <= yend)
        {
            return gap_endy - std::max(gap_stry, ystr);
        }
        overlap = 0;
    }
    return  overlap;
}
/*
 * Update value while keep the original sign bit 
 */
void _updateCordsStrEndValue(String<uint64_t> & cords_str,
                             String<uint64_t> & cords_end,
                             unsigned i,
                             uint64_t cord1,
                             uint64_t cord2,
                             uint64_t thd_tile_size)
{
    if (empty(cords_end))
    {
        resize (cords_end, length(cords_str));
        for (unsigned i = 0; i < length(cords_str); i++)
        {
            cords_end[i] = shift_cord (cords_str[i], thd_tile_size, thd_tile_size);
        }
    }
    set_cord_xy (cords_str[i], get_cord_x(cord1), get_cord_y(cord1));
    set_cord_xy (cords_end[i], get_cord_x(cord2), get_cord_y(cord2));
}
inline uint64_t get_tile_strand (uint64_t val)
{
    return _defaultTile.getStrand(val);
}
inline void set_tile_end (uint64_t & val)
{
    _defaultTile.setTileEnd(val);
}
 void set_tile_start (uint64_t & val)
{
    _defaultTile.setTileStart(val);
}
 void remove_tile_sgn (uint64_t & val)
{
    _defaultTile.removeTileSgn(val);
}
void copy_tile_sgn (uint64_t tile1, uint64_t & tile2)
{
    _defaultTile.copyTileSgn(tile1, tile2);
}
void remove_tile_sgn_start(uint64_t &val)
{
    _defaultTile.removeTileSgnStart(val);
}
void remove_tile_sgn_end(uint64_t & val)
{
    _defaultTile.removeTileSgnEnd(val);
}
bool is_tile_start(uint64_t val)
{
    return _defaultTile.isTileStart(val);
}
uint64_t is_tile_end(uint64_t val)
{
    return _defaultTile.isTileEnd(val);
}
bool is_tile_body(uint64_t val)
{
    return _defaultTile.isTileBody(val);
}
uint64_t shift_tile(uint64_t const & val, int64_t x, int64_t y)
{
    return shift_cord (val, x, y);
}
uint64_t get_tile_x (uint64_t val)
{
    return get_cord_x(val);
}
uint64_t get_tile_y (uint64_t val)
{
    return get_cord_y(val);
}
uint64_t get_tile_id(uint64_t val)
{
    return get_cord_id(val);
}
uint64_t create_tile (uint64_t id, uint64_t cordx, uint64_t cordy, uint64_t strand)
{
    return create_cord(id, cordx, cordy, strand);
}
void set_tile_strand(uint64_t & val)
{
    _defaultTile.setStrand(val);
}

void g_print_tile (uint64_t tile, CharString str)
{
    std::cout << str << " " 
              << get_cord_id(tile) << " " 
              << get_cord_strand(tile) << " " 
              << get_cord_x(tile) << " "
              << get_cord_y(tile) << " " 
              << get_cord_x(tile) - get_cord_y (tile) << "\n";    
}
void g_print_tiles_(String<uint64_t> & tiles, CharString str)
{
    for (unsigned i = 0; i < length(tiles); i++)
    {
        std::cout << i << " ";
        g_print_tile (tiles[i], str);
        if (is_tile_end(tiles[i]) || i == length(tiles) - 1)
        {
            std::cout << str << "end\n\n";
        }
    }
}

/*=============================================
=           Index free Map and clip           =
=============================================*/
/**
 * Part 2
 * NOTE: index free mapping for gaps 
 */
/**
 * g_hs_anchor: N/A[13]|strand[1]|anchorX[30]|cord_y[20]
 * @strand := shape strand of kmer in genome ^ shape strand of kmer in read
 * While the kmers are always picked up from the genome and read rather than
 * the reverse complement of the read. 
 * This is different from anchors used in chainning.
 * @anchor := n/a[..]|strand[1]|anchorX[30]
 * @anchorX = x - y + g_hs_anchor_zero. g_hs_anchor_zero to restrict @anchorX > 0. 
   (bits overflow otherwise) such that -g_hs_anchor_zero <= x - y < g_hs_anchor_zero
 */
uint64_t const g_hs_anchor_mask1 = (1ULL << 20) - 1;
uint64_t const g_hs_anchor_mask1_ = ~ g_hs_anchor_mask1;
uint64_t const g_hs_anchor_mask3 = (1ULL << 30) - 1;
uint64_t const g_hs_anchor_mask5 = (1ULL << 31) - 1;
uint64_t const g_hs_anchor_bit1 = 20;
uint64_t const g_hs_anchor_bit2 = 50;
uint64_t const g_hs_anchor_mask2 = ~(1ULL << 50);
uint64_t const g_hs_anchor_zero = 1ULL << (20);

uint64_t g_hs_anchor_getCord (uint64_t anchor)
{
    return anchor & g_hs_anchor_mask1;
}

uint64_t g_hs_anchor_getStrAnchor (uint64_t anchor) // return @anchor := strand + anchorx
{
    return ((anchor >> g_hs_anchor_bit1) & g_hs_anchor_mask5) - g_hs_anchor_zero;
}

uint64_t g_hs_anchor_getX (uint64_t val)
{
    return (((val >> g_hs_anchor_bit1)) & g_hs_anchor_mask3) - g_hs_anchor_zero + 
           (val & g_hs_anchor_mask1);
}

uint64_t g_hs_anchor_getY (uint64_t val)
{
    return val & g_hs_anchor_mask1;
}

uint64_t g_hs_anchor_get_strand(uint64_t val)
{
    return (val >> g_hs_anchor_bit2) & 1;
}

///g_hs: N/A[1]|xval[30]|type[2]|strand[1]|coordinate[30]
///type=0: from genome, type=1: from read
const uint64_t g_hs_bit1 = 30;
const uint64_t g_hs_bit2 = 31;
const uint64_t g_hs_bit3 = 33;
const uint64_t g_hs_mask2 = (1ULL << 30) - 1;
const uint64_t g_hs_mask3 = (1ULL << 32) - 1;

uint64_t g_hs_makeGhs_(uint64_t xval, 
                  uint64_t type, 
                  uint64_t strand, 
                  uint64_t coord)
{
    return (xval << 33) + (type<< 31) + (strand << 30) + coord;
}

int64_t g_hs_getCord(uint64_t & val)
{
    return int64_t(val & g_hs_mask2);
}
//given cord1, return g_hs type anchor = strand + anchor:= strand[1]|anchor[30]
uint64_t g_hs_Cord2StrAnchor(uint64_t cord)
{
    uint64_t strand = get_cord_strand(cord);
    return get_cord_x(cord) - get_cord_y(cord) + (strand << (g_hs_anchor_bit2 - g_hs_anchor_bit1));
}
void g_hs_setAnchor_(uint64_t & val, 
                     uint64_t const & hs1, /*genome*/
                     uint64_t const & hs2, /*read*/
                     uint64_t revscomp_const)
{
    uint64_t strand = ((hs1 ^ hs2) >> 30 ) & 1;
    uint64_t x = revscomp_const * strand - _nStrand(strand) * (hs2 & g_hs_mask2); 
    val = (((hs1 + g_hs_anchor_zero - x) & (g_hs_mask2))<< 20) +  x + (strand << g_hs_anchor_bit2);
}
//create anchors in clip, where strand is ommited
uint64_t c_2Anchor_(uint64_t const & hs1, uint64_t const & hs2)
{
    ///hs1 genome, hs2 read
    uint64_t x = hs2 & g_hs_mask2; 
    return (((hs1 - x + g_hs_anchor_zero) & (g_hs_mask2)) << g_hs_anchor_bit1) + x;
}
///get xvalue and type
uint64_t g_hs_getXT (uint64_t const & val)
{
    return (val >> 31) & g_hs_mask3;
}
uint64_t g_hs_getX (uint64_t const & val)
{
    uint64_t mask = ((1ULL << 30) - 1);
    return ((val >> 33) & mask) ;
}
uint64_t g_hs_anchor2Tile (uint64_t anchor)
{
    uint64_t strand = (anchor >> g_hs_anchor_bit2) & 1;
    /**
     * The cord of read read (y) is shown in its own direction.
     */
    uint64_t y = g_hs_anchor_getY(anchor);
    return (((anchor - (g_hs_anchor_zero << 20) + 
            ((anchor & g_hs_anchor_mask1)<< 20)) & 
              g_hs_anchor_mask2) & g_hs_anchor_mask1_) + y + (strand << 61);
}
int64_t tile_distance_x (uint64_t tile1, uint64_t tile2, uint64_t readlen)
{
    (void)readlen;
    return (int64_t)(get_cord_x(tile2)) - (int64_t)(get_cord_x(tile1));
}
int64_t tile_distance_y (uint64_t tile1, uint64_t tile2, uint64_t readlen)
{
    return  get_tile_strand (tile1 ^ tile2) ? get_tile_y(tile2) - readlen + 1 + get_tile_y(tile1) :
    (int64_t)(_defaultTile.getY(tile2)) - (int64_t)(_defaultTile.getY(tile1));
}

/*
 * shortcut to return if 2 cords have different anchors 
 * @cord1 and @cord2 are required to have same strand
 * @thd_dxy_min lower bound of dx dy
 * @thd_da_zero < is treated as 0.
 * @greater > 0 return true if anchor1 >> anchor2 (significantly larger)
 * @greater < 0 ...  anchor1 << anchor2
 * @greater = 0 ...  |anchor1 - anchor2| >> 0
 */
bool is_diff_anchor (uint64_t cord1, uint64_t cord2, int greater, int64_t thd_dxy_min, float thd_da_zero)
{
    int64_t dy = get_cord_y(cord2) - get_cord_y(cord1);
    int64_t dx = get_cord_x(cord2) - get_cord_x(cord1);
    int64_t dmax = std::max(std::abs(dx), std::abs(dy));
    return std::abs(dy - dx) > int64_t(std::max(thd_dxy_min, dmax) * thd_da_zero) && 
           dy - dx * greater >= 0;
}

/**
 * Shortcut to set main and recd flag for tiles
 * main and recd are sign of Cords.
 * tiles sgn will be cleared and replaced by cords sgn.
 */
void set_tiles_cords_sgns(String<uint64_t> & tiles, uint64_t sgn)
{
    for (int i = 0; i < length(tiles); i++)
    {
        remove_tile_sgn(tiles[i]);
        set_cord_gap(tiles[i]);
        set_cord_recd(tiles[i] , sgn);
    }
}

/**
 * collecting k-mers in 'seq' to 'g_hs'
 */
int g_mapHs_kmer_(String<Dna5> & seq, 
                   String<uint64_t> & g_hs, 
                   uint64_t str, 
                   uint64_t end, 
                   int shape_len,
                   int step,  
                   uint64_t type)
{
    LShape shape(shape_len);
    hashInit(shape, begin(seq) + str);
    int count = 0; 
    int i = 0; 
    uint64_t val = 0;
    for (uint64_t k = str; k < end; k++)
    {
        val = hashNextV(shape, begin(seq) + k);
        if (++count == step)  //collecting every step bases
        {
            //TODO: k - getT(shape)
            appendValue(g_hs, g_hs_makeGhs_(val, type, shape.strand, k));
            count = 0;
        }
    }
    return length(g_hs);
}

/**
 * Stream the block of @g_hs specified by @p1, @p2, @k and 
   keep anchors that with in [@acnhor_lower, anchor_upper).
 * When @direction is towards left, the function only collect anchors that can extend 
   from @gap_end to @gap_str.
 */
int g_mapHs_setAnchors_ (String<uint64_t> & g_hs, 
                         String<uint64_t> & g_anchor,
                         int p1, int p2, int k, 
                         uint64_t revscomp_const, int64_t anchor_lower, int64_t anchor_upper, 
                         uint64_t gap_str, uint64_t gap_end, int direction, GapParms & gap_parms)
{
    unsigned n = 0;
    if (direction == 0)
    {
        for (int i = p1; i < p2; i++) 
        {
            for (int j = p2; j < k; j++) 
            {
                uint64_t tmp_anchor;
                g_hs_setAnchor_(tmp_anchor, g_hs[i], g_hs[j], revscomp_const);
                int64_t tmp = g_hs_anchor_getStrAnchor(tmp_anchor);
                if (tmp < anchor_upper && tmp >= anchor_lower){
                    appendValue(g_anchor, tmp_anchor);
                }
            }   
        }
    }
    else if (direction < 0) //towards left : for mapping extend collect anchor of same strand and close to gap_str or gap_end, NOTE<red> anchors of different strands with gap_str or gap_end are skipped
    {
        int64_t y_end = get_cord_y(gap_end);
        int64_t anchor_base = g_hs_Cord2StrAnchor(gap_end);
        int64_t d_anchor = (1LL << 7) * gap_parms.thd_gmsa_d_anchor_rate;
        int64_t anchor_lower2, anchor_upper2;
        for (int i = p1; i < p2; i++) 
        {
            for (int j = p2; j < k; j++) 
            {
                uint64_t tmp_anchor;
                g_hs_setAnchor_(tmp_anchor, g_hs[i], g_hs[j], revscomp_const);
                int64_t tmp = g_hs_anchor_getStrAnchor(tmp_anchor);
                int64_t dy = y_end - g_hs_anchor_getY(tmp_anchor); 
                if (dy < 0 || (g_hs_anchor_get_strand(tmp_anchor) ^ get_cord_strand(gap_str))){
                    continue;
                }
                else{
                    int64_t d_anchor_acc = std::max((dy >> 7) * d_anchor, int64_t(50)); // d_anchor accumulated
                    anchor_lower2 = std::max(anchor_base - d_anchor_acc, int64_t(0)); // 1<<7 =128, every 128bp increase by d_anchor
                    anchor_upper2 = anchor_base + d_anchor_acc;
                }
                if (tmp < anchor_upper2 && tmp >= anchor_lower2){
                    appendValue(g_anchor, tmp_anchor);
                }
            }   
        }
    }
    else if (direction > 0) 
    {
        int64_t y_str = get_cord_y(gap_str);
        int64_t anchor_base = g_hs_Cord2StrAnchor(gap_str);
        int64_t d_anchor = (1LL << 7) * gap_parms.thd_gmsa_d_anchor_rate;
        int64_t anchor_lower2, anchor_upper2;
        for (int i = p1; i < p2; i++) 
        {
            for (int j = p2; j < k; j++) 
            {
                uint64_t tmp_anchor;
                g_hs_setAnchor_(tmp_anchor, g_hs[i], g_hs[j], revscomp_const);
                int64_t tmp = g_hs_anchor_getStrAnchor(tmp_anchor);
                int64_t dy = g_hs_anchor_getY(tmp_anchor) - y_str; 
                if (dy < 0 || (g_hs_anchor_get_strand(tmp_anchor) ^ get_cord_strand(gap_str))){
                    continue;
                }
                else{
                    int64_t d_anchor_acc = std::max((dy >> 7) * d_anchor, int64_t(50)); // d_anchor accumulated
                    anchor_lower2 = std::max(anchor_base - d_anchor_acc, int64_t(0)); // 1<<7 =128, every 128bp increase d_anchor
                    anchor_upper2 = anchor_base + d_anchor_acc;
                }
                if (tmp < anchor_upper2 && tmp >= anchor_lower2){
                    appendValue(g_anchor, tmp_anchor);
                }
            }   
        }
    }
    return 0;
}

/*----------  Section of MapAnchor2_: function, parm and wrapper  ----------*/
/**
  ::dcgx::simple X-drop by counting gap lens
  Chains in ascending order required : xi < xj && yi < yj provided i < j
*/
int dropChainGapX(String<uint64_t> & chains, uint64_t (*getX)(uint64_t), uint64_t(*getY)(uint64_t), int direction, bool f_erase, GapParms & gap_parms)
{
    if (direction == g_map_rght)
    {
        for (int i = 1; i < length(chains); i++) 
        {
            uint di = i + 1 >= gap_parms.thd_dcgx_window_size ? gap_parms.thd_dcgx_window_size : 1;
            if (getX(chains[i]) - getX(chains[i - 1]) > gap_parms.thd_dcgx_Xdrop_peak || 
                getX(chains[i]) - getX(chains[i + 1 - di]) > gap_parms.thd_dcgx_Xdrop_sum ||
                getY(chains[i]) - getY(chains[i - 1]) > gap_parms.thd_dcgx_Xdrop_peak || 
                getY(chains[i]) - getY(chains[i + 1 - di]) > gap_parms.thd_dcgx_Xdrop_sum)
            {
                if (f_erase)
                {
                    resize (chains, i);
                    return length(chains);
                }
                return i;
            }
        }
        return length(chains);
    }
    if (direction == g_map_left)
    {
        for (int i = length(chains) - 2; i > 0; i--) 
        {
            uint di = length(chains) - i >= gap_parms.thd_dcgx_window_size ? gap_parms.thd_dcgx_window_size : 1;
            if (getX(chains[i + 1]) - getX(chains[i]) > gap_parms.thd_dcgx_Xdrop_peak || 
                getX(chains[i + di - 1]) - getX(chains[i]) > gap_parms.thd_dcgx_Xdrop_sum ||
                getY(chains[i + 1]) - getY(chains[i]) > gap_parms.thd_dcgx_Xdrop_peak || 
                getY(chains[i + di - 1]) - getY(chains[i]) > gap_parms.thd_dcgx_Xdrop_sum)
            {
                if (f_erase)
                {
                    erase(chains, 0, i + 1);
                    return 0;
                }
                return i;
            }
        }
        return 0;
    }
}

unsigned _get_tile_f_ (uint64_t & tile,
                       StringSet<FeaturesDynamic> & f1,
                       StringSet<FeaturesDynamic> & f2)
{
   // uint64_t tile = shift_tile(k_tile, )
    uint thd_abort_score = UMAX;
    uint64_t tile_x = _defaultTile.getX(tile);
    uint64_t tile_y = _defaultTile.getY(tile);
    uint64_t n1 = get_tile_strand(tile);
    uint64_t n2 = get_tile_id(tile);
    unsigned fscore;
    if (n1 < length(f1) && n2 < length(f2))
    {
        fscore = _windowDist(f1[n1],  f2[n2],
                    _DefaultCord.cord2Cell(tile_y), 
                    _DefaultCord.cord2Cell(get_tile_x(tile)));
    }
    else
    {
        fscore = thd_abort_score;
    }
    return fscore;
}
/*
 * only when score < @thd_accept_score, @new_tile is meaningful
 *
unsigned _get_tile_f_tri_ (uint64_t & new_tile,
                           StringSet<FeaturesDynamic > & f1,
                           StringSet<FeaturesDynamic > & f2, 
                           unsigned thd_accept_score,
                           int thd_tile_size)
{
    int shift = thd_tile_size / 4;
    unsigned thd_abort_score = UMAX; // make sure thd_abort_score > thd_accept_score
    double t1 = sysTime();
    unsigned fscore =  _get_tile_f_ (new_tile, f1, f2) ;

    if (fscore <= thd_accept_score)
    {
        return fscore;
    }
    else 
    {
        uint64_t tile_l = shift_tile(new_tile, -shift, -shift);
        fscore = std::min(_get_tile_f_(tile_l, f1, f2), fscore);
        if (fscore <= thd_accept_score)
        {
            new_tile = tile_l;
            return fscore;
        }
        else
        {
            uint64_t tile_r = shift_tile(new_tile, shift, shift);
            new_tile = tile_r;
            fscore = std::min(_get_tile_f_(tile_r, f1, f2), fscore);
            return fscore;
        }
    }
    return fscore;
}
*/
/*
 * only when score < @thd_accept_score, @new_tile is meaningful
 */
unsigned _get_tile_f_tri_ (uint64_t & new_tile,
                           StringSet<FeaturesDynamic > & f1,
                           StringSet<FeaturesDynamic > & f2, 
                           unsigned thd_accept_score,
                           int thd_tile_size)
{
    int shift = thd_tile_size / 4;
    unsigned thd_abort_score = UMAX; // make sure thd_abort_score > thd_accept_score
    unsigned fscore1 =  _get_tile_f_ (new_tile, f1, f2) ;
    unsigned min_score = fscore1;
    uint64_t tile_l = shift_tile(new_tile, -shift, -shift);
    unsigned  fscore2 = _get_tile_f_(tile_l, f1, f2);

    if (fscore2 < fscore1)
    {
        new_tile = tile_l;
        min_score = fscore2;
    }
    uint64_t tile_r = shift_tile(new_tile, shift, shift);
    unsigned fscore3 = _get_tile_f_(tile_r, f1, f2);
    if (fscore3 < min_score)
    {
        new_tile = tile_r;
        min_score = fscore3;
    }
    return min_score;
}

int createTilesFromAnchors1_(String<uint64_t> & anchor, 
                             String<uint64_t> & tiles, 
                             StringSet<FeaturesDynamic> & f1,
                             StringSet<FeaturesDynamic> & f2,
                             uint64_t gap_str, 
                             uint64_t gap_end,
                             int anchor_end, 
                             int const & thd_tile_size,
                             float const & thd_err_rate,
                             int const & thd_pattern_in_window,
                             float const & thd_anchor_density,
                             int64_t const & thd_min_segment,
                             GapParms & gap_parms)
{
    int anchor_len = 0;
    std::sort (begin(anchor), begin(anchor) + anchor_end);
    anchor[anchor_end] = ~0;
    int prek = 0;
    for (int k = 0; k < anchor_end + 1; k++)
    {
        //TODO: handle thd_min_segment, anchor 
        int64_t d = std::abs((int64_t)g_hs_anchor_getY(anchor[k]) - (int64_t)g_hs_anchor_getY(anchor[prek]));
        if (g_hs_anchor_getStrAnchor(anchor[k]) - g_hs_anchor_getStrAnchor(anchor[prek]) > 
            thd_err_rate * std::max(thd_min_segment, d))
        {
            int thd_anchor_accpet = thd_anchor_density * 
            std::abs(int64_t(g_hs_anchor_getY(anchor[k - 1]) - 
                             g_hs_anchor_getY(anchor[prek])));
            thd_anchor_accpet = std::max (thd_anchor_accpet, 2);
            thd_anchor_accpet = std::min (g_thd_anchor, thd_anchor_accpet);
            if (anchor_len > thd_anchor_accpet) 
            {
                std::sort (begin(anchor) + prek, begin(anchor) + k, 
                           [](uint64_t & s1, uint64_t & s2)
                           {return g_hs_anchor_getX(s2) > g_hs_anchor_getX(s1);
                           });
//                g_CreateTilesFromChains_(anchor, tiles, f1, f2, gap_str, prek, k,
//                                &g_hs_anchor_getX, &g_hs_anchor_getY, &g_hs_anchor_get_strand, gap_parms);
            }
            prek = k;
            anchor_len = 0;
        }
        else
        {
            anchor_len++;
        }
    }
    return 0;
}

//ATTENTION::the Adjust @thd_abort_score if the function is changed
int getGapAnchorsChainScore(uint64_t const & anchor1, uint64_t const & anchor2, ChainScoreParms & chn_score_parms)
{
    int64_t dy = g_hs_anchor_getY(anchor1) - g_hs_anchor_getY(anchor2);
    int64_t dx = g_hs_anchor_getX(anchor1) - g_hs_anchor_getX(anchor2);
    if (dy < 0 || g_hs_anchor_get_strand(anchor1 ^ anchor2) || (std::abs(dx) < 8 && dx != dy)) //abort too close dx, such as dx == 0, dy == 100;
    {
        return -10000;
    }

    int64_t thd_min_dy = 50;
    int64_t da = std::abs(int64_t(g_hs_anchor_getStrAnchor(anchor2) - g_hs_anchor_getStrAnchor(anchor1)));
    int64_t derr =  (100 * da) / std::max(dy, thd_min_dy); // 1/100 = 0.01
    int score_derr;
    int score_dy;
    //d_err
    if (derr < 10)
    {
        score_derr = 0;
    }
    else if (derr < 15)
    {
        score_derr = 10 + 2 * derr ;
    }
    else 
    {
        score_derr =  derr * derr / 10 + 40;
    }

    //d_y
    if (dy < 100)
    {
        score_dy = dy / 4;
    }
    else if (dy < 200)
    {
        score_dy = dy / 3 - 9;
    }
    else 
    {
        score_dy = dy - 145;
    }
    return 100 - score_dy - score_derr ;
}
//chain compact anchors whose anchor value are very close
//supposed to use in extend existing anchor that might be called when mapping ins
//For 9mer:step1 = 5:step2 = 1
int getGapAnchorsChainScore2(uint64_t const & anchor1, uint64_t const & anchor2, ChainScoreParms & chn_score_parms)
{
    int64_t dy = g_hs_anchor_getY(anchor1) - g_hs_anchor_getY(anchor2);
    int64_t dx = g_hs_anchor_getX(anchor1) - g_hs_anchor_getX(anchor2);
    if (dy < 0 || g_hs_anchor_get_strand(anchor1 ^ anchor2) 
        || ((std::abs(dx) < 8 || std::abs(dy) < 8)&& dx != dy)) //abort too close dx, such as dx == 0, dy == 100;
    {
        return -10000;
    }

    int64_t thd_min_dy = 50;
    int64_t da = std::abs(int64_t(g_hs_anchor_getStrAnchor(anchor2) - g_hs_anchor_getStrAnchor(anchor1)));
    int64_t derr =  (100 * da) / std::max({dx, dy, thd_min_dy}); // 1/100 = 0.01
    int score_derr;
    int score_dy;
    //d_err
    if (derr < 5)
    {
        score_derr = 4 * derr;
    }
    else if (derr < 10)
    {
        score_derr = 6 * derr - 10;
    }
    else 
    {
        score_derr =  derr * derr - 5 * derr;
    }

    score_dy = dy * (dy + 300) / 300; 
    return 100 - score_dy - score_derr ;
}

//Warn::yellow > sychronize getApxChainScore3 of same logic if necessary when modifiy this function  
//Warn::red dup(dx < thd_min_dx) is not allowed in this score function.
int getGapBlocksChainScore2(uint64_t const & cord11, uint64_t const & cord12, uint64_t const & cord21, uint64_t const & cord22, uint64_t const & read_len, ChainScoreParms & chn_score_parms)
{
    int64_t thd_min_dy = -40;
    int64_t thd_min_dx = -40;
    int64_t dx, dy, da, d_err; 
    //int f_type = getForwardChainDxDy(cord11, cord12, cord21, cord22, read_len, dx, dy);

    int f_type = getChainBlockDxDy(cord11, cord12, cord21, cord22, read_len, chn_score_parms.chn_block_strand, dx, dy);
    
    int64_t thd_max_dy = 500; 
    int64_t thd_max_dx = 15000; //inv at end can be infinity
    int64_t thd_dup_trigger = -50;
    int64_t dx_ = std::abs(dx);
    int64_t dy_ = std::abs(dy);
    da = dx - dy;
    int score = 0;
    //if (dy < thd_min_dy || (f_type == 0 && dy > thd_max_dy) || dx_ > thd_max_dx)
    if (dx < thd_min_dx || dy < thd_min_dy)
    {
        score = INT_MIN;
        //score = INT_MIN;
    }
    else 
    {
        int64_t score_dy = dy_ > 300 ? dy_ / 4 - 25 : dy_ / 6;  
        int64_t score_dx = dx_ > 300 ? dx_ / 4 - 25 : dx_ / 6;  
        if (f_type == 1) //inv
        {
            score = 80 - score_dy; 
        }
        else if (da < -std::max(dx_ / 4, int64_t(50))) //1/4 = *0.25 , maximum sequence error_rate
        {
            if (dx > thd_dup_trigger) //ins
            {
                score = 80 - score_dx; // any large dy is theoretically allowed
            }
            else //dup
            {
                score = 40 - score_dy; // different from ins the dy of dup is suppoesd to be close enough
            }
        }
        else if (da > std::max(dy / 4, int64_t(50))) //del
        {
            score = 80 - score_dy;
        }
        else //normal 
        {
            score = 100 - score_dy;
        }
    }
    return score;
}

//chain blocks that are very close comppatly
//supposed to be used in extending existing anchor for ins/del
int getGapBlocksChainScore3(uint64_t const & cord11, uint64_t const & cord12, uint64_t const & cord21, uint64_t const & cord22, uint64_t const & read_len, ChainScoreParms & chn_score_parms)
{
    int64_t thd_min_dy = 0;
    int64_t thd_min_dx = 0;
    int64_t dx, dy, da, d_err; 
    //int f_type = getForwardChainDxDy(cord11, cord12, cord21, cord22, read_len, dx, dy);
    int f_type = getChainBlockDxDy(cord11, cord12, cord21, cord22, read_len, chn_score_parms.chn_block_strand, dx, dy);
    
    int64_t thd_max_dy = 500; 
    int64_t thd_max_dx = 15000; //inv at end can be infinity
    int64_t thd_dup_trigger = -50;
    int64_t dx_ = std::abs(dx);
    int64_t dy_ = std::abs(dy);
    da = dx - dy;
    int score = 0;
    if (dx < thd_min_dx || dy < thd_min_dy)
    {
        return INT_MIN;
        //score = INT_MIN;
    }
    int64_t score_dy = dy_ > 300 ? dy_ / 4 - 25 : dy_ / 6;  
    int64_t score_dist;  
    int64_t score_da;
    int64_t da_ratio;
    if (f_type == 1) //inv
    {
        score = 20 - score_dy; 
    }
    else 
    {
        da_ratio = 100 * std::abs(da) / std::max({dx_, dy_, int64_t(100)});
        if (da < 15)
        {
            score_da = da_ratio * (da_ratio + 20) / 40;
        }
        else if (da >= 15 && da < 30)
        {
            score_da = da_ratio * (da_ratio + 50) / 45;
        }
        else 
        {
            score_da = da_ratio * (da_ratio + 100) / 45;
        }

        /*
        if (da_ratio < 5)
        {
            score_da = 4 * da_ratio;
        }
        else if (da_ratio < 10)
        {
            score_da = 6 * da_ratio - 10;
        }
        else 
        {
            score_da = 10 * da_ratio - 50;
        }
        */
        int64_t max_dx_dy_ = std::max(dx_, dy_);
        score_dist = max_dx_dy_ * (max_dx_dy_ + 450)/2000;

        score = 100 - score_da - score_dist;
    }
    return score;
}

int chainTiles(String<uint64_t> & tiles, uint64_t read_len, uint64_t thd_gather_block_gap_size, GapParms & gap_parms)
{
    //insert(tiles, 0, 0);
    String<UPair> str_ends;
    String<UPair> str_ends_p;
    String<int> str_ends_p_score;
    gather_blocks_(tiles, str_ends, str_ends_p, 0, length(tiles), read_len, thd_gather_block_gap_size, 0, 0, &is_tile_end, &set_tile_end);
    
    //preFilterChains2(tiles, str_ends_p, &set_tile_end);
    //ChainScoreMetric chn_score(0, &getGapChainScore2);
    chainBlocksCords(tiles, str_ends_p, gap_parms.chn_score2, read_len, 64, gap_parms.thd_cts_major_limit, &remove_tile_sgn_end, &set_tile_end, 0);
    return 0;
}
int g_CreateChainsFromAnchors_(String<uint64_t> & anchors, String<uint64_t> & tiles,
                             uint64_t & gap_str, uint64_t & gap_end, uint64_t read_len, 
                             GapParms & gap_parms)
{
    uint64_t thd_anchor_gap_size = 100; //warn::not the thd_gap_size
    StringSet<String<uint64_t> > anchors_chains;
    String<int> anchors_chains_score;
    uint block_str = 0;
    uint thd_chain_depth = 20;
    uint64_t thd_chain_dx_depth = 80;
    std::sort(begin(anchors), end(anchors), [](uint64_t & a, uint64_t & b){return g_hs_anchor_getX(a) > g_hs_anchor_getX(b);});
    int thd_best_n = 20;
    chainAnchorsBase(anchors, anchors_chains, anchors_chains_score, 0, length(anchors), 
        thd_chain_depth, thd_chain_dx_depth, thd_best_n, gap_parms.chn_score1, &g_hs_anchor_getX);
    resize (tiles, lengthSum(anchors_chains)); 
    int it = 0;
    for (int i = 0; i < length(anchors_chains); i++)
    {
        for (int j = 0; j < length(anchors_chains[i]); j++)
        {
            tiles[it++] = g_hs_anchor2Tile(anchors_chains[i][j]);
        }
        set_tile_end(tiles[it - 1]);
    } 
    chainTiles(tiles, read_len, thd_anchor_gap_size, gap_parms);
    return 0;
}

/**
 *get the chain closed to gap_str (if direction = right) or gap_end (if ..left)
 *@f_erase_tiles, if remove tile in tmp_tiles that have been filtered out.
 */
std::pair<int, int> getClosestExtensionChain_(String<uint64_t> & tmp_tiles, uint64_t gap_str, uint64_t gap_end, bool f_erase_tiles, GapParms & gap_parms)
{
    int pre_i = 0;
    for (int i = 0; i < length(tmp_tiles); i++)
    {
        if (is_tile_end(tmp_tiles[i]))
        {
            int64_t danchor, dx, dy;
            if (gap_parms.direction < 0) 
            {
                dy = get_tile_y(gap_end) - get_tile_y(tmp_tiles[i]); 
                dx = get_tile_x(gap_end) - get_tile_x(tmp_tiles[i]); 
                danchor = dx - dy;
            }
            else if (gap_parms.direction > 0)
            {
                dy = get_tile_y(tmp_tiles[pre_i]) - get_tile_y(gap_str);
                dx = get_tile_x(tmp_tiles[pre_i]) - get_tile_x(gap_str);
                danchor = dx - dy;
            }
            if (std::abs(danchor) < gap_parms.thd_ctfas2_connect_danchor && 
                std::max(std::abs(dy), std::abs(dx)) < gap_parms.thd_ctfas2_connect_dy_dx)
            {
                if (f_erase_tiles == true)
                {
                    erase (tmp_tiles, 0, pre_i);
                    resize(tmp_tiles, i + 1 - pre_i); //remove [i + 1, end)
                    return std::pair<int, int>(0, length(tmp_tiles));
                }
                else
                {
                    return std::pair<int, int>(pre_i, i + 1);
                }
                break;
            } 
            pre_i = i + 1;
        }
    }  
    //failed to find closest chain
    if (f_erase_tiles)
    {
        clear(tmp_tiles);
    }
    return std::pair<int, int>(0, 0);
}
//Create tiles for the block of @chains within [@it_str, it_end).
//Note::chains are supposed to have one tile_end sign (one block) at most
//The new tiles created from [@chains[it_str], @chains[it_end]) are appended to the @tiles
int g_CreateTilesFromChains_ (String<uint64_t> & chains, 
                            String<uint64_t> & tiles, 
                            StringSet<FeaturesDynamic> & f1,
                            StringSet<FeaturesDynamic> & f2,
                            uint64_t gap_str, 
                            int it_str, 
                            int it_end, 
                            uint64_t(*get_x)(uint64_t), //get_x of chains rather than tile
                            uint64_t(*get_y)(uint64_t),
                            uint64_t(*get_strand)(uint64_t),
                            GapParms & gap_parms)
{
    if (it_end - it_str == 0)
    {
        return 0;
    }
    uint thd_fscore = getWindowThreshold(f1); // todo seqeunce error related 
    uint64_t pre_chain = chains[it_str];
    uint64_t pre_tile = 0;
    int64_t tmp_shift = gap_parms.thd_tile_size / 2;
    uint64_t step = gap_parms.thd_tile_size / 3;
    int kcount = 0; //count of kmers in range of each step
    int scan_str = it_str;
    int scan_end = it_str;
    for (int i = it_str; i <= it_end; i++) //i == it_end is out of anchor, this is suit the last one it_end - 1
    {
        if (i == it_end || get_strand(chains[i] ^ pre_chain) || get_x(chains[i]) > get_x(pre_chain) + step || 
                get_y(chains[i]) > get_y(pre_chain) + step)
        {
            if (i == it_end)
            {
                scan_end = it_end;
            }
            for (int j = scan_end - 1; j >= scan_str; j--)
            {
                uint64_t new_tile = create_tile(get_cord_id(gap_str), 
                                                get_x(chains[j]) - tmp_shift, 
                                                get_y(chains[j]) - tmp_shift,
                                                get_strand(chains[j]));
                //unsigned score = _get_tile_f_(new_tile, f1, f2);
                //g_print_tile(new_tile, "gs2");
                unsigned score =  _get_tile_f_tri_(new_tile, f1, f2, 
                    gap_parms.thd_ctfcs_accept_score, gap_parms.thd_tile_size);
                //g_print_tile(new_tile, "gs22");
                if (kcount >= gap_parms.thd_ctfcs_pattern_in_window && score <= 32 && 
                    get_tile_y(new_tile) > get_tile_y(pre_tile))
                {
                    if (empty (tiles) || is_tile_end(back(tiles)))
                    {
                        set_tile_start(new_tile);
                    }
                    appendValue (tiles, new_tile);
                    pre_tile = new_tile;
                    kcount = i - j;
                    pre_chain = chains[j];
                    break;
                }
            }
            scan_str = i;
            scan_end = i + 1;
        }
        else
        {
            scan_end++;
            kcount++;
        }
    }
    if (!empty(tiles))
    {
        set_tile_end(back(tiles)) ;
    }
    return 0;
}
/*
 * Create tiles for the block of @chains within [@it_str, it_end).
 * Note::chains are supposed to have one tile_end sign (one block) at most
 * This funtion requires @chains already clipped at the start and end of the chain,
   Thus the first elment of @tiles_str and last element of @tiles_end == @chains[0] 
   and back(@chains)
 * The function also generates the @tiles_end
 * @chains are required within [@gap_str, @gap_end)
 * @chains are required to be on one strand
 */
int g_CreateTilesFromChains_ (String<uint64_t> & chains, 
                              String<uint64_t> & tiles_str, 
                              String<uint64_t> & tiles_end,
                              StringSet<FeaturesDynamic> & f1,
                              StringSet<FeaturesDynamic> & f2,
                              uint64_t gap_str, 
                              uint64_t gap_end,
                              int it_str, 
                              int it_end, 
                              uint64_t(*get_x)(uint64_t), //get_x of chains rather than tile
                              uint64_t(*get_y)(uint64_t),
                              uint64_t(*get_strand)(uint64_t),
                              GapParms & gap_parms)
{
    (void)gap_str;
    (void)gap_end;
    //<<debug
    for (int i = it_str; i < it_end; i++)
    {
        g_print_tile(chains[i], "gctf3");
    }
    //>>debug
g_print_tile(gap_str, "gctf51");
g_print_tile(gap_end, "gctf52");

    int tiles_str_i = length(tiles_str);
    String<uint64_t> tiles_str_tmp;
    String<uint64_t> tiles_end_tmp;
    //g_print_tiles_(chains, "gctf1");
    g_CreateTilesFromChains_(chains, tiles_str_tmp, f1, f2, gap_str, it_str, it_end,
         get_x, get_y, get_strand, gap_parms);
    //std::cout << "gctf2" << it_str << it_end << length(tiles_str_tmp)<< empty(tiles_str_tmp) << "\n";
    if (empty (tiles_str_tmp))
    {
        std::cout << "gctf5" << empty(tiles_str_tmp) << "\n";
        return 0;
    }
    g_print_tiles_(tiles_str_tmp, "g35");
    int64_t tile_size = gap_parms.thd_tile_size;
    for (unsigned i = 0; i < length(tiles_str_tmp); i++)
    {
        int64_t dx1 = get_x(chains[it_str]) - get_tile_x(tiles_str_tmp[i]);
        int64_t dy1 = get_y(chains[it_str]) - get_tile_y(tiles_str_tmp[i]);
        //dout << "gctf61" << dx1 << dx2 << dy1 << dy2 << "\n";
        if (dx1 <= 0 && dy1 <= 0)
        {
            if (dx1 ==  0 && dy1 == 0)
            {
                break;
            }
            uint64_t new_head_str = chains[it_str];
            remove_tile_sgn(new_head_str);
            dout << "gctf33" << get_tile_x(new_head_str) << get_tile_y(new_head_str) << is_tile_end(tiles_str_tmp[i]) << is_tile_end(new_head_str) << "\n";
            if (i == 0)
            {
                insertValue(tiles_str_tmp, 0, new_head_str);
            }
            else
            {
                tiles_str_tmp[i - 1] = new_head_str;
                erase(tiles_str_tmp, 0, i - 1);
            }
            break;
        }
        if (i == length(tiles_str_tmp) - 1) //if not found such...
        {
            clear(tiles_str_tmp);
            appendValue(tiles_str_tmp, chains[it_str]);
        }
    }
    resize(tiles_end_tmp, length(tiles_str_tmp));
    for (unsigned i = 0; i < length(tiles_str_tmp); i++)
    {
        tiles_end_tmp[i] = shift_tile (tiles_str_tmp[i], tile_size, tile_size);
    }
    g_print_tiles_(tiles_str_tmp, "g41");
    g_print_tiles_(tiles_end_tmp, "g42");
    for (int i = int(length(tiles_end_tmp)) - 1; i >= 0; i--)
    {
        int64_t dx1 = get_x(chains[it_end - 1]) - get_tile_x(tiles_end_tmp[i]);
        int64_t dy1 = get_y(chains[it_end - 1]) - get_tile_y(tiles_end_tmp[i]);

        if (dx1 >= 0 && dy1 >= 0)
        {
            if (dx1 == 0 && dy1 == 0)
            {
                break;
            }
            erase (tiles_str_tmp, i + 1, length(tiles_str_tmp));
            erase (tiles_end_tmp, i + 1, length(tiles_end_tmp));
            uint64_t new_tail_end = chains[it_end - 1];
            uint64_t new_tail_str = shift_tile(new_tail_end, -tile_size, -tile_size);
            if (is_tile_end(tiles_str_tmp[i]))
            {
                remove_tile_sgn(tiles_str_tmp[i]);
                remove_tile_sgn(tiles_end_tmp[i]); 
                set_tile_end(new_tail_str);
                set_tile_end(new_tail_end);
            }
            dout << "gctf6" << is_tile_end(new_tail_str) << dx1 << dy1 << tile_size << i << length(tiles_str_tmp) << "\n";
            g_print_tile(new_tail_str, "gctf71");
            g_print_tile(new_tail_end, "gctf72");
            appendValue(tiles_str_tmp, new_tail_str);
            appendValue(tiles_end_tmp, new_tail_end);
            break;
        }
        if (i == 0) //if not found such.., then erase all except the first one
        {
            erase(tiles_str_tmp, 1, length(tiles_str_tmp));
            erase(tiles_end_tmp, 1, length(tiles_end_tmp));
            tiles_end_tmp[0] = shift_tile(tiles_end_tmp[0], dx1, dy1);
            //set_tile_end(tiles_str_tmp[0]);
            //set_tile_end(tiles_end_tmp[0]);
        }
    }
    append(tiles_str, tiles_str_tmp);
    append(tiles_end, tiles_end_tmp);
    g_print_tiles_(tiles_str, "gctf31");
    g_print_tiles_(tiles_end, "gctf32");
    return 0;
}
/*
int mapClipChains(String<uint64_t> & chain)
{
    for (int i = 0; i < length(tmp_tiles); i++)
    {
        if (is_tile_end(tmp_tiles[i]))
        {
            g_CreateTilesFromChains_(tmp_tiles, tiles, f1, f2, gap_str, pre_i, i + 1, &get_tile_x, &get_tile_y, &get_tile_strand, gap_parms);
            pre_i = i + 1;
        }
        else if (i < length(tmp_tiles) - 1 && 
            get_tile_strand(tmp_tiles[i] ^ tmp_tiles[i + 1]))
        {
            int len = length(tiles);
    ///extendClipInterval(ref, read, comstr, tmp_tiles, direction, gap_parms);
            g_CreateTilesFromChains_(tmp_tiles, tiles, f1, f2, gap_str, pre_i, i + 1, &get_tile_x, &get_tile_y, &get_tile_strand, gap_parms);
            if (len != length(tiles))
            {
                remove_tile_sgn_end(back(tiles));
            }
            pre_i = i + 1;    
        }
    }     
}
*/
int trimTiles(String<uint64_t> & tiles, 
              StringSet<FeaturesDynamic> & f1, StringSet<FeaturesDynamic> & f2,
              uint64_t gap_str,  uint64_t gap_end, uint64_t revscomp_const, int direction,
              GapParms & gap_parms)
{
/**
 * step1.Extend patch
 * extend window if there are gaps between tiles until the 
   coordinates x1 - x2 < window_size or the gap can't be extend any more
 * ATTENTION: This methods takes no account of the relation between y1 and y2.
 */
    int thd_gap_size = gap_parms.thd_tts_gap_size;
    uint64_t thd_tile_size = gap_parms.thd_tile_size;
    uint64_t thd_overlap_size = gap_parms.thd_tts_overlap_size;

    uint64_t cord_str = gap_str;
    int64_t shift_x = std::min(int64_t(get_cord_x(gap_end) - get_cord_x(gap_str)), int64_t(thd_tile_size));
    int64_t shift_y = std::min(int64_t(get_cord_y(gap_end) - get_cord_y(gap_str)), int64_t(thd_tile_size));
    uint64_t cord_end = shift_cord(gap_end, -shift_x, -shift_y);

    for (int i = 0; i < length(tiles); i++)
    {
        if (is_tile_start(tiles[i]) && direction >= 0)
        {
            int new_num = extendPatch(f1, f2, tiles, i, cord_str, tiles[i], revscomp_const, thd_overlap_size, thd_gap_size, gap_parms.thd_accept_score);
            if (new_num)
            {
                set_tile_start(tiles[i]);
                i += new_num;
                remove_tile_sgn_start(tiles[i]);
            }
        }
        if (is_tile_end(tiles[i]) && direction <= 0)
        {
            int new_num = extendPatch(f1, f2, tiles, i + 1, tiles[i], cord_end, revscomp_const, thd_overlap_size, thd_gap_size, gap_parms.thd_accept_score);   
            if (new_num)
            {
                remove_tile_sgn_end(tiles[i]);
                i += new_num;
                set_tile_end(tiles[i]);
            }
        }
        if (i >= 1 && !is_tile_end (tiles[i - 1]) && !is_tile_start(tiles[i]))
        {
            i += extendPatch(f1, f2, tiles, i, tiles[i - 1], tiles[i], revscomp_const, thd_overlap_size, thd_gap_size, gap_parms.thd_accept_score);   
        }
    }
    g_print_tiles_(tiles, "tms11");
    //g_print_tiles_(tiles, "tms12");
    //step2.Remove tiles out of bound.
    int64_t x_str = get_tile_x(gap_str);
    int64_t y_str = get_tile_y(gap_str);
    int64_t x_end = get_cord_x(gap_end);
    int64_t y_end = get_cord_y(gap_end);
    int di = 0;
    for (int i = 0; i < length(tiles); i++)
    {
        uint64_t x_t = get_tile_x(tiles[i]);
        uint64_t y_t = get_tile_strand(tiles[i] ^ gap_str) ?
                       revscomp_const - 1 - get_tile_y(tiles[i]) - thd_tile_size :
                       get_tile_y(tiles[i]);
        if (x_t < x_str || x_t + thd_tile_size > x_end || 
            y_t < y_str || y_t + thd_tile_size > y_end) //out of bound of [gap_str, gap_end)
        {
            if (is_tile_start (tiles[i]) && is_tile_end(tiles[i]))
            {
                //NONE
            }
            else if (is_tile_start(tiles[i]))
            {
                if (i + 1 < length(tiles)) 
                {
                    set_tile_start(tiles[i + 1]);
                }
            }
            else if (is_tile_end(tiles[i]))
            {
                if (i - di - 1 > 0)
                {
                    set_tile_end (tiles[i - di - 1]);
                }
            }
            else{
                //NONE
            }
            di++; 
        }
        else 
        {
            tiles[i - di] = tiles[i];
        }
    }
    if (di)
    {
        resize (tiles, length(tiles) - di);
    }

    return 0;
}

int g_create_anchors_ (String<uint64_t> & g_hs,
                       String<uint64_t> & g_hs_anchor,
                       int shape_len, 
                       int direction,
                       int64_t anchor_lower,
                       int64_t anchor_upper,
                       uint64_t rvcp_const,
                       uint64_t gap_str,
                       uint64_t gap_end,
                       GapParms & gap_parms)
{
    uint64_t mask = (1ULL << (2 * shape_len + g_hs_bit3)) - 1;
    std::sort (begin(g_hs), end(g_hs), [mask](uint64_t & a, uint64_t & b){return (a & mask) < (b & mask);});
    int p1 = 0, p2 = 0;
    for (int k = 1; k < length(g_hs); k++)
    {    
        switch (g_hs_getXT((g_hs[k] ^ g_hs[k - 1]) & mask))
        {
            case 0:       //x1 = x2 both from genome or read
                break;
            case 1:       //x1 = x2 one from genome the other from read
                p2 = k;
                break;
            default:      //anchor current block before process next block 
                g_mapHs_setAnchors_(g_hs, g_hs_anchor, p1, p2, k, rvcp_const, anchor_lower, anchor_upper, gap_str, gap_end, direction, gap_parms);
                p1 = k;
                p2 = k; 
        }
    }
    return 0;
}

int g_CreateExtendAnchorsPair_ (String<uint64_t> & g_hs,
                       String<uint64_t> & g_hs_anchor1,
                       String<uint64_t> & g_hs_anchor2,
                       int shape_len, 
                       uint64_t rvcp_const,
                       uint64_t gap_str1,
                       uint64_t gap_end1,
                       uint64_t gap_str2,
                       uint64_t gap_end2,
                       GapParms & gap_parms)
{
    uint64_t mask = (1ULL << (2 * shape_len + g_hs_bit3)) - 1;
    std::sort (begin(g_hs), end(g_hs), [mask](uint64_t & a, uint64_t & b){return (a & mask) < (b & mask);});
    int p1 = 0, p2 = 0;
    int direction1 = 1;
    int direction2 = -1;
    for (int k = 1; k < length(g_hs); k++)
    {    
        switch (g_hs_getXT((g_hs[k] ^ g_hs[k - 1]) & mask))
        {
            case 0:       //x1 = x2 both from genome or read
                break;
            case 1:       //x1 = x2 one from genome the other from read
                p2 = k;
                break;
            default:      //anchor current block before process next block 
                g_mapHs_setAnchors_(g_hs, g_hs_anchor1, p1, p2, k, rvcp_const, 0, 0, gap_str1, gap_end1, direction1, gap_parms);
                g_mapHs_setAnchors_(g_hs, g_hs_anchor2, p1, p2, k, rvcp_const, 0, 0, gap_str2, gap_end2, direction2, gap_parms);
                p1 = k;
                p2 = k; 
        }
    }
    return 0;
}

int g_stream_(String<Dna5> & seq1, //genome
              String<Dna5> & seq2, //read
              String<uint64_t> & g_hs,
              uint64_t gap_str,
              uint64_t gap_end, 
              unsigned shape_len,
              int step1,
              int step2,
              GapParms & gap_parms)
{
    //clear(g_hs);
    //resize(g_hs, 1ULL << 20);
    uint64_t gs_str = get_cord_x(gap_str);
    uint64_t gs_end = get_cord_x(gap_end);
    uint64_t gr_str = get_cord_y(gap_str);
    uint64_t gr_end = get_cord_y(gap_end);
    if (get_cord_strand(gap_str))
    {
        gr_str = length(seq2) - gr_str - 1;
        gr_end = length(seq2) - gr_end - 1;
        std::swap (gr_end, gr_str);
    }
    g_mapHs_kmer_(seq1, g_hs, gs_str, gs_end, shape_len, step1, 0);
    g_mapHs_kmer_(seq2, g_hs, gr_str, gr_end, shape_len, step2, 1);    
    return 0;
}
/*----------  Clip function  ----------*/
/**
 * stream seq creating hs
 */
int c_stream_(String<Dna5> & seq,String<uint64_t> & g_hs, 
              uint64_t sq_str, uint64_t sq_end, int step, int shape_len, uint64_t type)
{
    LShape shape(shape_len);
    hashInit_hs(shape, begin(seq) + sq_str, 0);
    int count = 0; 
    int i = 0; 
    uint64_t val = 0;

    for (uint64_t k = sq_str; k < sq_end; k++)
    {
        val = hashNext_hs(shape, begin(seq) + k);
        if (++count == step)  //collecting every step bases
        {
            //TODO: k - getT(shape)
            appendValue(g_hs, g_hs_makeGhs_(val, type, 0, k));
            count = 0;
        }
    }
    return length(g_hs);
}

//using de brujin sequence to calculate the clz and ctz of 4-mers
int const clzb_4_index_[8] = {0, 0, 3, 1, 3, 2, 2, 1}; // de brujin sequence table / 2
 int clzb_4__ (uint64_t a)
{
    uint64_t tmp = a & ((~a) + 1);
    return clzb_4_index_[(tmp - (tmp >> 4) - (tmp >> 5)) & 255]; 
}

short clzb_4_(uint64_t a)
{
    return (a)?__builtin_clz(unsigned (a)) / 2 - 12:4;
}

short ctzb_4_(uint64_t a)
{
    return (a)?__builtin_ctz(unsigned(a)) / 2:4;

}

/**
 * Stream the block of 'g_hs' within [p1,p2)x[p2,k), and convert the    
   production of cords to anchors with the restrictions of
   |candidates_anchor - 'anchor' | < band
 */
int c_createAnchorsBlocks_ (String<uint64_t> & g_hs, 
                             String<uint64_t> & g_anchor,
                             int p1, 
                             int p2, 
                             int k, 
                             int thd_band_level,  //dx >> band_level 
                             int thd_band_lower,  //band lower bound
                             int64_t anchor_x,
                             int64_t anchor_y,
                             int64_t x_lower = 0, //lower bound 
                             int64_t x_upper = 0) 
{
    int64_t dx_lower, dx_upper;
    if (x_lower == 0 && x_upper == 0)
    {
        dx_lower = ~0;
        dx_upper = (1LL << 63) - 1;
    }
    else 
    {
        dx_lower = x_lower - anchor_x;
        dx_upper = x_upper - anchor_x;
    }
    for (int i = p1; i < p2; i++) 
    {
        int dx = g_hs_getCord(g_hs[i]) - anchor_x;
        for (int j = p2; j < k; j++) 
        {
            int dy = g_hs_getCord(g_hs[j]) - anchor_y;
            int d_anchor = std::abs(dx - dy);
            if (d_anchor <= std::max(std::abs(dx) >> thd_band_level, thd_band_lower) 
                && dx < dx_upper && dx > dx_lower)
            {
                appendValue(g_anchor, c_2Anchor_(g_hs[i], g_hs[j]));
            }
        }   
    }
    return length(g_anchor);
}

int c_createAnchors (String<uint64_t> & g_hs, 
                     String<uint64_t> & g_anchors,
                     int g_hs_end,
                     int band_level,
                     int band_lower,
                     int64_t anchor_x,
                     int64_t anchor_y,
                     int64_t x_lower = 0,
                     int64_t x_upper = 0) 
{
    int p1 = 0, p2 = 0;
    std::sort (begin(g_hs), end(g_hs));
    for (int k = 1; k < g_hs_end; k++)
    {
        switch (g_hs_getXT(g_hs[k] ^ g_hs[k - 1]))
        {
            case 0:
                break;
            case 1:
                p2 = k;
                break;
            default:
                 c_createAnchorsBlocks_(
                    g_hs, g_anchors, 
                    p1, p2, k, 
                    band_level, band_lower, 
                    anchor_x, anchor_y, 
                    x_lower, x_upper);
                p1 = k;
                p2 = k; 
        }
    }
    return length(g_anchors);
}
//Create anchors within the given range
int c_createAnchors2 (String<uint64_t> & g_hs, 
                      String<uint64_t> & g_anchors,
                      int g_hs_end,
                      int64_t anchor_lower,
                      int64_t anchor_upper) 
{
    int p1 = 0, p2 = 0;
    std::sort (begin(g_hs), end(g_hs));
    for (int k = 1; k < g_hs_end; k++)
    {
        switch (g_hs_getXT(g_hs[k] ^ g_hs[k - 1]))
        {
            case 0:
                break;
            case 1:
                p2 = k;
                break;
            default:
                for (int i = p1; i < p2; i++) 
                {
                    int64_t x = g_hs_getCord(g_hs[i]);
                    for (int j = p2; j < k; j++) 
                    {
                        int64_t y = g_hs_getCord(g_hs[j]);
                        if (anchor_lower <= x - y && x - y < anchor_upper)
                        {
                            appendValue(g_anchors, c_2Anchor_(g_hs[i], g_hs[j]));
                        }
                    }   
                }
                p1 = k;
                p2 = k; 
        }
    }
    return length(g_anchors);
}

/**
 * []::f9
 * clip by anchors
 * @val1 length of match 
 * @val2 length of cluster
 * !!todo::tune thd_exp_err thd_min_len
 * ---------mmmmmmmmmm
 */
inline int64_t c_sc_(int val1, 
                     int val2, 
                     float thd_exp_err = 0.85)
{
    float rate = (float)val1 / val2;
    float thd_err1 = thd_exp_err;
    float thd_err2 = thd_exp_err; 
    int thd_min_len1 = 10; //20
    int thd_min_len2 = 10;
    if (val1 > 25)
    {
        thd_err1 -= 0.05;
        thd_err2 -= 0.15;
    }  
    else if (val1 > 15)
    {
        thd_err1 += 0.05;
        thd_err2 -= 0.05;
    }
    else
    {
        thd_err1 += 0.1;
        thd_err2 += 0.05;
    }
    if (rate > thd_err1 && val1 > thd_min_len1)
    {
        return val1 << 2; 
    }
    else if (rate > thd_err2 && val1 > thd_min_len2)
    {
        return val1 << 1;
    }
    else 
    {
        return val1;
    }
}

//small anchors


/**
 * Clip the anchors at the end of the leftmost (-1) or rightmost (1) anchor that is well extended. 
 * @clip_direction: -1 gap-match; 1 match-gap
 */
int64_t c_clip_anchors_ (String<uint64_t> & anchor, 
                         uint64_t clip_str,
                         uint64_t clip_end,
                         int shape_len, 
                         int thd_merge1, // thd of anchor
                         int thd_merge1_lower,
                         int thd_merge2, //thd of x
                         int clip_direction,
                         int thd_clip_sc = c_sc_(25, 30),
                         int thd_accept_score = c_sc_(c_shape_len + 3, (c_shape_len + 3) * 2)
                        )
{
    uint64_t gs_str = get_tile_x(clip_str);
    uint64_t gr_str = get_tile_y(clip_str);
    uint64_t gs_end = get_tile_x(clip_end);
    uint64_t gr_end = get_tile_y(clip_end);
    uint64_t genomeId = get_tile_id (clip_str);
    uint64_t gr_strand = get_tile_strand(clip_str);
    int direction = (clip_direction < 0) ? -1 : 1;
    int it = 0;
    int bit1 = 20;
    int bit2 = g_hs_anchor_bit1 + bit1;
    uint64_t mask = (1LL << bit1) - 1;
    uint64_t ct_conts = 0;
    appendValue(anchor, ~0);
    if (length(anchor) < 1)
    {
        return direction > 0 ? clip_str : clip_end;
    }
    std::sort (begin(anchor), end(anchor));
    int i_str = 0;
    for (int i = 0; i < length(anchor) - 1; i++)
    {
        if (g_hs_anchor_getY(anchor[i + 1] - anchor[i]) == 1 &&
            g_hs_anchor_getStrAnchor(anchor[i + 1]) - g_hs_anchor_getStrAnchor(anchor[i])== 0)
        {
            ct_conts++;
        }
        else
        {
            i_str = i - ct_conts ;
            uint64_t x = g_hs_anchor_getX(anchor[i_str]) - gs_str;
            uint64_t y = g_hs_anchor_getY(anchor[i_str]) - gr_str;
            anchor[it++] = (x << bit2) + ((ct_conts + 1) << g_hs_anchor_bit1) + y;
            ct_conts = 0;
            //collect continuos patterns (no gaps)
        }
    }
    i_str = length(anchor) - 1 - ct_conts;
    uint64_t x = g_hs_anchor_getX(anchor[i_str]) - gs_str;
    uint64_t y = g_hs_anchor_getY(anchor[i_str]) - gr_str;
    anchor[it++] = (x << bit2) + ((ct_conts + 1) << g_hs_anchor_bit1) + y;

    if (it < 1)
    {
        return direction > 0 ? clip_str : clip_end;
    }
    //!NOTE::Value of anchor has been changed to := x|ct_conts|y
    int64_t y1 = 0, y2 = 0;
    int64_t x1 = 0, x2 = 0;
    int64_t x1_end = 0, x2_end = 0;
    if (direction < 0)
    {
        int max_score = 0;
        uint64_t max_anchor = 0;
        std::sort(begin(anchor), begin(anchor) + it);
        for (int i = 0; i < it; i++) //extend anchor[i]
        {
            y1 = g_hs_anchor_getY(anchor[i]);
            x1 = (anchor[i] >> bit2) & mask;
            x1_end = x1 + ((anchor[i] >> bit1) & mask) + shape_len - 1;
            int score = c_sc_(x1_end - x1, x1_end - x1);
            int dj = 0;
            for (int j = i + 1; j < it; j++) 
            {
                //#anchor will be shrinked(overwrite anchor[j]) 
                //if anchor[j] can be merged to the
                //the chain starting from the anchor[i].
                y2 = g_hs_anchor_getY(anchor[j]);
                x2 = (anchor[j] >> bit2) & mask;
                x2_end = x2 + ((anchor[j] >> bit1) & mask) + shape_len - 1;
                int64_t da = x2 - x1 - y2 + y1;
                int thd_da_accept = std::max(int(x2 - x1) >> thd_merge1, 
                                             thd_merge1_lower);
                if (std::abs(da) < thd_da_accept && 
                    x2 - x1 < thd_merge2 && 
                    x1 < x2)
                {
                    score += c_sc_(x2_end - x2, x2_end - x1_end);
                    y1 = y2;
                    x1 = x2; 
                    x1_end = x2_end;
                    ++dj;
                }
                else
                {
                    anchor[j - dj] = anchor[j];
                }
            }
            if (score > thd_clip_sc)
            {
                int64_t rslt_x = (anchor[i] >> bit2) & mask; 
                int64_t rslt_y = (anchor[i] & mask);
                uint64_t clip = create_cord(genomeId, gs_str + rslt_x, gr_str + rslt_y, gr_strand);
                return clip;
            }
            else if (score > max_score && score > thd_accept_score)
            {
                max_score = score;
                max_anchor = anchor[i];
            }
            it -= dj;
        }
        if (max_score > 0)
        {
            int64_t rslt_x = (max_anchor >> bit2) & mask;
            int64_t rslt_y = (max_anchor & mask);
            uint64_t clip = create_cord(genomeId, gs_str + rslt_x, gr_str + rslt_y, gr_strand);
            return clip;
        } 
        else
        {
            return clip_end;
        }
    }
    else if (direction > 0)
    {
        int max_score = 0;
        uint64_t max_anchor = 0;
        std::sort(begin(anchor), begin(anchor) + it, std::greater<uint64_t>());
        for (int i = 0; i < it; ++i)
        {
            y1 = g_hs_anchor_getY(anchor[i]);
            x1 = (anchor[i] >> bit2) & mask;
            x1_end = x1 + ((anchor[i] >> bit1) & mask) + shape_len - 1;
            int score = c_sc_(x1_end - x1, x1_end - x1);
            int dj = 0;
            for (int j = i; j < it; ++j)
            {
                y2 = g_hs_anchor_getY(anchor[j]);
                x2 = (anchor[j] >> bit2) & mask;
                x2_end = x2 + ((anchor[j] >> bit1) & mask) + shape_len - 1;
                int64_t da = x2 - x1 - y2 + y1;
                int thd_da_accept = std::max(int(x1 - x2) >> thd_merge1, 
                                             thd_merge1_lower);
                if (std::abs(da) < thd_da_accept && 
                    x1 - x2_end < thd_merge2 && 
                    x1 > x2)
                {
                    score += c_sc_(x2_end - x2, x2_end - x1_end);
                    x1 = x2; 
                    y1 = y2;
                    x1_end = x2_end;
                    ++dj;
                }
                else
                {
                    anchor[j - dj] = anchor[j];
                }
            }
            if (score > thd_clip_sc)
            {
                int64_t rslt_x = (anchor[i] >> bit2) & mask;
                int64_t rslt_y = (anchor[i] & mask);
                uint64_t clip = create_cord(genomeId, gs_str + rslt_x, gr_str + rslt_y, gr_strand);
                return clip;
            }
            else if (score > max_score && score > thd_accept_score)
            {
                max_score = score;
                max_anchor = anchor[i];
            }
            it -= dj;
        } 
        if (max_score > 0)
        {
            int64_t rslt_x = (max_anchor >> bit2) & mask;
            int64_t rslt_y = (max_anchor & mask);
            uint64_t clip = create_cord(genomeId, gs_str + rslt_x, gr_str + rslt_y, gr_strand);
            return clip;

        }
        else
        {
            return clip_str;
        }
    }
    else
    {
        return clip_str;
    }
}

/**
* kmer of t1 is left to t2
*/
int c_isGapMatch_(uint64_t & dv, short& t1, short & t2, short & l1, short & l2, short k)
{
   /*
   if (dv == 0) {
       return 1; // match
   }
   if (t1 + l1 - k + 1 == 0){
       return 2; // mismatch
   }
   if (t2 + l1 - k == 0 && t2 != k && l1 != k){
       return 3; //del 
   }
   if (t1 + l2 - k + 1 == 0){
       return 4;  //ins
   }
   return 0;
   */
   return ((dv == 0) || (t1 + l1 - k + 1 == 0) || 
    (t2 + l1 - k == 0 && t2 != k && l1 != k) || (t1 + l2 - k + 1 == 0)) ? 1 : 0;
}

/***********************<Section: extend clip*************************/
//chain compact and small anchors 
//supposed to use in extendClip 5mer:step1 = 5:step2=1
int getExtendClipScore(uint64_t const & anchor1, uint64_t const & anchor2, ChainScoreParms & chn_score_parms)
{
    int64_t dy = g_hs_anchor_getY(anchor1) - g_hs_anchor_getY(anchor2);
    int64_t dx = g_hs_anchor_getX(anchor1) - g_hs_anchor_getX(anchor2);
    if (dy <= 0 || g_hs_anchor_get_strand(anchor1 ^ anchor2) 
        || ((std::abs(dx) < 3 || std::abs(dy) < 3) && dx != dy)) //abort too close dx, such as dx == 0, dy == 100;
    {
        return -10000;
    }

    int64_t thd_min_dy = 10;
    int64_t da = std::abs(int64_t(g_hs_anchor_getStrAnchor(anchor2) - g_hs_anchor_getStrAnchor(anchor1)));
    int64_t derr =  (100 * da) / std::max({dx, dy, thd_min_dy}); // 1/100 = 0.01
    int score_da;
    int score_dy;
    //d_err
    if (da == 0)
    {
        score_dy == 0;
    }
    if (da < 2)
    {
        score_da = 30 + 5 * da;
    }
    else if (da < 5)
    {
        score_da = 36 +  2 * da;
    }
    else
    {
        score_da = 41 + da;
    }
    score_dy = dy * (12 * dy + 650) / 450; 

    return 100 - score_dy - score_da ;
}
/*
 * Simple accumlated score of counting matches, taking less computational complexity
 * @shape_len : length of shape to create the @chain
 * @_getX : pass getY if accumulate y 
 * Always accumulate from left to right, clip direction is not considerd in the function(so do not use too large kmers to avoid introduced by length of kmers).
 */
int accumulateSimpleGapScore1(String<uint64_t> & chain, String<int> & gaps_score, int shape_len, uint64_t(*_getX)(uint64_t), GapParms & gap_parms)
{
    if (empty(chain))
    {
        return -1;
    }
    resize(gaps_score, length(chain), 0);
    uint64_t pre_x =  _getX(chain[0]);
    for (int i = 1; i < length(chain); i++)
    {
        uint64_t x_i = _getX(chain[i]);
        int new_gap = int(x_i - pre_x) > shape_len ? x_i - pre_x - shape_len : 0;
        gaps_score[i] += gaps_score[i - 1] + new_gap * gap_parms.int_precision;
        pre_x = x_i; 
    }    
    return 0;
}
/*
 * Find and Clip the chain at the breakpoint;
 * Method: ds(b) = max{s(b - dx) - s(b + dx)}
 * Namely clip at the point where the difference of score of two windows at the two sides of the point reaches the maximum.
 */
int clipChain_(String<uint64_t> & chain, 
               String<int> & gaps_score_x, 
               String<int> & gaps_score_y, 
               int direction, 
               bool f_clip,  
               uint64_t (*_get_x)(uint64_t), 
               uint64_t (*_get_y)(uint64_t), 
               GapParms & gap_parms)
{
    if (empty(chain))
    {
        return -1;
    }
    g_print_tiles_(chain, "cc1");
    int clip_i = isClipTowardsLeft(direction) ? - 1 : length(chain) - 1;
    int clip_x1, clip_x2, clip_y1, clip_y2;
    int thd_window_i_size = gap_parms.thd_ccps_window_size; //window_bps >= 
    //step1 * (thd_window_i_size - 1) + shape_len (when no gaps, equal)
    int max_d_clip = INT_MIN;
    int f_found_clip = 0;
    for (int i = 1; i < length(chain) - 1; i++)
    {
        int i_str = std::max(i - thd_window_i_size, 0);
        int i_end = std::min(i + thd_window_i_size, int(length(chain) - 1));
        int d1 = i - i_str;
        int d2 = i_end - i;
        clip_x1 = (gaps_score_x[i] - gaps_score_x[i_str]) / d1;
        clip_x2 = (gaps_score_x[i_end] - gaps_score_x[i]) / d2;
        clip_y1 = (gaps_score_y[i] - gaps_score_y[i_str]) / d1;
        clip_y2 = (gaps_score_y[i_end] - gaps_score_y[i]) / d2;

        if (isClipTowardsLeft(direction))
        {
            std::swap (clip_x1, clip_x2);
            std::swap (clip_y1, clip_y2);
        }
        int d_clip = clip_x2 - clip_x1 + clip_y2 - clip_y1;
        dout << "cc2" << get_cord_y(chain[i]) << d_clip << clip_x2 << clip_x1  << clip_y1 << clip_y2 << d1 << d2 << gap_parms.thd_ccps_clip1_upper << gap_parms.thd_ccps_clip2_lower << gaps_score_x[i] << gaps_score_x[i] << "\n";
        if (d_clip > max_d_clip && 
            clip_x1 < gap_parms.thd_ccps_clip1_upper && clip_y1 < gap_parms.thd_ccps_clip1_upper &&
            (clip_x2 > gap_parms.thd_ccps_clip2_lower || clip_y2 > gap_parms.thd_ccps_clip2_lower))
        {
            max_d_clip = d_clip;
            clip_i = i;
            dout << "cc3" << get_tile_y(chain[i]) << max_d_clip << "\n";
            f_found_clip = 1;
        }
    }
    if (f_clip && f_found_clip)
    {
        if (isClipTowardsLeft(direction))
        {
            erase(chain, 0, clip_i + 1);
        }
        else 
        {
            resize(chain, clip_i + 1);
        }
    }
    g_print_tiles_(chain, "cc4");
    return clip_i + 1;
}

int clipChain(String<uint64_t> & chain, int shape_len, int direction, bool f_clip,  uint64_t (*_get_x)(uint64_t), uint64_t (*_get_y)(uint64_t), GapParms & gap_parms)
{
    gap_parms.clipChainParms(shape_len, gap_parms.thd_err); //init clip parms

    String <int> gaps_score_x;
    String <int> gaps_score_y; 

    accumulateSimpleGapScore1(chain, gaps_score_x, shape_len, _get_x, gap_parms);
    accumulateSimpleGapScore1(chain, gaps_score_y, shape_len, _get_y, gap_parms);
    return clipChain_(chain, gaps_score_x, gaps_score_y, direction, f_clip, _get_x, _get_y, 
                    gap_parms);
}

/*
 * Generic function
   to filter records in @chain1 that located around the records of @chain2.
 * chain1 are required to be sorted by x in DESCENDING order already, 
   It's supposed to be the reversely sorted anchors before chaning
 * chain2 are required to be sorted by x in AESCENDING order already, 
   It's supposed to be the sorted chain after chaining.
 * @chain1 is the chain to be filtered, @chain2 is the main chain
 * Record in @chain1 r1 is sticked to the largest record in @chain2 r2 and r1 > r2,
    namely r1_x > r2_x;
 */
int stickMainChain(String<uint64_t> & chain1, 
                   String<uint64_t> & chain2, 
                   uint64_t(*getX1)(uint64_t), 
                   uint64_t(*getY1)(uint64_t), 
                   uint64_t(*getX2)(uint64_t), 
                   uint64_t(*getY2)(uint64_t), 
                   GapParms & gap_parms)
{
    if (empty(chain1) || empty(chain2))
    {
        return 0;
    }
    int di = 0, jj = length(chain2) - 1;
    uint64_t x1, x2 = getX2(chain2[jj]);
    for (int i = 0; i < length(chain1); i++)
    {
        x1 = getX1(chain1[i]);
        if (x1 < x2)
        {
            for (int j = jj - 1; j >= 0; j--)    
            {
                x2 = getX2(chain2[j]);
                if (x1 >= x2)
                {
                    jj = j;
                    break;
                }
            }
        }
        if (x1 < x2)//none such x2 exists in chain2 : x1 <= all in chains1
        {
            jj = 0; //replaced with the first element of chain2 that is the closet to x1.
        }
        int64_t anchor1 = x1 - getY1(chain1[i]); 
        int64_t anchor2 = getX2(chain2[jj]) - getY2(chain2[jj]);
        if (anchor1 >= anchor2 + gap_parms.thd_smcn_danchor ||  
            anchor1 < anchor2 - gap_parms.thd_smcn_danchor)
        {
            di++;
        }
        else
        {
            chain1[i - di] = chain1[i];
        }
    }
    resize (chain1, length(chain1) - di);
    return 0;
}
/*
 * Extend and clip within the range specified by @ext_str and @ext_end;
 * The @tiles_str and @tiles_end are empty string to store the result
   which is within ext_str<= .. <ext_end.
 * Note<red>::The function uses single strand hash, thus seq2 is required to
   be on the same strand of the @ext_str.
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
   ERROR::Don't use this function cause the c_stream doesn't contain strand.
   The newly inserted chain thus doesn't have strand.
 */
uint64_t extendClipRange(String<Dna5> & seq1, 
                         String<Dna5> & seq2, 
                         String<uint64_t> & tiles_str, 
                         String<uint64_t> & tiles_end, 
                         uint64_t ext_str, 
                         uint64_t ext_end, 
                         int direction, 
                         GapParms & gap_parms)
{
    //ext_str = shift_tile(ext_str, -50, -50);
    if (get_cord_strand(ext_str ^ ext_end))
    {
        return 1;
    }
    int thd_best_n = 2;
    int shape_len = 5;// gap_parms.thd_ecr_shape_len;
    int step1 = 3;
    int step2 = 1;
    String<uint64_t> g_hs; 
    String<uint64_t> g_hs_anchors;
    reserve(g_hs, 1024);
    reserve(g_hs_anchors, 1024);
    //int64_t anchor_lower = isClipTowardsLeft(direction) ? g_hs_Cord2StrAnchor(ext_end) - 50 : g_hs_Cord2StrAnchor(ext_str) - 50;
    //int64_t anchor_upper = isClipTowardsLeft(direction) ? g_hs_Cord2StrAnchor(ext_end) + 50 : g_hs_Cord2StrAnchor(ext_str) + 50;

    //create_g_anchors(seq1, seq2, g_hs, g_hs_anchors, ext_str, ext_end, anchor_lower, anchor_upper, shape_len, step1, step2);
    c_stream_(seq1, g_hs, get_tile_x(ext_str), get_tile_x(ext_end), step1, shape_len, 0);
    c_stream_(seq2, g_hs, get_tile_y(ext_str), get_tile_y(ext_end), step2, shape_len, 1);
    c_createAnchors(g_hs, g_hs_anchors, length(g_hs), 3, 50, get_tile_x(ext_str), get_tile_y(ext_str));
    //g_create_anchors_(g_hs, g_hs_anchors, shape_len, anchor_lower, anchor_upper, length(seq2) - 1);
    StringSet<String<uint64_t> > anchors_chains;
    String<int> anchors_chains_score;
    uint thd_chain_depth = 15;
    uint64_t thd_chain_dx_depth = 30;
    std::sort(begin(g_hs_anchors), end(g_hs_anchors), [](uint64_t & a, uint64_t & b){return g_hs_anchor_getX(a) > g_hs_anchor_getX(b);});
    chainAnchorsBase(g_hs_anchors, anchors_chains, anchors_chains_score, 0, length(g_hs_anchors), 
        thd_chain_depth, thd_chain_dx_depth, thd_best_n, gap_parms.chn_ext_clip_metric1, &g_hs_anchor_getX);
    if (empty(anchors_chains))
    {
        return 1;
    }
    //select the chain whose anchor of first (clip2right) or last tile(clip2left) is cloest to the init ext_str or ext_end as the final chain.
    int closest_i = 0;
    int score_i_str = 0;
    int len_sum = 0;
    int64_t min_da = LLMAX;
    for (int i = 0; i < length(anchors_chains) && !empty(anchors_chains[i]); i++)
    {
        uint64_t new_tile = isClipTowardsRight(direction) ? g_hs_anchor2Tile(anchors_chains[i][0]) :
            g_hs_anchor2Tile(back(anchors_chains[i]));
        uint64_t connect_tile = isClipTowardsRight(direction) ? ext_str : ext_end;
        int64_t da = std::abs(int64_t(get_tile_x(new_tile) - get_tile_x(connect_tile) - get_tile_y(new_tile) + get_tile_y(connect_tile)));
        if (da < min_da)
        {
            min_da = da;
            closest_i = i;
            score_i_str = len_sum;
        }
        len_sum += length(anchors_chains[i]);
    }
    if (min_da > gap_parms.thd_ecr_reject_da || length(anchors_chains[closest_i]) < 1)
    {
        return 2; 
    }
    String<int> chain_score;
    resize(chain_score, length(anchors_chains[closest_i]));
    for (int i = 0; i < length(anchors_chains[closest_i]); i++)
    {
        chain_score[i] = anchors_chains_score[score_i_str + i];
    }
    //insert new tiles:
    //start to create new tile from the breakpoint(anchor) towards the ext_str(when clip towards right) 
    //or ext_end (when clip towards left);
    int64_t thd_new_tile_step = gap_parms.thd_tile_size / 2; 

    if (isClipTowardsLeft(direction))
    {
        uint64_t pre_str_tile = 0;
        uint64_t ext_end_y = get_tile_y(ext_end);
        //int clip_i = clipChain (anchors_chains[closest_i], shape_len, direction, false, &g_hs_anchor_getX, &g_hs_anchor_getY, gap_parms);
        for (int i = 0; i < length(anchors_chains[closest_i]); i++)
        {
            uint64_t new_str_tile = g_hs_anchor2Tile(anchors_chains[closest_i][i]);
            if (i == 0 || get_tile_y(new_str_tile) > get_tile_y(pre_str_tile) + thd_new_tile_step)
            {
                appendValue(tiles_str, new_str_tile);
                if (get_tile_y(new_str_tile) + gap_parms.thd_tile_size < ext_end_y) // assure tiles_end_y < ext_end_y 
                {
                    appendValue(tiles_end, shift_tile(new_str_tile, gap_parms.thd_tile_size, gap_parms.thd_tile_size));
                }
                else //otherwise use the end of last anchor as the end of the tile
                {
                    uint64_t last_end_tile = shift_tile(g_hs_anchor2Tile(back(anchors_chains[closest_i])), shape_len, shape_len);
                    appendValue(tiles_end, last_end_tile);
                    break;
                }
                pre_str_tile = new_str_tile;
            }
        }
    }
    else if (isClipTowardsRight(direction)) 
    {
        uint64_t pre_end_tile = 0;
        uint64_t ext_str_y = get_tile_y(ext_str);
        //int clip_i = clipChain (anchors_chains[closest_i], shape_len, direction, false, &g_hs_anchor_getX, &g_hs_anchor_getY, gap_parms);
        for (int i = length(anchors_chains[closest_i]) - 1; i >= 0; i--)
        {
            uint64_t new_end_tile = shift_tile(g_hs_anchor2Tile(anchors_chains[closest_i][i]), shape_len, shape_len);
            if (i == length(anchors_chains[closest_i]) - 1 || 
                get_tile_y(pre_end_tile) > get_tile_y(new_end_tile) + thd_new_tile_step)
            {
                insertValue(tiles_end, 0, new_end_tile);
                if (get_tile_y(new_end_tile) >= ext_str_y + gap_parms.thd_tile_size)
                {
                    insertValue(tiles_str, 0, shift_tile(new_end_tile, -gap_parms.thd_tile_size, -gap_parms.thd_tile_size));
                }
                else
                {
                    uint64_t first_str_tile = g_hs_anchor2Tile(anchors_chains[closest_i][0]);
                    insertValue(tiles_str, 0, first_str_tile);
                    break;
                }
                pre_end_tile = new_end_tile;
            }
        }    
    }
    //Note::red, reomve this when using create_g_anchors; c_stream negelect strand, so add strand here.
    if (get_tile_strand(ext_str))
    {
        for (int i = 0; i < length(tiles_str); i++)
        {
            set_tile_strand (tiles_str[i]);
            set_tile_strand (tiles_end[i]);
        }
    }
    return 0;
}
/********************************************************************/
/**
 * Extend the region around the breakpoint and clip it by gapped pattern.
 * The clip function is splitted into two independent functoins for
 * two @clip_direction value {-1,1} to reduce branches of if and else.
 */
uint64_t c_clip_extend_( uint64_t & ex_d, // results
                    String<uint64_t> & hashs,
                    String<Dna5> & seq1, 
                    String<Dna5> & seq2, 
                    uint64_t extend_str,
                    uint64_t extend_end,
                    int thd_scan_radius,
                    int thd_error_level,
                    int thd_merge_anchor,
                    int thd_drop,
                    int clip_direction)
{    
    CmpInt64 g_cmpll;
    int thd_init_chain_da = 10; 
    int thd_init_chain_num = 6;    
    int thd_init_scan_radius = 20;

    int thd_da_upper = 3;
    int thd_da_lower = -3;
    uint64_t hs_len1 = get_cord_x(extend_end - extend_str);
    uint64_t hs_len2 = get_cord_y(extend_end - extend_str);
    unsigned shape_len = c_shape_len3;
    if (length(hashs) < hs_len1 ||
        hs_len1 < shape_len  || 
        hs_len2 < shape_len)
    {
        return 1;
    }
    int k_str = 0;
    int64_t j_str = 0;
    int64_t j_end = 0;
    int chain_init_len;
    float drop_count = 0;
    String<short> tzs; //trailing zero of each pattern
    String<short> lzs; //leading zero of each pattern
    String<short> chain_x;
    String<short> chain_y;

    Iterator<String<Dna5> >::Type it_str1 = begin(seq1) + get_cord_x(extend_str);
    Iterator<String<Dna5> >::Type it_str2 = begin(seq2) + get_cord_y(extend_str);
    LShape shape(shape_len);
    hashInit_hs(shape, it_str1, 0);
    for (int i = 0; i < hs_len1; i++)
    {
        hashs[i] = hashNext_hs(shape, it_str1 + i);
    }   
    if (isClipTowardsLeft(clip_direction)) //-Gap-Match-
    {
        appendValue(chain_x, hs_len1 - 1);
        appendValue(chain_y, hs_len2 - 1);
        chain_init_len = length(chain_y);
        hashInit_hs(shape, it_str2 + hs_len2 - shape.span, 1);
        for (int64_t i = hs_len2 - shape.span; i > 0; --i) //scan the read 
        {
            clear(tzs);
            clear(lzs);
            bool f_extend = false; 
            uint64_t hash_read = hashPre_hs(shape, it_str2 + i);  
            int64_t di = (hs_len2 - i) >> thd_error_level; 
                    di = std::max((int64_t)thd_scan_radius, di);
            int64_t dj = thd_scan_radius << 1;

//TODO::check boundary of j_str and j_end
            g_cmpll.min(j_end, i + di) << int64_t(chain_x[k_str]) 
                                       << int64_t(hs_len2 - 1)
                                       << int64_t(length(hashs));
            g_cmpll.max(j_str, i - di) >> int64_t(j_end - dj) >> 0 ;
            if (j_str > length(hashs) || j_end > length(hashs))
            {
                return 0;
            }
            //TODO::check boundary of j_str and j_end
            for (int64_t j = j_end; j > j_str; j--) //scan the genome
            {
                uint64_t dhash = hash_read ^ hashs[j];
                appendValue(tzs, ctzb_4_(dhash));
                appendValue(lzs, clzb_4_(dhash));
                short x = j;
                short y = i;
                int m = j_end - j - 1;
                if (m >= 0 && c_isGapMatch_(dhash, tzs[m], tzs[m + 1], lzs[m], lzs[m + 1], c_shape_len3))
                {
                    bool f_first = true;
                    if (hs_len2 - 1 - i < thd_init_chain_num)
                    {
                        if (std::abs(y - x) < thd_init_chain_da)
                        {
                            appendValue(chain_x, x);
                            appendValue(chain_y, y);
                            f_extend = true;
                            chain_init_len = length(chain_y);
                        }
                    }
                    else 
                    {
                        for (int k = k_str; k < length(chain_y); k++)
                        {
                            short dx = chain_x[k] - x;
                            short dy = chain_y[k] - y;
                            short da = dx - dy;
                            if (std::abs(dy) <= shape_len && f_first) 
                            {
                                k_str = k;
                                f_first = false;
                            }
                            if (da < thd_da_upper && da > thd_da_lower) 
                            {
                                appendValue(chain_x, x);
                                appendValue(chain_y, y);
                                f_extend = true;
                                break;
                            }
                        }
                    }
                }
            }
            if (!f_extend)
            {
                if (++drop_count > thd_drop)
                {
                    if (length(chain_x) > chain_init_len) //!empty
                    {
                        ex_d = shift_cord(extend_str, back(chain_x), back(chain_y));
                    }
                    else
                    {
                        ex_d = extend_end;
                    }
                    return ex_d;
                }
            }
            else
            {
                drop_count = std::max (drop_count - 0.5, 0.0);
            }
        }
    }
    else if (isClipTowardsRight(clip_direction)) //-Match-Gap
    {    
        appendValue(chain_x, 0);
        appendValue(chain_y, 0);
        chain_init_len = length(chain_y);
        hashInit_hs(shape, it_str2, 0);
        for (int64_t i = 0; i < hs_len2 - shape.span + 1; i++) //scan the read 
        {
            clear(tzs);
            clear(lzs);
            bool f_extend = false; 
            uint64_t hash_read = hashNext_hs(shape, it_str2 + i);  
            int64_t di = std::max (i >> thd_error_level, int64_t(thd_scan_radius));
            int64_t dj = thd_scan_radius << 1;
            g_cmpll.max(j_str, i - di) >> int64_t(chain_x[k_str]) >> 0 ;
            g_cmpll.min(j_end, i + di) << int64_t(j_str + dj) 
                                       << int64_t(hs_len2) 
                                       << length(hashs);
            if (j_str > length(hashs) || j_end > length(hashs))
            {
                return 0;
            }
            for (int64_t j = j_str; j < j_end; j++) //scan the genome
            {
                uint64_t dhash = hash_read ^ hashs[j];
                appendValue(tzs, ctzb_4_(dhash));
                appendValue(lzs, clzb_4_(dhash));
                short x = j;
                short y = i;
                int m = j - j_str - 1;
                if (m >= 0 && c_isGapMatch_(dhash, tzs[m], tzs[m + 1], lzs[m], lzs[m + 1], c_shape_len3))
                {
                    bool f_first = true;
                    if (i < thd_init_chain_num)
                    {
                        if (std::abs(y - x) < thd_init_chain_da)
                        {
                            appendValue(chain_x, x);
                            appendValue(chain_y, y);
                            chain_init_len = length(chain_y);
                            f_extend = true;
                        }
                    }
                    else 
                    {
                        for (int k = k_str; k < length(chain_y); k++)
                        {
                            short dx = x - chain_x[k];
                            short dy = y - chain_y[k];
                            short da = dx - dy;
                            if (std::abs(dy) <= shape_len && f_first) 
                            {
                                k_str = k;
                                f_first = false;
                            }
                            if (da < thd_da_upper && da > thd_da_lower) 
                            {
                                appendValue(chain_x, x);
                                appendValue(chain_y, y);
                                f_extend = true;
                                break;
                            }
                        }
                    }
                }
            }
            if (!f_extend)
            {
                if (++drop_count > thd_drop)
                {
                    if (length(chain_x) > chain_init_len)
                    {
                        ex_d = shift_cord (extend_str, back(chain_x), back(chain_y));
                        return ex_d;
                    }
                    else
                    {
                        ex_d = extend_str;
                        return ex_d;
                    }
                }
            }
            else
            {
                drop_count = std::max (drop_count - 0.5, 0.0);
            }
        }
    }
    return 0;
}

struct ParmClipExtend
{
    int thd_min_scan_delta;
    int thd_error_level;
    int thd_gap_shape;
    int thd_merge_anchor;
    int thd_merge_drop;   
};

/*
 *@clip_str, @clip_end required to have the equivalent strand
 */
uint64_t c_clip_(String<Dna5> & genome,  
                 String<Dna5> & read,
                 String<Dna5> & comstr,    //complement revers of read
                 uint64_t clip_str,
                 uint64_t clip_end,
                 float thd_band_ratio,
                 int clip_direction = 1)
{
    CmpInt64 g_cmpll;
    uint64_t gs_str = get_tile_x(clip_str);
    uint64_t gr_str = get_tile_y(clip_str);
    uint64_t gs_end = get_tile_x(clip_end);
    uint64_t gr_end = get_tile_y(clip_end);
    uint64_t genomeId = get_tile_id (clip_str);
    uint64_t gr_strand = get_tile_strand(clip_str);
    String<uint64_t> g_hs;
    String<uint64_t> g_anchor;
    reserve(g_hs, 1024);
    reserve(g_anchor, 1024);

    String<Dna5> & seq1 = genome;
    String<Dna5> & seq2 = (gr_strand) ? comstr : read;
    int band = (gs_end - gs_str) * thd_band_ratio;
    ///clip scaffold
 
//step1. extend anchor : currently aborted since the gap map is precies enough
    c_stream_(seq1, g_hs, gs_str, gs_end, 1, c_shape_len, 0);
    c_stream_(seq2, g_hs, gr_str, gr_end, 1, c_shape_len, 1);
    //std::sort (begin(g_hs), end(g_hs));
    int band_level = 3; //>>3 = /8 = * 12.5% error rate 
    c_createAnchors(g_hs, g_anchor, length(g_hs), band_level, band, gs_str, gr_str);

    int thd_merge1 = 3;
    int thd_merge1_lower = 10;
    int thd_merge2 = 20;
    int thd_width = 10;
    uint64_t clip = c_clip_anchors_(g_anchor, clip_str, clip_end, c_shape_len, 
                                    thd_merge1, thd_merge1_lower, thd_merge2,  clip_direction);
    //uint64_t clip = (clip_direction < 0 )? clip_end : clip_str;

//step2. clip_extend gap pattern further.
    int64_t extend_window = 100;
    int64_t thd_ovlp_shift = 10;
    int64_t thd_merge_anchor = 5;
    int64_t thd_merge_drop = 6;
    int64_t thd_error_level = 3;  // >>3 == * 0.125
    int64_t thd_scan_radius = 3;  //at least scan 5 elements in the genome for each kmer in the read
    int64_t dx = 0;
    int64_t dy = 0; 
    uint64_t extend_str;
    uint64_t extend_end;
    uint64_t breakpoints;
    int64_t shift;
    if (isClipTowardsLeft (clip_direction))
    {
        g_cmpll.min(dx, thd_ovlp_shift) << int64_t(get_cord_x(clip_end - clip) - 1);
        g_cmpll.min(dy, thd_ovlp_shift) << int64_t(get_cord_y(clip_end - clip) - 1);
        extend_end = shift_cord (clip, dx, dy);
        g_cmpll.min(shift, extend_window) 
                    << get_tile_x(extend_end - clip_str)
                    << get_tile_y(extend_end - clip_str);
        extend_str = shift_cord (extend_end, -shift, -shift);
    }
    else if (isClipTowardsRight (clip_direction))
    {
        g_cmpll.min(dx, thd_ovlp_shift) << int64_t(get_cord_y(clip - clip_str));
        g_cmpll.min(dy, thd_ovlp_shift) << int64_t(get_cord_y(clip - clip_str));
        extend_str = shift_cord(clip, -dx, -dy);
        g_cmpll.min(shift, extend_window) 
                    << get_tile_x(clip_end - extend_str) 
                    << get_tile_y(clip_end - extend_str);
        extend_end = shift_cord(extend_str, shift, shift);

    }

    c_clip_extend_(clip, g_hs, seq1, seq2, extend_str, extend_end,
                   thd_scan_radius, thd_error_level, thd_merge_anchor, thd_merge_drop,clip_direction);
    return clip;
}

/*
 * Clip exact breakpoints of the given tile at the front or end according to the clip direction.
 */
uint64_t clip_tile (String<Dna5> & seq1,
                    String<Dna5> & seq2,
                    String<Dna5> & comstr,
                    uint64_t tile, 
                    int sv_flag, 
                    int tile_size)
{
    int64_t thd_max_gap_size = tile_size;
    float thd_band_ratio = 0.5;

    CmpInt64 g_cmpll;
    uint64_t clip = EmptyClipConst;
    int64_t shift;
    int clip_direction;
    if (sv_flag & g_sv_r)
    {
        g_cmpll.min(shift) << int64_t(tile_size / 2)
                           << length(seq1) - 1 - get_tile_x(tile) 
                           << length(seq2) - 1 - get_tile_y(tile);
        uint64_t clip_str = shift_tile(tile, shift, shift);
        g_cmpll.min(shift, shift + thd_max_gap_size) 
                           << length(seq1) - 1 - get_tile_x(tile) 
                           << length(seq2) - 1 - get_tile_y(tile);
        uint64_t clip_end = shift_tile(tile, shift, shift);
        int clip_direction = 1;
        clip = c_clip_ (seq1, seq2, comstr, clip_str, clip_end, thd_band_ratio, clip_direction); 
        //clip = clip_str;
    }
    else if (sv_flag & g_sv_l)
    {
        g_cmpll.min(shift) << tile_size / 2 
                           << length(seq1) - 1 - get_tile_x(tile) 
                           << length(seq2) - 1 - get_tile_y(tile);
        uint64_t clip_end = shift_tile(tile, shift, shift);
        g_cmpll.min(shift) << tile_size / 2
                           << get_tile_x(tile)
                           << get_tile_y(tile);
        uint64_t clip_str = shift_tile(tile, -shift, -shift);
        clip_direction = -1;
        clip = c_clip_ (seq1, seq2, comstr, clip_str, clip_end, thd_band_ratio, clip_direction); 
        //clip = clip_end;
        //remove_tile_sgn(clip2);
    }
    remove_tile_sgn(clip);
    return clip;
}

/*----------  Reform tiles of gaps  ----------*/

int isTilesConsecutive(uint64_t & tile1, uint64_t tile2, uint64_t thd_cord_gap)
{
    return isCordsConsecutive_(tile1, tile2, thd_cord_gap);
}

/*
 * @tile_str and tile_end refer to the start and end of the same tile rather than two different cords
 * @thd_tile_size is the regular size of tile, 96x96 by default;
   However the size of the tile specified by the @tile_str and @tile_end is allowed to be smaller than that.
 */
uint64_t reform_tile_ (String<Dna5> & seq1,
                       String<Dna5> & seq2,
                       String<Dna5> & comstr,
                       uint64_t & tile_str,
                       uint64_t & tile_end,
                       int sv_flag, 
                       int thd_tile_size)
{
    int f_e = is_tile_end(tile_str);
    uint64_t clip = clip_tile (seq1, seq2, comstr, tile_str, sv_flag, thd_tile_size);
    int64_t shift;
    if (sv_flag & g_sv_r)
    {
        tile_end = clip;
    }
    else if (sv_flag & g_sv_l)
    {
        tile_str = clip;
    }
    if (f_e)
    {
        set_tile_end(tile_str);
        set_tile_end(tile_end);
    }

    return clip;
}
/*
 * reform, extend and clip tlies of @tiles_str[it] @tiles_end[it]
 * For simplicity, only one tile is allowed to be reform: @tiles_it[it]; 
   Hence in case of towards right, ext_str should >= @tiles_str[it] to make sure the newly reformed tile can be connectted to its predecessors
 * to do:: to reform a series of tiles, add for loop to compare the newly reformed tiles with the predecessors.
 */
int reformExtendClipTile (String<Dna5> & seq1, String<Dna5> & seq2, String<Dna5> & comstr,
                          String<uint64_t> & tiles_str, String<uint64_t> & tiles_end, 
                          int it, int sv_flag, GapParms gap_parms)
{
    int d_it = 0; 
    String<uint64_t> tmp_tiles_str;
    String<uint64_t> tmp_tiles_end;
    String<Dna5> & read = get_tile_strand(tiles_str[it]) ? comstr : seq2;

    if (sv_flag & g_sv_r)
    {
        uint64_t ext_str = tiles_str[it];
        uint64_t ext_end = shift_tile(tiles_end[it], gap_parms.thd_tile_size,  gap_parms.thd_tile_size);
        extendClipRange(seq1, read, tmp_tiles_str, tmp_tiles_end, ext_str, ext_end, g_clip_rght, gap_parms);
        if (!empty(tmp_tiles_str))
        {
            int ii = it;
            while (ii>= 0 && !(ii < it && is_tile_end(tiles_str[ii])) 
                && (get_tile_y(tiles_end[ii]) > get_tile_y(tmp_tiles_end[0]) || 
                    get_tile_x(tiles_end[ii]) > get_tile_x(tmp_tiles_end[0])))
            {ii--;}
            if (ii < it) //erase at least one  
            {
                tmp_tiles_str[0] = tiles_str[ii + 1];
                if (is_tile_end(tiles_str[it]))
                {
                    set_tile_end(back(tmp_tiles_str));
                    set_tile_end(back(tmp_tiles_end));
                }
            }
            //predecessor tiles_end out bound of newly clipped tiles_end.
            
            erase(tiles_str, ii + 1, it + 1);
            erase(tiles_end, ii + 1, it + 1);
            insert(tiles_str, ii + 1, tmp_tiles_str);
            insert(tiles_end, ii + 1, tmp_tiles_end); 
            d_it = length(tmp_tiles_str) - it + ii;
        } 
    }
    else if (sv_flag & g_sv_l)
    {
        uint64_t ext_str = shift_tile(tiles_str[it], -gap_parms.thd_tile_size, -gap_parms.thd_tile_size);
        uint64_t ext_end = tiles_end[it];    
        extendClipRange(seq1, read, tmp_tiles_str, tmp_tiles_end, ext_str, ext_end, g_clip_left, gap_parms);

        if (!empty(tmp_tiles_str))
        {
            int ii = it;
            while (ii < length(tiles_str) && !(ii > 0 && is_tile_end(tiles_str[ii - 1])) && 
                (get_tile_y(tiles_str[ii]) < get_tile_y(back(tmp_tiles_str)) || 
                 get_tile_x(tiles_str[ii]) < get_tile_x(back(tmp_tiles_str))))
            {ii++;}
            if (it < ii)
            {
                back(tmp_tiles_end) = tiles_end[ii - 1]; //keep the tiles_end[ii]
                if (is_tile_end(tiles_str[ii - 1]))
                {
                    set_tile_end(back(tmp_tiles_str));
                    set_tile_end(back(tmp_tiles_end));
                }
            }
            erase(tiles_str, it, ii);
            erase(tiles_end, it, ii);
            insert(tiles_str, it, tmp_tiles_str);
            insert(tiles_end, it, tmp_tiles_end);
            d_it = length(tmp_tiles_str) - std::min(ii - it, 1); //always points the original it, if the original it is erased then point to the last of tmp_tiles_str since the last_tmp_tiles_end is set as the original;
        }
    }
    return d_it;
}

/**
 * Scan @tiles_str[x] @tiles_end[x] of within x:[@pos_str, @pos_end) and reform if 1. head or taile, 2. gaps exists;
 * When (length(tiles_str) == 1), g_map_closed is not allowed
 * else the fist tile of each tile block in @tiles_str is clipped towards left only if @direction == g_map_left;
    The last tile ...towards right... only if == g_map_rght; 
 *  Others clipped towards both direction if gap exists.
 * @tiles_end required to be initied with same length of @tiles_str before call the  functoin
 */
int reform_tiles_(String<Dna5> & seq1, 
                  String<Dna5> & seq2/*read*/, 
                  String<Dna5> & comstr, //complement reverse of seq2
                  String<uint64_t> & tiles_str, 
                  String<uint64_t> & tiles_end,
                  String<uint64_t> & sp_tiles, //record tiles id that needs additional process(dups).
                  int direction, 
                  GapParms & gap_parms)
{

    if (empty (tiles_str))
    {
        return 0; 
    }
    int64_t thd_dxy_min = 80; //skip ins del < this value todo tune later
    float thd_da_zero = gap_parms.thd_err; 
    if (length(tiles_str) == 1 && direction == g_map_closed)
    {
        return 1;
    }
    //for (int i = pos_str; i < pos_end; i++)
    return 0;
    for (int it = 0; it < length(tiles_str); it++)
    {
        uint64_t thd_gap_size = gap_parms.thd_gap_len_min;
        if ((it == 0 || it > 0 && is_tile_end (tiles_str[it - 1])) && direction == g_map_left)
        {
            it += reformExtendClipTile (seq1, seq2, comstr, tiles_str, tiles_end, it, g_sv_l, gap_parms);
        }
        /*
        else if (it > 0)
        {
            uint64_t x11 = get_tile_x(tiles_str[it - 1]);
            uint64_t y11 = get_tile_y(tiles_str[it - 1]);
            uint64_t x12 = get_tile_x(tiles_end[it - 1]);
            uint64_t y12 = get_tile_y(tiles_end[it - 1]);
            uint64_t x21 = get_tile_x(tiles_str[it]);
            uint64_t y21 = get_tile_y(tiles_str[it]);

            if (!(x11 <= x21 && y11 <= y21 && x12 + thd_gap_size >= x21 && y12 + thd_gap_size >= y21 && 
                !get_tile_strand(tiles_str[it - 1] ^ tiles_str[it])))
            {
                it += reformExtendClipTile (seq1, seq2, comstr, tiles_str, tiles_end, it - 1, g_sv_r, gap_parms);
                if (is_diff_anchor (tiles_str[it - 1], tiles_str[it], 1, thd_dxy_min, thd_da_zero)) 
                {
                    insertClipStr(sp_tiles, it - 1);  //tiles[i - 1] and tiles[i] might be dups 
                    insertClipEnd(sp_tiles, it);
                }
                it += reformExtendClipTile (seq1, seq2, comstr, tiles_str, tiles_end, it, g_sv_l, gap_parms);
                //return 0;
            }
        }
        */ 
        if (direction == g_map_rght && (is_tile_end(tiles_str[it]) || it == length(tiles_str) - 1))
        {
            it += reformExtendClipTile (seq1, seq2, comstr, tiles_str, tiles_end, it, g_sv_r, gap_parms);
        }

    }

    return 0;
}

/*
 @direction == g_map_left @gap_str is ommited::will not append to tiles to be clipped
 @direction == g_map_rght @gap_end is ommited::will not append to tiles to be clipped
 */
int reform_tiles(String<Dna5> & seq1, 
                 String<Dna5> & seq2, 
                 String<Dna5> & comstr, 
                 String<uint64_t> & tiles_str, 
                 String<uint64_t> & tiles_end,
                 String<uint64_t> & sp_tiles,
                 uint64_t gap_str, 
                 uint64_t gap_end, 
                 int direction,
                 GapParms & gap_parms)
{
    //step.1 init tiles_str, tiles_end;
    //Insert head_tile and tail tile at the front and end for each block of tiles
    g_print_tiles_(tiles_str, "rft11") ;
    g_print_tiles_(tiles_end, "rft12") ;
    uint64_t head_tile = gap_str;
    uint64_t tail_tile = shift_tile(gap_end, -gap_parms.thd_tile_size, -gap_parms.thd_tile_size); 
    if (get_tile_x(gap_end) > gap_parms.thd_tile_size && get_tile_y (gap_end) > gap_parms.thd_tile_size)
    {
        tail_tile = shift_tile(gap_end, -gap_parms.thd_tile_size, -gap_parms.thd_tile_size);
    }
    else
    {
        tail_tile = head_tile;
    }
    remove_tile_sgn(head_tile);
    //remove_tile_sgn(tail_tile);
    set_tile_end(tail_tile);
    if (!empty(tiles_str))
    {
        copy_tile_sgn(back(tiles_str), tail_tile);
        copy_tile_sgn(tiles_str[0], head_tile);
        remove_tile_sgn(back(tiles_str));
        remove_tile_sgn(tiles_str[0]);
    }
    g_print_tiles_(tiles_str, "rft41") ;
    g_print_tiles_(tiles_end, "rft42") ;
    if (direction != g_map_left)
    {
        insertValue(tiles_str, 0, head_tile);

    }
    if (direction != g_map_rght)
    {
        appendValue(tiles_str, tail_tile);
    }
    g_print_tiles_(tiles_str, "rft21") ;
    g_print_tiles_(tiles_end, "rft22") ;
    int64_t d = shift_tile (0ULL, gap_parms.thd_tile_size, gap_parms.thd_tile_size);
    if (empty(tiles_end))
    {
        resize(tiles_end, length(tiles_str));
        for (int i = 0; i < length(tiles_str); i++)
        {
            tiles_end[i] = tiles_str[i] + d;
        }
    }
    else
    {
        if (direction != g_map_left)
        {
            insertValue(tiles_end, 0, shift_tile(head_tile, gap_parms.thd_tile_size, 
                gap_parms.thd_tile_size));
        }
        if (direction != g_map_rght)
        {
            appendValue(tiles_end, shift_tile(tail_tile, gap_parms.thd_tile_size,
                gap_parms.thd_tile_size));
        }
    }
    //step.2 reform tiles:clip and break block if necessary
    if (gap_parms.f_rfts_clip)
    {
        //reform_tiles_(seq1, seq2, comstr, tiles_str, tiles_end, sp_tiles, direction, 
        //            gap_parms);
    }
    dout << "rftl1" << direction << length(tiles_str) << length(tiles_end) << "\n";
    g_print_tiles_(tiles_str, "rft31") ;
    g_print_tiles_(tiles_end, "rft32") ;

    return 0;
}

/**
 * shortcut to insert @tiles at @cords[@pos] or @cords[@pos + 1] according to the @direction
 * if direction = g_map_left: erase back(@tiles) and insert(cords, @pos), 
 * if direction = g_map_closed erase first and last tiles then insert(cords, @pos)
 * if direction = g_map_rght: erase @tiles[0] and insert(cords, @pos + 1);
 * @tiles is required to have at least 2 tiles;
 * @tile[0] and back(@tiles) are the clipped cords of gap_str, gap_end
 * They will replace the two joint cords between which @tiles are inserted 
 */
int insert_tiles2Cords_(String<uint64_t> & cords, unsigned & pos, String<uint64_t> & tiles,
                        int direction, int thd_max_segs_num)
{
    if ((length(tiles) < 2 && direction == g_map_closed) || empty(tiles))
    {
        return 1;
    }
    int segs_num = 0;
    int return_type1 = 1 << 30;
    for (auto & tile : tiles)
    {
        if (is_tile_end(tile))
        {
            set_cord_end (tile);
            ++segs_num;
        }
    }

    if (segs_num > thd_max_segs_num)
    {
        return segs_num | return_type1;
    }
    String<uint64_t> tmp_tiles;
    if (direction == g_map_left) //insert at front of cords
    {
        uint64_t recd = get_cord_recd(cords[pos]);//set cord flag
        set_tiles_cords_sgns (tiles, recd);
        if (is_cord_block_end (cords[pos]))
        {
            set_cord_end (back(tiles));
        }
        else 
        {
            _DefaultHit.unsetBlockEnd(back(tiles));
        }
        cords[pos] = back(tiles);
        resize (tiles, length(tiles) - 1);
        insert(cords, pos, tiles);
        pos += length(tiles);
        clear(tiles);
    }
    else if (direction == g_map_rght) //insert at end
    {
        uint64_t recd = get_cord_recd(cords[pos]); 
        set_tiles_cords_sgns (tiles, recd);
        uint64_t cordtmp = cords[pos];
        cords[pos] = tiles[0];
        for (int i = 0; i < length(tiles) - 1; i++)
        {
            tiles[i] = tiles[i + 1];
        }
        resize(tiles, length(tiles) - 1);
        insert(cords, pos + 1, tiles);
        pos += length(tiles);
        if (is_cord_block_end(cordtmp))
        {
            _DefaultHit.setBlockEnd(cords[pos]);
        }
        else 
        { 
            _DefaultHit.unsetBlockEnd(cords[pos]);
        }
        clear(tiles);
    }
    else if (direction == g_map_closed)
    {
        uint64_t recd = get_cord_recd(cords[pos]);//set cord flag
        set_tiles_cords_sgns (tiles, recd);
        uint64_t cordtmp = cords[pos];
        cords[pos - 1] = tiles[0];
        cords[pos] = back(tiles);

        if (is_cord_block_end(cordtmp))
        {
            _DefaultHit.setBlockEnd(cords[pos]);
        }
        else
        {
            _DefaultHit.unsetBlockEnd(cords[pos]);
        }
        for (int i = 0; i < length(tiles) - 2; i++)
        {
            tiles[i] = tiles[i + 1];
        }
        resize (tiles, length(tiles) - 2);
        insert(cords, pos, tiles);
        pos += length(tiles);
        clear(tiles);
    }
    return 0;
}

int insert_tiles2Cords_(String<uint64_t> & cords_str, 
                        String<uint64_t> & cords_end,
                        unsigned & pos,
                        String<uint64_t> & tiles_str,
                        String<uint64_t> & tiles_end,
                        int direction,
                        int thd_cord_size,
                        int thd_max_segs_num)
{
    if (empty(cords_end))
    {
        resize (cords_end, length(cords_str));
        int64_t d = shift_cord (0ULL, int64_t(thd_cord_size), int64_t(thd_cord_size));
        for (int i = 0; i < length(cords_str); i++)
        {
            cords_end[i] = cords_str[i] + d;
        }
    }
    unsigned postmp = pos;
    int len = length(tiles_str);
    int len2 = length(cords_str);
    insert_tiles2Cords_(cords_str, pos, tiles_str, direction, thd_max_segs_num);
    insert_tiles2Cords_(cords_end, postmp, tiles_end, direction, thd_max_segs_num);
    return 0;
}
/********************* START: Extend mapping for INS/DEL **********************/
/* The section is the extendsion for the ins/del. 
   Since many ins are repeats. 
   This part is to customize indels.
*/ 
/*
 * Get the overalpped regions of @chain1 and @chain2 
 * NOTE:Resut [pair.first, len1) of chain1, [0, pair.second) of chain2  having overlaps within the region
 */
std::pair<int, int> getExtendsIntervalChainsOverlaps(String<uint64_t> & chain1, 
            String<uint64_t> & chain2, 
            uint64_t(*getX)(uint64_t), 
            uint64_t(*getY)(uint64_t), 
            GapParms & gap_parms)
{
    if (empty(chain1) || empty(chain2))
    {
        return std::pair<int, int>(length(chain1), 0);
    }
    uint64_t x2 = getX(chain2[0]);
    uint64_t y2 = getY(chain2[0]);
    x2 = (x2 > gap_parms.thd_dcomx_err_dx) ? x2 - gap_parms.thd_dcomx_err_dx : 0;
    y2 = (y2 > gap_parms.thd_dcomx_err_dy) ? y2 - gap_parms.thd_dcomx_err_dy : 0; 
    int i1 = 0;
    for (int i = length(chain1) - 1; i >= 0; i--)     
    {
        if (getX(chain1[i]) < x2 && getY(chain1[i]) < y2)
        {
            i1 = i + 1;
            break;
        }
    }
    uint64_t x1 = getX(back(chain1)) + gap_parms.thd_dcomx_err_dx;
    uint64_t y1 = getY(back(chain1)) + gap_parms.thd_dcomx_err_dy;
    x1 = (gap_parms.ref_len - x1 > gap_parms.thd_dcomx_err_dx) ? x1 + gap_parms.thd_dcomx_err_dx 
            : gap_parms.ref_len;
    y1 = (gap_parms.read_len - y1 > gap_parms.thd_dcomx_err_dy) ? y1 + gap_parms.thd_dcomx_err_dy 
            : gap_parms.read_len;
    int i2 = 0;
    for (int i = 0; i < length(chain2); i++)     
    {
        if (getX(chain2[i]) > x1 && getY(chain2[i]) > y1)
        {
            i2 = i;
            break;
        }
    }
    return std::pair<int, int> (i1, i2);
}
/* 
 * Generic function
 * Given @chains, map continously with smaller patterns along @chains
 * @map_str, @map_end are required to be on the same strand
 * @chains is on the single strand,(values of different strand of @chain aren't allowed) 
 * Map along @chains[i], where i within [@i_str, @i_end) 
   the result is written to the @tiles
 */
int mapAlongChain(String<Dna5> & seq1, 
                  String<Dna5> & seq2, 
                  String<uint64_t> & chains, 
                  String<uint64_t> & tiles,  
                  int i_str, int i_end, 
                  int shape_len, 
                  int step1, 
                  int step2, 
                  uint64_t(*getX)(uint64_t), 
                  uint64_t(*getY)(uint64_t), 
                  uint64_t(*getStrand)(uint64_t),
                  void    (*setStrand) (uint64_t &),
                  uint64_t(*mac_chain2Tile)(uint64_t), 
                  GapParms & gap_parms)
{
    dout << "macs2" << i_str << i_end << length(chains) << step1 << "\n";
    if (empty(chains) || i_str < 0 || i_end > length(chains) || i_end <= i_str)
    {
        return -1;
    }
    int thd_best_n = 1;
    g_print_tiles_(chains, "macs1");
    String<uint64_t> hs; 
    String<uint64_t> hs_anchors;
    reserve(hs, 1024);
    reserve(hs_anchors, 1024);
    int64_t anchor_str = getX(chains[i_str]) - getY(chains[i_str]);
    int64_t anchor_end = getX(chains[i_end - 1]) - getY(chains[i_end - 1]);
    c_stream_(seq1, hs, getX(chains[i_str]), getX(chains[i_end - 1]), step1, shape_len, 0);
    c_stream_(seq2, hs, getY(chains[i_str]), getY(chains[i_end - 1]), step2, shape_len, 1);
    c_createAnchors2(hs, hs_anchors, length(hs), std::min(anchor_str, anchor_end) - 30, std::max(anchor_str, anchor_end) + 30);
    std::sort(begin(hs_anchors), end(hs_anchors), [](uint64_t & a, uint64_t & b){return g_hs_anchor_getX(a) > g_hs_anchor_getX(b);});
    stickMainChain(hs_anchors, chains, &g_hs_anchor_getX, &g_hs_anchor_getY, getX, getY, gap_parms);
    StringSet<String<uint64_t> > anchors_chains;
    String<int> anchors_chains_score;
    uint thd_chain_depth = 15;
    uint64_t thd_chain_dx_depth = 30;

    chainAnchorsBase(hs_anchors, anchors_chains, anchors_chains_score, 0, length(hs_anchors), 
        thd_chain_depth, thd_chain_dx_depth, thd_best_n, gap_parms.chn_ext_clip_metric1, &g_hs_anchor_getX); 
    if (!empty (anchors_chains))
    {//choose the first chain
        int f_strand = getStrand(chains[0]);
        for (int i = 0; i < length(anchors_chains[0]); i++)
        {
            uint64_t new_tile = mac_chain2Tile(anchors_chains[0][i]);
            if (f_strand)
            {
                setStrand(new_tile);
            }
            appendValue (tiles, new_tile);
        }
    }
    return 0;
}
/*
 * Find and clip the breakpoint of two chains as ins/del.
   The function keeps x of the breakpoint two chains in case of ins and y of breakpoint two chains
   in case of dels closed to each other and finding the breakpoint such that the score is
   maxmized.
 * NOTE: getX and getY are abstract operations. Swap getX and getY to clip by y, namely call
   extendsIn..laps(chain1, chain2, getY, getX): choose y as the main coordinate in case of del
 */
int __extendsIntervalClipOverlapsInsDel_(String<uint64_t> & chain1, 
            String<uint64_t> & chain2, 
            int shape_len, 
            int step1, 
            int step2, 
            bool f_clip, 
            uint64_t(*getX)(uint64_t), 
            uint64_t(*getY)(uint64_t), 
            GapParms & gap_parms)
{
    if (empty(chain1) || empty(chain2))
    {
        return 0;
        //return std::pair<int, int>(length(chain1), 0);
    }
    String<int> gaps_score11;
    String<int> gaps_score12;
    String<int> gaps_score21; 
    String<int> gaps_score22;
    gap_parms.clipChainParms(shape_len, gap_parms.thd_err); //init clip parms
    accumulateSimpleGapScore1(chain1, gaps_score11, shape_len, getX, gap_parms);
    accumulateSimpleGapScore1(chain1, gaps_score12, shape_len, getY, gap_parms);
    accumulateSimpleGapScore1(chain2, gaps_score21, shape_len, getX, gap_parms);
    accumulateSimpleGapScore1(chain2, gaps_score22, shape_len, getY, gap_parms);
    clipChain_(chain1, gaps_score11, gaps_score12, g_map_rght, true, getX, getY, gap_parms);
    clipChain_(chain2, gaps_score21, gaps_score22, g_map_left, true, getX, getY, gap_parms);

    int j1 = 0, j2 = 0, i_clip = 0, j_clip = -1;  
    int j1_pre = 0, j2_pre = 0;
    int score11, score12, score21, score22, score1, score2;
    int min_score = INT_MAX;
    uint64_t x21 = getX(chain2[0]), x22 = getX(chain2[0]); //[x21, x22)
    for (int i = 0; i < length(chain1); i++)        
    {
        uint64_t x1 = getX(chain1[i]);
        uint64_t x2_lower = x1;
        uint64_t x2_upper = x1 + gap_parms.thd_eicos_clip_dxy;
        for (int j = j1_pre; j < length(chain2) && x21 < x2_lower; j++)
        {
            x21 = getX(chain2[j]); 
            j1 = j;
        }
        if (x21 > x2_upper) 
        {
            continue;
        }
        if (x21 < x2_lower)
        {
            break;
        }
        for (int j = j2_pre; j < length(chain2) && x22 <= x2_upper; j++)
        {
            x22 = getX(chain2[j]);
            j2 = j;
        }
        if (x22 < x2_lower)
        {
            break;
        }
        //x22 < x2_upper occur only when j2 == length(chain2) - 1, while this is allowed;
        //x2_lower<= x21 < x2_upper, x22 >= x2_upper &&  != pre(newly updated)
        //then search max score within [x21, x22), namely [j1, j2)
        if (j1 > j_clip || j2_pre != j2)
        {
            int ii1 = std::max(i - gap_parms.thd_eicos_window_size, int(0));
            //int ii2 = std::min(i + gap_parms.thd_eicos_window_size, int(length(chain1) - 1));
            //score11 = gaps_score11[i] - gaps_score11[ii1];
            score11 = gaps_score11[i];
            score12 = gaps_score12[i];
            for (int j = std::max(j1, j2_pre); j < j2; j++)
            {
                //int jj1 = std::max(j - gap_parms.thd_eicos_window_size + 1, int(0));
                int jj2 = std::min(j + gap_parms.thd_eicos_window_size, int(length(chain2) - 1));
                //score21 = gaps_score21[jj2] - gaps_score21[j]; 
                score21 = back(gaps_score21) - gaps_score21[j];
                score22 = back(gaps_score22) - gaps_score22[j];
                int score_connect = getX(chain2[j]) - getX(chain1[i]) > shape_len ? (getX(chain2[j]) - getX(chain1[i]) - shape_len) * gap_parms.int_precision : 0;
                int score = score11 + score12 + score21 + score22 + score_connect;
                if (score < min_score)
                {
                    min_score = score;
                    i_clip = i;
                    j_clip = j;
                }
            }
        }        
        j1_pre = j1;
        j2_pre = j2;
    }
    if (f_clip)
    {
        resize(chain1, i_clip); 
        j_clip = j_clip < 0 ? 0 : j_clip;
        erase(chain2, 0, j_clip);
        //return std::pair<int, int>(i_clip, 0);
        return 0;
    }
    else
    {
        //return std::pair<int, int> (i_clip, j_clip);
        return 0;
    }
    return 0;
}
/*
 * Wrapper to clip the overlaps of @chain1 and @chain2.
   The method is specifically for ins and dels.
 */
int extendsIntervalClipOverlapsInsDel_(String<uint64_t> & chain1, String<uint64_t> & chain2, 
    int shape_len, int step1, int step2, uint64_t(*getX)(uint64_t), uint64_t(*getY)(uint64_t),
     GapParms & gap_parms)
{
    if (empty(chain1) && empty(chain2))
    {
        return 0;
    }
    else if (empty(chain1))
    {
        clipChain (chain2, shape_len, g_map_left, true, getX, getY, gap_parms);
    }
    else if (empty(chain2))
    {
        clipChain (chain1, shape_len, g_map_rght, true, getX, getY, gap_parms);
    }
    else 
    {
        if (!gap_parms.thd_eicos_f_as_ins)
        {
            clipChain (chain1, shape_len, g_map_rght, true, getX, getY, gap_parms);
            clipChain (chain2, shape_len, g_map_left, true, getX, getY, gap_parms);
        }
        else  //clip as indel rather than dups
        {
            __extendsIntervalClipOverlapsInsDel_(chain1, chain2, shape_len, step1, step2, true, getX, getY, gap_parms);
        }
    }
    return 0;
}
/*
 * @direction < 0, find the record chain[i] in chain that x-x0 > dx andd y-y0> dy
   Then rechain the records from chain[0] to chain[i] and clip the new chain.
 *
int extendClipChainGaps(String<Dna5> & ref, 
                    String<Dna5> & read, 
                    String<Dna5> & comstr, 
                    String<uint64_t> & chain,
                    int i_str,
                    uint64_t dx,
                    uint64_t dy,
                    int direction,
                    GapParms & gap_parms)
{
    if (empty(chain))
    {
        return 0;
    }
    if (direction < 0) //left
    {
        uint64_t x0 = get_cord_x(chain[0]);
        uint64_t y0 = get_cord_y(chain[0]);
        for (int i = 0; i < length(chain); i++)
        {
            if (get_cord_x(chain[i]) - x0 >= dx && 
                get_cord_y(chain[i]) - y0 >= dy)
            {
                extendClipChain_(ref, read, comstr, chain, 0, i, direction, gap_parms);
                break;
            }
        }
    }
    else if (direction > 0)
    {
        uint64_t x1 = get_cord_x(back(chain));
        uint64_t y1 = get_cord_y(back(chain));
        for (int i = length(chain) - 1; i >= 0; i--)
        {
            if (get_cord_x(chain[i]) - x1 >= dx && 
                get_cord_y(chain[i]) - y1 >= dy)
            {
                extendClipChain_(ref, read, comstr, chain, i, length(chain), 
                    direction, gap_parms);
                break;
            }
        }
    }
    return 0;
}
*/

/*
 *Drop at the breakpoints of two overlapped chains as ins or dels such that the score is maximized
 *@chains1 direction = right, @chains direction = left
 */
int extendsIntervalMapOverlaps_(String<Dna5> & ref, 
                                String<Dna5> & read, 
                                String<Dna5> & comstr, 
                                String<uint64_t> & tiles1, 
                                String<uint64_t> & tiles2, 
                                uint64_t gap_str1, 
                                uint64_t gap_end1, 
                                uint64_t gap_str2, 
                                uint64_t gap_end2,
                                int shape_len, 
                                int step1, 
                                int step2, 
                                GapParms & gap_parms)
{
    dropChainGapX(tiles1, &get_tile_x, &get_tile_y, g_map_rght, true, gap_parms);
    dropChainGapX(tiles2, &get_tile_x, &get_tile_y, g_map_left, true, gap_parms);
    String<uint64_t> overlap_tiles1;
    String<uint64_t> overlap_tiles2; 
    std::pair<int, int> overlaps = getExtendsIntervalChainsOverlaps(tiles1, tiles2, &get_tile_x, 
        & get_tile_y, gap_parms);
    if (!empty(tiles1))
    {
        String<Dna5> & seq2 = get_tile_strand(tiles1[0]) ? comstr : read;
        mapAlongChain(ref, seq2, tiles1, overlap_tiles1, overlaps.first, length(tiles1), shape_len, 
            step1, step2, &get_tile_x, &get_tile_y, &get_tile_strand, &set_tile_strand,
             &g_hs_anchor2Tile, gap_parms);
    }
    if (!empty(tiles2))
    {
        String<Dna5> & seq2 = get_tile_strand(tiles2[0]) ? comstr : read;
        mapAlongChain(ref, seq2, tiles2, overlap_tiles2, 0,  overlaps.second, shape_len, step1, step2, 
            &get_tile_x, &get_tile_y, &get_tile_strand, &set_tile_strand,
             &g_hs_anchor2Tile, gap_parms);
    }
    if (get_tile_x(gap_str1) - get_tile_y(gap_str1) > get_tile_x(gap_end2) - get_tile_y(gap_end2))
    {
        extendsIntervalClipOverlapsInsDel_(overlap_tiles1, overlap_tiles2, shape_len, step1, step2, &get_tile_x, &get_tile_y, gap_parms); //ins, don't wrapper the shape_len here into gap_parms
    }
    else
    {
        extendsIntervalClipOverlapsInsDel_(overlap_tiles1, overlap_tiles2, shape_len, step1, step2, &get_tile_y, &get_tile_x, gap_parms); //del
    }
    dout << "ovlap" << length(overlap_tiles1) << length(overlap_tiles2) << get_cord_x(gap_str1) << get_cord_x(gap_end2) << "\n";
    resize(tiles1, overlaps.first);
    if (!empty(overlap_tiles1))
    {
        append(tiles1, overlap_tiles1);
    }
    erase(tiles2, 0, overlaps.second);
    if (!empty(overlap_tiles2))
    {
        insert(tiles2, 0, overlap_tiles2);
    }

    return 0;
}
/*
 * Supposed to call in @extendsInterval()
 * step1 map(create chain) for each end of the interval direction sperately
 * step2 remap of the overlap of the 2 chains of different direction : extendsIntervalMapOverlaps_
 * step3 clip the overlaps
 * step4 create tiles along the chains
 * NOTE:If clip is ins or del is according to @gap_str1, @gap_end1, @gap_str2, @gap_end2.
   Hence make sure @gap_str[1/2] and @gap_end[1/2] can correctly indicate ins/del.
 */
int extendsTilesFromAnchors (String<Dna5> & ref, 
                             String<Dna5> & read, 
                             String<Dna5> & comstr,
                             String<uint64_t> & anchors1, 
                             String<uint64_t> & anchors2, 
                             String<uint64_t> & tiles1, 
                             String<uint64_t> & tiles2,
                             StringSet<FeaturesDynamic> & f1, 
                             StringSet<FeaturesDynamic> & f2,
                             uint64_t gap_str1, 
                             uint64_t gap_end1, 
                             uint64_t gap_str2, 
                             uint64_t gap_end2,
                             uint64_t read_len,
                             GapParms & gap_parms)
{
    int original_direction = gap_parms.direction;
    int direction1 = g_map_rght;
    int direction2 = g_map_left;
    String<uint64_t> tmp_tiles1;
    String<uint64_t> tmp_tiles2;
    //!map right part
    gap_parms.direction = direction1;
    g_CreateChainsFromAnchors_(anchors1, tmp_tiles1, gap_str1, gap_end1, read_len, gap_parms);
    getClosestExtensionChain_(tmp_tiles1, gap_str1, gap_end1, true, gap_parms);
    //!map left part
    gap_parms.direction = direction2;
    g_CreateChainsFromAnchors_(anchors2, tmp_tiles2, gap_str2, gap_end2, read_len, gap_parms);
    getClosestExtensionChain_(tmp_tiles2, gap_str2, gap_end2, true, gap_parms);
    //!find and clip at the common breakpoint of the left and right chains
    int shape_len = gap_parms.thd_etfas_shape_len;
    int step1 = gap_parms.thd_etfas_step1;
    int step2 = gap_parms.thd_etfas_step2;
    extendsIntervalMapOverlaps_(ref, read, comstr, tmp_tiles1, tmp_tiles2, gap_str1, gap_end1, gap_str2, gap_end2,  shape_len, step1, step2, gap_parms);
    g_CreateTilesFromChains_(tmp_tiles1, tiles1, f1, f2, gap_str1, 0, length(tmp_tiles1), &get_tile_x, &get_tile_y, &get_tile_strand, gap_parms);    
    trimTiles(tiles1, f1, f2, gap_str1, gap_end2, read_len - 1, direction1, gap_parms);
    g_CreateTilesFromChains_(tmp_tiles2, tiles2, f1, f2, gap_str2, 0, length(tmp_tiles2), &get_tile_x, &get_tile_y, &get_tile_strand, gap_parms);  
    trimTiles(tiles2, f1, f2, gap_str1, gap_end2, read_len - 1, direction2, gap_parms);

    gap_parms.direction = original_direction;
    return 0;
}

/**
 * [@gap_str, @gap_end) to create a chain of tiles to extend the
   mapping area as long as possible.
 * @gap_str and @gap_end should have the same strand
 */
int extendsInterval(String<Dna5> & ref, //genome
                 String<Dna5> & read, //read
                 String<Dna5> & comstr,
                 String<uint64_t> & tiles1,    //results
                 String<uint64_t> & tiles2,    //results
                 StringSet<FeaturesDynamic > & fts_ref,  
                 StringSet<FeaturesDynamic > & fts_read,
                 uint64_t gap_str1,
                 uint64_t gap_end1, 
                 uint64_t gap_str2,
                 uint64_t gap_end2, 
                 GapParms & gap_parms) // extern parm
{
    if (get_cord_strand (gap_str1 ^ gap_end1) || get_cord_strand (gap_str2 ^ gap_end2) || get_cord_strand (gap_str1 ^ gap_str2))
    {
        return 1;
    }
    int original_direction = gap_parms.direction;
    int shape_len = gap_parms.thd_eis_shape_len; 
    int step1 = gap_parms.thd_eis_step1; //seq1 pattern step
    int step2 = gap_parms.thd_eis_step2; //seq2...
    String<uint64_t> g_hs;
    String<uint64_t> g_hs_anchors1;
    String<uint64_t> g_hs_anchors2;
    reserve(g_hs, 2048);
    reserve(g_hs_anchors1, 2048);
    reserve(g_hs_anchors2, 2048);

    int64_t thd_max_extend2 = 5000; //normal case 

    int direction1 = g_map_rght;
    int direction2 = g_map_left;
    uint64_t id = get_cord_id(gap_str1);
    uint64_t strand = get_cord_strand(gap_str1);
    uint64_t x1 = std::min(get_cord_x(gap_str1), get_cord_x(gap_str2));
    uint64_t y1 = std::min(get_cord_y(gap_str1), get_cord_y(gap_str2));
    uint64_t x2 = std::max(get_cord_x(gap_end1), get_cord_x(gap_end1));
    uint64_t y2 = std::max(get_cord_y(gap_end1), get_cord_y(gap_end2));
    uint64_t stream_str = create_cord(id, x1, y1, strand);
    uint64_t stream_end = create_cord(id, x2, y2, strand);
    double t1 = sysTime();
    g_stream_(ref, read, g_hs, stream_str, stream_end, shape_len, step1, step2, gap_parms);
    t1 = sysTime() - t1;
    double t2 = sysTime();
    g_CreateExtendAnchorsPair_(g_hs, g_hs_anchors1, g_hs_anchors2, shape_len, length(read) - 1, gap_str1, gap_end1, gap_str2, gap_end2, gap_parms);
    t2 = sysTime() - t2;
    double t3 = sysTime();
    extendsTilesFromAnchors(ref, read, comstr, g_hs_anchors1, g_hs_anchors2, tiles1, tiles2, fts_ref, fts_read, gap_str1, gap_end1, gap_str2, gap_end2, length(read), gap_parms);
    t3 = sysTime() - t3;
    //--direction = 1 part;
    /*
    gap_parms.direction = direction1;
    mapTilesFromAnchors (g_hs_anchors1, tiles1, fts_ref, fts_read, gap_str1, gap_end1, length(read) - 1, direction1, gap_parms);
    //--diection = -1 part;
    gap_parms.direction = direction2;
    mapTilesFromAnchors (g_hs_anchors2, tiles2, fts_ref, fts_read, gap_str2, gap_end2, length(read) - 1, direction2, gap_parms);
    */
    gap_parms.direction = original_direction;
    return 0;
}
/*
 * Remap the chain towards one direction with shorter patterns and clip the well mapped part
   @direction < 0: remap region [@chain[0], @chain[i_end]);
   @direction > 0: remap region [@chain[i_str], back(@chian));
 */
int remapChainOneEnd(String<Dna5> & ref, 
                     String<Dna5> & read, 
                     String<Dna5> & comstr, 
                     String<uint64_t> & chain, 
                     int shape_len, 
                     int step1, 
                     int step2, 
                     int remap_num,
                     int direction,
                     uint64_t (*getChainX) (uint64_t),
                     uint64_t (*getChainY) (uint64_t),
                     uint64_t (*getChainStrand) (uint64_t),
                     void     (*setChainStrand) (uint64_t &),
                     uint64_t (*anchor2Chain) (uint64_t),
                     GapParms & gap_parms)
{
    if (!direction || empty(chain))
    {
        return 0;
    }
    dout << "rcoe1" << length(chain) << direction << "\n";
    g_print_tiles_(chain, "rcoe11");
    //dropChainGapX(chain, getChainX, getChainY, direction, true, gap_parms);
    String<Dna5> & seq2 = getChainStrand(chain[0]) ? comstr : read;
    String<uint64_t> remap_chain; 
    int i_str, i_end;
    dout << "rcoe4" << direction << isClipTowardsLeft(direction) << isClipTowardsRight(direction) << "\n";
    if (isClipTowardsLeft(direction))
    {
        i_str = std::max(0, int(length(chain) - remap_num));
        i_end = length(chain);
        dout << "rcoe51" << i_end << "\n";
    }
    else if (isClipTowardsRight(direction))
    {
        i_str = 0;
        i_end = std::min(int(length(chain)), remap_num);
        dout << "rcoe52" << i_end << "\n";
    }
    dout << "rcoe4" << direction << i_end << length(chain) << isClipTowardsRight(direction) << isClipTowardsLeft(direction) << "\n";
    mapAlongChain(ref, seq2, chain, remap_chain, i_str, i_end, shape_len, step1,
             step2, getChainX, getChainY, getChainStrand, setChainStrand, anchor2Chain, gap_parms); 
    g_print_tiles_(remap_chain, "rcoe2");
    clipChain(remap_chain, shape_len, direction, true, getChainX, getChainY, gap_parms);
    g_print_tiles_(remap_chain, "rcoe3");
    if (isClipTowardsLeft(direction))
    {
        erase(chain, 0, i_end);
        if (!empty(remap_chain))
        {
            insert(chain, 0, remap_chain);
        }
    }
    else if (isClipTowardsRight(direction))
    {
        if (!empty(remap_chain))
        {
            resize(chain, i_str);
            append(chain, remap_chain);
        }
    }
    return 0;
}
/*
 * Wrapper to call function remapChainOneEnd at @chain[i_ptr_str] or @chain[i_ptr_end]
   if direction < 0  
   The region [@chain[@i_ptr_str]-lower, @chain[@i_ptr_str] + upper] will be remapped.
   if direction > 0
   The region [@chain[@i_ptr_end]-lower, @chain[@i_ptr_end] + upper] will be remapped.
   And the new chain will replace the original ones which is in the region.
 * The function returns the increased length of chain after re-extending.
   The return value >= -length(chain), when the chain is all erased, return value == -length(chain) 
   Hence i_ptr_end + reExtendChainOneSide >= -1,  
   Take care of  out of bound of i_ptr_end == -1 when iterating 
 */
int reExtendChainOneSide(String<Dna5> & ref, 
                        String<Dna5> & read, 
                        String<Dna5> & comstr, 
                        String<uint64_t> & chain, 
                        int i_ptr_str,
                        int i_ptr_end,
                        int lower,
                        int upper,
                        int shape_len, 
                        int step1, 
                        int step2, 
                        int direction,
                        uint64_t (*getChainX) (uint64_t),
                        uint64_t (*getChainY) (uint64_t),
                        uint64_t (*getChainStrand) (uint64_t),
                        void     (*setChainStrand) (uint64_t &),
                        uint64_t (*shiftChain)(uint64_t const &, int64_t, int64_t),
                        uint64_t (*anchor2Chain) (uint64_t),
                        GapParms & gap_parms)
{
    if (empty(chain) || i_ptr_str < 0 || i_ptr_end < 0)
    {
        return 0;
    }
    int ii, i_str, i_end;
    int len = length(chain);
    String <uint64_t> reextend_chain;
    g_print_tiles_(chain, "recos1");
    if (isClipTowardsLeft(direction))
    {
        int64_t d = -std::min({int64_t(get_cord_x(chain[i_ptr_str])), 
                              int64_t(get_tile_y(chain[i_ptr_str])),
                              int64_t(lower)});
        for (ii = i_ptr_str; ii < i_ptr_end; ii++)
        {
            if (get_tile_x(chain[ii]) - get_cord_x(chain[i_ptr_str]) >= upper) 
            {
                break;
            } 
        }
        //ii = std::min(int(i_ptr_end - 1), ii);
        resize (reextend_chain, ii - i_ptr_str + 2);
        reextend_chain[0] = shiftChain(chain[i_ptr_str], d, d); //insert the lower bound to extend to
        for (unsigned i = 0; i < ii - i_ptr_str + 1; i++)
        {
            reextend_chain[i + 1] = chain[i_ptr_str + i];
        }
        i_str = i_ptr_str;
        i_end = ii + 1;
        g_print_tiles_(reextend_chain, "recos21");
    }
    else if (isClipTowardsRight(direction))
    {
        int d = std::min({int64_t(length(ref) - get_cord_x(chain[i_ptr_end]) - 1),
                          int64_t(length(read) - get_cord_y(chain[i_ptr_end]) - 1), 
                          int64_t(upper)}); 
        for (ii = i_ptr_end; ii > i_ptr_str; ii--)
        {
            if (get_tile_x(chain[i_ptr_end]) - get_tile_x(chain[ii]) >= lower)
            {
                break;
            }
        }
        //ii = std::min(int(length(chain)) - 1, ii);
        resize (reextend_chain, i_ptr_end - ii + 2);
        for (unsigned i = 0; i < i_ptr_end - ii + 1; i++)
        {
            reextend_chain[i] = chain[ii + i];
        }
        back(reextend_chain) = shiftChain(chain[i_ptr_end], d, d);
        i_str = ii;
        i_end = i_ptr_end + 1;
        g_print_tiles_(reextend_chain, "recos22");
    }
    remapChainOneEnd(ref, read, comstr, reextend_chain, shape_len, step1, step2, 
            length(reextend_chain), direction,
            getChainX, getChainY, getChainStrand, setChainStrand, anchor2Chain, gap_parms);
    erase(chain, i_str, i_end);
    insert(chain, i_str, reextend_chain);
    
    dout << "ers1 " << i_ptr_str << i_ptr_end  << " s " << int(length(chain) - len) << " " << direction << i_str << i_end << length(reextend_chain) << "\n";
    return length(chain) - len;
}
/*
 * extends tiles towards one direction
 */
int extendTilesOneSide(String<Dna5> & ref, 
                       String<Dna5> & read, 
                       String<Dna5> & comstr,
                       String<uint64_t> & anchors, 
                       String<uint64_t> & tiles1, 
                       StringSet<FeaturesDynamic> & f1, 
                       StringSet<FeaturesDynamic> & f2,
                       uint64_t gap_str, 
                       uint64_t gap_end, 
                       uint64_t read_len,
                       int direction,
                       GapParms & gap_parms)
{
    int original_direction = gap_parms.direction;
    String<uint64_t> chain;
    gap_parms.direction = direction;
    g_CreateChainsFromAnchors_(anchors, chain, gap_str, gap_end, read_len, gap_parms);
    getClosestExtensionChain_(chain, gap_str, gap_end, true, gap_parms);

    //!find and clip at the common breakpoint of the left and right chains
    int shape_len = gap_parms.thd_etfas_shape_len;
    int step1 = gap_parms.thd_etfas_step1;
    int step2 = gap_parms.thd_etfas_step2;
    int remap_num = 50;
    remapChainOneEnd(ref, read, comstr, chain, shape_len, step1, step2, remap_num,
        direction, &get_tile_x, &get_tile_y, &get_tile_strand, &set_tile_strand,
         &g_hs_anchor2Tile, gap_parms);
    g_CreateTilesFromChains_(chain, tiles1, f1, f2, gap_str, 0, length(chain), &get_tile_x, &
        get_tile_y, &get_tile_strand, gap_parms);    
    trimTiles(tiles1, f1, f2, gap_str, gap_end, read_len - 1, direction, gap_parms);
    gap_parms.direction = original_direction;
    return 0;
}
int extendIntervalOneSide(String<Dna5> & ref, //genome
                 String<Dna5> & read, //read
                 String<Dna5> & comstr,
                 String<uint64_t> & tiles,    //results
                 StringSet<FeaturesDynamic > & fts_ref,  
                 StringSet<FeaturesDynamic > & fts_read,
                 uint64_t gap_str,
                 uint64_t gap_end,
                 int direction, 
                 GapParms & gap_parms) // extern parm
{
    if (get_cord_strand (gap_str ^ gap_end))
    {
        return 1;
    }
    int original_direction = gap_parms.direction;
    int shape_len = gap_parms.thd_eis_shape_len; 
    int step1 = gap_parms.thd_eis_step1; //seq1 pattern step
    int step2 = gap_parms.thd_eis_step2; //seq2...
    gap_parms.direction = direction;
    String<uint64_t> g_hs;
    String<uint64_t> g_hs_anchors;
    reserve(g_hs, 2048);
    reserve(g_hs_anchors, 2048);

    g_stream_(ref, read, g_hs, gap_str, gap_end, shape_len, step1, step2, gap_parms);
    g_create_anchors_(g_hs, g_hs_anchors, shape_len, direction, 0, 0, length(read) - 1, gap_str, gap_end, gap_parms);
    extendTilesOneSide(ref, read, comstr, g_hs_anchors, tiles, fts_ref, fts_read, 
        gap_str, gap_end, length(read), direction, gap_parms);
    gap_parms.direction = original_direction;
    return 0;
}

int mapExtendResultFilter_(String<uint64_t> & tiles_str, uint64_t gap_str, uint64_t gap_end, int direction, GapParms & gap_parms)
{
    if (isClipTowardsRight(direction))
    {
        uint64_t pre_tile = gap_str;
        for (int i = 0; i < length(tiles_str); i++)
        {
            int64_t dy = get_cord_y(tiles_str[i]) - get_tile_y(pre_tile);
            int64_t dx = get_cord_y(tiles_str[i]) - get_tile_x(pre_tile);
            if (dy > gap_parms.thd_me_reject_gap || dx > gap_parms.thd_me_reject_gap)
            {
                erase(tiles_str, i, length(tiles_str));
                break;
            }
            pre_tile = tiles_str[i];
        }
    }
    if (isClipTowardsLeft(direction))
    {
        uint64_t pre_tile = gap_end;
        for (int i = length(tiles_str) - 1; i >= 0; i--)
        {
            int64_t dy = get_cord_y(pre_tile) - get_tile_y(tiles_str[i]);
            int64_t dx = get_cord_y(pre_tile) - get_tile_x(tiles_str[i]);
            if (dy > gap_parms.thd_me_reject_gap || dx > gap_parms.thd_me_reject_gap)
            {
                erase(tiles_str, 0, i + 1);
                break;
            }
            pre_tile = tiles_str[i];
        }
    }

    return 0;
}
/*
 * map from @gap_str to @gap_end if direction > 0
   or from @gap_end to @gap_str if direction < 0  
 */
int mapExtend(StringSet<String<Dna5> > & seqs, 
              String<Dna5> & read, String<Dna5> & comstr,
              StringSet<FeaturesDynamic > & f1, 
              StringSet<FeaturesDynamic > & f2,
              String<uint64_t> & tiles_str, 
              String<uint64_t> & tiles_end, 
              uint64_t gap_str, 
              uint64_t gap_end, 
              int direction,
              GapParms & gap_parms)
{
    /*--Specify gap parms map extending--*/
    float d_anchor_rate_origin = gap_parms.thd_gmsa_d_anchor_rate;
    gap_parms.direction = direction;
    gap_parms.thd_ctfas2_connect_danchor = 50;
    gap_parms.thd_ctfas2_connect_dy_dx = 150;
    gap_parms.f_gmsa_direction = direction;
    gap_parms.thd_cts_major_limit = 3;
    gap_parms.f_me_map_extend = 1;
    gap_parms.thd_gmsa_d_anchor_rate = 0.25;
    print_cord(gap_str, "me1");
    print_cord(gap_end, "me2");
    String <Dna5> & ref = seqs[get_cord_id(gap_str)];
    String<uint64_t> sp_tiles; 
    //mapInterval(ref, read, comstr, tiles_str, f1, f2, gap_str, gap_end, 0, 0, direction, gap_parms);
    extendIntervalOneSide(ref, read, comstr, tiles_str, f1, f2, gap_str, gap_end,
                  direction, gap_parms);
    //filter out tiles of large gaps
    g_print_tiles_(tiles_str, "me3");
    mapExtendResultFilter_(tiles_str, gap_str, gap_end, direction, gap_parms);
    if (!empty(tiles_str) && isClipTowardsRight(direction))
    {
        remove_tile_sgn_end(back(tiles_str));
    }
    reform_tiles(ref, read, comstr, tiles_str, tiles_end, sp_tiles, 
                 gap_str, gap_end, direction, gap_parms);

    gap_parms.f_me_map_extend = 0;
    gap_parms.thd_gmsa_d_anchor_rate = d_anchor_rate_origin;
    return 0;
}
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
               GapParms & gap_parms)
{
    /*--Specify gap parms map extending--*/
    gap_parms.thd_ctfas2_connect_danchor = 50;
    gap_parms.thd_ctfas2_connect_dy_dx = 150;
    gap_parms.thd_cts_major_limit = 3;
    gap_parms.f_me_map_extend = 1;
    int original_direction = gap_parms.direction;
    int original_f_rfts_clip = gap_parms.f_rfts_clip;
    int direction1 = g_map_rght, direction2 = g_map_left;
    gap_parms.f_rfts_clip = 0; //disable clip in when reform tiles.
    String<Dna5> & ref = seqs[get_cord_id(gap_str1)];
    String<uint64_t> sp_tiles1; 
    String<uint64_t> sp_tiles2; 

    extendsInterval(ref, read, comstr, tiles_str1, tiles_str2, f1, f2, gap_str1, gap_end1, gap_str2, gap_end2, gap_parms);
    //direction = 1 part
    gap_parms.direction = direction1;
    mapExtendResultFilter_(tiles_str1, gap_str1, gap_end1, direction1, gap_parms);
    if (!empty(tiles_str1))
    {
        remove_tile_sgn_end(back(tiles_str1));
    }
    reform_tiles(ref, read, comstr, tiles_str1, tiles_end1, sp_tiles1, 
                 gap_str1, gap_end1, direction1, gap_parms);
    //<<debug
    if (!empty(tiles_end1))
    {
    //    back(tiles_end1) = shift_tile(back(tiles_end1), -95, -95);
    }
    //>>debug
    //direction = -1 part
    gap_parms.direction = direction2; 
    mapExtendResultFilter_(tiles_str2, gap_str2, gap_end2, direction2, gap_parms);
    reform_tiles(ref, read, comstr, tiles_str2, tiles_end2, sp_tiles2, 
                 gap_str2, gap_end2, direction2, gap_parms);
    //restore gap_parms

    gap_parms.direction = original_direction;
    gap_parms.f_rfts_clip = original_f_rfts_clip;
    gap_parms.f_me_map_extend = 0;
    return 0;    
}
/*---------------  Map of generic type  ---------------*/
/*
 * Wrapper of calling reExtendChainOneSide to clip
 * @extend_lower_cord, @extend_upper_cord is the bound where the chain extended to,
   usually gap_str or gap_end
 * In [i_ptr_str, i_ptr_end], the function is called in the closed domain
 */
int reExtendClipOneSide(String<Dna5> & ref, 
                        String<Dna5> & read, 
                        String<Dna5> & comstr, 
                        String<uint64_t> & chain, 
                        uint64_t extend_lower_cord,
                        uint64_t extend_upper_cord,
                        int i_ptr_str,
                        int i_ptr_end,
                        int direction,
                        GapParms & gap_parms)
{
    if (empty(chain) || i_ptr_str < 0 || i_ptr_end < 0)
    {
        return 0;
    }
    int lower = 60, upper = 60;
    int shape_len = gap_parms.thd_etfas_shape_len; 
    int step1 = gap_parms.thd_etfas_step1;
    int step2 = gap_parms.thd_etfas_step2;
    if (isClipTowardsLeft(direction))
    {
        int dx = get_tile_x(chain[i_ptr_str]) - get_tile_x(extend_lower_cord);
        int dy = get_tile_strand(chain[i_ptr_str]) ^ get_tile_strand(extend_lower_cord) ?
            get_tile_y(extend_upper_cord)  - length(read) + get_tile_y(chain[i_ptr_str]) :
            get_tile_y(chain[i_ptr_str]) - get_tile_y(extend_lower_cord);
        lower = std::min({dx, dy, lower});
    }
    else if (isClipTowardsRight(direction))
    {
        int dx = get_tile_x(extend_upper_cord) - 1 - get_tile_x(chain[i_ptr_end]);
        int dy = get_tile_strand(chain[i_ptr_end]) ^ get_tile_strand(extend_upper_cord) ?
            length(read) - 1 - get_tile_y(chain[i_ptr_end]) - get_tile_y(extend_lower_cord) :
            get_tile_y(extend_upper_cord) - get_tile_y(chain[i_ptr_end]);
        upper = std::min({dx, dy, upper});
    //dout << "luex" << lower << upper << direction << dx << dy  << "\n";   
    }
    return reExtendChainOneSide(ref, read, comstr, chain, i_ptr_str, i_ptr_end, lower, upper, 
                shape_len, step1, step2, direction, 
                &get_tile_x, &get_tile_y, &get_tile_strand, &set_tile_strand, &shift_tile, 
                &g_hs_anchor2Tile, 
                gap_parms);
}
int createTilesFromAnchors2_(String<Dna5> & ref,
                             String<Dna5> & read,
                             String<Dna5> & comstr,
                             String<uint64_t> & anchors, 
                             String<uint64_t> & tiles_str,
                             String<uint64_t> & tiles_end,
                             StringSet<FeaturesDynamic> & f1,
                             StringSet<FeaturesDynamic> & f2,
                             uint64_t gap_str,
                             uint64_t gap_end,
                             uint64_t read_len,
                             int direction,
                             GapParms & gap_parms)
{

    String<uint64_t> tmp_tiles;
    g_CreateChainsFromAnchors_(anchors, tmp_tiles, gap_str, gap_end, read_len, gap_parms);
    int pre_i = 0;
    g_print_tile(gap_str, "ctfa2");
    g_print_tile(gap_end, "ctfa2");
    g_print_tiles_(tmp_tiles, "ctfa2s");
    for (int i = 0; i < length(tmp_tiles); i++)
    {
        if (is_tile_end(tmp_tiles[i]))
        {
            uint64_t head_tile = tmp_tiles[pre_i];
            uint64_t tail_tile = tmp_tiles[i];
            std::cout << "cx0 " << pre_i << " " << i << "\n";
            i += reExtendClipOneSide(ref, read, comstr, tmp_tiles, gap_str, gap_end,
                                     pre_i, i, -1, gap_parms);
            g_print_tiles_(tmp_tiles, "ctfass1");
            std::cout << "cx1 " << pre_i << " " << i << "\n";
            g_print_tiles_(tmp_tiles, "ctfass2");
            i += reExtendClipOneSide(ref, read, comstr, tmp_tiles, gap_str, gap_end, 
                                     pre_i, i, 1, gap_parms);
            std::cout << "cx2 " << pre_i << " " << i << "\n";
            g_print_tiles_(tmp_tiles, "ctfass3");
            if (!(empty(tmp_tiles) || pre_i < 0 || i < 0))
            {
                copy_tile_sgn(head_tile, tmp_tiles[pre_i]);
                copy_tile_sgn(tail_tile, tmp_tiles[i]);

                g_CreateTilesFromChains_(tmp_tiles, tiles_str, tiles_end, f1, f2, 
                    gap_str, gap_end, pre_i, i + 1, &get_tile_x, &get_tile_y, 
                    &get_tile_strand, gap_parms);
            }
            pre_i = i + 1;
        }
        else if (i < length(tmp_tiles) - 1 && 
            get_tile_strand(tmp_tiles[i] ^ tmp_tiles[i + 1]))
        {

            int len = length(tiles_str);
            uint64_t head_tile = tmp_tiles[pre_i];
            uint64_t tail_tile = tmp_tiles[i];
            dout << "c33 " << i << length(tmp_tiles) << pre_i << "\n";
            //<<debug
            //if (!get_tile_strand(tmp_tiles[pre_i])) 
            //{
            //>>debug
            std::cout << "cx3 " << pre_i << i << "\n";
            i += reExtendClipOneSide(ref, read, comstr, tmp_tiles, gap_str, gap_end,
                                     pre_i, i, -1, gap_parms);
            std::cout << "cx4 " << pre_i << i << "\n";
            i += reExtendClipOneSide(ref, read, comstr, tmp_tiles, gap_str, gap_end, 
                                     pre_i, i, 1, gap_parms);
            std::cout << "cx5 " << pre_i << i << "\n";
           g_print_tiles_(tmp_tiles, "ctfa23");
            //}
            if (!(empty(tmp_tiles) || pre_i < 0 || i <0))
            {
                copy_tile_sgn(head_tile, tmp_tiles[pre_i]);
                copy_tile_sgn(tail_tile, tmp_tiles[i]);
                g_CreateTilesFromChains_(tmp_tiles, tiles_str, tiles_end, f1, f2, 
                    gap_str, gap_end, pre_i, i + 1, &get_tile_x, &get_tile_y, 
                    &get_tile_strand, gap_parms);
                if (len != length(tiles_str))
                {
                    remove_tile_sgn_end(back(tiles_str));
                    remove_tile_sgn_end(back(tiles_end));
                }
                dout << "ctfa25 " << pre_i << i << "\n";
           g_print_tiles_(tiles_str, "ctfa241");
           g_print_tiles_(tiles_end, "ctfa242");
            }
        //}
            pre_i = i + 1;    
            //<<debug
            //return 0;
            //>>debug 
        }
    }    
    //set_tile_end(back(tiles_str));
    //set_tile_end(back(tiles_end));
    g_print_tiles_(tiles_str, "ctfa3s1");
    g_print_tiles_(tiles_end, "ctfa3s2");
    return 0;
}
/**
 * ~Methods function of mapGAnchor2_
 * Map gaps of [gap_str, gap_end)
 * Cluster anchors and trim tiles for sv
 * Change thd_tile_size for differe size of window in the apx mapping
 */  
 int mapTilesFromAnchors (String<Dna5> & ref,
                     String<Dna5> & read,
                     String<Dna5> & comstr,
                     String<uint64_t> & anchors, 
                     String<uint64_t> & tiles_str, 
                     String<uint64_t> & tiles_end, 
                     StringSet<FeaturesDynamic> & f1,
                     StringSet<FeaturesDynamic> & f2,
                     uint64_t gap_str,
                     uint64_t gap_end,
                     uint64_t revscomp_const,
                     int direction,
                     GapParms & gap_parms)
{
    //method 1: sort to chain anchors o(nlg(n))
    //createTilesFromAnchors1_(anchor, tiles, f1, f2, gap_str, gap_end, anchor_end, thd_tile_size, thd_err_rate, thd_pattern_in_window, thd_anchor_density, thd_min_segment, gap_parms);

    //method 2: dp to chain anchors o(mn) ~ o(n)
    createTilesFromAnchors2_(ref, read, comstr, anchors, tiles_str, tiles_end, f1, f2, gap_str, 
        gap_end, revscomp_const, direction, gap_parms);
    g_print_tiles_(tiles_str, "mtfas1");
    g_print_tiles_(tiles_end, "mtfas2");
    //trimTiles(tiles_str, f1, f2, gap_str, gap_end, revscomp_const, direction, gap_parms);
    g_print_tiles_(tiles_str, "mtfas3");
    g_print_tiles_(tiles_end, "mtfas4");

    return 0;
}
/**
 * Map interval [@gap_str, @gap_end) to extend the
   mapped region as long as possible.
 * @gap_str and @gap_end should have the same strand
 */
int mapInterval(String<Dna5> & seq1, //genome
                 String<Dna5> & seq2, //read
                 String<Dna5> & comstr,
                 String<uint64_t> & tiles_str,    //results
                 String<uint64_t> & tiles_end,    //results
                 StringSet<FeaturesDynamic > & f1,  
                 StringSet<FeaturesDynamic > & f2,
                 uint64_t gap_str,
                 uint64_t gap_end, 
                 int64_t anchor_lower,
                 int64_t anchor_upper,
                 int direction,
                 GapParms & gap_parms) // extern parm
{
    if (get_cord_strand (gap_str ^ gap_end))
    {
        return 1;
    }
    int shape_len = 9; 
    int step1 = 5; //seq1 pattern step
    int step2 = 1; //seq2...
    String<uint64_t> g_hs;
    String<uint64_t> g_hs_anchors;
    reserve(g_hs, 2048);
    reserve(g_hs_anchors, 2048);
    g_print_tile(gap_str, "mi1");
    g_print_tile(gap_end, "mi2");
    g_stream_(seq1, seq2, g_hs, gap_str, gap_end, shape_len, step1, step2, gap_parms);
    g_create_anchors_(g_hs, g_hs_anchors, shape_len, direction, anchor_lower, anchor_upper, length(seq2) - 1, gap_str, gap_end, gap_parms);
    mapTilesFromAnchors (seq1, seq2, comstr, g_hs_anchors, tiles_str, tiles_end, f1, f2, gap_str, gap_end, length(seq2) - 1, direction, gap_parms);
    return 0;
}
/*
 * Map generic gap of [@gap_str, gap_end)
 * Output one best @tiles_str1 and @tiles_end1 
 * Map direction = 0 (closed)
 */
int mapGeneric(StringSet<String<Dna5> > & seqs, 
               String<Dna5> & read, 
               String<Dna5> & comstr,
               StringSet<FeaturesDynamic > & f1, 
               StringSet<FeaturesDynamic > & f2,
               String<uint64_t> & tiles_str, 
               String<uint64_t> & tiles_end,  
               uint64_t gap_str, 
               uint64_t gap_end, 
               GapParms & gap_parms)
{
    int t_direction = 0;
    uint64_t thd_gather_block_gap_size = 100; //warn::not the thd_gap_size
    String<uint64_t> sp_tiles_inv;
    int f_rfts_clip = gap_parms.f_rfts_clip;
    gap_parms.f_rfts_clip = 0; // createTileFromaAnchors2_ in mapInterval alredy clipped chain.
    g_print_tile(gap_str, "mg221");
    g_print_tile(gap_end, "mg222");
    mapInterval(seqs[get_tile_id(gap_str)], read, comstr, tiles_str, tiles_end, f1, f2,
                        gap_str, gap_end, LLMIN, LLMAX, t_direction, gap_parms);  
    //chainTiles(tiles_str1, length(read), thd_gather_block_gap_size, gap_parms);
    g_print_tiles_(tiles_str, "mg13");
    reform_tiles(seqs[get_tile_id(gap_str)], read, comstr, tiles_str, tiles_end, 
        sp_tiles_inv, gap_str, gap_end, t_direction, gap_parms);
    gap_parms.f_rfts_clip = f_rfts_clip;
    g_print_tiles_(tiles_str, "mg11");
    g_print_tiles_(tiles_end, "mg12");
    return 0;
}