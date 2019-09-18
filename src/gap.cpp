#include <utility> 
#include "base.h"
#include "shape_extend.h"
#include "cords.h"
#include "gap.h"

using std::endl;
/*=============================================
=               Glaobal Utilities             =
=============================================*/
/*
struct BitVector() //uint64_t bit operations extractor
{
    String<unsigned> bits;
    String<uint64_t> masks;
    //setValue (int)
    BitVector(String<int> & seg_lens)
}
BitVector::BitVector(String<int> & seg_lens)
{
    uint s = 0;
    for (int i = 0; i < length(seg_lens); i++)
    {
        s += seg_lens[i];
        if (s < 64)
        {
            appendValut(bits, s);
            appendValue(masks, (1ULL << seg_lens[i]) - 1);
        }
    }
}
*/

//NOTE::clip & map direction:: towards left < 0, right > 0, both 0
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

int const g_map_left = -1;
int const g_map_closed = 0;
int const g_map_rght = 1;

/*
 * NOTE! the following parameters are correlated.
 * Change them carefully 
 */
int const g_shape_len = 16;
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
    if (i << 1 < length(clips))
    {
        return clips[i << 1];
    }
    else
    {
        return 0;
    }
}

uint64_t getClipEnd(String<uint64_t> & clips, int i)
{
    if (i << 1 < length(clips) - 1)
    {
        return clips[(i << 1) + 1];
    }
    else
    {
        return 0;
    }
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
    return clip_direction < 0; 
}

int isClipTowardsRight (int clip_direction)
{
    return clip_direction > 0;
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
                         uint64_t bit = _defaultTileBase.strandBit
                        )
    {
        return (val >> bit) & 1;
    }
     uint64_t slide (uint64_t val, 
                           uint64_t x, 
                           uint64_t y,
                           uint64_t bit = _defaultTileBase.yBitLen)
    {
        return val + (x << bit) + y;
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
    bool isTileEnd(uint64_t & val, uint64_t bit = _defaultTileBase.sgnBit_end)
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
    void removeTileSgnStart(uint64_t & val,
        uint64_t bit = ~_defaultTileBase.sgnBit_str)
    {
        val &= bit;
    }
    void removeTileSgnEnd(uint64_t & val,
        uint64_t bit = ~_defaultTileBase.sgnBit_end)
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
inline uint16_t get_tile_strand (uint64_t val)
{
    return _defaultTile.getStrand(val);
}
inline void set_tile_end (uint64_t & val)
{
    _defaultTile.setTileEnd(val);
}
inline void set_tile_start (uint64_t & val)
{
    _defaultTile.setTileStart(val);
}
inline void remove_tile_sgn (uint64_t & val)
{
    _defaultTile.removeTileSgn(val);
}
inline void copy_tile_sgn (uint64_t tile1, uint64_t & tile2)
{
    _defaultTile.copyTileSgn(tile1, tile2);
}
inline void remove_tile_sgn_start(uint64_t &val)
{
    _defaultTile.removeTileSgnStart(val);
}
inline void remove_tile_sgn_end(uint64_t & val)
{
    _defaultTile.removeTileSgnEnd(val);
}
inline bool is_tile_start(uint64_t val)
{
    return _defaultTile.isTileStart(val);
}
inline bool is_tile_end(uint64_t val)
{
    return _defaultTile.isTileEnd(val);
}
inline bool is_tile_body(uint64_t val)
{
    return _defaultTile.isTileBody(val);
}
inline uint64_t shift_tile(uint64_t const & val, int64_t x, int64_t y)
{
    return shift_cord (val, x, y);
}
inline uint64_t get_tile_x (uint64_t val)
{
    return get_cord_x(val);
}
inline uint64_t get_tile_y (uint64_t val)
{
    return get_cord_y(val);
}
inline uint64_t get_tile_id(uint64_t val)
{
    return get_cord_id(val);
}
inline uint64_t create_tile (uint64_t id, uint64_t cordx, uint64_t cordy, uint64_t strand)
{
    return create_cord(id, cordx, cordy, strand);
}

/**
 * debug utils
 */
void g_print_tile (uint64_t tile, CharString str)
{
    std::cout << str << " " 
              << get_cord_id(tile) << " " 
              << get_cord_strand(tile) << " " 
              << get_cord_x(tile) << " "
              << get_cord_y(tile) << " " 
              << get_cord_x(tile) - get_cord_y (tile) << "\n";    
}
void g_print_tiles_(String<uint64_t> & tiles, CharString str = "print_tiles")
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

uint64_t g_hs_anchor_getAnchor (uint64_t anchor) // return @anchor := strand + anchorx
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
//debug util
void print_g_hs_anchor(String<uint64_t> & anchors,
                       int start,
                       int end,
                       CharString header = "print_g_hs_anchor")
{
    for (int i = start; i < end; i++) 
    {
    std::cout << header <<" " << i << " " 
              << g_hs_anchor_getX(anchors[i]) << " "
              << g_hs_anchor_getY(anchors[i]) << " "
              << g_hs_anchor_get_strand(anchors[i]) << " "
              << "\n";
    }
}
///g_hs: N/A[1]|xval[30]|type[2]strand[1]|coordinate[30]
///type=0: from genome, type=1: from read
const uint64_t g_hs_bit1 = 30;
const uint64_t g_hs_bit2 = 31;
const uint64_t g_hs_bit3 = 33;
const uint64_t g_hs_mask2 = (1ULL << 30) - 1;
const uint64_t g_hs_mask3 = (1ULL << 32) - 1;

void g_hs_setGhs_(uint64_t & val, 
                  uint64_t xval, 
                  uint64_t type, 
                  uint64_t strand, 
                  uint64_t coord)
{
    val = (xval << 33) + (type<< 31) + (strand << 30) + coord;
}

int64_t g_hs_getCord(uint64_t & val)
{
    return int64_t(val & g_hs_mask2);
}

void g_hs_setAnchor_(uint64_t & val, 
                     uint64_t const & hs1, /*genome*/
                     uint64_t const & hs2, /*read*/
                     uint64_t revscomp_const)
{
    uint64_t strand = ((hs1 ^ hs2) >> 30 ) & 1;
    uint64_t x = revscomp_const * strand - _nStrand(strand) * (hs2 & g_hs_mask2); 
    val = (((hs1 + g_hs_anchor_zero - x) & (g_hs_mask2))<< 20) + 
          x + (strand << g_hs_anchor_bit2);
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

uint64_t g_hs_anchor2Tile (uint64_t & anchor)
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
    if (get_tile_strand (tile1 ^ tile2))
    {
        return get_tile_y(tile2) - readlen + 1 + get_tile_y(tile1);
    }
    else
    {
        return (int64_t)(_defaultTile.getY(tile2)) - (int64_t)(_defaultTile.getY(tile1));
    }
}

/*
 * shortcut to return if 2 cords have different anchors 
 * @cord1d and @cord2 are required to have same strand
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
    int64_t da = dy - dx;
    int64_t dmax = std::max(std::abs(dx), std::abs(dy));
    return std::abs(da) > int64_t(std::max(thd_dxy_min, dmax) * thd_da_zero) && 
           da * greater >= 0;
}

/**
 * debug utility
 */
int _fscore (String<Dna5> & seq, String<Dna5> read, uint64_t gstart, uint64_t rstart)
{
    String<int> gfs;
    String<int> rfs;
    resize (gfs, 4);
    resize (rfs, 4);
    for (int i = 0; i < 4; i++)
    {
        gfs[i] = 0;
        rfs[i] = 0;
    }
    int fsc = 0;
    int script_len = 32;
    int delta = 100;
    for (delta; delta < 10000; delta+=32)
    {
        for (int d2 = 16; d2 < 5000; d2 += 10)
        {
            fsc = 0;
            for (int i = 0; i < window_size / script_len ; i++)
            {
                for (int j = 0; j < script_len; j++)
                {
                    gfs[ordValue(*(begin(seq) + gstart + i * script_len + j + delta))]++;
                    rfs[ordValue(*(begin(read) + rstart + i * script_len + j))]++;   
                }
                for (int j = 0; j < 4; j++)
                {
                    fsc += std::abs(gfs[j] - rfs[j]);
                    gfs[j] = rfs[j] = 0;
                }   
            
            }   
        }
        
    }
    return fsc;
}

/**
 * debug utility 
 */
int check_tiles_(String<uint64_t> & tiles, uint64_t g_start, uint64_t g_end)
{
    for (uint64_t k = 1; k < length(tiles); k++)
    {
        if (_defaultTile.getX(tiles[k] - tiles[k - 1]) > window_size)
        {
            return 1;
        }
    }
    
    if (length(tiles) > 0)
    {
        if (_defaultTile.getX(tiles[0]) -  g_start > window_size )//|| g_end - _defaultTile.getX(tiles[length(tiles) - 1]) > window_size)
        {
            return 1;
        }   
    }
    else 
    {
        if (g_end - g_start > window_size)
        {
            return 1;
        }
    }
    return 0;
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
                   uint64_t start, 
                   uint64_t end, 
                   int g_hs_start, 
                   int shape_len,
                   int step,  
                   uint64_t type)
{
    LShape shape(shape_len);
    hashInit(shape, begin(seq) + start);
    int count = 0; 
    int i = 0; 
    uint64_t val = 0;
    for (uint64_t k = start; k < end; k++)
    {
        val = hashNextV(shape, begin(seq) + k);
        if (++count == step)  //collecting every step bases
        {
            //TODO: k - getT(shape)
            //<<debug
            g_hs_setGhs_(g_hs[g_hs_start + i++], val, type, shape.strand, k);
            count = 0;
        }
    }
    return g_hs_start + i;
}

/**
 * Stream block of @g_hs specified by @p1, @p2, @k and 
   filter anchor that is in [@acnhor_lower, anchor_upper).
 */
int g_mapHs_setAnchors_ (String<uint64_t> & g_hs, 
                         String<uint64_t> & g_anchor,
                         int p1, 
                         int p2, 
                         int k, 
                         int g_anchor_end,
                         uint64_t revscomp_const,
                         int64_t anchor_lower,
                         int64_t anchor_upper)
{
    unsigned n = 0;
    for (int i = p1; i < p2; i++) 
    {
        for (int j = p2; j < k; j++) 
        {
            uint64_t tmp_anchor;
            g_hs_setAnchor_(tmp_anchor, g_hs[i], g_hs[j], revscomp_const);
            int64_t tmp = g_hs_anchor_getAnchor(tmp_anchor);
            if (tmp < anchor_upper && tmp >= anchor_lower)
            {
                g_anchor[g_anchor_end + n++] = tmp_anchor;
            }
        }   
    }
    return g_anchor_end + n;
}


/*----------  Section of MapAnchor2_: function, parm and wrapper  ----------*/

unsigned _get_tile_f_ (uint64_t & tile,
                       StringSet<FeaturesDynamic > & f1,
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
 */
unsigned _get_tile_f_tri_ (uint64_t & tile,
                           uint64_t & new_tile,
                           StringSet<FeaturesDynamic > & f1,
                           StringSet<FeaturesDynamic > & f2, 
                           unsigned thd_accept_score,
                           int thD_tile_size)
{
    int shift = thD_tile_size / 4;
    unsigned thd_abort_score = UMAX; // make sure thd_abort_score > thd_accept_score
    unsigned fscore =  _get_tile_f_ (tile, f1, f2) ;

    uint64_t x = get_tile_x(tile);
    uint64_t y = get_tile_y(tile);
    if (fscore < thd_accept_score)
    {
        new_tile = tile;
        return fscore;
    }
    else 
    {
        uint64_t tile_l = shift_tile(tile, -shift, -shift);
        uint64_t tile_r = shift_tile(tile, shift, shift);
        fscore = std::min(_get_tile_f_(tile_l, f1, f2), fscore);
        if (fscore < thd_accept_score)
        {
            new_tile = tile_l;
            return fscore;
        }
        else
        {
            new_tile = tile_r;
            fscore = std::min(_get_tile_f_(tile_r, f1, f2), fscore);
            if (fscore < thd_accept_score)
            {
                return fscore;
            }
            else
            {
                new_tile = tile;
                return UMAX;
            }
        }
    }
    return UMAX;
}

int apxCreateTilesFromAnchors_ (String<uint64_t> & anchor, 
                                String<uint64_t> & tiles, 
                                StringSet<FeaturesDynamic> & f1,
                                StringSet<FeaturesDynamic> & f2,
                                uint64_t gap_str, 
                                int prek, int k, 
                                int const & thD_tile_size,
                                int const & thd_pattern_in_window)
{
    uint thd_fscore = getWindowThreshold(f1);
    uint64_t new_tile = 0;
    int64_t tmp_shift = thD_tile_size / 2;
    int64_t prex = g_hs_anchor_getX(anchor[prek]);
    int64_t prey = g_hs_anchor_getY(anchor[prek]);
    int64_t centroid_x = 0;
    int64_t centroid_y = 0;
    int kcount = 0;
    int prej = prek;
    for (int j = prek; j < k; j++)
    {
        uint64_t x = g_hs_anchor_getX(anchor[j]);
        uint64_t y = g_hs_anchor_getY(anchor[j]);
        if ((x > prex + thD_tile_size || y > prey + thD_tile_size) || j == k - 1)
        {
            if (j == prej) //only when k == prek + 1
            {
                centroid_x = x;
                centroid_y = y;
            }   
            else 
            {
                centroid_x /= (j - prej);
                centroid_y /= (j - prej);
            }

            tmp_shift = std::min({int64_t(thD_tile_size) >> 1, centroid_x, centroid_y});
            uint64_t tmp_tile = create_tile(get_cord_id(gap_str), 
                                            centroid_x - tmp_shift,
                                            centroid_y - tmp_shift,
                                            g_hs_anchor_get_strand(anchor[prek]));
            unsigned score = 
            _get_tile_f_tri_(tmp_tile, new_tile, f1, f2, thd_fscore, thD_tile_size);
            if (kcount >= thd_pattern_in_window && score < thd_fscore)
            {
                if (empty (tiles) || is_tile_end(back(tiles)))
                {
                    set_tile_start(new_tile);
                    appendValue (tiles, new_tile);
                }
                else  
                {
                    if (get_tile_x (new_tile) > get_tile_x(back(tiles)) &&
                        get_tile_y (new_tile) > get_tile_y(back(tiles)))
                    {
                        appendValue (tiles, new_tile);
                    }
                }
            }
            if (j != k - 1)
            {
                prex = g_hs_anchor_getX(anchor[j - 1]);
                prey = g_hs_anchor_getY(anchor[j - 1]);
                prej = j;
                centroid_x = g_hs_anchor_getX(anchor[j]);
                centroid_y = g_hs_anchor_getY(anchor[j]);
                kcount=1;
            }
        }
        else
        {
            centroid_x += g_hs_anchor_getX(anchor[j]);
            centroid_y += g_hs_anchor_getY(anchor[j]);
            kcount++;
        }
    }
    if (!empty(tiles))
    {
        set_tile_end(back(tiles)) ;
    }
    return 0;
}

int createTilesFromAnchors1_(String<uint64_t> & anchor, 
                             String<uint64_t> & tiles, 
                             StringSet<FeaturesDynamic> & f1,
                             StringSet<FeaturesDynamic> & f2,
                             uint64_t gap_str, 
                             uint64_t gap_end,
                             int anchor_end, 
                             int const & thD_tile_size,
                             float const & thD_err_rate,
                             int const & thd_pattern_in_window,
                             float const & thd_anchor_density,
                             int64_t const & thd_min_segment)
{
    int anchor_len = 0;
    std::sort (begin(anchor), begin(anchor) + anchor_end);
    anchor[anchor_end] = ~0;
    int prek = 0;
    for (int k = 0; k < anchor_end + 1; k++)
    {
        //TODO: handle thd_min_segment, anchor 
        int64_t d = std::abs((int64_t)g_hs_anchor_getY(anchor[k]) - (int64_t)g_hs_anchor_getY(anchor[prek]));
        if (g_hs_anchor_getAnchor(anchor[k]) - g_hs_anchor_getAnchor(anchor[prek]) > 
            thD_err_rate * std::max(thd_min_segment, d))
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
                apxCreateTilesFromAnchors_(anchor, tiles, f1, f2, gap_str, prek, k, thD_tile_size, thd_pattern_in_window);
            }
            prek = k;
            anchor_len = 0;
        }
        else
        {
            anchor_len++;
        }
    }
}


//ATTENTION::the Adjust @thd_abort_score if the function is changed
int getGapChainScore(uint64_t const & anchor1, uint64_t const & anchor2)
{
    int dy = g_hs_anchor_getY(anchor1) - g_hs_anchor_getY(anchor2);
    if (dy < 0)
    {
        return -10000;
    }

    int thd_min_dy = 50;
    int da = std::abs(int64_t(g_hs_anchor_getAnchor(anchor2) - g_hs_anchor_getAnchor(anchor1)));
    int d_err =  (100 * da) / std::max(dy, thd_min_dy); // 1/100 = 0.01
    //d_err
    if (d_err < 10)
    {
        d_err = 0;
    }
    else if (d_err < 15)
    {
        d_err = 10 + 2 * d_err ;
    }
    else 
    {
        d_err =  d_err * d_err / 10 + 40;
    }

    //d_y
    if (dy < 50)
    {
        dy = 0;
    }
    else if (dy < 100)
    {
        dy = dy - 30;
    }
    else
    {
        dy = dy * dy / 200 + 20;
    }
    return 100 - dy - d_err ;
}

int createTilesFromAnchors2_(String<uint64_t> & anchor, 
                             String<uint64_t> & tiles,
                             StringSet<FeaturesDynamic> & f1,
                             StringSet<FeaturesDynamic> & f2,
                             uint64_t & gap_str,
                             uint64_t & gap_end,
                             int anchor_end,
                             int const & thD_tile_size,
                             int const & thD_err_rate,
                             int const & thd_pattern_in_window)
{
    StringSet<String<uint64_t> > chains;
    ChainScoreMetric chn_score1(50, &getGapChainScore);
    std::sort(begin(anchor), begin(anchor) + anchor_end, 
        [](uint64_t & a, uint64_t & b){return g_hs_anchor_getX(a) > g_hs_anchor_getX(b);});
    createChainsFromAnchors (chains, anchor, anchor_end, chn_score1);
    for (auto & chain : chains)
    {
        apxCreateTilesFromAnchors_(chain, tiles, f1, f2, gap_str, 0, length(chain), thD_tile_size, thd_pattern_in_window);
    }
    return 0;
}

/**
 * ~Methods function of mapGAnchor2_
 * Map gaps of [gap_str, gap_end)
 * Cluster anchors and trim tiles for sv
 * Change thD_tile_size for differe size of window in the apx mapping
 */  
 int map_g_anchor2_ (String<uint64_t> & anchor, 
                   String<uint64_t> & tiles, 
                   StringSet<FeaturesDynamic> & f1,
                   StringSet<FeaturesDynamic> & f2,
                   uint64_t gap_str,
                   uint64_t gap_end,
                   int anchor_end, 
                   int revscomp_const,
                   int direction,
                   int const & thD_tile_size,
                   float const & thD_err_rate,
                   int const & thd_pattern_in_window,
                   int const & thd_overlap_size,
                   int const & thd_gap_size,
                   float const & thd_overlap_tile,
                   float const & thd_swap_tile,
                   float const & thd_anchor_density,
                   int64_t const & thd_min_segment)
{
    //dout << "sv2\n";
    int thd_abort_tiles = 10;

    //step 1. 
    //sort to chain anchors o(nlg(n))
    //createTilesFromAnchors1_(anchor, tiles, f1, f2, gap_str, gap_end, anchor_end, thD_tile_size, thD_err_rate, thd_pattern_in_window, thd_anchor_density, thd_min_segment);

    //dp to chain anchors o(mn) ~ o(n)
    createTilesFromAnchors2_(anchor, tiles, f1, f2, gap_str, gap_end, anchor_end, thD_tile_size, thD_err_rate, thd_pattern_in_window);

    //step 2. merge check: if different segments in tiles can be megered; if cant then do nothing
    String<uint64_t> tmp_tiles = tiles;
    std::sort (begin(tmp_tiles), end(tmp_tiles),
    [](uint64_t & s1, uint64_t & s2)
    {  
        return get_tile_x(s1) < get_tile_x(s2);
    });
    int merge_flag = 1;
    for (int i = 1; i < length(tmp_tiles); ++i)
    {
        if (_defaultTile.getY(tmp_tiles[i]) < _defaultTile.getY(tmp_tiles[i - 1]) &&
            !_defaultTile.getStrand(tmp_tiles[i] ^ tmp_tiles[i - 1]))
        {
            merge_flag = 0;
            break;
        }
    }
    if (merge_flag)
    {
        for (int i = 0; i < length(tmp_tiles); ++i)
        {
            tiles[i] = tmp_tiles[i];
            remove_tile_sgn(tiles[i]); //remove start and end sign
        }
        if (!empty(tiles))
        {
            set_tile_start(tiles[0]);
            set_tile_end(back(tiles));
        }
    }
/**
 * step 3. extend patch
 * extend window if there are gaps between tiles until the 
   coordinates x1 - x2 < window_size or the gap can't be extend any more
 * ATTENTION: This methods takes no account of the relation between y1 and y2.
 */

    uint64_t cord_str = gap_str;
    int64_t shift_x = std::min(int64_t(get_cord_x(gap_end) - get_cord_x(gap_str)), int64_t(thD_tile_size));
    int64_t shift_y = std::min(int64_t(get_cord_y(gap_end) - get_cord_y(gap_str)), int64_t(thD_tile_size));
    uint64_t cord_end = shift_cord(gap_end, -shift_x, -shift_y);
    //uint64_t cord_end = gap_end;
    if (empty(tiles))
    {
        extendPatch(f1, f2, tiles, 0, cord_str, cord_end, revscomp_const, thd_overlap_size, thd_gap_size);
        if (!empty(tiles))
        {
            set_tile_start(tiles[0]);
            set_tile_end(back(tiles));
        }
        return 0;
    } 
    for (int i = 0; i < length(tiles); i++)
    {
        if (is_tile_start(tiles[i]))
        {
            int new_num = extendPatch(f1, f2, tiles, i, cord_str, tiles[i], revscomp_const, thd_overlap_size, thd_gap_size);
            if (new_num)
            {
                set_tile_start(tiles[i]);
                i += new_num;
                remove_tile_sgn_start(tiles[i]);
            }
        }
        if (is_tile_end(tiles[i]))
        {
            int new_num = extendPatch(f1, f2, tiles, i + 1, tiles[i], cord_end, revscomp_const, thd_overlap_size, thd_gap_size);   
            if (new_num)
            {
                remove_tile_sgn_end(tiles[i]);
                //g_print_tile(tiles[i], "esend");
                //g_print_tile(tiles[i + new_num], "esend2");
                i += new_num;
                set_tile_end(tiles[i]);
            }
        }
        if (i >= 1 && !is_tile_end (tiles[i - 1]) && !is_tile_start(tiles[i]))
        {
            //g_print_tile(tiles[i], "esmid");
            //g_print_tile(tiles[i + new_num], "esmid2");
            i += extendPatch(f1, f2, tiles, i, tiles[i - 1], tiles[i], revscomp_const, thd_overlap_size, thd_gap_size);   
        }
    }
    //Remove out of bound tiles.
    int64_t x_str = get_tile_x(gap_str);
    int64_t y_str = get_tile_y(gap_str);
    int64_t x_end = get_cord_x(gap_end);
    int64_t y_end = get_cord_y(gap_end);
    int di = 0;
    for (int i = 0; i < length(tiles); i++)
    {
        uint64_t x_t = get_tile_x(tiles[i]);
        uint64_t y_t = get_tile_strand(tiles[i] ^ gap_str) ?
                       revscomp_const - 1 - get_tile_y(tiles[i]) - thD_tile_size :
                       get_tile_y(tiles[i]);
        if (x_t < x_str || x_t + thD_tile_size > x_end || 
            y_t < y_str || y_t + thD_tile_size > y_end) //out of bound of [gap_str, gap_end)
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
                    //set_cord_end (tiles[i - di - 1]);
                }
            }
            else
            {
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

//!!TODO::[Important] tune the parm, especially the thD_err_rate [ing] 
//Parms of function map_g_Anchor2_
struct MapGAnchor2Parm_
{
    int thd_pattern_in_window;
    int thd_overlap_size;
    int thd_gap_size;
    float thd_overlap_tile;
    float thd_swap_tile;
    float thd_anchor_density;
    int64_t thd_min_segment;

    MapGAnchor2Parm_ (int thD_tile_size, float thD_err_rate) 
    {
        thd_pattern_in_window = 1;
        thd_overlap_size = 170;
        thd_gap_size = 180;
        thd_overlap_tile = thD_tile_size * 0.4;
        thd_swap_tile = thD_tile_size * 0.05;
        thd_anchor_density = 0.03;
        thd_min_segment = 100;
    }
};

//Function wrapper
int map_g_anchor2_ (String<uint64_t> & anchor, 
                    String<uint64_t> & tiles, 
                    StringSet<FeaturesDynamic> & f1,
                    StringSet<FeaturesDynamic> & f2,
                    uint64_t gap_str,
                    uint64_t gap_end,
                    int anchor_end, 
                    int revscomp_const,
                    int direction,
                    int thD_tile_size, //global parm
                    float thD_err_rate,  //global parm
                    MapGAnchor2Parm_ const & parm)
{
    return map_g_anchor2_ (anchor, tiles, f1, f2,
                           gap_str, gap_end, anchor_end, 
                           revscomp_const, direction,
                           thD_tile_size,
                           thD_err_rate,
                           parm.thd_pattern_in_window,
                           parm.thd_overlap_size,
                           parm.thd_gap_size,
                           parm.thd_overlap_tile,
                           parm.thd_swap_tile,
                           parm.thd_anchor_density,
                           parm.thd_min_segment
                           );
}
/**
 * Main Function Parm Wrapper 
 */
struct MapGAnchorParm
{
    MapGAnchor2Parm_ mapGAnchor2parm_;
    MapGAnchorParm(int thD_tile_size, float thD_err_rate):
        mapGAnchor2parm_(thD_tile_size, thD_err_rate)
    {};
};
/**
 * Main Function MapAnchor_ 
 */
int map_g_anchor (String<uint64_t> & anchor, 
                  String<uint64_t> & tiles, 
                  StringSet<FeaturesDynamic > & f1,
                  StringSet<FeaturesDynamic > & f2,
                  uint64_t cord_str,
                  uint64_t cord_end,
                  int anchor_end, 
                  int revscomp_const,
                  int direction,
                  int thD_tile_size,
                  float thD_err_rate,
                  MapGAnchorParm const & parm
                  )
{
    return map_g_anchor2_(anchor, tiles, f1, f2, cord_str, cord_end, anchor_end, revscomp_const, direction, thD_tile_size, thD_err_rate, parm.mapGAnchor2parm_);
}

int g_create_anchors_ (String<uint64_t> & g_hs,
                       String<uint64_t> & g_hs_anchor,
                       int & g_hs_end,
                       int & g_hs_anchor_end,
                       int shape_len, 
                       int64_t anchor_lower,
                       int64_t anchor_upper,
                       uint64_t rvcp_const)
{
    uint64_t mask = (1ULL << (2 * shape_len + g_hs_bit3)) - 1;
    std::sort (begin(g_hs), begin(g_hs) + g_hs_end, [mask](uint64_t & a, uint64_t & b){return (a & mask) < (b & mask);});
    int p1 = 0, p2 = 0;
    for (int k = 1; k < g_hs_end; k++)
    {    
        switch (g_hs_getXT((g_hs[k] ^ g_hs[k - 1]) & mask))
        {
            case 0:       //x1 = x2 both from genome or read
                break;
            case 1:       //x1 = x2 one from genome the other from read
                p2 = k;
                break;
            default:      //anchor current block before process next block 
                g_hs_anchor_end = g_mapHs_setAnchors_(g_hs, g_hs_anchor, p1, p2, k, g_hs_anchor_end, rvcp_const, anchor_lower, anchor_upper);
                p1 = k;
                p2 = k; 
        }
    }
    return 0;
}

/**
 * Map interval [@gap_str, @gap_end) to create a chain of tiles to extend the
   mapping area as long as possible.
 * @gap_str and @gap_end should have the same strand
 */
int map_interval(String<Dna5> & seq1, //genome
                 String<Dna5> & seq2, //read
                 String<Dna5> & comstr,
                 String<uint64_t> & g_hs,
                 String<uint64_t> & g_hs_anchor,
                 String<uint64_t> & g_hs_tile,    //results
                 StringSet<FeaturesDynamic > & f1,  
                 StringSet<FeaturesDynamic > & f2,
                 uint64_t gap_str,
                 uint64_t gap_end, 
                 int64_t cord_lower,
                 int64_t cord_upper,
                 int direction,
                 int thD_tile_size,   //WARNING 192 not allowed to change.
                 float thD_err_rate,
                 MapGAnchorParm const & parm1 // extern parm
                )
{
    if (get_cord_strand (gap_str ^ gap_end))
    {
        return 1;
    }
    int g_hs_end = 0;
    int g_hs_anchor_end = 0;
    uint64_t rvcp_const = length(seq2) - 1;
    uint64_t gs_str = get_cord_x(gap_str);
    uint64_t gs_end = get_cord_x(gap_end);
    uint64_t gr_str = get_cord_y(gap_str);
    uint64_t gr_end = get_cord_y(gap_end);
    if (get_cord_strand(gap_str))
    {
        gr_str = rvcp_const - gr_str;
        gr_end = rvcp_const - gr_end;
        std::swap (gr_end, gr_str);
    }

    g_hs_end = g_mapHs_kmer_(seq1, g_hs, gs_str, gs_end, g_hs_end, 8, 10, 0);
    g_hs_end = g_mapHs_kmer_(seq2, g_hs, gr_str, gr_end, g_hs_end, 8, 1, 1);

    g_create_anchors_(g_hs, g_hs_anchor, g_hs_end, g_hs_anchor_end, 8, cord_lower, cord_upper, rvcp_const);

    int f = map_g_anchor (g_hs_anchor, g_hs_tile, f1, f2, 
                         gap_str, gap_end,
                         g_hs_anchor_end, 
                         rvcp_const,
                         direction,
                         thD_tile_size,
                         thD_err_rate,
                         parm1
                        );

    return f;

}

/**
 * Patch for duplication (only) where additional mapping in (,gap_str] or [gap_end,) is necessary.
 */
int try_dup(String<Dna5> & seq,
            String<Dna5> & read,
            String<Dna5> & comstr,
            StringSet<FeaturesDynamic > & f1,
            StringSet<FeaturesDynamic > & f2,
            String<uint64_t> & g_hs,
            String<uint64_t> & g_anchor,
            String<uint64_t> & tiles, //result
            uint64_t gap_str,
            uint64_t gap_end,
            float thD_err_rate,
            uint64_t thD_tile_size,
            MapGAnchorParm const & parm1)
            //float band_ratio,
            //int direction)
{
    int64_t dy;
    uint64_t tile1 = gap_str;
    uint64_t tile2 = shift_tile(gap_end, -thD_tile_size, -thD_tile_size);
    if (get_tile_strand(tile1 ^ tile2))
    {
        dy = length(read) - 1 - get_tile_y(tile2) - get_tile_y(tile1);
    }
    else
    {
        dy = get_tile_y(tile2) - get_tile_y(tile1);
    }
    if (dy < 0)
    {
        return 1;
    }
    int64_t shift_y = dy; 
    int64_t shift_x = dy * (1 + thD_err_rate);

    //extend tile1 towards right
    uint64_t try_str = shift_tile(tile1, thD_tile_size / 2, thD_tile_size / 2);
    uint64_t try_end = shift_tile(try_str, shift_x, shift_y);
    int direction = 1;
    map_interval(seq, read, comstr,
                 g_hs, g_anchor, tiles, f1, f2,
                 try_str, try_end,
                 LLMIN, LLMAX,             
                 direction,
                 thD_tile_size,
                 thD_err_rate,
                 parm1
                 );
    //extend tile2 towards left
    try_end = shift_tile(tile2, thD_tile_size / 2, thD_tile_size / 2);
    try_str = shift_tile(try_end, -shift_x, -shift_y);
    direction = -1;
    map_interval(seq, read, comstr,
                 g_hs, g_anchor, tiles, f1, f2,
                 try_str, try_end,
                 LLMIN, LLMAX,
                 direction,
                 thD_tile_size,
                 thD_err_rate,
                 parm1
                 );
    return 0;
}

int try_tiles_dup(String<Dna5> & seq,
                  String<Dna5> & read,
                  String<Dna5> & comstr,
                  StringSet<FeaturesDynamic > & f1,
                  StringSet<FeaturesDynamic > & f2,
                  String<uint64_t> & g_hs,
                  String<uint64_t> & g_anchor,
                  String<uint64_t> & tiles, //result
                  uint64_t tile1,
                  uint64_t tile2,
                  float thD_err_rate,
                  uint64_t thD_tile_size,
                  MapGAnchorParm parm1)
{
    //shortcut to call dup for tiles
    uint64_t gap_str = tile1;
    uint64_t gap_end = shift_tile(tile2, thD_tile_size, thD_tile_size);
    return 0;
    int res = try_dup (seq, read, comstr,
                       f1, f2,
                       g_hs, g_anchor, tiles, //result
                       gap_str,
                       gap_end,
                       thD_err_rate,
                       thD_tile_size, 
                       parm1);
    return res;
}

int g_mapHs_(String<Dna5> & seq1, //genome
             String<Dna5> & seq2, //read
             String<Dna5> & comstr,
             String<uint64_t> & g_hs,
             String<uint64_t> & g_hs_anchor,
             String<uint64_t> & tiles,    //results
             StringSet<FeaturesDynamic > & f1,  
             StringSet<FeaturesDynamic > & f2,
             uint64_t gap_str,
             uint64_t gap_end, 
             uint64_t anchor_lower,
             uint64_t anchor_upper,
             float thD_err_rate,
             int thD_tile_size,   //WARNING 192 not allowed to change.
             int64_t thd_dxy_min,
             int direction, //= g_map_closed
             MapGAnchorParm const & parm1
             )
{
    clear(g_hs);
    clear(g_hs_anchor);
    resize(g_hs, 1ULL << 20);
    resize(g_hs_anchor, 1ULL << 20);

    float thd_da_zero = thD_err_rate; 
    map_interval(seq1, seq2, comstr, g_hs, g_hs_anchor, tiles, f1, f2, gap_str, gap_end, anchor_lower, anchor_upper, direction, thD_tile_size, thD_err_rate, parm1);
}

int g_mapHs_(String<Dna5> & seq1, //genome
             String<Dna5> & seq2, //read
             String<Dna5> & comstr,
             String<uint64_t> & g_hs,
             String<uint64_t> & g_hs_anchor,
             String<uint64_t> & tiles,    //results
             StringSet<FeaturesDynamic > & f1,  
             StringSet<FeaturesDynamic > & f2,
             uint64_t gap_str,
             uint64_t gap_end, 
             float thD_err_rate,
             int thD_tile_size,   //WARNING 192 not allowed to change.
             int64_t thd_dxy_min,
             int direction, //= g_map_closed
             MapGAnchorParm const & parm1
             )
{
    g_mapHs_(seq1, seq2, comstr,
             g_hs, g_hs_anchor, tiles, 
             f1, f2,
             gap_str,gap_end, 
             LLMIN,LLMAX,
             thD_err_rate,
             thD_tile_size,   //WARNING 192 not allowed to change.
             thd_dxy_min,
             direction, //= g_map_closed
             parm1
             );
}

/*----------  Clip function  ----------*/

int64_t g_anchor_dx_(uint64_t val1, uint64_t val2)
{
    return int64_t (g_hs_anchor_getX(val2) - g_hs_anchor_getX(val1));
}
int64_t g_anchor_da_(uint64_t val1, uint64_t val2)
{
    return int64_t (g_hs_anchor_getAnchor(val2) - g_hs_anchor_getAnchor(val1));
}

/**
 * stream seq creating hs
 */
 int c_stream_(String<Dna5> & seq,
               String<uint64_t> & g_hs, 
               uint64_t start, 
               uint64_t end, 
               int g_hs_start, 
               int step,  
               uint64_t type)
{
    LShape shape(c_shape_len);
    hashInit_hs(shape, begin(seq) + start, 0);
    int count = 0; 
    int i = 0; 
    uint64_t val = 0;

    for (uint64_t k = start; k < end; k++)
    {
        val = hashNext_hs(shape, begin(seq) + k);
        if (++count == step)  //collecting every step bases
        {
            //TODO: k - getT(shape)
            g_hs_setGhs_(g_hs[g_hs_start + i++], val, type, 0, k);
            count = 0;
        }
    }
    return g_hs_start + i;
}

 void c_2Anchor_(uint64_t & val, uint64_t const & hs1, uint64_t const & hs2)
{
    ///hs1 genome, hs2 read
    uint64_t x = hs2 & g_hs_mask2; 
    val = (((hs1 - x) & (g_hs_mask2)) << g_hs_anchor_bit1) + x;
}

 void c_2GapAnchor_(uint64_t & val, uint64_t const & hs1, uint64_t const & hs2, uint64_t gap_type)
{
    ///hs1 genome, hs2 read
    uint64_t x = hs2 & g_hs_mask2; 
    val = (((hs1 - x) & (g_hs_mask2)) << g_hs_anchor_bit1) + x + (gap_type << g_hs_anchor_bit2);
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
 int c_create_anchor_block_ (String<uint64_t> & g_hs, 
                             String<uint64_t> & g_anchor,
                             int p1, 
                             int p2, 
                             int k, 
                             int g_anchor_end,
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
                c_2Anchor_(g_anchor[g_anchor_end++], g_hs[i], g_hs[j]);
            }
        }   
    }
    return g_anchor_end;
}

 int c_create_anchors_ (String<uint64_t> & g_hs, 
                        String<uint64_t> & g_anchor,
                        int g_hs_end,
                        int band_level,
                        int band_lower,
                        int64_t anchor_x,
                        int64_t anchor_y,
                        int64_t x_lower = 0,
                        int64_t x_upper = 0
                       ) 
{
    int p1 = 0, p2 = 0;
    int g_anchor_end = 0;
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
                g_anchor_end = c_create_anchor_block_(
                    g_hs, g_anchor, 
                    p1, p2, k, 
                    g_anchor_end, 
                    band_level, band_lower, 
                    anchor_x, anchor_y, 
                    x_lower, x_upper);
                p1 = k;
                p2 = k; 
        }
    }
    return g_anchor_end;
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
                         int anchor_end,
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
    anchor[anchor_end] = ~0;
    if (anchor_end < 1)
    {
        return direction > 0 ? clip_str : clip_end;
    }
    std::sort (begin(anchor), begin(anchor) + anchor_end);
    int i_str = 0;
    for (int i = 0; i < anchor_end - 1; i++)
    {
        if (g_hs_anchor_getY(anchor[i + 1] - anchor[i]) == 1 &&
            g_anchor_da_(anchor[i], anchor[i + 1]) == 0)
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
    i_str = anchor_end - 1 - ct_conts;
    uint64_t x = g_hs_anchor_getX(anchor[i_str]) - gs_str;
    uint64_t y = g_hs_anchor_getY(anchor[i_str]) - gr_str;
    anchor[it++] = (x << bit2) + ((ct_conts + 1) << g_hs_anchor_bit1) + y;

    if (it < 1)
    {
        return direction > 0 ? clip_str : clip_end;
    }
    //!NOTE::Value of anchor has been changed to := x|ct_conts|y
    int64_t y1 = 0;
    int64_t y2 = 0;
    int64_t x1 = 0;
    int64_t x2 = 0;
    int64_t x1_end = 0;
    int64_t x2_end = 0;
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
    if (dv == 0)  
    {
        return 1; // match
    }
    if (t1 + l1 - k + 1 == 0)
    {
        return 2; // mismatch
    }
    if (t2 + l1 - k == 0 && t2 != k && l1 != k)
    {
        return 3; //del 
    }
    if (t1 + l2 - k + 1 == 0)
    {
        return 4;  //ins
    }
    return 0;
 }

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
    unsigned shape_len =c_shape_len3;
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
}

struct ParmClipExtend
{
    int thd_min_scan_delta;
    int thd_error_level;
    int thd_gap_shape;
    int thd_merge_anchor;
    int thd_merge_drop;   
};
/**
 * Wrapper 
 * @clip_direction:
 *  1 towards right
 * -1 towards left 
 */
int c_clip_extend_(uint64_t & clip, 
                   String<uint64_t> & hs, 
                   String<uint64_t> & anchors,
                   String<Dna5> & seq1,
                   String<Dna5> & seq2,
                   uint64_t extend_str,
                   uint64_t extend_end,
                   ParmClipExtend & parm,
                   int clip_direction = -1 
    )
{
    c_clip_extend_(clip, 
                   hs, 
                   seq1,
                   seq2,
                   extend_str,
                   extend_end,
                   parm.thd_min_scan_delta,
                   parm.thd_error_level,
                   parm.thd_merge_anchor,
                   parm.thd_merge_drop,
                   clip_direction 
    ); 
}

/*
 *@clip_str, @clip_end required to have the equivalent strand
 */
uint64_t c_clip_(String<Dna5> & genome,  
                 String<Dna5> & read,
                 String<Dna5> & comstr,    //complement revers of read
                 uint64_t clip_str,
                 uint64_t clip_end,
                 String<uint64_t> & g_hs,
                 String<uint64_t> & g_anchor,
                 float thd_band_ratio,
                 int clip_direction = 1
                )
{
    CmpInt64 g_cmpll;
    uint64_t gs_str = get_tile_x(clip_str);
    uint64_t gr_str = get_tile_y(clip_str);
    uint64_t gs_end = get_tile_x(clip_end);
    uint64_t gr_end = get_tile_y(clip_end);
    uint64_t genomeId = get_tile_id (clip_str);
    uint64_t gr_strand = get_tile_strand(clip_str);

//step1. extend anchor
    String<Dna5> & seq1 = genome;
    String<Dna5> * p = (gr_strand)?(&comstr):(&read);
    String<Dna5> & seq2 = *p;
    int band = (gs_end - gs_str) * thd_band_ratio;
    ///clip scaffold
 
    int g_hs_end = c_stream_(seq1, g_hs, gs_str, gs_end, 0, 1, 0);
        g_hs_end = c_stream_(seq2, g_hs, gr_str, gr_end, g_hs_end, 1, 1);
    std::sort (begin(g_hs), begin(g_hs) + g_hs_end);
    int band_level = 3; //>>3 = /8 = * 12.5% error rate 
    uint64_t g_anchor_end = c_create_anchors_(g_hs, 
                                              g_anchor, 
                                              g_hs_end, 
                                              band_level, 
                                              band, 
                                              gs_str, 
                                              gr_str);
    int thd_merge1 = 3;
    int thd_merge1_lower = 10;
    int thd_merge2 = 20;
    int thd_width = 10;
    uint64_t clip = c_clip_anchors_(g_anchor, 
                                    clip_str,
                                    clip_end,
                                    g_anchor_end, 
                                    c_shape_len, 
                                    thd_merge1, 
                                    thd_merge1_lower, 
                                    thd_merge2, 
                                    clip_direction);    

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

    c_clip_extend_(clip, 
                   g_hs,
                   seq1,
                   seq2,
                   extend_str,
                   extend_end,
                   thd_scan_radius,
                   thd_error_level,
                   thd_merge_anchor,
                   thd_merge_drop,
                   clip_direction
                  );
    return clip;
}
/*
 * Clip exact breakpoints of the given tile at the front or end according to the clip direction.
 */
uint64_t clip_tile (String<Dna5> & seq1,
                    String<Dna5> & seq2,
                    String<Dna5> & comstr,
                    String<uint64_t> & g_hs,
                    String<uint64_t> & g_hs_anchor,
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
        clip = c_clip_ (seq1, seq2, comstr, clip_str, clip_end, g_hs, g_hs_anchor, thd_band_ratio, clip_direction); 
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
        clip = c_clip_ (seq1, seq2, comstr, clip_str, clip_end, g_hs, g_hs_anchor, thd_band_ratio, clip_direction); 
        //remove_tile_sgn(clip2);
    }
    remove_tile_sgn(clip);
    return clip;
}

/*----------  Reform tiles of gaps  ----------*/

int isTilesConsecutive_(uint64_t & tile1, uint64_t tile2, uint64_t thd_cord_gap)
{
    return isCordsConsecutive_(tile1, tile2, thd_cord_gap);
}

/*
 * @tile_str and tile_end refer to the start and end of the same tile rather than two different cords
 * @thD_tile_size is the regular size of tile, normally 192x192;
   However the size of the tile specified by the @tile_str and @tile_end is is allowed to be smaller than that.
 */
uint64_t reform_tile_ (String<Dna5> & seq1,
                       String<Dna5> & seq2,
                       String<Dna5> & comstr,
                       String<uint64_t> & g_hs,
                       String<uint64_t> & g_hs_anchor,
                       uint64_t & tile_str, 
                       uint64_t & tile_end, 
                       int sv_flag, 
                       int thD_tile_size
                       )
{
    int f_e = is_tile_end(tile_str);
    uint64_t clip = clip_tile (seq1, seq2, comstr, g_hs, g_hs_anchor, tile_str, sv_flag, thD_tile_size);
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

/**
 * Scan and Clip tile in @tiles_str from [@pos_str, @pos_end) if there exists gap;
 * if (length(tiles_str) == 1)g_map_closed is not allowed
 * else
 *  The fist tile of each tile block in @tiles_str is clipped towards right only if @direction == g_map_lef;
    The last tile ...towards left... only if == g_map_rght; 
 *  Others clipped towards both direction if necessary (exists gap).
 */
int reform_tiles_(String<Dna5> & seq1,
                  String<Dna5> & seq2,   //read 
                  String<Dna5> & comstr, //complement reverse of the read (seq2)
                  String<uint64_t> & tiles_str,
                  String<uint64_t> & tiles_end,
                  String<uint64_t> & clips,
                  String<uint64_t> & g_hs,
                  String<uint64_t> & g_hs_anchor,
                  String<uint64_t> & sp_tiles, //record tiles id that needs additional process(dups).
                  int pos_str,
                  int pos_end,
                  int direction,
                  int thd_cord_gap,
                  int thD_tile_size,
                  float thD_err_rate)
{

    if (empty (tiles_str))
    {
        return 0; 
    }
    //Parms
    int64_t thd_dxy_min = 80; //skip ins del < this value todo tune later
    float thd_da_zero = thD_err_rate; 

    CmpInt64 g_cmpll;
    int sv_exists = 0;
    if (pos_end - pos_str == 1)
    {
        if (direction == g_map_rght)
        {
            reform_tile_ (seq1, seq2, comstr, g_hs, g_hs_anchor, tiles_str[pos_str], tiles_end[pos_str], g_sv_r, thD_tile_size);
            return 0;
        }
        else if (direction = g_map_left)
        {
            reform_tile_ (seq1, seq2, comstr, g_hs, g_hs_anchor, tiles_str[pos_str], tiles_end[pos_str], g_sv_l, thD_tile_size);
            return 0;
        }
        else //== g_map_closed
        {
            return 1;
        }
    }
    else
    {
        if (direction == g_map_left)
        {
            reform_tile_ (seq1, seq2, comstr, g_hs, g_hs_anchor, tiles_str[pos_str], tiles_end[pos_str], g_sv_l, thD_tile_size);
        }
        for (int i = pos_str + 1; i < pos_end; i++)
        {
            if (i > 0 && is_tile_end (tiles_str[i - 1]) && direction == g_map_left)
            {
                reform_tile_ (seq1, seq2, comstr, g_hs, g_hs_anchor, tiles_str[i], tiles_end[i], g_sv_l, thD_tile_size);
                continue;
            }
            int64_t dx = tile_distance_x(tiles_str[i - 1], tiles_str[i], uint64_t(length(seq2)));
            int64_t dy = tile_distance_y(tiles_str[i - 1], tiles_str[i], uint64_t(length(seq2))); 
            if (! isTilesConsecutive_(tiles_str[i - 1], tiles_str[i], thd_cord_gap)) //inv
            {
                if (is_diff_anchor (tiles_str[i - 1], tiles_str[i], 1, thd_dxy_min, thd_da_zero)) 
                {
                    insertClipStr(sp_tiles, i - 1);  //tiles[i - 1], and tiles[i] needs additional process.
                    insertClipEnd(sp_tiles, i);
                }

                reform_tile_ (seq1, seq2, comstr, g_hs, g_hs_anchor, tiles_str[i - 1], tiles_end[i - 1], g_sv_r, thD_tile_size);
                reform_tile_ (seq1, seq2, comstr, g_hs, g_hs_anchor, tiles_str[i], tiles_end[i], g_sv_l, thD_tile_size);
            }
            if (direction == g_map_rght && is_tile_end(tiles_str[i]))
            {
                reform_tile_ (seq1, seq2, comstr, g_hs, g_hs_anchor, tiles_str[i] , tiles_end[i], g_sv_r, thD_tile_size);
            }
        }
        return 0;
    }
}

/*
 @direction == g_map_left @gap_str is ommited::will not append to tiles to be clipped
 @direction == g_map_rght @gap_end is ommited::will not append to tiles to be clipped
 */
int reform_tiles(String<Dna5> & seq1,
                 String<Dna5> & seq2,
                 String<Dna5> & comstr, //complement reverse of the read (seq2)
                 String<uint64_t> & tiles_str,
                 String<uint64_t> & tiles_end,
                 String<uint64_t> & clips,
                 String<uint64_t> & g_hs,
                 String<uint64_t> & g_hs_anchor,
                 String<uint64_t> & sp_tiles,
                 uint64_t gap_str,
                 uint64_t gap_end, 
                 int direction,
                 int thd_cord_gap,
                 int thD_tile_size,
                 float thD_err_rate
                )
{
    //step.1 init tiles_str, tiles_end;
    //Insert head_tile and tail tile at the front and end for each block of tiles
    uint64_t head_tile = gap_str;
    uint64_t tail_tile = shift_tile(gap_end, -thD_tile_size, -thD_tile_size); 
    if (get_tile_x(gap_end) > thD_tile_size && get_tile_y (gap_end) > thD_tile_size)
    {
        tail_tile = shift_tile(gap_end, -thD_tile_size, -thD_tile_size);
    }
    else
    {
        tail_tile = head_tile;
    }
    remove_tile_sgn(head_tile);
    remove_tile_sgn(tail_tile);

    String<uint64_t> tiles_str_tmp;
    String<uint64_t> tiles_end_tmp;
    if (direction != g_map_left)
    {
        appendValue (tiles_str_tmp, head_tile);
    }
    if (empty(tiles_str))
    {
        if (direction != g_map_rght)
        {
            appendValue (tiles_str_tmp, tail_tile);
        }
    }
    else
    {
        for (int i = 0; i < length(tiles_str); i++)
        {
            appendValue(tiles_str_tmp, tiles_str[i]);
            if (is_tile_end(tiles_str[i]) || i == length(tiles_str) - 1)
            {
                if (direction != g_map_rght)
                {
                    copy_tile_sgn(back(tiles_str_tmp), tail_tile);
                    remove_tile_sgn(back(tiles_str_tmp));
                    appendValue(tiles_str_tmp, tail_tile);
                }
                if (i != length(tiles_str) - 1 && direction != g_map_left)
                {
                    copy_tile_sgn(tiles_str[i + 1], head_tile);
                    appendValue(tiles_str_tmp, head_tile);
                }
            }
        }
    }
    int64_t d = shift_tile (0ULL, thD_tile_size, thD_tile_size);
    resize(tiles_end_tmp, length(tiles_str_tmp));
    for (int i = 0; i < length(tiles_str_tmp); i++)
    {
        tiles_end_tmp[i] = tiles_str_tmp[i] + d;
    }

    //step.2 reform tiles:clip and break block if necessary
    reform_tiles_(seq1, seq2, comstr, 
                tiles_str_tmp,
                tiles_end_tmp,
                clips, 
                g_hs, 
                g_hs_anchor, 
                sp_tiles,
                0, length(tiles_str_tmp),
                direction, 
                thd_cord_gap, 
                thD_tile_size,
                thD_err_rate);
    tiles_str = tiles_str_tmp;
    tiles_end = tiles_end_tmp;
    return 0;
}

//return if two blocks can be chained according to last and first cords.
//@cord1, @cord2 are the last and first cord of each block
//return 1 if 2 blocks can be linked; else 0
int b_ins_link = 1;
int b_gap_link = 2; 
int b_inv_link = 4;
int isBlocksLinkable(uint64_t cord1, 
                     uint64_t cord2, 
                     uint64_t cmprevconst,
                     int64_t thd_max_chain_distance)
{
    if (!is_cord_block_end(cord1)) 
    {
        return 0;
    }
    if (get_cord_strand (cord1 ^ cord2))
    {
        //todo::check the logic of the reverse strand in this part
        int64_t x1 = get_cord_x(cord1);
        int64_t x2 = get_cord_y(cord2);
        int64_t y1 = get_cord_y(cord1);
        int64_t y2 = cmprevconst - get_cord_y(cord2);
        if (y1 < y2 && y1 + thd_max_chain_distance > y2 &&
            std::abs(x2 - x1) < thd_max_chain_distance) // ins or dup
        {
            return b_ins_link | b_inv_link;
        }
        else if (y1 < y2 && y1 + thd_max_chain_distance > y2 &&
                 x1 < x2 && x1 + thd_max_chain_distance > x2) //regular gap
        {
            return b_gap_link; 
        }
        else
        {
            return 0;
        }
        return 0;
    }
    else
    {
        int64_t x1 = get_cord_x(cord1);
        int64_t x2 = get_cord_x(cord2);
        int64_t y1 = get_cord_y(cord1);
        int64_t y2 = get_cord_y(cord2);
        if (y1 < y2 && y1 + thd_max_chain_distance > y2 &&
            std::abs(x2 - x1) < thd_max_chain_distance) //ins or dup
        {
            return b_ins_link;
        }
        else if (y1 < y2 && y1 + thd_max_chain_distance > y2 &&
                 x1 < x2 && x1 + thd_max_chain_distance > x2) //regular gap
        {
            return b_inv_link; 
        }
        else
        {
            return 0;
        }
    }
}

/**
 * identify svs which is across different blocks
 */
int try_blocks_sv_ (String<uint64_t> & cords, 
                    String<uint64_t> & clips,
                    uint64_t cmprevconst,
                    int64_t thd_max_chain_distance,
                    int64_t thD_tile_size)
{
    if (length(cords) < 2)
    {
        return 0;
    }
    for (int i = 1; i < length(cords); i++)
    {
        if (i != 1 && is_cord_block_end(cords[i - 1])) 
        {
            int link_type = isBlocksLinkable(cords[i - 1], cords[i], cmprevconst, thd_max_chain_distance);
                //insert
            if (link_type & b_inv_link)
            {

            }
            if (link_type & b_ins_link)
            {
                if (get_cord_x(cords[i - 1]) > get_cord_x(cords[i])) //dup
                {
                    insertClipStr(clips, cords[i]);
                    insertClipEnd(clips, shift_cord(cords[i - 1], thD_tile_size, thD_tile_size));
                }
                else //ins
                {
                    insertClipStr(clips, shift_tile(cords[i - 1], thD_tile_size, thD_tile_size));
                    insertClipEnd(clips, cords[i]);
                }
            }
            else if (link_type & b_gap_link)
            {

            }

        }
    }
    return 0;
}

/**
 * shortcut to insert @tiles at @cords[@pos] or @cords[@pos + 1] according to the @direction
 * if direction = g_map_left: erase back(@tiles) and insert(cords, @pos), 
 * if direction = g_map_closed then insert(cords, @pos)
 * if direction = g_map_rght: erase @tiles[0] and insert(cords, @pos + 1);
 * @tiles is required to have at least 2 tiles;
 * @tile[0] and back(@tiles) are the clipped cords of gap_str, gap_end
 * They will replace the two joint cords between which @tiles are inserted 
 */
int insert_tiles2Cords_(String<uint64_t> & cords, 
                        unsigned & pos,
                        String<uint64_t> & tiles,
                        int direction,
                        int thd_max_segs_num)
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
    insert_tiles2Cords_(cords_str, pos, tiles_str, direction, thd_max_segs_num);
    insert_tiles2Cords_(cords_end, postmp, tiles_end, direction, thd_max_segs_num);
    return 0;
}

/*----------------------  Gap main func ----------------------*/
int mapExtend_(StringSet<String<Dna5> > & seqs,
               String<Dna5> & read,
               String<Dna5> & comstr,
               uint64_t gap_str, 
               uint64_t gap_end, 
               uint64_t anchor_lower,
               uint64_t anchor_upper,
               String<uint64_t> & g_hs,
               String<uint64_t> & g_anchor,
               StringSet<FeaturesDynamic > & f1,
               StringSet<FeaturesDynamic > & f2,
               String<uint64_t> & tiles_str,     //results
               String<uint64_t> & tiles_end,     //results
               String<uint64_t> & clips,     //results 
               int direction,
               int thd_cord_gap, 
               int thD_tile_size,
               int thd_cord_remap, 
               float thD_err_rate,
               int64_t thd_dxy_min,
               MapGAnchorParm const & parm1)
{
    String <Dna5> & ref = seqs[get_cord_id(gap_str)];
    String<uint64_t> sp_tiles; 
    float thd_da_zero = thD_err_rate; 
    g_mapHs_(ref, read, comstr,
             g_hs, g_anchor, tiles_str, f1, f2,
             gap_str, gap_end, anchor_lower, anchor_upper,
             thd_da_zero, thD_tile_size, thd_dxy_min,
             g_map_rght, parm1
       );
    reform_tiles(ref, read, comstr, 
                 tiles_str, tiles_end,
                 clips, g_hs, g_anchor, sp_tiles, 
                 gap_str, gap_end, direction, 
                 thd_cord_gap, thD_tile_size, thD_err_rate);
}

/*
 * Map any possible anchors within the gap
 */
int mapGap_ (StringSet<String<Dna5> > & seqs,
             String<Dna5> & read,
             String<Dna5> & comstr,
             uint64_t gap_str, 
             uint64_t gap_end, 
             String<uint64_t> & g_hs,
             String<uint64_t> & g_anchor,
             StringSet<FeaturesDynamic > & f1,
             StringSet<FeaturesDynamic > & f2,
             String<uint64_t> & tiles_str,     //results
             String<uint64_t> & tiles_end,     //results
             String<uint64_t> & clips,         //results 
             int direction,
             int thd_cord_gap, 
             int thD_tile_size,
             int thd_cord_remap, 
             float thD_err_rate,
             int64_t thd_dxy_min,
             MapGAnchorParm const & parm1
            )
{

//>>debug
    CmpInt64 g_cmpll;
    float thd_da_zero = thD_err_rate; 
    clear(tiles_str);
    clear(tiles_end);
    _DefaultHit.unsetBlockEnd(gap_str); //remove cord sgn, format cord to tiles
    _DefaultHit.unsetBlockEnd(gap_end);

    String <Dna5> & ref = seqs[get_cord_id(gap_str)];
    int64_t x1 = get_cord_x(gap_str);
    int64_t x2 = get_cord_x(gap_end);
    int64_t y1 = get_cord_y(gap_str);
    int64_t y2 = get_cord_y(gap_end);
    int64_t da = x2 - x1 - y2 + y1;
    int64_t shift_x, shift_y;
    String<uint64_t> sp_tiles; 

    if (get_cord_strand(gap_str ^ gap_end))
    {
        if (direction != g_map_closed)
        {
            return -1; //this case is not allowed
        }
        String<uint64_t> tiles_str1;
        String<uint64_t> tiles_str2;
        String<uint64_t> tiles_end1;
        String<uint64_t> tiles_end2;
        String<uint64_t> clips1;
        String<uint64_t> clips2;
        int direction1 = g_map_rght;
        int direction2 = g_map_left;

        g_cmpll.min(shift_x, x2 - x1) << int64_t(length(ref) - 1 - get_cord_x(gap_str));
        g_cmpll.min(shift_y, (x2 - x1) * (1 + thD_err_rate)) << int64_t(length(read) - 1 - get_cord_y(gap_str));
        uint64_t gap_str1 = gap_str;
        uint64_t gap_end1 = shift_cord (gap_str, shift_x, shift_y);
        mapExtend_ (seqs, read, comstr, 
                    gap_str1, gap_end1, 
                    LLMIN, LLMAX,
                    g_hs, g_anchor, 
                    f1, f2,  
                    tiles_str1, tiles_end1, clips1, 
                    direction1,
                    thd_cord_gap, thD_tile_size, thd_cord_remap,
                    thD_err_rate, thd_dxy_min, parm1);

        g_cmpll.min(shift_x, x2 - x1) << int64_t(get_cord_x(gap_end));
        g_cmpll.min(shift_y, (x2 - x1) * (1 + thD_err_rate)) << int64_t(get_cord_y(gap_end));
        uint64_t gap_str2 = shift_cord (gap_end, -shift_x, -shift_y);
        uint64_t gap_end2 = gap_end;
        mapExtend_ (seqs, read, comstr, 
                    gap_str2, gap_end2, 
                    LLMIN, LLMAX,
                    g_hs, g_anchor, 
                    f1, f2,  
                    tiles_str2, tiles_end2, clips2, 
                    direction2,
                    thd_cord_gap, thD_tile_size, thd_cord_remap,
                    thD_err_rate, thd_dxy_min, parm1);

        if (!empty(tiles_str1))
        {
            append(tiles_str, tiles_str1);
            append(tiles_end, tiles_end1);
        }
        if (!empty(tiles_str2))
        {
            append(tiles_str, tiles_str2);
            append(tiles_end, tiles_end2);
        }
    }
    else if (get_cord_y(gap_end) - get_cord_y(gap_str) > (length(g_hs)) ||
             get_cord_x(gap_end) - get_cord_x(gap_str) > (length(g_hs)) ||
             get_cord_y(gap_end) > length(read) - 1 ||
             get_cord_y(gap_str) > length(read) - 1 ||
             get_cord_x(gap_end) > length(ref) - 1 ||
             get_cord_x(gap_str) > length(ref) - 1)
    {
        //return 1;
        //todo
    }
    else
    {
        g_mapHs_(ref, read, comstr,
                 g_hs, g_anchor, tiles_str, f1, f2,
                 gap_str, gap_end,
                 thd_da_zero, thD_tile_size, thd_dxy_min,
                 direction,
                 parm1
                );

        reform_tiles(ref, read, comstr, 
                     tiles_str, tiles_end,
                     clips, g_hs, g_anchor, sp_tiles,
                     gap_str, gap_end, direction, 
                     thd_cord_gap, thD_tile_size, thD_err_rate);
        //try dup for each sp_tiles element.
        for (int j = 0; j < getClipsLen(sp_tiles); j++)
        {
            String <uint64_t> dup_tiles_str;
            String <uint64_t> dup_tiles_end;
            String <uint64_t> dup_clips;
            String <uint64_t> dup_sp_tiles;
            int dup_direction = g_map_closed;
            uint64_t dup_str = tiles_str[getClipStr(sp_tiles, j)];
            uint64_t dup_end = tiles_end[getClipEnd(sp_tiles, j)];
            try_tiles_dup (ref, read, comstr, 
                           f1, f2,
                           g_hs, g_anchor, dup_tiles_str,
                           dup_str, dup_end, 
                           thD_err_rate, 
                           thD_tile_size,
                           parm1);
            reform_tiles(ref, read, comstr,
                         dup_tiles_str, dup_tiles_end,
                         dup_clips, g_hs, g_anchor, dup_sp_tiles,
                         dup_str, dup_end, dup_direction,
                         thd_cord_gap, thD_tile_size, thD_err_rate);
        }
    }

    return 0;
}

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
                             uint64_t thD_tile_size)
{
    if (empty(cords_end))
    {
        resize (cords_end, length(cords_str));
        for (unsigned i = 0; i < length(cords_str); i++)
        {
            cords_end[i] = shift_cord (cords_str[i], thD_tile_size, thD_tile_size);
        }
    }
    set_cord_xy (cords_str[i], get_cord_x(cord1), get_cord_y(cord1));
    set_cord_xy (cords_end[i], get_cord_x(cord2), get_cord_y(cord2));
}

/**
 * Re-map gaps in cords.
 * Gaps at the front or end of the cords are also remapped.
 * !!NOTE::Each block in the cords are required to be sorted according to the y value of the first cord, such that isBlocksLinkable can work appropriately.
 */
int mapGaps(StringSet<String<Dna5> > & seqs, 
            String<Dna5> & read, 
            String<Dna5> & comstr,
            String<uint64_t> & cords_str, 
            String<uint64_t> & cords_end, 
            String<uint64_t> & g_hs,
            String<uint64_t> & g_anchor,
            String<uint64_t> & clips, // string for clips cords
            String<UPair> & apx_gaps,
            StringSet<FeaturesDynamic > & f1,
            StringSet<FeaturesDynamic > & f2,
            int64_t thd_gap, 
            int64_t thD_tile_size,
            float thD_err_rate
           )
{
    CmpInt64 g_cmpll;
    if (length(cords_str) <= 1)
    {
        return 0;
    }
    String <uint64_t> tiles_str;
    String <uint64_t> tiles_end;
    String<UPair> str_ends;
    String<UPair> str_ends_p;

    int64_t shift_x;
    int64_t shift_y;

    int thd_max_segs_num = 1000; //max segs num allowed in each gap, gaps > this will abort all tiles

    uint64_t thd_max_extend = 2000;  //Important::tune::whether process in gap or pmpfiner.h
    uint64_t thd_max_extend_x = thd_max_extend;
    uint64_t thd_max_extend_y = thd_max_extend; //>large gaps are supposed to be handled during the apx part.
    int64_t thd_dxy_min = 80;
    int64_t thd_extend_xy = 3; //extend at first or last cord of seqs
    int64_t thd_max_chain_distance = thd_max_extend;
    int64_t block_size = thD_tile_size; 
    int64_t thd_cord_size = thD_tile_size; //NOTE::this is the cord size, may be different from thD_tile_size.
    int64_t thd_cord_remap = 100;
    int64_t thd_cord_gap = thd_gap + block_size;

    MapGAnchorParm parm1(thD_tile_size, thD_err_rate);

    clear(apx_gaps);
    gather_blocks_ (cords_str, str_ends, str_ends_p, length(read), thd_cord_gap, thd_cord_size, 0);
    gather_gaps_y_ (cords_str, str_ends, apx_gaps, length(read), thd_cord_gap);

    int count = 0;
    //NOTE cords_str[0] is the head, regular cords_str starts from 1
    for (unsigned i = 1; i < length(cords_str); i++)
    {
        unsigned sid = get_cord_id(cords_str[i]);
        uint64_t x1 = get_cord_x(cords_str[i - 1]);
        uint64_t y1 = get_cord_y(cords_str[i - 1]);
        uint64_t x2 = get_cord_x(cords_str[i]);
        uint64_t y2 = get_cord_y(cords_str[i]);
        int64_t dx = x2 - x1;
        int64_t dy = y2 - y1;
        int direction;
        if (_DefaultCord.isBlockEnd(cords_str[i - 1]))  //clip first cord
        {
            shift_x = std::min (int64_t(length(seqs[sid]) - 1 - get_cord_x(cords_str[i])),  block_size);
            shift_y = std::min (int64_t(length(read) - 1 - get_cord_y(cords_str[i])), block_size);
            uint64_t gap_end = shift_cord(cords_str[i], shift_x, shift_y);
            uint64_t gap_str;
            if (get_cord_y(gap_end) > thd_cord_gap) 
            {            
                shift_x = std::min(thd_max_extend, get_cord_x(gap_end));
                shift_y = std::min(thd_max_extend, get_cord_y(gap_end));
                shift_x = std::min(shift_x, shift_y * thd_extend_xy);
                uint64_t infi_cord = shift_cord(gap_end, -shift_x, -shift_y); 
                direction = g_map_left;
                gap_str = infi_cord;
                _DefaultHit.unsetBlockEnd(gap_str);
                _DefaultHit.unsetBlockEnd(gap_end);

                int max_gap_overlap_y = _getMaxGapsyOverlap(apx_gaps, gap_str, gap_end);
                if (max_gap_overlap_y > thd_cord_gap)
                {
                    mapGap_ (seqs, read, comstr, 
                             gap_str, gap_end, 
                             g_hs, g_anchor, 
                             f1, f2,  
                             tiles_str, 
                             tiles_end,
                             clips, 
                             direction,
                             thd_cord_gap, 
                             thD_tile_size,
                             thd_cord_remap,
                             thD_err_rate,
                             thd_dxy_min,
                             parm1);
                    insert_tiles2Cords_(cords_str, cords_end, i, tiles_str, tiles_end, direction, thd_cord_size, thd_max_segs_num);
                }
                else
                {
                    uint64_t tile_str = cords_str[i];
                    uint64_t tile_end = shift_cord(cords_str[i], block_size, block_size);

                    reform_tile_ (seqs[get_cord_id(cords_str[i])],
                          read,
                          comstr,
                          g_hs,
                          g_anchor,
                          tile_str, 
                          tile_end, 
                          g_sv_l, 
                          thD_tile_size
                        );
                    _updateCordsStrEndValue(cords_str, cords_end, i, tile_str, tile_end, block_size);
                    
                    int64_t anchor_base = get_cord_x(gap_end) - get_cord_y(gap_end);
                    int64_t anchor_lower = anchor_base - 100;
                    int64_t anchor_upper = anchor_base + 100;
                    mapExtend_ (seqs, read, comstr, 
                                gap_str, gap_end, 
                                anchor_lower, anchor_upper,
                                g_hs, g_anchor, 
                                f1, f2,  
                                tiles_str, 
                                tiles_end,
                                clips, 
                                direction,
                                thd_cord_gap, 
                                thD_tile_size,
                                thd_cord_remap,
                                thD_err_rate,
                                thd_dxy_min, parm1);
                    insert_tiles2Cords_(cords_str, cords_end, i, tiles_str, tiles_end, direction, thd_cord_size, thd_max_segs_num);
                }
            }
        }
        else if (!isCordsConsecutive_(cords_str[i - 1], cords_str[i], thd_cord_gap))
        {
            g_cmpll.min(shift_x, block_size) 
                    << int64_t(length(seqs[sid]) - 1 - get_cord_x(cords_str[i]));
            g_cmpll.min(shift_y, block_size)
                    << int64_t(length(read) - 1 - get_cord_y(cords_str[i]));
            uint64_t gap_end = shift_cord(cords_str[i], shift_x, shift_y);
            uint64_t gap_str = cords_str[i - 1]; 
            direction = g_map_closed;
            _DefaultHit.unsetBlockEnd(gap_str);
            _DefaultHit.unsetBlockEnd(gap_end);

            mapGap_(seqs, read, comstr, 
                    gap_str, gap_end, 
                    g_hs, g_anchor, 
                    f1, f2, 
                    tiles_str, 
                    tiles_end,
                    clips, 
                    direction,
                    thd_cord_gap, 
                    thD_tile_size,
                    thd_cord_remap,
                    thD_err_rate,
                    thd_dxy_min, 
                    parm1);
            insert_tiles2Cords_(cords_str, cords_end, i, tiles_str, tiles_end, direction, thd_cord_size, thd_max_segs_num);
            //<<debug
            //if (count++ == 1)
            //return 0;
            //>>debug
        }
        if (_DefaultHit.isBlockEnd(cords_str[i]))  ///right clip end cord
        {
            uint64_t gap_str = cords_str[i];
            if (length(read) - 1 - get_cord_y(gap_str) > thd_cord_gap)
            {
                shift_x = std::min(thd_max_extend_x, 
                                  length(seqs[sid]) - get_cord_x(gap_str) - 1);
                shift_y = std::min(thd_max_extend_y, 
                                   length(read) - get_cord_y(gap_str) - 1);
                shift_x = std::min(shift_x, shift_y * thd_extend_xy);
                uint64_t infi_cord = _DefaultCord.shift(gap_str, shift_x, shift_y);
                uint64_t gap_end = infi_cord;
                direction = g_map_rght;
                _DefaultHit.unsetBlockEnd(gap_str);
                _DefaultHit.unsetBlockEnd(gap_end);
                int max_gap_overlap_y = _getMaxGapsyOverlap(apx_gaps, gap_str, gap_end);

                if (max_gap_overlap_y > thd_cord_gap)
                {
                    mapGap_ (seqs, read, comstr, 
                             gap_str, gap_end, 
                             g_hs, g_anchor, 
                             f1, f2,  
                             tiles_str, 
                             tiles_end, 
                             clips, 
                             direction,
                             thd_cord_gap, 
                             thD_tile_size, 
                             thd_cord_remap,
                             thD_err_rate,
                             thd_dxy_min,
                             parm1);
                    insert_tiles2Cords_(cords_str, cords_end, i, tiles_str, tiles_end, direction, thd_cord_size, thd_max_segs_num);
                }
                else
                {
                    
                    uint64_t tile_str = cords_str[i];
                    uint64_t tile_end = shift_cord(cords_str[i], block_size, block_size);
                    reform_tile_ (seqs[get_cord_id(cords_str[i])],
                                  read,
                                  comstr,
                                  g_hs,
                                  g_anchor,
                                  tile_str, 
                                  tile_end, 
                                  g_sv_r, 
                                  thD_tile_size
                                  );
                    _updateCordsStrEndValue(cords_str, cords_end, i, tile_str, tile_end, block_size);

                    int64_t anchor_base = get_cord_x(gap_str) - get_cord_y(gap_end);
                    int64_t anchor_lower = anchor_base - 100;
                    int64_t anchor_upper = anchor_base + 100; 

                    mapExtend_ (seqs, read, comstr, 
                                gap_str, gap_end, 
                                anchor_lower, anchor_upper,
                                g_hs, g_anchor, 
                                f1, f2,  
                                tiles_str, 
                                tiles_end,
                                clips, 
                                direction,
                                thd_cord_gap, 
                                thD_tile_size,
                                thd_cord_remap,
                                thD_err_rate,
                                thd_dxy_min, 
                                parm1);
                    insert_tiles2Cords_(cords_str, cords_end, i, tiles_str, tiles_end, direction, thd_cord_size, thd_max_segs_num);
                }
            }
        }
    }
    //try_blocks_sv_(cords_str, clips, length(read) - 1, thd_max_chain_distance, thD_tile_size);
    return 0;
}

/*=====  End of Index free Map and clip  ======*/

