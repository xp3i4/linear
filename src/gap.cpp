//GNode: N/A[1]|xval[32]|strand[1]|coordinate[30]
#include <utility> 
#include "base.h"
#include "shape_extend.h"
#include "cords.h"
#include "gap.h"
 
/*=============================================
=             Interface function              =
=============================================*/
using std::endl;

/**
 * class ClipRecords:
 * based on struct tile!!
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
        //set_cord_id (back(clips), get_cord_id(clip));
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

int print_clips_gvf_(StringSet<String<uint64_t> > & clips, 
                     StringSet<CharString> & readsId, 
                     StringSet<CharString> & genomesId,
                     std::ofstream & of)
                     //std::string outputPrefix)
{
    //std::string file_path = outputPrefix + ".gvf";
    //std::cerr << "[]::filepath " << file_path << "\n";
    //of.open(toCString(file_path));
    of << "##gvf-version 1.10\n";
    std::string source = ".";
    std::string type = ".";
    for (unsigned i = 0; i < length(clips); i++)
    {
        //std::cout << "pcg1 " << i << " " << getClipsLen(clips[i]) << " " << length(clips) << "\n";
        //continue;
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

/*=======  End of interface function  =======*/


struct GNodeBase
{
    const unsigned xBitLen;
    const unsigned sBitLen;
    const unsigned cBitLen; 
    const uint64_t xmask;
    const uint64_t smask;
    const uint64_t cmask;

    GNodeBase():
        xBitLen(32),
        sBitLen(1),
        cBitLen(30),
        xmask((1ULL << xBitLen) - 1),
        smask((1ULL << sBitLen) - 1),
        cmask((1ULL << cBitLen) - 1)
        {}
}_defaultGNodeBase;

struct GNode
{
     void setValue (uint64_t & val, uint64_t const & xval, 
                   uint64_t const & strand, uint64_t const & coordinate, 
                   uint64_t const & sbit = _defaultGNodeBase.sBitLen,
                   uint64_t const & cbit = _defaultGNodeBase.cBitLen)
    {
        val = (xval << (cbit + sbit)) + (strand << cbit) + coordinate;
    }
     uint64_t makeValue (uint64_t const & xval, uint64_t & strand, uint64_t const & coordinate, 
                        uint64_t const & sbit = _defaultGNodeBase.sBitLen,
                        uint64_t const & cbit = _defaultGNodeBase.cBitLen)
    {
        return (xval << (cbit + sbit)) + (strand << cbit) + coordinate;
    }
     uint64_t getXValue (uint64_t const & xval, 
                        uint64_t const & bit = _defaultGNodeBase.cBitLen + _defaultGNodeBase.sBitLen, 
                        uint64_t const & xmask = _defaultGNodeBase.xmask)
    {
        return (xval >> bit) & xmask;
    }
     uint64_t getStrand (uint64_t const & xval, 
                        uint64_t const & bit = _defaultGNodeBase.cBitLen, 
                        uint64_t const & mask = _defaultGNodeBase.smask)
    {
        return (xval >> bit) & mask;
    }
     uint64_t getCoord (uint64_t const & xval, 
                        uint64_t const & mask = _defaultGNodeBase.cmask)
    {
        return (xval & mask);
    }
}_defaultGNode;

/*
 *NOTE! the following parameters highly correlated.
 *Do not change them independently 
 */
int const g_shape_len = 8;
int const g_thd_anchor = 6;
float const g_thd_anchor_density = 0.03;
float const g_thd_error_percent = 0.2;

struct GIndex
{
    String <uint64_t> g_hs;
    String <uint64_t> g_dir;
    LShape shape;
    GIndex():
    shape(g_shape_len){};
    GIndex(unsigned shape_len):
    shape(shape_len)
    {};
};

int g_createDir_(String<Dna5> & seq, uint64_t gs_start, uint64_t gs_end, 
                 String<uint64_t> & g_hs, String<uint64_t> & g_dir, 
                 LShape & shape)
{
    hashInit(shape, begin(seq) + gs_start);
    unsigned count = 0; 
    unsigned const step = 10;
    unsigned countx = 0;
    uint64_t preX = 0;
    clear(g_hs);
    for (uint64_t k = gs_start; k < gs_end; k++)
    {
        if (++count == step)  //collecting every 10 bases
        {
            hashNext(shape, begin(seq) + k);
            if (shape.XValue == preX && countx < shape.span - shape.weight)
            {
                ++countx;
            }
            else
            {
                //TODO: k - getT(shape)
                appendValue(g_hs, _defaultGNode.makeValue(shape.XValue, shape.strand, k));
                preX = shape.XValue;
                countx = 0;
            }
            count = 0;
        }
        else
        {
            hashNexth(shape, begin(seq) + k);
        }
    }
    appendValue(g_hs, ~0);
    std::sort (begin(g_hs), end(g_hs) - 1);
    count = 0;
    for (uint64_t k = 0; k < length(g_hs) - 1; k++)
    {
        if (_defaultGNode.getXValue(g_hs[k] ^ g_hs[k + 1])) //g_hs[k].xval != g_hs[k+1].xval
        {
            g_dir[_defaultGNode.getXValue(g_hs[k])] = k - count;
            count = 0;
        }
        else
        {
            ++count; 
        }
    }
    return 0;
}    

void g_createDir(String<Dna5> & seq, uint64_t gs_start, uint64_t gs_end, GIndex & g_index) 
{
    clear(g_index.g_dir);
    resize(g_index.g_dir,  1 << (g_index.shape.weight * 2));
    g_createDir_(seq, gs_start, gs_end, g_index.g_hs, g_index.g_dir, g_index.shape);
}

//ExtendSa gs[36]|rs[14]|gb[14]
//gs: starting coordinate in genome
//rs: starting coordinate in read
//gb: number of blocks in the gap
//Gap no larger than gb x blockSize
/*
struct gapbase_
{
    const unsigned gBitLen = 36;
    const unsigned rBitLen = 14;
    const unsigned bBitLen = 14;
    const unsigned bBit = 6; // blockSize = 1 << bBitLen
    const unsigned gmask = (1ULL << gBitLen) - 1;
    const unsigned rmask = (1ULL << rBitLen) - 1;
    const unsigned bmask = (1ULL << bBitLen) - 1;
} _defaultgapbase_;

struct Gap_
{
     uint64_t getGStart(uint64_t & val, 
                       uint64_t const & bit = _defaultgapbase_.rBitLen + _defaultgapbase_.bBitLen, 
                       uint64_t const & mask = _defaultgapbase_.gmask, 
                       uint64_t const & bbit = _defaultgapbase_.bBit)
    {
        return ((val >> bit) & mask) << bbit;
    }
     uint64_t getRStart(uint64_t & val, 
                       uint64_t const & bit = _defaultgapbase_.bBitLen, 
                       uint64_t const & mask = _defaultgapbase_.rmask,
                       uint64_t const & bbit = _defaultgapbase_.bBit)

    {
        return ((val >> bit) & mask) << bbit;
    }
     uint64_t getGEnd(uint64_t & val, 
                    uint64_t const & bit = _defaultgapbase_.gBitLen, 
                    uint64_t const & smask = _defaultgapbase_.smask,
                    uint64_t const & gmask = _defaultgapbase_.gmask)
    {
        return ((val >> bit) & smask) + (val & gmask);
    }
     uint64_t makeValue (uint64_t start, uint64_t len, 
                    uint64_t const & bit = _defaultgapbase_.gBitLen)
    {
       return (start << bit) + len;
    }
}_defaultGap_;
*/

//gap id[10]|start[30]|gap length[24]

struct gapbase_
{
    const unsigned iBitLen = 10;
    const unsigned cBitLen = 30;
    const unsigned bBitLen = 24;
    const unsigned imask = (1ULL << iBitLen) - 1;
    const unsigned cmask = (1ULL << cBitLen) - 1;
    const unsigned bmask = (1ULL << bBitLen) - 1;
} _defaultgapbase_;

struct Gap_
{
     uint64_t getId (uint64_t & val, 
                           uint64_t const & bit = _defaultgapbase_.cBitLen + _defaultgapbase_.bBitLen,
                           uint64_t const & mask = _defaultgapbase_.bmask)
    {
        return (val >> bit) & mask;
    }
     uint64_t getStart(uint64_t & val, 
                      uint64_t const & bit = _defaultgapbase_.bBitLen, 
                      uint64_t const & mask = _defaultgapbase_.cmask)
    {
        return (val >> bit) & mask;
    }
     uint64_t getEnd(uint64_t & val, 
                    uint64_t const & bit = _defaultgapbase_.bBitLen, 
                    uint64_t const & cmask = _defaultgapbase_.cmask,
                    uint64_t const & bmask = _defaultgapbase_.bmask)
    {
        return ((val >> bit) & cmask) + (val & bmask);
    }
     uint64_t getLength (uint64_t & val, uint64_t const & mask = _defaultgapbase_.bmask)
    {
        return val & mask;
    }
     uint64_t makeValue (uint64_t start, uint64_t len, 
                        uint64_t const & bit = _defaultgapbase_.bBitLen)
    {
        return (start << bit) + len;
    }
}_defaultGap_;

struct Gap
{
    uint64_t gap1; //seq
    uint64_t gap2; //read
     uint64_t getId1()
    {
        return _defaultGap_.getId(gap1);
    }
     uint64_t getId2()
    {
        return _defaultGap_.getId(gap2);
    }
     uint64_t getStart1 ()
    {
        return _defaultGap_.getStart (gap1);
    }
     uint64_t getStart2 ()
    {
        return _defaultGap_.getStart (gap2);
    }
     uint64_t getEnd1 ()
    {
        return _defaultGap_.getEnd (gap1);
    }
     uint64_t getEnd2 ()
    {
        return _defaultGap_.getEnd (gap2);
    }
    Gap()
    {
        
    }
    Gap(uint64_t start1, uint64_t end1, uint64_t start2, uint64_t end2)
    {
        gap1 = _defaultGap_.makeValue (start1, end1 - start1);
        gap2 = _defaultGap_.makeValue (start2, end2 - start2);   
    }
     void setValue (uint64_t start1, uint64_t end1, uint64_t start2, uint64_t end2)
    {    
        gap1 = _defaultGap_.makeValue (start1, end1 - start1);
        gap2 = _defaultGap_.makeValue (start2, end2 - start2);   
    }
};

 Gap makeGap (uint64_t start1, uint64_t end1, uint64_t start2, uint64_t end2)
{    
    Gap gap (start1, end1, start2, end2);
    return gap;
}

/*
 * Merge gaps in the reference  
 */
/*
int mergeGap(String <std::pair<uint64_t, uint64_t> > & gaps)
{
    std::sort (begin(gaps), end(gaps), [](std::pair<uint64_t, uint64_t> & p1, std::pair<uint64_t, uint64_t> & p2) -> bool
    {return p1.first < p2.first}); // !Note: ascending rather than descending 
    unsigned pt = 0;
    for (unsigned k = 0; k < length(gaps) - 1; k++)
    {
        if (_defaultGap_.getEnd(gaps[k]) >  _defaultGap_.getStart(gaps[k + 1]))
        {
            len += _defaultGap_.getLength(gaps[k + 1])
        }
        else 
        {
            gaps[pt] = _defaultGap_.makeValue(_defaultGap_.getStart(gaps[pt]), 
                                                _defaultGap_.getEnd(gaps[k]));
            ++pt;
        }
    }
    resize (gaps, pt);
}
*/

//N/A[2]strand[1]|N/A[1]|anchor[40]|coord[20]
//strand = 1 or 0
struct ACoordBase
{
    unsigned const sBitLen = 1;
    unsigned const aBitLen = 40;
    unsigned const cBitLen = 20;
    unsigned const sBit = 61;
    
    uint64_t const amask = (1ULL << aBitLen) - 1;
    uint64_t const cmask = (1ULL << cBitLen) - 1;
    
}_defaultACoordBase;

struct ACoord
{
     uint64_t makeValue(uint64_t cf, uint64_t cr,  uint64_t strand,
                       uint64_t const & sbit = _defaultACoordBase.sBit,
                       uint64_t const & bit = _defaultACoordBase.cBitLen) 
    {
        return (strand << sbit) + ((cf + _nStrand(strand) * cr) << bit) + cr;
    }
     uint64_t reverseAnchor(uint64_t & anchor, uint64_t const & mask = _defaultACoordBase.amask)
    {
        return (-anchor) & mask;
    }
     uint64_t getAnchor (uint64_t val, 
                        uint64_t const & mask = _defaultACoordBase.amask, 
                        uint64_t const & bit = _defaultACoordBase.cBitLen,
                        uint64_t const & mask2 = (1ULL << _defaultACoordBase.sBit)) 
    {
        return (val & mask2) + ((val >> bit) & mask);
    }
     uint64_t getCoord (uint64_t val, 
                       uint64_t const & mask = _defaultACoordBase.cmask)
    {
        return val & mask;
    }
     uint64_t getX (uint64_t val, 
                   uint64_t const & bit = _defaultACoordBase.cBitLen, 
                   uint64_t const & mask = _defaultACoordBase.amask, 
                   uint64_t const & bit2 = _defaultACoordBase.sBit,
                   uint64_t const & mask3 = _defaultACoordBase.cmask
                  )
    {
        return ((val >> bit) & mask) + _nStrand((val >> bit2) & 1) * (val & mask3);
    }
     uint64_t getY (uint64_t val)
    {
        return getCoord(val);
    }
    
}_defaultACoord;

//Format::Tile tile_sign[2]|strand[1]|tileEnd[1]|x[40]|y[20]
//tile_sign:=1 start, 2 end, 0 body;
//0-61 bits same as the Format::Cord
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
void g_print_tiles_(String<uint64_t> & tiles, CharString str = "print_tiles")
{
    for (unsigned i = 0; i < length(tiles); i++)
    {
        std::cout << str << " " 
                  << i << " " 
                  << tiles[i] << " " 
                  << get_cord_strand(tiles[i]) << " " 
                  << get_cord_y(tiles[i]) << " " 
                  <<  get_cord_id(tiles[i]) << " " 
                  << get_cord_x(tiles[i]) << "\n";
        if (is_tile_end(tiles[i]))
        {
            std::cout << str << "\n\n";
        }
    }
}

 uint64_t acoord2Tile(uint64_t val, 
                         uint64_t const & bit = _defaultACoordBase.cBitLen,
                         uint64_t const & bit2 = _defaultACoordBase.sBit,
                         uint64_t const & mask = _defaultACoordBase.cmask)
{
    return val - ((_nStrand((val >> bit2) & 1) * (val & mask)) << bit);
}

int mapGap_(GIndex & g_index, String <Dna5> & read,  uint64_t start2, uint64_t end2, 
            String<uint64_t> & tile, uint64_t const thd_tileSize)
{
    //TODO: this parameter needs discussion.
    float thd_error_percent = 0.2; 
    uint64_t thd_min_segment = 100;
    //double time = sysTime();
    unsigned count = 0, step = 1;
    uint64_t preX = 0;
    String <uint64_t> anchor;
    hashInit(g_index.shape, begin(read) + start2);
    for (unsigned k = start2; k < end2; k++)
    {
        if (++count == step)
        {
            hashNext(g_index.shape, begin(read) + k);
            if (g_index.shape.XValue != preX)
            {
                uint64_t hsStart = g_index.g_dir[g_index.shape.XValue]; 
                //NOTE: For efficiency, the index.dir is uninitialized. 
                //So it contains pointers point to the memory out of the boundary of index.hs. 
                //So before searching the value, the boundary must be checked.
                if (hsStart < length(g_index.g_hs) )
                {
                    while (_defaultGNode.getXValue(g_index.g_hs[hsStart]) == g_index.shape.XValue)
                    {
                        uint64_t strand = _defaultGNode.getStrand(g_index.g_hs[hsStart]) ^ g_index.shape.strand;
                        //appendValue(anchor, ((_defaultGNode.getCoord(g_index.g_hs[hsStart]) + _nStrand(strand) * k) << 20) | k | (strand << 63));
                        appendValue(anchor, _defaultACoord.makeValue(_defaultGNode.getCoord(g_index.g_hs[hsStart]), k, strand));
                        ++hsStart;
                        // if strand == 1,  different strands
                        // + k
                        // if strand == 0,  the same strand
                        // - k
                    }   
                }
                preX = g_index.shape.XValue;
            }
            count = 0;
        }
        else
        {
            hashNexth(g_index.shape, begin(read) + k);
        }
    }
    std::sort (begin(anchor), end(anchor));
    appendValue (anchor, ~0);
    uint64_t prek = 0;
    unsigned anchor_len = 0, max_anchor_len = 0, max_prek = 0, max_k = 0;
    for (unsigned k = 0; k < length(anchor); k++)
    {
        //TODO: handle thd_min_segment, anchor 
        if (_defaultACoord.getAnchor(anchor[k]) - _defaultACoord.getAnchor(anchor[prek]) > 
            thd_error_percent * std::max(thd_min_segment, _defaultACoord.getCoord(anchor[k] - anchor[prek])))
        {
            if (anchor_len > max_anchor_len)
            {
                max_anchor_len = anchor_len;
                max_prek = prek;
                max_k = k;
                anchor_len = 0;
            }
            prek = k;
        }
        else
        {
            anchor_len++;
        }
    }
    std::cerr << "[]::mapGap_ " << max_prek << " " << max_k << "\n";
            
            std::sort (begin(anchor) + max_prek, begin(anchor) + max_k, 
                       [](uint64_t & s1, uint64_t & s2){
                return _defaultACoord.getCoord(s2) > _defaultACoord.getCoord(s1);
            });
            appendValue(tile, acoord2Tile(anchor[max_prek]));
            for (unsigned j = max_prek + 1; j < max_k; j++)
            {
                //add a tile if anchor[j] is not in the last tile
                if (_defaultACoord.getX(anchor[j]) > _defaultTile.getX(back(tile)) + thd_tileSize
                 || _defaultACoord.getY(anchor[j]) > _defaultTile.getY(back(tile)) + thd_tileSize)
                {
                    appendValue (tile, acoord2Tile(anchor[j - 1]));
                }
            }
        //TODO:anchor2Tile
    return 0;
}

/**
 * NOTE: Functions mapping gaps by the index.
 * Aborted due to reason regarding the speed.
 * Replaced by index-free methods.
 */ 
/*
int mapGap(String <Dna5> & seq, String <Dna5> & read, Gap & gap, 
           String<uint64_t> & tile, uint64_t const thd_tileSize)
{
    GIndex g_index;
    g_createDir(seq, gap.getStart1(), gap.getEnd1(), g_index);
    std::cerr << "[]::mg2 " << gap.getStart1() << " " << gap.getEnd1() << "\n";
    mapGap_ (g_index, read, gap.getStart2(), gap.getEnd2(), tile, thd_tileSize);
    //mapGap_ (g_index, _reverse(read), length(read) - end2 - 1, length(read) - start2 - 1);
    //optimizeGap();
    //TODO: map _reverse and do optimize
    return 0;
}
*/
/*
 * Extract gaps from cords then re-map gaps.
 * use index to map
 * []::mg3
 *
int mapGaps(StringSet<String<Dna5> > & seqs, String<Dna5> & read, String<uint64_t> & cords, 
            unsigned const thd_gap, unsigned const thd_tileSize)
{
    String <uint64_t> tile;
    Gap gap;
    //uint64_t thd_cordGap = _DefaultCord.createCord(thd_gap, thd_gap);
    uint64_t delta = (uint64_t)thd_gap / 2;
    unsigned count = 0;
    for (unsigned k = 2; k < length(cords); k++)
    {
        if (_DefaultCord.getCordX(cords[k] - cords[k - 1]) > thd_gap && 
            get_cord_y(cords[k] - cords[k - 1]) > thd_gap &&
            !_DefaultHit.isBlockEnd(cords[k - 1]))
        {
            clear(tile);
            gap = makeGap(_DefaultCord.getCordX(cords[k - 1]) + delta,
                          _DefaultCord.getCordX(cords[k]) + delta,
                          get_cord_y(cords[k - 1]) + delta, 
                          get_cord_y(cords[k]) + delta);
            //mapGap(seqs[_getSA_i1(_DefaultCord.getCordX(cords[k - 1]))], read, gap, tile, thd_tileSize);
            count += _DefaultCord.getCordX(cords[k] - cords[k - 1]);
            insert(cords, k, tile);
            k += length(tile);
        }
    }
    //TODO: tile -> cords

    return count;
}
*/

/*=============================================
=           Index free Map and clip           =
=============================================*/
/**
 * Part 2
 * NOTE: index free local mapping for gaps 
 */
/**
 * g_hs_anchor: N/A[13]|strand[1]|anchor[30]|cord_y[20]
 * @strand := shape strand of kmer in genome ^ shape strand of kmer in read
 * While the kmers are always picked up from the genome and read rather than
 * the reverse complement of the read. 
 * This is different from anchors used in chainning.
 * @anchor = x - y + g_hs_anchor_zero. g_hs_anchor_zero to restrict @anchor > 0.
   -g_hs_anchor_zero <= x - y < g_hs_anchor_zero
 */
uint64_t const g_hs_anchor_mask1 = (1ULL << 20) - 1;
uint64_t const g_hs_anchor_mask1_ = ~ g_hs_anchor_mask1;
uint64_t const g_hs_anchor_mask3 = (1ULL << 30) - 1;
uint64_t const g_hs_anchor_mask5 = (1ULL << 31) - 1;
uint64_t const g_hs_anchor_bit1 = 20;
uint64_t const g_hs_anchor_bit2 = 50;
uint64_t const g_hs_anchor_mask2 = ~(1ULL << 50);
uint64_t const g_hs_anchor_zero = 1ULL << (20);

int const g_sv_noe = 0;     //none
int const g_sv_inv = 1;     //inversion
int const g_sv_ins = 2;     //insertion
int const g_sv_del = 4;     //deletion
int const g_sv_trs = 8;     //translocation
int const g_sv_dup = 16;    //duplication
int const g_sv_l = 32;      //clip towards left
int const g_sv_r = 64;      //clip towards right
int const g_sv_gap = 128;   //gap unmapped:inversion duplication translocation.
int const g_sv_shrink = 256; // -->|  |<-- clip norml to the gap  (gap shrink).
int const g_sv_extend = 512; // |<--  -->| clip the gap to normal (gap extend)

int const g_align_left = -1;
int const g_align_closed = 0;
int const g_align_right = 1;
uint64_t g_hs_anchor_getCord (uint64_t anchor)
{
    return anchor & g_hs_anchor_mask1;
}

uint64_t g_hs_anchor_getAnchor (uint64_t anchor)
{
    return ((anchor >> g_hs_anchor_bit1) - g_hs_anchor_zero)& g_hs_anchor_mask5;
}

uint64_t g_hs_anchor_getX (uint64_t val)
{
    return (((val >> g_hs_anchor_bit1) - g_hs_anchor_zero) & g_hs_anchor_mask3) + 
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
static const uint64_t g_hs_mask2 = (1ULL << 30) - 1;
static const uint64_t g_hs_mask3 = (1ULL << 32) - 1;

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
    val = (((hs1 - x + g_hs_anchor_zero) & (g_hs_mask2))<< 20) + 
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

uint64_t g_hs_anchor_2Tile (uint64_t & anchor, /*uint64_t main_strand,*/ uint64_t revscomp_const)
{
    uint64_t strand = (anchor >> g_hs_anchor_bit2) & 1;
    /**
     * The cord of read (y) is shown in the same direction of the main strand
     * more readable 
     */
    //uint64_t y = (main_strand ^ strand) * revscomp_const - _nStrand(main_strand ^ strand) * g_hs_anchor_getY (anchor);
    /**
     * The cord of read read (y) is shown in its own direction.
     * easy for following processing (align)
     */
    uint64_t y = g_hs_anchor_getY(anchor);
	return (((anchor - (g_hs_anchor_zero << 20) + 
            ((anchor & g_hs_anchor_mask1)<< 20)) & 
              g_hs_anchor_mask2) & g_hs_anchor_mask1_) + y + (strand << 61);
}

int64_t tile_distance_x (uint64_t tile1, uint64_t tile2)
{
    return (int64_t)(get_cord_x(tile2)) - (int64_t)(get_cord_x(tile1));
}

int64_t tile_distance_y (uint64_t tile1, uint64_t tile2)
{
    return (int64_t)(_defaultTile.getY(tile2)) - (int64_t)(_defaultTile.getY(tile1));
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
 */
void set_tiles_flags_(String<uint64_t> & tiles, uint64_t sgn)
{
    for (int i = 0; i < length(tiles); i++)
    {
        set_cord_gap(tiles[i]);
        set_cord_recd(tiles[i] , sgn);
    }
}
/*
 * collecting minimizer
 *
 int g_mapHs_kmer_(String<Dna5> & seq, 
                         String<uint64_t> & g_hs, 
                         uint64_t start, 
                         uint64_t end, 
                         int g_hs_start, 
                         int step,  
                         uint64_t type)
{
    LShape(g_shape_len)  shape;
    hashInit(shape, begin(seq) + start);
    int count = 0; 
    int countx = 0;
    int i = 0; 
    uint64_t preX = 0;
    for (uint64_t k = start; k < end; k++)
    {
        if (++count == step)  //collecting every 10 bases
        {
            hashNext(shape, begin(seq) + k);
            if (shape.XValue == preX && countx < shape.span - shape.weight)
            {
                ++countx;
            }
            else
            {
                //TODO: k - getT(shape)
                g_hs_setGhs_(g_hs[g_hs_start + i++], shape.XValue, type, shape.strand, k);
                preX = shape.XValue;
                countx = 0;
            }
            count = 0;
        }
        else
        {
            hashNexth(shape, begin(seq) + k);
        }
    }
    return g_hs_start + i;
}
*/

/**
 * collecting k-mers in 'seq' to 'g_hs'
 */
 int g_mapHs_kmer_(String<Dna5> & seq, 
                   String<uint64_t> & g_hs, 
                   uint64_t start, 
                   uint64_t end, 
                   int g_hs_start, 
                   int step,  
                   uint64_t type)
{
    LShape shape(g_shape_len);
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
 * Stream part of 'g_hs' and convert the elements into anchors
 */
 int g_mapHs_setAnchors_ (String<uint64_t> & g_hs, 
                            String<uint64_t> & g_anchor,
                            int p1, 
                            int p2, 
                            int k, 
                            uint64_t revscomp_const,
                            int g_anchor_end) 
{
    unsigned n = 0;
    for (int i = p1; i < p2; i++) 
    {
        for (int j = p2; j < k; j++) 
        {
            g_hs_setAnchor_(g_anchor[g_anchor_end + n++], g_hs[i], g_hs[j], revscomp_const);
        }   
    }
    return g_anchor_end + n;
}
/*
 * select the longest anchor: faster but lose accuracy
 void g_mapHs_anchor_ (String<uint64_t> & anchor, 
                             String<uint64_t> & tile, 
                             int anchor_end, 
                             int thd_tileSize,
                             uint64_t  main_strand, 
                             int revscomp_const
                            )
{
    int64_t thd_min_segment = 100;
    int64_t prek = 0;
    int64_t prex = 0;
    int64_t prey = 0; 
    int anchor_len = 0, max_anchor_len = 0, max_prek = 0, max_k = 0;
    float thd_error_percent = 0.2;
    std::sort (begin(anchor), begin(anchor) + anchor_end);
    anchor[anchor_end] = ~0;
    for (int k = 0; k < anchor_end + 1; k++)
    {
        //TODO: handle thd_min_segment, anchor 
        //int64_t d = std::abs(_defaultACoord.getCoord(anchor[k]) - _defaultACoord.getCoord(anchor[prek]));
        int64_t d = std::abs((int64_t)g_hs_anchor_getY(anchor[k]) - (int64_t)g_hs_anchor_getY(anchor[prek]));
        //std::cout << "[]::g_mapHs_anchor_ 3 " << g_hs_anchor_getAnchor(anchor[k]) << " " << g_hs_anchor_getX(anchor[k]) << " " << g_hs_anchor_getY(anchor[k]) << " " << anchor_len << "\n";
        if (g_hs_anchor_getAnchor(anchor[k] - anchor[prek]) > 
            //150)
            thd_error_percent * std::max(thd_min_segment, d))
        {
            if (anchor_len > max_anchor_len)
            {
                max_anchor_len = anchor_len;
                max_prek = prek;
                max_k = k;
            }
            anchor_len = 0;
            prek = k;
        }
        else
        {
            anchor_len++;
        }
    }
        if (anchor_len > max_anchor_len)
        {
            max_prek = prek;
            max_k = anchor_end - 1;
        }
   
    
    std::sort (begin(anchor) + max_prek, 
            begin(anchor) + max_k, 
            [](uint64_t & s1, uint64_t & s2)
            {
                //return g_hs_anchor_getCord(s2) > g_hs_anchor_getCord(s1);
                return g_hs_anchor_getX(s2) > g_hs_anchor_getX(s1);
            });
    //appendValue(tile, g_hs_anchor_2Tile(anchor[max_prek], main_strand, revscomp_const));
    std::cout << "[]::g_mapHs_anchor_ " << max_prek + 1 << " " << max_k << " " << anchor_len << "\n";
    for (int j = max_prek; j < max_k; j++)
    {
        //TODO: change for invs, y is in descending order 
        if (g_hs_anchor_getX(anchor[j]) > prex + thd_tileSize
            || g_hs_anchor_getY(anchor[j]) > prey + thd_tileSize)
        {
        //std::cout << "[]::map_anchors " << j - 1<< " " << 10870 - g_hs_anchor_getY(anchor[j]) << " " << g_hs_anchor_getX(anchor[j]) << "\n";
            prex = g_hs_anchor_getX(anchor[j - 1]);
            prey = g_hs_anchor_getY(anchor[j - 1]);
            appendValue (tile, g_hs_anchor_2Tile(anchor[j - 1], main_strand, revscomp_const));
        }
    }
    if (prey != g_hs_anchor_getY(anchor[max_k -1]))
    {
        appendValue (tile, g_hs_anchor_2Tile(anchor[max_k - 1], main_strand, revscomp_const));
    }
}
*/

/*
 * cluster all the anchor candidates within the gap: 
 */
 void g_mapHs_anchor_ (String<uint64_t> & anchor, 
                       String<uint64_t> & tile, 
                       int anchor_end, 
                       int thd_tileSize,
                       uint64_t  main_strand, 
                       int revscomp_const
                       )
{
    int64_t thd_min_segment = 100;
    int64_t prek = 0;
    int64_t prex = -1;
    int64_t prey = -1; 
    int anchor_len = 0, max_anchor_len = 0, max_prek = 0, max_k = 0;
    int thd_k_in_window = 1;
    float thd_error_percent = 0.6;
    std::sort (begin(anchor), begin(anchor) + anchor_end);
    anchor[anchor_end] = ~0;
    for (int k = 0; k < anchor_end + 1; k++)
    {
        //TODO: handle thd_min_segment, anchor 
        int64_t d = std::abs((int64_t)g_hs_anchor_getY(anchor[k]) - (int64_t)g_hs_anchor_getY(anchor[prek]));
        if (g_hs_anchor_getAnchor(anchor[k]) - g_hs_anchor_getAnchor(anchor[prek]) > 
            thd_error_percent * std::max(thd_min_segment, d))
        {
            if ((std::abs(anchor_len / (float)(g_hs_anchor_getY(anchor[k - 1]) - g_hs_anchor_getY(anchor[prek]))) > g_thd_anchor_density && anchor_len > 2) || anchor_len > g_thd_anchor)
            {
                std::sort (begin(anchor) + prek, 
                    begin(anchor) + k, 
                    [](uint64_t & s1, uint64_t & s2)
                    {
                        return g_hs_anchor_getX(s2) > g_hs_anchor_getX(s1);
                    });
                prex = prey  = -1;
                for (int j = prek + 1; j < k; j++)
                {
                //TODO: change for invs, y is in descending order 
                    if ((g_hs_anchor_getX(anchor[j]) > prex + thd_tileSize
                        || g_hs_anchor_getY(anchor[j]) > prey + thd_tileSize))
                    {
                        prex = g_hs_anchor_getX(anchor[j - 1]);
                        prey = g_hs_anchor_getY(anchor[j - 1]);
                        appendValue (tile, g_hs_anchor_2Tile(anchor[j - 1],  revscomp_const));
                    }
                }
                appendValue (tile, g_hs_anchor_2Tile(anchor[k - 1], revscomp_const));
            }
            prek = k;
            anchor_len = 0;
        }
        else
        {
            anchor_len++;
        }
    }

    std::sort (begin(tile), 
               end(tile),
               [](uint64_t & s1, uint64_t & s2)
               {
                   return _defaultTile.getX(s1) < _defaultTile.getX(s2);
                });
}

/**
 * cluster all anchors and trim tiles for sv
 * Will conduct additional processing. 
 * Don't call it in the pipeline of approximate mapping.
 * 
 * ATTENTION: gr_start and gr_end is generated according to strand of the reference 
 * which is always regarded as the forward strand (strand = 0) rather than the main_strand
 */
 void g_mapHs_anchor_sv1_ (String<uint64_t> & anchor, 
                           String<uint64_t> & tiles, 
                           StringSet<FeaturesDynamic> & f1,
                           StringSet<FeaturesDynamic> & f2,
                           uint64_t gs_start,
                           uint64_t gs_end,
                           uint64_t gr_start,
                           uint64_t gr_end,
                           uint64_t  main_strand, 
                           uint64_t genomeId,
                           int anchor_end, 
                           int thd_tileSize,
                           int revscomp_const,
                           int direction
                           )
{
    int64_t thd_min_segment = 100;
    int thd_k_in_window = 1;
    int thd_fscore = 45;
    float thd_overlap_tile = thd_tileSize * 0.4;
    float thd_err_rate = 0.6;
    float thd_swap_tile = thd_tileSize * 0.05;
    
    /**
     * cluster anchors
     */
    int64_t prek = 0;
    int64_t prex = -1;
    int64_t prey = -1; 
    int anchor_len = 0, max_anchor_len = 0, max_prek = 0, max_k = 0;
    std::sort (begin(anchor), begin(anchor) + anchor_end);
    anchor[anchor_end] = ~0;
    String<int> score;

    for (int k = 0; k < anchor_end + 1; k++)
    {
        //TODO: handle thd_min_segment, anchor 
        
        int64_t d = std::abs((int64_t)g_hs_anchor_getY(anchor[k]) - (int64_t)g_hs_anchor_getY(anchor[prek]));
        if (g_hs_anchor_getAnchor(anchor[k]) - g_hs_anchor_getAnchor(anchor[prek]) > 
            thd_err_rate * std::max(thd_min_segment, d))
        {
            if ((std::abs(anchor_len / (float)(g_hs_anchor_getY(anchor[k - 1]) - g_hs_anchor_getY(anchor[prek]))) > g_thd_anchor_density && anchor_len > 2) || anchor_len > g_thd_anchor)
            {
                std::sort (begin(anchor) + prek, 
                    begin(anchor) + k, 
                    [](uint64_t & s1, uint64_t & s2)
                    {
                        return g_hs_anchor_getX(s2) > g_hs_anchor_getX(s1);
                    });
                prex = prey = -1;
                int kcount = 0;
                uint64_t first_low_bound_x = g_hs_anchor_getX(anchor[prek]) + window_size;
                uint64_t first_low_bound_y = g_hs_anchor_getY(anchor[prek]) + window_size;
                for (int i = prek + 1; g_hs_anchor_getX(anchor[i]) < first_low_bound_x && g_hs_anchor_getY(anchor[i]) < first_low_bound_y; i++)
                {
                    kcount++;
                }
                for (int j = prek + 1; j < k; j++)
                {
                    if ((g_hs_anchor_getX(anchor[j]) > prex + thd_tileSize ||  
                         g_hs_anchor_getY(anchor[j]) > prey + thd_tileSize))
                    {
                        prex = g_hs_anchor_getX(anchor[j - 1]);
                        prey = g_hs_anchor_getY(anchor[j - 1]);
                        appendValue (tiles, g_hs_anchor_2Tile(anchor[j - 1], 
                                     revscomp_const));
                        appendValue (score, kcount);
                        kcount=0;
                    }
                    else
                    {
                        kcount++;
                    }
                }
                appendValue (tiles, g_hs_anchor_2Tile(anchor[k - 1], 
                             revscomp_const));
                appendValue (score, kcount);
            }
            prek = k;
            anchor_len = 0;
        }
        else
        {
            anchor_len++;
        }
    }
    /**
     * remove poorly anchored tile: score < thd_tileSize 
     */
    float tz = thd_tileSize / 2;
    int prep = 0;
    for (int i = 0; i < length(score); i++)
    {
        
        if (score[i] < thd_k_in_window) 
        {
            continue;
        }
        else
        {
            uint64_t tile_x = _defaultTile.getX(tiles[i]);
            uint64_t tile_y = _defaultTile.getY(tiles[i]);
            unsigned fscore = _windowDist(f1[_defaultTile.getStrand(tiles[i])],
                                          f2[get_tile_id(tiles[i])],
                                          _DefaultCord.cord2Cell(tile_y), 
                                          _DefaultCord.cord2Cell(get_tile_x(tiles[i])));
            //if (fscore < windowThreshold)
            if (fscore < thd_fscore)
            {
                tiles[prep] = tiles[i];
                score[prep] = score[i];
                prep++;  
            }
        }
    }
    /**
     * sort cords according to the x, y and strand
     * using a weight function as
     * y + (x << w1) + (strand << w2);
     * or defined as y + 8 * x + 512 * strand
     * NOTE flip y according to the strand (for inversions)
     * NOTE strand is the relative strand to the main_strand.
     * The function first cluster according to the strand, then the reference direction, at last the read direction
     * For an example
     * Given s1 = 0, s2 = 1, y1 = 1, y2 = 10;
     * then:
     * when x1 - x2 > 512/8=64(bases), the cord1 > cord2 even if s1 < s2.  
     * when x2 < x1 < x2 + 64, the cord1 < cord2. This case is regarded as  x1 are not significantly bigger than x2, so the cord is sorted according to the strand.
     * The similar case for y 
     */
    std::sort (begin(tiles), 
               begin(tiles) + prep,
               [main_strand, revscomp_const](uint64_t & s1, uint64_t & s2)
               {
                    uint64_t strand1 = _defaultTile.getStrand(s1) ^ main_strand; // main_strand as the 0 strand
                    uint64_t strand2 = _defaultTile.getStrand(s2) ^ main_strand;
                    // _flip y to the main_strand if strand1 == 1, otherwise do nothing
                    uint64_t y1 = _flipCoord(_defaultTile.getY(s1), revscomp_const, strand1); 
                    uint64_t y2 = _flipCoord(_defaultTile.getY(s2), revscomp_const, strand2);
                   
                   return  y1 + (_defaultTile.getX(s1) << 3) + (strand1 << 9) < 
                           y2 + (_defaultTile.getX(s2) << 3) + (strand2 << 9);  

            });
    g_print_tiles_ (tiles, "gmas2");
    /**
     * remove tiles overlap with its adjacent tiles except the first and last tiles
     */
    int prep2 = 0;
    for (int i = 1; i < prep - 1; i++)
    {
        
        /**
         * the if condition can't be or '||', since for dels and invs edges of one side can be very close.
         * ---------------------- x
         *         /t1/\ t2\ 
         *        -----------     y
         * Example of ins above has large distance of y while x of two tiles are very close.
         */
        if (std::abs(int64_t(_defaultTile.getX(tiles[i + 1]) - _defaultTile.getX(tiles[prep2]))) < thd_overlap_tile && 
            std::abs(int64_t(_defaultTile.getY(tiles[i + 1]) - _defaultTile.getY(tiles[prep2]))) < thd_overlap_tile)
        {
            continue;
        }
        else
        {
            tiles[prep2] = tiles[i];
            prep2++;
        }
    }
    resize (tiles, prep2);
    
    /**
     * extend window if there are gaps between tiles until the horizontal coordinates x1 - x2 < window_size or the gap can't be extend any more
     * ATTENTION: relation between y1 and y2 currently are not considered.
     */
    ///extend the middle tiles
    for (int i = 1; i < length(tiles); i++)
    {
        i += extendPatch(f1, f2, tiles, i, tiles[i - 1], tiles[i], revscomp_const);   
    }
    
    ///extend the last and first tiles
    ///flip the coordinates from the direction of the reference genome to the direction of the data structure 'Cord'.
    uint64_t gr_start_flip = _flipCoord(gr_start, revscomp_const, main_strand);
    uint64_t gr_end_flip = _flipCoord(gr_end, revscomp_const, main_strand);
    if (main_strand)
    {
        std::swap (gr_start_flip, gr_end_flip);
    }
    uint64_t startCord = create_cord(genomeId, gs_start, gr_start_flip, main_strand);
    uint64_t endCord = create_cord(genomeId, gs_end, gr_end_flip, main_strand);
    if (empty(tiles))
    {
        extendPatch(f1, f2, tiles, 0, startCord, endCord, revscomp_const);
    }
    else
    {
        extendPatch(f1, f2, tiles, 0, startCord, tiles[0], revscomp_const);
        extendPatch(f1, f2, tiles, length(tiles), back(tiles), endCord, revscomp_const);   
    }
    
    //g_print_tiles_(tiles, f1, f2);
}

unsigned _get_tile_f_ (uint64_t & tile,
                  StringSet<FeaturesDynamic > & f1,
                  StringSet<FeaturesDynamic> & f2)
{

   // uint64_t tile = shift_tile(k_tile, )
    uint64_t tile_x = _defaultTile.getX(tile);
    uint64_t tile_y = _defaultTile.getY(tile);
    unsigned fscore = 
        _windowDist(f1[get_tile_strand(tile)], 
                    f2[get_tile_id(tile)],
                    _DefaultCord.cord2Cell(tile_y), 
                    _DefaultCord.cord2Cell(get_tile_x(tile)));
    return fscore;
}

unsigned _get_tile_f_tri_ (uint64_t & tile,
                  uint64_t & new_tile,
                  StringSet<FeaturesDynamic > & f1,
                  StringSet<FeaturesDynamic > & f2, 
                  unsigned thd_accept_score,
                  int block_size = window_size
                  )
{
    int shift = window_size / 4;
    uint64_t tile_l = shift_tile(tile, -shift, -shift);
    uint64_t tile_r = shift_tile(tile, shift, shift);
    unsigned fscore =  _get_tile_f_ (tile, f1, f2) ;
    if (fscore < thd_accept_score)
    {
        new_tile = tile;
        return fscore;
    }
    else
    {
        fscore = std::min(_get_tile_f_(tile_l, f1, f2), fscore);
        if (fscore < thd_accept_score)
        {
            new_tile = tile_l;
            return fscore;
        }
        else
        {
            new_tile = tile_r;
            return std::min(_get_tile_f_(tile_r, f1, f2), fscore);

        }
    }
}
//!!TODO::[Important] tune the parm, especially the thd_err_rate 
struct MapAnchorParm
{

};
/**
 * cluster all anchors and trim tiles for sv
 * Will conduct additional processing. 
 * Don't call it in the pipeline of approximate mapping.
 */ 
 int g_mapHs_anchor_sv2_ (String<uint64_t> & anchor, 
                          String<uint64_t> & tiles, 
                          StringSet<FeaturesDynamic> & f1,
                          StringSet<FeaturesDynamic> & f2,
                          uint64_t gap_str,
                          uint64_t gap_end,
                          int anchor_end, 
                          int thd_tileSize,
                          int revscomp_const,
                          int direction
                          )
{
    CmpInt64 g_cmpll;
    int block_size = window_size;
    int64_t thd_min_segment = 100;
    int thd_pattern_in_window = 1;
    unsigned thd_fscore = get_windowThreshold(f1);
    float thd_overlap_tile = thd_tileSize * 0.4;
    float thd_err_rate = 0.2;
    float thd_swap_tile = thd_tileSize * 0.05;
    uint64_t main_strand = get_cord_strand(gap_str);
    
    //step 1. cluster anchors cord_end
    int64_t prek = 0;
    int64_t prex = -1;
    int64_t prey = -1; 
    int64_t tmp_shift = block_size >> 1;
    int anchor_len = 0;
    std::sort (begin(anchor), begin(anchor) + anchor_end);
    anchor[anchor_end] = ~0;
    int pre_tile_end = 0;

    for (int k = 0; k < anchor_end + 1; k++)
    {
        //TODO: handle thd_min_segment, anchor 
        int64_t d = std::abs((int64_t)g_hs_anchor_getY(anchor[k]) - (int64_t)g_hs_anchor_getY(anchor[prek]));
        if (g_hs_anchor_getAnchor(anchor[k]) - g_hs_anchor_getAnchor(anchor[prek]) > 
            thd_err_rate * std::max(thd_min_segment, d))
        {
            int thd_anchor_accpet = g_thd_anchor_density * 
            std::abs(int64_t(g_hs_anchor_getY(anchor[k - 1]) - 
                             g_hs_anchor_getY(anchor[prek])));
            thd_anchor_accpet = std::max (thd_anchor_accpet, 2);
            thd_anchor_accpet = std::min (g_thd_anchor, thd_anchor_accpet);

            if (anchor_len > thd_anchor_accpet) 
            {
                std::sort (begin(anchor) + prek, 
                           begin(anchor) + k, 
                [](uint64_t & s1, uint64_t & s2)
                {return g_hs_anchor_getX(s2) > g_hs_anchor_getX(s1);
                });
                //prex = prey = - 1;
                prex = g_hs_anchor_getX(anchor[prek]);
                prey = g_hs_anchor_getY(anchor[prek]);
                int kcount = 0;

                uint64_t upper_x = g_hs_anchor_getX(anchor[prek]) + thd_tileSize;
                uint64_t upper_y = g_hs_anchor_getY(anchor[prek]) + thd_tileSize;
                int64_t centroid_x = 0;
                int64_t centroid_y = 0;
                int prej = prek;
                for (int j = prek; j < k; j++)
                {
                    uint64_t x = g_hs_anchor_getX(anchor[j]);
                    uint64_t y = g_hs_anchor_getY(anchor[j]);
                    if ((x > prex + thd_tileSize || y > prey + thd_tileSize)||
                         j == k - 1)
                    {
                        if (j != k - 1)
                        {
                            centroid_x /= (j - prej);
                            centroid_y /= (j - prej);
                        }
                        else
                        {
                            //centroid_x += g_hs_anchor_getX(anchor[j]);
                            //centroid_y += g_hs_anchor_getY(anchor[j]);
                            centroid_x /= (j - prej);
                            centroid_y /= (j - prej);
                        }
                        g_cmpll.min(tmp_shift, block_size >> 1) << centroid_x 
                                                                << centroid_y; 
                        int64_t tmp_tile_x = centroid_x - tmp_shift ;
                        int64_t tmp_tile_y = centroid_y - tmp_shift ;
                        uint64_t tmp_tile = create_tile(0, tmp_tile_x, tmp_tile_y, 
                                            g_hs_anchor_get_strand(anchor[prek]));
                        uint64_t new_tile;
                        unsigned score = _get_tile_f_tri_(tmp_tile, new_tile, f1, f2, thd_fscore, thd_tileSize);
                        if (kcount >= thd_pattern_in_window && 
                            score < thd_fscore)
                        {
                            if (empty(tiles) || is_tile_end(back(tiles)))
                            {
                                set_tile_start(new_tile);
                                appendValue (tiles, new_tile);
                            }
                            else
                            {
                                appendValue (tiles, new_tile);
                            }
                        }
                        prex = g_hs_anchor_getX(anchor[j - 1]);
                        prey = g_hs_anchor_getY(anchor[j - 1]);
                        prej = j;
                        centroid_x = g_hs_anchor_getX(anchor[j]);
                        centroid_y = g_hs_anchor_getY(anchor[j]);
                        kcount=1;
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
                pre_tile_end = length(tiles);
            }
            prek = k;
            anchor_len = 0;
        }
        else
        {
            anchor_len++;
        }
    }
    //step 2. clip and merge 
    String<uint64_t> tmp_tiles = tiles;
    std::sort (begin(tmp_tiles), end(tmp_tiles),
    [](uint64_t & s1, uint64_t & s2)
    {  

      //  return _defaultTile.getX(s1) < _defaultTile.getX(s2);
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
 * ATTENTION: relation between y1 and y2 currently are not considered.
 */
    int thd_overlap_size = 170;
    int thd_gap_size = 180;
    uint64_t cord_str = gap_str;
    uint64_t cord_end = gap_end;
    if (main_strand)
    {
        cmpRevCord(cord_str, cord_end, cord_str, cord_end, revscomp_const);
    }
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
                i += new_num;
                set_tile_end(tiles[i]);
            }
        }
        if (i > 1 && !is_tile_end (tiles[i - 1]) && !is_tile_start(tiles[i]))
        {
            i += extendPatch(f1, f2, tiles, i, tiles[i - 1], tiles[i], revscomp_const, thd_overlap_size, thd_gap_size);   
        }
    }
    int ct_tile_block = 0;
    for (int i = 0; i < length(tiles); i++)
    {
        if (is_tile_end(tiles[i]))
        {
            ++ct_tile_block;
        }
    }
    if (ct_tile_block > 1)
    {
        for (int i = 0; i < length(tiles); i++)
        {
            if (is_tile_end(tiles[i]))
            {
                set_cord_end(tiles[i]);
            }
        }
    }
    return 0;
}
/**
 * wrapper
 */
 void g_mapHs_anchor_sv_ (String<uint64_t> & anchor, 
                                String<uint64_t> & tiles, 
                                StringSet<FeaturesDynamic > & f1,
                                StringSet<FeaturesDynamic > & f2,
                                uint64_t cord_str,
                                uint64_t cord_end,
                                int anchor_end, 
                                int thd_tileSize,
                                int revscomp_const,
                                int direction
                                )
{
    g_mapHs_anchor_sv2_ (anchor, tiles, f1, f2, cord_str, cord_end, anchor_end, thd_tileSize, revscomp_const, direction);
}

/**
 * Map gaps specified by the @cord_str and @cord_end and create the chain of tiles 
 * to cover the gap as long as possible.
 * @cord_str will be extended towards the right side, @cord_end will 
 * be extended towards the left side.
 * Cords between them will be extended towards both sides.
 * @f1 and @f2 are frequency vector for approximate mapping
 */
 int g_mapHs_(String<Dna5> & seq1, //genome
              String<Dna5> & seq2, //read
              String<Dna5> & comstr,
              String<uint64_t> & g_hs,
              String<uint64_t> & g_hs_anchor,
              String<uint64_t> & g_hs_tile,    //results
              StringSet<FeaturesDynamic > & f1,  
              StringSet<FeaturesDynamic >& f2,
              uint64_t cord_str,
              uint64_t cord_end, 
              int thd_tileSize,   //WARNING 192 not allowed to change.
              int direction = g_align_closed
             )
{
    int g_hs_end = 0;
    int g_hs_anchor_end = 0;
    uint64_t rvcp_const = length(seq2) - 1;
    uint64_t gs_str = get_cord_x(cord_str);
    uint64_t gs_end = get_cord_x(cord_end);
    uint64_t gr_str = get_cord_y(cord_str);
    uint64_t gr_end = get_cord_y(cord_end);

    if (get_cord_strand(cord_str))
    {
        gr_str = rvcp_const - gr_str;
        gr_end = rvcp_const - gr_end;
        std::swap (gr_end, gr_str);
    }
    g_hs_end = g_mapHs_kmer_(seq1, g_hs, gs_str, gs_end, g_hs_end, 10, 0);
    g_hs_end = g_mapHs_kmer_(seq2, g_hs, gr_str, gr_end, g_hs_end, 1, 1);

    std::sort (begin(g_hs), begin(g_hs) + g_hs_end);

    int p1 = 0, p2 = 0;
    for (int k = 1; k < g_hs_end; k++)
    {    
        switch (g_hs_getXT(g_hs[k] ^ g_hs[k - 1]))
        {
            case 0:       //x1 = x2 both from genome or read
                break;
            case 1:       //x1 = x2 one from genome the other from read
                p2 = k;
                break;
            default:      //anchor current block before process next block 
                g_hs_anchor_end = g_mapHs_setAnchors_(g_hs, g_hs_anchor, p1, p2, k, rvcp_const, g_hs_anchor_end);
                p1 = k;
                p2 = k; 
        }
    }
    g_mapHs_anchor_sv_(g_hs_anchor, g_hs_tile, f1, f2, 
                       cord_str, cord_end,
                       g_hs_anchor_end, 
                       thd_tileSize,
                       rvcp_const,
                       direction
                      );
}

 int64_t g_anchor_dx_(uint64_t val1, uint64_t val2)
{
    return int64_t (g_hs_anchor_getX(val2) - g_hs_anchor_getX(val1));
}
 int64_t g_anchor_da_(uint64_t val1, uint64_t val2)
{
    return int64_t (g_hs_anchor_getAnchor(val2) - g_hs_anchor_getAnchor(val1));
}

///c_ functions to clip breakpoints by counting kmers
const unsigned c_shape_len = 8; //anchor shape
const unsigned c_shape_len2 = 4; //base-level clipping shape
const unsigned c_shape_len3 = 4; //base-level clipping gap shape

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
 * Stream a block of 'g_hs' within [p1,p2)x[p2,k), and convert the combination of elements into anchors with the following restrictions
 * |candidates_anchor - 'anchor' | < band
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
 * ---------mmmmmmmmmm
 */
inline int64_t c_sc_(int val1, int val2, 
                     float thd_exp_err = 0.85)
{
    float rate = (float)val1 / val2;
    float thd_err1 = thd_exp_err;
    float thd_err2 = thd_exp_err; 
    int thd_min_len = 10;
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
    if (rate > thd_err1 && val1 > thd_min_len)
    {
        return val1 << 2; 
    }
    else if (rate > thd_err2 && val1 > thd_min_len)
    {
        return val1 << 1;
    }
    else 
    {
        return val1;
    }
}
/**
 * 
 * Clip the anchors when it reach the first anchor that can be well extended. 
 * @clip_direction: -1 gap-match; 1 match-gap
 */
int64_t c_clip_anchors_ (String<uint64_t> & anchor, 
                         uint64_t gs_start,
                         uint64_t gr_start,
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
    if (anchor_end < 1)
    {
        return 0;
    }
    int direction = (clip_direction < 0) ? -1 : 1;
    int it = 0;
    int bit1 = 20;
    int bit2 = g_hs_anchor_bit1 + bit1;
    uint64_t mask = (1LL << bit1) - 1;
    uint64_t ct_conts = 0;
    anchor[anchor_end] = ~0;
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
            uint64_t x = g_hs_anchor_getX(anchor[i_str]) - gs_start;
            uint64_t y = g_hs_anchor_getY(anchor[i_str]) - gr_start;
            anchor[it++] = (x << bit2) + ((ct_conts + 1) << g_hs_anchor_bit1) + y;
            ct_conts = 0;
            //collect continuos patterns (no gaps)
        }
    }
    i_str = anchor_end - 1 - ct_conts;
    uint64_t x = g_hs_anchor_getX(anchor[i_str]) - gs_start;
    uint64_t y = g_hs_anchor_getY(anchor[i_str]) - gr_start;
    anchor[it++] = (x << bit2) + ((ct_conts + 1) << g_hs_anchor_bit1) + y;

    if (it < 1)
    {
        return 0;
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
                return (rslt_x << 32) + rslt_y;
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
            return (rslt_x << 32) + rslt_y;
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
                return (rslt_x << 32) + rslt_y;
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
            return (rslt_x << 32) + rslt_y;
        }
    }

    return 0;
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
int c_isGapMatch_2anchor_(uint64_t & anchor, uint64_t & x, uint64_t & y, uint64_t & dv, short& t1, short & t2, short & l1, short & l2, short k)
{
    uint64_t type = c_isGapMatch_(dv, t1, t2, l1, l2, k);
    if (type)
    {
        return type;
    }
    else
    {
        c_2GapAnchor_(anchor, x, y, 3ULL);
        return type;
    }
}
/**
 * Aborted
 */
uint64_t c_clip_anchors_precise(String<uint64_t> & anchor, 
                             uint64_t gs_start,
                             uint64_t gr_start,
                             int anchor_end,
                             unsigned shape_len, 
                             int thd_merge1, // thd of anchor
                             int thd_merge1_lower,
                             int thd_merge2, //thd of x
                             int thd_width, //WARNING: elements to search. Not band of bases
                             int thd_clip_sc = c_sc_(27, 30)
                             //uint64_t  main_strand, 
                             //int revscomp_const
                            )
{
    int last_k = 0, start_k = 0;
    uint64_t ct_conts = 0;
    uint64_t thd_break_a = 2; ///WARNING::tune 
    anchor[anchor_end] = ~0;
    int it = 0;
    int bit1 = 20;
    int bit = g_hs_anchor_bit1 + bit1;
    int64_t mask = (1LL << bit) - 1;
    int64_t mask1 = (1LL << bit1) - 1;
    int64_t mask2 = (mask1 << bit1);
    std::sort (begin(anchor), begin(anchor) + anchor_end);
    for (int k = 0; k < anchor_end; k++)
    {
        uint64_t dx = std::abs(g_anchor_dx_(anchor[k], anchor[k + 1]));
        uint64_t dy = g_hs_anchor_getY(anchor[k + 1] - anchor[k]);
        int64_t da = g_anchor_da_(anchor[k], anchor[k + 1]);
        if (da == 0 && dy == 1) 
        {
            ct_conts++;
        }
        else
        {
            last_k = k - ct_conts ;
            uint64_t x = g_hs_anchor_getX(anchor[last_k]) - gs_start;
            uint64_t y = g_hs_anchor_getY(anchor[last_k]) - gr_start;
            anchor[it++] = (x << bit) + ((ct_conts + 1) << g_hs_anchor_bit1) + y;
            start_k = k + 1;
            ct_conts = 0;
        }
    }
    std::sort(begin(anchor), begin(anchor) + it, std::greater<uint64_t>());

    int64_t itgenome = 0;
    int64_t itread = 0;
    for (int i = 0; i < it - 1; i++)
    {
        int64_t y0 = g_hs_anchor_getY(anchor[i]);
        int64_t x0 = (anchor[i] >> bit) & mask1;
        int c_conts = ((anchor[i] >> bit1) & mask1) + shape_len - 1;
        int nxt_c_conts = c_conts;
        int sc = 0;
        int64_t nxt_end = x0 + c_conts;
        int sum_conts = c_conts;
        int seg_len = 0;
        int dj = 0;
        int width_count = 0;
        int64_t x1;
        int64_t y1;
        int64_t d_anc;
        int64_t x_end = x0;
        int64_t y_end = y0;
        for (int j = i + 1; j < it - 1 && width_count < thd_width; j++)
        {
            y1 = g_hs_anchor_getY(anchor[j]);
            x1 = (anchor[j] >> bit) & mask1;
            d_anc = int64_t(x0 - x1 - y0 + y1);
            c_conts = ((anchor[j] >> bit1) & mask1) + shape_len - 1;
            nxt_end = x1 + c_conts;
            nxt_c_conts = c_conts;
            //if (std::abs(d_anc) < std::max(int(x1 - x0) >> thd_merge1, thd_merge1_lower) && (x1 - pre_end) < thd_merge2 && x1 > x0)
            if (std::abs(d_anc) < std::max(int(x0 - x1) >> thd_merge1, thd_merge1_lower) && (x0 - nxt_end) < thd_merge2 && x0 > x1)
            {
                sum_conts += c_conts;
                seg_len += x0 - x1;
                sc += c_sc_(nxt_c_conts, x0 - x1);
                x0 = x1; y0 = y1;
                x_end = x1;
                y_end = y1;
                ++dj;
            }
            else
            {
                anchor[j - dj] = anchor[j];
                ++width_count = 0;
            }
        }
        sc += c_sc_(nxt_c_conts, nxt_c_conts);
        it -= dj;
        if (sum_conts > 15)
        {
            itgenome = x_end;
            itread = y_end;
        }
    } 
    return (itgenome << 32) + itread;
//    return 0;
    ///clip by sliding windows
    /*
    int w_len = 5;
    int w_ct_l = 0;
    int w_ct_r = 0;
    float max_d = 0;
    uint64_t clip = 0;
    for (int i = max_prek; i < max_prek + w_len; i++)
    {
        w_ct_l += g_anchor_dx_ (anchor[i], anchor[i + 1]);
        w_ct_r += g_anchor_dx_ (anchor[i + w_len], anchor[i + w_len + 1]);
    }
    for (int i = max_prek + w_len; i < max_k - w_len - 1; i++)
    {
        w_ct_l += g_anchor_dx_(anchor[i], anchor[i + 1]) 
                - g_anchor_dx_(anchor[i - w_len], anchor[i - w_len + 1]);
        w_ct_r += g_anchor_dx_(anchor[i + w_len], anchor[i + w_len + 1])
                - g_anchor_dx_(anchor[i], anchor[i + 1]);
        float d = ((float)w_ct_l - w_ct_r) / (w_ct_l * w_ct_r);
        if (d > max_d)
        {
            max_d = d;
            clip = anchor[i];
        }
    }
    return g_hs_anchor_getAnchor(max_anchor_median);
    */
}

String<Dna5> & int2str(String<Dna5> & s, uint64_t s1, int n)
{
    int mask = 3;
    for (int i = 0; i < n; i++)
    {
        int t = (s1 >> (i * 2)) & mask;
        if (t == 0)
        {
            s[n - 1 - i] = 'A';
        }
        if (t == 1)
        {
            s[n - 1 - i] = 'C';
        }
        if (t == 2)
        {
            s[n - 1 - i] = 'G';
        }
        if (t == 3)
        {
            s[n - 1 - i] = 'T';
        }
    }
    return s;
}

//<<debug util
String<char>  hash2kmer(uint64_t val)
{
    String<char> s;
    for (int n = 64; n >= 1; n/=4)
    {
        int c = val / n;
        val -= c * n; 
        char a;
        if (c == 0)
            a = 'A';
        if (c == 1)
            a = 'C';
        if (c == 2)
            a = 'G';
        if (c == 3)
            a = 'T';
        appendValue(s, a);
        //std::cerr << val << " " << c << " " << n<< "\n";
    }
        //std::cerr << s << "\n";
    return s;
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
    uint64_t g_anchor_val = 0;
    g_anchor_val = c_clip_anchors_(g_anchor, 
                                   gs_str, 
                                   gr_str, 
                                   g_anchor_end, 
                                   c_shape_len, 
                                   thd_merge1, 
                                   thd_merge1_lower, 
                                   thd_merge2, 
                                   clip_direction);    
    
    int64_t dx = (g_anchor_val >> 32);
    int64_t dy = (g_anchor_val & ((1ULL << 32) - 1));
    uint64_t clip = create_cord(genomeId, gs_str + dx, gr_str + dy, gr_strand);

//step2. clip_extend gap pattern further.
    int extend_window = 100;
    int thd_ovlp_shift = 10;
    int thd_merge_anchor = 5;
    int thd_merge_drop = 6;
    int thd_error_level = 3;  // >>3 == * 0.125
    int thd_scan_radius = 3;  //at least scan 5 elements in the genome for each kmer in the read
    uint64_t extend_str;
    uint64_t extend_end;
    uint64_t breakpoints;
    int64_t shift;
    if (isClipTowardsLeft (clip_direction))
    {
        //<<<debug
        g_cmpll.min(dx, dx + thd_ovlp_shift) << get_cord_x(clip_end - clip_str) - 1;
        g_cmpll.min(dy, dy + thd_ovlp_shift) << get_cord_y(clip_end - clip_str) - 1;
        extend_end = shift_cord (clip_str, dx, dy);
        //>>>debug
  //      extend_end = shift_cord (clip_str, dx, dy);
        g_cmpll.min(shift, extend_window) 
                    << get_tile_x(extend_end - clip_str)
                    << get_tile_y(extend_end - clip_str);
        extend_str = shift_cord (extend_end, -shift, -shift);
    }
    else if (isClipTowardsRight (clip_direction))
    {
        g_cmpll.max(dx, dx - thd_ovlp_shift) >> 0;
        g_cmpll.max(dy, dy - thd_ovlp_shift) >> 0;
        extend_str = shift_cord(clip_str, dx, dy);
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

/**
 * cord_str and cord_end required to have same strand.
 * TODO::process multiple tile_start and tile_end. Currently it handles tiles having one block.
 */
int g_extend_clip_(String<Dna5> & seq1,
                   String<Dna5> & seq2,
                   String<Dna5> & comstr, //complement reverse of the read (seq2)
                   String<uint64_t> & tiles,
                   String<uint64_t> & clips,
                   String<uint64_t> & g_hs,
                   String<uint64_t> & g_hs_anchor,
                   uint64_t cord_str,
                   uint64_t cord_end, 
                   int direction,
                   int thd_cord_gap
                  )
{
    CmpInt64 g_cmpll;
    if (get_cord_strand(cord_str ^ cord_end))
    {
        return 0;
    }
    uint64_t head_tile = cord_str;
    uint64_t tail_tile = cord_end; 
    uint64_t main_strand = get_cord_strand(cord_str);
    int64_t tile_size = window_size;
    int64_t thd_max_gap_size = tile_size;
    float thd_band_ratio = 0.5;
    insert (tiles, 0, head_tile);
    appendValue (tiles, tail_tile);
    String<int> sv_flags;
    resize (sv_flags, length(tiles));
    for (int i = 0; i < length(sv_flags); i++)
    {
        sv_flags[i] = 0;
    }    
    int sv_exists = 0;
    for (int i = 1; i < length(tiles); i++)
    {
        int64_t dx = tile_distance_x(tiles[i - 1], tiles[i]);
        int64_t dy = tile_distance_y(tiles[i - 1], tiles[i]); 
        ///check sv type
        if (_defaultTile.getStrand(tiles[i] ^ tiles[i - 1])) //inv
        {
            if (get_tile_strand(tiles[i]) ^ main_strand)
            {
                sv_flags[i] |= g_sv_inv + g_sv_l + g_sv_extend;
            }
            else
            {
                sv_flags[i - 1] |= g_sv_inv + g_sv_r + g_sv_extend;
            }
            sv_exists = 1;
        }
        else if (dx - dy > tile_size)   //del
        {
            sv_flags[i - 1] |= g_sv_ins + g_sv_r + g_sv_shrink;
            sv_flags[i] |= g_sv_ins + g_sv_l + g_sv_shrink;
            sv_exists = 1;
        }
        else if (dy - dx > tile_size)  //ins 
        {
            sv_flags[i - 1] |= g_sv_del + g_sv_r + g_sv_shrink;
            sv_flags[i] |= g_sv_del + g_sv_l + g_sv_shrink;
            sv_exists = 1;
        }
        else if (dy > thd_cord_gap && dx > thd_cord_gap &&
                 dx > 0 && dy > 0) //gap, inv
        {
            sv_flags[i - 1] += g_sv_gap + g_sv_r + g_sv_shrink;
            sv_flags[i] += g_sv_gap + g_sv_l + g_sv_shrink;
            sv_exists = 1;
        }
    }

    if (direction == g_align_left) 
    {
        if (get_tile_y(tiles[0]) < thd_cord_gap && 
            ! get_tile_strand(tiles[0] ^ tiles[1]))
        {
            sv_flags[1] &= ~(g_sv_l);
        }
        sv_flags[0] = 0;
    }
    else if (direction == g_align_right)
    {
        uint64_t tile1 = back(tiles);
        uint64_t tile2 = tiles[length(tiles) - 2];
        if (get_tile_y(tile1 - tile2) < thd_cord_gap &&
            !get_tile_strand(tile1 ^ tile2))
        {
            sv_flags[length(sv_flags) - 2] &= ~(g_sv_r);
        }
        back(sv_flags) = 0;
    }
    
    uint64_t clip = -1;
    uint64_t cgstart, cgend;
    uint64_t crstart, crend;
    uint64_t crstrand, band;
    uint64_t delta;
    uint64_t clip_str, clip_end;
    int64_t shift;
    String<uint64_t> clips_tmp;
    int clip_direction;
    if (sv_exists)
    {
        for (int i = 0; i < length(sv_flags); i++)
        {
            uint64_t tile = tiles[i];
            uint64_t clip1 = EmptyClipConst;
            uint64_t clip2 = EmptyClipConst;
            if (sv_flags[i] & g_sv_r)
            {
                g_cmpll.min(shift) << int64_t(tile_size / 2)
                                   << length(seq1) - 1 - get_tile_x(tile) 
                                   << length(seq2) - 1 - get_tile_y(tile);
                clip_str = shift_tile(tile, shift, shift);
                g_cmpll.min(shift, shift + thd_max_gap_size) 
                                   << length(seq1) - 1 - get_tile_x(tile) 
                                   << length(seq2) - 1 - get_tile_y(tile);
                clip_end = shift_tile(tile, shift, shift);
                clip_direction = 1;
                thd_band_ratio = 0.5;
                clip = c_clip_ (seq1, seq2, comstr, clip_str, clip_end, g_hs, g_hs_anchor, thd_band_ratio, clip_direction); 
                clip1 = clip;
                remove_tile_sgn(clip1);
                appendValue(clips_tmp, clip1);
            }
            if (sv_flags[i] & g_sv_l)
            {
                g_cmpll.min(shift) << tile_size / 2 
                                   << length(seq1) - 1 - get_tile_x(tile) 
                                   << length(seq2) - 1 - get_tile_y(tile);
                clip_end = shift_tile(tile, shift, shift);
                g_cmpll.min(shift) << tile_size / 2
                                   << get_tile_x(tile)
                                   << get_tile_y(tile);
                clip_str = shift_tile(tile, -shift, -shift);
                clip_direction = -1;
                thd_band_ratio = 0.5;
                clip = c_clip_ (seq1, seq2, comstr, clip_str, clip_end, g_hs, g_hs_anchor, thd_band_ratio, clip_direction); 
                clip2 = clip;
                remove_tile_sgn(clip2);
                appendValue(clips_tmp, clip2);
            }
        }
        //Insert clip and str and end to the clips records list.
        //Using the brackets matching methods for clips of same strands as the following
        //[{(< >)}]   [] {} () <> indicates pairs of clip str and end.
        std::sort (begin(clips_tmp), end(clips_tmp));
        int end0 = 0;
        for (int i = 0; i < length(clips_tmp); i++)
        {
            if (!get_tile_strand(clips_tmp[i]))
            {
                ++end0;
            }
        }
        if (length(clips_tmp) / 2 * 2 == length(clips_tmp))
        {
        //TODO::change the if 
        //WARNING::it's incorrect, just for temporarily use.
        for (int i = 0; i < end0 / 2; i++)
        {
            insertClipStr(clips, clips_tmp[i]);
            insertClipEnd(clips, clips_tmp[end0 - 1 - i]);
        }
        for (int i = end0; i < (end0 + length(clips_tmp)) / 2 ; i++)
        {
            insertClipStr(clips, clips_tmp[i]);
            insertClipEnd(clips, clips_tmp[length(clips_tmp) + end0 - 1 - i]);
        }
        }
    }
    if (length(tiles) > 2)
    {
        for (int i = 0; i < length(tiles) - 2; i++)
        {
            tiles[i] = tiles[i + 1];
        }   
        resize (tiles, length(tiles) - 2);
    }
    else
    {
        clear(tiles);
    }
    return 0;
    /// remove the inserted head and tail tiles;
}


/**
 * 1.Map gaps between given cords [cord1, cord2] as follows
 *  |cord1---------cord2|
 *  |--------gap--------|
 *  Both cord1 and cord2 will be extend towards right or left at each side.
 * 
 * 2. ATTENTION cord1 and cord2 are required to have the same strand.
 * case 1:
 *   seg1  gap  seg2
 *  -----xxxxxxx-----  correct
 * 
 * case 2:
 *   seg1  gap  seg2
 *  +++++xxxxxxx-----  incorrect : different strands of seg1 and seg2
 */
 int mapGap_ (StringSet<String<Dna5> > & seqs,
              String<Dna5> & read,
              String<Dna5> & comstr,
              uint64_t cord1, 
              uint64_t cord2, 
              String<uint64_t> & g_hs,
              String<uint64_t> & g_anchor,
              StringSet<FeaturesDynamic > & f1,
              StringSet<FeaturesDynamic >& f2,
              String<uint64_t> & tiles,     //results
              String<uint64_t> & clips,     //results 
              int direction,
              int const thd_cord_gap, 
              int const thd_tileSize,
              int thd_cord_remap
             )
{
    clear(tiles);
    if (_DefaultCord.getCordStrand(cord1 ^ cord2))
    {
        return -1;
    }
    uint64_t strand = _DefaultCord.getCordStrand (cord1);
/// WARNING: the main strand is defined as the strand of the cord1
/// ATTENTION: gr_start and gr_end are flipped !!! if the main strand is complement reversed
    uint64_t genomeId = get_cord_id(cord1);
    if (get_cord_x(cord2 - cord1) > thd_cord_remap && 
        get_cord_y(cord2 - cord1) > thd_cord_remap)
    {
        g_mapHs_(seqs[genomeId], read, comstr,
                 g_hs, g_anchor, tiles, f1, f2,
                 cord1, cord2,
                 thd_tileSize,
                 direction
               );

        g_extend_clip_(seqs[genomeId], read, comstr, tiles, clips, g_hs, g_anchor, cord1, cord2, direction, thd_cord_gap);
    }
    /*
    else 
    {
        appendValue(tiles, cord1);
        appendValue(tiles, cord2);
        g_extend_clip_(seqs[genomeId], read, comstr, tiles, clips, g_hs, g_anchor, gs_start, gs_end, gr_start, gr_end, strand, genomeId, direction, thd_cord_gap);
    }
    */
    return length(tiles);
}

/**
 * Re-map gaps in cords.
 * Gaps at the beginning or end are also remapped.
 */
int mapGaps(StringSet<String<Dna5> > & seqs, 
            String<Dna5> & read, 
            String<Dna5> & comstr,
            String<uint64_t> & cords, 
            String<uint64_t> & g_hs,
            String<uint64_t> & g_anchor,
            String<uint64_t> & clips, // string for clips cords
            StringSet<FeaturesDynamic > & f1,
            StringSet<FeaturesDynamic >& f2,
            int const thd_gap, 
            int const thd_tileSize
           )
{
    CmpInt64 g_cmpll;
    if  (length(cords) <= 1)
    {
        return 0;
    }
    String <uint64_t> tiles;
    int64_t shift_x;
    int64_t shift_y;
    uint64_t genomeId;
    uint64_t gap_len = 0;
    uint64_t thd_max_extend = 10000;
    int block_size = window_size;
    int thd_cord_remap = 100;
    int thd_cord_gap = thd_gap + block_size;
    ///NOTE cords[0] is the head cord, so starts from 1
    for (unsigned i = 1; i < length(cords); i++)
    {
        unsigned sid = get_cord_id(cords[i]);
        uint64_t cordx1 = get_cord_x(cords[i - 1]);
        uint64_t cordy1 = get_cord_x(cords[i - 1]);
        uint64_t cordx2 = get_cord_x(cords[i]);
        uint64_t cordy2 = get_cord_x(cords[i]);
        int64_t dcordx = cordx2 - cordx1;
        int64_t dcordy = cordy2 - cordy1;
        if (_DefaultCord.isBlockEnd(cords[i - 1]))  ///left clip first cord
        {
            uint64_t cord2 = cords[i];
            if (get_cord_y(cord2) > thd_cord_gap) 
            {            
                shift_x = std::min(thd_max_extend, get_cord_x(cord2));
                shift_y = std::min(thd_max_extend, get_cord_y(cord2));
                uint64_t infi_cord = shift_cord(cord2, -shift_x, -shift_y);   
                uint64_t cord1 = infi_cord;
                mapGap_ (seqs, read, comstr, 
                         cord1, cord2, 
                         g_hs, g_anchor, f1, f2,  
                         tiles, clips, 
                         g_align_left, 
                         thd_cord_gap, 
                         thd_tileSize,
                         thd_cord_remap);
                if (length(tiles) > 0)
                {
                    uint64_t recd = get_cord_recd(cords[i]);//set cord flag
                    //dout << "recd" << get_cord_y(tiles[0]) << get_cord_x(tiles[0]) << recd << "\n";
                    set_tiles_flags_ (tiles, recd);
                    insert(cords, i, tiles);
                    i += length(tiles);
                    clear(tiles);
                }
                gap_len += get_cord_y(cord2);
            }
            continue;
        }
        ///clip closed interval for middle cords
        else if ((cordx2 - cordx1 > thd_cord_gap  ||
                  cordy2 - cordy1 > thd_cord_gap) && 
                  cordx1 < cordx2 &&
                  cordy1 < cordy2 &&
                  dcordx < thd_max_extend && 
                  dcordy < thd_max_extend
                )        
        {
            g_cmpll.min(shift_x, block_size / 2) 
                    << int64_t(length(seqs[sid]) - 1 - get_cord_x(cords[i - 1]));
            g_cmpll.min(shift_y, block_size / 2)
                    << int64_t(length(read) - 1 - get_cord_y(cords[i - 1]));
            uint64_t cord2 = shift_cord(cords[i], shift_x, shift_y);
            g_cmpll.min(shift_x, shift_x) << get_cord_x(cords[i] - cords[i - 1]); //ins
            g_cmpll.min(shift_y, shift_y) << get_cord_y(cords[i] - cords[i - 1]); //del
            //uint64_t cord1 = shift_cord(cords[i - 1], shift_x, shift_y);
            uint64_t cord1 = cords[i - 1]; 
            cord2 = cords[i];
            mapGap_(seqs, read, comstr, 
                    cord1, cord2, 
                    g_hs, g_anchor, f1, f2, 
                    tiles, clips, 
                    g_align_closed,
                    thd_cord_gap, 
                    thd_tileSize,
                    thd_cord_remap);
            if (!empty(tiles))
            {
                uint64_t recd = get_cord_recd(cords[i]);//set cord flag
                set_tiles_flags_ (tiles, recd);
                insert(cords, i, tiles);
                i += length(tiles);
                clear(tiles);
            }   
            gap_len += std::max(get_cord_x(cord2 - cord1), 
                                get_cord_y(cord2 - cord1));
        }
        if (_DefaultHit.isBlockEnd(cords[i]))  ///right clip end cord
        {
            uint64_t cord1 = cords[i];
            if (length(read) - 1 - get_cord_y(cord1) > thd_cord_gap)
            {
                shift_x = std::min(thd_max_extend, 
                                  length(seqs[sid]) - get_cord_x(cord1));
                shift_y = std::min(thd_max_extend, 
                                   length(read) - get_cord_y(cord1) - 1);
                uint64_t infi_cord = _DefaultCord.shift(cord1, shift_x, shift_y);
                uint64_t cord2 = infi_cord;
                mapGap_ (seqs, read, comstr, 
                         cord1, cord2, 
                         g_hs, g_anchor, f1, f2,  
                         tiles, clips, 
                         g_align_right,
                         thd_cord_gap, 
                         thd_tileSize, 
                         thd_cord_remap);
                if (!empty(tiles))
                {
                    uint64_t recd = get_cord_recd(cords[i]);//set cord flag
                    set_tiles_flags_ (tiles, recd);
                    insert(cords, i + 1, tiles);
                    _DefaultHit.unsetBlockEnd(cords[i]);
                    i += length(tiles);
                    _DefaultHit.setBlockEnd(cords[i]);
                    clear(tiles);
                }   
                gap_len += length(read) - get_cord_y(cord1) - window_size;
            }
        }
        
    }
    return gap_len;
}
/*=====  End of Index free Map and clip  ======*/

