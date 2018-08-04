//GNode: N/A[1]|xval[32]|strand[1]|coordinate[30]
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
        {
        }
}_defaultGNodeBase;

struct GNode
{
    inline void setValue (uint64_t & val, uint64_t const & xval, 
                   uint64_t const & strand, uint64_t const & coordinate, 
                   uint64_t const & sbit = _defaultGNodeBase.sBitLen,
                   uint64_t const & cbit = _defaultGNodeBase.cBitLen)
    {
        val = (xval << (cbit + sbit)) + (strand << cbit) + coordinate;
    }
    inline uint64_t makeValue (uint64_t const & xval, uint64_t & strand, uint64_t const & coordinate, 
                        uint64_t const & sbit = _defaultGNodeBase.sBitLen,
                        uint64_t const & cbit = _defaultGNodeBase.cBitLen)
    {
        return (xval << (cbit + sbit)) + (strand << cbit) + coordinate;
    }
    inline uint64_t getXValue (uint64_t const & xval, 
                        uint64_t const & bit = _defaultGNodeBase.cBitLen + _defaultGNodeBase.sBitLen, 
                        uint64_t const & xmask = _defaultGNodeBase.xmask)
    {
        return (xval >> bit) & xmask;
    }
    inline uint64_t getStrand (uint64_t const & xval, 
                        uint64_t const & bit = _defaultGNodeBase.cBitLen, 
                        uint64_t const & mask = _defaultGNodeBase.smask)
    {
        return (xval >> bit) & mask;
    }
    inline uint64_t getCoord (uint64_t const & xval, 
                        uint64_t const & mask = _defaultGNodeBase.cmask)
    {
        return (xval & mask);
    }
}_defaultGNode;

int const g_shape_len = 15;

struct GIndex
{
    String <uint64_t> g_hs;
    String <uint64_t> g_dir;
    Shape<Dna5, Minimizer<g_shape_len> > shape;
};


int g_createDir_(String<Dna5> & seq, uint64_t gs_start, uint64_t gs_end, 
                String<uint64_t> & g_hs, String<uint64_t> & g_dir, Shape<Dna5, Minimizer<g_shape_len> > & shape)
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
    inline uint64_t getGStart(uint64_t & val, 
                       uint64_t const & bit = _defaultgapbase_.rBitLen + _defaultgapbase_.bBitLen, 
                       uint64_t const & mask = _defaultgapbase_.gmask, 
                       uint64_t const & bbit = _defaultgapbase_.bBit)
    {
        return ((val >> bit) & mask) << bbit;
    }
    inline uint64_t getRStart(uint64_t & val, 
                       uint64_t const & bit = _defaultgapbase_.bBitLen, 
                       uint64_t const & mask = _defaultgapbase_.rmask,
                       uint64_t const & bbit = _defaultgapbase_.bBit)

    {
        return ((val >> bit) & mask) << bbit;
    }
    inline uint64_t getGEnd(uint64_t & val, 
                    uint64_t const & bit = _defaultgapbase_.gBitLen, 
                    uint64_t const & smask = _defaultgapbase_.smask,
                    uint64_t const & gmask = _defaultgapbase_.gmask)
    {
        return ((val >> bit) & smask) + (val & gmask);
    }
    inline uint64_t makeValue (uint64_t start, uint64_t len, 
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
    inline uint64_t getId (uint64_t & val, 
                           uint64_t const & bit = _defaultgapbase_.cBitLen + _defaultgapbase_.bBitLen,
                           uint64_t const & mask = _defaultgapbase_.bmask)
    {
        return (val >> bit) & mask;
    }
    inline uint64_t getStart(uint64_t & val, 
                      uint64_t const & bit = _defaultgapbase_.bBitLen, 
                      uint64_t const & mask = _defaultgapbase_.cmask)
    {
        return (val >> bit) & mask;
    }
    inline uint64_t getEnd(uint64_t & val, 
                    uint64_t const & bit = _defaultgapbase_.bBitLen, 
                    uint64_t const & cmask = _defaultgapbase_.cmask,
                    uint64_t const & bmask = _defaultgapbase_.bmask)
    {
        return ((val >> bit) & cmask) + (val & bmask);
    }
    inline uint64_t getLength (uint64_t & val, uint64_t const & mask = _defaultgapbase_.bmask)
    {
        return val & mask;
    }
    inline uint64_t makeValue (uint64_t start, uint64_t len, 
                        uint64_t const & bit = _defaultgapbase_.bBitLen)
    {
        return (start << bit) + len;
    }
}_defaultGap_;

struct Gap
{
    uint64_t gap1; //seq
    uint64_t gap2; //read
    inline uint64_t getId1()
    {
        return _defaultGap_.getId(gap1);
    }
    inline uint64_t getId2()
    {
        return _defaultGap_.getId(gap2);
    }
    inline uint64_t getStart1 ()
    {
        return _defaultGap_.getStart (gap1);
    }
    inline uint64_t getStart2 ()
    {
        return _defaultGap_.getStart (gap2);
    }
    inline uint64_t getEnd1 ()
    {
        return _defaultGap_.getEnd (gap1);
    }
    inline uint64_t getEnd2 ()
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
    inline void setValue (uint64_t start1, uint64_t end1, uint64_t start2, uint64_t end2)
    {    
        gap1 = _defaultGap_.makeValue (start1, end1 - start1);
        gap2 = _defaultGap_.makeValue (start2, end2 - start2);   
    }
};

inline Gap makeGap (uint64_t start1, uint64_t end1, uint64_t start2, uint64_t end2)
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

/*
 * flip strand from 0, 1 to -1, 1;
 * strand = 0, 1, other values is not allowed
 * return -1 , 1
 */
inline uint64_t _nStrand(uint64_t strand)
{
    return (strand << 1) - 1;
}

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
    inline uint64_t makeValue(uint64_t cf, uint64_t cr,  uint64_t strand,
                       uint64_t const & sbit = _defaultACoordBase.sBit,
                       uint64_t const & bit = _defaultACoordBase.cBitLen) 
    {
        return (strand << sbit) + ((cf + _nStrand(strand) * cr) << bit) + cr;
    }
    inline uint64_t reverseAnchor(uint64_t & anchor, uint64_t const & mask = _defaultACoordBase.amask)
    {
        return (-anchor) & mask;
    }
    inline uint64_t getAnchor (uint64_t val, 
                        uint64_t const & mask = _defaultACoordBase.amask, 
                        uint64_t const & bit = _defaultACoordBase.cBitLen,
                        uint64_t const & mask2 = (1ULL << _defaultACoordBase.sBit)) 
    {
        return (val & mask2) + ((val >> bit) & mask);
    }
    inline uint64_t getCoord (uint64_t val, 
                       uint64_t const & mask = _defaultACoordBase.cmask)
    {
        return val & mask;
    }
    inline uint64_t getX (uint64_t val, 
                   uint64_t const & bit = _defaultACoordBase.cBitLen, 
                   uint64_t const & mask = _defaultACoordBase.amask, 
                   uint64_t const & bit2 = _defaultACoordBase.sBit,
                   uint64_t const & mask3 = _defaultACoordBase.cmask
                  )
    {
        return ((val >> bit) & mask) + _nStrand((val >> bit2) & 1) * (val & mask3);
    }
    inline uint64_t getY (uint64_t val)
    {
        return getCoord(val);
    }
    
}_defaultACoord;

//Tile N/A[2]|strand[1]|tileEnd[1]|x[40]|y[20]
//same format as cord
struct TileBase
{
    uint64_t const xBitLen = 40;
    uint64_t const yBitLen = 20;
    uint64_t const xmask = (1ULL << xBitLen) - 1;
    uint64_t const ymask = (1ULL << yBitLen) - 1;
}_defaultTileBase;

struct Tile
{
    inline uint64_t getX (uint64_t val, 
                   uint64_t const & bit = _defaultTileBase.yBitLen,  
                   uint64_t const & mask = _defaultTileBase.xmask)
    {
        return (val >> bit) & mask;
    }
    inline uint64_t getY (uint64_t val, 
                   uint64_t const & mask = _defaultTileBase.ymask)
    {
        return val & mask;
    }
    inline uint64_t makeValue(uint64_t x, uint64_t y, uint64_t const & bit = _defaultTileBase.yBitLen)
    {
        return (x << bit) + y;
    }
    
}_defaultTile;

inline uint64_t acoord2Tile(uint64_t val, 
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
    std::cout << "[]::mapGap_ " << start2 << " " << end2 << "\n";
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
    std::cout << "[]::len anchor" << length(anchor) << "\n";
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


//========================
//g_hs_anchor: N/A[13]|strand[1]|anchor[30]|cord[20]
static const uint64_t g_hs_anchor_mask1 = (1ULL << 20) - 1;
static const uint64_t g_hs_anchor_mask3 = (1ULL << 30) - 1;
static const uint64_t g_hs_anchor_bit1 = 20;
static const uint64_t g_hs_anchor_bit2 = 50;
static const uint64_t g_hs_anchor_mask2 = ~(1ULL << 50);

inline uint64_t g_hs_anchor_getCord (uint64_t & anchor)
{
    return anchor & g_hs_anchor_mask1;
}

inline uint64_t g_hs_anchor_getAnchor (uint64_t anchor)
{
    return (anchor >> g_hs_anchor_bit1) & g_hs_anchor_mask3;
}

inline uint64_t g_hs_anchor_getX (uint64_t & val)
{
    return ((val >> g_hs_anchor_bit1) & g_hs_anchor_mask3) + (val & g_hs_anchor_mask1);
}

inline uint64_t g_hs_anchor_getY (uint64_t & val)
{
    return val & g_hs_anchor_mask1;
}


//g_hs: N/A[1]|xval[30]|type[2]strand[1]|coordinate[30]
static const uint64_t g_hs_mask1 = (1ULL << 30);
static const uint64_t g_hs_mask2 = (1ULL << 30) - 1;
static const uint64_t g_hs_mask3 = (1ULL << 32) - 1;

inline void g_hs_setGhs_(uint64_t & val, 
                         uint64_t xval, 
                         uint64_t type, 
                         uint64_t strand, 
                         uint64_t coord)
{
    val = (xval << 33) + (type<< 31) + (strand << 30) + coord;
}

inline void g_hs_setAnchor_(uint64_t & val, 
                            uint64_t const & hs1, /*genome*/
                            uint64_t const & hs2, /*read*/
                            uint64_t revscomp_const)
{
    uint64_t strand = ((hs1 ^ hs2) >> 30 ) & 1;
    uint64_t x = revscomp_const * strand - _nStrand(strand) * (hs2 & g_hs_mask2); 
    //std::cout << "[]::g_hs_setAnchor_ " << revscomp_const << " " << strand << " " << x << "\n";
    val = (((hs1 - x) & (g_hs_mask2)) << 20) + x + (strand << g_hs_anchor_bit2);
}

inline uint64_t g_hs_getXT (uint64_t const & val)
{
    return (val >> 31) & g_hs_mask3;
}

inline uint64_t g_hs_anchor_2Tile (uint64_t & anchor)
{
    uint64_t strand = (anchor >> g_hs_anchor_bit2) & 1;
	return ((anchor + ((anchor & g_hs_anchor_mask1) << 20)) & g_hs_anchor_mask2) + (strand << 61);
}

/*
 * collecting minimizer
 *
inline int g_mapHs_kmer_(String<Dna5> & seq, 
                         String<uint64_t> & g_hs, 
                         uint64_t start, 
                         uint64_t end, 
                         int g_hs_start, 
                         int step,  
                         uint64_t type)
{
    Shape<Dna5, Minimizer<g_shape_len> >  shape;
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
/*
 * collecting kmers
 */
inline int g_mapHs_kmer_(String<Dna5> & seq, 
                         String<uint64_t> & g_hs, 
                         uint64_t start, 
                         uint64_t end, 
                         int g_hs_start, 
                         int step,  
                         uint64_t type)
{
    Shape<Dna5, Minimizer<g_shape_len> >  shape;
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
            g_hs_setGhs_(g_hs[g_hs_start + i++], val, type, shape.strand, k);
            count = 0;
        }
    }
    return g_hs_start + i;
}

inline int g_mapHs_setAnchors_ (String<uint64_t> & g_hs, 
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

inline void g_mapHs_anchor_ (String<uint64_t> & anchor, 
                             String<uint64_t> & tile, 
                             int anchor_end, 
                             int thd_tileSize
                            )
{
    int64_t thd_min_segment = 100;
    int64_t prek = 0;
    int anchor_len = 0, max_anchor_len = 0, max_prek = 0, max_k = 0;
    float thd_error_percent = 0.2;
    std::sort (begin(anchor), begin(anchor) + anchor_end);
    anchor[anchor_end] = ~0;
    for (int k = 0; k < anchor_end + 1; k++)
    {
        //TODO: handle thd_min_segment, anchor 
        //int64_t d = std::abs(_defaultACoord.getCoord(anchor[k]) - _defaultACoord.getCoord(anchor[prek]));
        int64_t d = std::abs((int64_t)g_hs_anchor_getY(anchor[k]) - (int64_t)g_hs_anchor_getY(anchor[prek]));
        //std::cout << "[]::g_mapHs_anchor_ 3 " << std::max(thd_min_segment, d) << "\n";
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

    std::sort (begin(anchor) + max_prek, 
            begin(anchor) + max_k, 
            [](uint64_t & s1, uint64_t & s2)
            {
                return g_hs_anchor_getCord(s2) > g_hs_anchor_getCord(s1);
            });
    appendValue(tile, g_hs_anchor_2Tile(anchor[max_prek]));
    for (int j = max_prek + 1; j < max_k; j++)
    {
        if (g_hs_anchor_getX(anchor[j]) > _defaultTile.getX(back(tile)) + thd_tileSize
            || g_hs_anchor_getY(anchor[j]) > _defaultTile.getY(back(tile)) + thd_tileSize)
        {
            appendValue (tile, g_hs_anchor_2Tile(anchor[j - 1]));
        }
    }
     
}

inline int g_mapHs_(String<Dna5> & seq, 
                    String<Dna5> & read,
                    uint64_t gs_start, 
                    uint64_t gs_end, 
                    uint64_t gr_start,
                    uint64_t gr_end,
                    String<uint64_t> & g_hs,
                    String<uint64_t> & g_hs_anchor,
                    String<uint64_t> & g_hs_tile,
                    int thd_tileSize)
{
    int g_hs_end = 0;
    int g_hs_anchor_end = 0;
    g_hs_end = g_mapHs_kmer_(seq, g_hs, gs_start, gs_end, g_hs_end, 10, 0);
    g_hs_end = g_mapHs_kmer_(read, g_hs, gr_start, gr_end, g_hs_end, 1, 1);
    std::sort (begin(g_hs), begin(g_hs) + g_hs_end);
    int p1 = 0, p2 = 0;

    for (int k = 0; k < g_hs_end; k++)
    {
        switch (g_hs_getXT(g_hs[k] ^ g_hs[k - 1]))
        {
            case 0:
                break;
            case 1:
                p2 = k;
                break;
            default:
                g_hs_anchor_end = g_mapHs_setAnchors_(g_hs, g_hs_anchor, p1, p2, k, length(read) - 1, g_hs_anchor_end);
                p1 = k;
                p2 = k; 
        }
    }
    g_mapHs_anchor_(g_hs_anchor, g_hs_tile, g_hs_anchor_end, thd_tileSize);
}

/*
 * []::mg2 
 */
int mapGap(String <Dna5> & seq, String <Dna5> & read, Gap & gap, 
           String<uint64_t> & tile, uint64_t const thd_tileSize)
{
    GIndex g_index;
    g_createDir(seq, gap.getStart1(), gap.getEnd1(), g_index);
    std::cerr << "[]::mg2 " << gap.getStart1() << " " << gap.getEnd1() << "\n";
    mapGap_ (g_index, read, gap.getStart2(), gap.getEnd2(), tile, thd_tileSize);
    //mapGap_ (g_index, _reverse(read), length(read) - end2 - 1, length(read) - start2 - 1);
    //optimizeGap();
    //TODO: map _reverse and do optimiz
    return 0;
}

/*
 * Extract gaps from cords then re-map gaps.
 * use index to map
 * []::mg3
 */
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
            _DefaultCord.getCordY(cords[k] - cords[k - 1]) > thd_gap &&
            !_DefaultHit.isBlockEnd(cords[k - 1]))
        {
            clear(tile);
            gap = makeGap(_DefaultCord.getCordX(cords[k - 1]) + delta,
                          _DefaultCord.getCordX(cords[k]) + delta,
                          _DefaultCord.getCordY(cords[k - 1]) + delta, 
                          _DefaultCord.getCordY(cords[k]) + delta);
            //mapGap(seqs[_getSA_i1(_DefaultCord.getCordX(cords[k - 1]))], read, gap, tile, thd_tileSize);
            count += _DefaultCord.getCordX(cords[k] - cords[k - 1]);
            insert(cords, k, tile);
            k += length(tile);
        }
    }
 //   for (unsigned j = 1; j < length(cords); j++)
 //       std::cout << "[]::mg3 " << length(tile) << " " << gap.getId2() << " " << gap.getStart2() << " " << gap.getEnd2() << " " << _DefaultCord.getCordY(cords[j])  << " " << _DefaultCord.getCordX(cords[j])<< "\n";
    //TODO: tile -> cords

    return count;
}

/*
 * Extract gaps from cords then re-map gaps.
 * use table to map
 * []::mg3
 * start and end of the gap in this function are supposed to have the same strand. Otherwise(like inverse) need more process.
 */
int mapGaps(StringSet<String<Dna5> > & seqs, 
            String<Dna5> & read, 
            String<uint64_t> & cords, 
            String<uint64_t> & g_hs,
            String<uint64_t> & g_anchor,
            int const thd_gap, 
            int const thd_tileSize)
{
    String <uint64_t> tile;
    Gap gap;
    uint64_t delta = (uint64_t)thd_gap / 2;
    unsigned count = 0;
    for (unsigned k = 2; k < length(cords); k++)
    {
        if (_DefaultCord.getCordX(cords[k] - cords[k - 1]) > thd_gap && 
            _DefaultCord.getCordY(cords[k] - cords[k - 1]) > thd_gap &&
            !_DefaultHit.isBlockEnd(cords[k - 1]))
        {
            clear(tile);
            //strand of cord[k - 1] are supposed to be eqaul to strand of cord[k]
            uint64_t strand = _DefaultCord.getCordStrand (cords[k - 1]);
            uint64_t gr_start = (length(read) - 1) * strand - 
            		(_nStrand(strand) * (_DefaultCord.getCordY(cords[k - 1]) + delta));
            uint64_t gr_end = (length(read) - 1) * strand - 
            		(_nStrand(strand) * (_DefaultCord.getCordY(cords[k])));
            if (strand)
                std::swap(gr_start, gr_end);
            g_mapHs_(seqs[_getSA_i1(_DefaultCord.getCordX(cords[k - 1]))], 
                     read, 
                     _getSA_i2(_DefaultCord.getCordX(cords[k - 1])) + delta,
                     _getSA_i2(_DefaultCord.getCordX(cords[k])) + delta,
                     gr_start,
                     gr_end,
                     //_DefaultCord.getCordY(cords[k -1]) + delta,
                     //_DefaultCord.getCordY(cords[k]) + delta,
                     g_hs,
                     g_anchor,
                     tile,
                     thd_tileSize
                    );
            count += _DefaultCord.getCordX(cords[k] - cords[k - 1]);
            insert(cords, k, tile);
            k += length(tile);
        }
    }
    //TODO: tile -> cords
    return count;
}
int mapGaps(StringSet<String<Dna5> > & seqs, StringSet<String<Dna5> > & reads, 
            StringSet<String<uint64_t> > & cords, unsigned const thd_gap, unsigned const thd_tileSize)
{
    for (unsigned k = 0; k < length(reads); k++)
    {
        mapGaps(seqs, reads[k], cords[k], thd_gap, thd_tileSize);
    }
    return 0;
}

int mapGaps(StringSet<String<Dna5> > & seqs, 
            StringSet<String<Dna5> >& reads, 
            StringSet<String<uint64_t> > & cords, 
            int start1,
            int start2,
            int len,
            int const thd_gap, 
            int const thd_tileSize)
{
    
    String<uint64_t>  g_hs;
    String<uint64_t>  g_anchor;
    resize (g_hs, 1ULL << 20);
    resize (g_anchor, 1ULL<<20);
    for (unsigned k = 0; k < len; k++)
    {
        mapGaps(seqs, reads[start1 + k], cords[start2 + k], g_hs, g_anchor, thd_gap, thd_tileSize);
    }
}











