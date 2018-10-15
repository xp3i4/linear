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

/*
 *NOTE! the following parameters highly related.
 * Do not change them independently 
 */
int const g_shape_len = 8;
int const g_thd_anchor = 6;
float const g_thd_anchor_density = 0.03;
float const g_thd_error_percent = 0.2;

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

//Format::Tile N/A[2]|strand[1]|tileEnd[1]|x[40]|y[20]
//same as the Format::Cord
struct TileBase
{
    uint64_t const xBitLen = 40;
    uint64_t const yBitLen = 20;
    uint64_t const strandBit = 61;
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
    inline uint16_t getStrand (uint64_t val,
                            uint64_t bit = _defaultTileBase.strandBit
                              )
    {
        return (val >> bit) & 1;
    }
    inline uint64_t slide (uint64_t val, 
                           uint64_t x, 
                           uint64_t y,
                           uint64_t bit = _defaultTileBase.yBitLen)
    {
        return val + (x << bit) + y;
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
    //std::cout << "[]::mapGap_ " << start2 << " " << end2 << "\n";
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
static const uint64_t g_hs_anchor_mask1_ = ~ g_hs_anchor_mask1;
static const uint64_t g_hs_anchor_mask3 = (1ULL << 30) - 1;
static const uint64_t g_hs_anchor_mask5 = (1ULL << 31) - 1;
static const uint64_t g_hs_anchor_bit1 = 20;
static const uint64_t g_hs_anchor_bit2 = 50;
static const uint64_t g_hs_anchor_mask2 = ~(1ULL << 50);

static const int g_sv_inv = 1;    //inversion
static const int g_sv_ins = 2;    //insertion
static const int g_sv_del = 4;    //deletion
static const int g_sv_trans = 8;  //translocation
static const int g_sv_dup = 16;   //duplication


inline uint64_t g_hs_anchor_getCord (uint64_t & anchor)
{
    return anchor & g_hs_anchor_mask1;
}

inline uint64_t g_hs_anchor_getAnchor (uint64_t anchor)
{
    return (anchor >> g_hs_anchor_bit1) & g_hs_anchor_mask5;
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
    //std::cout << "[]::g_hs_setAnchor_ " << hs2 << " " << (g_hs_anchor_getY(val)) << " " << (hs2 & g_hs_mask2) << " "<< revscomp_const << " " << strand << "\n";
}

inline uint64_t g_hs_getXT (uint64_t const & val)
{
    return (val >> 31) & g_hs_mask3;
}

inline uint64_t g_hs_anchor_2Tile (uint64_t & anchor, uint64_t & main_strand, uint64_t revscomp_const)
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
    //std::cerr << "[]::y " << strand <<  " " << g_hs_anchor_getY(anchor) << "\n";
	return (((anchor + ((anchor & g_hs_anchor_mask1) << 20)) & g_hs_anchor_mask2) & g_hs_anchor_mask1_) + y + (strand << 61);
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

/**
 * Keep k-mers from the 'seq' into 'g_hs'
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

/**
 * Stream part of 'g_hs' and convert the elements into anchors
 */
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
/*
 * select the longest anchor: faster but lose accuracy
inline void g_mapHs_anchor_ (String<uint64_t> & anchor, 
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
inline void g_mapHs_anchor_ (String<uint64_t> & anchor, 
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
        if (g_hs_anchor_getAnchor(anchor[k] - anchor[prek]) > thd_error_percent * std::max(thd_min_segment, d))
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
                        appendValue (tile, g_hs_anchor_2Tile(anchor[j - 1], main_strand, revscomp_const));
                    }
                }
                appendValue (tile, g_hs_anchor_2Tile(anchor[k - 1], main_strand, revscomp_const));
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
 */
inline void g_mapHs_anchor_sv_ (String<uint64_t> & anchor, 
                                String<uint64_t> & tile, 
                                StringSet<String<short> > & f1,
                                StringSet<String<short> >& f2,
                                uint64_t gs_start,
                                uint64_t gs_end,
                                uint64_t gr_start,
                                uint64_t gr_end,
                                int anchor_end, 
                                int thd_tileSize,
                                uint64_t  main_strand, 
                                int revscomp_const
                                )
{
    int64_t thd_min_segment = 100;
    int thd_k_in_window = 1;
    float thd_overlap_tile = thd_tileSize * 0.4;
    float thd_error_percent = 0.6;
    
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
        if (g_hs_anchor_getAnchor(anchor[k] - anchor[prek]) > thd_error_percent * std::max(thd_min_segment, d))
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
                int kcount = 0;
                uint64_t first_low_bound_x = g_hs_anchor_getX(anchor[prek]) + window_size;
                uint64_t first_low_bound_y = g_hs_anchor_getY(anchor[prek]) + window_size;
                for (int i = prek + 1; g_hs_anchor_getX(anchor[i]) < first_low_bound_x && g_hs_anchor_getY(anchor[i]) < first_low_bound_y; i++)
                {
                    kcount++;
                }
                for (int j = prek + 1; j < k; j++)
                {
                    if ((g_hs_anchor_getX(anchor[j]) > prex + thd_tileSize
                        || g_hs_anchor_getY(anchor[j]) > prey + thd_tileSize))
                    {
                        prex = g_hs_anchor_getX(anchor[j - 1]);
                        prey = g_hs_anchor_getY(anchor[j - 1]);
                        appendValue (tile, g_hs_anchor_2Tile(anchor[j - 1], main_strand, revscomp_const));
                        appendValue (score, kcount);
                        kcount=0;
                    }
                    else
                    {
                        kcount++;
                    }
                }
                appendValue (tile, g_hs_anchor_2Tile(anchor[k - 1], main_strand, revscomp_const));
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
    std::cout << "g_mapHs_anchor_sv_ thd_tileSize " << thd_tileSize << "\n";
    float tz = thd_tileSize / 2;
    int prep = 0;
    for (int i = 0; i < length(score); i++)
    {
        //std::cout << "score " << i << " " << score[i] << "\n";
        
        if (score[i] < thd_k_in_window) 
        {
            continue;
        }
        else
        {
            uint64_t tile_x = _defaultTile.getX(tile[i]);
            uint64_t tile_y = _defaultTile.getY(tile[i]);
            unsigned fscore = _windowDist(begin(f1[_defaultTile.getStrand(tile[i])]) + _DefaultCord.cord2Cell(tile_y), 
                                        begin(f2[_getSA_i1(tile_x)]) + _DefaultCord.cord2Cell(_getSA_i2(tile_x)));
            std::cout << "[]::g_mapHs_anchor_sv_ fscore " << i << " " << fscore << " " << tile_y << " "<< length(score)<< "\n";
            if (fscore < windowThreshold)
            {
                tile[prep] = tile[i];
                score[prep] = score[i];
                prep++;  
            }
        }
    }
    
    std::sort (begin(tile), 
               begin(tile) + prep,
               [](uint64_t & s1, uint64_t & s2)
               {
                   return _defaultTile.getX(s1) < _defaultTile.getX(s2);
            });
    /**
     * remove tiles overlap with its adjacent tiles except the first and last tiles
     */
    int prep2 = 0;
    for (int i = 1; i < prep - 1; i++)
    {
        
        /// the if condition can't be or '||', since for dels and invs edges of one side can be very close.
        /// ---------------------- x
        ///         /t1/\ t2\ 
        ///        -----------     y
        /// Above tile t1 and t2, y of two tiles is far away while x of two tiles are almost overlapped 
        
        if (_defaultTile.getX(tile[i + 1] - tile[prep2]) < thd_overlap_tile && 
            _defaultTile.getY(tile[i + 1] - tile[prep2]) < thd_overlap_tile)
        {
            continue;
        }
        else
        {
            tile[prep2] = tile[i];
            prep2++;
        }
    }
    resize (tile, prep2);
    
    /**
     * extend window if there are gaps between tiles
     */
    //extend the middle tiles
    uint64_t genomeId = _getSA_i1(_defaultTile.getX(tile[0]));
    for (int i = 1; i < length(tile); i++)
    {

        std::cout << "[]::g_mapHs_anchor_sv_::patch " << _defaultTile.getY(tile[i - 1]) << " " << _defaultTile.getY(tile[i])<< "\n";
        i += extendPatch(f1[_defaultTile.getStrand(tile[i])], f2[genomeId], tile, i, tile[i - 1], tile[i]);   
    }
    //extend the last and first tiles
    uint64_t startCord = _DefaultCord.createCord(_createSANode(genomeId, gs_end), gr_end, _defaultTile.getStrand(tile[0]));
    extendPatch(f1[_defaultTile.getStrand(tile[0])], f2[genomeId], tile, 0, startCord, tile[0]);
    uint64_t endCord = _DefaultCord.createCord(_createSANode(genomeId, gs_end), gr_end, _defaultTile.getStrand(back(tile)));
    extendPatch(f1[_defaultTile.getStrand(back(tile))], f2[genomeId], tile, length(tile), back(tile), endCord);   
}
/**
* NOTE:!!! this is the part test for pairwise alignment in the tile
*/
inline int g_align_gap_(String<Dna5> & seq,
                        String<Dna5> & read,
                        String<Dna5> & comstr, /*complement reverse of the read*/
                        String<uint64_t> & tiles,
                        uint64_t g_start,
                        uint64_t g_end, 
                        uint64_t r_start,
                        uint64_t r_end,
                        uint64_t  main_strand
                        )
{
    int thd_overlap = 50;
    String<int> score;
    String<int> score_mark;
    String<int> patch_cord;
    resize(score, length(tiles));
    resize(score_mark, length(tiles));
    /**
     * check multiple path
     * \
     *  \ ---
     *   \   \
     *    \   \
     *         \
     *   path1 path2
     */
   /* 
    for (int k = 0; k < length(tiles); k++)
    {
        if ()
    }
    */
    /**
     * scan for sv and clip the break point
     * discontinuous tiles: possible ins or del
     * strand flip: invs
     */ 
    int clip = -1;
    int sv_flag = 0;
    for (int k = 1; k < length(tiles); k++)
    {
        /**
         * decide sv type
         */
        if (_defaultTile.getStrand(tiles[k - 1] ^ tiles[k]))
        {
            sv_flag |= g_sv_inv;
        }
        if (_getSA_i2(_defaultTile.getX(tiles[k] - tiles[k - 1])) > window_size)  //window: 192x192
        {
            sv_flag |= g_sv_ins;
        }
        if (_defaultTile.getY(tiles[k] - tiles[k - 1]) > window_size)
        {
            sv_flag |= g_sv_del;
        }
        
        /**
         * clip sv
         */
        
        clip = -1;
        switch (sv_flag)
        {
            case g_sv_inv: //1
                break;
            case g_sv_ins: //2
                break;
            case g_sv_inv + g_sv_del : //5
                break;
                /*
            {
                uint64_t mid_cord = _defaultTile.slide(tiles[k - 1], window_size / 2, window_size / 2);
                std::cout << "xxxxxxxif0\n";
                clip = clip_window(seq, read, comstr, mid_cord);
                if (clip < 0)
                {
                    std::cout << "xxxxxxxif1\n";
                    clip = clip_window(seq, read, comstr, tiles[k-1]);
                    if (clip < 0)
                    {
                    std::cout << "xxxxxxxif2\n";
                        clip = clip_window(seq, read, comstr, tiles[k]);
                    }
                }
            }*/
         //   default:
         //       std::cout << "[]::g_align_gap_: none \n";
        }
        std::cout << "[]::g_align_gap_ : sv_flag " << sv_flag << "\n";
        sv_flag = 0;
        
    }

    //std::cout << "[]::g_align_gap_ " << length
}

inline int g_mapHs_(String<Dna5> & seq, 
                    String<Dna5> & read,
                    String<Dna5> & comstr,
                    uint64_t gs_start, 
                    uint64_t gs_end, 
                    uint64_t gr_start,
                    uint64_t gr_end,
                    uint64_t main_strand,
                    String<uint64_t> & g_hs,
                    String<uint64_t> & g_hs_anchor,
                    String<uint64_t> & g_hs_tile,
                    StringSet<String<short> > & f1,
                    StringSet<String<short> >& f2,
                    int thd_tileSize)
{
    int g_hs_end = 0;
    int g_hs_anchor_end = 0;

    g_hs_end = g_mapHs_kmer_(seq, g_hs, gs_start, gs_end, g_hs_end, 10, 0);
    g_hs_end = g_mapHs_kmer_(read, g_hs, gr_start, gr_end, g_hs_end, 1, 1);

    //std::cout << "[g_mapHs_]::gs,r start " << gs_start << " " << gs_end << " " << gr_start << " " << gr_end << "\n";
    std::sort (begin(g_hs), begin(g_hs) + g_hs_end);
    int p1 = 0, p2 = 0;
    for (int k = 0; k < g_hs_end; k++)
    {
        //std::cout << "[g_mapHs_] " << k << " " << ((g_hs[k] >> 31) & 3) << " " << (g_hs[k] >> 33) << " " << (g_hs[k] & ((1ULL << 30) - 1)) << "\n";
        switch (g_hs_getXT(g_hs[k] ^ g_hs[k - 1]))
        {
            case 0:
                //std::cout << "case 0\n";
                break;
            case 1:
                //std::cout <<"case 1 \n";
                p2 = k;
                break;
            default:
                //std::cout << "case 2 " << p1 << " " << p2 << " " << k << "\n";
                g_hs_anchor_end = g_mapHs_setAnchors_(g_hs, g_hs_anchor, p1, p2, k, length(read) - 1, g_hs_anchor_end);
                p1 = k;
                p2 = k; 
        }
    }
    /*
    for (int k = 0; k < g_hs_anchor_end; k++)
    {
        std::cout << "[3]::g_mapHs_::anchor " << g_hs_anchor_getAnchor(g_hs_anchor[k]) << " " << g_hs_anchor_getY(g_hs_anchor[k]) << "\n";
    }
    */
    //g_mapHs_anchor_(g_hs_anchor, g_hs_tile, g_hs_anchor_end, thd_tileSize, main_strand, length(read) - 1);
    g_mapHs_anchor_sv_(g_hs_anchor, 
                       g_hs_tile, 
                       f1, 
                       f2, 
                       gs_start, 
                       gs_end, 
                       gr_start, 
                       gr_end, 
                       g_hs_anchor_end, 
                       thd_tileSize, 
                       main_strand, 
                       length(read) - 1);
    g_align_gap_(seq, read, comstr, g_hs_tile, gs_start, gs_end, gr_start, gr_end, main_strand);
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

/*
 * Extract gaps from cords then re-map gaps.
 * use table to map
 * []::mg3
 * gaps at the begin or end of the read are also processed.
 */
int mapGaps(StringSet<String<Dna5> > & seqs, 
            String<Dna5> & read, 
            String<Dna5> & comstr,
            String<uint64_t> & cords, 
            String<uint64_t> & g_hs,
            String<uint64_t> & g_anchor,
            StringSet<String<short> > & f1,
            StringSet<String<short> >& f2,
            int const thd_gap, 
            int const thd_tileSize)
{
    String <uint64_t> tile;
    Gap gap;
    uint64_t delta = (uint64_t)thd_tileSize / 2;
    unsigned count = 0;
    /*
    if (_DefaultCord.getX(cord[0]) > thd_gap )
    {
        uint64_t strand = _DefaultCord.getCordStrand (cords[0]);
        uint64_t gr_start = (length(read) - 1) * strand
        uint64_t gr_end = (length(read) - 1) * strand - (_nStrand(strand) * (_DefaultCord.getCordY(cords[k])));
        if (strand)
            std::swap (gr_start, gr_end);
                    g_mapHs_(seqs[_getSA_i1(_DefaultCord.getCordX(cords[0]))], 
                     read, 
                     _getSA_i2(_DefaultCord.getCordX(cords[k - 1])),
                     _getSA_i2(_DefaultCord.getCordX(cords[0]) + delta),
                     gr_start,
                     gr_end,
                     strand,
                     //_DefaultCord.getCordY(cords[k -1]) + delta,
                     //_DefaultCord.getCordY(cords[k]) + delta,
                     g_hs,
                     g_anchor,
                     tile,
                     thd_tileSize
                    );
    }
    */
    /*
     * TODO the first and last cord
     */
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
            		(_nStrand(strand) * (_DefaultCord.getCordY(cords[k - 1])));
            uint64_t gr_end = (length(read) - 1) * strand - 
            		(_nStrand(strand) * (_DefaultCord.getCordY(cords[k])));
            if (strand)
                std::swap(gr_start, gr_end);
            std::cout << "[]::mapGaps " << gr_start << " " << gr_end << "\n";
            g_mapHs_(seqs[_getSA_i1(_DefaultCord.getCordX(cords[k - 1]))], 
                     read, 
                     comstr,
                     _getSA_i2(_DefaultCord.getCordX(cords[k - 1])),
                     _getSA_i2(_DefaultCord.getCordX(cords[k])),
                     gr_start,
                     gr_end,
                     strand,
                     //_DefaultCord.getCordY(cords[k -1]) + delta,
                     //_DefaultCord.getCordY(cords[k]) + delta,
                     g_hs,
                     g_anchor,
                     tile,
                     f1,
                     f2,
                     thd_tileSize
                    );
            //g_align_gap_(seqs, read, comstr, tile);
            //align(seqs, read, comstr, tile);
            std::cout << "[]::mapGaps:check_tiles_ " << check_tiles_(tile, _getSA_i2(_DefaultCord.getCordX(cords[k - 1])), _getSA_i2(_DefaultCord.getCordX(cords[k]))) << "\n";
            count += _DefaultCord.getCordX(cords[k] - cords[k - 1]);
            insert(cords, k, tile);
            k += length(tile);
        }
    }
    //TODO: tile -> cords
    return count;
}


/*
 * Extract gaps from cords then re-map gaps.
 * use table to map
 * []::mg3
 * gaps at the begin or end of the read are note processed.
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
    uint64_t delta = (uint64_t)thd_tileSize / 2;
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
            		(_nStrand(strand) * (_DefaultCord.getCordY(cords[k - 1])));
            uint64_t gr_end = (length(read) - 1) * strand - 
            		(_nStrand(strand) * (_DefaultCord.getCordY(cords[k])));
            if (strand)
                std::swap(gr_start, gr_end);
            g_mapHs_(seqs[_getSA_i1(_DefaultCord.getCordX(cords[k - 1]))], 
                     read, 
                     _getSA_i2(_DefaultCord.getCordX(cords[k - 1])),
                     _getSA_i2(_DefaultCord.getCordX(cords[k])),
                     gr_start,
                     gr_end,
                     strand,
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
*/

int mapGaps(StringSet<String<Dna5> > & seqs, StringSet<String<Dna5> > & reads, 
            StringSet<String<uint64_t> > & cords, unsigned const thd_gap, unsigned const thd_tileSize)
{
    for (unsigned k = 0; k < length(reads); k++)
    {
        mapGaps(seqs, reads[k], cords[k], thd_gap, thd_tileSize);
    }
    return 0;
}










