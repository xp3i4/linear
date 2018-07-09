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

int const g_shape_len = 18;

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

//gap start[40]|gap length[24]
struct gapbase_
{
    const unsigned cBitLen = 40;
    const unsigned bBitLen = 24;
    const unsigned cmask = (1ULL << cBitLen) - 1;
    const unsigned bmask = (1ULL << bBitLen) - 1;
} _defaultgapbase_;

struct Gap_
{
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
//N/A[2]|strand[1]|N/A[1]|anchor[40]|coord[20]
//strand = 1 or -1 
struct ACoordBase
{
    unsigned const sBitLen = 1;
    unsigned const aBitLen = 40;
    unsigned const cBitlen = 20;
    unsigned const sBit = 61;
    
    uint64_t const amask = (1ULL << aBitLen) - 1;
    
}_defaultACoordBase;

struct ACoord
{
    makeValue(uint64_t strand, uint64_t anchor, uint64_t coord) 
    {
        return (((anchor * strand) & _defaultACoordBase.amask) << abit) + coord;
    }
    uint64_t reverseAnchor(uint64_t & anchor, uint64_t const & mask = _defaultACoordBase.amask)
    {
        return (-anchor) & mask
    }
    
}_defaultACoord;

/*
 * flip strand from 0, 1 to -1, 1;
 * strand = 0, 1, other values is not allowed
 * return -1 , 1
 */
inline uint64_t _nStrand(uint64_t strand)
{
    return (strand << 1) - 1;
}

int mapGap_(GIndex & g_index, String <Dna5> & read,  uint64_t start2, uint64_t end2)
{
    float thd_error_percent = 0.2; 
    unsigned thd_min_segment = 100;
    double time = sysTime();
    unsigned count = 0, step = 5;
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
                while (_defaultGNode.getXValue(g_index.g_hs[hsStart]) == g_index.shape.XValue)
                {
                    uint64_t strand = _defaultGNode.getStrand(g_index.g_hs[hsStart]) ^ g_index.shape.strand;
                    //appendValue(anchor, ((_defaultGNode.getCoord(g_index.g_hs[hsStart]) + _nStrand(strand) * k) << 20) | k | (strand << 63));
                    appendValue(anchor,  _defaultACoord.makeValue(_defaultGNode.getCoord(g_index.g_hs[hsStart]) +
                                _nStrand(strand) * k, k, strand));
                    ++hsStart;
                    // if strand == 1, on different strands
                    // + k
                    // if strand == 0,   on the same strand
                    // - k
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
    for (unsigned k = 0; k < length(anchor); k++)
    {
        //TODO: handle thd_min_segment, anchor 
        //if (anchor[k] - anchor[prek] > thd_segment * anchor[] )
        if (anchor[k] - anchor[prek] > (thd_error_percent * std::max(thd_min_segment, anchor) << 20) )
        {
            std::sort (begin(anchor) + prek, begin(anchor) + k, 
                       [](uint64_t & s1, uint64_t & s2){
              //  return _defaultAnchorCoord.getCoord(s2 - s1) > 0;
                return (s2 & ((1ULL << 20) - 1)) > (s1 & ((1ULL << 20) - 1)) ;
            });
            /*
            for (unsigned j = prek; j < k; j++)
            {
                if ()
            }
            */
            prek = k;
        //TODO:transform seed to tiles
        }
    }
    for (unsigned k = 0; k < length (anchor); k++)
    {
        std::cout << "[]::mapGap_ " << (anchor[k] >> 20) << " " << (anchor[k] & ((1 << 20) - 1)) << "\n";
    }
}

int mapGap(String <Dna5> & seq, String <Dna5> & read, Gap & gap)
{
    GIndex g_index;
    g_createDir(seq, gap.getStart1(), gap.getEnd1(), g_index);
    mapGap_ (g_index, read, gap.getStart2(), gap.getEnd2());
    //mapGap_ (g_index, _reverse(read), length(read) - end2 - 1, length(read) - start2 - 1);
    //optimizeGap();
    return 0;
}

int mapGaps(String<Dna5> & seq, String<Dna5> & read, String <Gap> & gaps)
{
    //collectGaps(cords, gaps);
    //mergeGaps(gaps);
    for (unsigned k = 0; k < length(gaps); k++)
    {
        mapGap(seq, read, gaps[k]);
    }
    return 0;
}

int g_test1 (String<Dna5> & seq)
{
    double time = sysTime ();
    GIndex g_index;
    for (unsigned k = 0; k < 100000; k++)
    {
        g_createDir(seq, k, k+1000, g_index);
        hashInit(g_index.shape, begin(seq) + k);
        unsigned count = 0;
        for (unsigned j = k; j < k + 1000; j++)
        {
            if (++count == 10)
            {
                hashNext(g_index.shape, begin(seq) + j);
                if (_defaultGNode.getXValue(g_index.g_hs[g_index.g_dir[g_index.shape.XValue]]) != g_index.shape.XValue)
                {
                    std::cerr << "[]::g_test::error " << k << " " << j << "\n";
                    return 1;
                }
                count = 0;
            }
            else
            {
                hashNexth(g_index.shape, begin(seq) + j);
            }
        }
    }
    std::cerr << "[]::g_test " << sysTime () - time << "\n";
    return 0;
}

int g_test (String<Dna5> & seq)
{
    std::cerr << "[]::g_test \n";
    double time = sysTime ();
    Gap gap (100, 700, 200, 500);
    for (unsigned k = 0; k <1; k++)
    {
        gap.setValue(k, k + 1000, k, k + 1000);
        mapGap(seq, seq, gap);
    }
    std::cerr << "[]::g_test " << sysTime () - time << "\n";
    return 0;
}






