#ifndef SEQAN_HEADER_INDEX_UTIL_H
#define SEQAN_HEADER_INDEX_UTIL_H
#include <seqan/sequence.h>
#include <seqan/index.h>
#include "shape_extend.h"

using namespace seqan;
typedef uint64_t HType1;
typedef uint64_t HType2;

extern const uint64_t _dirEmpty;
extern const uint64_t _bitEmpty;
extern const uint64_t _bitLength;
extern const uint64_t _bitValue;
extern const unsigned _bitLength_END;
extern const unsigned _bitValue_END;
extern const unsigned _HeadValue_bits;
extern const uint64_t _HeadType_code;
extern const uint64_t _HeadTypeVtl_code;
extern const uint64_t _HeadTypeHVl_code;
extern const unsigned _BodyValue_bits;
extern const unsigned _BodyType_bits;
extern const uint64_t _BodyType_code;
extern const uint64_t _BodyTypeEnd_code;
extern const uint64_t _BodyType_key;
extern const uint64_t _getBody;
extern const uint64_t _bitEmptyType;
extern const uint64_t _bitCode;                                
extern const uint64_t _bitEmptyCode;
extern const uint64_t _bitBodyCode;
extern const uint64_t _bitHeadCode;
extern const uint64_t _bitVtlHeadCode;
extern uint64_t _BaseNum_bits;    
extern uint64_t _SeqNum_bits;    
extern uint64_t _BaseNum_code;
extern uint64_t _BaseNum_SeqMask;
extern const uint64_t _Empty_Dir_;
extern const unsigned blocklimit;
extern const unsigned index_shape_len;
extern const float def_alpha;

//========================================================
//The is the section to optimize 25-mer index for mapping
//Structt Hs: String<uint64_t>
//Types of node in Hs: 1.head node and 2.body node
//Head: Headflag[1] = 0|sortFlag[1]|N/A[2]|Pointer[20]| xvalue[40]
//Body: bodyflag[1] = 1|N/A[2]|yvalue[20] |typeCode[1]|sa[40]
static const unsigned XValueBit = 40;

struct HsBase
{
    const unsigned bit;
    const unsigned bodyYBit; 
    const unsigned bodyYBitLen; 
    const unsigned bodyYMask; 
    const unsigned bodyCodeBit;
    const unsigned pointerBit;
    const unsigned pointerBitLen;
    const uint64_t mask;
    const uint64_t pointerMask;
    const uint64_t maxPointer;
    const uint64_t headTypeFlag;
    const uint64_t typeFlag;
    const uint64_t typeFlag2;
    const uint64_t typeMask;
    const uint64_t bodyCodeFlag;
    
    HsBase(bool cerr);
};
extern HsBase _DefaultHsBase;

struct Hs
{
    typedef HType1 ValueType; 
    typedef HType2 ValueBodyType;
   
    bool isHead(uint64_t const &, 
                uint64_t const & = _DefaultHsBase.typeFlag);
    uint64_t MinusX(uint64_t const &, uint64_t const &, 
                    uint64_t const & = _DefaultHsBase.mask);
    void setHsHead(uint64_t &, uint64_t const &, uint64_t const &, 
                   uint64_t const & bit = _DefaultHsBase.pointerBit, 
                   uint64_t const & typeFlag = _DefaultHsBase.typeMask);
    uint64_t makeHsHead(uint64_t const &, uint64_t const &, 
                   uint64_t const & bit = _DefaultHsBase.pointerBit, 
                   uint64_t const & typeFlag = _DefaultHsBase.typeMask);
    uint64_t getHeadX(uint64_t const &, 
                      uint64_t const & = _DefaultHsBase.mask);
    uint64_t getHeadPtr(uint64_t const &, 
                        uint64_t const & = _DefaultHsBase.pointerBit, 
                        uint64_t const & = _DefaultHsBase.pointerMask);
    void setHsBody(uint64_t &, uint64_t const &,  uint64_t const & id, uint64_t const & pos,
                   uint64_t const & typeFlag = _DefaultHsBase.typeFlag);
    uint64_t makeHsBody(uint64_t const &,  uint64_t const & id, uint64_t const & pos,
                   uint64_t const & typeFlag = _DefaultHsBase.typeFlag);
    uint64_t getHsBodyY(uint64_t const &,
        uint64_t const & = _DefaultHsBase.bodyYBit, 
        uint64_t const & = _DefaultHsBase.bodyYMask);
    uint64_t getHsBodyS(uint64_t const & val, 
        uint64_t const & mask = _DefaultHsBase.mask); 
    void setHsHeadPtr(uint64_t &, uint64_t const &, 
                      uint64_t const & = _DefaultHsBase.bit, 
                      uint64_t const & = _DefaultHsBase.mask);
    bool isBody(uint64_t const & val, uint64_t const & flag = _DefaultHsBase.typeFlag);
    bool isBodyYEqual(uint64_t const & hval, uint64_t const & yval, 
                    uint64_t const & bit = _DefaultHsBase.bodyYBit,
                    uint64_t const & flag = _DefaultHsBase.typeFlag2);
    void setHsBodyReverseStrand(uint64_t & val);
    bool isHsBodyReverseStrand(uint64_t & val);
    void setHsBodyY(uint64_t & val, uint64_t y, uint64_t const & bit = _DefaultHsBase.bodyYBit, uint64_t const & mask = _DefaultHsBase.bodyYMask);
    
};
extern Hs _DefaultHs;

struct XNode
{
    typedef uint64_t TypeV1;
    typedef unsigned TypeV2;
    typedef uint64_t TypeV2L;
    
    TypeV1 val1;
    TypeV2 val2;
};

struct XString
{
    String<XNode> xstring;
    uint64_t mask;
    
    XString(){};
    XString(uint64_t const & seqlen);
    uint64_t _fullSize(uint64_t const & seqlen, float const & alpha = def_alpha);
    void clear();
};

class HIndex
{
public:

    String<HType1>     ysa;
    XString            xstr;
    LShape             shape;
    double             alpha;    
    uint64_t emptyDir;

    HIndex();
    HIndex(StringSet<String<Dna5> > const & text);
    void clear();
}; 

typedef HIndex LIndex;

uint64_t _getSA_i1(uint64_t const & node);
uint64_t _getSA_i2(uint64_t const & node);
uint64_t getXDir(HIndex const & index, uint64_t const & xval, uint64_t const & yval);
uint64_t getXYDir(HIndex const & index, uint64_t const & xval, uint64_t const & yval);
bool createHIndex(StringSet<String<Dna5> > & seq, LIndex & index, unsigned & threads, bool efficient);

// uint64_t getXDir(LIndex const & index, uint64_t const & xval, uint64_t const & yval)

#endif