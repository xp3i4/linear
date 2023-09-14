#ifndef LINEAR_HEADER_INDEX_UTIL_H
#define LINEAR_HEADER_INDEX_UTIL_H
#include <seqan/parallel.h>
#include <seqan/sequence.h>
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

/*=============================================
=            StringSet for dindex           =
=============================================*/
struct CreatIndexParms : public Parms
{
    int64_t thd_shape_len;
    int64_t thd_min_step;
    int64_t thd_max_step;
    int64_t thd_omit_block;
};
struct CreateSIndexParms : public CreatIndexParms
{
    CreateSIndexParms();
};
struct CreateDIndexParms : public CreatIndexParms
{
    CreateDIndexParms();
};
struct CreateHIndexParms : public CreatIndexParms
{
    CreateHIndexParms();
};

class SIndex
{
    std::vector<std::vector<int64_t> > hs; 
    LShape shape;
public:
    int cas_key;
    String<bool> cas_keys;
    SIndex();
    SIndex(unsigned); //shape_len
    std::vector<std::vector<int64_t> > & getHs();
    LShape & getShape();
    int fullSize();
    int getShapeLen();
    int64_t getVal(int64_t, int64_t);
    int printStatistics();
    int64_t queryHsStr(int64_t);
    int64_t queryHsEnd(int64_t);
    int64_t queryHsBlockLen(int64_t);
    void clear();
};
int createSIndex(StringSet<String<Dna5> > & seqs, 
                 SIndex & index, 
                 unsigned gstr,
                 unsigned gend,
                 int64_t thd_min_step, 
                 int64_t thd_max_step,
                 int64_t thd_omit_block,
                 unsigned threads
                );

/*=============================================
=  short mer direct index for error rate > 0.2 =
=============================================*/

class DIndex 
{
    String <int> dir;
    String <uint64_t> hs;
    void * pt_dir[2];
    void * pt_hs[2];
    int pt;
    uint64_t rs;
    uint64_t fs;
    LShape shape;
public:
    Iterator<String<Dna5> >::Type tmp_it; 
    DIndex();
    DIndex(unsigned); //shape_len
    String<int> & getDir();
    String<uint64_t> & getHs();
    uint64_t val2Anchor(int64_t &, uint64_t & y, uint64_t & read_len, LShape & shape);
    LShape & getShape();
    int fullSize();
    int getShapeLen();
    void clear();
}; 
uint64_t shape2DIndexCordy(LShape & shape);
uint64_t getDIndexCordy(uint64_t index_val);
int createDIndex(StringSet<String<Dna5> > & seqs, 
                 DIndex & index, 
                 unsigned gstr,
                 unsigned gend,
                 int64_t thd_min_step, 
                 int64_t thd_max_step,
                 int64_t thd_omit_block,
                 unsigned threads
                );
int64_t queryHsStr(DIndex & index, uint64_t & xval);
int64_t queryHsEnd(DIndex & index, uint64_t & xval);
void queryHsStrEnd(DIndex & index, uint64_t & xval, std::pair<int64_t, int64_t> & str_end);

/*=============================================
=        Section to optimize 25-mer index     =
=============================================*/
//Structt Hs: String<uint64_t>
//Types of node in Hs: 1.head node and 2.body node
//Head: Headflag[1] = 0|sortFlag[1]|N/A[2]|Pointer[20]| xvalue[40]
//Body: bodyflag[1] = 1|N/A[2]|yvalue[20] |typeCode[1]|sa[40]

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
    void setHsBodyY(uint64_t & val, uint64_t y, uint64_t const & bit = _DefaultHsBase.bodyYBit, uint64_t const & mask = _DefaultHsBase.bodyYMask);
    uint64_t getHsBodyStrand(uint64_t & val);
    uint64_t getLength(String<uint64_t> & ); //length equals num of head of ptr != 0 and body excludes the addition info empty Head nodes ptr==0 at the end of the array.
    
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

    String<HType1>                          ysa;
    String<HType1>                          ysa_str_end;
    String<std::pair<uint64_t, uint64_t> >  ysa_sorted_records;
    XString            xstr;
    LShape             shape;
    double             alpha;    
    uint64_t           emptyDir;

    HIndex();
    HIndex(unsigned shape_len, float index_alpha);
    bool isEmptyDir(uint64_t);
    bool ifYsaSorted(uint64_t g_str, uint64_t g_end);
    bool insertYsaSortedRecord(uint64_t g_str, uint64_t g_end);
    void clear();
    LShape & getShape();
}; 

typedef HIndex LIndex; 
uint64_t _getSA_i1(uint64_t const & node);
uint64_t _getSA_i2(uint64_t const & node);
uint64_t getXDir(HIndex const & index, uint64_t const & xval, uint64_t const & yval);
uint64_t getXYDir(HIndex const & index, uint64_t const & xval, uint64_t const & yval);

bool createHIndex(StringSet<String<Dna5> > & seq, 
                  LIndex & index, 
                  unsigned gstr,
                  unsigned gend,
                  unsigned & threads, 
                  unsigned thd_step, 
                  bool efficient);

/*===============================================
=            Section of IndexDynamic            =
===============================================*/

extern int const typeDIx;
extern int const typeSIx;
extern int const typeHIx;
extern int const typeMIx;
extern int const typeFIx;
struct IndexDynamic 
{
    HIndex hindex;
    DIndex dindex;
    SIndex sindex;
    int typeIx;
    int isHIndex();
    int isSIndex();
    int isDIndex();
    int isMIndex();
    void setHIndex();
    void setDIndex();
    void setSIndex();
    void setMIndex();   //mapper index
    void setMHIndex();  //mapper::HIndx
    void setMDIndex();  //mapper::DIndex
    void setFIndex();   //filter index
    void setFDIndex();  //filter::Dindex
    void setFHIndex();  //filter::HIndex;
    void setIndexType(int);
    void clearIndex();
    IndexDynamic(StringSet<String<Dna5> > &);
};

bool createIndexDynamic(StringSet<String<Dna5> > & seq, 
                        IndexDynamic & index, 
                        unsigned gstr,
                        unsigned gend,
                        unsigned threads, 
                        bool efficient);
// uint64_t getXDir(LIndex const & index, uint64_t const & xval, uint64_t const & yval)

/*=====  End of Section comment block  ======*/
#endif