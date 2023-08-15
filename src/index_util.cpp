#include "base.h"
#include "cords.h"
#include "shape_extend.h"
#include "index_util.h"

using namespace seqan; 

/**
 * definition of types of elements in dir (uint64_t)
 * bodyNode=
 * length1[4] length2[38] YValue[20] code[2]: code={1}
 * 
 * headNode=
 * hvalue(/XValue)[62] code[2]: code={2,3}, 2:head, 3:virtual head
 * 
 * code  0 empty
 * code  2 head 10
 * code  3 virtual head 11
 * code  1 body  01
 */
const uint64_t _dirEmpty = (uint64_t)1 << 63;
const uint64_t _bitEmpty = 0;
const uint64_t _bitLength = ((uint64_t)1 << 58) - 1;
const uint64_t _bitValue = ((uint64_t)1 << 56) - 1;
const unsigned _bitLength_END = 60;
const unsigned _bitValue_END = 2;
//HeadNode(H):
//H :=  (h/X)Value[62]|HeadType[2]: 0:=empty 2:=head 3:=virtual head
const unsigned _HeadValue_bits = 3;
const uint64_t _HeadType_code = 2;
const uint64_t _HeadTypeVtl_code = 3;
const uint64_t _HeadTypeHVl_code = 4;
// BodyNode(B):
// B := YValue[23]|BodyType[1]: |counth[40]
// occ = counth[n+1] - count[n] = (B[n+1] - B[n]) & bit, bit = 00..0011...11
const unsigned _BodyValue_bits = 41;
const unsigned _BodyType_bits = 40;
const uint64_t _BodyType_code = 1;
const uint64_t _BodyTypeEnd_code = 0;
const uint64_t _BodyType_key = ~((uint64_t)1 << _BodyType_bits);
const uint64_t _getBody = ((uint64_t)1 << _BodyType_bits) - 1;
const uint64_t _bitEmptyType = 0;
const uint64_t _bitCode = 3;    
const uint64_t _bitEmptyCode = 0;
const uint64_t _bitBodyCode = 1;
const uint64_t _bitHeadCode = 2;
const uint64_t _bitVtlHeadCode = 3;
//SA node:= seq num i1[10]| base num i2[30]
uint64_t _BaseNum_bits = 30 ;    
uint64_t _SeqNum_bits = _BodyType_bits - _BaseNum_bits;    
uint64_t _BaseNum_code = ((uint64_t)1 << _BaseNum_bits) - 1;
uint64_t _BaseNum_SeqMask = (1ULL << _SeqNum_bits) - 1;
const uint64_t _Empty_Dir_ = -1;
const unsigned blocklimit = 32;
    
 uint64_t _makeHeadNode(uint64_t code)
{
    return (code << _HeadValue_bits) + _HeadType_code;
} 
 uint64_t _makeVtlHeadNode(uint64_t code)
{
    return (code << _HeadValue_bits) + _HeadTypeVtl_code;
}
 uint64_t _makeHVlHeadNode(uint64_t code)
{
    return (code << _HeadValue_bits) + _HeadTypeHVl_code; 
}
 void _setHVlHeadNode(uint64_t & headNode, uint64_t const & hValue)
{
    headNode = (hValue << _HeadValue_bits) + _HeadTypeHVl_code; 
}
 void _setHeadNode(uint64_t & headNode, uint64_t const & hValue)
{
    headNode = (hValue << _HeadValue_bits) + _HeadType_code;
}
 uint64_t _makeEmptyNode(uint64_t code)
{
    return (code << _HeadValue_bits) + _bitEmptyType;
}
 void _setBodyType_Begin(uint64_t & code){
    code &= _BodyType_key; 
}
 uint64_t _ifBodyType(uint64_t code){
    return code & (~_BodyType_key);
}
 uint64_t _getHeadValue(uint64_t  code)
{
    return code >> _HeadValue_bits;  
}
 void _setBodyNode(uint64_t & bodyNode, uint64_t const & YValue, uint64_t const & type, uint64_t const & counth)
{
    bodyNode = (YValue << _BodyValue_bits) + (type << _BodyType_bits) + counth;
}
 uint64_t _getBodyValue(uint64_t code)
{
    return code >> _BodyValue_bits;
}
 uint64_t _getBodyCounth(uint64_t & code)
{
    return code & _getBody;
}
 uint64_t _createSANode(uint64_t const & i1, uint64_t const & i2)
{
    return (i1 << _BaseNum_bits) + i2;
}
 void _setSANode(uint64_t & node, uint64_t const & i1, uint64_t const & i2)
{
    node = (i1 << _BaseNum_bits) + i2;
}
 uint64_t _getSA_i1(uint64_t const & node)
{
    return (node >> _BaseNum_bits) & _BaseNum_SeqMask;
}
 uint64_t _getSA_i2(uint64_t const & node)
{
    return node & _BaseNum_code;
}

/*=====================================================
=            Ehanced open addressing index            =
=======================================================
This is the optimized version of 25-mer index deriverd 
from  the original index.
This version index is specifically optimized for calling 
Approximate mapping. 
Sufficient benchmarks are required before wrapping.  
Do not use it for any other purpose!
*/
//Structt Hs: String<uint64_t>
//Types of node in Hs: 1.head node and 2.body node
//Head: Headflag[1] = 0|sortFlag[1]|N/A[2]|Pointer[20]| xvalue[40]
//Body: bodyflag[1] = 1|N/A[2]|yvalue[20] |typeCode[1]|sa[40]
//WARN::ptr==0 only when it's at the end of the arry which are for additional 
//infos and exclued from the normal head and body nodes.
const unsigned XValueBit = 40;
HsBase::HsBase(bool cerr):
        bit(XValueBit),
        bodyYBit(_BodyValue_bits), // 41
        bodyYBitLen(20),
        bodyYMask((1ULL << bodyYBitLen) - 1),
        bodyCodeBit(_BodyType_bits), //40
        pointerBit(bit), 
        pointerBitLen(23),
        mask((1ULL << bit) - 1),
        pointerMask((1ULL << (pointerBitLen)) - 1),
        maxPointer(1 << pointerBitLen),
        headTypeFlag(0),
        typeFlag(1ULL << 63),
        typeFlag2(1ULL << (63 - bodyYBit)),
        typeMask(typeFlag - 1),
        bodyCodeFlag(1ULL << bodyCodeBit)    // complement reverse strand flag
        {
            if (cerr)
                std::cerr << "HsBase::pointerBit " << pointerBit << std::endl;
        }
HsBase _DefaultHsBase(false);
/**
 * XNode = struct v1[64],v2[32]: .v1: hashvalue; .v2:pointer to YNode
 * .v1: v2 type[2]|value[60]|v1 type[2]
 * Types of .v1 including: 1.empty node{00} 2.xvalue head{01} 3.hvalue{10} head 4 virtual ad{11}
 *  0 empty
 *  1 xvalue head
 *  2 head 10
 *  3 virtual head 
 */
struct XNodeBase   //define dirNode
{
    typedef uint64_t NodeType;

    uint64_t bit;
    uint64_t mask;
    uint64_t bit2;
    uint64_t mask2;
    
    NodeType emptyNode;
    NodeType xHead;
    NodeType xHead1;
    NodeType virtualHead;
    
    uint64_t returnDir;
    uint64_t returnSa;
   
    NodeType _Empty_Dir_;
    
    XNodeBase():
        bit(2),
        mask((1<<bit) - 1),
        bit2(62),
        mask2((1 << bit2) - 1),
        emptyNode(0),
        xHead(1),
        xHead1(2),
        virtualHead(3),
        returnDir(0),
        returnSa(1),
        _Empty_Dir_(~0)
    {}
}_DefaultXNodeBase;

struct XNodeFunc
{
    uint64_t getAddY();
    uint64_t hash(uint64_t const &);
    void setXNode(XNode &, XNode::TypeV1 const & val1, XNode::TypeV2 const &, 
                  uint64_t const &, uint64_t const &, 
                  uint64_t const & = _DefaultXNodeBase.bit,
                  uint64_t const & = _DefaultXNodeBase.bit2
                 );
    XNode::TypeV1 makeYXKey(typename Hs::ValueBodyType const &, XNode::TypeV1 const &, 
                            uint64_t const & = _DefaultHsBase.bodyYMask << _DefaultHsBase.bodyYBit, 
                            uint64_t const & = _DefaultHsBase.mask);
    XNode::TypeV1 collision(XNode::TypeV1 const &, XNode::TypeV1 const &, uint64_t const & = _DefaultXNodeBase.mask2);
    XNode::TypeV2L makeReturnVal(XNode const &, uint64_t const & = _DefaultXNodeBase.mask2);
    
}_DefaultXNodeFunc;

XString::XString(uint64_t const & seqlen)
{
    _fullSize(seqlen);
}
uint64_t XString::_fullSize(uint64_t const & seqlen, float const & alpha)
{
    uint64_t len = 1ULL; 
    while ((len) < seqlen * alpha)
        len <<=1;
    //std::cout << "seql" << len << "\n";
    resize(xstring, len);
    mask = len - 1;
    for (uint64_t k = 0; k < len; k++)
        xstring[k].val1 = 0;
    return len;
}
void XString::clear()
{
    seqan::clear(xstring);
    shrinkToFit(xstring);
}
bool Hs::isHead(uint64_t const & val, uint64_t const & flag)
{
    return (val & flag) ^ flag;
}
void Hs::setHsHead(uint64_t & head, uint64_t const & ptr, uint64_t const & xval, uint64_t const & bit, uint64_t const & mask)
{
    head = ((ptr << bit) + xval) & mask;
}
uint64_t Hs::makeHsHead(uint64_t const & ptr, uint64_t const & xval, uint64_t const & bit, uint64_t const & mask)
{
    return ((ptr << bit) + xval) & mask;
}
uint64_t Hs::MinusX(uint64_t const & value1, uint64_t const & value2, uint64_t const & mask)
{
    return ((value1 - value2) & mask);
}
uint64_t Hs::getHeadX(uint64_t const & value, uint64_t const & mask)
{
    return value & mask;
}  
uint64_t Hs::getHeadPtr(uint64_t const & val, uint64_t const & bit, uint64_t const & mask)
{
    return (val >> bit) & mask;
}
void Hs::setHsBody(uint64_t & val, uint64_t const & yval, uint64_t const & id, uint64_t const & pos, uint64_t const & typeFlag)
{
    val = ((yval << _BodyValue_bits)|typeFlag) + (id << _BaseNum_bits) + (pos);
    
}
uint64_t Hs::makeHsBody(uint64_t const & yval, uint64_t const & id, uint64_t const & pos, uint64_t const & typeFlag)
{
    return ((yval << _BodyValue_bits)|typeFlag) + (id << _BaseNum_bits) + (pos);
    
}
uint64_t Hs::getHsBodyY(uint64_t const & val, uint64_t const & bit, uint64_t const & mask)
{
    return (val >> bit) & mask;
}
void Hs::setHsHeadPtr(uint64_t & val, uint64_t const & ptr,  uint64_t const & bit, uint64_t const & mask)
{
    val = (val & mask) + (ptr << bit);
}
uint64_t Hs::getHsBodyS(uint64_t const & val, uint64_t const & mask )
{
    return val & mask;
}
bool Hs::isBody(uint64_t const & val, uint64_t const & flag)
{
    return val & flag;
}
bool Hs::isBodyYEqual(uint64_t const & hval, uint64_t const & yval, uint64_t const & bit, uint64_t const & flag)
{
    return ((hval >> bit) ^ yval) == flag;
}
void Hs::setHsBodyReverseStrand(uint64_t & val)
{
    val |= _DefaultHsBase.bodyCodeFlag;
}
void Hs::setHsBodyY(uint64_t & val, uint64_t y, uint64_t const & bit, uint64_t const & mask)
{
    val = (val & (~(mask << bit))) | (y << bit);
}
uint64_t Hs::getHsBodyStrand(uint64_t & val)
{
    return (val >> _DefaultHsBase.bodyCodeBit) & 1ULL;
}
uint64_t Hs::getLength(String<uint64_t> & hs)
{
    uint64_t real_len = length(hs);
    while (real_len > 0 && (isHead(hs[real_len - 1]) && !getHeadPtr(hs[real_len - 1])))
    {real_len--;}
    return real_len;
}

Hs _DefaultHs;
const unsigned index_shape_len = 25; 
const float def_alpha = 1.6;
HIndex::HIndex():
    shape(index_shape_len),
    alpha(def_alpha)
    {}
HIndex::HIndex(unsigned shape_len, float val):
    shape(shape_len),
    alpha(val)
{
}
LShape & HIndex::getShape()
{
    return shape; 
}
bool HIndex::isEmptyDir(uint64_t pos)
{
    return pos == emptyDir;
}
void HIndex::clear()
{
    seqan::clear(ysa);
    shrinkToFit(ysa);
    xstr.clear();
}
/* Chenck if region in ysa contianing genomes from @g_str to @g_end [@g_str, g_end) 
   have been already sorted;
 * New records are required always append at the end of @ysa_sort_records. So the function checks starting from 
   the end of the @ysa_sort_records to check the latestd record that is ovelapped with [@g_str, g_end).
 * Return true only if [@g_str, g_end) is exactly in the @ysa_sorted_records !
 */
bool HIndex::ifYsaSorted(uint64_t g_str, uint64_t g_end)
{
    for (int i = length(ysa_sorted_records) - 1; i >= 0; i--)
    {
        if (std::max(g_str, ysa_sorted_records[i].first) < std::min(g_end, ysa_sorted_records[i].second)) //overlap
        {
            return (g_str != ysa_sorted_records[i].first || g_end != ysa_sorted_records[i].second) ?
                    false : true;
        }
    }
    return false;
}

bool HIndex::insertYsaSortedRecord(uint64_t g_str, uint64_t g_end)
{
    bool f_return = ifYsaSorted(g_str, g_end);
    if (!f_return)
    {
        appendValue (ysa_sorted_records, std::pair<uint64_t, uint64_t> (g_str, g_end));
    }
    return f_return;
}

/*
 * serial sort hs
 * bucket[]+
 */
 bool _hsSortX(Iterator<String<uint64_t> >::Type const & begin, 
                     Iterator<String<uint64_t> >::Type const & end, 
                     unsigned const & xValBitLen)
{
    if (xValBitLen <34 || xValBitLen > 42)
    {
        std::cerr << "[Error]: _dirSortX " << xValBitLen << "\n";
        return false;
    }
    
    unsigned bit[18] = {9,4,9,4,9,4,8,5,8,5,8,5,8,5,7,6,7,6}; //xValueBitLen 34 - 42;
    
    unsigned p_bit = bit[(xValBitLen - 34) << 1];
    unsigned l =  bit[((xValBitLen - 34) << 1) + 1];
    //std::cerr << p_bit << " " << l << std::endl;
    unsigned  l_move = 64, r_move = 64 - p_bit;
    uint64_t count[513]; // 2^max(bit[18]) + 1
    //int count[1024];
    
    String<uint64_t> output;
    resize(output, end - begin);
    //std::cerr << "end - begin " <<end - begin << std::endl;
    for (uint64_t j = 0; j < l; j++)
    {
        l_move -= p_bit;
        for (int k = 0; k < (1<<p_bit); k++)
            count[k]=0;
        for (int64_t k = 0; k < end - begin; k += _DefaultHs.getHeadPtr(*(begin + k)))
        {
            count[(*(begin + k) << l_move >> r_move) + 1] += _DefaultHs.getHeadPtr(*(begin + k));
        }
        for (int k = 1; k < (1 << p_bit); k++)
        {
            count[k] += count[k - 1];
        }
        for (int64_t k = 0; k < end - begin; k++)
        {
        
            if (_DefaultHs.isHead(*(begin + k)))
            {
                uint64_t x = *(begin + k) << l_move >> r_move;
                uint64_t ptr = _DefaultHs.getHeadPtr(*(begin + k));
                for (uint64_t it = 0; it < ptr; it++)
                {
                    output[count[x] + it] = *(begin + k + it);
                }
                count[x] += ptr;
            }
            
        }
        for (int64_t k = 0; k < end - begin; k++)
            *(begin + k) = output[k];
    }
    return true;
}
/*
 * parallel sort hs for index either collecting kmers or collecting minimizers
 * bucket[]+
 * However requires larger memory footprint.
 * xValBitLen = shape.weight * 2
 */
 bool _hsSortX_1(Iterator<String<uint64_t> >::Type const & begin, 
                       Iterator<String<uint64_t> >::Type const & end, 
                       unsigned const & xValBitLen, 
                       unsigned threads)
{

    const unsigned lowerBound = 10;
    const unsigned upperBound = 42;
    if (xValBitLen < lowerBound || xValBitLen > upperBound)
    {
        std::cerr << "[Error]: _dirSortX " << xValBitLen << "\n";
       // return false;
    }
    
    //unsigned const bit[18] = {9,4,9,4,9,4,8,5,8,5,8,5,8,5,7,6,7,6}; //xValueBitLen 34 - 42;
    unsigned const bit[upperBound - lowerBound + 2] 
        = {5,2,4,3,7,2,4,4,9,2,10,2,11,2,12,2,7,4,7,4,10,3,8,4,9,4,9,4,8,5,8,5,7,6}; //xValueBitLen 34 - 42;
    unsigned const p_bit = bit[(xValBitLen - lowerBound + 1) >> 1 << 1];
    unsigned const l =  bit[((xValBitLen - lowerBound + 1) >> 1 << 1) + 1];
    unsigned const r_move = 64 - p_bit;
    unsigned l_move = 64;
    uint64_t const mask = (1 << p_bit) - 1;
    uint64_t size = (end - begin) / threads;
    unsigned thd1 = end - begin - size * threads;
    uint64_t thd_n1 = (size + 1) * thd1;
    std::vector<std::vector<uint64_t> > ctd(threads, std::vector<uint64_t>((1<<p_bit) + 1, 0));
    std::vector<std::vector<std::vector<uint64_t> > > next(threads, 
                    std::vector<std::vector<uint64_t> >(threads, std::vector<uint64_t>((1 << p_bit) + 1, 0)));
    String<uint64_t> output;
    resize(output, end - begin);
    
    #pragma omp parallel 
    {
        unsigned thd_id = omp_get_thread_num();
        #pragma omp for
        for (int64_t k = 0; k < end - begin; k++)
        {
            if (_DefaultHs.isHead(*(begin + k)))
            {
                uint64_t x = *(begin + k) & mask;
                uint64_t ptr = _DefaultHs.getHeadPtr(*(begin + k));
                if (thd_id == threads - 1)
                    ctd[0][x + 1] += ptr;    
                    //ctd[x+1][0] += ptr;
                else
                    ctd[thd_id + 1][x] += ptr;
                    //ctd[x][thd_id + 1] += ptr;
            }
        }
    }

    unsigned PBit = 1 << p_bit;
    for (uint64_t j = 0; j < l; j++)
    {
        std::cerr << "=>Index::SortHash             "  << (float) j / l * 100 << "%\r";
        l_move -= p_bit;
        for (unsigned m = 1; m < threads; m++)
            ctd[m][0] += ctd[m - 1][0];
            //ctd[0][m] += ctd[0][m - 1];
        for (unsigned n = 1; n < PBit; n++)
        {
            ctd[0][n] += ctd[threads - 1][n - 1];
            //ctd[n][0] += ctd[n-1][threads - 1];
            for (unsigned m=1; m < threads; m++)
            {
                ctd[m][n] += ctd[m-1][n];
                //ctd[n][m] += ctd[n][m-1];
            }
        }
        #pragma omp parallel 
        {
            unsigned thd_id = omp_get_thread_num();
            #pragma omp for
            for (int64_t k = 0; k < end - begin; k++)
            {
                if (_DefaultHs.isHead(*(begin + k)))
                {
                    uint64_t x = *(begin + k) << l_move >> r_move;
                    uint64_t ptr = _DefaultHs.getHeadPtr(*(begin + k));
                    for (uint64_t it = 0; it < ptr; it++)
                    {
                        output[ctd[thd_id][x] + it] = *(begin + k + it);
                        //output[ctd[x][thd_id] + it] = *(begin + k + it);
                    }
                    unsigned thd_num = (ctd[thd_id][x] < thd_n1) ? 
                                       (ctd[thd_id][x]) / (size + 1) : 
                                       (ctd[thd_id][x] - thd_n1) / size + thd1;
                    ctd[thd_id][x] += ptr;
                    x = *(begin + k) << (l_move - p_bit)>> r_move;
                    if (thd_num == threads - 1)
                    {
                        next[thd_id][0][x + 1] += ptr;
                    }
                    else
                    {
                        next[thd_id][thd_num + 1][x] += ptr;
                    }
                }
            }
        }
        if (j < l - 1)
        {       
                    #pragma omp parallel for
                for (unsigned k=0; k < threads; k++)
                {
                    //std::fill(ctd[k].begin(), ctd[k].end(), 0);
                    for (unsigned n = 0; n < PBit; n++)
                    {
                        ctd[k][n] = 0;
                        //ctd[n][k] = 0;
                        for (unsigned m = 0; m < threads; m++)
                        {
                            ctd[k][n] += next[m][k][n];
                            //ctd[n][k] += next[m][k][n];
                            next[m][k][n] = 0;
                        }   
                    }
                        
                }
        }
        #pragma omp parallel for
        for(int64_t k = 0; k < end - begin; k++)
        {
            *(begin + k) = output[k];
        }
        
    }

    return true;
}

/*
 * Interface to sort x
 */
 bool _hsSortX(Iterator<String<uint64_t> >::Type const & begin, 
                     Iterator<String<uint64_t> >::Type const & end, 
                     unsigned const & xValBitLen, 
                     unsigned threads)
{
    //_hsSortX(begin, end, xValBitLen);
    _hsSortX_1(begin, end, xValBitLen, threads);
    //_hsSortX_2(begin, end, xValBitLen, threads);
    return true;
}

 void insertSort(Iterator<String<uint64_t> >::Type const & begin, 
                       Iterator<String<uint64_t> >::Type const & end)
{
    typedef Iterator<String<uint64_t> >::Type TIter;
    typename Value<TIter>::Type key;
    
    for (int j = 1; j < end - begin; j++)
    {
        key = *(begin + j); 
        int k = j - 1;
        while (k >= 0)
        {
            if (*(begin + k) < key)
                *(begin+k+1) = *(begin+k);
            else
            {   
                break;
            }   
            k--;
        }        
        *(begin+k+1) = key;
    }
}

 bool _sort_YSA_Block(Iterator<String<uint64_t> >::Type const & begin, 
                            Iterator<String<uint64_t> >::Type const & end, 
                            unsigned const & sortModeThr = 60) // sort y and sa
{
    typedef typename Iterator<String<uint64_t> >::Type TIter;
    typedef typename Value<TIter>::Type ValueType;
    if (end - begin< sortModeThr)
        insertSort(begin, end);
    else
        std::sort(begin, end, std::greater<ValueType>());
    return true;
}

 bool _hsSortY_SA(Iterator<String<uint64_t> >::Type const & begin, 
                  Iterator<String<uint64_t> >::Type const & end) 
{
    typedef typename Iterator<String<uint64_t> >::Type TIter;
    typedef typename Value<TIter>::Type ValueType;
    uint64_t k = 0, ptr;
    unsigned sortModeThr = 20;
    while ((int64_t)k < end - begin)
    {
        ptr = _DefaultHs.getHeadPtr(*(begin + k));
        if (ptr< sortModeThr)
            insertSort(begin + k + 1, begin + k + ptr);
        else
            std::sort(begin + k + 1, begin + k + ptr, std::greater<ValueType>());
        k += ptr;
    }
    return true;
}

 void _hsSort(Iterator<String<uint64_t> >::Type const & begin, 
                    Iterator<String<uint64_t> >::Type const & end, 
                    unsigned const & shapeWeight)
{
    _hsSortX(begin, end, shapeWeight << 1);
}

 void _hsSort(Iterator<String<uint64_t> >::Type const & begin, 
                    Iterator<String<uint64_t> >::Type const & end, 
                    unsigned const & shapeWeight, unsigned & threads)
{
    _hsSortX(begin, end, shapeWeight << 1, threads);
}

/*
 * serial creat hash array
 */
 bool _createHsArray(StringSet<String<Dna5> > & seq, 
                     String<uint64_t> & hs, 
                     LShape & shape,
                     unsigned gstr,
                     unsigned gend
                     )
{
    double time = sysTime();
    uint64_t preX = ~0;
    //-k
    //int64_t ptr = 0, count = 1;
    int64_t ptr = 0, count = -1;
    resize(hs, lengthSum(seq) << 1);
    //-k
    //unsigned start = 2;
    uint64_t thd_step = 3;
    hs[0] = hs[1] = 0;
    for(uint64_t j = gstr; j < gend; j++)
    {
        hashInit(shape, begin(seq[j]));
        for (uint64_t k =0; k < length(seq[j]) - shape.span + 1; k++)
        {
            if(ordValue(*(begin(seq[j]) + k + shape.span - 1)) == 4)
            {
                k += hashInit(shape, begin(seq[j]) + k);
                if(k > length(seq[j]) - shape.span + 1)
                    break;
            }
            
            hashNext(shape, begin(seq[j]) + k);
            if (k % thd_step != 0)
            {
                if (shape.XValue ^ preX)
                {
                    _DefaultHs.setHsHead(hs[++count - ptr], ptr, preX);
                    ptr = 2;
                    preX = shape.XValue; 
                }
                else
                {
                    ++ptr;
                }
                _DefaultHs.setHsBody(hs[++count], shape.YValue, j, k); 
                //std::cout << "hs[count] " << (hs[count] & ((1ULL << 30) - 1))<< std::endl;

            }
            
        }
    }
    _DefaultHs.setHsHead(hs[++count - ptr], ptr, shape.XValue);
    _DefaultHs.setHsHead(hs[count], 0, 0);
    std::cerr << "      init Time[s]" << sysTime() - time << " " << std::endl;
//-k
    resize(hs, count + 1);
    //shrinkToFit(hs);
    _hsSort(begin(hs), begin(hs) + count, shape.weight);
    //_hsSort(begin(hs) + start, begin(hs) + count, shape.weight);
    
    std::cerr << "      End createHsArray " << std::endl;
    return true;
}



/*
 * parallel creat hash array
 * creating index only collecting mini hash value [minindex]
 * state::warnning. for seq contains 'N', error. since the k in openmp doesn't change correctly
 */
 uint64_t __createHsArray(StringSet<String<Dna5> > & seq, 
                     String<uint64_t> & hs, 
                     LShape & shape, 
                     unsigned gstr,
                     unsigned gend,
                     unsigned threads,
                     unsigned thd_step)
{
    uint64_t hsRealEnd = 0;
    resize (hs, lengthSum(seq) * 2 / thd_step + 1000);
    std::vector<int64_t> hsRealSize(threads, 0);
    std::vector<int64_t> seqChunkSize(threads, 0);
    std::vector<int64_t> hss(threads, 0);
    //uint64_t seqChunkSize[4] = {0};
    for(uint64_t j = gstr; j < gend; j++)
    {
        std::cerr << "=>Index::Initiate             " << (float)j / length(seq) * 100<< "%           \r";
        uint64_t thd_count = 0; // count number of elements in hs[] for each thread
        #pragma omp parallel reduction(+: thd_count)
        {
            LShape tshape = shape; 
            uint64_t preX = ~0;
            int64_t ptr = 0;
            uint64_t size2 = (length(seq[j]) - tshape.span + 1) / threads;
            uint64_t start;
            uint64_t hsStart; 
            unsigned thd_id = omp_get_thread_num();
            //unsigned ct_step = 2;
            if (thd_id < (length(seq[j]) - tshape.span + 1) - size2 * threads)
            {
                seqChunkSize[thd_id] = size2 + 1;
                start = (size2 + 1) * thd_id;
                hsStart = hsRealEnd + (start << 1) / thd_step + thd_id * 10;
                hss[thd_id] = hsStart;
            }
            else
            {
                seqChunkSize[thd_id] = size2;
                start =  length(seq[j]) + 1 - tshape.span - size2 * (threads - thd_id);
                hsStart = hsRealEnd + (start << 1) / thd_step + thd_id * 10;
                hss[thd_id] = hsStart;
            }
 
            hashInit(tshape, begin(seq[j]) + start);
            //dout << "start" << start << seqChunkSize[thd_id] << "\n";
            for (uint64_t k = start; k < start + seqChunkSize[thd_id]; k++)
            {
                //std::cerr << k << "\n";
                if(ordValue(*(begin(seq[j]) + k + tshape.span - 1)) == 4)
                {
                    k += hashInit(tshape, begin(seq[j]) + k);

                    if (k > seqChunkSize[thd_id] - tshape.span + 1 + start)
                    {
                        k = seqChunkSize[thd_id] - (seqChunkSize[thd_id] + start) % thd_step + thd_step + start;
                    }
                }
                hashNext(tshape, begin(seq[j]) + k);
                if (k % thd_step == 0)
                {
                    if (tshape.XValue ^ preX)
                    {
                        _DefaultHs.setHsHead(hs[hsStart + thd_count - ptr], ptr, preX);
                        _DefaultHs.setHsBody(hs[hsStart + ++thd_count], tshape.YValue, j, k); 
                        if (tshape.strand)
                        {
                            _DefaultHs.setHsBodyReverseStrand(hs[hsStart + thd_count]);
                        }
                        preX = tshape.XValue; 
                        ++thd_count;
                        ptr = 2;
                    }
                }
            }
            _DefaultHs.setHsHead(hs[hsStart + thd_count - ptr], ptr, tshape.XValue);
            hsRealSize[thd_id] = thd_count;
        }
        for (unsigned k = 1; k < threads; k++)
        {
            hsRealSize[k] += hsRealSize[k - 1];
            seqChunkSize[k] += seqChunkSize[k - 1];
        }
        for (unsigned j = 1; j < threads; j++)
        {
            //uint64_t it = hsRealEnd + (seqChunkSize[j - 1] << 1);
            uint64_t it = hss[j];
            for (uint64_t k = hsRealEnd + hsRealSize[j - 1]; k < hsRealEnd + hsRealSize[j]; k++)
            {
                hs[k] = hs[it];
                ++it;
            }   
        }   

        hsRealEnd += thd_count;
    }
    resize (hs, hsRealEnd + 1);
    _DefaultHs.setHsHead(hs[hsRealEnd], 0, 0);
    //std::cout << "hs_size " << length(hs) << "\n";
    return hsRealEnd;
}

bool _createHsArray(StringSet<String<Dna5> > & seq, 
                    String<uint64_t> & hs, 
                    LShape & shape, 
                    unsigned gstr,
                    unsigned gend,
                    unsigned threads,
                    unsigned thd_step, 
                    bool efficient = true)
{
    double time = sysTime();
    uint64_t hsRealEnd = __createHsArray(seq, hs, shape, gstr, gend, threads, thd_step);
    if (efficient)
    {
        clear(seq);
        shrinkToFit(seq);   
    }
    std::cerr << "--Index::Initiate             Elapsed Time[s] " << sysTime() - time << "\n";
    std::cerr << "=>Index::SortHash                                                   \r";
    time = sysTime();
    _hsSort(begin(hs), begin(hs) + hsRealEnd, shape.weight, threads);
    std::cerr << "  Index::SortHash             Elapsed Time[s] " << sysTime() - time << std::endl;
    return true;
}

/*
 * parallel creat hash array
 * creating index only collecting mini hash value [minindex]
 * genomes will be detroyed during the functoin to reduce memory consumption
 * appendvalue instead of resize
 * state::debug succ for seq without 'N', seq containing 'N' not tested 
 */
 bool _createHsArray2_MF(StringSet<String<Dna5> >  & seq, String<uint64_t> & hs, LShape & shape, unsigned & threads)
{
    std::cerr << "[prallel_createHsArray2_MF]\n";
    double time = sysTime();
    unsigned const step = 10;
    for(uint64_t j = 0; j < length(seq); j++)
    {
        #pragma omp parallel
        {
            LShape tshape = shape; 
            String<uint64_t> hsTmp;
            clear(hsTmp);
            uint64_t preX = ~0;
            int64_t ptr = 2;
            int64_t kn = 0;
            bool flag = true; // if it is the first k for each thread then true;
            #pragma omp for
            for (uint64_t k = 0; k < length(seq[j]); k++)
            {
                //printf("id %d %d %d\n", omp_get_thread_num(), k, length(seq[j]));
                if ((int64_t)k >= kn)
                {
                    if(ordValue(*(begin(seq[j]) + k + tshape.span - 1)) == 4 || flag)
                    {
                        kn = hashInit(tshape, begin(seq[j]) + k) + k;
                        flag = false;
                    }
                    else
                    {
                        hashNext(tshape, begin(seq[j]) + k);
//!Note in openmp, the value of k can't be changed. 
//  so here define the kn to store the value of k
//function above using ompenmp need to be modified!
                        
                        if (k % step == 0)
                        {
                            if (tshape.XValue ^ preX)
                            {
                                appendValue(hsTmp, _DefaultHs.makeHsHead(ptr, tshape.XValue));
                                appendValue(hsTmp, _DefaultHs.makeHsBody(tshape.YValue, j, k));
                                if (tshape.strand)
                                {
                                    _DefaultHs.setHsBodyReverseStrand(back(hsTmp));
                                }
                                preX = tshape.XValue; 
                            }
                        }
                    }
                }
            }
            #pragma omp for ordered
            for(unsigned tj = 0; tj < threads; tj++)
            #pragma omp ordered
            {
                append(hs, hsTmp);
            }
        }
    }
    clear(seq);
    shrinkToFit(seq);
    appendValue(hs, _DefaultHs.makeHsHead(0, 0));
    std::cerr << "[debug] length of hs " << length(hs)  << " lengthsum " << lengthSum(seq) * 2 / step << "\n";
    std::cerr << "      init Time[s]" << sysTime() - time << " " << std::endl;
    
    _hsSort(begin(hs), end(hs) - 1, shape.weight, threads);
    
    std::cerr << "      End createHsArray " << std::endl;
    return true;
}

bool checkHsSort(String<uint64_t> const & hs)
{
    uint64_t preX = _DefaultHs.getHeadX(hs[0]);
    uint64_t preY = hs[1];
    uint64_t k = 0;
    uint64_t count = 0;
    uint64_t countT = 0;
    double time = sysTime();
    while (_DefaultHs.getHeadPtr(hs[k]))
    {
        if (_DefaultHs.getHeadX(hs[k]) < preX)
        {
            std::cerr << "sort x error " << k;
            return false;
        }
        preY = hs[k + 1];
        for (unsigned j = k + 1; j < k + _DefaultHs.getHeadPtr(hs[k]); j++)
        {
            //std::cout << "hs[j]" << hs[j] << std::endl;
            if (hs[j] > preY)
            {
                //std::cerr << "sort y error " << k << " " << j;
                return false;
            }
            preY = hs[j];
        }
        //std::cout << std::endl;
        preX = _DefaultHs.getHeadX(hs[k]);
        k += _DefaultHs.getHeadPtr(hs[k]);
    }
    //std::cerr << "last " << k << std::endl;
    
    k = 0;
    time = sysTime();
    preX = _DefaultHs.getHeadX(hs[0]);
    while (_DefaultHs.getHeadPtr(hs[k]))
    {
        if (preX != _DefaultHs.getHeadX(hs[k]))
        {
            count++;
            if (preX >= (3ULL<30))
                countT++;
        }
        preX = _DefaultHs.getHeadX(hs[k]);
        k += _DefaultHs.getHeadPtr(hs[k]);
    }
    std::cerr << "checkHsSort() " << count << " " << (float)countT / count << " " << sysTime() - time << std::endl;
    return true;
}

 uint64_t XNodeFunc::hash(uint64_t const & val)
{
    uint64_t key = val;
    key = (~key) + (key << 21); 
    key = key ^ (key >> 24);
    key = (key + (key << 3)) + (key << 8);
    key = key ^ (key >> 14);
    key = (key + (key << 2)) + (key << 4); 
    key = key ^ (key >> 28);
    key = key + (key << 31);
    return key;        
}

 void XNodeFunc::setXNode(XNode & val, XNode::TypeV1 const & val1, XNode::TypeV2 const &val2, uint64_t const & t, uint64_t const & r, uint64_t const & bit, uint64_t const & bit2)
{
    val.val1 = (val1 << bit) + t + (r << bit2);
    val.val2 = val2;
}

 XNode::TypeV1 XNodeFunc::makeYXKey(typename Hs::ValueBodyType const & hval, XNode::TypeV1 const & xval, uint64_t const & masky, uint64_t const & maskx)
{
    return (hval & masky) + (xval & maskx);
}

 XNode::TypeV1 XNodeFunc::collision(XNode::TypeV1 const & val, XNode::TypeV1 const & xnode_val1, uint64_t const & mask)
{
    return (val ^ xnode_val1) & mask;
}

 XNode::TypeV2L XNodeFunc::makeReturnVal(XNode const & xnode, uint64_t const & mask)
{
    return (xnode.val1 & (~mask)) + xnode.val2;
    //return xnode.val2;
}

 uint64_t requestXNode_noCollision (XString & xstr, uint64_t const & xval, 
                uint64_t const & val2,  uint64_t const & nodeType, uint64_t const returnType)
{
    uint64_t h1 = _DefaultXNodeFunc.hash(xval) & xstr.mask;
    uint64_t delta = 0;
    while (xstr.xstring[h1].val1) //0 stands for empty
    {
        h1 = (h1 + delta +1) & xstr.mask;
        delta++;
    }
    _DefaultXNodeFunc.setXNode(xstr.xstring[h1], xval, val2, nodeType, returnType);
    return h1;
}

 uint64_t requestXNode_noCollision_Atomic (XString & xstr, uint64_t const & xval, 
                uint64_t const & val2,  uint64_t const & nodeType, uint64_t const returnType)
{
    uint64_t h1 = _DefaultXNodeFunc.hash(xval) & xstr.mask;
    uint64_t delta = 0;
    //while (xstr.xstring[h1].val1) //0 stands for empty
    while (atomicCas<uint64_t>(xstr.xstring[h1].val1, 0, 1))
    {
        h1 = (h1 + delta +1) & xstr.mask;
        delta++;
        //printf("[debug]::delta %d\n", delta);
    }
    _DefaultXNodeFunc.setXNode(xstr.xstring[h1], xval, val2, nodeType, returnType);
    return h1;
}

uint64_t getXDir(XString const & xstr, uint64_t const & xval, uint64_t const & yval)
{
        uint64_t val, delta = 0;
        uint64_t h1 = _DefaultXNodeFunc.hash(xval) & xstr.mask;
        
        //_setHeadNode(val, xval);
//!!!!! need to modify;
        val = (xval << 2) + _DefaultXNodeBase.xHead;
        while (xstr.xstring[h1].val1)
        {
            //switch (xstr.xstring[h1].val1 ^ val) 
            switch(_DefaultXNodeFunc.collision(xstr.xstring[h1].val1, val))
            {
                case 0:
                    //return _DefaultXNodeFunc.makeReturnVal(xstr.xstring[h1]);
                    return xstr.xstring[h1].val2;
                case 2:
//!!!!! need to modify;
                    val = (yval << 42) + (xval << 2) + _DefaultXNodeBase.xHead;
                    h1 = _DefaultXNodeFunc.hash((yval << 40) + xval) & xstr.mask;
                    delta = 0;
                    //std::cerr << "case2\n" ;
                    break;
                //case 3:
                //    return ( ^ shape.yvalue)?xstr[h1].val2:_DefaultXNodeBase._Empty_Dir_;
                default:
                    //std::cerr << "case4\n" ;
                    h1 = (h1 + delta + 1) & xstr.mask;
                    delta++;
            }
        }
        //return _DefaultXNodeBase._Empty_Dir_;
    return 0;
}

uint64_t getXDir(HIndex const & index, uint64_t const & xval, uint64_t const & yval)
{
    uint64_t val, delta = 0;
    uint64_t h1 = _DefaultXNodeFunc.hash(xval) & index.xstr.mask;
    val = (xval << 2) + _DefaultXNodeBase.xHead;
    while (index.xstr.xstring[h1].val1)
    {
        switch(_DefaultXNodeFunc.collision(index.xstr.xstring[h1].val1, val))
        {
            case 0:
                return index.xstr.xstring[h1].val2;
            case 2:
                val = (yval << 42) + (xval << 2) + _DefaultXNodeBase.xHead;
                h1 = _DefaultXNodeFunc.hash((yval << 40) + xval) & index.xstr.mask;
                delta = 0;
                break;
            default:
                h1 = (h1 + delta + 1) & index.xstr.mask;
                delta++;
        }
    }
    return index.emptyDir;
}

 uint64_t getXYDir(HIndex const & index, uint64_t const & xval, uint64_t const & yval)
{
    uint64_t val, delta = 0;
    uint64_t h1 = _DefaultXNodeFunc.hash(xval) & index.xstr.mask;
    uint64_t pos;
    
    //_setHeadNode(val, xval);
//!!!!! need to modify;
    val = (xval << 2) + _DefaultXNodeBase.xHead;
    //unsigned count = 0;
    while (index.xstr.xstring[h1].val1)
    {
        //count++;
        //printf("[debug]::getXYDir %d\n", count);
        //switch (index.xstr.xstring[h1].val1 ^ val) 
        switch(_DefaultXNodeFunc.collision(index.xstr.xstring[h1].val1, val))
        {
            case 0:
                //std::cerr << "case1\n";
                pos = index.xstr.xstring[h1].val2;
                do{
                    if (yval == _DefaultHs.getHsBodyY(index.ysa[pos]))
                        return pos;
                    if (yval > _DefaultHs.getHsBodyY(index.ysa[pos]))
                        return index.emptyDir;
                }while(_DefaultHs.isBody(index.ysa[++pos]));
                //return _DefaultXNodeFunc.makeReturnVal(index.xstr.xstring[h1]);
            case 2:
//!!!!! need to modify;
                val = (yval << 42) + (xval << 2) + _DefaultXNodeBase.xHead;
                h1 = _DefaultXNodeFunc.hash((yval << 40) + xval) & index.xstr.mask;
                delta = 0;
                //std::cerr << "case2\n" ;
                break;
            //case 3:
            //    return ( ^ shape.yvalue)?index.xstr[h1].val2:_DefaultXNodeBase._Empty_Dir_;
            default:
                //std::cerr << "case4\n" ;
                h1 = (h1 + delta + 1) & index.xstr.mask;
                delta++;
        }
    }
    return index.emptyDir;
}

 uint64_t gethDir(XString const & xstr, uint64_t const & hval)
{        
    uint64_t  val = (hval << 2) + _DefaultXNodeBase.xHead;
    uint64_t  h1 = _DefaultXNodeFunc.hash(hval) & xstr.mask;
    uint64_t delta = 0;
        
    //!!!!! need to modify;
    while (xstr.xstring[h1].val1)
    {
        switch(_DefaultXNodeFunc.collision(xstr.xstring[h1].val1, val))
        {
            case 0:
                //return _DefaultXNodeFunc.makeReturnVal(xstr.xstring[h1]);
                return (uint64_t)xstr.xstring[h1].val2;
            default:
                h1 = (h1 + delta + 1) & xstr.mask;
                delta++;
        }
    }
    //return _DefaultXNodeBase._Empty_Dir_;
    return 0;
}
//#define debug_c_ysa
/*  
 * serial sort ysa 
 */
bool _createYSA(String<uint64_t> & hs, XString & xstr, uint64_t & indexEmptyDir, uint64_t thd_blocklimit)
{
//-k
    uint64_t k = _DefaultHs.getHeadPtr(hs[0]);
    uint64_t preX = _DefaultHs.getHeadX(hs[0]);
    //uint64_t k = 2 + _DefaultHs.getHeadPtr(hs[2]);
    //uint64_t preX = _DefaultHs.getHeadX(hs[2]);
    uint64_t ptr = k;
    uint64_t block_size = ptr;
    uint64_t countMove = 0, prek = 0;
    double time = sysTime();
    //remove duplicate hs header
    while(_DefaultHs.getHeadPtr(hs[k]))
    {
        ptr = _DefaultHs.getHeadPtr(hs[k]);
        if (preX != _DefaultHs.getHeadX(hs[k]))
        {
            
            hs[k - countMove] = hs[k];
            _DefaultHs.setHsHeadPtr(hs[prek], block_size);
            prek = k - countMove;
            block_size = ptr;
            preX = _DefaultHs.getHeadX(hs[k]);
        }
        else
        {
            countMove++;
            block_size += ptr - 1;
        }
        for (uint64_t j = k + 1; j < k + ptr; j++)
        {
            hs[j - countMove]= hs[j];       
        }
           
        k += ptr;
    }
    _DefaultHs.setHsHeadPtr(hs[prek], block_size);
    //hs[k - countMove] = 0;
    //hs[k - countMove + 1] = 0;
    _DefaultHs.setHsHead(hs[k - countMove], 0, 0);
    _DefaultHs.setHsHead(hs[k - countMove + 1], 0, 0);

    resize(hs, k + 2 - countMove);
    shrinkToFit(hs);
    indexEmptyDir = k - countMove;
//-k 
    k=0;
    preX = _DefaultHs.getHeadX(hs[0]);
    //k=2;
    //preX = _DefaultHs.getHeadX(hs[2]);

    ptr = 0;
    block_size = 0;

    std::cerr << "      preprocess sort y " << sysTime() - time << std::endl;
    time = sysTime();
    while(_DefaultHs.getHeadPtr(hs[k]))
    {
        ptr = _DefaultHs.getHeadPtr(hs[k]);
        _sort_YSA_Block(begin(hs) + k + 1, begin(hs) + k + ptr);
        k += ptr;
    }
    std::cerr << "      sort y " << sysTime() - time << std::endl;
//-k
    ptr = 0; k = 0;

    uint64_t count = 0; 
    
//    ptr = 0; k = 2;
    time = sysTime();
    while(_DefaultHs.getHeadPtr(hs[k]))
    {
        ptr = _DefaultHs.getHeadPtr(hs[k]);
        if (ptr < thd_blocklimit)
        {
            ++count;
        }
        else
        {   
            for (unsigned j = k + 1; j < k + ptr; j++)
            {
                if(_DefaultHs.getHsBodyY(hs[j] ^ hs[j - 1]))
                {
                    ++count;
                }
            }
            ++count;
        }
        k += ptr;
    }
    xstr._fullSize(count);
    ptr = 0; k = 0;
    std::cerr << "preprocess2 resize xstr " << sysTime() - time << std::endl;
    time = sysTime();
    while(_DefaultHs.getHeadPtr(hs[k]))
    {
        ptr = _DefaultHs.getHeadPtr(hs[k]);
        if (ptr < thd_blocklimit)
        {
            requestXNode_noCollision(xstr, _DefaultHs.getHeadX(hs[k]), 
                    k + 1, _DefaultXNodeBase.xHead, _DefaultXNodeBase.returnDir);
        }
        else
        {
            uint64_t xval = _DefaultHs.getHeadX(hs[k]);
            requestXNode_noCollision(xstr, xval, 
                    ~1, _DefaultXNodeBase.virtualHead, _DefaultXNodeBase.returnDir);
            for (unsigned j = k + 1; j < k + ptr; j++)
            {
                if(_DefaultHs.getHsBodyY(hs[j] ^ hs[j - 1]))
                {
                    requestXNode_noCollision(xstr, (xval + ((hs[j] & ((1ULL<<61) - (1ULL<<41))) >>1)), 
                            j, _DefaultXNodeBase.xHead, _DefaultXNodeBase.returnDir);
                }
                   
            }
        }
        k += ptr;
    }
    std::cerr << "RequestDir " << sysTime() - time << std::endl;
    return true;
}

/* parallel sort ysa
 * this function is to create index of minihash value (xval), yval is omitted![minindex]
 * create ysa within [@hs_str, hs_end) of @hs
 * WARN::be sure @hs[hs_str] and @hs[hs_end] are Hs::headNode type
 */
bool _createYSA(String<uint64_t> & hs, XString & xstr, uint64_t hs_str, uint64_t hs_end, uint64_t & indexEmptyDir, bool f_shrink_hs, bool f_ysa_sorted, unsigned threads, uint64_t thd_blocklimit)
{
    unused(f_ysa_sorted);
    std::cerr << "=>Index::SortYSA                                                  \r";
    double time = sysTime();
    uint64_t ptr  = _DefaultHs.getHeadPtr(hs[hs_str]);
    uint64_t preX = _DefaultHs.getHeadX  (hs[hs_str]);
    uint64_t prek = hs_str;
    uint64_t k     = hs_str + ptr;
    uint64_t block_size = ptr;
    uint64_t countMove = 0;

    while(k < hs_end && _DefaultHs.getHeadPtr(hs[k]))
    {
        ptr = _DefaultHs.getHeadPtr(hs[k]);
        if (preX != _DefaultHs.getHeadX(hs[k]))
        {
            hs[k - countMove] = hs[k];
            _DefaultHs.setHsHeadPtr(hs[prek], block_size);
            prek = k - countMove;
            block_size = ptr;
            preX = _DefaultHs.getHeadX(hs[k]);
        }
        else
        {
            countMove++;
            block_size += ptr - 1;
        }
        for (uint64_t j = k + 1; j < k + ptr; j++)
        {
            hs[j - countMove]= hs[j];       
        }
        k += ptr;
    }
    uint64_t hs_str_modified;
    uint64_t hs_end_modified;
    if (countMove > 2) //for most cases, countMove > 2
    {
        hs_str_modified = hs_str;
        hs_end_modified = k -countMove;
        _DefaultHs.setHsHeadPtr(hs[prek], block_size);
        _DefaultHs.setHsHead(hs[k - countMove], 0, 0); 
        _DefaultHs.setHsHead(hs[k - countMove + 1], 0, 0); 
        indexEmptyDir = k - countMove;
    }
    else //abort the last block, it nearly have no impact on large genomes
    {
        hs_str_modified = hs_str;
        hs_end_modified = prek;
        _DefaultHs.setHsHead(hs[prek], 0, 0);
        _DefaultHs.setHsHead(hs[prek + 1], 0, 0); 
        indexEmptyDir = prek; 
    }

    if (f_shrink_hs)
    {
        resize(hs, k + 2 - countMove);
        shrinkToFit(hs);
    }

    #pragma omp parallel
    {
        uint64_t ptr = 0;
        #pragma omp for
        for (uint64_t k = hs_str_modified; k < hs_end_modified; k++)
        { 
            if(_DefaultHs.isHead(hs[k]))
            {
                ptr = _DefaultHs.getHeadPtr(hs[k]);
                _sort_YSA_Block(begin(hs) + k + 1, begin(hs) + k + ptr);
            }   
        }
    }

    std::cerr << "  Index::SortYSA              Elapsed Time[s] " << sysTime() - time << std::endl;
    std::cerr << "=>Index::ResizeXstr                                            \r" ;
    k = hs_str_modified;
    uint64_t count = 0; 
    time = sysTime();
    float c1 = 0, c2 = 0;
    while(_DefaultHs.getHeadPtr(hs[k]) && k < hs_end_modified)
    {
        ptr = _DefaultHs.getHeadPtr(hs[k]);
        if (ptr < thd_blocklimit)
        {
            ++count;
            c1+=ptr;
        }
        else
        {   
            for (unsigned j = k + 1; j < k + ptr; j++)
            {
                if(_DefaultHs.getHsBodyY(hs[j] ^ hs[j - 1]))
                {
                    ++count;
                }
            }
            ++count;
            c2+=ptr;
        }
        k += ptr;
    }
    xstr._fullSize(count);
    std::cerr << "  Index::ResizeXstr           Elapsed Time[s] " << sysTime() - time << std::endl;
    std::cerr << "=>Index::RequestDir                                                  \r";
    time = sysTime();
    //std::cout << "lbks " << c1 / (c1 + c2) << "\n";

    #pragma omp parallel 
    {
        uint64_t ptr = 0;
        #pragma omp for 
        for (uint64_t i = hs_str_modified; i < hs_end_modified; i++)
        {

            if (_DefaultHs.isHead(hs[i]) && _DefaultHs.getHeadPtr(hs[i]))
            {
                ptr = _DefaultHs.getHeadPtr(hs[i]);
                if (ptr < thd_blocklimit)
                {
                    for (unsigned j = i + 1; j < i + ptr; j++)
                    {
                        _DefaultHs.setHsBodyY(hs[j], 0);
                    }
                    requestXNode_noCollision_Atomic(xstr, _DefaultHs.getHeadX(hs[i]), 
                           i + 1, _DefaultXNodeBase.xHead, _DefaultXNodeBase.returnDir);   
                }
                else
                {
                    uint64_t xval = _DefaultHs.getHeadX(hs[i]);
                    requestXNode_noCollision(xstr, xval, 
                        ~1, _DefaultXNodeBase.virtualHead, _DefaultXNodeBase.returnDir);
                    for (unsigned j = i + 1; j < i + ptr; j++)
                    {
                        if(_DefaultHs.getHsBodyY(hs[j] ^ hs[j - 1]))
                        {
                            requestXNode_noCollision(xstr, (xval + ((hs[j] & ((1ULL<<61) - (1ULL<<41))) >>1)), 
                                j, _DefaultXNodeBase.xHead, _DefaultXNodeBase.returnDir);
                        }
                    }
                }
            }
        }
    }

    std::cerr << "  Index::RequestDir           Elapsed Time[s] " << sysTime() - time << std::endl;
    (void) threads;
    return true;
}

bool _createYSA(String<uint64_t> & hs, XString & xstr, uint64_t & indexEmptyDir, unsigned threads, uint64_t thd_blocklimit)
{
    return _createYSA(hs, xstr, 0, _DefaultHs.getLength(hs), indexEmptyDir, true, false, threads, thd_blocklimit);
}

/*
 * free geonme sequence during creating, for raw map
 */
 bool _createQGramIndexDirSA_parallel(StringSet<String<Dna5> > & seq, 
        XString & xstr, String<uint64_t> & hs,  
        LShape & shape, 
        uint64_t & indexEmptyDir, 
        unsigned gstr,
        unsigned gend,
        unsigned & threads, 
        unsigned thd_step,
        uint64_t thd_blocklimit,
        bool efficient)    
{
    double time = sysTime();
    //_createHsArray2_MF(seq, hs, shape, threads);
    _createHsArray(seq, hs, shape, gstr, gend, threads, thd_step, efficient);
    _createYSA(hs, xstr, indexEmptyDir, threads, thd_blocklimit);
    std::cerr << "  End creating Index Time[s]:" << sysTime() - time << " \n";
    return true; 
}

bool createHIndex(StringSet<String<Dna5> > & seq, HIndex & index, unsigned gstr, unsigned gend, unsigned & threads, unsigned thd_step, uint64_t thd_blocklimit, bool efficient)
{
    return _createQGramIndexDirSA_parallel(seq, index.xstr, index.ysa, index.shape, index.emptyDir, gstr, gend, threads, thd_step, thd_blocklimit, efficient);
}

/*----------  DIndex and wrapper  ----------*/

int const typeDIx = 1;
int const typeHIx = 2;
int const typeSIx = 4;
int const typeMIx = 8;
int const typeFIx = 16;

unsigned dshape_len = 21; //!!WARN::only odd number, even is not allowed
DIndex::DIndex():
    rs( _DefaultCordBase.flag_strand),
    fs(~_DefaultCordBase.flag_strand),
    shape(dshape_len)
{}
DIndex::DIndex (unsigned len):
    shape(len)
{}
LShape & DIndex::getShape()
{
    return shape;
}
int DIndex::fullSize()
{
    return (1 << shape.weight << shape.weight) + 1;
}
String<int> & DIndex::getDir()
{
    return dir;
}
String<uint64_t> & DIndex::getHs()
{
    return hs;
}
uint64_t DIndex::val2Anchor(int64_t & i, uint64_t  & y, uint64_t & read_len, LShape & shape)
{
    if (get_cord_strand(hs[i]) ^ shape.strand)
    {
        uint64_t cordy = read_len - 1 - y;
        return (hs[i] - (cordy << 20) + cordy) | rs;
    }
    else 
    {
        return (hs[i] - (y << 20) + y) & fs;
    }    
}
void DIndex::clear()
{
    seqan::clear(dir);
    shrinkToFit(dir);
    seqan::clear(hs);
    shrinkToFit(hs);
}
int createDIndex_serial(StringSet<String<Dna5> > & seqs, 
                        DIndex & index, 
                        int64_t thd_min_step, 
                        int64_t thd_max_step,
                        unsigned gstr,
                        unsigned gend)
{
    unused(gstr);
    unused(gend);
    LShape & shape = index.getShape();
    String<int> & dir = index.getDir();
    String<uint64_t> & hs = index.getHs();
    resize (index.getDir(), index.fullSize(), 0);
    uint64_t preVal = ~0;
    uint64_t last_j = 0;
    for (uint64_t i = 0; i < length(seqs); i++)
    {
        int64_t count = 0;
        hashInit (shape, begin(seqs[i]));
        for (uint64_t j = 0; j < length(seqs[i]) - shape.span; j++)
        {
            hashNexth(shape, begin(seqs[i]) + j);
            if (++count > thd_min_step)
            {
                hashNextX(shape, begin(seqs[i]) + j);
                if (preVal != shape.XValue || j - last_j > (uint64_t)thd_max_step)
                {
                    ++dir[shape.XValue];
                    preVal = shape.XValue;
                    last_j = j;
                }
                count = 0;
            }
        }
    }
    int64_t sum = 0;
    for (uint64_t i = 0; i < length(dir); i++)
    {
        sum += dir[i];
        dir[i] = sum - dir[i];
    }
    last_j = 0;
    int64_t EmptyVal = create_cord(length(seqs),0,0,0); 
    //make sure genomeid >= length(seqs) and cord y be 0! y points to next empty.
    resize (hs, sum, EmptyVal);
    for (uint64_t i = 0; i < length(seqs); i++)
    {
        int64_t count = 0; 
        hashInit (shape, begin(seqs[i]));
        for (uint64_t j = 0; j < length(seqs[i]) - shape.span; j++)
        {
            hashNexth (shape, begin(seqs[i]) + j);
            if (++count > thd_min_step)
            {
                hashNextX (shape, begin(seqs[i]) + j);
                if (preVal != shape.XValue || j - last_j > (uint64_t)thd_max_step)
                {
                    int k = dir[shape.XValue];
                    k += get_cord_y (hs[k]);
                    hs[k] = create_cord(i, j, 0, shape.strand);
                    hs[dir[shape.XValue]] = shift_cord(hs[dir[shape.XValue]], 0, 1); //y++;
                    //std::cout << shape.XValue << " " << k <<"xxx3\n";

                    //hs[k]++;
                    preVal = shape.XValue;
                    last_j = j;

                } 
                count = 0;
            }
        }  
    }
    return 0;
}

uint64_t const DINDEXY_BITS2 = 5; //shape.YValue bits
uint64_t const DINDEXY_BITS1 = 20 - DINDEXY_BITS2; //block share pointer bits
uint64_t const DINDEXY_MASK1 = (1ULL << DINDEXY_BITS1) - 1;
uint64_t shape2DIndexCordy(LShape & shape)
{
    //The DIndex stores a part of shape.YValue only
    //get lower DINDEXY_BITS2 bits of yvalue and left shift DINDEXY_BITS1 bits as the cordy.
    return (shape.YValue & ((1ULL << (DINDEXY_BITS2)) - 1) << DINDEXY_BITS1);
}
uint64_t getDIndexCordy(uint64_t index_val) 
{
    return get_cord_y(index_val) & (~DINDEXY_MASK1);
}
uint64_t getDIndexBlockPointer_(int64_t index_val)
{
    return get_cord_y(index_val) & (DINDEXY_MASK1);
}

/**
 * create index on [@gstr, @gend)th genomes 
 * @hs[i] use cord struct. while it's yvalue is differnt from the cord
 * y in hs[i]:= a part of shape.yvalue[DINDEXY_BITS2 bits]|
                pointer to block.if first in the block: pointer to the last, 
                otherwise 0[DINDEXY_BITS1 bits]
 */
int createDIndex(StringSet<String<Dna5> > & seqs, 
                 DIndex & index, 
                 int64_t thd_min_step, 
                 int64_t thd_max_step,
                 int64_t thd_omit_block,
                 unsigned gstr,
                 unsigned gend,
                 unsigned threads)
{
    //std::cerr << std::fixed << std::setprecision(2);
    //the maximum block size is (1 << DINDEXY_BITS1) - 1
    serr.print_message("=>Index::Initiating ", 0, 2, std::cerr);
    double t = sysTime();
    LShape & t_shape = index.getShape();
    String<int> & dir = index.getDir();
    String<uint64_t> & hs = index.getHs();
    resize (dir, index.fullSize(), 0);
    unsigned genome_lens_sum = 0;
    float current_lens_ratio = 0;
    float finished_ratio = 0;
    for (unsigned i = gstr; i < gend; i++)
    {
        genome_lens_sum += length(seqs[i]);
    }
    for (unsigned i = gstr; i < gend; i++)
    {
        String<int64_t> t_blocks;
        current_lens_ratio = float(length(seqs[i])) / genome_lens_sum;
        for (unsigned j = 0; j < threads; j++)
        {
            appendValue(t_blocks, length(seqs[i]) / threads * j); 
        }
        appendValue (t_blocks, length(seqs[i]) - t_shape.span);
        #pragma omp parallel
        {
            unsigned t_id = omp_get_thread_num();
            int64_t t_str = t_blocks[t_id];
            int64_t t_end = t_blocks[t_id + 1];
            int64_t last_j = t_str - 1;
            int64_t count = 0;
            uint64_t preVal = ~0;
            LShape shape = t_shape;
            hashInit (shape, begin(seqs[i]) + t_str);
            //dout << "cdx1 " << t_str<< t_end <<"\n";
            int64_t t_percent_cerr = (t_end - t_str) * 2 / 100; //cerr percent every 2%
            int64_t j_count = 0;
            for (int64_t j = t_str; j < t_end; j++)
            {
                hashNexth(shape, begin(seqs[i]) + j);
                if (++count > thd_min_step)
                {
                    hashNextX(shape, begin(seqs[i]) + j);
                    if (preVal != shape.XValue || j - last_j > thd_max_step)
                    {
                        atomicInc(dir[shape.XValue]);
                        preVal = shape.XValue;
                        last_j = j;
                    }
                    count = 0;
                }
                if (t_id == 0 && j_count == t_percent_cerr)
                {
                    j_count = 0; 
                    finished_ratio += float(t_percent_cerr) / (t_end - t_str) * 100 * current_lens_ratio;

                    serr.print_message("=>Index::Initiating [", 0, 0, std::cerr);
                    serr.print_message(int(finished_ratio), 0, 0, std::cerr);
                    serr.print_message("%]", 0, 2, std::cerr);
                }
                ++j_count;
            }
        }
    }
    //double t4 = sysTime();
    int64_t sum = 0;
    for (uint64_t i = 0; i < length(dir); i++)
    {
        if (dir[i] > thd_omit_block)
        {
            dir[i] = 0;
        }
        /*
        if (dir[i] != 0)
        {
            dout << "cdx" << dir[i] << "\n";
        }
        */
        sum += dir[i];
        dir[i] = sum - dir[i];
    }
    int64_t EmptyVal = create_cord(length(seqs),0,0,0); 
    //make sure genomeid >= length(seqs) and cord y be 0! y points to next empty.
    resize (hs, sum, EmptyVal);
    serr.print_message("--Index::Initiate[100%]   ", 0, 1, std::cerr);
    serr.print_message("=>Index::Hashing", 0, 2, std::cerr);

    current_lens_ratio = 0;
    finished_ratio = 0;
    for (uint64_t i = 0; i < length(seqs); i++)
    {
        String<int64_t> t_blocks;
        current_lens_ratio = float(length(seqs[i])) / genome_lens_sum;
        for (unsigned j = 0; j < threads; j++)
        {
            appendValue(t_blocks, length(seqs[i]) / threads * j); 
        }
        appendValue (t_blocks, length(seqs[i]) - t_shape.span); 
        //dout << "cx22"<<t_blocks << "\n";
        #pragma omp parallel
        {   
            unsigned t_id = omp_get_thread_num();
            int64_t t_str = t_blocks[t_id];
            int64_t t_end = t_blocks[t_id + 1];
            int64_t last_j = t_str - 1;
            int64_t count = 0;
            uint64_t preVal = ~0;
            LShape shape = t_shape;
            hashInit (shape, begin(seqs[i]) + t_str);
            //dout << "strd " << t_str << t_end << "\n";
            int64_t t_percent_cerr = (t_end - t_str) * 1 / 50;
            int64_t j_count = 0;
            for (int64_t j = t_str; j < t_end; j++)
            {
                hashNexth (shape, begin(seqs[i]) + j);
                if (++count > thd_min_step)
                {
                    hashNextX (shape, begin(seqs[i]) + j);
                    if (preVal != shape.XValue || j - last_j > thd_max_step)
                    {
                        if (dir[shape.XValue + 1] - dir[shape.XValue])
                        {
                            int64_t slot_str = dir[shape.XValue];
                            int64_t k = slot_str + getDIndexBlockPointer_(atomic_inc_cord_y(hs[slot_str])) - 1;
                            //int64_t new_cord = create_cord(i, j, 0, shape.strand); 
                            uint64_t new_cord = create_cord(i, j + const_anchor_zero, 0, shape.strand); //be sure lower bits of new_cord_y for share pointer == 0 
                            if (k == slot_str) 
                            {   //atomic creation the first cord which cotains shared pointer
                                new_cord -= EmptyVal;
                                atomicAdd(hs[k], new_cord); //original hs[k] = EmptyVal + pointer
                            }
                            else
                            {
                                hs[k] = new_cord;
                            }
                            preVal = shape.XValue;
                            last_j = j;
                        }
                    } 
                    count = 0;
                }
                if (t_id == 0 && j_count == t_percent_cerr)
                {
                    j_count = 0; 
                    finished_ratio += float(t_percent_cerr) / (t_end - t_str) * 100 * current_lens_ratio;
                    serr.print_message("=>Index::Hashing [", 0, 0, std::cerr);
                    serr.print_message(int(finished_ratio), 0, 0, std::cerr);
                    serr.print_message("%]", 0, 2, std::cerr);
                }
                ++j_count;
            }  
        }
    }
    //dout << "size" << float(length(dir)) * 8 / 1024/1024/1024 << float(length(hs)) /128/1024/1024 << "\n";
    //std::cout << "createDIndex " << sysTime() - t << " " << sysTime() - t2 << "\n";
    serr.print_message("Index::Hash    [100%]              ", 2, 1, std::cerr);
    serr.print_message("End creating index ", 2, 0, std::cerr);
    serr.print_message("Elapsed time[s] ", 2, 0, std::cerr);
    serr.print_message(sysTime() - t, 2, 1, std::cerr);
    return 0;
}
/*----------  SIndex  ----------*/
SIndex::SIndex():
    shape(21),
    cas_key(0)
{}
SIndex::SIndex (unsigned len):
    shape(len),
    cas_key(0)
{}
LShape & SIndex::getShape()
{
    return shape;
}
int SIndex::fullSize()
{
    return (1 << shape.weight << shape.weight) + 1;
}
std::vector<std::vector<int64_t> > & SIndex::getHs()
{
    return hs;
}
void SIndex::clear()
{
    //seqan::clear(hs);
    //shrinkToFit(hs);
    hs.clear();
    hs.shrink_to_fit();
}
int64_t SIndex::getVal(int64_t key, int64_t pos)
{
    return hs[key][pos];
}
int SIndex::printStatistics()
{
    if (empty(hs))
    {
        return 0;
    }
    float mean = 0 ;
    String<float> lens;
    resize(lens, hs.size(), 0);
    for (unsigned i = 0; i < hs.size(); i++)
    {
        lens[i] += hs[i].size();
    }
    std::sort(begin(lens), end(lens));
    int c0 = 0; //count of 0 elments
    for (unsigned i = 0; i < length(lens); i++) 
    {
        if (lens[i] != 0)
        {
            c0 = i - 1;
            break;
        }
    }
    //count none 0 elments only
    erase (lens, 0, c0);
    int count_large_block = 0;
    int count_large_block_key = 0;
    for (unsigned i = 0; i < length(lens); i++) 
    {
        if (lens[i] > 1024)
        {
            count_large_block++;
            count_large_block_key += lens[i];
        }
        mean += lens[i];
    }
    mean /= (length(lens));

    float var = 0;
    for (unsigned i = 0; i < length(lens); i++) 
    {
        var += (lens[i] - mean) * (lens[i] - mean);
    }
    var /= (length(lens) - 1);
    dout << "sindx_statistics" << 
             mean << 
             var  << 
             length(lens) << 
             "m" <<
             lens[0] << 
             lens[unsigned(length(lens) * 0.25)] << 
             lens[unsigned(length(lens) * 0.5)]
        << lens[unsigned(length(lens) * 0.75)]
        <<lens[unsigned(length(lens) * 0.95)] << back(lens) << mean * length(lens) << count_large_block 
        << count_large_block_key / float(mean * length(lens))
        << "\n";
    return 0;
}
int64_t SIndex::queryHsStr(int64_t xval)
{
    (void) xval;
    return 0;
}
int64_t SIndex::queryHsEnd(int64_t xval)
{
    return length(hs[xval]);
}
int64_t SIndex::queryHsBlockLen(int64_t xval)
{
    return queryHsEnd(xval);
}
/*
 * N
 */
/*
int _createSIndexHs11Thread(String<Dna5> & seq, 
                 SIndex & index,
                 uint64_t genome_n, 
                 int64_t thd_min_step, 
                 int64_t thd_max_step,
                 int64_t thd_omit_block,
                 unsigned seq_str,
                 unsigned seq_end,
                 unsigned thread_id)
{
    LShape & t_shape = index.getShape();
    StringSet<String<int64_t> > & hs = index.getHs()[thread_id];
    std::cerr << "csh1 " << thread_id << " " << seq_str << "\n";
    if (empty(hs))
    {
        resize (hs, index.fullSize()); 
    }
    int64_t count = 0;
    LShape shape = t_shape;
    hashInit (shape, begin(seq[genome_n]) + seq_str);
    for (int64_t i = seq_str; i < seq_end; i++)
    {
        hashNexth(shape, begin(seq) + i);
        if (++count > thd_min_step)
        {
            hashNextX(shape, begin(seq) + i);
            //if (preVal != shape.XValue || j - last_j > thd_max_step)
            //{
                //atomicInc(dir[shape.XValue]);
                //if (index.queryHashBlockSize(shape.XValue) < thd_omit_block)
                //{
                    appendValue(hs[shape.XValue], create_cord(genome_n, i, 0, shape.strand));
                //}
            //}
            count = 0;
        }
    }
    return 0;
}
int _createSIndexHs12Thread(String<Dna5> & seq, 
                            SIndex & index,
                            uint64_t genome_n, 
                            int64_t thd_min_step, 
                            int64_t thd_max_step,
                            int64_t thd_omit_block,
                            unsigned seq_str,
                            unsigned seq_end,
                            unsigned thread_id)
{
    LShape & t_shape = index.getShape();
    StringSet<String<int64_t> > & hs = index.getHs()[thread_id];
    if (empty(hs))
    {
        std::cout << "csh21 " << thread_id << " " << seq_str << " " << index.fullSize() << "\n";
        resize (hs, index.fullSize()); 
    }
    int64_t thd_rehash = thd_omit_block / 2;
    int64_t last_i = 0;
    LShape shape = t_shape;
        std::cout << "csh22 " << thread_id << " " << seq_str << " " << index.fullSize() << "\n";
    hashInit (shape, begin(seq[genome_n]) + seq_str);
    for (int64_t i = seq_str; i < seq_end; i++)
    {
        hashNexth(shape, begin(seq) + i);
        if (i - last_i > thd_min_step)
        {
            hashNextX(shape, begin(seq) + i);
            if (index.queryHashBlockSize(shape.XValue) < thd_rehash || i - last_i > thd_max_step)
            {
                appendValue(hs[shape.XValue], create_cord(genome_n, i, 0, shape.strand));
                last_i = i; 
            }
        }
    }
    return 0;
}
*/

void _initSIndexHs(SIndex & index, unsigned threads)
{
    //resize (index.getHs(), 1);
    //resize (index.getHs()[0], index.fullSize()); //hs[0] is always non empty no matter how many threads
    //resize (index.cas_keys, index.fullSize(), 0);
    index.getHs().resize(index.fullSize()); //hs[0] is always non empty no matter how many threads
    resize (index.cas_keys, index.fullSize(), 0);
    (void) threads;
}
/**
 * All threads access the same dir which is a StringSet
 * More memory efficient than indepent StringSet for each thread
 * Thread safe not tested 
int _createSIndexHsThreadUnit0(String<Dna5> & seq, 
                 SIndex & index,
                 uint64_t genome_n, 
                 int64_t thd_min_step, 
                 int64_t thd_max_step,
                 int64_t thd_omit_block,
                 unsigned seq_str,
                 unsigned seq_end,
                 unsigned thread_id)
{
    LShape & t_shape = index.getShape();
    std::vector<std::vector<int64_t> > & hs = index.getHs();
    
    int64_t thd_rehash = thd_omit_block / 2;
    int64_t last_i = seq_str;
    int64_t pre_key = -1;
    LShape shape = t_shape;
    String<std::tuple<int64_t, int64_t, int64_t> > key_val_lens; 
    hashInit (shape, begin(seq[genome_n]) + seq_str);
    for (int64_t i = seq_str; i < seq_end; i++)
    {
        hashNexth(shape, begin(seq) + i);
        if (i - last_i > thd_min_step)
        {
            hashNextX(shape, begin(seq) + i);
            appendValue(key_val_lens, std::make_tuple(
                        int64_t(shape.XValue), 
                        int64_t(create_cord(genome_n, i, 0, shape.strand)), 
                        int64_t(index.queryHsBlockLen(shape.XValue))));
            if (i - last_i >= thd_max_step)
            {
                unsigned min_j = length(key_val_lens) - 1;
                int64_t new_key = std::get<0>(back(key_val_lens));
                int64_t new_val = std::get<1>(back(key_val_lens));
                int64_t min_len = std::get<2>(back(key_val_lens));
                if (std::get<2>(back(key_val_lens)) > thd_rehash)
                {
                    for (unsigned j = 0; j < length(key_val_lens); j++)
                    {
                        if (std::get<2>(key_val_lens[j]) < min_len)
                        {
                            min_len = std::get<2>(key_val_lens[j]);
                            min_j = j;
                        }
                    }
                    new_key = std::get<0>(key_val_lens[min_j]);
                    new_val = std::get<1>(key_val_lens[min_j]);
                }
                while (1)
                {
                    //dout << "sindx2" << i << shape.XValue << thread_id << "\n";
                    if (!atomicCas(index.cas_keys[new_key], false, true))
                    {
                        hs[new_key].push_back(new_val);
                        atomicCas(index.cas_keys[new_key], true, false);
                        pre_key = new_key;
                        last_i = get_cord_x(new_val);
                        break;
                    }
                }
                erase(key_val_lens, 0, std::min(unsigned(min_j + 1 + thd_min_step), unsigned(length(key_val_lens))));
            }
        }
    }
    return 0;
}
*/
/**
 * All threads access the same dir which is a StringSet
 * More memory efficient than indepent StringSet for each thread
 * Thread safe not tested 
 * @finished_ratio can be modified by one thread only
 */
int _createSIndexHsThreadUnit(String<Dna5> & seq, 
                 SIndex & index,
                 uint64_t genome_n, 
                 int64_t thd_min_step, 
                 int64_t thd_max_step,
                 int64_t thd_omit_block,
                 unsigned seq_str,
                 unsigned seq_end,
                 unsigned thread_id,
                 unsigned genome_lens_sum,
                 float & finished_ratio)
{
    LShape & t_shape = index.getShape();
    std::vector<std::vector<int64_t> > & hs = index.getHs();
    
    int64_t thd_rehash = thd_omit_block / 2;
    int64_t last_i = seq_str;
    LShape shape = t_shape;
    //Iterator<String<Dna5> > it = begin(seq[genome_n]) + seq_str;
    //hashInit (shape, it);
    hashInit(shape, begin(seq) + seq_str);
    int64_t preVal = -1;
    unsigned j_count = 0;
    int64_t t_percent_cerr = (seq_end - seq_str) * 2 / 100; //cerr percent every 2%
    float current_lens_ratio = 0;
    for (int64_t i = seq_str; i < seq_end; i++)
    {
        current_lens_ratio = float(length(seq)) / genome_lens_sum;
        hashNexth(shape, begin(seq) + i);
        //hashNexth_hpc(shape, it, begin(seq) + seq_end);
        if (i - last_i > thd_min_step)
        {
            hashNextX(shape, begin(seq) + i);
            //hashNextX(shape, it);
            if ((preVal != int64_t(shape.XValue)) || i - last_i > thd_max_step)
            //if ((preVal != shape.XValue && index.queryHsBlockLen(shape.XValue) < thd_rehash) || i - last_i > thd_max_step)
            {
                while (1)
                {
                    if (!atomicCas(index.cas_keys[shape.XValue], false, true))
                    {
                        hs[shape.XValue].push_back(create_cord(genome_n, i, 0, shape.strand));
                        atomicCas(index.cas_keys[shape.XValue], true, false);
                        last_i = i; 
                        preVal = shape.XValue;
                        break;
                    }
                }
            }
        }
        if (thread_id == 0 && j_count == t_percent_cerr)
        {
            j_count = 0; 
            finished_ratio += float(t_percent_cerr) / (seq_end - seq_str) * 100 * current_lens_ratio;
            serr.print_message("=>Index::Hashing [", 0, 0, std::cerr);
            serr.print_message(int(finished_ratio), 0, 0, std::cerr);
            serr.print_message("%]", 0, 2, std::cerr);
        }
        ++j_count;
    }
    (void)thd_rehash;
    return 0;
}
int createSIndex(StringSet<String<Dna5> > & seqs, 
                 SIndex & index, 
                 int64_t thd_min_step, 
                 int64_t thd_max_step,
                 int64_t thd_omit_block,
                 unsigned gstr,
                 unsigned gend,
                 unsigned threads)
{
    double t0 = sysTime();
    serr.print_message("=>SIndex::Initiating ", 0, 2, std::cerr);
    LShape & t_shape = index.getShape();
    _initSIndexHs(index, threads);
    serr.print_message("--SIndex::Initiate[100%]   ", 0, 1, std::cerr);
    serr.print_message("=>SIndex::Hashing", 0, 2, std::cerr);
    double t1 = sysTime();
    unsigned genome_lens_sum = 0;
    float finished_ratio = 0;
    for (unsigned i = gstr; i < gend; i++)
    {
        genome_lens_sum += length(seqs[i]);
    }
    for (int64_t i = gstr; i < gend; i++)
    {
        String<int64_t> t_blocks;
        for (unsigned j = 0; j < threads; j++)
        {
            appendValue(t_blocks, length(seqs[i]) / threads * j); 
        }
        appendValue (t_blocks, length(seqs[i]) - t_shape.span);
        #pragma omp parallel
        {
            unsigned t_id = omp_get_thread_num();
            int64_t t_str = t_blocks[t_id];
            int64_t t_end = t_blocks[t_id + 1];
            _createSIndexHsThreadUnit (seqs[i], index, i, 
                thd_min_step, thd_max_step, thd_omit_block, t_str, t_end, t_id, genome_lens_sum, finished_ratio);
        }
    }
    t1 = sysTime() - t1;
    //index.printStatistics();
    #pragma omp for 
    for (uint64_t i = 0; i < length(index.getHs()); i++)
    {
        if (length(index.getHs()[i]) > (uint64_t) thd_omit_block)
        {
            clear(index.getHs()[i]);
            index.getHs()[i].shrink_to_fit();
        }
    }
    clear(index.cas_keys);
    shrinkToFit(index.cas_keys);

    serr.print_message("SIndex::Hash[100%]                    ", 2, 1, std::cerr);
    serr.print_message("End creating index ", 2, 0, std::cerr);
    serr.print_message("Elapsed time[s] ", 2, 0, std::cerr);
    serr.print_message(sysTime() - t0, 2, 1, std::cerr);
    return 0;
}
CreateSIndexParms::CreateSIndexParms()
{
    thd_shape_len = 24;
    thd_min_step = 8; 
    thd_max_step = 10;
    thd_omit_block = 1024;
}

/*----------  MDindex  ----------*/
int _createDIndexFromHs(String<uint64_t> & hs, String<uint64_t> & hs_str_end, DIndex & index, int64_t thd_omit_block, unsigned threads)
{
    String<int> & dir = index.getDir();
    String<uint64_t> & d_hs = index.getHs();
    resize (dir, index.fullSize(), 0);
    //double t2 = sysTime();
    //dout << "idx2" << length(dir) << t_shape.weight << thd_min_step << thd_max_step << thd_omit_block<< "\n";
    //dout << threads << "threads\n"; 
    String<int64_t> t_blocks;
    //assign task region within the @hs to each threads
    uint64_t t_str = 0;
    for (unsigned i = 0; i < threads; i++)
    {
        while (t_str < _DefaultHs.getLength(hs) && !_DefaultHs.isHead(hs[t_str]))
        {
            t_str++;
        }
        appendValue(t_blocks, t_str); 
        t_str += _DefaultHs.getLength(hs) / threads;
    }
    appendValue (t_blocks, _DefaultHs.getLength(hs));
    //init @dir
    #pragma omp parallel
    {
        unsigned t_id = omp_get_thread_num();
        int64_t t_str = t_blocks[t_id];
        int64_t t_end = t_blocks[t_id + 1];
        uint64_t xval = 0;
        uint64_t pre_g_id = (t_str == 0) ? ~0 : _getSA_i1(_DefaultHs.getHsBodyS(hs[t_str - 1]));
        for (int64_t i = t_str; i < t_end; i++)
        {
            if (_DefaultHs.isHead(hs[i]))
            {
                xval = _DefaultHs.getHeadX(hs[i]);
                uint64_t current_g_id = _getSA_i1(_DefaultHs.getHsBodyS(hs[i + 1]));
                if (current_g_id != pre_g_id)
                {
                    //dout << "ci" << current_g_id << i << "\n";
                    hs_str_end[current_g_id] = i;
                    pre_g_id = current_g_id; 
                }
            }
            else{
                atomicInc(dir[xval]);
            }
        }
    }
    if (!empty(hs_str_end))
    {
        back(hs_str_end) = _DefaultHs.getLength(hs);
    }
    int64_t sum = 0;
    for (uint64_t i = 0; i < length(dir); i++)
    {
        if (dir[i] > thd_omit_block)
        {
            dir[i] = 0;
        }
        sum += dir[i];
        dir[i] = sum - dir[i];
    }
    //dout << "dt" << sum << back(dir) << "\n";
    int64_t max_seq_num = (1LL << 20) - 1; //ids in cords occupies 20 bit
    int64_t EmptyVal = create_cord(max_seq_num, 0, 0, 0); 
    //make sure genomeid >= length(seqs) and cord y be 0! y points to next empty.
    resize (d_hs, sum, EmptyVal);
    serr.print_message("--Index::inite   ", 0, 1, std::cerr);
    serr.print_message("=>Index::creating dindex", 0, 2, std::cerr);

    #pragma omp parallel
    {   
        unsigned t_id = omp_get_thread_num();
        int64_t t_str = t_blocks[t_id];
        int64_t t_end = t_blocks[t_id + 1];
        uint64_t xval = 0;
        for (int64_t j = t_str; j < t_end; j++)
        {
            if (_DefaultHs.isHead(hs[j]))
            {
                xval = _DefaultHs.getHeadX(hs[j]);
            }
            else if (dir[xval + 1] - dir[xval]) //too large blocks > thd_omit_block should (has been) be ommited in the @dir
            {
                int64_t slot_str = dir[xval];
                int64_t k = slot_str + get_cord_y(atomic_inc_cord_y(d_hs[slot_str])) - 1;
                int64_t new_cord = create_cord(_getSA_i1(_DefaultHs.getHsBodyS(hs[j])), 
                                               _getSA_i2(_DefaultHs.getHsBodyS(hs[j])), 0, 
                                               _DefaultHs.getHsBodyStrand(hs[j])); //be sure new_cord_y == 0 
                if (k == slot_str) 
                {   //atomic creating the first cord which cotains shared pointer
                    new_cord -= EmptyVal;
                    atomicAdd(d_hs[k], new_cord); //original hs[k] = EmptyVal + pointer
                }
                else
                {
                    d_hs[k] = new_cord;
                }
            }
            else {/*NONE*/}
        }
    }

    return 0;
}

int createMDIndex(StringSet<String<Dna5> > & seqs,
                  IndexDynamic & index, 
                  unsigned gstr, unsigned gend, 
                  int64_t thd_omit_block, unsigned threads, unsigned thd_step)
{
    unused(gstr);
    unused(gend);
    unused(thd_step);
    //t = sysTime();
    double t = sysTime();
    serr.print_message("=>Index::creating MD index          ", 0, 2, std::cerr);
    resize (index.hindex.ysa_str_end, length(seqs) + 1, 0);
    _createDIndexFromHs(index.hindex.ysa, index.hindex.ysa_str_end, index.dindex, thd_omit_block, threads);
    serr.print_message("Index::hash        ", 2, 1, std::cerr);
    serr.print_message("End creating index ", 2, 0, std::cerr);
    serr.print_message(sysTime() - t, 2, 1, std::cerr);
    return 0;
}

/*----------  MHindx  ----------*/

int _createHIndexFromHs(String<uint64_t> & hs, 
                        XString & xstr,
                        LShape & shape,
                        uint64_t & indexEmptyDir,
                        uint64_t g_hs_str, //[@g_hs_str, @g_hs_end) in @hs that belongs to the genome.
                        uint64_t g_hs_end,
                        bool     f_ysa_sorted, //flag indicating if ysa has already been sorted
                        uint64_t thd_blocklimit,
                        unsigned threads)
{
    //dout << "chssort" << g_hs_str << length(hs) << g_hs_end << "\n";
    _hsSort(begin(hs) + g_hs_str, begin(hs) + g_hs_end, shape.weight, threads);
    _createYSA(hs, xstr, g_hs_str, g_hs_end, indexEmptyDir, false, f_ysa_sorted, threads, thd_blocklimit);
    return 0;
}

int createMHIndex(IndexDynamic & index, uint64_t g_str, uint64_t g_end, uint64_t thd_blocklimit, unsigned threads)
{
    HIndex & hindex = index.hindex;
    bool f_ysa_sorted = hindex.insertYsaSortedRecord(g_str, g_end);
    //dout << "css1" << g_str << g_end << f_ysa_sorted << "\n";
    _createHIndexFromHs(hindex.ysa, hindex.xstr, hindex.getShape(), hindex.emptyDir,
                        hindex.ysa_str_end[g_str], hindex.ysa_str_end[g_end], f_ysa_sorted, thd_blocklimit, threads);
    return 0;
}

/*----------  DynamicIndex   ----------*/

int64_t queryHsStr(DIndex & index, uint64_t & xval)
{
    return index.getDir()[xval];
}
int64_t queryHsEnd(DIndex & index, uint64_t & xval)
{
    return index.getDir()[xval + 1];
}
void queryHsStrEnd(DIndex & index, uint64_t & xval, std::pair<int64_t, int64_t> & str_end)
{
    str_end.first  = index.getDir()[xval];
    str_end.second = index.getDir()[xval + 1];
}
int IndexDynamic::isHIndex()
{
    return typeIx & typeHIx;
}
int IndexDynamic::isSIndex()
{
    return typeIx & typeSIx;
}
int IndexDynamic::isDIndex()
{
    return typeIx & typeDIx;
}
int IndexDynamic::isMIndex()
{
    return typeIx & typeMIx;
}
void IndexDynamic::setHIndex()
{
    typeIx &= ~typeDIx;
    typeIx &= ~typeSIx;
    typeIx |=  typeHIx;
}
void IndexDynamic::setDIndex()
{
    typeIx &= ~typeHIx;
    typeIx &= ~typeSIx;
    typeIx |=  typeDIx;
}
void IndexDynamic::setSIndex()
{
    typeIx &= ~typeHIx;
    typeIx &= ~typeDIx;
    typeIx |=  typeSIx;
}
void IndexDynamic::setMIndex()
{
    typeIx &= ~typeFIx;
    typeIx |=  typeMIx;
}
void IndexDynamic::setMDIndex()
{
    setMIndex();
    setDIndex();
}
void IndexDynamic::setMHIndex()
{
    setMIndex();
    setHIndex();
}
void IndexDynamic::setFIndex()
{
    typeIx &= ~typeMIx;
    typeIx |=  typeFIx;
}
void IndexDynamic::setFDIndex()
{
    setFIndex();
    setDIndex();
}
void IndexDynamic::setFHIndex()
{
    setFIndex();
    setHIndex();
}
void IndexDynamic::setIndexType(int ix_type)
{
    if (1 << (ix_type - 1)== typeHIx)     {setHIndex();}
    else if (1 << (ix_type - 1) == typeDIx){setDIndex();}
    else if (1 << (ix_type - 1) == typeMIx){setMIndex();}
    else if (1 << (ix_type - 1) == typeFIx){setFIndex();}
    else if (1 << (ix_type - 1) == typeSIx){setSIndex();}
    else if (1 << (ix_type - 1) == typeHIx + typeMIx){setMHIndex();}
    else if (1 << (ix_type - 1) == typeDIx + typeMIx){setMDIndex();}
    else if (1 << (ix_type - 1) == typeHIx + typeFIx){setFHIndex();}
    else if (1 << (ix_type - 1) == typeDIx + typeFIx){setFDIndex();}
    else {setHIndex();}
}

void IndexDynamic::clearIndex()
{
    if (isHIndex())
    {
       hindex.clear();
    }
    else if (isDIndex())
    {
        dindex.clear();
    }
}
IndexDynamic::IndexDynamic(StringSet<String<Dna5> > & seqs) : hindex(), dindex(), typeIx(typeHIx)
{unused(seqs);}

int createMHsArray(StringSet<String<Dna5> > & seqs, IndexDynamic & index, uint64_t gstr, uint64_t gend, unsigned threads, uint64_t thd_step, bool f_recreate)
{
    if (empty(index.hindex.ysa) || f_recreate) 
    {
        double t = sysTime();
        serr.print_message("=>Index::creating hash array          ", 0, 2, std::cerr);
        __createHsArray(seqs, index.hindex.ysa, index.hindex.getShape(),  gstr, gend, threads, thd_step);
        serr.print_message("End creating hash array              ", 2, 0, std::cerr);
        serr.print_message(sysTime() - t, 2, 1, std::cerr);   
    }
    return 0;
}

bool createIndexDynamic(StringSet<String<Dna5> > & seqs, IndexDynamic & index, unsigned gstr, unsigned gend, unsigned threads, bool efficient)
{
    //dout << "idx" << index.typeIx << "\n";
    if (index.isMIndex())    
    {
        unsigned thd_shape_len = 23;
        int64_t thd_step = 10;
        index.dindex.getShape().init_shape_parm(thd_shape_len);
        index.hindex.getShape().init_shape_parm(thd_shape_len);
        if (index.isDIndex())
        {
            int64_t thd_omit_block = 50; 
            createMHsArray(seqs, index, gstr, gend, threads, thd_step, true);
            createMDIndex (seqs, index, gstr, gend, thd_omit_block, threads, thd_step);
        }
        else if (index.isHIndex())
        {
            //dout << "idx2" << index.typeIx << "\n";
            uint64_t thd_blocklimit = 32;
            createMHsArray(seqs, index, gstr, gend, threads, thd_step, false);
            createMHIndex (index, gstr, gend, thd_blocklimit, threads);
        }
        else{/*NONE*/}
    }
    else
    {
        if (index.isSIndex())
        //if (1)
        {
            unsigned thd_shape_len = 21;

            int64_t thd_min_step = 8 ;
            int64_t thd_max_step = 10;//thd_shape_len + 5;
            int64_t thd_omit_block = 200;//1024000; 

            index.sindex.getShape().init_shape_parm(thd_shape_len);
            createSIndex(seqs, index.sindex, thd_min_step, thd_max_step, thd_omit_block,
                                gstr, gend, threads);
            //test
            /*
            LShape shape = index.sindex.getShape();
            String<Dna5> & seq = seqs[0];
            unsigned count = 0;
            unsigned last_i = 0;
            std::pair<int, int> it;
            
            hashInit(shape, begin(seq));
            for (unsigned i = 0; i < length(seq); i++)
            {
                hashNexth(shape, begin(seq) + i);
                if (i - last_i > thd_min_step)
                {
                    int f_found = 0;
                    hashNextX(shape, begin(seq) + i); 
                    for (unsigned j = index.sindex.queryHsStr(shape.XValue); j < index.sindex.queryHsEnd(shape.XValue); j++)
                    {
                        if (i == unsigned(get_cord_x(index.sindex.getVal(shape.XValue, j))))
                        {
                            f_found = 1;
                            break;
                        }
                    }
                    std::cout << "sinx2 " << f_found << "\n";
                    last_i = i;
                }
            }
            */
            //index.sindex.printStatistics();
            //index.sindex.clear();
            //return createDIndex_serial(seqs, index.dindex, 4, 10);
        }
        else if (index.isDIndex())
        {
            int64_t thd_min_step = 8;
            int64_t thd_max_step = 10; 
            int64_t thd_omit_block = 400; //std::min(1024, (1 << DINDEXY_BITS1) - 1);  
            unsigned thd_shape_len = 21;
            index.dindex.getShape().init_shape_parm(thd_shape_len);
            return createDIndex(seqs, index.dindex, thd_min_step, thd_max_step, thd_omit_block,
                                gstr, gend, threads);
            //return createDIndex_serial(seqs, index.dindex, 4, 10);
        }
        else if (index.isHIndex())
        {
            //chr1: step=5, shape_len=25 | 0.077% blocks > 32 (thd_block_limit), 5.7% > kmers in these blocks
            //chr1: 5,  21 | 0.70% > 32 |0.22% > 64 
            //chr1: 10, 21 | 0.49% > 32 |0.16% > 64 4.8% kmers in...
            //grch37 1-22,x,y : 10, 25, 32| 1.3% > 32, 27.02% kmers in|0.6% > 64, 23.6%  
            unsigned thd_step = 8;
            unsigned thd_shape_len = 17; //WARN only odd number is allowed due to double strand hash
            uint64_t thd_blocklimit = 1024;
            float alpha = 1.6;
            index.hindex.shape.init_shape_parm(thd_shape_len / 2 * 2 + 1);
            index.hindex.alpha = alpha;
            return createHIndex(seqs, index.hindex, 
                                gstr, gend, 
                                threads, thd_step, thd_blocklimit, efficient);
        }
        else {/*NONE*/}
    }
    return false;
}

