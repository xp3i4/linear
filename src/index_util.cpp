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
This is the simplified version of optimized 25-mer index 
deriverd from  the original genric index.
This version index is specifically optimized for the use 
of Approximate mapping. 
Interface is not fully wrappered due to concerns regarding
performance issue. 
Sufficient benchmarks are required before wrapping.  
Please do not use it for other purpose!
*/
//Structt Hs: String<uint64_t>
//Types of node in Hs: 1.head node and 2.body node
//Head: Headflag[1] = 0|sortFlag[1]|N/A[2]|Pointer[20]| xvalue[40]
//Body: bodyflag[1] = 1|N/A[2]|yvalue[20] |typeCode[1]|sa[40]
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
 bool Hs::isHsBodyReverseStrand(uint64_t & val)
{
    return val & (_DefaultHsBase.bodyCodeFlag);
}
 void Hs::setHsBodyY(uint64_t & val, uint64_t y, uint64_t const & bit, uint64_t const & mask)
{
    val = (val & (~(mask << bit))) | (y << bit);
}

Hs _DefaultHs;
const unsigned index_shape_len = 25; 
const float def_alpha = 1.6;
HIndex::HIndex():
    alpha(def_alpha), 
    shape(index_shape_len)
    {}
HIndex::HIndex(unsigned shape_len, float val):
    alpha(val), 
    shape(shape_len)
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
 */
 bool _hsSortX_1(Iterator<String<uint64_t> >::Type const & begin, 
                       Iterator<String<uint64_t> >::Type const & end, 
                       unsigned const & xValBitLen, 
                       unsigned threads)
{
    const unsigned lowerBound = 20;
    const unsigned upperBound = 42;
    if (xValBitLen < lowerBound || xValBitLen > upperBound)
    {
        std::cerr << "[Error]: _dirSortX " << xValBitLen << "\n";
       // return false;
    }
    
    //unsigned const bit[18] = {9,4,9,4,9,4,8,5,8,5,8,5,8,5,7,6,7,6}; //xValueBitLen 34 - 42;
    unsigned const bit[upperBound - lowerBound + 1] 
        = {10,2,11,2,12,2,7,4,7,4,10,3,9,4,9,4,8,5,8,5,7,6}; //xValueBitLen 34 - 42;
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
        std::cerr << "=>Index::sortHash             " << std::setprecision(1) << std::fixed << (float) j / l * 100 << "%\r";
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
    while (k < end - begin)
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
            for (uint64_t k = start; k < start + seqChunkSize[thd_id]; k++)
            {

                if(ordValue(*(begin(seq[j]) + k + tshape.span - 1)) == 4)
                {
                    k += hashInit(tshape, begin(seq[j]) + k);

                    if (k > seqChunkSize[thd_id] - tshape.span + 1 + start)
                    {
                        k = seqChunkSize[thd_id] - (seqChunkSize[thd_id] + start) % thd_step + thd_step + start;
                    }
                }
                hashNext(tshape, begin(seq[j]) + k);
                //if (++ct_step != step)
                if (k % thd_step == 0)
                {
                    if (tshape.XValue ^ preX)
                    {
                        //if (ptr != 2)
                        //    printf("[debug]::ptr\n");
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
               // else 
               // {
               //     ct_step = 0;
               // }
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
    std::cout << "hs_size " << length(hs) << "\n";
    if (efficient)
    {
        clear(seq);
        shrinkToFit(seq);   
    }
    _DefaultHs.setHsHead(hs[hsRealEnd], 0, 0);
    std::cerr << "--Index::Initiate             Elapsed Time[s] " << sysTime() - time << " " << std::endl;
    std::cerr << "=>Index::SortHash                                                   \r";
    time = sysTime();
    _hsSort(begin(hs), begin(hs) + hsRealEnd, shape.weight, threads);
    std::cerr << "  Index::sortHash             Elapsed Time[s] " << sysTime() - time << std::endl;
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
                if (k >= kn)
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
                std::cerr << "sort y error " << k << " " << j;
                return false;
            }
            preY = hs[j];
        }
        //std::cout << std::endl;
        preX = _DefaultHs.getHeadX(hs[k]);
        k += _DefaultHs.getHeadPtr(hs[k]);
    }
    std::cerr << "last " << k << std::endl;
    
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
    //unsigned count = 0;
    //_setHeadNode(val, val);
//!!!!! need to modify;
    val = (xval << 2) + _DefaultXNodeBase.xHead;
    while (index.xstr.xstring[h1].val1)
    {
     //   ++count;
        //switch (index.xstr.xstring[h1].val1 ^ val) 
        switch(_DefaultXNodeFunc.collision(index.xstr.xstring[h1].val1, val))
        {
            case 0:
                //std::cerr << "case1\n";
                //return _DefaultXNodeFunc.makeReturnVal(index.xstr.xstring[h1]);
                return index.xstr.xstring[h1].val2;
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
    //std::cout << count << "\n";
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
    std::cerr << "request dir " << sysTime() - time << std::endl;
    return true;
}

/*
 * parallel sort ysa
 * this function is for index only collecting minihash value [minindex]
 */
 bool _createYSA(String<uint64_t> & hs, XString & xstr, uint64_t & indexEmptyDir, unsigned threads, uint64_t thd_blocklimit)
{

    uint64_t k = _DefaultHs.getHeadPtr(hs[0]);
    uint64_t preX = _DefaultHs.getHeadX(hs[0]);
    uint64_t ptr = k;
    uint64_t block_size = ptr;
    uint64_t countMove = 0, prek = 0;
    double time = sysTime();
    uint64_t countx = 0;
    std::vector<uint64_t> thd_hsStart(threads + 1, 0);
    std::cerr << "=>Index::SortYSA                                                  \r";
    while(_DefaultHs.getHeadPtr(hs[k]))
    {
        ptr = _DefaultHs.getHeadPtr(hs[k]);
        if (preX != _DefaultHs.getHeadX(hs[k]))
        {
            
            hs[k - countMove] = hs[k];
            _DefaultHs.setHsHeadPtr(hs[prek], block_size);
            //printf("[debug]::block %d\n", block_size);
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
        ++countx;
        k += ptr;
    }
    _DefaultHs.setHsHeadPtr(hs[prek], block_size);
    _DefaultHs.setHsHead(hs[k - countMove], 0, 0);
    _DefaultHs.setHsHead(hs[k - countMove + 1], 0, 0);

    thd_hsStart[threads] = prek + 1;
    resize(hs, k + 2 - countMove);
    shrinkToFit(hs);
    indexEmptyDir = k - countMove ;
    k=0;
    preX = _DefaultHs.getHeadX(hs[0]);
    ptr = 0;
    block_size = 0;

    #pragma omp parallel
    {
        uint64_t ptr = 0;
        #pragma omp for
        for (uint64_t k = 0; k < length(hs) - 2; k++)
        { 
            if(_DefaultHs.isHead(hs[k]))
            {
                ptr = _DefaultHs.getHeadPtr(hs[k]);
                _sort_YSA_Block(begin(hs) + k + 1, begin(hs) + k + ptr);
            }   
        }
    
    }
    std::cerr << "  Index::SortYSA              Elapsed Time[s] " << sysTime() - time << std::endl;
    std::cerr << "=>Index::resize xstr                                            \r" ;
    ptr = 0; k = 0;
    uint64_t count = 0; 
    
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
    std::cerr << "  Index::resize xstr          Elapsed Time[s] " << sysTime() - time << std::endl;
    std::cerr << "=>Index::request dir                                                  \r";
    time = sysTime();

    #pragma omp parallel 
    {
        //uint64_t thd_id = omp_get_thread_num();
        uint64_t ptr = 0;

    #pragma omp for 
        for (uint64_t m = 0; m < length(hs); m++)
        {

            if (_DefaultHs.isHead(hs[m]) && _DefaultHs.getHeadPtr(hs[m]))
            {
                ptr = _DefaultHs.getHeadPtr(hs[m]);
                //std::cout << "thd_blocklimit " << ptr << "\n";
                if (ptr < thd_blocklimit)
                {
                    for (unsigned j = m + 1; j < m + ptr; j++)
                    {
                        _DefaultHs.setHsBodyY(hs[j], 0);
                    }
                    requestXNode_noCollision_Atomic(xstr, _DefaultHs.getHeadX(hs[m]), 
                           m + 1, _DefaultXNodeBase.xHead, _DefaultXNodeBase.returnDir);   
                }
                else
                {
                    uint64_t xval = _DefaultHs.getHeadX(hs[m]);
                    //printf ("[debug]::hash %" PRIu64 " %d\n", xval, ptr);
                    
                    requestXNode_noCollision(xstr, xval, 
                        ~1, _DefaultXNodeBase.virtualHead, _DefaultXNodeBase.returnDir);
                    for (unsigned j = m + 1; j < m + ptr; j++)
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

    std::cerr << "  Index::request dir          Elapsed Time[s] " << sysTime() - time << std::endl;
    (void) threads;
    return true;
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

unsigned dshape_len = 21; //!!WARN::only odd number, even is not allowed
DIndex::DIndex():
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
String<int64_t> & DIndex::getHs()
{
    return hs;
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
                        unsigned gend
                        )
{
    double t = sysTime();
    LShape & shape = index.getShape();
    String<int> & dir = index.getDir();
    String<int64_t> & hs = index.getHs();
    resize (index.getDir(), index.fullSize(), 0);
    double t2 = sysTime();
    int64_t preVal = ~0;
    int64_t last_j = 0;
    for (int64_t i = 0; i < length(seqs); i++)
    {
        int64_t count = 0;
        hashInit (shape, begin(seqs[i]));
        for (int64_t j = 0; j < length(seqs[i]) - shape.span; j++)
        {
            hashNexth(shape, begin(seqs[i]) + j);
            if (++count > thd_min_step)
            {
                hashNextX(shape, begin(seqs[i]) + j);
                if (preVal != shape.XValue || j - last_j > thd_max_step)
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
    for (int64_t i = 0; i < length(dir); i++)
    {
        sum += dir[i];
        dir[i] = sum - dir[i];
    }
    last_j = 0;
    int64_t EmptyVal = create_cord(length(seqs),0,0,0); 
    //make sure genomeid >= length(seqs) and cord y be 0! y points to next empty.
    resize (hs, sum, EmptyVal);
    for (int64_t i = 0; i < length(seqs); i++)
    {
        int64_t count = 0; 
        hashInit (shape, begin(seqs[i]));
        for (int64_t j = 0; j < length(seqs[i]) - shape.span; j++)
        {
            hashNexth (shape, begin(seqs[i]) + j);
            if (++count > thd_min_step)
            {
                hashNextX (shape, begin(seqs[i]) + j);
                if (preVal != shape.XValue || j - last_j > thd_max_step)
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
    std::cout << "createDIndex " << sysTime() - t << " " << sysTime() - t2 << "\n";
}

/**
 * create index on [@gstr, @gend)th genomes 
 */
int createDIndex(StringSet<String<Dna5> > & seqs, 
                 DIndex & index, 
                 int64_t thd_min_step, 
                 int64_t thd_max_step,
                 int64_t thd_omit_block,
                 unsigned gstr,
                 unsigned gend,
                 unsigned threads
                )
{
    serr.print_message("=>Index::initing ", 0, 2, std::cerr);
    double t = sysTime();
    LShape & t_shape = index.getShape();
    String<int> & dir = index.getDir();
    String<int64_t> & hs = index.getHs();
    resize (dir, index.fullSize(), 0);
    double t2 = sysTime();
    dout << "idx2" << length(dir) << t_shape.weight << thd_min_step << thd_max_step << thd_omit_block<< "\n";
    dout << threads << "threads\n"; 
    for (int64_t i = gstr; i < gend; i++)
    {
        String<int64_t> t_blocks;
        for (int j = 0; j < threads; j++)
        {
            appendValue(t_blocks, length(seqs[i]) / threads * j); 
        }
        appendValue (t_blocks, length(seqs[i]) - t_shape.span);
        #pragma omp parallel
        {
            unsigned t_id = omp_get_thread_num();
            dout << "idx3 " << t_id << "\n";
            int64_t t_str = t_blocks[t_id];
            int64_t t_end = t_blocks[t_id + 1];
            int64_t preVal = ~0;
            int64_t last_j = t_str - 1;
            int64_t count = 0;
            LShape shape = t_shape;
            hashInit (shape, begin(seqs[i]) + t_str);
            //dout << "cdx1 " << t_str<< t_end <<"\n";
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
            }
        }
    }
    //double t4 = sysTime();
    int64_t sum = 0;
    for (int64_t i = 0; i < length(dir); i++)
    {
        if (dir[i] > thd_omit_block)
        {
            dir[i] = 0;
        }
        sum += dir[i];
        dir[i] = sum - dir[i];
    }
    dout << "dt" << sum << back(dir) << "\n";
    //dout << sysTime() - t4 << "x5\n";
    int64_t EmptyVal = create_cord(length(seqs),0,0,0); 
    //make sure genomeid >= length(seqs) and cord y be 0! y points to next empty.
    resize (hs, sum, EmptyVal);
    serr.print_message("--Index::inite   ", 0, 1, std::cerr);
    serr.print_message("=>Index::hashing", 0, 2, std::cerr);
    int64_t maxinfi = LLMAX;
    for (int64_t i = 0; i < length(seqs); i++)
    {
        String<int64_t> t_blocks;
        for (int j = 0; j < threads; j++)
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
            int64_t preVal = ~0;
            int64_t last_j = t_str - 1;
            int64_t count = 0;
            LShape shape = t_shape;
            hashInit (shape, begin(seqs[i]) + t_str);
            //dout << "strd " << t_str << t_end << "\n";
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
                            int64_t k = slot_str + get_cord_y(atomic_inc_cord_y(hs[slot_str])) - 1;
                            int64_t new_cord = create_cord(i, j, 0, shape.strand); //be sure new_cord_y == 0 
                            if (k == slot_str) 
                            {           //atomic creation the first cord which cotains shared pointer
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
            }  
        }
    }
    std::cout << "createDIndex " << sysTime() - t << " " << sysTime() - t2 << "\n";
    serr.print_message("Index::hash        ", 2, 1, std::cerr);
    serr.print_message("End createing index ", 2, 0, std::cerr);
    serr.print_message(sysTime() - t, 2, 1, std::cerr);
    return 0;
}

int64_t queryHsStr(DIndex & index, int64_t xval)
{
    return index.getDir()[xval];
}
int64_t queryHsEnd(DIndex & index, int64_t xval)
{
    return index.getDir()[xval + 1];
}

int IndexDynamic::isHIndex()
{
    return typeIx == typeHIx;
}

int IndexDynamic::isDIndex()
{
    return typeIx == typeDIx;
}

void IndexDynamic::setHIndex()
{
    typeIx = typeHIx;
}

void IndexDynamic::setDIndex()
{
    typeIx = typeDIx;
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
IndexDynamic::IndexDynamic(StringSet<String<Dna5> > & seqs):hindex()
{}

bool createIndexDynamic(StringSet<String<Dna5> > & seqs, IndexDynamic & index, unsigned gstr, unsigned gend, unsigned threads, bool efficient)
{

    if (index.isDIndex())
    {
        int64_t thd_min_step = 4;
        int64_t thd_max_step = 10;
        int64_t thd_omit_block = 50; 
        unsigned thd_shape_len = 21;
        index.hindex.shape.init_shape_parm(thd_shape_len);
        std::cout << "cidx\n";
        //TODO::parm wrapping 
        return createDIndex(seqs, index.dindex, 
                            thd_min_step, 
                            thd_max_step, 
                            thd_omit_block,
                            gstr, gend,
                            threads);
//       return createDIndex_serial(seqs, index.dindex, 4, 10);
    }
    else if (index.isHIndex())
    {
        //chr1: step=5, shape_len=25, block_limt =32|  0.077% > 32(block size) 
        //chr1, 5,  21, 32 | 0.7% > 32 0.22% > 64
        //chr1, 10, 21, 32 | 0.49% > 32  0.16% > 64
        unsigned thd_step = 10;
        unsigned thd_shape_len = 25; //WARN only odd number is allowed due to double strand hash
        uint64_t thd_blocklimit = 32;
        float alpha = 1.6;
        index.hindex.shape.init_shape_parm(thd_shape_len / 2 * 2 + 1);
        dout << "span" << index.hindex.shape.span << "\n";
        index.hindex.alpha = alpha;
        return createHIndex(seqs, index.hindex, 
                            gstr, gend, 
                            threads, thd_step, thd_blocklimit, efficient);
    }
    /*
    else if (index.isMHIndex()) //Hindex part of Mix of DIndex and HIndex, super! 
    {

    }
    else if (index.isMDIndex()) //DIndex part of Mix of DIndex and HIndex, !
    {

    }
    */
}

