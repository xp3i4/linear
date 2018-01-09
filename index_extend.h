// ==========================================================================
//                          Mappeing SMRT reads
// ==========================================================================
// Copyright (c) 2006-2016, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: cxpan <chenxu.pan@fu-berlin.de>
// ==========================================================================

#ifndef SEQAN_HEADER_INDEX_PM_H
#define SEQAN_HEADER_INDEX_PM_H

namespace seqan{

//x-begin: min shape open index 
 template < typename TObject, unsigned TSPAN, unsigned TWEIGHT>
    class Index<TObject, IndexQGram<Minimizer<TSPAN, TWEIGHT>, OpenAddressing> >
    {
    private:
        static const double defaultAlpha;
    public:
        typedef typename Member<Index, QGramText>::Type     TTextMember;
        typedef typename Fibre<Index, QGramText>::Type        TText;
        typedef typename Fibre<Index, QGramSA>::Type        TSA;
        typedef typename Fibre<Index, QGramDir>::Type        TDir;
        typedef typename Fibre<Index, QGramCounts>::Type    TCounts;
        typedef typename Fibre<Index, QGramCountsDir>::Type    TCountsDir;
        typedef typename Fibre<Index, QGramShape>::Type        TShape;
        typedef typename Fibre<Index, QGramBucketMap>::Type    TBucketMap;
        typedef typename Cargo<Index>::Type                    TCargo;
        typedef typename Size<Index>::Type                    TSize;

        TTextMember     text;        // underlying text
//        TSA                sa;            // suffix array sorted by the first q chars
        String<uint64_t> sa;
        TDir            dir;        // bucket directory
        TCounts            counts;        // counts each q-gram per sequence
        TCountsDir        countsDir;    // directory for count buckets
        TShape            shape;        // underlying shape
        TCargo            cargo;        // user-defined cargo
        TBucketMap        bucketMap;    // bucketMap table (used by open-addressing index)
        TSize            stepSize;    // store every <stepSize>'th q-gram in the index

        double            alpha;        // for m entries the hash map has at least size alpha*m
    //x-begin//
        uint64_t           start;          // 
        uint64_t         _Empty_Dir_ = start - 2;  
    //x-end//
        Index():
            stepSize(1),
            alpha(defaultAlpha) {}

        Index(Index &other):
            text(other.text),
            sa(other.sa),
            dir(other.dir),
            counts(other.counts),
            countsDir(other.countsDir),
            shape(other.shape),
            cargo(other.cargo),
            bucketMap(other.bucketMap),
            stepSize(1),
            alpha(defaultAlpha) {}

        Index(Index const &other):
            text(other.text),
            sa(other.sa),
            dir(other.dir),
            counts(other.counts),
            countsDir(other.countsDir),
            shape(other.shape),
            cargo(other.cargo),
            bucketMap(other.bucketMap),
            stepSize(1),
            alpha(defaultAlpha) {}

        template <typename TText_>
        Index(TText_ &_text):
            text(_text),
            stepSize(1),
            alpha(defaultAlpha) {}

        template <typename TText_>
        Index(TText_ const &_text):
            text(_text),
            stepSize(1),
            alpha(defaultAlpha) {}

        template <typename TText_, typename TShape_>
        Index(TText_ &_text, TShape_ const &_shape):
            text(_text),
            shape(_shape),
            stepSize(1),
            alpha(defaultAlpha) {}

        template <typename TText_, typename TShape_>
        Index(TText_ const &_text, TShape_ const &_shape):
            text(_text),
            shape(_shape),
            stepSize(1),
            alpha(defaultAlpha) {}
    };
template <typename TObject, unsigned TSPAN, unsigned TWEIGHT>
 const double Index<TObject, IndexQGram<Minimizer<TSPAN, TWEIGHT>, OpenAddressing> >::defaultAlpha = 1.6;


//=============================================================
// definition of some types for mimizer index
    struct hPair
    {
        uint64_t i1;
        uint64_t i2;
        inline hPair & operator = (hPair & b)
        {
            i1 = b.i1;
            i2 = b.i2;
            return *this;
        }
    };
//==============================================================

//=============================================================
//definition of types of elements in dir (uint64_t)
    //bodyNode=
    //length1[4] length2[38] YValue[20] code[2]: code={1}
    //
    //headNode=
    //hvalue(/XValue)[62] code[2]: code={2,3}, 2:head, 3:virtual head
    //
    //code  0 empty
    //code  2 head 10
    //code  3 virtual head 11
    //code  1 body  01

    static const uint64_t _dirEmpty = (uint64_t)1 << 63;
    static const uint64_t _bitEmpty = 0;
    static const uint64_t _bitLength = ((uint64_t)1 << 58) - 1;
    static const uint64_t _bitValue = ((uint64_t)1 << 56) - 1;
    static const unsigned _bitLength_END = 60;
    static const unsigned _bitValue_END = 2;
    
    //HeadNode(H):
    //H :=  (h/X)Value[62]|HeadType[2]: 0:=empty 2:=head 3:=virtual head
    static const unsigned _HeadValue_bits = 3;
    static const uint64_t _HeadType_code = 2;
    static const uint64_t _HeadTypeVtl_code = 3;
    static const uint64_t _HeadTypeHVl_code = 4;
    // BodyNode(B):
    // B := YValue[23]|BodyType[1]: |counth[40]
    // occ = counth[n+1] - count[n] = (B[n+1] - B[n]) & bit, bit = 00..0011...11
    static const unsigned _BodyValue_bits = 41;
    static const unsigned _BodyType_bits = 40;
    //static const unsigned _BodyValue_bits = 2;
    static const uint64_t _BodyType_code = 1;
    static const uint64_t _BodyTypeEnd_code = 0;
    static const uint64_t _BodyType_key = ~((uint64_t)1 << _BodyType_bits);
    static const uint64_t _getBody = ((uint64_t)1 << _BodyType_bits) - 1;

    static const uint64_t _bitEmptyType = 0;

    static const uint64_t _bitCode = 3;                                 // node(value) & _bitCode to acquire the type of the node 
    //static const uint64_t _bitValue2 = ((uint64_t)1 << 22) - 1;
    static const uint64_t _bitEmptyCode = 0;
    static const uint64_t _bitBodyCode = 1;
    static const uint64_t _bitHeadCode = 2;
    static const uint64_t _bitVtlHeadCode = 3;
  
    //SA node:= seq num i1[10]| base num i2[30]
    static uint64_t _BaseNum_bits = 30 ;    
    static uint64_t _SeqNum_bits = _BodyType_bits - _BaseNum_bits;    
    static uint64_t _BaseNum_code = ((uint64_t)1 << _BaseNum_bits) - 1;
    static uint64_t _BaseNum_SeqMask = (1ULL << _SeqNum_bits) - 1;
        
 
    static const uint64_t _Empty_Dir_ = -1;

    static const unsigned blocklimit = 32;
    
//==============================================================
    template <typename HValue>
    inline HValue _makeHeadNode(HValue code)
    {
        return (code << _HeadValue_bits) + _HeadType_code;
    } 
    template <typename HValue>
    inline HValue _makeVtlHeadNode(HValue code)
    {
        return (code << _HeadValue_bits) + _HeadTypeVtl_code;
    }
    template <typename HValue>
    inline HValue _makeHVlHeadNode(HValue code)
    {
        return (code << _HeadValue_bits) + _HeadTypeHVl_code; 
    }
    template <typename HValue>
    inline void _setHVlHeadNode(HValue & headNode, HValue const & hValue)
    {
        headNode = (hValue << _HeadValue_bits) + _HeadTypeHVl_code; 
    }
    template <typename HValue>
    inline void _setHeadNode(HValue & headNode, HValue const & hValue)
    {
        headNode = (hValue << _HeadValue_bits) + _HeadType_code;
    }
    
    template <typename HValue>
    inline HValue _makeEmptyNode(HValue code)
    {
        return (code << _HeadValue_bits) + _bitEmptyType;
    }
    template <typename HValue, unsigned TSPAN, unsigned TWEIGHT>
    inline HValue _getDirStart(Index<StringSet<DnaString>, IndexQGram<Minimizer<TSPAN, TWEIGHT>, OpenAddressing> >  & index)
    {
        return index.start;
    }
    template <typename HValue>
    inline void _setBodyType_Begin(HValue & code){
        code &= _BodyType_key; 
    }

    template <typename HValue>
    inline HValue _ifBodyType(HValue code){
        return code & (~_BodyType_key);
    }
    
    template <typename HValue>
    inline HValue _getHeadValue(HValue  code)
    {
        return code >> _HeadValue_bits;  
    }
    template <typename HValue>
    inline void _setBodyNode(HValue & bodyNode, HValue const & YValue, HValue const & type, HValue const & counth)
    {
        bodyNode = (YValue << _BodyValue_bits) + (type << _BodyType_bits) + counth;
    }
    template <typename HValue>
    inline HValue _getBodyValue(HValue code)
    {
        return code >> _BodyValue_bits;
    }
    template <typename HValue>
    inline HValue _getBodyCounth(HValue & code)
    {
        return code & _getBody;
    }

    template <typename HValue>
    inline HValue _createSANode(HValue const & i1, HValue const & i2)
    {
        return (i1 << _BaseNum_bits) + i2;
    }

    template <typename HValue>
    inline void _setSANode(HValue & node, HValue const & i1, HValue const & i2)
    {
        node = (i1 << _BaseNum_bits) + i2;
    }
    

    template <typename HValue> 
    inline HValue _getSA_i1(HValue const & node)
    {
        return (node >> _BaseNum_bits) & _BaseNum_SeqMask;
    }
    template <typename HValue>
    inline HValue _getSA_i2(HValue const & node)
    {
        return node & _BaseNum_code;
    }
    
//x-end: min shape open index 


   
    //x1-begin
    //template <typename TIndex, typename THashValue, typename TParallelTag>

    template < typename TBucketMap, typename TValue >
    inline TValue
    _hashFunction1(TBucketMap const &, TValue val)
    {
    uint64_t key = val;
          key = (~key) + (key << 21); // key = (key << 21) - key - 1;
  key = key ^ (key >> 24);
  key = (key + (key << 3)) + (key << 8); // key * 265
  key = key ^ (key >> 14);
  key = (key + (key << 2)) + (key << 4); // key * 21
  key = key ^ (key >> 28);
  key = key + (key << 31);
  return key;        
    }

    template<typename TDir, typename THashValue>
    inline THashValue
    requestDir(TDir & dir, THashValue hlen, THashValue code, THashValue code1)//, Tag<TParallelTag> parallelTag)
    //code:headNode, code1:pointerTobody
    {
        //std::cerr << code << " " << len << std::endl;
        typedef unsigned long TSize;
        //TSize hlen = 536870913 ;
        //TSize hlen = 2147483649;
        if (hlen == 0ul) return code;
        
        TSize h1 = _hashFunction1(dir, _getHeadValue(code));
#ifdef SEQAN_OPENADDRESSING_COMPACT
        --hlen;
        h1 %= hlen;
#else
        hlen -= 2;
        h1 &= hlen;
#endif
        TSize delta = 0;
        (void)delta;
        while(dir[h1] | dir[h1+1]) 
        {
            switch(code ^ dir[h1]){
                case 0:
                    return h1;
                case 1:
                    return h1;
                default:
                    h1 = (h1 + delta + 1) & hlen;
                    delta++;
            }
        }
        dir[h1] = code;
        dir[h1 + 1] = code1;
        return h1;
    }

    //x-end2:
    
    template <typename TObject, unsigned TSPAN, unsigned TWEIGHT, typename TValue, typename TSpec>
    inline typename Value< Shape<TValue, Minimizer<TSPAN, TWEIGHT, TSpec> > >::Type 
    getDir(Index<TObject, IndexQGram<Minimizer<TSPAN, TWEIGHT>, OpenAddressing> > const & index, Shape<TValue, Minimizer<TSPAN, TWEIGHT, TSpec> > const & shape)
    {
        typedef unsigned long TSize;
        typedef typename Value< Shape<TValue, Minimizer<TSPAN, TWEIGHT, TSpec> > >::Type THashValue;
        // get size of the index

        // check whether bucket map is disabled and
        // where the hash should be found if no collision took place before
    
        THashValue key, it;
        TSize hlen = index.start - 2;
        TSize h1 = _hashFunction1(index.dir, shape.XValue) & hlen;
        //if (hlen == 0ul) return index._Empty_Dir_;

//#ifdef SEQAN_OPENADDRESSING_COMPACT
//        --hlen;
//        h1 %= hlen;
//#else
        //hlen -= 2;
        //h1 &= hlen;
//#endif
        TSize delta = 0;
        (void)delta;
        _setHeadNode(key,shape.XValue);
        while (index.dir[h1] | index.dir[h1+1])
        {
            switch (index.dir[h1] ^ key) 
            {
                case 0:
                    it = _getHeadValue(index.dir[h1+1]);
                    do{
                    
                        if (shape.YValue ==  _getBodyValue(index.dir[it]))
                        {    
                            return it;
                        } 
                    }while(_ifBodyType(index.dir[++it])); //until the begin of next block
                    return index._Empty_Dir_ ;
                case 1:
                    _setHVlHeadNode(key, shape.hValue);
                    h1 = _hashFunction1(index.dir, shape.hValue) & hlen;
                    delta = 0;
                    break;
                default:
                    h1 = (h1 + delta + 1) & hlen;
                    delta++;
            }
        }
        return index._Empty_Dir_; 
    }

    
    template <typename TObject, unsigned TSPAN, unsigned TWEIGHT>
    inline __int64 _fullDirLength(Index<TObject, IndexQGram<Minimizer<TSPAN, TWEIGHT>, OpenAddressing> > const &index)
    {
        typedef Index<TObject, IndexQGram<Minimizer<TSPAN, TWEIGHT>, OpenAddressing> >    TIndex;
        //typedef typename Fibre<TIndex, FibreShape>::Type                    TShape;
        //typedef typename Host<TShape>::Type                                    TTextValue;

        double num_qgrams = _qgramQGramCount(index) * index.alpha;
        //double max_qgrams = 2*pow((double)ValueSize<TTextValue>::VALUE, (double)length(indexShape(index)));
        __int64 qgrams;

        qgrams = (__int64)ceil(num_qgrams);
#ifndef SEQAN_OPENADDRESSING_COMPACT
        __int64 power2 = 1;
        while (power2 < qgrams)
            power2 <<= 1;
        qgrams = power2;
    #endif
            resize(const_cast<TIndex &>(index).bucketMap.qgramCode, qgrams + 1, Exact());
        return qgrams + 1;
    }

template <typename TObj, unsigned TSpan, unsigned TWeight>
void _qgramClearDir(Index<StringSet<String<TObj> >, IndexQGram<Minimizer<TSpan, TWeight>, OpenAddressing> > & index)
{
    typedef Shape<Dna, Minimizer<TSpan, TWeight> > TM_Shape;
    typedef typename Value<TM_Shape>::Type HValue;
    if (length(indexText(index)) > (1ULL<<_SeqNum_bits))
    {
        _BaseNum_bits = 17;    
        _SeqNum_bits = _BodyType_bits - _BaseNum_bits;
        _BaseNum_code = ((uint64_t)1 << _BaseNum_bits) - 1;

    }
    resize (indexDir(index), _fullDirLength(index) + lengthSum(indexText(index)) + 2);
    index.start = _fullDirLength(index);
    index._Empty_Dir_ = 0;
    for (HValue k = 0; k < length(index.dir); k++) 
    {
        index.dir[k] = _bitEmpty;
    }
    std::cerr << "    _qgramClearDir():" << std::endl;
    std::cerr << "        _fullDirLength(index) = " << _fullDirLength(index) << std::endl;
    std::cerr << "        lengh(index.dir) = " << length(index.dir) << std::endl;
    std::cerr << "        End _qgramClearDir()" << std::endl;
}
/*
template <unsigned TSpan, unsigned TWeight>
void _qgramCountQGrams2(Index<StringSet<DnaString>, IndexQGram<Minimizer<TSpan, TWeight>, OpenAddressing > > & index)
{
    typedef Shape<Dna, Minimizer<TSpan, TWeight> > TM_Shape;
    typedef Iterator<String<Dna> >::Type TIter;
    typedef typename Value<TM_Shape>::Type HValue;
    typedef std::tuple<HValue, HValue, HValue, HValue> HTuple;
    typedef String<HTuple> Stringtuple;

    TM_Shape shape;
    Stringtuple hs, hs1;
    HValue  m = 0, sum = 0;

    double time = sysTime();
    resize(hs, lengthSum(indexText(index)) - length(indexText(index)) * (shape.span - 1) + 1);

    std::cerr << "        _qgramCountQGrams() sysTime(): " << sysTime() - time << std::endl;
    std::cerr << "            lengthSum(StringSet) = " << lengthSum(indexText(index)) << std::endl;
    for(HValue k = 0; k < length(indexText(index)); k++)
    {
        TIter it = begin(indexText(index)[k]);
        hashInit(shape, it);
        for (HValue j = 0; j < length(indexText(index)[k]) - shape.span + 1; j++)
        {
            hashNext(shape, it + j);
            hs[m++] = std::make_tuple(shape.XValue, shape.hValue, shape.YValue, _createSANode(k, j));
        }
    }
    std::cerr << "            make_tuple sysTime(): " << sysTime() - time << std::endl;
    //hs[length(hs) - 1] = std::make_tuple((HValue)0, (HValue)0, (HValue)0, (HValue)1);
    std::sort(begin(hs), end(hs) - 1,
        [](const HTuple &a, const HTuple &b)
        {return (std::get<0>(a) > std::get<0>(b)||(std::get<0>(a) == std::get<0>(b) && std::get<1>(a) > std::get<1>(b)));});
    hs[length(hs) - 1] = std::make_tuple((HValue)1, (HValue)1, (HValue)0, (HValue)1);
    std::cerr << "            sort sysTime(): " << sysTime() - time << " " << std::get<0>(hs[length(hs) - 2]) << " " << std::get<1>(hs[length(hs)- 2]) << std::endl;
    HValue countx = 1, counth = 1, tmp = 0, countdh = 0, countb = 0, hk = 0;
    resize(index.sa, length(hs) - 1);
    for (HValue k = 1; k < length(hs); k++)
    {
        index.sa[k - 1] = std::get<3>(hs[k-1]);
        if (std::get<1>(hs[k]) != std::get<1>(hs[k - 1]))
        { _setBodyNode(index.dir[index.start + hk], std::get<2>(hs[k-1]), _BodyType_code, tmp);
        if (std::get<0>(hs[k - 1]) == 0 && std::get<2>(hs[k - 1]) == 441204)
            std::cerr <<"_getBodyValue = " << _getBodyValue(index.dir[index.start + hk]) << std::endl;
            hk++;
            countb++;
            countdh++;
            tmp = counth;
        }
        //else
        counth++;

        if (std::get<0>(hs[k]) != std::get<0>(hs[k - 1]))
        {
            if (countb < blocklimit)
            {
                requestDir(index.dir, index.start, _makeHeadNode(std::get<0>(hs[k-1])), _makeEmptyNode(index.start + hk - countb));
                for (HValue j = 0; j < countb; j++)
                    _setBodyType_Begin(index.dir[index.start + hk - countb]);
            }
            else
            {
                hk -= countb;
                requestDir(index.dir, index.start, _makeVtlHeadNode(std::get<0>(hs[k-1])), _makeEmptyNode(index.start + hk));
                for (HValue j = k - countx; j < k; j++)
                    if (std::get<1>(hs[j]) != std::get<1>(hs[j + 1]))
                    {
                        requestDir(index.dir, index.start, _makeHVlHeadNode(std::get<1>(hs[j])), _makeEmptyNode(index.start+hk));
                        _setBodyType_Begin(index.dir[index.start + hk]);
                        hk++;
                    }
            }
            countb = 0;
            countx = 1;
        }
        else
        {
            countx++;
        }
    }
    std::cerr << std::endl;
    std::cerr << counth << std::endl;
    resize(index.dir, index.start + countdh + 10);
    _setBodyNode(index.dir[index.start + countdh], _bitEmpty, _BodyType_code, counth - 1); 
    _setBodyType_Begin(index.dir[index.start + countdh]);
    index._Empty_Dir_ = index.start + countdh + 1;
    //_setBodyNode(index.dir[index.start + countdh], _bitEmpty, _BodyTypeEnd_code, counth - 1); 
    //_setBodyType_Begin(index.dir[index.start + countdh]);
    //index._Empty_Dir_ = index.start + countdh;
    //_setBodyNode(index.dir[index.start + countdh + 1], _bitEmpty, _BodyTypeEnd_code, counth - 1); 
    //_setBodyType_Begin(index.dir[index.start + countdh + 1]);
    std::cerr << "            End _qgramCountQGrams() sysTime(): " << sysTime() - time << std::endl;
}



template <unsigned TSpan, unsigned TWeight>
void _qgramCountQGrams3(Index<StringSet<DnaString  >, IndexQGram<Minimizer<TSpan, TWeight>, OpenAddressing > > & index)
{
    typedef Shape<Dna, Minimizer<TSpan, TWeight> > TM_Shape;
    typedef Iterator<String<Dna> >::Type TIter;
    typedef typename Value<TM_Shape>::Type HValue;
    typedef std::tuple<HValue, HValue, HValue> HTuple;
    typedef String<HTuple> Stringtuple;

    TM_Shape shape;
    Stringtuple hs, hs1;
    HValue  m = 0, sum = 0;

    double time = sysTime();
    resize(hs, lengthSum(indexText(index)) - length(indexText(index)) * (shape.span - 1) + 1);

    std::cerr << "        _qgramCountQGrams() sysTime(): " << sysTime() - time << std::endl;
    std::cerr << "            lengthSum(StringSet) = " << lengthSum(indexText(index)) << std::endl;
    for(HValue k = 0; k < length(indexText(index)); k++)
    {
        TIter it = begin(indexText(index)[k]);
        hashInit(shape, it);
        for (HValue j = 0; j < length(indexText(index)[k]) - shape.span + 1; j++)
        {
            hashNext(shape, it + j);
            hs[m++] = std::make_tuple(shape.XValue, shape.YValue, _createSANode(k, j));
        }
    }
    std::cerr << "            make_tuple sysTime(): " << sysTime() - time << std::endl;
    //hs[length(hs) - 1] = std::make_tuple((HValue)0, (HValue)0, (HValue)0, (HValue)1);
    std::sort(begin(hs), end(hs) - 1,
        [](const HTuple &a, const HTuple &b)
        {return (std::get<0>(a) < std::get<0>(b)||(std::get<0>(a) == std::get<0>(b) && std::get<1>(a) > std::get<1>(b)));});
    hs[length(hs) - 1] = std::make_tuple((HValue)1, (HValue)0, (HValue)1);
    std::cerr << "            sort sysTime(): " << sysTime() - time << std::endl;
    HValue countx = 1, counth = 1, tmp = 0, countdh = 0, countb = 1, hk = 0;
    resize(index.sa, length(hs) - 1);
    for (HValue k = 1; k < length(hs); k++)
    {
        if (std::get<0>(hs[k]) ^ std::get<0>(hs[k - 1])|std::get<1>(hs[k]) ^ std::get<1>(hs[k - 1]))
        { _setBodyNode(index.dir[index.start + hk], std::get<1>(hs[k-1]), _BodyType_code, tmp);
            hk++;
            countb++;
            countdh++;
            tmp = counth;
        }
        if (std::get<0>(hs[k]) ^ std::get<0>(hs[k - 1]))
        {
            
            if (countb < blocklimit)
            {
                requestDir(index.dir, index.start, _makeHeadNode(std::get<0>(hs[k-1])), _makeEmptyNode(index.start + hk - countb));
                for (HValue j = 0; j < countb; j++)
                    _setBodyType_Begin(index.dir[index.start + hk - countb]);
            }
            else
            {
                hk -= countb;
                requestDir(index.dir, index.start, _makeVtlHeadNode(std::get<0>(hs[k-1])), _makeEmptyNode(index.start + hk));
                for (HValue j = k - countx; j < k; j++)
                    if ((std::get<0>(hs[j]) ^ std::get<0>(hs[j + 1])) | (std::get<1>(hs[j]) ^ std::get<1>(hs[j + 1])))
                    {
                        requestDir(index.dir, index.start, _makeHVlHeadNode(xy2h(shape, std::get<0>(hs[j]),std::get<1>(hs[j]))), _makeEmptyNode(index.start+hk));
                        _setBodyType_Begin(index.dir[index.start + hk]);
                        hk++;
                    }
            }
            countb = 0;
            countx = 1;
        }
        else
        {
            countx++;
        }
        
        index.sa[k - 1] = std::get<2>(hs[k-1]);
        counth++;
    }
    std::cerr << std::endl;
    std::cerr << counth << std::endl;
    resize(index.dir, index.start + countdh + 10);
    _setBodyNode(index.dir[index.start + countdh], _bitEmpty, _BodyType_code, counth - 1); 
    _setBodyType_Begin(index.dir[index.start + countdh]);
    index._Empty_Dir_ = index.start + countdh + 1;
    //_setBodyNode(index.dir[index.start + countdh], _bitEmpty, _BodyTypeEnd_code, counth - 1); 
    //_setBodyType_Begin(index.dir[index.start + countdh]);
    //index._Empty_Dir_ = index.start + countdh;
    //_setBodyNode(index.dir[index.start + countdh + 1], _bitEmpty, _BodyTypeEnd_code, counth - 1); 
    //_setBodyType_Begin(index.dir[index.start + countdh + 1]);
    std::cerr << "            End _qgramCountQGrams() sysTime(): " << sysTime() - time << std::endl;
}
*/

//template <typename TObj, typename Compare>
template <typename TIter>
//inline void _mergeSort(String<Iterator<TObj>::Type> & begin, String<Iterator<TObj>::Type> & end, Compare cmp)
inline void _mergeSort(TIter const & it, String<unsigned> begin, String<unsigned> end)
{
    //typedef typename Iterator<TObj>::Type TIter;
    //typedef typename Value<TObj>::Type _Value;
    String<Pair<uint64_t, uint64_t> >tmp;
    uint64_t max;
    resize(tmp, end[length(end) -1 ] - begin[0]);
    unsigned maxk = 0;
    for (unsigned j = 0; j < length(tmp); j++)
    { 
        max = 0;
        for (unsigned k = 0; k < length(begin); k++) 
        {
            if (begin[k] != end[k])
                if((it + begin[k])->i1 > max)
                {
                    max = (it + begin[k])->i1;
                    maxk = k;
                }
        }
        tmp[j] = *(it+begin[maxk]);
        begin[maxk] += 1;
    } 
}
/*
//template <typename TObj, typename Compare>
template <typename TIter>
inline void 
//_radixSort(Iterator<TObj>::Type const & begin, Iterator<TObj>::Type const & end, Compare compare)
_radixSort(TIter const & begin,  TIter const & end, 
            unsigned const p_bit, unsigned const & l)
{
    unsigned  l_move = 64, r_move = 64 - p_bit;
    //uint64_t count[1<<p_bit];
    uint64_t count[1024];
    //int count[1024];
    String<Pair<uint64_t, uint64_t> > output;
    resize(output, end - begin);
    for (uint64_t j = 0; j < l; j++)
    {
        l_move -= p_bit;
        for (int k = 0; k< (1<<p_bit); k++)
            count[k]=0;
        for (int64_t k = 0; k < end - begin; k++)
            count[(begin + k)->i1 << l_move >> r_move]++;
        for (int k = 1; k < (1 << p_bit); k++)
            count[k] += count[k - 1];
        for (int64_t k = end - begin - 1; k >=0; k-- )
            output[--count[(begin + k)->i1 << l_move >> r_move]] = *(begin + k);
        for (int64_t k = 0; k < end - begin; k++)
            *(begin + k)  = output[k];
    }
}
*/

template <typename TIter>
inline void 
//_radixSort(Iterator<TObj>::Type const & begin, Iterator<TObj>::Type const & end, Compare compare)
_radixSort(TIter const & begin,  TIter const & end, 
            unsigned const p_bit, unsigned const & l)
{
    unsigned  l_move = 64, r_move = 64 - p_bit;
    //uint64_t count[1<<p_bit];
    uint64_t count[1024];
    //int count[1024];
    String<Pair<uint64_t, uint64_t> > output;
    resize(output, end - begin);
    TIter begin1 = begin, end1 = end;
    for (uint64_t j = 0; j < l; j++)
    {
        l_move -= p_bit;
        for (int k = 0; k< (1<<p_bit); k++)
            count[k]=0;
        for (int64_t k = 0; k < end - begin; k++)
            count[(begin + k)->i1 << l_move >> r_move]++;
        for (int k = 1; k < (1 << p_bit); k++)
            count[k] += count[k - 1];
        for (int64_t k = end - begin - 1; k >=0; k-- )
            output[--count[(begin + k)->i1 << l_move >> r_move]] = *(begin + k);
        //for (int64_t k = 0; k < end - begin; k++)
        //    *(begin + k)  = output[k];
        std::swap(begin, begin1);
        std::swap(end, end1);
    }
}



template <typename TIter>
inline void RMSort(TIter const & begin, TIter const & end)
{
    unsigned _radixBlock = 32 << 20, numBlock = ceil((end - begin)/_radixBlock);
    String<unsigned> segb, sege;
    resize(segb, numBlock);
    resize(sege, numBlock);
    for (unsigned k = 0; k < numBlock; k++)
    {
        segb[k] = k * _radixBlock;
        if (k != numBlock)
            sege[k] = segb[k] + _radixBlock;
        else
            sege[k] = end - begin;
        _radixSort(begin + segb[k], begin + sege[k], 9, 5);
    }
    _mergeSort(begin, segb, sege);
}

template <typename TIt>
inline void _insertSort(TIt const & begin, TIt const & end )
{
    Pair<uint64_t, uint64_t> key;
    for (int j = 1; j < end - begin; j++)
    {
        key = *(begin + j); 
        int k = j - 1;
        while (k >= 0)
        {
            if (((begin + k)->i2 < key.i2))
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

template <typename TIt>
inline void _sort3(TIt const & begin, const TIt & end, unsigned const & p_bit, unsigned const & l)
{
    unsigned  l_move = 64, r_move = 64 - p_bit;
    //uint64_t count[1<<p_bit];
    uint64_t count[1024];
    //int count[1024];
    String<Pair<uint64_t, uint64_t> > output;
    resize(output, end - begin);
    for (uint64_t j = 0; j < l; j++)
    {
        l_move -= p_bit;
        for (int k = 0; k< (1<<p_bit); k++)
            count[k]=0;
        for (int64_t k = 0; k < end - begin; k++)
            count[(begin + k)->i1 << l_move >> r_move]++;
        for (int k = 1; k < (1 << p_bit); k++)
            count[k] += count[k - 1];
        for (int64_t k = end - begin - 1; k >=0; k-- )
            output[--count[(begin + k)->i1 << l_move >> r_move]] = *(begin + k);
        for (int64_t k = 0; k < end - begin; k++)
            *(begin + k)  = output[k];
    }
}

template <typename TIt>
inline void _sort3_i2_(TIt const & begin, const TIt & end, unsigned const & p_bit, unsigned const & l)
{
    unsigned  l_move = 64, r_move = 64 - p_bit;
    //uint64_t count[1<<p_bit];
    uint64_t count[1024];
    //int count[1024];
    String<Pair<uint64_t, uint64_t> > output;
    resize(output, end - begin);
    for (uint64_t j = 0; j < l; j++)
    {
        l_move -= p_bit;
        for (int k = 0; k< (1<<p_bit); k++)
            count[k]=0;
        for (int64_t k = 0; k < end - begin; k++)
            count[(begin + k)->i2 << l_move >> r_move]++;
        for (int k = 1; k < (1 << p_bit); k++)
            count[k] += count[k - 1];
        for (int64_t k = end - begin - 1; k >=0; k-- )
            output[--count[(begin + k)->i2 << l_move >> r_move]] = *(begin + k);
        for (int64_t k = 0; k < end - begin; k++)
            *(begin + k)  = output[k];
    }
}



template <unsigned TSPAN, unsigned TWEIGHT>
void _createValueArray2(StringSet<DnaString> & reads, String<Pair<uint64_t, uint64_t> > & hs, Shape<Dna, Minimizer<TSPAN, TWEIGHT> > & shape, int step, int l)
{
    typedef String<Pair<uint64_t, uint64_t> > StringPair;

    //TShape shape;
    StringPair tmp;
    String<uint64_t> tmp3;
    uint64_t p = 0, q=0, c = 0, n = -1, pre = ~0, count = 0, mask = ((uint64_t)1 << 63);
    //std::cerr << "    _createValueArray() " << std::endl;
    double time = sysTime();
    resize(tmp, lengthSum(reads) - length(reads) * (shape.span - 1));
    resize(tmp3, lengthSum(reads) - length(reads) * (shape.span - 1));
    for(uint64_t j = 0; j < length(reads); j++)
    {
        hashInit(shape, begin(reads[j]));
        for (uint64_t k = 0; k < length(reads[j]) - shape.span + 1; k++)
        {
            hashNext(shape, begin(reads[j]) + k);
            _setBodyNode(tmp3[p], shape.YValue, (uint64_t)1, _createSANode(j, k));
            if (pre ^ shape.XValue)
            {
                tmp[q].i1 = shape.XValue;
                tmp[q].i2 = p;
                pre = shape.XValue;
                tmp3[p] |= mask;
                q++;
            }
            p++;
        }
    }
    resize(tmp, q);
    *(end(tmp3)) |= (mask);
    std::cerr << "        loading time " << sysTime() - time << std::endl;
    c = tmp[0].i2;
    p = q = count = 0;
    n = -1;
    _sort3(begin(tmp), end(tmp), step, l);         // sort parameters 1
    std::cerr << "        sort xvalue time " << sysTime() - time << std::endl;
    c = tmp[0].i2;
    for (uint64_t q = 0;  q < length(hs) - 1; q++)
    {
        if (tmp3[c] & mask)
        {
            c = tmp[++n].i2;
        }
        hs[q].i1 = tmp[n].i1;
        hs[q].i2 = tmp3[c] & (~mask);
        c++;
        //std::cerr << c << " " << length(tmp3) << " " << q << " " << length(hs) << std::endl;
    }
    std::cerr << "        xvalue expand " << sysTime() - time << std::endl;
    hs[length(hs)-1].i2 |= mask;
    pre = hs[0].i1;
    for (uint64_t k = 0; k < length(hs); k++)
    {
        if (pre ^ hs[k].i1)
        {
            pre = hs[k].i1;
            if (count < 20)                   // sort parameters
                _insertSort(begin(hs) + k - count, begin(hs) + k);
            else
                //std::sort(begin(hs) + k -count, begin(hs) + k, [](Pair<uint64_t, uint64_t> & a,
                //Pair<uint64_t, uint64_t> & b){return a.i2 > b.i2;});
    
                //std::stable_sort(begin(hs) + k -count, begin(hs) + k, comp);
                _sort3_i2_(begin(hs) + k - count, begin(hs) + k,8,8);

            count = 0;
        }
        count++;
   }


   std::cerr << "        End sort sysTime(): " <<  sysTime() - time << std::endl;
}

template <unsigned TSPAN, unsigned TWEIGHT>
void _createValueArray2(StringSet<String<Dna5> > & reads, String<Pair<uint64_t, uint64_t> > & hs, Shape<Dna5, Minimizer<TSPAN, TWEIGHT> > & shape, int step, int l)
{
    typedef String<Pair<uint64_t, uint64_t> > StringPair;

    //TShape shape;
    StringPair tmp;
    String<uint64_t> tmp3;
    uint64_t p = 0, q=0, c = 0, n = -1, pre = ~0, count = 0, mask = ((uint64_t)1 << 63);
    std::cerr << "    _createValueArray() String<Dna5> " << std::endl;
    double time = sysTime();
    resize(tmp, lengthSum(reads) - length(reads) * (shape.span - 1));
    resize(tmp3, lengthSum(reads) - length(reads) * (shape.span - 1)+1);
    for(uint64_t j = 0; j < length(reads); j++)
    {
        hashInit(shape, begin(reads[j]));
        for (uint64_t k =0; k < length(reads[j]) - shape.span + 1; k++)
        {
            if(ordValue(*(begin(reads[j]) + k + shape.span - 1)) == 4)
            {
                k += hashInit(shape, begin(reads[j]) + k);
                if(k >  length(reads[j]) - shape.span + 1)
                    break;
            }
            hashNext(shape, begin(reads[j]) + k);
            //std::cerr << k << " " << shape.hValue << std::endl;
            _setBodyNode(tmp3[p], shape.YValue, (uint64_t)1, _createSANode(j, k));
            if (pre ^ shape.XValue)
            {
                tmp[q].i1 = shape.XValue;
                tmp[q].i2 = p;
                pre = shape.XValue;
                tmp3[p] |= mask;
                q++;
            }
            p++;
        }
    }
    //back(tmp3) |= mask;
    resize(tmp, q);
    //*(end(tmp3)) |= (mask);
    back(tmp3) |= mask;

    std::cerr << "        loading time " << sysTime() - time << std::endl;
    c = tmp[0].i2;
    p = q = count = 0;
    n = -1;
    _sort3(begin(tmp), end(tmp), step, l);         // sort parameters 1
    std::cerr << "        sort xvalue time " << sysTime() - time << std::endl;
    c = tmp[0].i2;
    for (uint64_t q = 0;  q < length(hs) - 1; q++)
    {
        if (tmp3[c] & mask)
        {
            c = tmp[++n].i2;
        }
        hs[q].i1 = tmp[n].i1;
        hs[q].i2 = tmp3[c] & (~mask);
        c++;
    //std::cerr << " c " << c << " tmp " << length(tmp3) << " length " << length(hs) << " q " << q << std::endl;
    }
    std::cerr << "        xvalue expand " << sysTime() - time << std::endl;
    hs[length(hs)-1].i2 |= mask;
    pre = hs[0].i1;
    for (uint64_t k = 0; k < length(hs); k++)
    {
        if (pre ^ hs[k].i1)
        {
            pre = hs[k].i1;
            if (count < 20)                   // sort parameters
                _insertSort(begin(hs) + k - count, begin(hs) + k);
            else
                //std::sort(begin(hs) + k -count, begin(hs) + k, [](Pair<uint64_t, uint64_t> & a,
                //Pair<uint64_t, uint64_t> & b){return a.i2 > b.i2;});
    
                //std::stable_sort(begin(hs) + k -count, begin(hs) + k, comp);
                _sort3_i2_(begin(hs) + k - count, begin(hs) + k,8,8);

            count = 0;
        }
        count++;
   }


   std::cerr << "        End sort sysTime(): " <<  sysTime() - time << std::endl;
}

template <typename TObj, unsigned TSpan, unsigned TWeight>
void _qgramCountQGrams(Index<StringSet<String<TObj> >, IndexQGram<Minimizer<TSpan, TWeight>, OpenAddressing > > & index)
{
    typedef Shape<TObj, Minimizer<TSpan, TWeight> > TM_Shape;
    //typedef Iterator<String<Dna> >::Type TIter;
    typedef typename Value<TM_Shape>::Type HValue;
    //typedef std::tuple<HValue, HValue, HValue> HTuple;
    typedef Pair<uint64_t, uint64_t> PairH;
    typedef String<PairH> StringPairH;
    //typedef String<HTuple> StringTuple;
    //StringSet<DnaString> reads;

    TM_Shape shape;
    StringPairH hs, hs1;
    HValue  m = 0, sum = 0;

    double time = sysTime();

    resize(hs, lengthSum(indexText(index)) - length(indexText(index)) * (shape.span - 1) + 1);

    //std::cerr << "        _qgramCountQGrams() sysTime(): " << sysTime() - time << std::endl;
    //std::cerr << "            lengthSum(StringSet) = " << lengthSum(indexText(index)) << std::endl;

     _createValueArray2(indexText(index), hs, shape, 9, 5);
    hs[length(hs) - 1].i1 = 1;
    _setBodyNode(hs[length(hs) - 1].i2,  (HValue)0, (HValue)0, (HValue)1);//std::make_tuple((HValue)1, (HValue)0, (HValue)1);
    std::cerr << "            sort sysTime(): " << sysTime() - time << std::endl;
    HValue countx = 1, counth = 1, tmp = 0, countdh = 0, countb = 0, hk = 0;
    resize(index.sa, length(hs) - 1);
    for (HValue k = 1; k < length(hs); k++)
    {
        //if (std::get<0>(hs[k]) ^ std::get<0>(hs[k - 1])|std::get<1>(hs[k]) ^ std::get<1>(hs[k - 1]))
        if ((hs[k].i1 ^ hs[k - 1].i1 )|(_getBodyValue(hs[k].i2 ^ hs[k - 1].i2)))
        { _setBodyNode(index.dir[index.start + hk], _getBodyValue(hs[k-1].i2), _BodyType_code, tmp);
            hk++;
            countb++;
            countdh++;
            tmp = counth;
        }
        if (hs[k].i1 ^ hs[k - 1].i1)
        {
            
            if (countb < blocklimit)
            {
                requestDir(index.dir, index.start, _makeHeadNode(hs[k-1].i1), _makeEmptyNode(index.start + hk - countb));
                for (HValue j = 0; j < countb; j++)
                    _setBodyType_Begin(index.dir[index.start + hk - countb]);
            }
            else
            {
                hk -= countb;
                requestDir(index.dir, index.start, _makeVtlHeadNode(hs[k-1].i1), _makeEmptyNode(index.start + hk));
                for (HValue j = k - countx; j < k; j++)
                    if ((hs[j].i1 ^ hs[j + 1].i1) | _getBodyValue(hs[j].i2 ^ hs[j + 1].i2))
                    {
                        requestDir(index.dir, index.start, _makeHVlHeadNode(xy2h(shape, hs[j].i1, _getBodyValue(hs[j].i2))), _makeEmptyNode(index.start+hk));
                        _setBodyType_Begin(index.dir[index.start + hk]);
                        hk++;
                    }
            }
            countb = 0;
            countx = 1;
        }
        else
        {
            countx++;
        }
        
        index.sa[k - 1] = _getBodyCounth(hs[k-1].i2);
        counth++;
    }
    //std::cerr << std::endl;
    //std::cerr << counth << std::endl;
    resize(index.dir, index.start + countdh + 10);
    _setBodyNode(index.dir[index.start + countdh], _bitEmpty, _BodyType_code, counth - 1); 
    _setBodyType_Begin(index.dir[index.start + countdh]);
    index._Empty_Dir_ = index.start + countdh + 1;
    //_setBodyNode(index.dir[index.start + countdh], _bitEmpty, _BodyTypeEnd_code, counth - 1); 
    //_setBodyType_Begin(index.dir[index.start + countdh]);
    //index._Empty_Dir_ = index.start + countdh;
    //_setBodyNode(index.dir[index.start + countdh + 1], _bitEmpty, _BodyTypeEnd_code, counth - 1); 
    //_setBodyType_Begin(index.dir[index.start + countdh + 1]);
    std::cerr << "            End _qgramCountQGrams() sysTime(): " << sysTime() - time << std::endl;
}

template <typename TObj, unsigned TSpan, unsigned TWeight>
void _createQGramIndex(Index<StringSet<String<TObj> >, IndexQGram<Minimizer<TSpan, TWeight>, OpenAddressing > >& index, StringSet<String<TObj> > & seq)
{
    double time = sysTime(); 
    //std::cerr << "    createQGramIndexDirOnly() sysTime(): " << std::endl;
    _qgramClearDir(index);
    _qgramCountQGrams(index);
    //std::cerr << "        End createQGramIndexDirOnly() sysTime(): " << sysTime() - time << std::endl;
    std::cerr << "        index.dir " << (float)length(index.dir) /1024/1024/1024 *8 << " GB index.sa " 
                << (float) length(index.sa) /1024/1024/128 << " GB" << std::endl;
    std::cerr << "    End creating index. Time[s] " << sysTime() - time << std::endl;
}

template <typename TObj, unsigned TSpan, unsigned TWeight>
void _createQGramIndex(Index<StringSet<String<TObj> >, IndexQGram<Minimizer<TSpan, TWeight>, OpenAddressing > >& index)
{
    double time = sysTime(); 
    //std::cerr << "    createQGramIndexDirOnly() sysTime(): " << std::endl;
    _qgramClearDir(index);
    _qgramCountQGrams(index);
    //std::cerr << "        End createQGramIndexDirOnly() sysTime(): " << sysTime() - time << std::endl;
    std::cerr << "        index.dir " << (float)length(index.dir) /1024/1024/1024 *8 << " GB index.sa " 
                << (float) length(index.sa) /1024/1024/128 << " GB" << std::endl;
    std::cerr << "    End creating index. Time[s] " << sysTime() - time << std::endl;
}

//========================================================
//Begin(P2):This section is to optimize 25mer for mapping

//Hs: String<uint64_t>
//types of node in Hs including: 1.head node and 2.body node
//head: Headflag[1] = 0|Pointer[23]| xvalue[40]
//body: bodyflag[1] = 1|N/A[2]|yvalue[20] |typeCode[1]|sa[40]
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
    
    HsBase(bool cerr):
        bit(XValueBit),
        bodyYBit(_BodyValue_bits),
        bodyYBitLen(20),
        bodyYMask((1ULL << bodyYBitLen) - 1),
        bodyCodeBit(_BodyType_bits),
        pointerBit(bit), 
        pointerBitLen(23),
        mask((1ULL << bit) - 1),
        pointerMask((1ULL << (pointerBitLen)) - 1),
        maxPointer(1 << pointerBitLen),
        headTypeFlag(0),
        typeFlag(1ULL << 63),
        typeFlag2(1ULL << (63 - bodyYBit)),
        typeMask(typeFlag - 1),
        bodyCodeFlag(1ULL << bodyCodeBit)

        {
            if (cerr)
                std::cerr << "HsBase::pointerBit " << pointerBit << std::endl;
        }
        
}_DefaultHsBase(false);

struct Hs
{
    typedef uint64_t ValueType; 
    typedef uint64_t ValueBodyType;
   
    bool isHead(uint64_t const &, 
                uint64_t const & = _DefaultHsBase.typeFlag);
    uint64_t MinusX(uint64_t const &, uint64_t const &, 
                    uint64_t const & = _DefaultHsBase.mask);
    void setHsHead(uint64_t &, uint64_t const &, uint64_t const &, 
                   uint64_t const & bit = _DefaultHsBase.pointerBit, 
                   uint64_t const & typeFlag = _DefaultHsBase.typeMask);
    uint64_t getHeadX(uint64_t const &, 
                      uint64_t const & = _DefaultHsBase.mask);
    uint64_t getHeadPtr(uint64_t const &, 
                        uint64_t const & = _DefaultHsBase.pointerBit, 
                        uint64_t const & = _DefaultHsBase.pointerMask);
    void setHsBody(uint64_t &, uint64_t const &,  uint64_t const & id, uint64_t const & pos,
                   uint64_t const & typeFlag = _DefaultHsBase.typeFlag
                  );
    uint64_t getHsBodyY(uint64_t const &,
        uint64_t const & = _DefaultHsBase.bodyYBit, 
        uint64_t const & = _DefaultHsBase.bodyYMask
    );

    uint64_t getHsBodyS(uint64_t const & val, 
        uint64_t const & mask = _DefaultHsBase.mask
    )
    {return val & mask;}
    void setHsHeadPtr(uint64_t &, uint64_t const &, 
                      uint64_t const & = _DefaultHsBase.bit, 
                      uint64_t const & = _DefaultHsBase.mask);
    bool isBody(uint64_t const & val, uint64_t const & flag = _DefaultHsBase.typeFlag)
    {return val & flag;}
    bool isBodyYEqual(uint64_t const & hval, uint64_t const & yval, 
                    uint64_t const & bit = _DefaultHsBase.bodyYBit,
                    uint64_t const & flag = _DefaultHsBase.typeFlag2
    )
    //return if hval is body and if yvalue of hval euqals to yval
        {return ((hval >> bit) ^ yval) == flag;}
    
}_DefaultHs;

//XNode = struct v1[64],v2[32]: .v1: hashvalue; .v2:pointer to YNode
//.v1: v2 type[2]|value[60]|v1 type[2]
//Types of .v1 including: 1.empty node{00} 2.xvalue head{01} 3.hvalue{10} head 4 virtual head{11}
//  0 empty
//  1 xvalue head
//  2 head 10
//  3 virtual head 


struct XNodeBase   //define dirNode
{
    typedef uint64_t NodeType;
    typedef uint64_t ReturnType;
    typedef uint64_t Bit;
    typedef uint64_t Mask;
    
    Bit bit;
    Mask mask;
    Bit bit2;
    Mask mask2;
    
    NodeType emptyNode;
    NodeType xHead;
    NodeType xHead1;
    NodeType virtualHead;
    
    ReturnType returnDir;
    ReturnType returnSa;
    
   
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

struct XNode
{
    typedef uint64_t TypeV1;
    typedef unsigned TypeV2;
    typedef uint64_t TypeV2L;
    
    TypeV1 val1;
    TypeV2 val2;
};

struct XNodeFunc
{
    uint64_t getAddY();
    uint64_t hash(uint64_t const &);
    void setXNode(XNode &, XNode::TypeV1 const & val1, XNode::TypeV2 const &, 
                XNodeBase::NodeType const &, XNodeBase::ReturnType const &, 
                XNodeBase::Bit const & = _DefaultXNodeBase.bit,
                XNodeBase::Bit const & = _DefaultXNodeBase.bit2
                 );
    XNode::TypeV1 makeYXKey(typename Hs::ValueBodyType const &, XNode::TypeV1 const &, 
                            XNodeBase::Mask const & = _DefaultHsBase.bodyYMask << _DefaultHsBase.bodyYBit, 
                            XNodeBase::Mask const & = _DefaultHsBase.mask);
    XNode::TypeV1 collision(XNode::TypeV1 const &, XNode::TypeV1 const &, XNodeBase::Mask const & = _DefaultXNodeBase.mask2);
    XNode::TypeV2L makeReturnVal(XNode const &, XNodeBase::Mask const & = _DefaultXNodeBase.mask2);
    
}_DefaultXNodeFunc;



struct XString
{
    String<XNode> xstring;
    uint64_t mask;
    
    XString(){};
    XString(uint64_t const & seqlen);
    uint64_t _fullSize(uint64_t const & seqlen, float const & alpha = 1.6);
};

template <unsigned TSPAN>
struct HIndexBase
{
    typedef StringSet<String<Dna5> > Text;
    typedef String<typename Hs::ValueType> YSA;
    typedef XString XStr;
    typedef Shape<Dna5, Minimizer<TSPAN> > TShape;
    
    static const double defaultAlpha;
};

template <unsigned TSPAN> 
const double HIndexBase<TSPAN>::defaultAlpha(1.6);
 
template <unsigned TSPAN>
class HIndex
{
    
    public:
        typedef typename HIndexBase<TSPAN>::TShape TShape;
        typename HIndexBase<TSPAN>::YSA             ysa;        
        typename HIndexBase<TSPAN>::XStr            xstr;       
        typename HIndexBase<TSPAN>::TShape          shape;
        double   alpha;    
        uint64_t emptyDir;
        
        HIndex():
            alpha(HIndexBase<TSPAN>::defaultAlpha) 
            {}
        HIndex(typename HIndexBase<TSPAN>::Text const & text):
            alpha(HIndexBase<TSPAN>::defaultAlpha) 
        {
            (void) text;
        }
        
};


XString::XString(uint64_t const & seqlen)
{
    std::cout << "1";
    _fullSize(seqlen);
    std::cout << "2";
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


inline bool Hs::isHead(uint64_t const & val, uint64_t const & flag)
{
    return (val & flag) ^ flag;
}

inline void Hs::setHsHead(uint64_t & head, uint64_t const & ptr, uint64_t const & xval, uint64_t const & bit, uint64_t const & mask)
{
    head = ((ptr << bit) + xval) & mask;
}

inline uint64_t Hs::MinusX(uint64_t const & value1, uint64_t const & value2, uint64_t const & mask)
{
    return ((value1 - value2) & mask);
}

inline uint64_t Hs::getHeadX(uint64_t const & value, uint64_t const & mask)
{
    return value & mask;
}  

inline uint64_t Hs::getHeadPtr(uint64_t const & val, uint64_t const & bit, uint64_t const & mask)
{
    return (val >> bit) & mask;
}

inline void Hs::setHsBody(uint64_t & val, uint64_t const & yval, uint64_t const & id, uint64_t const & pos, uint64_t const & typeFlag)
{
    val = ((yval << _BodyValue_bits)|typeFlag) + (id << _BaseNum_bits) + (pos);
    
}

inline uint64_t Hs::getHsBodyY(uint64_t const & val, uint64_t const & bit, uint64_t const & mask)
{
    return (val >> bit) & mask;
}

inline void Hs::setHsHeadPtr(uint64_t & val, uint64_t const & ptr,  uint64_t const & bit, uint64_t const & mask)
{
    val = (val & mask) + (ptr << bit);
}
/*
 * bucket[]-
 */
/*
template <typename TIt>
inline bool _hsSortX(TIt const & begin, TIt const & end, unsigned const & xValBitLen)
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
    uint64_t count[512];
    //int count[1024];
    
    String<uint64_t> output;
    resize(output, end - begin);
    //std::cerr << "end - begin " <<end - begin << std::endl;
    for (uint64_t j = 0; j < l; j++)
    {
        l_move -= p_bit;
        for (int k = 0; k< (1<<p_bit); k++)
            count[k]=0;
        for (int64_t k = 0; k < end - begin; k += _DefaultHs.getHeadPtr(*(begin + k)))
        {
            count[*(begin + k) << l_move >> r_move] += _DefaultHs.getHeadPtr(*(begin + k));
        }
        for (int k = 1; k < (1 << p_bit); k++)
        {
            count[k] += count[k - 1];
        }
        for (int64_t k = end - begin - 1;  k >=0; k--)
        {
        
            if (_DefaultHs.isHead(*(begin + k)))
            {
                uint64_t x = *(begin + k) << l_move >> r_move;
                uint64_t ptr = _DefaultHs.getHeadPtr(*(begin + k));
                count[x] -= ptr;
                for (uint64_t it = 0; it < ptr; it++)
                {
                    output[count[x] + it] = *(begin + k + it);
                }
            }
            
        }
        for (int64_t k = 0; k < end - begin; k++)
            *(begin + k) = output[k];
    }
    return true;
}
*/

/*
 * parallel sort hs
 * bucket[]+
 */
static unsigned const maxThread = 4;
static unsigned const maxBucket = 513;

template <typename TIt>
inline bool _hsSortX(TIt const & begin, TIt const & end, unsigned const & xValBitLen)
{
    if (xValBitLen <34 || xValBitLen > 42)
    {
        std::cerr << "[Error]: _dirSortX " << xValBitLen << "\n";
        return false;
    }
    unsigned const bit[18] = {9,4,9,4,9,4,8,5,8,5,8,5,8,5,7,6,7,6}; //xValueBitLen 34 - 42;
    unsigned const p_bit = bit[(xValBitLen - 34) << 1];
    unsigned const l =  bit[((xValBitLen - 34) << 1) + 1];
    unsigned const r_move = 64 - p_bit;
    unsigned l_move = 64;
    unsigned threads = 4;//omp_get_num_threads();
    if (threads > maxThread)
        threads = maxThread;
    uint64_t const mask = (1 << p_bit) - 1;
    uint64_t size = (end - begin) / threads;
    unsigned thd1 = end - begin - size * threads;
    uint64_t thd_n1 = (size + 1) * thd1;
    omp_set_num_threads(threads);
    std::cerr << "[hssort] " << threads << "\n";
    std::vector<std::vector<uint64_t> > ctd(threads, std::vector<uint64_t>((1<<p_bit) + 1, 0));
    std::vector<std::vector<std::vector<uint64_t> > > next(threads, 
                    std::vector<std::vector<uint64_t> >(threads, std::vector<uint64_t>((1 << p_bit) + 1, 0)));
    String<uint64_t> output;
    resize(output, end - begin);
    
    //Initialize ctd[][] 
    
    //std::cerr << "[hssort2] " << threads << "\n";
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
                //printf("[originHs] %d %d %d %d %d\n", x, ptr, k, size ,end - begin);
            }
        }
    }

    //std::cerr << "[hssort3] " << threads << "\n";
    unsigned PBit = 1 << p_bit;
    for (uint64_t j = 0; j < l; j++)
    {
        //std::cerr << "[j] " << j << "\n";
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
                    unsigned thd_num = (ctd[thd_id][x] < thd_n1)?ctd[thd_id][x]/(size + 1):(ctd[thd_id][x] - thd_n1) / size + thd1;
                    //unsigned thd_num = (ctd[x][thd_id] < thd_n1)?ctd[x][thd_id]/(size + 1):(ctd[x][thd_id] - thd_n1) / size + thd1;
                    //printf("[thd_num] %d %d %d\n", thd_id, x, ctd[thd_id][x]);
                    ctd[thd_id][x] += ptr;
                    //ctd[x][thd_id] += ptr;
                    x = *(begin + k) << (l_move - p_bit)>> r_move;
                    if (thd_num == threads - 1)
                        next[thd_id][0][x + 1] += ptr;
                    else
                        next[thd_id][thd_num + 1][x] += ptr;
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
            
           //     for (unsigned n = 0; n < (1 << p_bit); n++)
           //     {
           //         //std::fill(ctd[k].begin(), ctd[k].end(), 0);
           //         #pragma omp parallel for
           //         for (unsigned k=0; k < threads; k++)
           //         {
           //             ctd[n][k] = 0;
           //             for (unsigned m = 0; m < threads; m++)
           //             {
           //                 //ctd[k][n] += next[m][k][n];
           //                 ctd[n][k] += next[m][k][n];
           //                 next[m][k][n] = 0;
           //             }   
           //         }
           //             
           //     }
        
        
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
 * serial sort hs
 * bucket[]+
 */
/*
template <typename TIt>
inline bool _hsSortX(TIt const & begin, TIt const & end, unsigned const & xValBitLen)
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
 */

template <typename TIter>//, typename Comp>
void insertSort(TIter const & begin, TIter const & end)//, Comp const & comp)
{
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

template <typename TIt>
inline bool _sort_YSA_Block(TIt const & begin, TIt const & end, unsigned const & sortModeThr = 60) // sort y and sa
{
    typedef typename Value<TIt>::Type ValueType;
    if (end - begin< sortModeThr)
        insertSort(begin, end);
    else
        std::sort(begin, end, std::greater<ValueType>());
    return true;
}

template <typename TIt>
inline bool _hsSortY_SA(TIt const & begin, TIt const & end) // sort y and sa
{
    typedef typename Value<TIt>::Type ValueType;
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

template <typename TIter>
inline void _hsSort(TIter const & begin, TIter const & end, unsigned const & shapeWeight)
{
    std::cerr << "      sorting xstr \n";
    double time = sysTime();
    //uint64_t mask = (1ULL << 9) - 1;
    _hsSortX(begin, end, shapeWeight << 1);
    //std::sort(begin, end);
    std::cerr << "      _dirSortX Time[s]" << sysTime() - time << std::endl;
    //_hsSortY_SA(begin, end);
    //std::cerr << " _dirSortY " << sysTime() - time << "\n";
}

/*
 * serial creat hash array
 */
template <unsigned SHAPELEN>
bool _createHsArray(StringSet<String<Dna5> > const & seq, String<uint64_t> & hs, Shape<Dna5, Minimizer<SHAPELEN> > & shape)
{
    double time = sysTime();
    uint64_t preX = ~0;
    //-k
    //int64_t ptr = 0, count = 1;
    int64_t ptr = 0, count = -1;
    resize(hs, lengthSum(seq) << 1);
    //-k
    //unsigned start = 2;
    hs[0] = hs[1] = 0;
    for(uint64_t j = 0; j < length(seq); j++)
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
            if (k % 3 != 0)
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
    shrinkToFit(hs);
    _hsSort(begin(hs), begin(hs) + count, shape.weight);
    //_hsSort(begin(hs) + start, begin(hs) + count, shape.weight);
    
    std::cerr << "      End createHsArray " << std::endl;
    return true;
}

/*
 * parallel creat hash array
 */
/*
template <unsigned SHAPELEN>
bool _createHsArray(StringSet<String<Dna5> > const & seq, String<uint64_t> & hs, Shape<Dna5, Minimizer<SHAPELEN> > & shape)
{
    double time = sysTime();
    int64_t sum = 0;
    unsigned const step = 3;
    omp_set_num_threads(4);
    for(uint64_t j = 0; j < length(seq); j++)
    {
        #pragma omp parallel reduction(+: sum)
        {
            uint64_t preX = ~0;
            int64_t ptr = 0, count = -1;
            uint64_t k;
            unsigned size2 = length(seq[j]) / omp_get_num_threads();
            unsigned ChunkSize = size2;
            if (omp_get_thread_num() < length(seq[j]) - size2 * omp_get_num_threads())
            {
                ChunkSize = size2 + 1;
            }
            printf("createhs %d %d\n", ChunkSize, omp_get_thread_num());
            Shape<Dna5, Minimizer<SHAPELEN> > tshape = shape; 
            String<uint64_t> tmpHs;
            resize(tmpHs, ChunkSize << 1);
            tmpHs[0] = tmpHs[1] = 0;
            hashInit(tshape, begin(seq[j]));
            #pragma omp for
            for (k = 0; k < length(seq[j]) - tshape.span + 1; k++)
            {
                if(ordValue(*(begin(seq[j]) + k + tshape.span - 1)) == 4)
                {
                    k += hashInit(tshape, begin(seq[j]) + k);
                    if (k > ChunkSize - tshape.span + 1)
                    {
                        k = ChunkSize - ChunkSize % step + step;
                    }
                }
                    
                hashNext(tshape, begin(seq[j]) + k);
                if (k % step != 0)
                {
                    if (tshape.XValue ^ preX)
                    {
                        _DefaultHs.setHsHead(tmpHs[++count - ptr], ptr, preX);
                        ptr = 2;
                        preX = tshape.XValue; 
                    }
                    else
                    {
                        ++ptr;
                    }
                    _DefaultHs.setHsBody(tmpHs[++count], tshape.YValue, j, k); 
                    //std::cout << "hs[count] " << (hs[count] & ((1ULL << 30) - 1))<< std::endl;
                }
            }
            _DefaultHs.setHsHead(tmpHs[++count - ptr], ptr, tshape.XValue);
            resize(tmpHs, count);
            sum += count;
            #pragma omp for ordered
            for (unsigned k = 0; k < omp_get_num_threads(); k++)
                #pragma omp ordered
                {
                    append (hs, tmpHs);
                }
        }
    }
    _DefaultHs.setHsHead(hs[sum], 0, 0);
    std::cerr << "      init Time[s]" << sysTime() - time << " " << std::endl;
//-k
    resize(hs, sum + 1);
    shrinkToFit(hs);
    _hsSort(begin(hs), begin(hs) + sum, shape.weight);
    //_hsSort(begin(hs) + start, begin(hs) + count, shape.weight);
    
    std::cerr << "      End createHsArray " << std::endl;
    return true;
}
*/
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

inline uint64_t XNodeFunc::hash(uint64_t const & val)
{
    uint64_t key = val;
    key = (~key) + (key << 21); // key = (key << 21) - key - 1;
    key = key ^ (key >> 24);
    key = (key + (key << 3)) + (key << 8); // key * 265
    key = key ^ (key >> 14);
    key = (key + (key << 2)) + (key << 4); // key * 21
    key = key ^ (key >> 28);
    key = key + (key << 31);
    return key;        
}

inline void XNodeFunc::setXNode(XNode & val, XNode::TypeV1 const & val1, XNode::TypeV2 const &val2, 
                            XNodeBase::NodeType const & t, XNodeBase::ReturnType const & r, 
                            XNodeBase::Bit const & bit, XNodeBase::Bit const & bit2)
{
    val.val1 = (val1 << bit) + t + (r << bit2);
    val.val2 = val2;
}

inline XNode::TypeV1 XNodeFunc::makeYXKey(typename Hs::ValueBodyType const & hval, XNode::TypeV1 const & xval, 
                            XNodeBase::Mask const & masky, XNodeBase::Mask const & maskx)
{
    return (hval & masky) + (xval & maskx);
}

inline XNode::TypeV1 XNodeFunc::collision(XNode::TypeV1 const & val, XNode::TypeV1 const & xnode_val1, XNodeBase::Mask const & mask)
{
    return (val ^ xnode_val1) & mask;
}

inline XNode::TypeV2L XNodeFunc::makeReturnVal(XNode const & xnode, XNodeBase::Mask const & mask)
{
    return (xnode.val1 & (~mask)) + xnode.val2;
    //return xnode.val2;
}

inline uint64_t requestXNode_noCollision (XString & xstr, uint64_t const & xval, 
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

inline uint64_t requestXNode_noCollision_Atomic (XString & xstr, uint64_t const & xval, 
                uint64_t const & val2,  uint64_t const & nodeType, uint64_t const returnType)
{
    uint64_t h1 = _DefaultXNodeFunc.hash(xval) & xstr.mask;
    uint64_t delta = 0;
    //while (xstr.xstring[h1].val1) //0 stands for empty
    while (atomicCas<uint64_t>(xstr.xstring[h1].val1, 0, 1))
    {
        h1 = (h1 + delta +1) & xstr.mask;
        delta++;
    }
    _DefaultXNodeFunc.setXNode(xstr.xstring[h1], xval, val2, nodeType, returnType);
    return h1;
}

inline uint64_t getXDir(XString const & xstr, uint64_t const & xval, uint64_t const & yval, bool & vflag)
{
        uint64_t val, delta = 0;
        uint64_t h1 = _DefaultXNodeFunc.hash(xval) & xstr.mask;
        
        //_setHeadNode(val, xval);
//!!!!! need to modify;
        val = (xval << 2) + _DefaultXNodeBase.xHead;
        vflag = false;
        while (xstr.xstring[h1].val1)
        {
            //switch (xstr.xstring[h1].val1 ^ val) 
            switch(_DefaultXNodeFunc.collision(xstr.xstring[h1].val1, val))
            {
                case 0:
                    //std::cerr << "case1\n";
                    //return _DefaultXNodeFunc.makeReturnVal(xstr.xstring[h1]);
                    return xstr.xstring[h1].val2;
                case 2:
//!!!!! need to modify;
                    vflag = true;
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
        //std::cout << xval << " " << yval << " error \n";
        //return _DefaultXNodeBase._Empty_Dir_;
    return 0;
}


template <unsigned span>
inline uint64_t getXDir(HIndex<span> const & index, uint64_t const & xval, uint64_t const & yval, bool & vflag)
{
    uint64_t val, delta = 0;
    uint64_t h1 = _DefaultXNodeFunc.hash(xval) & index.xstr.mask;
    
    //_setHeadNode(val, val);
//!!!!! need to modify;
    val = (xval << 2) + _DefaultXNodeBase.xHead;
    vflag = false;
    while (index.xstr.xstring[h1].val1)
    {
        //switch (index.xstr.xstring[h1].val1 ^ val) 
        switch(_DefaultXNodeFunc.collision(index.xstr.xstring[h1].val1, val))
        {
            case 0:
                //std::cerr << "case1\n";
                //return _DefaultXNodeFunc.makeReturnVal(index.xstr.xstring[h1]);
                return index.xstr.xstring[h1].val2;
            case 2:
//!!!!! need to modify;
                vflag = true;
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

template <unsigned span>
inline uint64_t getXYDir(HIndex<span> const & index, uint64_t const & xval, uint64_t const & yval)
{
    uint64_t val, delta = 0;
    uint64_t h1 = _DefaultXNodeFunc.hash(xval) & index.xstr.mask;
    uint64_t pos;
    
    //_setHeadNode(val, xval);
//!!!!! need to modify;
    val = (xval << 2) + _DefaultXNodeBase.xHead;
    while (index.xstr.xstring[h1].val1)
    {
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

/*
template <unsigned span>
inline uint64_t getNextXYDir(HIndex<span> const & index, HShape<span> const & shape)
{
    if (shape.preX ^ shape.XValue)
    {
        pos = getXYDir(index, shape.XValue, shape.YValue, shape.vflag);
        shape.preX = shape.XValue;
    }
    else 
    {
        if(shape.vflag)
        {
            pos = getXYDir(index, shape.XValue + (shape.YValue << 40), 0, shape.vflag);
            shape.vflag = true;
        }
        else 
            do{
                if (shape.YValue == _DefaultHs.getHsBodyY(index.ystr[pos]))
                    return pos;
                if (shapeYValue > _DefaultHs.getHsBodyY(index.ystr[pos]))
                    return index.emptyDir;
            }while(_DefaultHs.isBody(index.ystr[++pos]));
    }
    return
}
*/

inline uint64_t gethDir(XString const & xstr, uint64_t const & hval)
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
/*
template <unsigned TSPAN, unsigned TWEIGHT>
bool _createYSA(String<uint64_t> & hs, XString & xstr, uint64_t & indexEmptyDir)
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
        if (ptr < blocklimit)
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
        if (ptr < blocklimit)
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

*/
/*
 * parallel sort ysa
 */

template <unsigned TSPAN, unsigned TWEIGHT>
bool _createYSA(String<uint64_t> & hs, XString & xstr, uint64_t & indexEmptyDir)
{
//-k
    omp_set_num_threads(4);
    uint64_t k = _DefaultHs.getHeadPtr(hs[0]);
    uint64_t preX = _DefaultHs.getHeadX(hs[0]);
    //uint64_t k = 2 + _DefaultHs.getHeadPtr(hs[2]);
    //uint64_t preX = _DefaultHs.getHeadX(hs[2]);
    uint64_t ptr = k;
    uint64_t block_size = ptr;
    uint64_t countMove = 0, prek = 0;
    double time = sysTime();
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
#pragma omp parallel
{
    uint64_t ptr = 0;
    #pragma omp for
    for (uint64_t k = 0; k < length(hs) - 2; k++)
    { 
//        printf("[sorty] %d\n", k);
        if(_DefaultHs.isHead(hs[k]))
        {
            ptr = _DefaultHs.getHeadPtr(hs[k]);
            _sort_YSA_Block(begin(hs) + k + 1, begin(hs) + k + ptr);
            k += ptr - 1;
        }   
    }
    
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
        if (ptr < blocklimit)
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
#pragma omp parallel 
{
    uint64_t ptr = 0;
    #pragma omp for
    for (k =0; k < length(hs); ++k)
    {
        //printf("[requestdir atomic] %d %d %d\n", omp_get_thread_num(), length(hs),k);
        if (_DefaultHs.isHead(hs[k]) && _DefaultHs.getHeadPtr(hs[k]))
        {
            ptr = _DefaultHs.getHeadPtr(hs[k]);
            if (ptr < blocklimit)
           {
               requestXNode_noCollision_Atomic(xstr, _DefaultHs.getHeadX(hs[k]), 
                       k + 1, _DefaultXNodeBase.xHead, _DefaultXNodeBase.returnDir);
           }
            else
            {
                uint64_t xval = _DefaultHs.getHeadX(hs[k]);
                requestXNode_noCollision_Atomic(xstr, xval, 
                        ~1, _DefaultXNodeBase.virtualHead, _DefaultXNodeBase.returnDir);
                for (unsigned j = k + 1; j < k + ptr; j++)
                {
                    if(_DefaultHs.getHsBodyY(hs[j] ^ hs[j - 1]))
                    {
                        requestXNode_noCollision_Atomic(xstr, (xval + ((hs[j] & ((1ULL<<61) - (1ULL<<41))) >>1)), 
                                j, _DefaultXNodeBase.xHead, _DefaultXNodeBase.returnDir);
                    }
                    
                }
            }
            k += ptr - 1;    
        }
    }

}

    std::cerr << "request dir " << sysTime() - time << std::endl;
    return true;
}


template <unsigned SHAPELEN>
bool _createQGramIndexDirSA(StringSet<String<Dna5> > const & seq, XString & xstr, String<uint64_t> & hs,  Shape<Dna5, Minimizer<SHAPELEN> > & shape, uint64_t & indexEmptyDir, bool Efficient)    
{
    
    typedef Shape<Dna5, Minimizer<SHAPELEN> > ShapeType;
    double time = sysTime();
    if (Efficient)
        _createHsArray(seq, hs, shape);
    //time = sysTime();
    _createYSA<LENGTH<ShapeType>::VALUE, WGHT<ShapeType>::VALUE>(hs, xstr, indexEmptyDir);
    std::cerr << "  End creating Index Time[s]:" << sysTime() - time << " \n";
    return true; 
    
}
 
template <typename TDna, unsigned span>
bool createHIndex(StringSet<String<TDna> > & seq, HIndex<span> & index)
{
    return _createQGramIndexDirSA(seq, index.xstr, index.ysa, index.shape, index.emptyDir, true);
}

template <typename TDna, unsigned TSpan>
bool _createQGramIndex(HIndex<TSpan> & index, StringSet<String<TDna> > & seq)
{
    return _createQGramIndexDirSA(seq, index.xstr, index.ysa, index.shape, index.emptyDir, true);
}

/*
template <unsigned span, typename Result>
bool streamSeq(HIndex<span> & index, StringSet<String<Dna5> > & seqs, Result & result, void (*call_back)(uint64_t const & , Result & ))    
{
    uint64_t preX = 0, pos = 0;
    bool vflag;
    
        for(uint64_t j = 0; j < length(seqs); j++)
        {
            hashInit(index.shape, begin(seqs[j]));
            for (uint64_t k =0; k < length(seqs[j]) - index.shape.span + 1; k++)
            {
                if(ordValue(*(begin(seqs[j]) + k + index.shape.span - 1)) == 4)
                {
                    k += hashInit(index.shape, begin(seqs[j]) + k);
                    
                }
                hashNext(index.shape, begin(seqs[j]) + k);
                
                if (preX ^ index.shape.XValue)
                {
                    pos = getXDir(index.xstr, index.shape.XValue, index.shape.YValue, vflag);
                    preX = index.shape.XValue;
                }
                else 
                {
                    if(vflag)
                    {
                //        pos = getXDir(index.xstr, index.shape.XValue, index.shape.YValue, vflag);
                        pos = getXDir(index.xstr, index.shape.XValue + (index.shape.YValue << 40), 0, vflag);
                        vflag = true;
                    }
                }
                uint64_t m = pos;
                do
                {
                    if(_DefaultHs.getHsBodyY(index.ysa[m]) == index.shape.YValue)
                    {
                        //sum ^= _DefaultHs.getHsBodyS(index.ysa[m]);
                        call_back(_DefaultHs.getHsBodyS(index.ysa[m]), result);
                        continue;
                    }
                    if (_DefaultHs.getHsBodyY(index.ysa[m]) < index.shape.YValue)
                        break;
                }
                while (_DefaultHs.isBody(index.ysa[++m]));
            }
        }
}


inline void streamCall_back(uint64_t const & sa,  uint64_t & result)
{
    result ^= sa;
}
*/
//End(P2)
//=========================================================================


/*
template <unsigned TSpan, unsigned TWeight>
void createQGramIndexDirOnly2(Index<StringSet<DnaString>, IndexQGram<Minimizer<TSpan, TWeight>, OpenAddressing > >& index)
{
    double time = sysTime(); 
    std::cerr << "    createQGramIndexDirOnly() sysTime(): " << std::endl;
    _qgramClearDir(index);
    _qgramCountQGrams2(index);
    std::cerr << "        index.dir[index._Empty_Dir_] = " << index.dir[index._Empty_Dir_] << std::endl << "        length(index.dir) = " << length(index.dir) << std::endl;
    std::cerr << "        End createQGramIndexDirOnly() sysTime(): " << sysTime() - time << std::endl;
}
*/
// ----------------------------------------------------------------------------
// Function open()
// ----------------------------------------------------------------------------
}
#endif //#ifndef SEQAN_HEADER_...

