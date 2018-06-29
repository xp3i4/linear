// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2015, Knut Reinert, FU Berlin
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
// Author: David Weese <david.weese@fu-berlin.de>
// ==========================================================================

#ifndef SEQAN_HEADER_INDEX_QGRAM_OPENADRESSING_H
#define SEQAN_HEADER_INDEX_QGRAM_OPENADRESSING_H
namespace SEQAN_NAMESPACE_MAIN
{

    struct OpenAddressing_;
    typedef Tag<OpenAddressing_> OpenAddressing;

    template <typename THashValue>
    struct BucketMap
    {
        static const THashValue EMPTY;
        String<THashValue> qgramCode;
    };

    template <typename THashValue>
    const THashValue BucketMap<THashValue>::EMPTY = (THashValue)-1;

    // use the index value type as shape value type
    template < typename TObject, typename TShapeSpec >
    struct Fibre< Index<TObject, IndexQGram<TShapeSpec, OpenAddressing> >, FibreBucketMap>
    {
        typedef typename Fibre< Index<TObject, IndexQGram<TShapeSpec, OpenAddressing> >, FibreShape>::Type TShape;
        typedef typename Value<TShape>::Type    THashValue;
        typedef BucketMap<THashValue>            Type;
    };

/*!
 * @class OpenAddressingQGramIndex
 * @extends IndexQGram
 * @headerfile <seqan/index.h>
 * @brief A <i>q</i>-gram that uses open addressing hashing instead of an array.
 *
 * @signature template <typename TIndex, typename TShapeSpec>
 *            class Index<TText, IndexQGram<TShapeSpec, OpenAddressing> >;
 *
 * @tparam TText      The @link TextConcept text type @endlink.
 * @tparam TShapeSpec The @link Shape @endlink specialization type.
 *
 * This index uses a non-trivial hashing for mapping q-gram hash values to buckets.  This reduces the sizes of bucket
 * directories (QGramDir, QGramCountsDir fibres) from &Sigma;<i><sup>q</sup></i> to min(<i>&alpha; &middot; n</i>,
 * \Sigma<i><sup>q</sup></i>), for a load factor <i>&alpha; &gt; 1</i>.  A bucket still stores occurrences (or counts)
 * of the same <i>q</i>-gram, but in contrast to the @link IndexQGram @endlink index, buckets are in random order due to
 * the hashing.
 *
 * @var double OpenAddressingQGramIndex::alpha
 * @brief Load factor.  Controls space/time-tradeoff and must be greater 1.  Default value is 1.6.
 */
#ifdef PLATFORM_WINDOWS_VS
#pragma warning( push )
// Disable warning C4521 locally (multiple copy constructors).
#pragma warning( disable: 4521 )
#endif  // PLATFORM_WINDOWS_VS

    template < typename TObject, typename TShapeSpec >
    class Index<TObject, IndexQGram<TShapeSpec, OpenAddressing> >
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
        TSA                sa;            // suffix array sorted by the first q chars
        TDir            dir;        // bucket directory
        TCounts            counts;        // counts each q-gram per sequence
        TCountsDir        countsDir;    // directory for count buckets
        TShape            shape;        // underlying shape
        TCargo            cargo;        // user-defined cargo
        TBucketMap        bucketMap;    // bucketMap table (used by open-addressing index)
        TSize            stepSize;    // store every <stepSize>'th q-gram in the index

        double            alpha;        // for m entries the hash map has at least size alpha*m

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
#ifdef PLATFORM_WINDOWS_VS
// Enable warning C4521 again (multiple copy operators).
#pragma warning( pop )
#endif  // PLATFORM_WINDOWS_VS


//x-begin: min shape open index 


 template < typename TObject, unsigned TSPAN, unsigned TWEIGHT>
    class Index<TObject, IndexQGram<MinimizerShape<TSPAN, TWEIGHT>, OpenAddressing> >
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
    // B := YValue[26]|BodyType[1]: |counth[37]
    // occ = counth[n+1] - count[n] = (B[n+1] - B[n]) & bit, bit = 00..0011...11
    static const unsigned _BodyValue_bits = 38;
    static const unsigned _BodyType_bits = 37;
    //static const unsigned _BodyValue_bits = 2;
    static const uint64_t _BodyType_code = 1;
    static const uint64_t _BodyTypeEnd_code = 0;
    static const uint64_t _BodyType_key = ~((uint64_t)1 << _BodyType_bits);
    static const uint64_t _getBody = ((uint64_t)1 << _BodyType_bits) - 1;

    static const uint64_t _bitEmptyType = 0;

    static const uint64_t _bitCode = 3;                                 // node(value) & _bitCode to acquire the type of the node 
    static const uint64_t _bitValue2 = ((uint64_t)1 << 22) - 1;
    static const uint64_t _bitEmptyCode = 0;
    static const uint64_t _bitBodyCode = 1;
    static const uint64_t _bitHeadCode = 2;
    static const uint64_t _bitVtlHeadCode = 3;
  
    //SA node:= seq num i1[40]| base num i2[16]
    static const uint64_t _BaseNum_bits = 32;    
    static const uint64_t _SeqNum_bits = 30;    
    static const uint64_t _BaseNum_code = ((uint64_t)1 << _BaseNum_bits) - 1;
 
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
    inline HValue _getDirStart(Index<StringSet<DnaString>, IndexQGram<MinimizerShape<TSPAN, TWEIGHT>, OpenAddressing> >  & index)
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
    inline HValue _getSA_i2(HValue const & node)
    {
        return node & _BaseNum_code;
    }
    
//x-end: min shape open index 


    template < typename TObject, typename TShapeSpec >
    const double Index<TObject, IndexQGram<TShapeSpec, OpenAddressing> >::defaultAlpha = 1.6;
//x-begin
    template < typename TObject,unsigned TSPAN, unsigned TWEIGHT>
    const double Index<TObject, IndexQGram<MinimizerShape<TSPAN, TWEIGHT>, OpenAddressing> >::defaultAlpha = 1.6;
//x-end
    //////////////////////////////////////////////////////////////////////////////
    // Counting sort - Step 1: Clear directory
    template < typename TDir, typename THashValue, typename TParallelTag >
    inline void _qgramClearDir(TDir &dir, BucketMap<THashValue> &bucketMap, Tag<TParallelTag> parallelTag)
    {
        typedef BucketMap<THashValue> TBucketMap;
        if (!empty(dir))
            arrayFill(begin(dir, Standard()), end(dir, Standard()), 0, parallelTag);
        if (!empty(bucketMap.qgramCode))
            arrayFill(begin(bucketMap.qgramCode, Standard()), end(bucketMap.qgramCode, Standard()), TBucketMap::EMPTY, parallelTag);
    }

    template < typename TBucketMap, typename TValue >
    inline TValue
    _hashFunction(TBucketMap const &, TValue val)
    {
        // WARNING:
        // As SSE4.2 is not always available, the bucket order is platform dependent
#ifdef __SSE4_2__
        return _mm_crc32_u64(0ul, val);
#else
        return ((val * 43) ^ (val >> 20)) + val;
#endif
    }
    
    //x1-begin
    //template <typename TIndex, typename THashValue, typename TParallelTag>
/*
    template<typename TDir, typename THashValue>
    inline THashValue
    requestBucket(TDir & dir, THashValue code, unsigned len)//, Tag<TParallelTag> parallelTag)
    {
        //std::cout << code << " " << len << std::endl;
        typedef unsigned long TSize;
        // get size of the index

        // check whether bucket map is disabled and
        // where the hash should be found if no collision took place before
        TSize hlen = length(dir);
        if (hlen == 0ul) return code;
        
        TSize h1 = _hashFunction(dir, code >> _bitValue_END);
#ifdef SEQAN_OPENADDRESSING_COMPACT
        --hlen;
        h1 %= hlen;
#else
        hlen -= 2;
        h1 &= hlen;
#endif
        // was the entry empty or occupied by our code?
        
        // if not we have a collision -> probe for our code or an empty entry
        //
        // do linear probing if we need to save memory (when SEQAN_OPENADDRESSING_COMPACT is defined)
        // otherwise do quadratic probing to avoid clustering (Cormen 1998)
        TSize delta = 0;
        (void)delta;
        //std::cout << "done\n";
        unsigned count = 0;
        //return h1;
        for (auto k = 0; k < len; k++)
        {
            //if (code == dir[h1 + k])
            //{
            //    std::cout << count<< std::endl;
            //    return h1;
            //}
            //std::cout << (dir[h1+k] & _bitCode) << std::endl;
            switch (dir[h1+k] & _bitCode){
                case 1:
                    h1 = (0 + h1 + (dir[h1 + k] >> _bitLength_END) + delta) &  hlen ;
                    //std::cout << len << " " << (dir[h1+k] >> _bitLength_END) << std::endl;
                    std::cout << k << std::endl;
                    ++delta; k=0;
                    count++;
                    break;
                case 2:  
                    h1 = (0 + h1 + (dir[h1 + k + 1] >> _bitLength_END) + delta) &  hlen ;
                    std::cout << k << std::endl;
                    //std::cout << len << " " << (dir[h1 + k + 1] >> _bitLength_END) << std::endl;
                    ++delta; k=0;
                    count++ ;
                    break;
                case 3:
                    std::cout << k << std::endl;
                    h1 = (h1 + 1 + delta) &  hlen ;
                    ++delta; k=0;
                    count++;
                    break;
            }
        } 
        std::cout << len << " count " << count << std::endl;
        return h1;
    }
*/
    //x1-end
/*
    template<typename TDir, typename THashValue>
    inline THashValue
    requestBucket(TDir & dir, THashValue code, unsigned len)//, Tag<TParallelTag> parallelTag)
    {
        //std::cout << code << " " << len << std::endl;
        typedef unsigned long TSize;
        // get size of the index

        // check whether bucket map is disabled and
        // where the hash should be found if no collision took place before
        TSize hlen = length(dir);
        if (hlen == 0ul) return code;
        
        TSize h1 = _hashFunction(dir, code >> _bitValue_END);
#ifdef SEQAN_OPENADDRESSING_COMPACT
        --hlen;
        h1 %= hlen;
#else
        hlen -= 2;
        h1 &= hlen;
#endif
        // was the entry empty or occupied by our code?
        
        // if not we have a collision -> probe for our code or an empty entry
        //
        // do linear probing if we need to save memory (when SEQAN_OPENADDRESSING_COMPACT is defined)
        // otherwise do quadratic probing to avoid clustering (Cormen 1998)
        TSize delta = 0;
        (void)delta;
        //std::cout << "done\n";
        unsigned count = 0;
        //return h1;
        for (auto k = 0; k < len; k++)
        {
            switch (dir[h1+k] & _bitCode){
                case 1:
                    h1 = (0 + h1 + (dir[h1 + k] >> _bitLength_END) + delta) &  hlen ;
                    //std::cout << len << " " << (dir[h1+k] >> _bitLength_END) << std::endl;
                    //std::cout << k << std::endl;
                    ++delta; k=0;
                    count++;
                    break;
                case 2:  
                    h1 = (0 + h1 + (dir[h1 + k + 1] >> _bitLength_END) + delta) &  hlen ;
                    std::cout << k << std::endl;
                    //std::cout << len << " " << (dir[h1 + k + 1] >> _bitLength_END) << std::endl;
                    ++delta; k=0;
                    count++ ;
                    break;
                case 3:
                    //std::cout << k << std::endl;
                    h1 = (h1 + 1 + delta) &  hlen ;
                    ++delta; k=0;
                    count++;
                    break;
            }
        } 
        std::cout << len << " count " << count << std::endl;
        return h1;
    }
*/

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
        //std::cout << code << " " << len << std::endl;
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

/*
    template < typename THashValue, typename THashValue2, typename TParallelTag >
    inline THashValue
    requestBucket(BucketMap<THashValue> &bucketMap, THashValue2 code, Tag<TParallelTag> parallelTag)
    {
        //std::cout << "requestBucket ";
        typedef BucketMap<THashValue> TBucketMap;
        typedef unsigned long TSize;
        // get size of the index

        // check whether bucket map is disabled and
        // where the hash should be found if no collision took place before
        TSize hlen = length(bucketMap.qgramCode);
        if (hlen == 0ul) return code;
        
        TSize h1 = _hashFunction(bucketMap, code);
        //TSize h1 = code;
#ifdef SEQAN_OPENADDRESSING_COMPACT
        --hlen;
        h1 %= hlen;
#else
        hlen -= 2;
        h1 &= hlen;
        //h1 %= hlen;
#endif
        // was the entry empty or occupied by our code?
        THashValue currentCode = atomicCas(bucketMap.qgramCode[h1], TBucketMap::EMPTY, code, parallelTag);
        if (currentCode == TBucketMap::EMPTY || currentCode == code)
            //return h1;
            return 1;

        // if not we have a collision -> probe for our code or an empty entry
        //
        // do linear probing if we need to save memory (when SEQAN_OPENADDRESSING_COMPACT is defined)
        // otherwise do quadratic probing to avoid clustering (Cormen 1998)
        TSize delta = 0;
        (void)delta;
        unsigned count = 1;
        do {
#ifdef SEQAN_OPENADDRESSING_COMPACT
            h1 = (h1 + 1) % hlen;               // linear probing guarantees that all entries are visited
#else
            h1 = (h1 + delta + 1) & hlen;       // for power2-sized tables the (i*i+i)/2 probing guarantees the same
            //h1 = (h1 + delta + 1) % hlen; 
            ++delta;
#endif
            currentCode = atomicCas(bucketMap.qgramCode[h1], TBucketMap::EMPTY, code, parallelTag);
            count++;
            //std::cout << h1 << " " << bucketMap.qgramCode[h1] << " " << code << std::endl;
        } while (currentCode != TBucketMap::EMPTY && currentCode != code);
        //std::cout << code << " " << bucketMap.qgramCode[h1] << std::endl;
        //std::cout << count << std::endl;
        //return h1;
        return count;
    }
 */
    template < typename THashValue, typename THashValue2, typename TParallelTag >
    inline THashValue
    requestBucket(BucketMap<THashValue> &bucketMap, THashValue2 code, Tag<TParallelTag> parallelTag)
    {
        typedef BucketMap<THashValue> TBucketMap;
        typedef unsigned long TSize;
        // get size of the index

        // check whether bucket map is disabled and
        // where the hash should be found if no collision took place before
        TSize hlen = length(bucketMap.qgramCode);
        if (hlen == 0ul) return code;

        TSize h1 = _hashFunction(bucketMap, code);
#ifdef SEQAN_OPENADDRESSING_COMPACT
        --hlen;
        h1 %= hlen;
#else
        hlen -= 2;
        h1 &= hlen;
#endif
        // was the entry empty or occupied by our code?
        THashValue currentCode = atomicCas(bucketMap.qgramCode[h1], TBucketMap::EMPTY, code, parallelTag);
        if (currentCode == TBucketMap::EMPTY || currentCode == code)
            return h1;

        // if not we have a collision -> probe for our code or an empty entry
        //
        // do linear probing if we need to save memory (when SEQAN_OPENADDRESSING_COMPACT is defined)
        // otherwise do quadratic probing to avoid clustering (Cormen 1998)
        TSize delta = 0;
        (void)delta;
        do {
#ifdef SEQAN_OPENADDRESSING_COMPACT
            h1 = (h1 + 1) % hlen;               // linear probing guarantees that all entries are visited
#else
            h1 = (h1 + delta + 1) & hlen;       // for power2-sized tables the (i*i+i)/2 probing guarantees the same
            ++delta;
#endif
            currentCode = atomicCas(bucketMap.qgramCode[h1], TBucketMap::EMPTY, code, parallelTag);
        } while (currentCode != TBucketMap::EMPTY && currentCode != code);
        return h1;
    }
   
    //x-begin2:
/*
    template <typename TDir, typename THashValue>
    inline THashValue
    getBucket(TDir const &dir,THashValue v1, THashValue v2, THashValue v3) 
    {
        typedef unsigned long TSize;
        // get size of the index

        // check whether bucket map is disabled and
        // where the hash should be found if no collision took place before
        TSize hlen = length(dir);
        if (hlen == 0ul) return v2;

        TSize h1 = _hashFunction(dir, _getHeadValue(v1));
#ifdef SEQAN_OPENADDRESSING_COMPACT
        --hlen;
        h1 %= hlen;
#else
        hlen -= 2;
        h1 &= hlen;
#endif
        TSize delta = 0;
        (void)delta;

        bool flag = false;
        unsigned k = 1;
        uint64_t len, start;
        while(dir[h1] | dir[h1+1])          // if dir[h1]==dir[h+1]==0 then stop;
        {
            switch(v1 ^ (dir[h1])){
                case 0: 
                    //len=_getBlockLength(dir[h1]);
                    start=_getHeadValue(dir[h1]);
                    do
                    {
                        if (v2 == _getBodyValue(dir[start + k]))
                            return start + k;
                    }while(_BlockEnd(dir[h1+k++]))
                    return _EMPTY_NODE;
                case 1:
                    v1 = v3;
                    h1 = _hashFunction(dir, _getHeadValue(v1));
                    break;
                default:
                    h1 = (h1 + delta + 1) & hlen;
                    delta++;
            }
        }
        return _EMPTY_NODE;
    }    
*/
/*
//    getBucket(TDir const &dir,THashValue v1, THashValue v2, THashValue v3) 
//    {
//        typedef unsigned long TSize;
//        // get size of the index
//
//        // check whether bucket map is disabled and
//        // where the hash should be found if no collision took place before
//        TSize hlen = length(dir);
//        if (hlen == 0ul) return v2;
//
//        TSize h1 = _hashFunction(dir, v1 >> _bitValue_END);
//#ifdef SEQAN_OPENADDRESSING_COMPACT
//        --hlen;
//        h1 %= hlen;
//#else
//        hlen -= 2;
//        h1 &= hlen;
//#endif
//        TSize delta = 0;
//        (void)delta;
//
//        bool flag = false;
//        unsigned k = 1;
//        uint64_t len = 1;
//        while(dir[h1] ^ _bitEmpty)
//        {
//            switch(v1 ^ (dir[h1]) & _bitValue_END){
//                case 0: 
//                    len=dir[h1+1] >> _bitLength;
//                    while(!flag)
//                    {
//                        if(v2 == (dir[h1 + k] & _bitValue2)||(k > len))
//                            return h1 + k;
//                        k++;
//                    }
//                    return h1 + 1;
//                case 1:
//                    v1 = v3;
//                    h1 = _hashFunction(dir, v1 >> _bitValue_END);
//                    flag = true;
//                    break;
//                default:
//                    h1 = (h1 + (dir[h1] >> _bitLength_END) +delta) & hlen;
//                    delta++;
//            }
//        }
//        return h1;
//    }
*/
    //x-end2:
    
    template <typename TObject, unsigned TSPAN, unsigned TWEIGHT, typename TValue, typename TSpec>
    inline typename Value< Shape<TValue, MinimizerShape<TSPAN, TWEIGHT, TSpec> > >::Type 
    getDir(Index<TObject, IndexQGram<MinimizerShape<TSPAN, TWEIGHT>, OpenAddressing> > const & index, Shape<TValue, MinimizerShape<TSPAN, TWEIGHT, TSpec> > const & shape)
    {
        typedef unsigned long TSize;
        typedef typename Value< Shape<TValue, MinimizerShape<TSPAN, TWEIGHT, TSpec> > >::Type THashValue;
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

    template < typename THashValue, typename THashValue2 >
    inline THashValue
    getBucket(BucketMap<THashValue> const &bucketMap, THashValue2 code)
    {
        typedef BucketMap<THashValue> TBucketMap;
        typedef unsigned long TSize;
        // get size of the index

        // check whether bucket map is disabled and
        // where the hash should be found if no collision took place before
        TSize hlen = length(bucketMap.qgramCode);
        if (hlen == 0ul) return code;

        TSize h1 = _hashFunction(bucketMap, code);
#ifdef SEQAN_OPENADDRESSING_COMPACT
        --hlen;
        h1 %= hlen;
#else
        hlen -= 2;
        h1 &= hlen;
#endif

        // probe for our code or an empty entry
        //
        // do linear probing if we need to save memory (when SEQAN_OPENADDRESSING_COMPACT is defined)
        // otherwise do quadratic probing to avoid clustering (Cormen 1998)
        TSize delta = 0;
        (void)delta;
        while (bucketMap.qgramCode[h1] != code && bucketMap.qgramCode[h1] != TBucketMap::EMPTY)
        {
#ifdef SEQAN_OPENADDRESSING_COMPACT
            h1 = (h1 + 1) % hlen;               // linear probing guarantees that all entries are visited
#else
            h1 = (h1 + delta + 1) & hlen;       // for power2-sized tables the (i*i+i)/2 probing guarantees the same
            ++delta;
#endif
        }
        return h1;
    }
 

    template <typename TBucketMap>
    inline bool _emptyBucketMap(TBucketMap const &bucketMap)
    {
        return empty(bucketMap);
    }

    inline bool _emptyBucketMap(Nothing const &)
    {
        return false;
    }

    template <typename TObject, typename TShapeSpec>
    inline __int64 _fullDirLength(Index<TObject, IndexQGram<TShapeSpec, OpenAddressing> > const &index)
    {
        typedef Index<TObject, IndexQGram<TShapeSpec, OpenAddressing> >    TIndex;
        typedef typename Fibre<TIndex, QGramDir>::Type                        TDir;
        typedef typename Fibre<TIndex, FibreShape>::Type                    TShape;
        typedef typename Host<TShape>::Type                                    TTextValue;
        typedef typename Value<TDir>::Type                                    TDirValue;
        typedef typename Value<TShape>::Type                                THashValue;

        double num_qgrams = _qgramQGramCount(index) * index.alpha;
        double max_qgrams = pow((double)ValueSize<TTextValue>::VALUE, (double)weight(indexShape(index)));
        __int64 qgrams;

        // compare size of open adressing with 1-1 mapping and use the smaller one
        if (num_qgrams * (sizeof(TDirValue) + sizeof(THashValue)) < max_qgrams * sizeof(TDirValue))
        {
            qgrams = (__int64)ceil(num_qgrams);
#ifndef SEQAN_OPENADDRESSING_COMPACT
            __int64 power2 = 1;
            while (power2 < qgrams)
                power2 <<= 1;
            qgrams = power2;
#endif
            resize(const_cast<TIndex &>(index).bucketMap.qgramCode, qgrams + 1, Exact());
        } else
        {
            qgrams = (__int64)ceil(max_qgrams);
            clear(const_cast<TIndex &>(index).bucketMap.qgramCode);    // 1-1 mapping, no bucket map needed
        }
        return qgrams + 1;
    }

    template <typename TObject, unsigned TSPAN, unsigned TWEIGHT>
    inline __int64 _fullDirLength(Index<TObject, IndexQGram<MinimizerShape<TSPAN, TWEIGHT>, OpenAddressing> > const &index)
    {
        typedef Index<TObject, IndexQGram<MinimizerShape<TSPAN, TWEIGHT>, OpenAddressing> >    TIndex;
        typedef typename Fibre<TIndex, FibreShape>::Type                    TShape;
        typedef typename Host<TShape>::Type                                    TTextValue;

        double num_qgrams = _qgramQGramCount(index) * index.alpha;
        double max_qgrams = 2*pow((double)ValueSize<TTextValue>::VALUE, (double)length(indexShape(index)));
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

template <unsigned TSpan, unsigned TWeight>
void _qgramClearDir(Index<StringSet<DnaString>, IndexQGram<MinimizerShape<TSpan, TWeight>, OpenAddressing> > & index)
{
    typedef Shape<Dna, MinimizerShape<TSpan, TWeight> > TM_Shape;
    typedef typename Value<TM_Shape>::Type HValue;
    resize (indexDir(index), _fullDirLength(index) + lengthSum(indexText(index)) + 2);
    index.start = _fullDirLength(index);
    index._Empty_Dir_ = 0;
    for (HValue k = 0; k < length(index.dir); k++) 
    {
        index.dir[k] = _bitEmpty;
    }
    std::cout << "        _qgramClearDir():" << std::endl;
    std::cout << "            _fullDirLength(index) = " << _fullDirLength(index) << std::endl;
    std::cout << "            lengh(index.dir) = " << length(index.dir) << std::endl;
    std::cout << "            End _qgramClearDir()" << std::endl;
}

template <unsigned TSpan, unsigned TWeight>
void _qgramCountQGrams2(Index<StringSet<DnaString>, IndexQGram<MinimizerShape<TSpan, TWeight>, OpenAddressing > > & index)
{
    typedef Shape<Dna, MinimizerShape<TSpan, TWeight> > TM_Shape;
    typedef Iterator<String<Dna> >::Type TIter;
    typedef typename Value<TM_Shape>::Type HValue;
    typedef std::tuple<HValue, HValue, HValue, HValue> HTuple;
    typedef String<HTuple> Stringtuple;

    TM_Shape shape;
    Stringtuple hs, hs1;
    HValue  m = 0, sum = 0;

    double time = sysTime();
    resize(hs, lengthSum(indexText(index)) - length(indexText(index)) * (shape.span - 1) + 1);

    std::cout << "        _qgramCountQGrams() sysTime(): " << sysTime() - time << std::endl;
    std::cout << "            lengthSum(StringSet) = " << lengthSum(indexText(index)) << std::endl;
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
    std::cout << "            make_tuple sysTime(): " << sysTime() - time << std::endl;
    //hs[length(hs) - 1] = std::make_tuple((HValue)0, (HValue)0, (HValue)0, (HValue)1);
    std::sort(begin(hs), end(hs) - 1,
        [](const HTuple &a, const HTuple &b)
        {return (std::get<0>(a) > std::get<0>(b)||(std::get<0>(a) == std::get<0>(b) && std::get<1>(a) > std::get<1>(b)));});
    hs[length(hs) - 1] = std::make_tuple((HValue)1, (HValue)1, (HValue)0, (HValue)1);
    std::cout << "            sort sysTime(): " << sysTime() - time << " " << std::get<0>(hs[length(hs) - 2]) << " " << std::get<1>(hs[length(hs)- 2]) << std::endl;
    HValue countx = 1, counth = 1, tmp = 0, countdh = 0, countb = 0, hk = 0;
    resize(index.sa, length(hs) - 1);
    for (HValue k = 1; k < length(hs); k++)
    {
        index.sa[k - 1] = std::get<3>(hs[k-1]);
        if (std::get<1>(hs[k]) != std::get<1>(hs[k - 1]))
        { _setBodyNode(index.dir[index.start + hk], std::get<2>(hs[k-1]), _BodyType_code, tmp);
        if (std::get<0>(hs[k - 1]) == 0 && std::get<2>(hs[k - 1]) == 441204)
            std::cout <<"_getBodyValue = " << _getBodyValue(index.dir[index.start + hk]) << std::endl;
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
    std::cout << std::endl;
    std::cout << counth << std::endl;
    resize(index.dir, index.start + countdh + 10);
    _setBodyNode(index.dir[index.start + countdh], _bitEmpty, _BodyType_code, counth - 1); 
    _setBodyType_Begin(index.dir[index.start + countdh]);
    index._Empty_Dir_ = index.start + countdh + 1;
    //_setBodyNode(index.dir[index.start + countdh], _bitEmpty, _BodyTypeEnd_code, counth - 1); 
    //_setBodyType_Begin(index.dir[index.start + countdh]);
    //index._Empty_Dir_ = index.start + countdh;
    //_setBodyNode(index.dir[index.start + countdh + 1], _bitEmpty, _BodyTypeEnd_code, counth - 1); 
    //_setBodyType_Begin(index.dir[index.start + countdh + 1]);
    std::cout << "            End _qgramCountQGrams() sysTime(): " << sysTime() - time << std::endl;
}



template <unsigned TSpan, unsigned TWeight>
void _qgramCountQGrams3(Index<StringSet<DnaString>, IndexQGram<MinimizerShape<TSpan, TWeight>, OpenAddressing > > & index)
{
    typedef Shape<Dna, MinimizerShape<TSpan, TWeight> > TM_Shape;
    typedef Iterator<String<Dna> >::Type TIter;
    typedef typename Value<TM_Shape>::Type HValue;
    typedef std::tuple<HValue, HValue, HValue> HTuple;
    typedef String<HTuple> Stringtuple;

    TM_Shape shape;
    Stringtuple hs, hs1;
    HValue  m = 0, sum = 0;

    double time = sysTime();
    resize(hs, lengthSum(indexText(index)) - length(indexText(index)) * (shape.span - 1) + 1);

    std::cout << "        _qgramCountQGrams() sysTime(): " << sysTime() - time << std::endl;
    std::cout << "            lengthSum(StringSet) = " << lengthSum(indexText(index)) << std::endl;
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
    std::cout << "            make_tuple sysTime(): " << sysTime() - time << std::endl;
    //hs[length(hs) - 1] = std::make_tuple((HValue)0, (HValue)0, (HValue)0, (HValue)1);
    std::sort(begin(hs), end(hs) - 1,
        [](const HTuple &a, const HTuple &b)
        {return (std::get<0>(a) < std::get<0>(b)||(std::get<0>(a) == std::get<0>(b) && std::get<1>(a) > std::get<1>(b)));});
    hs[length(hs) - 1] = std::make_tuple((HValue)1, (HValue)0, (HValue)1);
    std::cout << "            sort sysTime(): " << sysTime() - time << std::endl;
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
    std::cout << std::endl;
    std::cout << counth << std::endl;
    resize(index.dir, index.start + countdh + 10);
    _setBodyNode(index.dir[index.start + countdh], _bitEmpty, _BodyType_code, counth - 1); 
    _setBodyType_Begin(index.dir[index.start + countdh]);
    index._Empty_Dir_ = index.start + countdh + 1;
    //_setBodyNode(index.dir[index.start + countdh], _bitEmpty, _BodyTypeEnd_code, counth - 1); 
    //_setBodyType_Begin(index.dir[index.start + countdh]);
    //index._Empty_Dir_ = index.start + countdh;
    //_setBodyNode(index.dir[index.start + countdh + 1], _bitEmpty, _BodyTypeEnd_code, counth - 1); 
    //_setBodyType_Begin(index.dir[index.start + countdh + 1]);
    std::cout << "            End _qgramCountQGrams() sysTime(): " << sysTime() - time << std::endl;
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

static String<Pair<unsigned, unsigned> > _SORT_PARA = {std::make_pair(8,4)};//{(8,4),(9,4),(9,4),(8,5),(8,5),(7,6),(8,6),(8,6),(8,6),(9,6)}; //weight:[16,26)

inline Pair<unsigned, unsigned> _getSortPara(uint64_t shapeWeight)
{
    return _SORT_PARA[shapeWeight - 16];
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

template <unsigned TSPAN, unsigned TWEIGHT>
void _createValueArray2(StringSet<DnaString> & reads, String<Pair<uint64_t, uint64_t> > & hs, Shape<Dna, MinimizerShape<TSPAN, TWEIGHT> > & shape, int step, int l)
{
    typedef String<Pair<uint64_t, uint64_t> > StringPair;

    //TShape shape;
    StringPair tmp;
    String<uint64_t> tmp3;
    uint64_t p = 0, q=0, c = 0, n = -1, pre = ~0, count = 0, mask = ((uint64_t)1 << 63);
    //std::cout << "    _createValueArray() " << std::endl;
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

    std::cout << "        loading time " << sysTime() - time << std::endl;
    c = tmp[0].i2;
 p = q = count = 0;
    n = -1;
    //_sort3(begin(tmp), end(tmp), step, l);         // sort parameters 1
    _sort3(begin(tmp), end(tmp), _getSortPara(shape.weight).i1, 
            _getSortPara(shape.weight).i2);         // sort parameters 1
    std::cout << "        sort xvalue time " << sysTime() - time << std::endl;
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
    }
    std::cout << "        xvalue expand " << sysTime() - time << std::endl;
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
                std::sort(begin(hs) + k -count, begin(hs) + k, [](Pair<uint64_t, uint64_t> & a,
                Pair<uint64_t, uint64_t> & b){return a.i2 > b.i2;});
            count = 0;
        }
        count++;
   }


   std::cout << "        End sort sysTime(): " <<  sysTime() - time << std::endl;
}




template <unsigned TSpan, unsigned TWeight>
void _qgramCountQGrams(Index<StringSet<DnaString>, IndexQGram<MinimizerShape<TSpan, TWeight>, OpenAddressing > > & index)
{
    typedef Shape<Dna, MinimizerShape<TSpan, TWeight> > TM_Shape;
    typedef Iterator<String<Dna> >::Type TIter;
    typedef typename Value<TM_Shape>::Type HValue;
    typedef std::tuple<HValue, HValue, HValue> HTuple;
    typedef Pair<uint64_t, uint64_t> PairH;
    typedef String<PairH> StringPairH;
    typedef String<HTuple> StringTuple;
    StringSet<DnaString> reads;

    TM_Shape shape;
    StringPairH hs, hs1;
    HValue  m = 0, sum = 0;

    double time = sysTime();

    resize(hs, lengthSum(indexText(index)) - length(indexText(index)) * (shape.span - 1) + 1);

    //std::cout << "        _qgramCountQGrams() sysTime(): " << sysTime() - time << std::endl;
    //std::cout << "            lengthSum(StringSet) = " << lengthSum(indexText(index)) << std::endl;

     _createValueArray2(indexText(index), hs, shape, 9, 5);
    hs[length(hs) - 1].i1 = 1;
    _setBodyNode(hs[length(hs) - 1].i2,  (HValue)0, (HValue)0, (HValue)1);//std::make_tuple((HValue)1, (HValue)0, (HValue)1);
    std::cout << "            sort sysTime(): " << sysTime() - time << std::endl;
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
    //std::cout << std::endl;
    //std::cout << counth << std::endl;
    resize(index.dir, index.start + countdh + 10);
    _setBodyNode(index.dir[index.start + countdh], _bitEmpty, _BodyType_code, counth - 1); 
    _setBodyType_Begin(index.dir[index.start + countdh]);
    index._Empty_Dir_ = index.start + countdh + 1;
    //_setBodyNode(index.dir[index.start + countdh], _bitEmpty, _BodyTypeEnd_code, counth - 1); 
    //_setBodyType_Begin(index.dir[index.start + countdh]);
    //index._Empty_Dir_ = index.start + countdh;
    //_setBodyNode(index.dir[index.start + countdh + 1], _bitEmpty, _BodyTypeEnd_code, counth - 1); 
    //_setBodyType_Begin(index.dir[index.start + countdh + 1]);
    std::cout << "            End _qgramCountQGrams() sysTime(): " << sysTime() - time << std::endl;
}

template <unsigned TSpan, unsigned TWeight>
void createQGramIndexDirOnly(Index<StringSet<DnaString>, IndexQGram<MinimizerShape<TSpan, TWeight>, OpenAddressing > >& index)
{
    double time = sysTime(); 
    //std::cout << "    createQGramIndexDirOnly() sysTime(): " << std::endl;
    _qgramClearDir(index);
    _qgramCountQGrams(index);
    //std::cout << "        End createQGramIndexDirOnly() sysTime(): " << sysTime() - time << std::endl;
    std::cout << sysTime() - time << std::endl;
}
/*
template <unsigned TSpan, unsigned TWeight>
void createQGramIndexDirOnly2(Index<StringSet<DnaString>, IndexQGram<MinimizerShape<TSpan, TWeight>, OpenAddressing > >& index)
{
    double time = sysTime(); 
    std::cout << "    createQGramIndexDirOnly() sysTime(): " << std::endl;
    _qgramClearDir(index);
    _qgramCountQGrams2(index);
    std::cout << "        index.dir[index._Empty_Dir_] = " << index.dir[index._Empty_Dir_] << std::endl << "        length(index.dir) = " << length(index.dir) << std::endl;
    std::cout << "        End createQGramIndexDirOnly() sysTime(): " << sysTime() - time << std::endl;
}
*/
// ----------------------------------------------------------------------------
// Function open()
// ----------------------------------------------------------------------------

template < typename TText, typename TShapeSpec>
inline bool open(Index<TText, IndexQGram<TShapeSpec, OpenAddressing> > &index, const char *fileName, int openMode)
{
    String<char> name;
    
    name = fileName;    append(name, ".txt");
    if (!open(getFibre(index, QGramText()), toCString(name), openMode)) return false;
    
    name = fileName;    append(name, ".sa");
    if (!open(getFibre(index, QGramSA()), toCString(name), openMode)) return false;
    
    name = fileName;    append(name, ".dir");
    if (!open(getFibre(index, QGramDir()), toCString(name), openMode)) return false;

    name = fileName;    append(name, ".bkt");
    if (!open(getFibre(index, QGramBucketMap()).qgramCode, toCString(name), openMode)) return false;

    return true;
}

template <typename TText, typename TShapeSpec>
inline bool open(Index<TText, IndexQGram<TShapeSpec, OpenAddressing> > &index, const char *fileName)
{
    return open(index, fileName, OPEN_RDONLY);
}

// ----------------------------------------------------------------------------
// Function save()
// ----------------------------------------------------------------------------

template <typename TText, typename TShapeSpec>
inline bool save(Index<TText, IndexQGram<TShapeSpec, OpenAddressing> > &index, const char *fileName, int openMode)
{
    String<char> name;
    
    name = fileName;    append(name, ".txt");
    if (!save(getFibre(index, QGramText()), toCString(name), openMode)) return false;
    
    name = fileName;    append(name, ".sa");
    if (!save(getFibre(index, QGramSA()), toCString(name), openMode)) return false;
    
    name = fileName;    append(name, ".dir");
    if (!save(getFibre(index, QGramDir()), toCString(name), openMode)) return false;

    name = fileName;    append(name, ".bkt");
    if (!save(getFibre(index, QGramBucketMap()).qgramCode, toCString(name), openMode)) return false;

    return true;
}

template <typename TText, typename TShapeSpec>
inline bool save(Index<TText, IndexQGram<TShapeSpec, OpenAddressing> > &index, const char *fileName)
{
    return save(index, fileName, OPEN_WRONLY | OPEN_CREATE);
}




}

#endif //#ifndef SEQAN_HEADER_...
