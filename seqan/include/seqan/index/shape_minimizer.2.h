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

#ifndef SEQAN_HEADER_SHAPE_MINIMIZER_H
#define SEQAN_HEADER_SHAPE_MINIMIZER_H

namespace seqan
{

const float _boundAlpha = 0.8;
struct ReverseComplement_;
typedef Tag<ReverseComplement_> const   ReverseComplement;

// ----------------------------------------------------------------------------
// Struct MinimizerShape
// ----------------------------------------------------------------------------

template <unsigned TSPAN, unsigned TWEIGHT, typename TSpec = void>
struct MinimizerShape;

// ----------------------------------------------------------------------------
// Class Shape<MinimizerShape>
// ----------------------------------------------------------------------------

template <typename TValue, unsigned TSPAN, unsigned TWEIGHT, typename TSpec>
class Shape<TValue, MinimizerShape<TSPAN, TWEIGHT, TSpec> >
{
public:
    typedef typename Value<Shape>::Type THashValue;

    unsigned span;
    unsigned weight;
    int first;
    int bound;
    THashValue hValue;      //minimizer hash
    THashValue XValue;
    THashValue YValue;
    static const THashValue leftFactor = Power<ValueSize<TValue>::VALUE, TSPAN - 1>::VALUE;
    static const THashValue m_leftFactor = Power<ValueSize<TValue>::VALUE, TWEIGHT - 1>::VALUE;
    TValue  leftChar;

    Shape():
        span(TSPAN),
        weight(TWEIGHT),
        first(-1),
        bound((unsigned)(TWEIGHT * _boundAlpha))
    {}
};

// ----------------------------------------------------------------------------
// Metafunction LENGTH
// ----------------------------------------------------------------------------

template <typename TValue, unsigned TSPAN, unsigned TWEIGHT, typename TSpec>
struct LENGTH<Shape<TValue, MinimizerShape<TSPAN, TWEIGHT, TSpec> > >
{
    enum { VALUE = TSPAN };
};

// ----------------------------------------------------------------------------
// Metafunction WEIGHT
// ----------------------------------------------------------------------------

template <typename TValue, unsigned TSPAN, unsigned TWEIGHT, typename TSpec>
struct WEIGHT<Shape<TValue, MinimizerShape<TSPAN, TWEIGHT, TSpec> > >
{
    enum { VALUE = TWEIGHT - TSPAN};
};

// ----------------------------------------------------------------------------
// Metafunction DELTA=LENGTH-WEIGHT
// ----------------------------------------------------------------------------

template <typename T>
struct WINDOW_BIT_SIZE;

template <typename TValue, unsigned TSPAN, unsigned TWEIGHT, typename TSpec>
struct WINDOW_BIT_SIZE<Shape<TValue, MinimizerShape<TSPAN, TWEIGHT, TSpec> > >
{
    enum { VALUE = 2 * (TSPAN - TWEIGHT)};
};

// ----------------------------------------------------------------------------
// Function weight()
// ----------------------------------------------------------------------------

template <typename TValue, unsigned TSPAN, unsigned TWEIGHT, typename TSpec>
inline SEQAN_HOST_DEVICE
typename Size< Shape<TValue, MinimizerShape<TSPAN, TWEIGHT, TSpec> > >::Type
weight(Shape<TValue, MinimizerShape<TSPAN, TWEIGHT, TSpec> > const &me)
{
    return me.weight;
}

// ----------------------------------------------------------------------------
// Function _minHash()
// ----------------------------------------------------------------------------
// return lexicographically smaller hash as the minimizer

template <typename TValue, unsigned TSPAN, unsigned TWEIGHT, typename TSpec, typename TIter>
inline typename Value<Shape<TValue, MinimizerShape<TSPAN, TWEIGHT, TSpec> > >::Type
hash(Shape<TValue, MinimizerShape<TSPAN, TWEIGHT, TSpec> > &me, TIter const &it)
{
    Range<TIter> range(it, it + length(me));
    Shape<TValue, UngappedShape<TSPAN> > u_tmpShape;

    return me.hValue;
}

// ----------------------------------------------------------------------------
// Function hashNext()
// ----------------------------------------------------------------------------

uint64_t hash_key=1;
template <typename TValue, unsigned TSPAN, unsigned TWEIGHT, typename TSpec, typename TIter>
inline void hashInit(Shape<TValue, MinimizerShape<TSPAN, TWEIGHT, TSpec> > &me, TIter const &it)
{
    typedef typename Size< Shape<TValue, SimpleShape> >::Type    TSize;

        SEQAN_ASSERT_GT((unsigned)me.span, 0u);

        me.leftChar = 0;
        me.hValue = ordValue(*it);
        hash_key=1;
        hash_key = (hash_key << (TSPAN*2 -2 )) - 1;
        for(TSize i = 2; i < me.span; ++i) {
            //me.hValue = me.hValue * ValueSize<TValue>::VALUE + ordValue((TValue)*(it +i-2));
            me.hValue = (me.hValue << 2) + ordValue((TValue)*(it +i-2));

        }
}

static const int MultiplyDeBruijnBitPosition[64] =
{
    0,  1,  2, 53,  3,  7, 54, 27,
    4, 38, 41,  8, 34, 55, 48, 28,
   62,  5, 39, 46, 44, 42, 22,  9,
   24, 35, 59, 56, 49, 18, 29, 11,
   63, 52,  6, 26, 37, 40, 33, 47,
   61, 45, 43, 21, 23, 58, 17, 10,
   51, 25, 36, 32, 60, 20, 57, 16,
   50, 31, 19, 15, 30, 14, 13, 12,
};
    template <typename TValue, unsigned TSPAN, unsigned TWEIGHT, typename TSpec, typename TIter>
    inline typename Value< Shape<TValue, MinimizerShape<TSPAN, TWEIGHT, TSpec> > >::Type
    hashNext(Shape<TValue, MinimizerShape<TSPAN, TWEIGHT, TSpec> > &me, TIter const &it)
    {
    SEQAN_CHECKPOINT
        // remove first, shift left, and add next character
        //typedef typename Value< Shape<TValue, TSpec> >::Type    THValue;
        typedef typename Size< Shape<TValue, TSpec> >::Type        TSize;
        unsigned t;
        SEQAN_ASSERT_GT((unsigned)me.span, 0u);

        //me.hValue =
        //    (me.hValue - ordValue(me.leftChar) * (THValue)me.leftFactor) * ValueSize<TValue>::VALUE
        //    + ordValue((TValue)*(it + ((TSize)me.span - 1)));
        //me.leftChar = *it;
        me.hValue=((me.hValue & hash_key)<<2)+ordValue((TValue)*(it + ((TSize)me.span - 1)));
        //t = MultiplyDeBruijnBitPosition[(((me.hValue & -me.hValue) * 0x022fdd63cc95386d)) >> 58];
        t = __builtin_clzll(me.hValue);
        //t = (__builtin_clz(me.hValue) < WINDOW_BIT_SIZE<Shape<TValue, MinimizerShape<TSPAN, TWEIGHT, TSpec> > >::VALUE)?__builtin_clz(me.hValue):WINDOW_BIT_SIZE<Shape<TValue, MinimizerShape<TSPAN, TWEIGHT, TSpec> > >::VALUE;
        t += (WINDOW_BIT_SIZE<Shape<TValue, MinimizerShape<TSPAN, TWEIGHT, TSpec> > >::VALUE - t) & ((WINDOW_BIT_SIZE<Shape<TValue, MinimizerShape<TSPAN, TWEIGHT, TSpec> > >::VALUE  - t) >> (sizeof(unsigned) * CHAR_BIT - 1));
        me.XValue=me.hValue << t >> WINDOW_BIT_SIZE<Shape<TValue, MinimizerShape<TSPAN, TWEIGHT, TSpec> > >::VALUE;
        me.YValue=(me.hValue >> (t-1) << (64-t))+(me.hValue <<t>>t)+t;
        return me.XValue;
    }

}


// namespace seqan



#endif
