// ==========================================================================
//                           Mapping SMRT reads 
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
// Autor: cxpan <chenxu.pan@fu-berlin.de>
// ==========================================================================

#ifndef SEQAN_HEADER_SHAPE_PM_H
#define SEQAN_HEADER_SHAPE_PM_H
//#include <seqan/basic.h>
//#include <seqan/index.h>

namespace seqan{

const float _boundAlpha = 0.8;
//struct ReverseComplement_;
//typedef Tag<ReverseComplement_> const   ReverseComplement;

// ----------------------------------------------------------------------------
// Struct Minimizer
// ----------------------------------------------------------------------------


struct MiniValueBit
{
    enum{VALUEBIT = 64};
};

struct MiniHEX{               
    enum {HEX = 8 };
};

template <unsigned shapeLength>
struct MiniWeight{
    enum{ WEIGHT = shapeLength - MiniHEX::HEX};
};


struct MiniYBit{
    enum{ YBIT = (MiniHEX::HEX << 1) + Log2<1 + MiniHEX::HEX>::VALUE};
};

template <unsigned TSPAN, unsigned TWEIGHT = MiniWeight<TSPAN>::WEIGHT, typename TSpec = void>
struct Minimizer;
typedef Minimizer<0, 0> SimpleMShape;

// ----------------------------------------------------------------------------
// Class Shape<Minimizer>
// ----------------------------------------------------------------------------

//template <typename TValue, typename TSpec>>
//typedef seqan::Shape<TValue, TSpec> Shape;

template <typename TValue, unsigned TSPAN, unsigned TWEIGHT, typename TSpec>
class Shape<TValue, Minimizer<TSPAN, TWEIGHT, TSpec> >
{
public:
    typedef typename seqan::Value<Shape>::Type THashValue;

    unsigned span;
    unsigned weight;
    THashValue hValue;     //hash value 
    THashValue crhValue;    //ReverseComplement
    THashValue XValue;     //minimizer 
    THashValue YValue;     //Y(h,x)
    THashValue strand;
    int  leftChar;
    int x;

    Shape():
        span(TSPAN),
        weight(TWEIGHT),
        hValue(0),
        crhValue(0),
        XValue(0),
        YValue(0),
        strand(0),
        leftChar(0),
        x(0)
    {}
};

// ----------------------------------------------------------------------------
// Metafunction LENGTH
// ----------------------------------------------------------------------------

template <typename TValue, unsigned TSPAN, unsigned TWEIGHT, typename TSpec>
struct LENGTH<Shape<TValue, Minimizer<TSPAN, TWEIGHT, TSpec> > >
{
    enum { VALUE = TSPAN };
};

// ----------------------------------------------------------------------------
// Metafunction WEIGHT
// ----------------------------------------------------------------------------

template <typename TValue, unsigned TSPAN, unsigned TWEIGHT, typename TSpec>
struct WEIGHT<Shape<TValue, Minimizer<TSPAN, TWEIGHT, TSpec> > >
{
    enum { VALUE = TWEIGHT - TSPAN};
};

template <typename T>
struct WGHT;

template <typename TValue, unsigned TSPAN, unsigned TWEIGHT, typename TSpec>
struct WGHT<Shape<TValue, Minimizer<TSPAN, TWEIGHT, TSpec> > >
{
    enum { VALUE = TWEIGHT - TSPAN};
};

// ----------------------------------------------------------------------------
// Metafunction DELTA=LENGTH-WEIGHT
// ----------------------------------------------------------------------------

template <typename T>
struct WINDOW_BIT_SIZE;

template <typename TValue, unsigned TSPAN, unsigned TWEIGHT, typename TSpec>
struct WINDOW_BIT_SIZE<Shape<TValue, Minimizer<TSPAN, TWEIGHT, TSpec> > >
{
    enum { VALUE = 2 * (TSPAN - TWEIGHT)};
};

template <typename TValue, unsigned TSPAN, unsigned TWEIGHT, typename TSpec>
 SEQAN_HOST_DEVICE
typename Size< Shape<TValue, Minimizer<TSPAN, TWEIGHT, TSpec> > >::Type
length(Shape<TValue, Minimizer<TSPAN, TWEIGHT, TSpec> > const &me)
{
    return me.span;
}

// ----------------------------------------------------------------------------
// Function weight()
// ----------------------------------------------------------------------------

template <typename TValue, unsigned TSPAN, unsigned TWEIGHT, typename TSpec>
 SEQAN_HOST_DEVICE
typename Size< Shape<TValue, Minimizer<TSPAN, TWEIGHT, TSpec> > >::Type
weight(Shape<TValue, Minimizer<TSPAN, TWEIGHT, TSpec> > const &me)
{
    return me.weight;
}

//-----------------------------------------------------------------------------
// Function resize()
//-----------------------------------------------------------------------------
template <typename TValue>
 void resize(Shape<TValue, SimpleMShape> & me, unsigned new_span, unsigned new_weight)
{   
    //typedef typename Value< Shape<TValue, SimpleShape> >::Type    THValue;
    //me.leftFactor = _intPow((THValue)ValueSize<TValue>::VALUE, new_span - 1); 
    //me.leftFactor2 = (_intPow((THValue)ValueSize<TValue>::VALUE, new_span) - 1) / (ValueSize<TValue>::VALUE - 1); 
    me.span = new_span;
    me.weight = new_weight;
}   

template <typename TValue, typename TSize, unsigned TSPAN, unsigned TWEIGHT>
 void resize(Shape<TValue, Minimizer<TSPAN, TWEIGHT> > & me, TSize new_span, TSize new_weight)
{   
    //typedef typename Value< Shape<TValue, SimpleShape> >::Type    THValue;
    //me.leftFactor = _intPow((THValue)ValueSize<TValue>::VALUE, new_span - 1); 
    //me.leftFactor2 = (_intPow((THValue)ValueSize<TValue>::VALUE, new_span) - 1) / (ValueSize<TValue>::VALUE - 1); 
    me.span = new_span;
    me.weight = new_weight;
}

// ----------------------------------------------------------------------------
// Function _minHash()
// ----------------------------------------------------------------------------
// return lexicographically smaller hash as the minimizer

template <typename TValue, unsigned TSPAN, unsigned TWEIGHT, typename TSpec, typename TIter>
 typename Value<Shape<TValue, Minimizer<TSPAN, TWEIGHT, TSpec> > >::Type
hash(Shape<TValue, Minimizer<TSPAN, TWEIGHT, TSpec> > &me, TIter const &it)
{
      typedef typename Size< Shape<TValue, SimpleShape> >::Type    TSize;
      me.leftChar = 0;
        me.hValue = 0;
        for (TSize i = 0; i < me.span; ++i)
        {
            me.hValue = (me.hValue << 2) + ordValue((TValue)*(it + i));
        }
    return me.hValue;
}

// ----------------------------------------------------------------------------
// Function hashNext()
// ----------------------------------------------------------------------------


template <unsigned span> 
struct MASK
{
    static const uint64_t VALUE = (1ULL << span) - 1;
};

static const uint64_t COMP4 = 3;
static const int  ordC = 3;

template <typename TValue>
 void phi(TValue & h)
{
    h+=h<<31;
}

template <typename TValue, unsigned TSPAN, unsigned TWEIGHT, typename TSpec, typename TIter>
 void hashInit(Shape<TValue, Minimizer<TSPAN, TWEIGHT, TSpec> > &me, TIter const &it)
{
    typedef typename Size< Shape<TValue, SimpleShape> >::Type    TSize;

        SEQAN_ASSERT_GT((unsigned)me.span, 0u);

        me.leftChar = 0;
        //me.hValue = ordValue(*it);
        me.hValue = 0;
        //for(TSize i = 2; i < me.span; ++i) {
        //    //me.hValue = me.hValue * ValueSize<TValue>::VALUE + ordValue((TValue)*(it +i-2));
        //    me.hValue = (me.hValue << 2) + ordValue((TValue)*(it +i-2));
        for (TSize i = 0; i < me.span - 1; ++i)
        {
            me.hValue = (me.hValue << 2) + ordValue((TValue)*(it + i));
        }
        me.x = 0;
        //}
}

template <unsigned TSPAN, unsigned TWEIGHT, typename TSpec, typename TIter>
 uint64_t hashInit(Shape<Dna5, Minimizer<TSPAN, TWEIGHT, TSpec> > &me, TIter const &it)
{
    std::cout << "hashInitxxx\n";
    SEQAN_ASSERT_GT((unsigned)me.span, 0u);

    me.leftChar = 0;
    //me.hValue = ordValue(*it);
    me.hValue = 0;
    me.crhValue = 0;
    //hash_key = ((uint64_t)1 << (me.span*2 -2 )) - 1;
    //hash_key3 = ((uint64_t)1 << (me.span*2 )) - 1;
    me.leftChar = 0;
    me.x = me.leftChar- ordC;
    uint64_t k =0, count = 0; //COMP for complemnet value;

    while (count < me.span)
    {
        if (ordValue(*(it + k + count)) == 4)
        {
            k += count + 1;
            count = 0;
        }
        else
            count++;
    }
    unsigned bit = 2;
    for (unsigned i = 0; i < me.span - 1; ++i)
    {
        uint64_t val = ordValue (*(it + k + i));
        me.x += (val << 1) - ordC;
        me.hValue = (me.hValue << 2) + val;
        me.crhValue += ((COMP4 - val) << bit);
        bit += 2;
        
    }
    //}
    return k;
}
//
/**
 *init for hashNexthS 
 */
template <unsigned TSPAN, unsigned TWEIGHT, typename TSpec, typename TIter>
 uint64_t hashInit_hs(Shape<Dna5, Minimizer<TSPAN, TWEIGHT, TSpec> > &me, TIter const &it, int d = 0)
{
        me.hValue = 0;
        for (unsigned i = d; i < me.span - 1 + d; ++i)
        {
            me.hValue = (me.hValue << 2) + ordValue (*(it + i));
        }
        me.hValue <<= (d << 1);
        return 0;
}

/*
 * this hashNext function is for index only collect mini hash value [minindex]
 */ 
template <typename TValue, unsigned TSPAN, unsigned TWEIGHT, typename TSpec, typename TIter>
 typename Value< Shape<TValue, Minimizer<TSPAN, TWEIGHT, TSpec> > >::Type
hashNext(Shape<TValue, Minimizer<TSPAN, TWEIGHT, TSpec> > &me, TIter const &it)
{
    std::cout << "hashNextxxx\n";
    //typedef typename Size< Shape<TValue, TSpec> >::Type  TSize;
    SEQAN_ASSERT_GT((unsigned)me.span, 0u);
    uint64_t v1;
    unsigned t, span = TSPAN << 1, weight = TWEIGHT << 1;
    uint64_t v2 = ordValue((TValue)*(it + me.span - 1));
    me.hValue=((me.hValue & MASK<TSPAN * 2 - 2>::VALUE)<<2)+ v2;
    me.crhValue=((me.crhValue >> 2) & MASK<TSPAN * 2 - 2>::VALUE) + 
                ((COMP4 - v2) << (span - 2));
    me.XValue = MASK<TSPAN * 2>::VALUE; 
    me.x += (v2 - me.leftChar) << 1;
    me.leftChar = ordValue(*(it));
    //printf("[debug]::hash %d\n", me.x);
    //v2 = (me.x > 0)?me.hValue:me.crhValue;
    if (me.x > 0)
    {
        v2 =me.hValue; 
        me.strand = 0;
    }
    else 
    {
        v2 = me.crhValue;
        me.strand = 1;
    }
    for (unsigned k = 64-span; k <= 64 - weight; k+=2)
    {
        v1 = v2 << k >> (64-weight);
        if(me.XValue > v1)
        {
            me.XValue=v1;
            t = k;
        }
    } 
    
    me.YValue = (v2 >> (64-t) << (64-t-weight)) +
            (v2 & ((1ULL<<(64-t-weight)) - 1)) + 
            (t << (span - weight - 1));
    //me.YValue = 0;
    return me.XValue; 
}

/*
 * this hashNext function is for index only collect mini hash value [minindex]
 * calculate hValue;
 */ 
template <typename TValue, unsigned TSPAN, unsigned TWEIGHT, typename TSpec, typename TIter>
 typename Value< Shape<TValue, Minimizer<TSPAN, TWEIGHT, TSpec> > >::Type
hashNexth(Shape<TValue, Minimizer<TSPAN, TWEIGHT, TSpec> > &me, TIter const &it)
{
    //typedef typename Size< Shape<TValue, TSpec> >::Type  TSize;
    SEQAN_ASSERT_GT((unsigned)me.span, 0u);
    uint64_t  v2 = ordValue((TValue)*(it + me.span - 1 ));
    me.hValue=((me.hValue & MASK<TSPAN * 2 - 2>::VALUE)<< 2)+ v2;
    me.crhValue=((me.crhValue >> 2) & MASK<TSPAN * 2 - 2>::VALUE) + 
                ((COMP4 - v2) << (TSPAN * 2 - 2));
    me.x += (v2 - me.leftChar) << 1;
    me.leftChar = ordValue(*(it));
    return me.x; 
}
/**
 * only calculate hash value for single strand
 * calculate hValue;
 */ 
template <typename TValue, unsigned TSPAN, unsigned TWEIGHT, typename TSpec, typename TIter>
 typename Value<Shape<TValue, Minimizer<TSPAN, TWEIGHT, TSpec> > >::Type
hashNext_hs(Shape<TValue, Minimizer<TSPAN, TWEIGHT, TSpec> > &me, TIter const &it)
{
    uint64_t v2 = ordValue((TValue)*(it + me.span - 1 ));
    me.hValue=((me.hValue & MASK<TSPAN * 2 - 2>::VALUE)<< 2)+ v2;
    return me.hValue; 
}


template <typename TValue, unsigned TSPAN, unsigned TWEIGHT, typename TSpec, typename TIter>
 typename Value<Shape<TValue, Minimizer<TSPAN, TWEIGHT, TSpec> > >::Type
hashPre_hs(Shape<TValue, Minimizer<TSPAN, TWEIGHT, TSpec> > &me, TIter const &it)
{
    uint64_t v2 = ordValue((TValue)*(it)) << ((TSPAN << 1)  - 2);
    me.hValue=((me.hValue >> 2) & MASK<TSPAN * 2 - 2>::VALUE)+ v2;
    return me.hValue; 
}
/*
 * this hashNext function is for index only collect mini hash value [minindex]
 * calculate hValue;
 */ 
template <typename TValue, unsigned TSPAN, unsigned TWEIGHT, typename TSpec, typename TIter>
 typename Value< Shape<TValue, Minimizer<TSPAN, TWEIGHT, TSpec> > >::Type
hashNextV(Shape<TValue, Minimizer<TSPAN, TWEIGHT, TSpec> > &me, TIter const &it)
{
    //typedef typename Size< Shape<TValue, TSpec> >::Type  TSize;
    SEQAN_ASSERT_GT((unsigned)me.span, 0u);
    uint64_t  v2 = ordValue((TValue)*(it + me.span - 1 ));
    me.hValue=((me.hValue & MASK<TSPAN * 2 - 2>::VALUE)<< 2)+ v2;
    me.crhValue=((me.crhValue >> 2) & MASK<TSPAN * 2 - 2>::VALUE) + 
                ((COMP4 - v2) << (TSPAN * 2 - 2));
    me.x += (v2 - me.leftChar) << 1;
    me.leftChar = ordValue(*(it));
    me.strand = (me.x >> 63) & 1; //Note: me.x type is uint64_t
    return (me.x > 0)?me.hValue:me.crhValue; 
}
/*
 * this hashNext function is for index only collect mini hash value [minindex]
 * calculate XValue
 */ 
template <typename TValue, unsigned TSPAN, unsigned TWEIGHT, typename TSpec, typename TIter>
 typename Value< Shape<TValue, Minimizer<TSPAN, TWEIGHT, TSpec> > >::Type
hashNextX(Shape<TValue, Minimizer<TSPAN, TWEIGHT, TSpec> > &me, TIter const &it)
{

    std::cout << "hashNextxxx\n";
    //typedef typename Size< Shape<TValue, TSpec> >::Type  TSize;
    SEQAN_ASSERT_GT((unsigned)me.span, 0u);
    uint64_t v1;
    unsigned span = TSPAN << 1, weight = TWEIGHT << 1;
    uint64_t t, v2;
    if (me.x > 0)
    {
        v2 =me.hValue; 
        me.strand = 0;
    }
    else 
    {
        v2 = me.crhValue;
        me.strand = 1;
    }
    me.XValue = MASK<TSPAN * 2>::VALUE;
    for (unsigned k = 64-span; k <= 64 - weight; k+=2)
    {
        v1 = v2 << k >> (64-weight);
        if(me.XValue > v1)
        {
            me.XValue=v1;
            t = k;
        }
    } 
    me.YValue = (v2 >> (64-t) << (64-t-weight)) +
                (v2 & ((1ULL<<(64-t-weight)) - 1)) + 
                (t << (span - weight - 1));
                
    (void)it;
    return me.XValue; 
}

template <typename TValue, unsigned TSPAN, unsigned TWEIGHT, typename TSpec>
 uint64_t getT(Shape<TValue, Minimizer<TSPAN, TWEIGHT, TSpec> > & me)
{
    return (me.YValue >> ((TSPAN - TWEIGHT) << 1));
}

template <typename TValue, unsigned TSPAN, unsigned TWEIGHT, typename TSpec>
 uint64_t h2y(Shape<TValue, Minimizer<TSPAN, TWEIGHT, TSpec> > &me,  uint64_t h)
{
    uint64_t x = -1, v1, t=0;
    for (unsigned k = 64-(me.span << 1) ; k <= 64 - (me.weight << 1); k+=2)
    {
        v1 = h << k >> (64-(me.weight<<1));
        if(x > v1)
        { 
            x=v1;
            t=k;
        }
    } 
    return (h>> (64-t) << (64-t-(me.weight<<1)))+(h& (((uint64_t)1<<(64-t-(me.weight<<1))) - 1))+(t<<(((me.span - me.weight) << 1) - 1));

}

template <typename TValue, unsigned TSPAN, unsigned TWEIGHT, typename TSpec>
 uint64_t xy2h(Shape<TValue, Minimizer<TSPAN, TWEIGHT, TSpec> > &me,  uint64_t x, uint64_t y)
{
    
    unsigned span = me.span << 1, weight = TWEIGHT << 1;
    uint64_t t = y >> (span - weight - 1);
    uint64_t t1 = 64 - weight -t;
    uint64_t mask = (1 << t1 ) - 1, mask1 = ((1 << (span - weight)) - 1);
    
    return ((y &(~mask) & mask1)<< weight) + (x << (64 - t - weight)) + (y & mask);
}

template <unsigned TSPAN, unsigned TWEIGHT>
 uint64_t xy2h(uint64_t x, uint64_t y)
{
    
    unsigned span = TSPAN << 1, weight = TWEIGHT << 1;
    uint64_t t = y >> (span - weight - 1);
    uint64_t t1 = 64 - weight -t;
    uint64_t mask = (1 << t1 ) - 1, mask1 = ((1 << (span - weight)) - 1);
    
    return ((y &(~mask) & mask1)<< weight) + (x << (64 - t - weight)) + (y & mask);
}

}

#endif

