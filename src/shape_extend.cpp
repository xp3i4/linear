#include "shape_extend.h"

using namespace seqan;

typedef Dna5 ShapeType;

LShape::LShape(unsigned span):
        span(span),
        weight(span - 8),
        hValue(0),
        crhValue(0),
        XValue(0),
        YValue(0),
        strand(0),
        leftChar(0),
        x(0)
{
    
}

void resize(LShape & me, unsigned new_span, unsigned new_weight)
{   
    me.span = new_span;
    me.weight = new_weight;
}   

static const uint64_t COMP4 = 3;
static const int  ordC = 3;

uint64_t getMask(unsigned bit)
{
    return (1ULL << bit) - 1;
}

uint64_t hashInit(LShape & me, TIterS it)
{
    SEQAN_ASSERT_GT((unsigned)me.span, 0u);

    me.leftChar = 0;
    me.hValue = 0;
    me.crhValue = 0;
    me.leftChar = 0;
    me.x = me.leftChar- ordC;
    uint64_t k =0, count = 0; //COMP for complemnet value;
    while (count < me.span)
    {
        if (ordValue((ShapeType)*(it + k + count)) == 4)
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
        uint64_t val = ordValue ((ShapeType)*(it + k + i));
        me.x += (val << 1) - ordC;
        me.hValue = (me.hValue << 2) + val;
        me.crhValue += ((COMP4 - val) << bit);
        bit += 2;
    }
    return k;
}
/**
 *init for hashNext_hs and hashPre_hs
 *@d=0 init for hashNext_hs
 *@d=1 init for hashPre_hs 
 */
uint64_t hashInit_hs(LShape & me, TIterS it, int d)
{
    me.hValue = 0;
    for (unsigned i = d; i < me.span - 1 + d; ++i)
    {
        me.hValue = (me.hValue << 2) + ordValue ((ShapeType)*(it + i));
    }
    me.hValue <<= (d << 1);
    return 0;
}
uint64_t hashNext(LShape & me, TIterS it)
{
    SEQAN_ASSERT_GT((unsigned)me.span, 0u);
    uint64_t v1;
    unsigned t, span = me.span << 1, weight = me.weight << 1;
    uint64_t v2 = ordValue((ShapeType)*(it + me.span - 1));
    uint64_t mask = getMask(span - 2);
    me.hValue=((me.hValue & mask)<<2)+ v2;
    me.crhValue=((me.crhValue >> 2) & mask) + 
                ((COMP4 - v2) << (span - 2));
    me.XValue = getMask(span); 
    me.x += (v2 - me.leftChar) << 1;
    me.leftChar = ordValue((ShapeType)*(it));
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
    me.YValue = (v2 >> (64-t) << (64 - t - weight)) +
            (v2 & ((1ULL<<(64 - t - weight)) - 1)) + 
            (t << (span - weight - 1));
    return me.XValue; 
}
/*
 * this hashNext function is for index only collect mini hash value [minindex]
 * calculate hValue;
 */ 
uint64_t hashNexth(LShape & me, TIterS it)
{
    SEQAN_ASSERT_GT((unsigned)me.span, 0u);
    uint64_t v2 = ordValue((ShapeType)*(it + me.span - 1 ));
    uint64_t mask = getMask((me.span << 1) - 2);
    me.hValue = ((me.hValue & mask) << 2) + v2;
    me.crhValue = ((me.crhValue >> 2) & mask) + 
                  ((COMP4 - v2) << ((me.span << 1) - 2));
    me.x += (v2 - me.leftChar) << 1;
    me.leftChar = ordValue((ShapeType)*(it));
    return me.x; 
}
/**
 * only calculate hash value for single strand
 * calculate hValue;
 */ 
uint64_t hashNext_hs(LShape & me, TIterS it)
{
    uint64_t v2 = ordValue((ShapeType)*(it + me.span - 1 ));
    uint64_t mask = getMask((me.span << 1) - 2);
    me.hValue = ((me.hValue & mask) << 2) + v2;
    return me.hValue; 
}
uint64_t hashPre_hs(LShape & me, TIterS it)
{
    uint64_t v2 = ordValue((ShapeType)*(it)) << ((me.span << 1)  - 2);
    uint64_t mask = getMask((me.span << 1) - 2);
    me.hValue = ((me.hValue >> 2) & mask) + v2;
    return me.hValue; 
}
/*
 * this hashNext function is for index only collect mini hash value [minindex]
 * calculate hValue;
 */ 
uint64_t hashNextV(LShape & me, TIterS it)
{
    SEQAN_ASSERT_GT((unsigned)me.span, 0u);
    uint64_t v2 = ordValue((ShapeType)*(it + me.span - 1 ));
    uint64_t mask = getMask((me.span << 1) - 2);
    me.hValue=((me.hValue & mask) << 2) + v2;
    me.crhValue=((me.crhValue >> 2) & mask) + 
                ((COMP4 - v2) << ((me.span <<1) - 2));
    me.x += (v2 - me.leftChar) << 1;
    me.leftChar = ordValue((ShapeType)*(it));
    me.strand = (me.x >> 63) & 1; //Note: me.x type is uint64_t
    return (me.x > 0)?me.hValue:me.crhValue; 
}
/*
 * this hashNext function is for index only collect mini hash value [minindex]
 * calculate XValue
 */ 
uint64_t hashNextX(LShape & me, TIterS it)
{
    SEQAN_ASSERT_GT((unsigned)me.span, 0u);
    uint64_t v1;
    unsigned span = me.span << 1, weight = me.weight << 1;
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
    me.XValue = getMask(me.span << 1);
    for (unsigned k = 64 - span; k <= 64 - weight; k += 2)
    {
        v1 = v2 << k >> (64 - weight);
        if(me.XValue > v1)
        {
            me.XValue = v1;
            t = k;
        }
    } 
    me.YValue = (v2 >> (64 - t) << (64 - t - weight)) +
                (v2 & ((1ULL<<(64 - t - weight)) - 1)) + 
                (t << (span - weight - 1));
    (void)it;
    return me.XValue; 
}


