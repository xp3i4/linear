#include <seqan/seq_io.h>
#include "shape_extend.h"

using namespace seqan;

//WARN!!:: Only odd Shape len is allowed if call the hash due to the double strand hash value
//All even shape len will be converted to len + 1
typedef Dna5 ShapeType;

int lexicoHash2kmer(uint64_t val, uint64_t k, String<char> & kmer)
{
    char s1[5] = {'A','C','G','T','N'};
    clear(kmer);
    for (uint64_t i = 0; i < k; i++)
    {
        uint64_t vi = (val >> (2 * (k - i - 1))) & (3ULL);
        //std::cout << "vi" << vi << s1[vi]<< k << " " << val <<"\n";
        appendValue(kmer, s1[vi]);
    }
    return 0;
}
int print_minimizer(Iterator<String<Dna5> >::Type it, uint64_t x, uint64_t y, uint64_t strand, uint64_t span, uint64_t weight, String<char> header)
{
    char s1[5] = {'A','C','G','T','N'};
    char s2[5] = {'T','G','C','A','N'};
    String<char> kmer;
    String<char> pre;
    String<char> suf;
    String<char> m1;
    String<char> m2;
    lexicoHash2kmer(x, weight, m1);
    lexicoHash2kmer(y, 4, m2);
    for (uint64_t i = 0; i < span + 4; i++)
    {
        if (strand == 0) 
        {
            if (i == span)
            {
                appendValue(kmer,'|');
            }
            appendValue(kmer,s1[seqan::ordValue((Dna5)*(it + i))]);
        }
        else if (strand == 1)
        {
            if (i == span)
            {
                appendValue(kmer,'|');
            }
            appendValue(kmer,s2[seqan::ordValue((Dna5)*(it + span - 1 - i ))]);
        }
    }
    std::cout << header << " " << strand << " " << kmer << " " << m1 << " " << m2 << " " << y << " | "<< "\n";
    return 0;
}
void LShape::init_shape_parm (unsigned shape_span)
{
    span = shape_span ;
    weight = span - 8;
}
LShape::LShape(unsigned shape_span):
        hValue(0),
        crhValue(0),
        XValue(0),
        YValue(0),
        strand(0),
        leftChar(0),
        x(0)
{
    init_shape_parm(shape_span);
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
    me.x = me.leftChar - ordC;
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
        me.x += (int(val) << 1) - ordC;
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
    unsigned t = 0, span = me.span << 1, weight = me.weight << 1;
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
    for (unsigned k = 64 - span; k <= 64 - weight; k += 2)
    {
        v1 = v2 << k >> (64-weight);
        if(me.XValue > v1)
        {
            me.XValue = v1;
            t = k;
        }
    } 
    me.YValue = (v2 >> (64 - t) << (64 - t - weight)) +
            (v2 & ((1ULL << (64 - t - weight)) - 1)) + 
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
    uint64_t mask = getMask((me.span << 1) - 2);
    int v2 = ordValue((ShapeType)*(it + me.span - 1 ));
    me.hValue = ((me.hValue & mask) << 2) + v2;
    me.crhValue = ((me.crhValue >> 2) & mask) + 
                  ((COMP4 - v2) << ((me.span << 1) - 2));
    me.x += (v2 - me.leftChar) << 1;
    me.leftChar = ordValue((ShapeType)*(it));
    return me.x < 0 ? me.hValue : me.crhValue;
}

/*
uint64_t hashNexth_hpc(LShape & me, TIterS & it, TIterS it_end)
{
    SEQAN_ASSERT_GT((unsigned)me.span, 0u);
    uint64_t mask = getMask((me.span << 1) - 2);
    int v2 = ordValue((ShapeType)*(it + me.span - 1 ));
    me.hValue = ((me.hValue & mask) << 2) + v2;
    me.crhValue = ((me.crhValue >> 2) & mask) + 
                  ((COMP4 - v2) << ((me.span << 1) - 2));
    me.x += (v2 - me.leftChar) << 1;
    me.leftChar = ordValue((ShapeType)*(it));
    TIterS it_origin = it + me.span - 1
    TIterS it_next = it + me.spand;

    while (it_next < it_end && *it_next == *it_origin)
    {
        ++it_next;
    }
    it = it_next < it_end ? it_next - me.span + 1 : it_end - me.span + 1;
    return me.x < 0 ? me.hValue : me.crhValue;
}
*/

/**
 *  calculate hash value for single strand
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
    uint64_t mask = getMask((me.span << 1) - 2);
    int v2 = ordValue((ShapeType)*(it + me.span - 1 ));
    me.hValue=((me.hValue & mask) << 2) + v2;
    me.crhValue=((me.crhValue >> 2) & mask) + 
                ((COMP4 - v2) << ((me.span <<1) - 2));
    me.x += (v2 - me.leftChar) << 1;
    me.leftChar = ordValue((ShapeType)*(it));
    me.strand = me.x < 0 ? 1 : 0;
    return me.strand ? me.crhValue : me.hValue;
}

uint64_t hashNextXX(LShape & me, TIterS it, uint64_t & v2, uint64_t & t)
{
    uint64_t v1;
    unsigned span = me.span << 1, weight = me.weight << 1;
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
    (void)it;
    return me.XValue;
}
uint64_t hashNextXY(LShape & me, TIterS it, uint64_t & v2, uint64_t & t)
{
    unsigned span = me.span << 1, weight = me.weight << 1;
    me.YValue = (v2 >> (64 - t) << (64 - t - weight)) +
                (v2 & ((1ULL<<(64 - t - weight)) - 1)) + 
                (t << (span - weight - 1));
    (void)it;

    return me.YValue;
}
uint64_t hashNextXY2(LShape & me, TIterS it, uint64_t & v2, uint64_t & t)
{
    me.YValue = 0;
    int64_t n = 4;
    //<<debug
    String<char> kmer1;
    String<char> kmer2;
    //>>debu
    if (me.x > 0)
    {
        int64_t d_it = (t >> 1) + me.span + me.weight - 32;
//        for (int64_t i = d_it; i < d_it + (me.span - me.weight); i++)
        for (int64_t i = d_it; i < d_it + n; i++)
        {
            int64_t val = ordValue((ShapeType)*(it + i));
            me.YValue = val > 3 ? (me.YValue << 2) : (me.YValue << 2) + val;
            //<<debug
            lexicoHash2kmer(me.XValue, me.weight, kmer1);
            lexicoHash2kmer(me.YValue, n, kmer2);
            //std::cout << "hnxy2 " << kmer1 << " " << kmer2 << " " << t << " " << d_it << "\n";
            //>>debug
        }
    }
    else
    {
        int64_t d_it = -(t >> 1) - me.weight + 31;
        //for (int64_t i = d_it; i < d_it + (me.span - me.weight); i++)
        //std::cout << "hnxy35 " << d_it << " " << d_it - n <<  "\n";
        for (int64_t i = d_it; i > d_it - n; i--)
        {
            int64_t val = COMP4 - ordValue((ShapeType)*(it + i));
 //           std::cout << "hnxy34 " << me.YValue << "\n"; 
            if (val < 0)
            {
                me.YValue = me.YValue << 2;
  //              std::cout << "hnxy33 " << me.YValue << " " << *(it + i) << "\n";
            }
            else 
            {
                me.YValue = (me.YValue << 2) + val;
   //             std::cout << "hnxy32 " << me.YValue << " " << *(it + i) << "\n";
            }
            //me.YValue = val < 0 ? (me.YValue << 2) : (me.YValue << 2) + val;
        } 
         //std::cout << "hnxy3 " << me.YValue <<" " << d_it << *(it + d_it) << *(it -1+d_it) << *(it -2+d_it) << *(it -3+d_it) << "\n";
    }
    /*
    if (me.YValue > 4)
    {
        std::cout << "hashNextXY2 " << me.YValue << "\n";
    }
    */
    (void) v2;
    return me.YValue;
}
/*
 * this hashNext function is for index only collect mini hash value [minindex]
 * calculate XValue
 */ 
uint64_t hashNextX(LShape & me, TIterS it)
{
    SEQAN_ASSERT_GT((unsigned)me.span, 0u);
    uint64_t t = 0, v2;
    hashNextXX(me, it, v2, t);
    hashNextXY2(me, it, v2, t);
    return me.XValue; 
}

uint64_t hashNextX2(LShape & me, TIterS it)
{
    SEQAN_ASSERT_GT((unsigned)me.span, 0u);
    uint64_t t = 0, v2;
    hashNextXX(me, it, v2, t);
    hashNextXY2(me, it, v2, t);
    return me.XValue; 
}
