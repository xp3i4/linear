#ifndef LINEAR_HEADER_SHAPE_EXTEND_H
#define LINEAR_HEADER_SHAPE_EXTEND_H
#include <seqan/sequence.h>
using namespace seqan;

typedef typename Iterator<String<Dna5> >::Type TIterS;
class LShape
{
public:
    unsigned span;
    unsigned weight;
    uint64_t hValue;    //hashValue
    uint64_t crhValue;  //inverse hash
    uint64_t XValue;    
    uint64_t YValue;
    uint64_t strand;
    int     leftChar;
    int     x;

    void init_shape_parm (unsigned shape_span);
    LShape(unsigned span);
};

void resize(LShape & me, unsigned new_span, unsigned new_weight);
uint64_t getMask(unsigned bit);
uint64_t hashInit(LShape & me, TIterS it);
uint64_t hashInit_hs(LShape & me, TIterS it, int d);
uint64_t hashNext(LShape & me, TIterS it);
uint64_t hashNexth(LShape & me, TIterS it);
uint64_t hashNext_hs(LShape & me, TIterS it);
uint64_t hashPre_hs(LShape & me, TIterS it);
uint64_t hashNextV(LShape & me, TIterS it);
uint64_t hashNextX(LShape & me, TIterS it);

#endif

