#ifndef SEQAN_HEADER_PMP_FINDER_H
#define SEQAN_HEADER_PMP_FINDER_H
#include "base.h"
#include "index_util.h"

using namespace seqan;

typedef Iterator <String <Dna5> >::Type TIter5;

struct CordBase
{
    //Cord(C): coordinates of the vertex of sliding windows
    //=|N/A[2]|strand[1]|cordEnd[1] genomeCord [40] |readCord [20bits]
    //cell [4] is the minimum length the window is allowed to slide in the alignment matrix.
    //genomeCord(gC or xC): = position in the genome >> cell_bit << cell_bit. the last cell_bit bits maybe set to 0
    //gC:= SA node = Seq num i1 [10] | Base num i2 [30]  
    //readCord(rC or yC): ~= position in the read >> cell_bit << cell_bit. the last cell_bit bits maybe set to 0 during process.
    //rC:= Base num [20]
    
    typedef unsigned Bit;
    typedef uint64_t Mask;
    typedef uint64_t Flag;
    typedef unsigned Size;
    
    Bit bit;
    uint64_t flagEnd;
    Mask mask;
    Mask maskx;
    Mask valueMask;
    Bit flag_bit;
    Flag flag_strand;
    Flag flag_end;
    Bit cell_bit;
    Size cell_size;
    Mask headFlag;
    Mask valueMask_dstr;
    unsigned bit_id;
    CordBase();
};
extern CordBase _DefaultCordBase;

struct Cord
{
    uint64_t getCordX(uint64_t const &, unsigned const & = _DefaultCordBase.bit, uint64_t const & = _DefaultCordBase.maskx) const;
    uint64_t getCordY(uint64_t const & cord, uint64_t const & mask = _DefaultCordBase.mask) const;
    uint64_t createCord(uint64_t const & x, 
                 uint64_t const & y, 
                 uint64_t const & strand,
                 unsigned const & bit = _DefaultCordBase.bit, 
                 unsigned const & bit2 = _DefaultCordBase.flag_bit) const;
    uint64_t hit2Cord(uint64_t const & hit, 
               unsigned const & bit = _DefaultCordBase.bit, 
               uint64_t const & mask = _DefaultCordBase.mask,
               uint64_t const & mask2 = _DefaultCordBase.valueMask
              ) const;
    uint64_t hit2Cord_dstr(uint64_t const & hit, 
               unsigned const & bit = _DefaultCordBase.bit, 
               uint64_t const & mask = _DefaultCordBase.mask,
               uint64_t const & mask2 = _DefaultCordBase.valueMask_dstr
              ) const;
    uint64_t cord2Cell(uint64_t const & cord, 
                unsigned const & bit = _DefaultCordBase.cell_bit) const;
    uint64_t cell2Cord(uint64_t const & cell, 
                unsigned const & bit = _DefaultCordBase.cell_bit) const;
    void setCordEnd(uint64_t & cord,
            typename CordBase::Flag const & strand = _DefaultCordBase.flag_strand,
            typename CordBase::Flag const & end = _DefaultCordBase.flag_end);
    uint64_t getCordStrand(uint64_t const & cord,
            unsigned const & strand = _DefaultCordBase.flag_bit) const;
    uint64_t isCordEnd(uint64_t const & cord,
                typename CordBase::Flag const & end = _DefaultCordBase.flag_end) const;
    void setMaxLen(String<uint64_t> &, uint64_t const &, uint64_t const & = _DefaultCordBase.mask);
    uint64_t getMaxLen(String<uint64_t> const &, uint64_t const & = _DefaultCordBase.mask);
    uint64_t shift(uint64_t const & val, int64_t x, int64_t y, unsigned const & = _DefaultCordBase.bit); //add x and y

    bool isCordsOverlap(uint64_t & val1, uint64_t & val2, int64_t thd);
    bool isBlockEnd(uint64_t &, uint64_t const & = _DefaultCordBase.flagEnd);
};
extern Cord _DefaultCord; 
uint64_t get_cord_x (uint64_t val);
uint64_t get_cord_y (uint64_t val); 
uint64_t get_cord_strand (uint64_t val);
uint64_t get_cord_id (uint64_t val);
uint64_t create_id_x (uint64_t id, uint64_t x);
uint64_t create_cord (uint64_t id, uint64_t cordx, uint64_t cordy, uint64_t strand);

void cmpRevCord(uint64_t val1, uint64_t val2, uint64_t & cr_val1, uint64_t & cr_val2, uint64_t read_len);
uint64_t set_cord_xy (uint64_t val, uint64_t x, uint64_t y);

unsigned _windowDist(Iterator<String<short> >::Type const & it1, 
                        Iterator<String<short> >::Type const & it2);
//WARNING:The length of read should be < 1MB;
extern const float band_width;
extern const unsigned cmask;
extern const unsigned cell_size;
extern const unsigned cell_num;
extern const unsigned window_size; //16*12
extern const unsigned window_delta;
extern const unsigned sup;
extern const unsigned med;
extern const unsigned inf;

extern const unsigned initx; 
extern const unsigned inity;

extern const unsigned scriptStep;
extern const unsigned scriptBit;
extern const unsigned scriptWindow; //script_length = 2^scriptWindow
extern const unsigned scriptWindow2;
extern const unsigned scriptWindow3;
extern const int scriptCount[5];
extern const int scriptMask;
extern const int scriptMask2;
extern const int scriptMask3;

extern const uint64_t hmask;
/**
 * ATTENTION TODO parameter needs tuning: will affect speed, gap extension, clip
 */
extern const unsigned windowThreshold; // 36;
/**
 *   struct hit:
 *   extend the structure Cord;
 *   NA[2]|strand[1]|head[1]|genome pos[40]|read pos[20]
 *   NodeType: 1 Head, 0 Body
 */
void createFeatures(TIter5 const &, TIter5 const &, String<short> & );
void createFeatures(StringSet<String<Dna5> > &, StringSet<String<short> > &, unsigned);


struct HitBase
{
    uint64_t bit;
    uint64_t bit2;
    uint64_t flag;
    uint64_t flag2;
    uint64_t mask;
    HitBase();
};
extern HitBase _DefaultHitBase;

struct Hit
{
    void setBlockStart(uint64_t &, uint64_t const & = _DefaultHitBase.flag);
    void setBlockBody(uint64_t &, uint64_t const & = _DefaultHitBase.flag);
    bool isBlockStart(uint64_t &, uint64_t const & = _DefaultHitBase.flag);
    void setBlockEnd(uint64_t &, uint64_t const & = _DefaultHitBase.flag);
    void unsetBlockEnd(uint64_t &, uint64_t const & = _DefaultHitBase.flag);
    void setBlockStrand(uint64_t &, uint64_t const &, 
                     uint64_t const & = _DefaultHitBase.flag2);
    bool isBlockEnd(uint64_t &, uint64_t const & = _DefaultHitBase.flag);
    unsigned getStrand(uint64_t const &, uint64_t const & = _DefaultHitBase.flag2);
};
extern Hit _DefaultHit;

uint64_t mnMapReadList( LIndex  & index,
                        String<Dna5> & read,
                        Anchors & anchors,
                        MapParm & mapParm,
                        String<uint64_t> & hit);
bool path_dst(
             typename Iterator<String<uint64_t> >::Type hitBegin, 
             typename Iterator<String<uint64_t> >::Type hitEnd, 
             StringSet<String<short> > & f1,
             StringSet<String<short> > & f2, 
             String<uint64_t> & cords,
             float const & cordLenThr
            );

int extendPatch(StringSet<String<short> > & f1, 
                StringSet<String<short> > & f2, 
                  String<uint64_t> & cords,
                  int k,
                  uint64_t cord1,
                  uint64_t cord2,
                  int revscomp_const
                );

int getsai1(int k);
#endif