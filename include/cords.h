#ifndef LINEAR_HEADER_CORDS_H
#define LINEAR_HEADER_CORDS_H

#include <seqan/sequence.h>

using seqan::String;
using seqan::StringSet;
using seqan::CharString;
using namespace seqan;
/* @anchorX = x - y + g_hs_anchor_zero. g_hs_anchor_zero to restrict @anchorX > 0. 
   (bits overflow otherwise) such that -g_hs_anchor_zero <= x - y < g_hs_anchor_zero
*/
extern uint64_t const_anchor_zero;
extern uint64_t const_cordx_max; 
extern uint64_t FORWARD_STRAND;
extern uint64_t REVERSE_STRAND; 
extern uint64_t INFI_CORD;
extern uint64_t MAX_CORD_ID;
extern uint64_t MAX_CORD_X ;
extern uint64_t MAX_CORD_Y ;
extern bool     f_debug;

/*
 * Cord(C): pair of bitwised coordinates of two sequences;
 * :=|main[1]|record[1]|strand[1]|blockEnd[1]|gC [40] |rC [20]
 * @main: cords that generated by the apx mapping. They are the main chain of the 
    mapping. Other cords are gap cords that generated by the gap mapping extened from the main cords.
    1:main, 0:gap 
 * @Record: a groub of blocks from the same template (alignment/mapping). 
    One record is composed of main_cords + gap_cords.
    In the String of cords, adjacent records are marked by different values as 
    00000000|11111111
    records1|records2 
 * @genomeCord(gC or xval):= position in the genome >> cell_bit << cell_bit. the 
    last ll_bit bits maybe set to 0
 * @blockEnd: mark end of block which is a group cords of same anchors.
 * @gC(genome coordinate):= SA node := Seq_id[10] | base_x i2[30]  
 * @rC:= base_y[20]. Read coordinate
 */
typedef uint64_t CordType;
typedef String<CordType> CordsType;
typedef StringSet<CordsType > CordsSetType;
struct CordsParms
{
    int64_t const thd_rcdxdy_min_dx = -20;
    int64_t const thd_rcdxdy_min_dy = -20;
    int64_t const thd_rcdxdy_merge_min_dx = 10;
    int64_t const thd_rcdxdy_merge_min_dy = 10;

    CordsParms();
};
struct CordBase
{
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
    uint64_t f_main;
    uint64_t f_recd;
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
    uint64_t get_hit_strx(uint64_t const & hit, 
               unsigned const & bit = _DefaultCordBase.bit, 
               uint64_t const & mask = _DefaultCordBase.mask,
               uint64_t const & mask2 = _DefaultCordBase.valueMask_dstr
              ) const;
    uint64_t cord2Cell(uint64_t const & cord, 
                unsigned const & bit = _DefaultCordBase.cell_bit) const;
    uint64_t cell2Cord(uint64_t const & cell, 
                unsigned const & bit = _DefaultCordBase.cell_bit) const;
    void setCordEnd(uint64_t & cord,
            typename CordBase::Flag const & end = _DefaultCordBase.flag_end);
    uint64_t getCordStrand(uint64_t const & cord,
            unsigned const & strand = _DefaultCordBase.flag_bit) const;
    void setMaxLen(String<uint64_t> &, uint64_t const &, uint64_t const & = _DefaultCordBase.mask);
    uint64_t getMaxLen(String<uint64_t> const &, uint64_t const & = _DefaultCordBase.mask);
    uint64_t shift(uint64_t const & val, int64_t x, int64_t y, unsigned const & = _DefaultCordBase.bit); //add x and y
    bool isCordsOverlap(uint64_t & val1, uint64_t & val2, int64_t thd);
    bool isBlockEnd(uint64_t &, uint64_t const & = _DefaultCordBase.flagEnd);
    uint64_t makeBlockEndVal(uint64_t, uint64_t const & = _DefaultCordBase.flagEnd);
};
extern Cord _DefaultCord; 
int initCords (String<uint64_t> &);
int initHits (String<uint64_t> &);
int initHitsScore (String<int> & hit_score);
int isHitsEmpty(String<uint64_t> & hits);
int isFirstHit(Iterator<String<uint64_t> >::Type & it );
int isLastHit(Iterator<String<uint64_t> >::Type & it);
Iterator<String<uint64_t> >::Type beginHits (String<uint64_t> & hits);
Iterator<String<uint64_t> >::Type endHits(String<uint64_t> & hits);
/**
 *  struct hit:
 *  extend the structure Cord;
 *  NA[1]|Long pattern[1]|strand[1]|head[1]|genome pos[40]|read pos[20]
 *  NodeType: 1 Head, 0 Body
 *  long pattern: the matched pattern is a long pattern
 */
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
    void setBlockStrand(uint64_t &, uint64_t const &, uint64_t const & = _DefaultHitBase.flag2);
    bool isBlockEnd(uint64_t &, uint64_t const & = _DefaultHitBase.flag);
    unsigned getStrand(uint64_t const &, uint64_t const & = _DefaultHitBase.flag2);
    uint64_t getAnchor(uint64_t const &);
    void setLongPattern(uint64_t &);
    void unsetLongPattern(uint64_t &);
    uint64_t isLongPattern(uint64_t &);
};
extern Hit _DefaultHit;

void cmpRevCord(uint64_t, uint64_t, uint64_t &, uint64_t &, uint64_t);
uint64_t get_cord_x (uint64_t);
uint64_t get_cord_y (uint64_t); 
uint64_t get_cord_strand (uint64_t);
uint64_t get_cord_id (uint64_t);
uint64_t get_cord_recd (uint64_t cord);
uint64_t get_cord_xy (uint64_t val);
uint64_t shift_cord(uint64_t const &, int64_t, int64_t);
uint64_t create_id_x (uint64_t, uint64_t);
uint64_t create_cord (uint64_t, uint64_t, uint64_t, uint64_t);
uint64_t new_xy_cord(uint64_t, uint64_t, uint64_t);
void set_cord_xy (uint64_t & val, uint64_t x, uint64_t y);
void set_cord_id (uint64_t & val, uint64_t id);
void set_cord_end (uint64_t &); 
void unset_cord_end(uint64_t &);
void set_cord_block_end(uint64_t & val);
void print_cord(uint64_t, CharString = "");
void print_cords(String<uint64_t> &, CharString = "", bool f_print = f_debug);
void set_cord_main (uint64_t & cord);
void set_cord_gap (uint64_t & cord);
void set_cord_recd(uint64_t & cord, uint64_t sgn);
uint64_t is_cord_main(uint64_t cord);
int64_t atomic_inc_cord_y (uint64_t & cord); // atomic cord++, return the new cord
uint64_t is_cord_block_end(uint64_t);
int isCordsConsecutive_(uint64_t & cord1, uint64_t cord2, uint64_t thd_cord_gap);
uint64_t make_anchor(uint64_t id, uint64_t x, uint64_t y, uint64_t strand);

//return if range [x11, x12) and [x21, x22) has overlap
bool _isRangeOverLap(uint64_t x11, uint64_t x12, uint64_t x21, uint64_t x22);
bool _isCordyOverLap(uint64_t cord11, uint64_t cord12, uint64_t cord21, uint64_t cord22, uint64_t read_len);
uint64_t isDiffCordsStrand(uint64_t & cord1, uint64_t & cord2);


uint64_t getAnchorX(uint64_t anchor);

typedef std::pair<uint64_t, uint64_t> UPair;
UPair getUPForwardy(UPair str_end, uint64_t read_len);

void printAnchors(String<uint64_t> & anchors, CharString header = "", bool f_print = f_debug);
int reformCords(String<uint64_t> & cords_str,
                String<uint64_t> & cords_end,
                int (*reformCordFunc) (String<uint64_t> &, String<uint64_t> &, String<int> &,
                    String<int> &, String<int> &, String<int> &, unsigned &, CordsParms &),
                CordsParms & cords_parms);

int reformCordsDxDy1(String<uint64_t> & cords_str,
                     String<uint64_t> & cords_end,
                     String<int> & bands11,
                     String<int> & bands12,
                     String<int> & bands21,
                     String<int> & bands22,
                     unsigned & it,
                     CordsParms & cords_parms);

struct CordInfo
{
    float score; 

    CordInfo(float score);
};
#endif