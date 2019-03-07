#ifndef SEQAN_HEADER_PMP_FINDER_H
#define SEQAN_HEADER_PMP_FINDER_H

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

}_DefaultCordBase;

struct Cord
{
    uint64_t getCordX(uint64_t const &, unsigned const &, uint64_t const &) const;
    uint64_t getCordY(uint64_t const &, uint64_t const &) const;
    uint64_t createCord(uint64_t const &, uint64_t const &, uint64_t const &, unsigned const &, unsigned const &) const ;
    uint64_t hit2Cord(uint64_t const &, unsigned const &, uint64_t const &, uint64_t const &) const;
    uint64_t hit2Cord_dstr(uint64_t const &, unsigned const &, uint64_t const &, uint64_t const &) const;
    uint64_t cord2Cell(uint64_t const &, unsigned const &) const;
    uint64_t cell2Cord(uint64_t const &, unsigned const &) const;
    void setCordEnd(uint64_t &, uint64_t const &, uint64_t const &);
    uint64_t getCordStrand(uint64_t const &, unsigned const &) const;
    uint64_t isCordEnd(uint64_t const &, uint64_t const &)const;
    void setMaxLen(String<uint64_t> &, uint64_t const &, uint64_t const & = _DefaultCordBase.mask);
    uint64_t getMaxLen(String<uint64_t> const &, uint64_t const & = _DefaultCordBase.mask);
    uint64_t shift(uint64_t const & val, int64_t x, int64_t y, unsigned const & = _DefaultCordBase.bit); //add x and y

    bool isCordsOverlap(uint64_t & val1, uint64_t & val2, int64_t thd);
    bool isBlockEnd(uint64_t &, uint64_t const & = _DefaultCordBase.flagEnd);
}_DefaultCord; 
uint64_t get_cord_x (uint64_t val);
uint64_t get_cord_y (uint64_t val); 
uint64_t get_cord_strand (uint64_t val);
uint64_t get_cord_id (uint64_t val);
uint64_t create_id_x (uint64_t id, uint64_t x);
uint64_t create_cord (uint64_t id, uint64_t cordx, uint64_t cordy, uint64_t strand);

void cmpRevCord(uint64_t val1, uint64_t val2, uint64_t & cr_val1, uint64_t & cr_val2, uint64_t read_len);
uint64_t set_cord_xy (uint64_t val, uint64_t x, uint64_t y);

//WARNING:The length of read should be < 1MB;
static const float band_width = 0.25;
static const unsigned cmask = ((uint64_t)1<<20) - 1;
static const unsigned cell_size = 16;
static const unsigned cell_num = 12;
static const unsigned window_size = cell_size * cell_num; //16*12
static const unsigned window_delta = window_size * (1 - 2 * band_width);
static const unsigned sup = cell_num;
static const unsigned med =ceil((1 - band_width) * cell_num);
static const unsigned inf = ceil((1 - 2 * band_width) * cell_num);

static const unsigned initx = 5; 
static const unsigned inity = 5;

//======================================================
static const unsigned scriptStep=16;
static const unsigned scriptBit=4;
static const unsigned scriptWindow=5; //script_length = 2^scriptWindow
static const unsigned scriptWindow2 = scriptWindow << 1;
static const unsigned scriptWindow3 = scriptWindow2 + scriptWindow;
static const int scriptCount[5] = {1, 1<<scriptWindow, 1 <<(scriptWindow * 2), 0, 0};
static const int scriptMask = (1 << scriptWindow) - 1;
static const int scriptMask2 = scriptMask << scriptWindow;
static const int scriptMask3 = scriptMask2 << scriptWindow;

static const uint64_t hmask = (1ULL << 20) - 1;
/**
 * ATTENTION TODO parameter needs tuning: will affect speed, gap extension, clip
 */
static const unsigned windowThreshold = 36; // 36;
/**
 *   struct hit:
 *   extend the structure Cord;
 *   NA[2]|strand[1]|head[1]|genome pos[40]|read pos[20]
 *   NodeType: 1 Head, 0 Body
 */

struct HitBase
{
    uint64_t bit;
    uint64_t bit2;
    uint64_t flag;
    uint64_t flag2;
    uint64_t mask;
    
    HitBase():
        bit(60),
        bit2(61),
        flag(1ULL<<bit),
        flag2(1ULL<<bit2),
        mask(flag - 1)
        {}
}_DefaultHitBase;

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

}_DefaultHit;

#endif