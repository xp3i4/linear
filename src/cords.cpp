#include <seqan/stream.h>
#include <seqan/parallel.h>
#include "cords.h"
 
using namespace seqan;

uint64_t const_anchor_zero = (1ULL << 20); // make sure y in cord < this  
uint64_t const_cordx_max = (1ULL << 30) - 1; //cordx <= this
uint64_t FORWARD_STRAND = 0;
uint64_t REVERSE_STRAND = 1; 
uint64_t INFI_CORD = shift_cord (0, (1ULL << 40) - 1,  (1ULL << 20) - 1);
uint64_t MAX_CORD_ID = (1ULL << 10) - 1;
uint64_t MAX_CORD_X  = (1ULL << 30) - 1;
uint64_t MAX_CORD_Y  = (1ULL << 20) - 1;
 
CordBase::CordBase():
        bit(20),
        flagEnd(1ULL << 60),
        mask(0xfffff),
        maskx(0xffffffffff),
        valueMask((1ULL<< 60) - 1),
        flag_bit(61),
        flag_strand(1ULL << flag_bit),
        flag_end(0x1000000000000000),
        cell_bit(4),
        cell_size(16),
        headFlag((1ULL<<63)),
        valueMask_dstr(valueMask | flag_strand), 
        bit_id (30),
        f_main(1ULL << 63),
        f_recd(1ULL << 62)
{}
CordBase _DefaultCordBase;  
Cord _DefaultCord;
HitBase::HitBase():
        bit(60),
        bit2(61),
        flag(1ULL<<bit),
        flag2(1ULL<<bit2),
        mask(flag - 1)
{}
HitBase _DefaultHitBase;
Hit _DefaultHit;
//TODO::re-wrapper the constatns!.
uint64_t Cord::getCordX(uint64_t const & cord, 
               unsigned const & bit,
               uint64_t const & mask) const
{
    return (cord >> bit) & mask; 
}
uint64_t Cord::getCordY(uint64_t const & cord, 
               uint64_t const & mask) const 
{
    return cord & mask;
}
uint64_t Cord::createCord(uint64_t const & x, 
                 uint64_t const & y, 
                 uint64_t const & strand,
                 unsigned const & bit, 
                 unsigned const & bit2) const
{
    return (x << bit) + y + (strand << bit2);
}

uint64_t Cord::hit2Cord(uint64_t const & hit, 
               unsigned const & bit, 
               uint64_t const & mask,
               uint64_t const & mask2
              ) const
{
    uint64_t new_cord = (hit + ((hit & mask) << bit) - (const_anchor_zero << bit)) & mask2;
    _DefaultHit.unsetLongPattern(new_cord);
    return new_cord;
}

uint64_t Cord::hit2Cord_dstr(uint64_t const & hit, 
               unsigned const & bit, 
               uint64_t const & mask,
               uint64_t const & mask2
              ) const
{
    uint64_t new_cord = (hit + ((hit & mask) << bit)  - (const_anchor_zero << bit)) & mask2;
    _DefaultHit.unsetLongPattern(new_cord);
    return new_cord;
}

uint64_t Cord::get_hit_strx(uint64_t const & hit, unsigned const & bit, 
  uint64_t const & mask1, uint64_t const & mask2) const
{ 
  return ((hit & mask2 ) >> bit) + (hit & mask1) - const_anchor_zero;
}

uint64_t Cord::cord2Cell(uint64_t const & cord, 
                unsigned const & bit) const
{
    return cord >> bit;
}

 uint64_t Cord::cell2Cord(uint64_t const & cell, 
                unsigned const & bit) const
{
    return cell << bit;
}

void Cord::setCordEnd(uint64_t & cord,
            typename CordBase::Flag const & end) //set cord block end flag
{
    cord |= end;
}

typename CordBase::Flag Cord::getCordStrand(uint64_t const & cord,
            unsigned const & strand) const
{
    return (cord >> strand) & 1ULL;
}

void Cord::setMaxLen(String<uint64_t> & cord, uint64_t const & len, uint64_t const & mask)
{
    if (len > (cord[0] & mask))
        cord[0] = len + ((cord[0]) & (~mask));
}

uint64_t Cord::getMaxLen(String<uint64_t> const & cord, uint64_t const & mask)
{
    if (empty(cord))
        return 0;
    return cord[0] & mask;
}
 
uint64_t Cord::shift(uint64_t const & val, int64_t x, int64_t y, unsigned const & bit) //add x and y
{
    if (x < 0)
        return val - ((-x) << bit) + y;
    else
        return val + (x << bit) + y;
}

bool Cord::isCordsOverlap(uint64_t & val1, uint64_t & val2, int64_t thd)
{
    int64_t dx = _DefaultCord.getCordX(val2 - val1);
    int64_t dy = get_cord_y(val2 - val1);
    return (dx >= 0) && (dx < thd) && (dy >= 0) && (dy < thd);
}

bool Cord::isBlockEnd(uint64_t & val, uint64_t const & flag)
{
    return val & flag;
}
uint64_t Cord::makeBlockEndVal(uint64_t val, uint64_t const & flag)
{
    return val | flag;
}
//todo::re wrapper !!!
uint64_t get_cord_x (uint64_t val) 
{
    return (val >> 20) & ((1ULL << 30) - 1);
}
uint64_t get_cord_y (uint64_t val) {return _DefaultCord.getCordY(val);}
uint64_t get_cord_strand (uint64_t val) {return _DefaultCord.getCordStrand(val);}
//TODO::re wrapper 
uint64_t get_cord_id (uint64_t val) 
{
    return (val >> 50 ) & ((1ULL << 10 ) - 1);
}
uint64_t shift_cord(uint64_t const & val, int64_t x, int64_t y)
{
    return _DefaultCord.shift(val, x, y);
}
void set_cord_end (uint64_t & val) {_DefaultCord.setCordEnd(val);} 
void unset_cord_end (uint64_t & val) {_DefaultHit.unsetBlockEnd(val);}
void set_cord_id (uint64_t & val, uint64_t id)
{
    val -= ((get_cord_id(val) - id) << 50);
}
uint64_t create_id_x(uint64_t const id, uint64_t const x)
{
    return (id << _DefaultCordBase.bit_id) + x;
}
void print_cord(uint64_t cord, CharString header)
{
    std::cout << header << " " 
              << get_cord_strand(cord) << " " 
              << get_cord_x(cord) << " "
              << get_cord_y(cord) << "\n";
}

uint64_t create_cord (uint64_t id, uint64_t cordx, uint64_t cordy, uint64_t strand)
{
    return _DefaultCord.createCord(create_id_x (id, cordx), cordy, strand);
}
int64_t atomic_inc_cord_y (int64_t & cord) 
{
    return atomicInc(cord);
}
void set_cord_block_end(uint64_t & val)
{
    val |= _DefaultCordBase.flagEnd;
}

void cmpRevCord(uint64_t val1, 
                uint64_t val2,
                uint64_t & cr_val1,
                uint64_t & cr_val2,
                uint64_t read_len)
{
    cr_val1 = (val1 - get_cord_y(val1) + read_len - get_cord_y(val2) - 1) ^ _DefaultCordBase.flag_strand;
    cr_val2 = (val2 - get_cord_y(val1) + read_len - get_cord_y(val2) - 1) ^ _DefaultCordBase.flag_strand;
}
uint64_t new_xy_cord (uint64_t val, uint64_t x, uint64_t y)
{
    return (val & (~_DefaultCordBase.valueMask)) + (x << _DefaultCordBase.bit) + y;
}
void set_cord_xy (uint64_t & val, uint64_t x, uint64_t y)
{
    val &= val & (~_DefaultCordBase.valueMask);
    val += (x << _DefaultCordBase.bit) + y;
}
uint64_t get_cord_xy (uint64_t val)
{
    return (val & _DefaultCordBase.valueMask);
}
void set_cord_main (uint64_t & cord)
{
    cord |= _DefaultCordBase.f_main;
}
void set_cord_gap (uint64_t & cord)
{
    cord &= ~_DefaultCordBase.f_main;
}

uint64_t is_cord_main(uint64_t cord)
{
    return cord & _DefaultCordBase.f_main;
}

void set_cord_recd(uint64_t & cord, uint64_t sgn)
{
    if (sgn != 0)
    {
        cord |= _DefaultCordBase.f_recd;
    }
    else
    {
        cord &= ~_DefaultCordBase.f_recd;
    }
}
uint64_t get_cord_recd (uint64_t cord)
{
    return cord & _DefaultCordBase.f_recd;
}

uint64_t is_cord_block_str (String<uint64_t> & cords, uint it)
{
    return it == 0 || is_cord_block_end(cords[it - 1]);
}

uint64_t is_cord_block_end (uint64_t cord)
{
    return _DefaultCord.isBlockEnd(cord);
}

void print_cords(String<uint64_t> & cords, CharString header)
{
    std::cout << header << "_cords_header \n";
    int64_t dx, dy, prex = 0, prey = 0;
    for (uint i = 0; i < length(cords); i++)
    {
        std::cout << header << " " << i << " "
                  << get_cord_id(cords[i]) << " "
                  << get_cord_strand(cords[i]) << " "
                  << get_cord_x(cords[i])  << " "
                  << get_cord_y (cords[i]) << " " 
                  << get_cord_x(cords[i]) - prex << " "
                  << get_cord_y(cords[i]) - prey << " "
                  << length(cords) << "\n";
        prex = get_cord_x(cords[i]);
        prey = get_cord_y(cords[i]);
        if (is_cord_block_end(cords[i]))
        {
            prex = 0; 
            prey = 0;
            std::cout << header << " end\n\n";
        }
    }
}
//Gap condition
int isCordsConsecutive_(uint64_t & cord1, uint64_t cord2, uint64_t thd_cord_gap)
{
    uint64_t x1 = get_cord_x (cord1);
    uint64_t x2 = get_cord_x (cord2);
    uint64_t y1 = get_cord_y (cord1);
    uint64_t y2 = get_cord_y (cord2);

    int f = !get_cord_strand(cord1 ^ cord2) &&
           x1 <= x2 && y1 <= y2 && 
           x2 - x1 < thd_cord_gap && y2 - y1 < thd_cord_gap;
    return f;
}

uint64_t make_anchor(uint64_t id, uint64_t x, uint64_t y, uint64_t strand)
{
    return create_cord (id, x - y + const_anchor_zero, y, strand);
}
/*----------  For simplicity hits and cords are not fully wrapped, use them with caution  ----------*/

int initCords (String<uint64_t> & cords)
{
    clear(cords);
    appendValue(cords, 0);
    _DefaultHit.setBlockEnd(cords[0]);
    return length(cords);
}
int initHits (String<uint64_t> & hits)
{
    return initCords(hits);
}
int initHitsScore (String<int> & hit_score)
{
    clear (hit_score);
    appendValue(hit_score, 0);
    return length(hit_score);
}
int isHitsEmpty(String<uint64_t> & hits)
{
    return length(hits) < 2;
}
int isFirstHit(Iterator<String<uint64_t> >::Type & it)
{
    return _DefaultHit.isBlockEnd(*(it - 1));
}
int isLastHit(Iterator<String<uint64_t> >::Type & it)
{
    return _DefaultHit.isBlockEnd(*it);
}
Iterator<String<uint64_t> >::Type beginHits (String<uint64_t> & hits)
{
    return begin(hits) + 1;
}
Iterator<String<uint64_t> >::Type endHits(String<uint64_t> & hits)
{
    return end(hits);
}

/*----------  Hits shortcuts  ----------*/

void Hit::setBlockStart(uint64_t & val, uint64_t const & flag)
{
    val |= flag;
}
void Hit::setBlockBody(uint64_t & val, uint64_t const & flag)
{
    val &= (~flag);
}
bool Hit::isBlockStart(uint64_t & val, uint64_t const & flag)
{
    return val & flag;
}
void Hit::setBlockEnd(uint64_t & val, uint64_t const & flag)
{
    val |= flag;
}
void Hit::unsetBlockEnd(uint64_t & val, uint64_t const & flag)
{
    val &= ~flag;
}
void Hit::setBlockStrand(uint64_t & val, uint64_t const & strand, uint64_t const & flag)
{
    if (strand)
        val |= flag;
    else
        val &= ~flag;
}
bool Hit::isBlockEnd(uint64_t & val, uint64_t const & flag)
{
    return val & flag;
}
unsigned Hit::getStrand(uint64_t const & val, uint64_t const & flag)
{
    return (val & flag)?1:0;
}
uint64_t Hit::getAnchor(uint64_t const & val)
{
    return (val >> 20 & ((1ULL << 42) - 1) & (~(1ULL << 41)));
}
void Hit::setLongPattern(uint64_t & val)
{
    val |= (1ULL << 62);
}
void Hit::unsetLongPattern(uint64_t & val)
{
    val &= ~(1ULL << 62);
}
uint64_t Hit::isLongPattern(uint64_t & val)
{
    return val & (1ULL << 62);
}
void _printHit(unsigned j, unsigned id1, unsigned id2, String<uint64_t> & hit, unsigned len)
{
    unsigned end;
    for (unsigned k = 0; k < length(hit); k++)
    {
        if (_DefaultHit.isBlockEnd(hit[k]))
            end = 1;
        else
            end = 0;
        printf("[printhit] %d %d %d %d %d\n", j, id1, id2, len, end);
    }
}

void _printHit(String<uint64_t>  & hit, CharString header)
{
    for (unsigned k = 0; k < length(hit); k++)
    {
        std::cout << "[P]::_printHit() " 
              << get_cord_id(_DefaultCord.hit2Cord(hit[k])) << " " 
              << get_cord_x(_DefaultCord.hit2Cord(hit[k])) << " " 
              << get_cord_y(hit[k]) << "\n";
        if (_DefaultHit.isBlockEnd(hit[k]))
        {
            std::cout << "[P]::_printHit() end\n";

        }
    }
}

//return if range [x11, x12) and [x21, x22) has overlap
bool _isRangeOverLap(uint64_t x11, uint64_t x12, uint64_t x21, uint64_t x22)
{
    return std::max(x11, x21) < std::min(x21, x22);
}

bool _isCordyOverLap(uint64_t cord11, uint64_t cord12, uint64_t cord21, uint64_t cord22, uint64_t read_len)
{
    return get_cord_strand(cord11 ^ cord21) ? 
                _isRangeOverLap(get_cord_y(cord11), get_cord_y(cord12),
                                read_len - 1 - get_cord_y(cord21), read_len - 1 - get_cord_y(cord22)) :
                _isRangeOverLap(get_cord_y(cord11), get_cord_y(cord12), get_cord_y(cord21), get_cord_y(cord22));
}

uint64_t getAnchorX(uint64_t anchor)
{
    return get_cord_x(_DefaultCord.hit2Cord_dstr(anchor));
}

//shortcut to extract y value of cords and convert them to y value(stry, endy) on the forward strand (y projection)
UPair getUPForwardy(UPair str_end, uint64_t read_len)
{
    if (get_cord_strand(str_end.first))
    {
        return UPair(read_len - get_cord_y(str_end.second) - 1,
                     read_len - get_cord_y(str_end.first) - 1);
    }
    else
    {
        return UPair(get_cord_y(str_end.first),
                     get_cord_y(str_end.second));
    }
}