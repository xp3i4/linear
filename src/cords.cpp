#include <seqan/stream.h>
#include <seqan/parallel.h>
#include "cords.h"
 
using namespace seqan;

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
    return (hit + ((hit & mask) << bit)) & mask2;
}

uint64_t Cord::hit2Cord_dstr(uint64_t const & hit, 
               unsigned const & bit, 
               uint64_t const & mask,
               uint64_t const & mask2
              ) const
{
    return (hit + ((hit & mask) << bit)) & mask2;
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
            typename CordBase::Flag const & end)
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
    int64_t tmp = atomicInc(cord);
    return tmp;
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

uint64_t is_cord_block_end (uint64_t cord)
{
    return _DefaultCord.isBlockEnd(cord);
}

void print_cords(String<uint64_t> & cords, CharString header)
{
    std::cout << header << "_cords_header \n";
    for (int i = 1; i < length(cords); i++)
    {
        std::cout << header << " " 
                  << get_cord_y (cords[i]) << " "  
                  << get_cord_x(cords[i]) << " "
                  << get_cord_strand(cords[i])
                  << is_cord_block_end(cords[i]) << "\n";
        if (is_cord_block_end(cords[i]))
        {
            std::cout << header << " end\n\n";
        }
    }
}
//Used as gap condition
int isCordsConsecutive_(uint64_t & cord1, uint64_t cord2, uint64_t thd_cord_gap)
{
    uint64_t x1 = get_cord_x (cord1);
    uint64_t x2 = get_cord_x (cord2);
    uint64_t y1 = get_cord_y (cord1);
    uint64_t y2 = get_cord_y (cord2);

    int f = !get_cord_strand(cord1 ^ cord2) &&
           x1 < x2 && y1 < y2 &&
           x2 - x1 < thd_cord_gap &&
           y2 - y1 < thd_cord_gap;
    if (!f)
    {
      std::cout << "iscc " << x1 << " " << x2 << " " << y1 << " " << y2 << " "<< get_cord_strand(cord1 ^ cord2) << "\n";
    }
    return f;
}

