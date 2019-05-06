#include <iostream>
#include <seqan/sequence.h>
#include <seqan/stream.h>
#include "base.h"
#include "pmpfinder.h"
#include "chain_map.h"

using namespace seqan;
const float band_width = 0.25;
const unsigned cmask = ((uint64_t)1<<20) - 1;
const unsigned cell_size = 16;
const unsigned cell_num = 12;
const unsigned window_size = cell_size * cell_num; //16*12
const unsigned window_delta = window_size * (1 - 2 * band_width);
const unsigned sup = cell_num;
const unsigned med =ceil((1 - band_width) * cell_num);
const unsigned inf = ceil((1 - 2 * band_width) * cell_num);

const unsigned initx = 5; 
const unsigned inity = 5;

const unsigned scpt_step=16;
const unsigned scpt_bit=4;
const unsigned scpt_len=5; 
const unsigned scpt_len2 = scpt_len << 1;
const unsigned scpt_len3 = scpt_len2 + scpt_len;
const int scriptMask = (1 << scpt_len) - 1;
const int scriptMask2 = scriptMask << scpt_len;
const int scriptMask3 = scriptMask2 << scpt_len;
const uint64_t hmask = (1ULL << 20) - 1;
const unsigned windowThreshold = 72; // 36;

int const typeFeatures1_32 = 1;
int const typeFeatures2_48 = 2;
int FeaturesDynamic::isFs1_32()
{
    return fs_type == typeFeatures1_32;
}
int FeaturesDynamic::isFs2_48()
{
    return fs_type == typeFeatures2_48;
}
void FeaturesDynamic::setFs1_32()
{
    fs_type = typeFeatures1_32;
}
void FeaturesDynamic::setFs2_48()
{
    fs_type = typeFeatures2_48;
}
FeaturesDynamic::FeaturesDynamic(int type)
{
    if (type == typeFeatures1_32)
    {
        fs_type = typeFeatures1_32;
    }
    else
    {
        fs_type = typeFeatures2_48;
    }
}

void decInt96 (int96 & val, int96 & dval) // val -= dval
{
    val[0] -= dval[0];
    val[1] -= dval[1];
    val[2] -= dval[2];
}
void setInt96 (int96 & val, int96 & dval) // val = dval
{
    val[0] = dval[0];
    val[1] = dval[1];
    val[2] = dval[2];
}
void incInt96 (int96 & val, int96 & dval) // val += dval
{
    val[0] += dval[0];
    val[1] += dval[1];
    val[2] += dval[2];
}
void printInt96(int96 val, CharString header)
{
    std::cout << header << " ";
    int sum = 0;
    for (int i = 0; i < 3; i++)
    {
        int v = val[i];
        for (int ii = 0; ii < 5; ii++)
        {
            sum += v & 63;
            std::cout << (v & 63) << " ";
            v >>= 6;
        }
        std::cout << "| ";
    }
    std::cout << sum << "\n";
}

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
        bit_id (40)
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

typename CordBase::Flag Cord::isCordEnd(uint64_t const & cord,
                typename CordBase::Flag const & end) const
{
    return cord & end;
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

uint64_t get_cord_x (uint64_t val) {return _getSA_i2(_DefaultCord.getCordX(val));}
uint64_t get_cord_y (uint64_t val) {return _DefaultCord.getCordY(val);}
uint64_t get_cord_strand (uint64_t val) {return _DefaultCord.getCordStrand(val);}
uint64_t get_cord_id (uint64_t val) {return _getSA_i1(_DefaultCord.getCordX(val));}
uint64_t shift_cord(uint64_t const & val, int64_t x, int64_t y)
{
    return _DefaultCord.shift(val, x, y);
}
void set_cord_end (uint64_t & val) {_DefaultCord.setCordEnd(val);}
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

void cmpRevCord(uint64_t val1, 
                uint64_t val2,
                uint64_t & cr_val1,
                uint64_t & cr_val2,
                uint64_t read_len)
{
    cr_val1 = (val1 - get_cord_y(val1) + read_len - get_cord_y(val2) - 1) ^ _DefaultCordBase.flag_strand;
    cr_val2 = (val2 - get_cord_y(val1) + read_len - get_cord_y(val2) - 1) ^ _DefaultCordBase.flag_strand;
}
uint64_t set_cord_xy (uint64_t val, uint64_t x, uint64_t y)
{
    return (val & (~_DefaultCordBase.valueMask)) + (x << _DefaultCordBase.bit) + y;
}

//Map "N" to a larger number > 64, such that (1LL << ordN_(N)) = 0.
inline int64_t ordN_(Dna5 & a)
{
    int64_t res = ordValue(a);
    return (res == 4) ? 100LL : res;
}

//======HIndex getIndexMatch()
//<<debug util
int printScript(int64_t & val, CharString header)
{
    int sum = 0;
    std::cout << header << " "; 
    for (int i = 0; i < 16; i++)
    {

        int v = (val >> (i << 2) & 15);
        std::cout << (val >> (i << 2) & 15) << " ";
        sum += v;
    }
    std::cout << " sum=" << sum << "\n";
    return sum;
}
/*----------  Script encoding type I  ----------*/
const int scptCount[5] = {1, 1<<scpt_len, 1 <<(scpt_len * 2), 0, 0};
int _scriptDist1_32(int const & s1, int const & s2)
{
    int res = std::abs((s1 & scriptMask) - (s2 & scriptMask)) +
              std::abs(((s1 >> scpt_len) & scriptMask) -
                       ((s2 >> scpt_len) & scriptMask)) + 
              std::abs((s1>> scpt_len2) - (s2 >> scpt_len2));
    return res;
}
unsigned _windowDist1_32(Iterator<String<short> >::Type const & it1, 
                         Iterator<String<short> >::Type const & it2)
{
    return _scriptDist1_32(*it1, *it2) 
         + _scriptDist1_32(*(it1 + 2), *(it2 + 2)) 
         + _scriptDist1_32(*(it1 + 4), *(it2 + 4)) 
         + _scriptDist1_32(*(it1 + 6), *(it2 + 6)) 
         + _scriptDist1_32(*(it1 + 8), *(it2 + 8)) 
         + _scriptDist1_32(*(it1 + 10), *(it2 + 10));
}
void createFeatures1_32(TIter5 const & itBegin, TIter5 const & itEnd, String<short> & f)
{
    unsigned next = 1;
    unsigned window = 1 << scpt_len;
    resize (f, ((itEnd - itBegin -window) >> scpt_bit) + 1);
    f[0] = 0;
    for (unsigned k = 0; k < window; k++)
    {
        f[0] += scptCount[ordValue(*(itBegin + k))];
    }
    for (unsigned k = scpt_step; k < itEnd - itBegin - window ; k+=scpt_step) 
    {
        f[next] = f[next - 1];
        for (unsigned j = k - scpt_step; j < k; j++)
        {
            f[next] += scptCount[ordValue(*(itBegin + j + window))] - scptCount[ordValue(*(itBegin + j))];
        }
        next++;
    }
}
uint64_t parallelParm_Static(uint64_t range, 
                             unsigned threads, 
                             unsigned & thd_id, 
                             uint64_t & thd_begin, 
                             uint64_t & thd_end)
{
    uint64_t ChunkSize = range / threads;
    unsigned id = range - ChunkSize * threads;
    if (thd_id < id)
    {
        thd_begin = ++ChunkSize * thd_id; 
    }
    else
    {
       thd_begin = id + thd_id * ChunkSize;
    }
    thd_end = thd_begin + ChunkSize; 
    return ChunkSize;
}
void createFeatures1_32(TIter5 const & itBegin, TIter5 const & itEnd, String<short> & f, unsigned threads)
{
    unsigned window = 1 << scpt_len;
    resize (f, ((itEnd - itBegin -window) >> scpt_bit) + 1);
#pragma omp parallel
{
    unsigned thd_id = omp_get_thread_num();
    uint64_t thd_begin;
    uint64_t thd_end;
    uint64_t range = (itEnd - itBegin - window - scpt_step) / scpt_step;
    parallelParm_Static(range, threads, 
                        thd_id,  thd_begin, thd_end);
    uint64_t next = thd_begin;
    thd_begin *= scpt_step;
    thd_end *= scpt_step;
    f[next] = 0;
    for (unsigned k = thd_begin; k < thd_begin + window; k++)
    {
        f[next] += scptCount[ordValue(*(itBegin + k))];
    }
    next++;
    for (unsigned k = thd_begin + scpt_step; k < thd_end ; k += scpt_step) 
    {
        f[next] = f[next - 1];
        for (unsigned j = k - scpt_step; j < k; j++)
        {
            f[next] += scptCount[ordValue(*(itBegin + j + window))] - scptCount[ordValue(*(itBegin + j))];
        }
        next++;
    }
}
}

/*----------  Script encoding type II  ----------*/
/* TG TC TA GT GG int1  \
   GC GA CT CG CC int2  -- int96 for 48 bases script; Each 2mer has 6 bits
   CA AT AG AC AA int3  /
 */
//TODO!!! for 2mer > 32 occurrences out of boundary

int const window48 = 48;
int const max31 = 31;
int const mxu31 = (max31 << 24) + (max31 << 18) + (max31 << 12) +
                  (max31 << 6) + max31;
int64_t __scriptDist63_31(int const & s1, int const & s2)
{
    int d = s1 + mxu31 - s2;
    int mask = 63; 
    return std::abs((d >> 24 & mask) - max31) +
           std::abs((d >> 18 & mask) - max31) +
           std::abs((d >> 12 & mask) - max31) +
           std::abs((d >> 6 & mask) - max31) +
           std::abs((d & mask) - max31);
}
int64_t _scriptDist63_31(int const & s11, int const & s12, int const & s13,
                         int const & s21, int const & s22, int const & s23)
{
    return __scriptDist63_31(s11, s21) +
           __scriptDist63_31(s12, s22) +
           __scriptDist63_31(s13, s23);
}
//wrapper into 96 bit integer.
int64_t _scriptDist63_31(int96 & s1, int96 & s2)
{
    return  _scriptDist63_31(s1[0], s1[1], s1[2],
                             s2[0], s2[1], s2[2]);
}
int64_t _windowDist48_4(Iterator<String<int96> >::Type it1,  //feature string iterator
                        Iterator<String<int96> >::Type it2)
{
    /*
    printInt96(*it1, "4841");
    printInt96(*(it1 + 1), "4841");
    printInt96(*(it1 + 2), "4841");
    printInt96(*(it1 + 3), "4841");
    printInt96(*it2, "4842");
    printInt96(*(it2 + 1), "4842");
    printInt96(*(it2 + 2), "4842");
    printInt96(*(it2 + 3), "4842");
    */
    return _scriptDist63_31(*it1, *it2) + 
           _scriptDist63_31(*(it1 + 1), *(it2 + 1)) +
           _scriptDist63_31(*(it1 + 2), *(it2 + 2)) + 
           _scriptDist63_31(*(it1 + 3), *(it2 + 3));
    
}
/**
 * 1.units[n] = (i << 8 + k) maps n to bits of int96 (ith int, kth bit);
 * 2.NOTE::In gcc, x << y == x << (y|31) due to optimization concerning.
 *   In other words if y > 31, it does nothing. 
 *   Beacause of that the infiN is set to 31 for TT, N* and *N 
 *   They are add to the 31bits and cut off the 30bits by & infi_mask30
 */
short infiN = 31; //for base '*N' 'N*' and 'TT' 
unsigned infi_mask30 = (1 << 31) - 1;
short units[25] = {
/*A**/ 0,             6,             12,            18,            infiN, 
/*C**/ 24,            (1 << 8) + 0,  (1 << 8) + 6,  (1 << 8) + 12, infiN,
/*G**/ (1 << 8) + 18, (1 << 8) + 24, (2 << 8) + 0,  (2 << 8) + 6,  infiN,
/*T**/ (2 << 8) + 12, (2 << 8) + 18, (2 << 8) +24,  infiN,         infiN,  
/*N**/ infiN,         infiN,         infiN,         infiN,         infiN};
void add2merInt96 (int96 & val, TIter5 it)
{
    unsigned ordV = ordValue(*it) * 5 + ordValue(*(it + 1)); // 5*ordV_i + ordV_{i+1}
    unsigned i = units[ordV] >> 8;
    unsigned addVal = (1 << (units[ordV] & 255)) & infi_mask30;
    val[i] += addVal;
    //std::cout << "a2i96 " << *it << *(it+1) << " " << ordV << " " << i << " " << addVal << " " << (1 << 255)<< "\n";
}
int createFeatures48(TIter5 it_str, TIter5 it_end, String<int96> & f)
{
    double t = sysTime();
    int addMod3[3] = {1, 2, 0};   // adMod3[i] = ++ i % 3;
    int96 zero96 = {0, 0, 0};
    std::vector<int96> buffer(3, zero96); //buffer of 3 cells in one script
    resize (f, (it_end - it_str - window48) / scpt_step + 1); 
    setInt96(f[0], zero96);
    //std::cout << "cf48 " << length(f) << "\n";
    for (int i = 0; i < 3; i++) //init f[0]
    {
        for (int j = i << scpt_bit; j < (i << scpt_bit) + scpt_step; j++)
        {
            add2merInt96(buffer[i], it_str + j);
        }
        //printInt96(buffer[i], "cf3");
        incInt96(f[0], buffer[i]);
    }
    //std::cerr << " cf48 time " << sysTime() - t << "\n";
    int next = 1; //stream f[next]
    int ii = 0;
    for (int i = scpt_step; i < it_end - it_str - window48 - 1; i += scpt_step) 
    {
        setInt96(f[next], f[next - 1]);
        decInt96(f[next], buffer[ii]);
        setInt96(buffer[ii], zero96); //clear to 0;
        for (int j = i - scpt_step + window48; j < i + window48; j++)
        {
            add2merInt96(buffer[ii], it_str + j);
            //std::cout << *(it_str + j);
        }
        //printInt96(buffer[ii], "cf2");
        incInt96(f[next], buffer[ii]);
        ii = addMod3[ii];
        next++;
    }
    return 0;
}
int createFeatures48(TIter5 it_str, TIter5 it_end, String<int96> & f, unsigned threads)
{
    int window = window48;
    if (it_end - it_str < window)
    {
        return 0;
    }
    resize (f, ((it_end - it_str -window) >> scpt_bit) + 1);
    int64_t range = (it_end - it_str - window) / scpt_step + 1; //numer of windows 
    if (range < threads)
    {
        createFeatures48(it_str, it_end, f);
        return 0;
    }
#pragma omp parallel
{
    unsigned thd_id = omp_get_thread_num();
    int64_t chunk_size = range / threads; 
    int64_t thd_begin = thd_id * (chunk_size + 1);
    unsigned id1 = range - chunk_size * threads;
    if (thd_id >= id1)
    {
        thd_begin = id1 + chunk_size * thd_id;
    }
    else
    {
        ++chunk_size;
    }
    int64_t thd_end = thd_begin + chunk_size;
    int64_t next = thd_begin;
    thd_begin *= scpt_step;
    thd_end *= scpt_step;
    //std::cout << "cfs2 " << thd_id << " " << thd_begin << " " << thd_end << "\n";

    int addMod3[3] = {1, 2, 0};   // adMod3[i] = ++ i % 3;
    int96 zero96 = {0, 0, 0};
    std::vector<int96> buffer(3, zero96); //buffer of 3 cells in one script
    setInt96(f[next], zero96);
    //std::cout << "cf48 " << length(f) << "\n";
    for (int i = 0; i < 3; i++) //init buffer and f[next]
    {
        int tmp = thd_begin + (i << scpt_bit); 
        for (int j = tmp; j < tmp + scpt_step; j++)
        {
            add2merInt96(buffer[i], it_str + j);
        }
        //printInt96(buffer[i], "cfp3");
        //std::cout << "cfp3 " << tmp << '\n';
        incInt96(f[next], buffer[i]);
    }
    int ii = 0;
    next++; //stream f[next]
    for (int i = thd_begin + scpt_step; i < thd_end; i += scpt_step) 
    {
        setInt96(f[next], f[next - 1]);
        decInt96(f[next], buffer[ii]);
        setInt96(buffer[ii], zero96); //clear to 0;
        for (int j = i - scpt_step + window; j < i + window; j++)
        {
            add2merInt96(buffer[ii], it_str + j);
            //std::cout << *(it_str + j);
        }
        //printInt96(buffer[ii], "cf2");
        incInt96(f[next], buffer[ii]);
        ii = addMod3[ii];
        next++;
    }
}
}

/*----------  Script encoding wrapper  ----------*/
unsigned _windowDist(Iterator<String<short> >::Type const & it1, 
                     Iterator<String<short> >::Type const & it2)
{
    return _windowDist1_32(it1, it2);
}
unsigned _windowDist(Iterator<String<int96> >::Type const & it1, 
                     Iterator<String<int96> >::Type const & it2)
{
    return _windowDist48_4(it1, it2);
}
int createFeatures(TIter5 it_str, TIter5 it_end, FeaturesDynamic & f)
{
    if (f.isFs1_32())
    {
        createFeatures1_32(it_str, it_end, f.fs1_32);
    }
    else if (f.isFs2_48())
    {
        createFeatures48(it_str, it_end, f.fs2_48);
    }
}
int createFeatures(TIter5 it_str, TIter5 it_end, FeaturesDynamic & f, unsigned threads)
{
    if (f.isFs1_32())
    {
        createFeatures1_32(it_str, it_end, f.fs1_32, threads);
    }
    else if (f.isFs2_48())
    {
        createFeatures48(it_str, it_end, f.fs2_48, threads);
    }
}
int createFeatures(StringSet<String<Dna5> > & seq, 
                   StringSet<FeaturesDynamic> & f)
{
    resize(f, length(seq));
    for (unsigned k = 0; k < length(seq); k++)
    {
        createFeatures(begin(seq[k]), end(seq[k]), f[k]);
    }
}
int createFeatures(StringSet<String<Dna5> > & seq, 
                   StringSet<FeaturesDynamic> & f, 
                   unsigned threads)
{
    resize(f, length(seq));
    for (unsigned k = 0; k < length(seq); k++)
    {
        createFeatures(begin(seq[k]), end(seq[k]), f[k], threads);
    }
}

/*----------  Dynamic programmign of extending path (tiles)  ----------*/
bool initCord(String<uint64_t> & hit, unsigned & currentIt, String<uint64_t> & cord)
{
    currentIt = 0;
    if (empty(hit))
        return false;
    else
        appendValue(cord, _DefaultCord.hit2Cord(hit[0]));
    return true;
}
/*----------  previouse & next window(tile)  ----------*/
//Template function is not used to speed up the compilation 
//So there code for short are copied to int96.
uint64_t previousWindow(String<short> & f1, 
                        String<short> & f2, 
                        uint64_t cord,
                        float & score,
                        unsigned window_threshold = windowThreshold)
{
    uint64_t genomeId = _getSA_i1(_DefaultCord.getCordX(cord));
    uint64_t strand = get_cord_strand(cord);
    uint64_t x_suf = _DefaultCord.cord2Cell(get_cord_x(cord));
    uint64_t y_suf = _DefaultCord.cord2Cell(get_cord_y(cord));
    uint64_t x_min = 0;
    uint64_t y;
    uint64_t new_cord = 0;
    
    if (y_suf < med || x_suf < sup)
        return 0;
    else 
        y = y_suf - med;

    unsigned min = ~0;
    for (uint64_t x = x_suf - sup; x < x_suf - inf; x += 1) 
    {
        unsigned tmp = _windowDist(begin(f1) + y, begin(f2) + x);
        if (tmp < min)
        {
            min = tmp;
            x_min = x;
        }
    }
    if (min > window_threshold)
        return 0;    
    else 
    {
        if ( x_suf - x_min > med)
        {
            new_cord = _DefaultCord.createCord(create_id_x(genomeId, _DefaultCord.cell2Cord(x_suf - med)),  _DefaultCord.cell2Cord(x_suf - x_min - med + y), strand);
        }
        else
        {
            new_cord = _DefaultCord.createCord(create_id_x(genomeId, _DefaultCord.cell2Cord(x_min)), _DefaultCord.cell2Cord(y), strand);
        } 
    }
    score += min;
    return new_cord;
}
uint64_t previousWindow(String<int96> & f1, 
                        String<int96> & f2, 
                        uint64_t cord,
                        float & score,
                        unsigned window_threshold = windowThreshold)
{
    uint64_t genomeId = _getSA_i1(_DefaultCord.getCordX(cord));
    uint64_t strand = get_cord_strand(cord);
    uint64_t x_suf = _DefaultCord.cord2Cell(get_cord_x(cord));
    uint64_t y_suf = _DefaultCord.cord2Cell(get_cord_y(cord));
    uint64_t x_min = 0;
    uint64_t y;
    uint64_t new_cord = 0;
    
    if (y_suf < med || x_suf < sup)
        return 0;
    else 
        y = y_suf - med;

    unsigned min = ~0;
    for (uint64_t x = x_suf - sup; x < x_suf - inf; x += 1) 
    {
        unsigned tmp = _windowDist(begin(f1) + y, begin(f2) + x);
        //<<debug
        std::cout << "pw1 " << tmp << " " << x * 16 << " " << get_cord_y(cord) << "\n";
        //>>debug
        if (tmp < min)
        {
            min = tmp;
            x_min = x;
        }
    }
    if (min > window_threshold)
        return 0;    
    else 
    {
        if ( x_suf - x_min > med)
        {
            new_cord = _DefaultCord.createCord(create_id_x(genomeId, _DefaultCord.cell2Cord(x_suf - med)),  _DefaultCord.cell2Cord(x_suf - x_min - med + y), strand);
        }
        else
        {
            new_cord = _DefaultCord.createCord(create_id_x(genomeId, _DefaultCord.cell2Cord(x_min)), _DefaultCord.cell2Cord(y), strand);
        } 
    }
    score += min;
    return new_cord;
}
uint64_t nextWindow(String<short> & f1, 
                    String<short> & f2, 
                    uint64_t cord,
                    float & score,
                    unsigned window_threshold = windowThreshold
                    )
{
    uint64_t genomeId = get_cord_id(cord);
    uint64_t strand = get_cord_strand(cord);
    uint64_t x_pre = _DefaultCord.cord2Cell(get_cord_x(cord));
    uint64_t y_pre = _DefaultCord.cord2Cell(get_cord_y(cord));
    uint64_t x_min = 0;
    uint64_t y;
    uint64_t new_cord = 0;
    unsigned min = ~0;
    
    if (y_pre + sup * 2 > length(f1) || x_pre + sup * 2> length(f2))
        return 0;
    else 
        y = y_pre + med;
    
    for (uint64_t x = x_pre + inf; x < x_pre + sup; x += 1) 
    {
        unsigned tmp = _windowDist(begin(f1) + y, begin(f2) + x);
        if (tmp < min)
        {
            min = tmp;
            x_min = x;
        }
    }
    if (min > window_threshold)
    {
       return 0;
    }
    else 
    {
        if ( x_min - x_pre > med)
        {
            new_cord = _DefaultCord.createCord(create_id_x(genomeId, _DefaultCord.cell2Cord(x_pre + med)),  _DefaultCord.cell2Cord(x_pre + med - x_min + y), strand);
        }
        else
        {
            new_cord = _DefaultCord.createCord(create_id_x(genomeId, _DefaultCord.cell2Cord(x_min)), _DefaultCord.cell2Cord(y), strand);
        }
    }
    score += min;
    return new_cord;
}
uint64_t nextWindow(String<int96> & f1, 
                    String<int96> & f2, 
                    uint64_t cord,
                    float & score,
                    unsigned window_threshold = windowThreshold
                    )
{
    uint64_t genomeId = get_cord_id(cord);
    uint64_t strand = get_cord_strand(cord);
    uint64_t x_pre = _DefaultCord.cord2Cell(get_cord_x(cord));
    uint64_t y_pre = _DefaultCord.cord2Cell(get_cord_y(cord));
    uint64_t x_min = 0;
    uint64_t y;
    uint64_t new_cord = 0;
    unsigned min = ~0;
    
    if (y_pre + sup * 2 > length(f1) || x_pre + sup * 2> length(f2))
        return 0;
    else 
        y = y_pre + med;
    
    for (uint64_t x = x_pre + inf; x < x_pre + sup; x += 1) 
    {
        unsigned tmp = _windowDist(begin(f1) + y, begin(f2) + x);
        std::cout << "nw1 " << tmp << " " << x * 16 << "\n";
        if (tmp < min)
        {
            min = tmp;
            x_min = x;
        }
    }
    if (min > window_threshold)
    {
       return 0;
    }
    else 
    {
        if ( x_min - x_pre > med)
        {
            new_cord = _DefaultCord.createCord(create_id_x(genomeId, _DefaultCord.cell2Cord(x_pre + med)),  _DefaultCord.cell2Cord(x_pre + med - x_min + y), strand);
        }
        else
        {
            new_cord = _DefaultCord.createCord(create_id_x(genomeId, _DefaultCord.cell2Cord(x_min)), _DefaultCord.cell2Cord(y), strand);
        }
    }
    score += min;
    return new_cord;
}

void checkPath(StringSet<String<Dna5> > & cords, StringSet<String<Dna5> > const & reads)
{
    unsigned count = 0;
    Iterator<StringSet<String<Dna5> > >::Type it = begin(cords);
    for (auto && read : reads)
    {
        if(empty(*it))
            count++;
        else
            if (get_cord_y(back(*it)) + window_size * 2 < length(read))
            {
                count++;
            }
        it++;
    }
}

/*=============================================
=            Mapping and anchoring            =
=============================================*/
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

void _printHit(String<uint64_t>  & hit)
{
    for (unsigned k = 0; k < length(hit); k++)
    {
        std::cout << "[P]::_printHit() " 
              << _getSA_i1(_DefaultCord.getCordX(_DefaultCord.hit2Cord(hit[k]))) << " " 
              << _getSA_i2(_DefaultCord.getCordX(_DefaultCord.hit2Cord(hit[k]))) << " " 
              << get_cord_y(hit[k]) << "\n";
        if (_DefaultHit.isBlockEnd(hit[k]))
        {
            std::cout << "[P]::_printHit() end\n";
        }
    }
}

/**
 * Search double strand pattern in the index and
 * append to anchors
 */
 unsigned getIndexMatchAll(LIndex & index,
                           String<Dna5> & read,
                           String<uint64_t> & set,
                           MapParm & mapParm)
{   
    int dt = 0;
    LShape shape(index.shape);
    uint64_t xpre = 0;
    hashInit(shape, begin(read));
    for (unsigned k = 0; k < length(read); k++)
    {
        hashNexth(shape, begin(read) + k);
        uint64_t pre = ~0;
        if (++dt == mapParm.alpha)
        {
            dt = 0;
            if(hashNextX(shape, begin(read) + k) ^ xpre)
            {
                xpre = shape.XValue;
                uint64_t pos = getXDir(index, shape.XValue, shape.YValue);
                uint64_t ptr = _DefaultHs.getHeadPtr(index.ysa[pos-1]);
                if (index.isEmptyDir(pos) || ptr >= mapParm.delta)
                {
                    dt = 0;
                    continue;
                }
                while ((_DefaultHs.getHsBodyY(index.ysa[pos]) == shape.YValue || 
                        _DefaultHs.getHsBodyY(index.ysa[pos]) == 0))
                {
                    if (_DefaultHs.getHsBodyS(pre - index.ysa[pos]) > mapParm.kmerStep)
                    {
                        if (((index.ysa[pos] & _DefaultHsBase.bodyCodeFlag) >>_DefaultHsBase.bodyCodeBit) ^ shape.strand)
                        {
                            
                            uint64_t cordy = length(read) - 1 - k;
                            appendValue(set, (((_DefaultHs.getHsBodyS(index.ysa[pos]) - cordy) << 20) 
                                + cordy) | _DefaultHitBase.flag2);
                        }
                        else
                        {    
                            uint64_t cordy = k;
                            appendValue(set, ((_DefaultHs.getHsBodyS(index.ysa[pos]) - cordy) << 20) | cordy);
                        }
                        pre = index.ysa[pos];
                    }
                    if (++pos > length(index.ysa) - 1)
                    {
                        break;
                    }
                }
            }
        }
    }
    return 0;
}

 uint64_t getAnchorMatchAll(Anchors & anchors, unsigned const & readLen, MapParm & mapParm, String<uint64_t> & hit)
{
    uint64_t ak, maxAnchor = 0;
    unsigned c_b=mapParm.shapeLen, sb=0, maxStart=0, maxEnd=0;
    anchors[0] = anchors[1];
    ak=anchors[0];
    anchors.sort(anchors.begin(), anchors.end());
    for (unsigned k = 1; k <= anchors.length(); k++)
    {
        //if (anchors[k] - ak >= AnchorBase::AnchorValue)
        //if (anchors[k] - ak > mapParm.anchorDeltaThr)
        if (anchors[k] - ak > (((unsigned)(0.1 * readLen)) << 20))//mapParm.anchorDeltaThr)
        {
            if (c_b > mapParm.anchorLenThr * readLen)
            {
                anchors.sortPos2(anchors.begin() + sb, anchors.begin() + k);
                for (unsigned m = sb; m < k; m++)
                {
                    appendValue(hit, anchors[m]);
                }
                _DefaultHit.setBlockEnd(back(hit));
            }
            if (maxAnchor < c_b)
            {
                maxAnchor = c_b;
                maxStart = sb; 
                maxEnd = k;
            }
            sb = k;
            ak = anchors[k];
            c_b = mapParm.shapeLen;
        }
        else 
        {
            if(anchors.deltaPos2(k, k - 1) >  mapParm.shapeLen)
                c_b += mapParm.shapeLen;
            else
                c_b += anchors.deltaPos2(k, k - 1); 
        }
    }
    if (empty(hit) && maxEnd ) // maxStart < maxEnd
    {
        for (unsigned k = maxStart; k < maxEnd; k++)
        {
            appendValue(hit, anchors[k]);
        }
        _DefaultHit.setBlockEnd(back(hit));

    }
    return maxAnchor;
}

 uint64_t getAnchorMatchFirst(Anchors & anchors, unsigned const & readLen, MapParm & mapParm, String<uint64_t> & hit)
{
    uint64_t ak, maxAnchor = 0;
    unsigned c_b=mapParm.shapeLen, sb=0, maxStart=0, maxEnd=0;
    anchors[0] = anchors[1];
    ak=anchors[0];
    anchors.sort(anchors.begin(), anchors.end());
    for (unsigned k = 1; k <= anchors.length(); k++)
    {
        if (anchors[k] - ak > (((unsigned)(0.1 * readLen)) << 20))//mapParm.anchorDeltaThr)
        {
            if (maxAnchor < c_b)
            {
                maxAnchor = c_b;
                maxStart = sb; 
                maxEnd = k;
            }
            sb = k;
            ak = anchors[k];
            c_b = mapParm.shapeLen;
        }
        else 
        {
            if(anchors.deltaPos2(k, k - 1) >  mapParm.shapeLen)
                c_b += mapParm.shapeLen;
            else
                c_b += anchors.deltaPos2(k, k - 1); 
        }
    }
    if (maxEnd) // maxStart < maxEnd
    {
        anchors.sortPos2(anchors.begin() + maxStart, anchors.begin() + maxEnd);
        for (unsigned k = maxStart; k < maxEnd; k++)
        {
            appendValue(hit, anchors[k]);
        }
        _DefaultHit.setBlockEnd(back(hit));
    }
    return maxAnchor;
}

 uint64_t getAnchorMatchList(Anchors & anchors, unsigned const & readLen, MapParm & mapParm, String<uint64_t> & hit)
{
    uint64_t ak;
    uint64_t c_b=mapParm.shapeLen, sb=0, sc = 0;
    anchors[0] = anchors[1];
    ak=anchors[0];
    anchors.sort(anchors.begin(), anchors.end());
    uint64_t mask = (1ULL << 20) - 1;
    int64_t list[200000]={0};
    unsigned lcount = 0;
    //printf("[debug]::getAnchorMatchList %d\n", anchors.length());
    for (unsigned k = 1; k < anchors.length(); k++)
    {
//TODO 0.1 should be replaced by a cutoff in marpparm
        if (anchors[k] - ak > (((unsigned)(0.1 * readLen)) << 20))//mapParm.anchorDeltaThr)
        {
            if (c_b > mapParm.anchorLenThr * readLen)
            {
                anchors.sortPos2(anchors.begin() + sb, anchors.begin() + k);
                list[lcount++] = (c_b << 40) + (sb << 20) + k;
            }
            sb = k;
            ak = anchors[k];
            c_b = mapParm.shapeLen;
        }
        else 
        {
            //cb += std::min(anchors.deltaPos2(k, k - 1), mapParm.shapeLen);
            if(anchors.deltaPos2(k, k - 1) >  mapParm.shapeLen)
                c_b += mapParm.shapeLen;
            else
                c_b += anchors.deltaPos2(k, k - 1); 
        }
    }
    /*
    //<<debug
        for (int ii = 0; ii < anchors.length(); ii++)
        {
            uint64_t mask1 = (1ULL << 20) - 1;
            uint64_t mask2 = (1ULL << 40) - 1;
            uint64_t tmp_cord = _DefaultCord.hit2Cord(anchors[ii]);

            std::cout << "gaml1 " << ii << " " << get_cord_y(tmp_cord) << " " << get_cord_x(anchors[ii]) << " " << get_cord_x(tmp_cord) << "\n";
        }
    //>>debug
    */
    if (list[0])
    {
        std::sort (list, list + lcount, std::greater<uint64_t>());
        for (unsigned k = 0; k < mapParm.listN; k++)
        {
          if (((list[0] / 10) < list[k])  && list[k])
          {
              sb = ((list[k] >> 20) & mask);
              sc = list[k] & mask;
              for (unsigned n = sb; n < sc; n++)
              {
                  appendValue(hit, anchors[n]);
              }   
              _DefaultHit.setBlockEnd(back(hit));
          }
          else
              break;
        }
        std::cout << "gaml2 " << anchors.length() << " " << length(hit) << "\n";
        return (list[0] >> 40);   
    }
    else
    {
       return 0;
    }
}

 uint64_t mnMapReadAll( LIndex  & index,
                           String<Dna5> & read,
                          Anchors & anchors,
                          MapParm & mapParm,
                          String<uint64_t> & hit  
                         )
{
    getIndexMatchAll(index, read, anchors.set, mapParm);    
    return getAnchorMatchAll(anchors, length(read), mapParm, hit);
}

 uint64_t mnMapReadFirst( LIndex  & index,
                           String<Dna5> & read,
                          Anchors & anchors,
                          MapParm & mapParm,
                          String<uint64_t> & hit  
                         )
{
    getIndexMatchAll(index, read, anchors.set,  mapParm);    
    return getAnchorMatchFirst(anchors, length(read), mapParm, hit);
}

 uint64_t mnMapReadList( LIndex  & index,
                           String<Dna5> & read,
                          Anchors & anchors,
                          MapParm & mapParm,
                          String<uint64_t> & hit)
{
    getIndexMatchAll(index, read, anchors.set, mapParm);    
    //printf("done getinxmatchall\n");
    return getAnchorMatchList(anchors, length(read), mapParm, hit);
}
/*=====  End of Mapping and anchoring  ======*/

/*
 * this is initCord for double strand index(with flag in cord value)
 */
 bool initCord(typename Iterator<String<uint64_t> >::Type & it, 
                     typename Iterator<String<uint64_t> >::Type & hitEnd,
                     unsigned & preCordStart,
                     String<uint64_t> & cord
                    )
{
        
    if (empty(cord))
    {
        appendValue(cord, 0);
        _DefaultHit.setBlockEnd(cord[0]);
    }

    if (it == hitEnd)
        return false;
    else
    {
        //hit2Cord_dstr keeps the strand flag in cord value
        appendValue(cord, _DefaultCord.hit2Cord_dstr(*(it)));
        ++it;
        preCordStart = length(cord) - 1;   
    }
    return true;
}

/*
 * endCord for double strand index
 */
 bool endCord( String<uint64_t> & cord,
                     unsigned & preCordStart,
                     float const & cordThr,
                     float & score)
{
    _DefaultCord.setMaxLen(cord, length(cord) - preCordStart);   
    if (length(cord) - preCordStart > cordThr)// > std::max(score/25, cordThr))
    {
        _DefaultHit.setBlockEnd(back(cord));
    }
    else
    {
        erase(cord, preCordStart, length(cord));
    }
    (void)score;
    return true;
}

/*
 * nextCord for double strand sequence
 */
 bool nextCord(typename Iterator<String<uint64_t> >::Type & it, 
               typename Iterator<String<uint64_t> >::Type const & hitEnd, 
               unsigned & preCordStart,
               String<uint64_t> & cord,
               float const & cordThr,
               float & score
               )
{
//TODO: add maxlen of anchor to the first node in cord;
    if (it >= hitEnd)
        return false;
    
    while(!_DefaultHit.isBlockEnd(*(it - 1)))
    {
        if(get_cord_y(*(it))>get_cord_y(back(cord)) +  window_delta) 
        {
            appendValue(cord, _DefaultCord.hit2Cord_dstr(*(it)));
            ++it;
            std::cout << "nc1 " << get_cord_y(back(cord)) << "\n";
            return true;
        }
        ++it;
    }
    _DefaultHit.setBlockEnd(back(cord));
    
    if(it < hitEnd)
    {
        if (length(cord) - preCordStart < std::max(score/25,  cordThr) )//|| std::abs(b1 - b2) < 10)
        {
            erase(cord, preCordStart, length(cord));
            //std::cerr << "[]::nextCord erase cord\n";
        }
        else
        {
            _DefaultCord.setMaxLen(cord, length(cord) - preCordStart);
    //printf("[debug]::nextcord new block %f %d %f\n", (float)score/(length(cord) - preCordStart), length(cord) - preCordStart, cordThr);
        }
        preCordStart = length(cord);
        appendValue(cord, _DefaultCord.hit2Cord_dstr(*(it)));
        score = 0;
        ++it;
        return true;
    }
    else
    {
        score = 0;
        return false;
    }
}

bool extendWindow(String<short> &f1, 
                  String<short> & f2, 
                  String<uint64_t> & cords, 
                  float & score, 
                  uint64_t & strand)
{
    uint64_t pre_cord_y = (_DefaultHit.isBlockEnd(cords[length(cords) - 2]))?
    0:get_cord_y(cords[length(cords) - 2]) + window_delta;
    unsigned len = length(cords) - 1;
    uint64_t new_cord;
    while (pre_cord_y<= get_cord_y(back(cords)))
    {
        new_cord = previousWindow(f1, f2, back(cords), score, windowThreshold);
        if (new_cord && get_cord_y(new_cord) > pre_cord_y)
        {
            appendValue(cords, new_cord);
        }
        else
        {
            break;
        }
    } 
    for (unsigned k = len; k < ((length(cords) + len) >> 1); k++) 
    {
        std::swap(cords[k], cords[length(cords) - k + len - 1]);
    }
    while (true)
    {
        new_cord = nextWindow(f1, f2, back(cords), score, windowThreshold);
        if (new_cord)
        {
            appendValue(cords, new_cord);
        }
        else
        {
            break;
        }
    }
    return true;
}
bool extendWindow(String<int96> & f1, 
                  String<int96> & f2, 
                  String<uint64_t> & cords, 
                  float & score, 
                   uint64_t & strand)
{
    uint64_t pre_cord_y = (_DefaultHit.isBlockEnd(cords[length(cords) - 2]))?
    0:get_cord_y(cords[length(cords) - 2]) + window_delta;
    unsigned len = length(cords) - 1;
    uint64_t new_cord;
    std::cout << "ew1 " << pre_cord_y << " " << get_cord_y(back(cords)) << "\n";
    while (pre_cord_y < get_cord_y(back(cords)))
    {
        new_cord = previousWindow(f1, f2, back(cords), score, windowThreshold);
        if (new_cord && get_cord_y(new_cord) > pre_cord_y)
        {
            appendValue(cords, new_cord);
        }
        else
        {
            break;
        }
    } 
    for (unsigned k = len; k < ((length(cords) + len) >> 1); k++) 
    {
        std::swap(cords[k], cords[length(cords) - k + len - 1]);
    }
    while (true)
    {
        new_cord = nextWindow(f1, f2, back(cords), score, windowThreshold);
        if (new_cord)
        {
            appendValue(cords, new_cord);
        }
        else
        {
            break;
        }
    }
    return true;
}
bool extendWindow(FeaturesDynamic & f1, 
                  FeaturesDynamic & f2, 
                  String<uint64_t> & cords, 
                  float & score, 
                   uint64_t & strand)
{
    if (f1.isFs1_32())
    {
        extendWindow (f1.fs1_32, f2.fs1_32, cords, score, strand);
    }
    else if (f1.isFs2_48())
    {
        extendWindow (f1.fs2_48, f2.fs2_48, cords, score, strand);
    }
}
bool isOverlap (uint64_t cord1, uint64_t cord2, 
                int revscomp_const, 
                int overlap_size
                )
{
    uint64_t strand1 = _DefaultCord.getCordStrand(cord1);
    uint64_t strand2 = _DefaultCord.getCordStrand(cord2);
    int64_t x1 = _DefaultCord.getCordX(cord1);
    //int64_t y1 = revscomp_const * strand1 - _nStrand(strand1) * get_cord_y(cord1);
    int64_t y1 = _DefaultCord.getCordY (cord1);
    int64_t x2 = _DefaultCord.getCordX(cord2);
    //int64_t y2 = revscomp_const * strand2 - _nStrand(strand2) * get_cord_y(cord2);
    int64_t y2 = _DefaultCord.getCordY (cord2);
    (void) revscomp_const;
    return std::abs(x1 - x2) < overlap_size && 
           std::abs(y1 - y2) < overlap_size && (!(strand1 ^ strand2));
}
/**
 * cord1 is predecessor of cord2 and they are not overlapped
 */
 bool isPreGap (uint64_t cord1, uint64_t cord2, 
                int revscomp_const, 
                int gap_size
                )
{

    int64_t x1 = _DefaultCord.getCordX(cord1);
    int64_t x2 = _DefaultCord.getCordX(cord2);
    (void) revscomp_const;
    return (x1 + gap_size <= x2);
}
/**
 * cord1 is successor of cord2 and they are not overlapped
 */
 bool isSucGap (uint64_t cord1, uint64_t cord2, 
                int revscomp_const,
                int gap_size
                )
{
    return isPreGap (cord2, cord1, revscomp_const, gap_size);
}
/**
 * Extend windows between cord1, and cord2 if they are not overlapped,
 * and insert the windows to the k th element of the cords. 
 * if cord1 and cord2 have the same strand 
 * then call previousWindow for cord1 and nextWindow for cord2 until x1 + window_size < x2
 * if cord1 and cord2 have different strand 
 * then call nextWindow for cord1 and previousWindow for cord2 along each own strand until it can't be extended any more.
 * 
 */         
 int extendPatch(StringSet<String<FeatureType> > & f1, 
                 StringSet<String<FeatureType> > & f2, 
                 String<uint64_t> & cords,
                 int kk,
                 uint64_t cord1,
                 uint64_t cord2,
                 int revscomp_const,
                 int overlap_size,
                 int gap_size
                )
{
    unsigned window_threshold = 30;
    float score = 0;
    std::cout << "eP dg1_1 " << get_cord_y(cord1) << " " << get_cord_y(cord2) << "\n";
    if (isOverlap(cord1, cord2, revscomp_const, overlap_size))
    {
        return 0;
    }
    uint64_t pcord = cord1;
    uint64_t scord = cord2;
    /*
    if (isSucGap(cord1, cord2, revscomp_const))
    {
        pcord = cord2;
        scord = cord1;
    }
    */
    uint64_t strand1 = _DefaultCord.getCordStrand(pcord);
    uint64_t strand2 = _DefaultCord.getCordStrand(scord);
    uint64_t genomeId1 = get_cord_id(pcord);
    uint64_t genomeId2 = get_cord_id(scord);
    int len = 0;
    uint64_t cord = pcord;
    String<uint64_t> tmp;
    std::cout << "dg1_ " << get_cord_y(cord) << " " << get_cord_y(scord) << "\n";
    while (isPreGap(cord, scord, revscomp_const, gap_size))
    {
        cord = nextWindow (f1[strand1], f2[genomeId1], cord, score, window_threshold);
        std::cout << "dg1_ " << get_cord_y(cord) << "\n";
        if (cord)
        {
            appendValue (tmp, cord);
        }
        else
        {
            break;
        }
    }
    uint64_t nw = pcord;
    if (!empty (tmp))
    {
        len += length(tmp);
        nw = back(tmp);
        insert(cords, kk, tmp);
        clear(tmp);
    }
    
    cord = scord;
    while (isSucGap(cord, nw, revscomp_const, gap_size))
    {
        cord = previousWindow(f1[strand2], f2[genomeId2], cord, score, window_threshold);
        if (cord)
        {
            //TODO do another round previousWindow if cord = 0.
            appendValue (tmp, cord);
        }
        else
        {
            break;
        }
    }
    if (!empty (tmp))
    {
        for (unsigned i = 0; i < length(tmp) / 2; i++)
        {
            std::swap (tmp[i], tmp[length(tmp) - i - 1]);
        }
        insert(cords, kk + len, tmp);
        len += length(tmp);
    }
    return len;
}

/*
 * path for double strand 
 */
 bool path_dst(typename Iterator<String<uint64_t> >::Type hitBegin, 
               typename Iterator<String<uint64_t> >::Type hitEnd, 
               StringSet<FeaturesDynamic> & f1,
               StringSet<FeaturesDynamic> & f2, 
               String<uint64_t> & cords,
               float const & cordLenThr
               )
{
    typename Iterator<String<uint64_t> >::Type it = hitBegin;
    unsigned preBlockPtr;
    float score = 0;
    if(initCord(it, hitEnd, preBlockPtr, cords))
    {
        do{
            uint64_t strand = _DefaultCord.getCordStrand(back(cords));
            extendWindow(f1[strand], f2[get_cord_id(back(cords))], cords, score, strand);
        }
        while (nextCord(it, hitEnd, preBlockPtr, cords, cordLenThr, score));
        return endCord(cords, preBlockPtr, cordLenThr, score);   
    }
    //std::cout << "[]::path_dist::cord " 
    return false;
}

//End all mapper module
//============================================
