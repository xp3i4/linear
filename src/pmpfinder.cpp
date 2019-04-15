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

const unsigned scriptStep=16;
const unsigned scriptBit=4;
const unsigned scriptWindow=5; //script_length = 2^scriptWindow
const unsigned scriptWindow2 = scriptWindow << 1;
const unsigned scriptWindow3 = scriptWindow2 + scriptWindow;
const int scriptCount[5] = {1, 1<<scriptWindow, 1 <<(scriptWindow * 2), 0, 0};
const int scriptMask = (1 << scriptWindow) - 1;
const int scriptMask2 = scriptMask << scriptWindow;
const int scriptMask3 = scriptMask2 << scriptWindow;

const uint64_t hmask = (1ULL << 20) - 1;
const unsigned windowThreshold = 36; // 36;

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


//======HIndex getIndexMatch()

int _scriptDist(int const & s1, int const & s2)
{
    int res = std::abs((s1 & scriptMask) - (s2 & scriptMask)) +
              std::abs(((s1 >> scriptWindow) & scriptMask) -
                       ((s2 >> scriptWindow) & scriptMask)) + 
              std::abs((s1>>scriptWindow2) - (s2>>scriptWindow2));
    return res;
}

void createFeatures(TIter5 const & itBegin, TIter5 const & itEnd, String<short> & f)
{
    unsigned next = 1;
    unsigned window = 1 << scriptWindow;
    resize (f, ((itEnd - itBegin -window) >> scriptBit) + 1);
    f[0] = 0;
    for (unsigned k = 0; k < window; k++)
    {
        f[0] += scriptCount[ordValue(*(itBegin + k))];
    }
    for (unsigned k = scriptStep; k < itEnd - itBegin - window ; k+=scriptStep) 
    {
        f[next] = f[next - 1];
        for (unsigned j = k - scriptStep; j < k; j++)
        {
            f[next] += scriptCount[ordValue(*(itBegin + j + window))] - scriptCount[ordValue(*(itBegin + j))];
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

/**
 * Parallel 
 */
 void createFeatures(TIter5 const & itBegin, TIter5 const & itEnd, String<short> & f, unsigned threads)
{
    unsigned window = 1 << scriptWindow;
    resize (f, ((itEnd - itBegin -window) >> scriptBit) + 1);
#pragma omp parallel
{
    unsigned thd_id = omp_get_thread_num();
    uint64_t thd_begin;
    uint64_t thd_end;
    uint64_t range = (itEnd - itBegin - window - scriptStep) / scriptStep;
    parallelParm_Static(range, threads, 
                        thd_id,  thd_begin, thd_end);
    uint64_t next = thd_begin;
    thd_begin *= scriptStep;
    thd_end *= scriptStep;
    f[next] = 0;
    for (unsigned k = thd_begin; k < thd_begin + window; k++)
    {
        f[next] += scriptCount[ordValue(*(itBegin + k))];
    }
    next++;
    for (unsigned k = thd_begin + scriptStep; k < thd_end ; k+=scriptStep) 
    {
        f[next] = f[next - 1];
        for (unsigned j = k - scriptStep; j < k; j++)
        {
            f[next] += scriptCount[ordValue(*(itBegin + j + window))] - scriptCount[ordValue(*(itBegin + j))];
        }
        next++;
    }
}
}

/**
 * Parallel
 */
 void createFeatures(StringSet<String<Dna5> > & seq, StringSet<String<short> > & f, unsigned threads)
{
    resize(f, length(seq));
    for (unsigned k = 0; k < length(seq); k++)
        createFeatures(begin(seq[k]), end(seq[k]), f[k], threads);
}

 void createFeatures(StringSet<String<Dna5> > & seq, StringSet<String<short> > & f)
{
    resize(f, length(seq));
    for (unsigned k = 0; k < length(seq); k++)
        createFeatures(begin(seq[k]), end(seq[k]), f[k]);
}

 unsigned _windowDist(Iterator<String<short> >::Type const & it1, 
                      Iterator<String<short> >::Type const & it2)
{
    return _scriptDist(*it1, *it2) 
         + _scriptDist(*(it1 + 2), *(it2 + 2)) 
         + _scriptDist(*(it1 + 4), *(it2 + 4)) 
         + _scriptDist(*(it1 + 6), *(it2 + 6)) 
         + _scriptDist(*(it1 + 8), *(it2 + 8)) 
         + _scriptDist(*(it1 + 10), *(it2 + 10));
}

 bool nextCord(String<uint64_t> & hit, unsigned & currentIt, String<uint64_t> & cord)
{
    uint64_t cordLY = get_cord_y(back(cord));
    while (++currentIt < length(hit)) 
    {
        uint64_t tmpCord = _DefaultCord.hit2Cord(hit[currentIt]);
        if(get_cord_y(tmpCord) > cordLY + window_delta)
        {
            appendValue(cord, tmpCord);
            return true;
        }
    }
    return false;
}

 bool initCord(String<uint64_t> & hit, unsigned & currentIt, String<uint64_t> & cord)
{
    currentIt = 0;
    if (empty(hit))
        return false;
    else
        appendValue(cord, _DefaultCord.hit2Cord(hit[0]));
    return true;
}

 bool previousWindow(String<short> & f1, 
                           String<short> & f2, 
                           String<uint64_t> & cord, 
                           uint64_t & strand)
{
    typedef uint64_t CordType;
    CordType genomeId = get_cord_id(back(cord));
    CordType x_suf = _DefaultCord.cord2Cell(get_cord_x(back(cord)));
    CordType y_suf = _DefaultCord.cord2Cell(get_cord_y(back(cord)));
    CordType x_min = 0;
    CordType y;
    if (y_suf < med || x_suf < sup)
        return false;
    else 
        y = y_suf - med;

    unsigned min = ~0;
    for (CordType x = x_suf - sup; x < x_suf - inf; x += 1) 
    {
        unsigned tmp = _windowDist(begin(f1) + y, begin(f2) + x);
        
        if (tmp < min)
        {
            min = tmp;
            x_min = x;
        }
    }
    if (min > windowThreshold)
        return false;
    else 
    {
        if ( x_suf - x_min > med)
        {
            appendValue(cord, _DefaultCord.createCord(create_id_x(genomeId, _DefaultCord.cell2Cord(x_suf - med)),  _DefaultCord.cell2Cord(x_suf - x_min - med + y), strand));
        }
        else
            appendValue(cord, _DefaultCord.createCord(create_id_x(genomeId, _DefaultCord.cell2Cord(x_min)), _DefaultCord.cell2Cord(y), strand));
    }
    return true;
}

 bool previousWindow(String<short> & f1, 
                           String<short> & f2, 
                           String<uint64_t> & cord, 
                           float & score, 
                           uint64_t & strand)
{
    typedef uint64_t CordType;
    CordType genomeId = get_cord_id(back(cord));
    CordType x_suf = _DefaultCord.cord2Cell(get_cord_x(back(cord)));
    CordType y_suf = _DefaultCord.cord2Cell(get_cord_y(back(cord)));
    CordType x_min = 0;
    CordType y;
    if (y_suf < med || x_suf < sup)
        return false;
    else 
       y = y_suf - med;

    unsigned min = ~0;
    for (CordType x = x_suf - sup; x < x_suf - inf; x += 1) 
    {
        unsigned tmp = _windowDist(begin(f1) + y, begin(f2) + x);
        if (tmp < min)
        {
            min = tmp;
            x_min = x;
        }
    }
    if (min > windowThreshold)
        return false;    
    else 
    {
        if ( x_suf - x_min > med)
        {
            appendValue(cord, _DefaultCord.createCord(create_id_x(genomeId, _DefaultCord.cell2Cord(x_suf - med)),  _DefaultCord.cell2Cord(x_suf - x_min - med + y), strand));
        }
        else
            appendValue(cord, _DefaultCord.createCord(create_id_x(genomeId, _DefaultCord.cell2Cord(x_min)), _DefaultCord.cell2Cord(y), strand));
    }
    score += min;
    return true;
}

 uint64_t previousWindow(String<short> & f1, 
                               String<short> & f2, 
                               uint64_t cordx,
                               uint64_t cordy,
                               uint64_t strand,
                               unsigned window_threshold = windowThreshold )
{
    typedef uint64_t CordType;
    CordType genomeId = _getSA_i1(cordx);
    CordType x_suf = _DefaultCord.cord2Cell(_getSA_i2(cordx));
    CordType y_suf = _DefaultCord.cord2Cell(cordy);
    CordType x_min = 0;
    CordType y;
    
    if (y_suf < med || x_suf < sup)
        return 0;
    else 
        y = y_suf - med;

    unsigned min = ~0;
    for (CordType x = x_suf - sup; x < x_suf - inf; x += 1) 
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
            return _DefaultCord.createCord(create_id_x(genomeId, _DefaultCord.cell2Cord(x_suf - med)),  _DefaultCord.cell2Cord(x_suf - x_min - med + y), strand);
        }
        else
        {
            return _DefaultCord.createCord(create_id_x(genomeId, _DefaultCord.cell2Cord(x_min)), _DefaultCord.cell2Cord(y), strand);
        } 
    }
    return 0;
}

 uint64_t previousWindow(String<short> & f1, String<short> & f2, uint64_t cord, unsigned window_threshold)
{
    return previousWindow(f1, 
                          f2, 
                          _DefaultCord.getCordX(cord),
                          get_cord_y(cord), 
                          _DefaultCord.getCordStrand(cord), window_threshold);
}

 bool nextWindow(String<short> &f1, 
                       String<short> & f2, 
                       String<uint64_t> & cord, 
                       uint64_t & strand)
{
    typedef uint64_t CordType;
    CordType genomeId = get_cord_id(back(cord));
    CordType x_pre = _DefaultCord.cord2Cell(get_cord_x(back(cord)));
    CordType y_pre = _DefaultCord.cord2Cell(get_cord_y(back(cord)));
    CordType x_min = 0;
    CordType y;
    if (y_pre + sup * 2> length(f1) || x_pre + sup * 2 > length(f2))
        return false;
    else 
        y = y_pre + med;
    
    unsigned min = ~0;
    for (CordType x = x_pre + inf; x < x_pre + sup; x += 1) 
    {

        unsigned tmp = _windowDist(begin(f1) + y, begin(f2) + x);
        if (tmp < min)
        {
            min = tmp;
            x_min = x;
        }
    }

    if (min > windowThreshold)
       return false;
    else 
        if ( x_min - x_pre > med)
        {
            appendValue(cord, _DefaultCord.createCord(create_id_x(genomeId, _DefaultCord.cell2Cord(x_pre + med)),  _DefaultCord.cell2Cord(x_pre + med - x_min + y), strand));
        }
        else
        {
            appendValue(cord, _DefaultCord.createCord(create_id_x(genomeId, _DefaultCord.cell2Cord(x_min)), _DefaultCord.cell2Cord(y), strand));
        }
    return true;
}

 bool nextWindow(String<short> &f1, 
                       String<short> & f2, 
                       String<uint64_t> & cord, 
                       float & score, 
                       uint64_t & strand)
{
    typedef uint64_t CordType;
    CordType genomeId = get_cord_id(back(cord));
    CordType x_pre = _DefaultCord.cord2Cell(get_cord_x(back(cord)));
    CordType y_pre = _DefaultCord.cord2Cell(get_cord_y(back(cord)));
    CordType x_min = 0;
    CordType y;
    unsigned min = ~0;
    
    if (y_pre + sup * 2 > length(f1) || x_pre + sup * 2> length(f2))
        return false;
    else 
        y = y_pre + med;
    
    for (CordType x = x_pre + inf; x < x_pre + sup; x += 1) 
    {
        unsigned tmp = _windowDist(begin(f1) + y, begin(f2) + x);
        if (tmp < min)
        {
            min = tmp;
            x_min = x;
        }
    }
    if (min > windowThreshold)
       return false;
    else 
        if ( x_min - x_pre > med)
        {
            appendValue(cord, _DefaultCord.createCord(create_id_x(genomeId, _DefaultCord.cell2Cord(x_pre + med)),  _DefaultCord.cell2Cord(x_pre + med - x_min + y), strand));
        }
        else
        {
            appendValue(cord, _DefaultCord.createCord(create_id_x(genomeId, _DefaultCord.cell2Cord(x_min)), _DefaultCord.cell2Cord(y), strand));
        }
    score += min;
    return true;
}

 uint64_t nextWindow(String<short> & f1, 
                           String<short> & f2, 
                           uint64_t cordx,
                           uint64_t cordy,
                           uint64_t strand,
                           unsigned window_threshold = windowThreshold
                          )
{
    typedef uint64_t CordType;
    CordType genomeId = _getSA_i1(cordx);
    CordType x_pre = _DefaultCord.cord2Cell(_getSA_i2(cordx));
    CordType y_pre = _DefaultCord.cord2Cell(cordy);
    CordType x_min = 0;
    CordType y;
    unsigned min = ~0;
    
    if (y_pre + sup * 2 > length(f1) || x_pre + sup * 2> length(f2))
        return 0;
    else 
        y = y_pre + med;
    
    for (CordType x = x_pre + inf; x < x_pre + sup; x += 1) 
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
        if ( x_min - x_pre > med)
        {
            return _DefaultCord.createCord(create_id_x(genomeId, _DefaultCord.cell2Cord(x_pre + med)),  _DefaultCord.cell2Cord(x_pre + med - x_min + y), strand);
        }
        else
        {
            return _DefaultCord.createCord(create_id_x(genomeId, _DefaultCord.cell2Cord(x_min)), _DefaultCord.cell2Cord(y), strand);
        }
    return 0;
}

 uint64_t nextWindow(String<short> & f1, String<short> & f2, uint64_t cord, unsigned window_threshold)
{
    return nextWindow(f1, f2, _DefaultCord.getCordX(cord), get_cord_y(cord), _DefaultCord.getCordStrand(cord), window_threshold);
}

 bool extendWindow(String<short> &f1, String<short> & f2, String<uint64_t> & cord, uint64_t strand)
{
    uint64_t preCord = (length(cord)==1)?0:back(cord);
    unsigned len = length(cord) - 1;
    while (preCord <= back(cord) && previousWindow(f1, f2, cord, strand)){}
    for (unsigned k = len; k < ((length(cord) + len) >> 1); k++) 
    {
        std::swap(cord[k], cord[length(cord) - k + len - 1]);
    }
    while (nextWindow(f1, f2, cord, strand)){}
    return true;
}

 bool path(String<Dna5> & read, String<uint64_t> hit, StringSet<String<short> > & f2, String<uint64_t> & cords)
{
    String<short> f1;
    unsigned currentIt = 0;
    if(!initCord(hit, currentIt, cords))
        return false;
    createFeatures(begin(read), end(read), f1);
    unsigned genomeId = get_cord_id(cords[0]);
    while (get_cord_y(back(cords)) < length(read) - window_size)
    {
        extendWindow(f1, f2[genomeId], cords, _DefaultCord.getCordStrand(back(cords)));
        if(!nextCord(hit, currentIt, cords))
                return false;
    }
    return true;
}

void path(String<uint64_t> & hits, StringSet<String<Dna5> > & reads, StringSet<String<Dna5> > & genomes, StringSet<String<uint64_t> > & cords)
{
    StringSet<String<short> > f2;
    createFeatures(genomes, f2);
    std::cerr << "raw mapping... " << std::endl;
    for (unsigned k = 0; k < length(reads); k++)
    {
        path(reads[k], hits[k], f2, cords[k]);
    }
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

/**================================================================
 *  The following part implements different method of mapping 
 */
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

//===!Note:Need to put this parameterin the mapper threshold
/*
template <typename TDna, typename TSpec>
 unsigned getIndexMatchAll(typename LIndex & index,
                              typename String<Dna5> const & read,
                              uint64_t* const set,
                              unsigned & len,
                              MapParm & mapParm
                             )
{    
    unsigned block = (mapParm.blockSize < length(read))?mapParm.blockSize:length(read);
    unsigned dt = block * (mapParm.alpha / (1 - mapParm.alpha));
    hashInit(index.shape, begin(read));
    for (unsigned h=0; h <= length(read) - block; h += dt)
    {
        hashInit(index.shape, begin(read) + h);
        for (unsigned k = h; k < h + block; k++)
        {
            hashNext(index.shape, begin(read) + k);
            uint64_t pre = ~0;
            uint64_t pos = getXYDir(index, index.shape.XValue, index.shape.YValue);
            if (_DefaultHs.getHsBodyY(index.ysa[std::min(pos + mapParm.delta, length(index.ysa) - 1)]) ^ index.shape.YValue)
            {
                while (_DefaultHs.isBodyYEqual(index.ysa[pos], index.shape.YValue))
                {
//!Note: needs change
//!Note: the sa is in reverse order in hindex. this is different from the generic index
                    if (pre - index.ysa[pos] > mapParm.kmerStep)
                    {
                        set[len++] = (_DefaultHs.getHsBodyS(index.ysa[pos]-k) << 20)+k;
                        pre = index.ysa[pos];
                    }
                    ++pos;
                }
            }
        }
    }
    return 0;
}
*/
/*
template <typename TDna, typename TSpec>
 unsigned getIndexMatchAll(typename LIndex & index,
                              typename String<Dna5> const & read,
                              uint64_t* const set,
                              unsigned & len,
                              MapParm & mapParm
                             )
{   
    typedef typename LIndex TIndex;
    typedef typename TIndex::TShape PShape;
    
    unsigned block = (mapParm.blockSize < length(read))?mapParm.blockSize:length(read);
    unsigned dt = block * (mapParm.alpha / (1 - mapParm.alpha));
    PShape shape;
    hashInit(shape, begin(read));
    for (unsigned h=0; h <= length(read) - block; h += dt)
    {
        hashInit(shape, begin(read) + h);
        for (unsigned k = h; k < h + block; k++)
        {
            hashNext(shape, begin(read) + k);
            uint64_t pre = ~0;
            uint64_t pos = getXYDir(index, shape.XValue, shape.YValue);
            if (_DefaultHs.getHsBodyY(index.ysa[std::min(pos + mapParm.delta, length(index.ysa) - 1)]) ^ shape.YValue)
            {
                while (_DefaultHs.isBodyYEqual(index.ysa[pos], shape.YValue))
                {
//!Note: needs change
//!Note: the sa is in reverse order in hindex. this is different from the generic index
                    if (pre - index.ysa[pos] > mapParm.kmerStep)
                    {
                        set[len++] = (_DefaultHs.getHsBodyS(index.ysa[pos]-k) << 20)+k;
                        pre = index.ysa[pos];
                    }
                    ++pos;
                }
            }
        }
    }
    return 0;
}
*/
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
    //<<debug
        for (int ii = 0; ii < anchors.length(); ii++)
        {
            uint64_t mask1 = (1ULL << 20) - 1;
            uint64_t mask2 = (1ULL << 40) - 1;
            uint64_t tmp_cord = _DefaultCord.hit2Cord(anchors[ii]);

            std::cout << "gaml1 " << ii << " " << get_cord_y(tmp_cord) << " " << get_cord_x(anchors[ii]) << " " << get_cord_x(tmp_cord) << "\n";
        }
    //>>debug
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

 bool endCord( String<uint64_t> & cord,
                     unsigned & preCordStart
                   )
{
    _DefaultHit.setBlockEnd(back(cord));
    _DefaultCord.setMaxLen(cord, length(cord) - preCordStart);
    return true;
}
/*
 * endCord for single strand index
 *
 bool endCord( String<uint64_t> & cord,
                     unsigned & preCordStart,
                     float const & cordThr,
                     uint64_t const & strand,
                     float & score
                   )
{
    _DefaultCord.setMaxLen(cord, length(cord) - preCordStart);   
    if (length(cord) - preCordStart > cordThr)
    {

        _DefaultHit.setBlockEnd(back(cord));
        _DefaultHit.setBlockStrand(cord[preCordStart], strand);
    }
    else
    {
        erase(cord, preCordStart, length(cord));
    }
    score = 0;
    return true;
}
*/

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
 bool nextCord(typename Iterator<String<uint64_t> >::Type & it, 
                     typename Iterator<String<uint64_t> >::Type const & hitEnd,
                     unsigned & preCordStart,
                     String<uint64_t> & cord)
{
//TODO: add maxlen of anchor to the first node in cord;
    if (it >= hitEnd)
        return false;
    while(!_DefaultHit.isBlockEnd(*(it - 1)))
    {
        if(get_cord_y(*(it))>get_cord_y(back(cord)) +  window_delta) 
        {
            appendValue(cord, _DefaultCord.hit2Cord(*(it)));
            ++it;
            return true;
        }
        ++it;
    }
    _DefaultHit.setBlockEnd(back(cord));
    if(it < hitEnd)
    {
        _DefaultCord.setMaxLen(cord, length(cord) - preCordStart);
        preCordStart = length(cord);
        appendValue(cord, _DefaultCord.hit2Cord(*(it)));
        ++it;
        return true;
    }
    else
        return false;
}
*/

/*
 * nextCord for single strand sequence
 *
 bool nextCord(typename Iterator<String<uint64_t> >::Type & it, 
                     typename Iterator<String<uint64_t> >::Type const & hitEnd, 
                     unsigned & preCordStart,
                     String<uint64_t> & cord,
                     float const & cordThr,
                     uint64_t const & strand,
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
            appendValue(cord, _DefaultCord.hit2Cord(*(it)));
            ++it;
            return true;
        }
        ++it;
    }
    _DefaultHit.setBlockEnd(back(cord));
    _DefaultHit.setBlockStrand(cord[preCordStart], strand);
    if(it < hitEnd)
    {
        if (length(cord) - preCordStart < cordThr )//|| std::abs(b1 - b2) < 10)
        {
            erase(cord, preCordStart, length(cord));
        }
        else
        {
            _DefaultCord.setMaxLen(cord, length(cord) - preCordStart);
        }
        preCordStart = length(cord);
        appendValue(cord, _DefaultCord.hit2Cord(*(it)));
        ++it;
        return true;
    }
    else
        return false;
    (void) score;
}*/

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

 bool extendWindowAll(String<short> &f1, 
                            String<short> & f2, 
                            String<uint64_t> & cord, 
                            float & score, 
                            uint64_t & strand)
{
    uint64_t preCordY = (_DefaultHit.isBlockEnd(cord[length(cord) - 2]))?
    0:get_cord_y(back(cord)) + window_delta;
    unsigned len = length(cord) - 1;
    while (preCordY<= get_cord_y(back(cord)) && 
        previousWindow(f1, f2, cord, score, strand)){}
    for (unsigned k = len; k < ((length(cord) + len) >> 1); k++) 
    {
        std::swap(cord[k], cord[length(cord) - k + len - 1]);
    }
    while (nextWindow(f1, f2, cord, score, strand)){}
    //std::cout << "[]::extendWindow " << k << "\n";
    return true;
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
 int extendPatch(StringSet<String<short> > & f1, 
                 StringSet<String<short> > & f2, 
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
        cord = nextWindow (f1[strand1], f2[genomeId1], cord, window_threshold);
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
        cord = previousWindow(f1[strand2], f2[genomeId2], cord, window_threshold);
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
 bool pathAll(String<Dna5> & read, 
                 typename Iterator<String<uint64_t> >::Type hitBegin, 
                 typename Iterator<String<uint64_t> >::Type hitEnd, 
                 StringSet<String<int> > & f2, 
                 String<uint64_t> & cords
                   )
{
    String<int> f1;
    typename Iterator<String<uint64_t> >::Type it = hitBegin;
    unsigned preBlockPtr;
    createFeatures(begin(read), end(read), f1);
    float score;
    if(initCord(it, hitEnd, preBlockPtr, cords))
    {
        do{
            extendWindowAll(f1, f2[get_cord_id(back(cords))], cords, score;
        }
        while (nextCord(it, hitEnd, preBlockPtr, cords));
        return endCord(cords, preBlockPtr);   
    }
    return false;
}

 bool pathAll(String<Dna5> & read, 
                 typename Iterator<String<uint64_t> >::Type hitBegin, 
                 typename Iterator<String<uint64_t> >::Type hitEnd, 
                 StringSet<String<int> > & f2, 
                 String<uint64_t> & cords,
                 float const & cordThr, 
                 uint64_t const & strand
                   )
{
    String<int> f1;
    typename Iterator<String<uint64_t> >::Type it = hitBegin;
    unsigned preBlockPtr;
    float cordLenThr = length(read) * cordThr / window_size;
    float score = 0;
    createFeatures(begin(read), end(read), f1);
    if(initCord(it, hitEnd, preBlockPtr, cords))
    {
        do{
            extendWindowAll(f1, f2[get_cord_id(back(cords))], cords, score;
        }
        while (nextCord(it, hitEnd, preBlockPtr, cords, cordLenThr, strand, score));
        return endCord(cords, preBlockPtr, cordLenThr, strand, score);   
    }
    return false;
}
*/

/*
 * path functoin for single strand
 bool pathAll(String<Dna5> & read, 
                 typename Iterator<String<uint64_t> >::Type hitBegin, 
                 typename Iterator<String<uint64_t> >::Type hitEnd, 
                 StringSet<String<int> > & f2, 
                 String<uint64_t> & cords,
                 float const & cordThr, 
                 uint64_t const & strand
                   )
{
    //std::cerr << "[pathall]\n";
    String<int> f1;
    typename Iterator<String<uint64_t> >::Type it = hitBegin;
    unsigned preBlockPtr;
    float cordLenThr = length(read) * cordThr / window_size;
    float score = 0;
    createFeatures(begin(read), end(read), f1);
    if(initCord(it, hitEnd, preBlockPtr, cords))
    {
        do{
            extendWindowAll(f1, f2[get_cord_id(back(cords))], cords, score;
        }
        while (nextCord(it, hitEnd, preBlockPtr, cords, cordLenThr, strand, score));
        return endCord(cords, preBlockPtr, cordLenThr, strand, score);   
    }
    return false;
}
*/

/*
 * path for double strand 
 */
 bool path_dst(
                 typename Iterator<String<uint64_t> >::Type hitBegin, 
                 typename Iterator<String<uint64_t> >::Type hitEnd, 
                 StringSet<String<short> > & f1,
                 StringSet<String<short> > & f2, 
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
            extendWindowAll(f1[strand], f2[get_cord_id(back(cords))], cords, score, strand);
        }
        while (nextCord(it, hitEnd, preBlockPtr, cords, cordLenThr, score));
        return endCord(cords, preBlockPtr, cordLenThr, score);   
    }
    //std::cout << "[]::path_dist::cord " 
    return false;
}

/*
 * this mapping function use the hash value of double strand 
 * that allow to stream the sequence of double strand for only once
 */
int rawMap_dst( LIndex   & index,
             StringSet<String<Dna5> > & reads,
             StringSet<String<Dna5> > & genomes, 
            MapParm & mapParm,
            StringSet<String<uint64_t> > & cords,
            unsigned & threads)
{
  
    typedef String<Dna5> Seq;
    double time=sysTime();
    float senThr = mapParm.senThr / window_size;
    float cordThr = mapParm.cordThr / window_size;
    MapParm complexParm = mapParm;
    complexParm.alpha = complexParm.alpha2;
    complexParm.listN = complexParm.listN2;
    std::cerr << "[rawMap_dst] Raw Mapping \n"; 

    
    StringSet<String<short> > f2;
    double time2 = sysTime();
    createFeatures(genomes, f2, threads);
    std::cerr << "init1 " << sysTime() - time2 << "\n";
    
    
#pragma omp parallel
{
    unsigned size2 = length(reads) / threads;
    unsigned ChunkSize = size2;
    Seq comStr;
    //Anchors anchors(Const_::_LLTMax, AnchorBase::size);
    Anchors anchors;
    String<uint64_t> crhit;
    StringSet<String<uint64_t> >  cordsTmp;
    StringSet< String<short> > f1;
    //printf ("done\n");
    unsigned thd_id =  omp_get_thread_num();
    if (thd_id < length(reads) - size2 * threads)
    {
        ChunkSize = size2 + 1;
    }
    resize(cordsTmp, ChunkSize);
    resize(f1, 2);
    unsigned c = 0;

    //printf ("done\n");
    #pragma omp for
    for (unsigned j = 0; j < length(reads); j++)
    //for (unsigned j = 14456; j < 14457; j++)
    {
        //printf ("[debug] id %d %d %d\n", thd_id, j, length(reads[j]));
        if (length(reads[j]) >= mapParm.minReadLen)
        {
            float cordLenThr = length(reads[j]) * cordThr;
            _compltRvseStr(reads[j], comStr);
            createFeatures(begin(reads[j]), end(reads[j]), f1[0]);
            createFeatures(begin(comStr), end(comStr), f1[1]);
            anchors.init(1);
            clear(crhit);
            mnMapReadList(index, reads[j], anchors, mapParm, crhit);
            path_dst(begin(crhit), end(crhit), f1, f2, cordsTmp[c], cordLenThr);
            //printf("done1\n");
            if (_DefaultCord.getMaxLen(cordsTmp[c]) < length(reads[j]) * senThr)// && 
               //_DefaultCord.getMaxLen(cordsTmp[c]) > 0)
            {
                clear(cordsTmp[c]);
                anchors.init(1);
                clear(crhit);
                mnMapReadList(index, reads[j], anchors, complexParm, crhit);
                path_dst(begin(crhit), end(crhit), f1, f2, cordsTmp[c], cordLenThr);
            }   
        }
        //printf("#done\n");
        c += 1;
    } 
    #pragma omp for ordered
    for (unsigned j = 0; j < threads; j++)
        #pragma omp ordered
        {
            append(cords, cordsTmp);
        }
}
    std::cerr << "    End raw mapping. Time[s]: " << sysTime() - time << std::endl;
    return 0;
}


//End all mapper module
//============================================
