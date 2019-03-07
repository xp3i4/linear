// ==========================================================================
//                          Mapping SMRT reads
// ==========================================================================
// Copyright (c) 2006-2016, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: cxpan <chenxu.pan@fu-berlin.de>
// ==========================================================================
#include <iostream>
#include <seqan/sequence.h>
#include <seqan/stream.h>
#include "base.h"
#include "index_util.h"
#include "pmpfinder.h"

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
        valueMask_dstr(valueMask | flag_strand)
{}
    
inline uint64_t 
Cord::getCordX(uint64_t const & cord, 
               unsigned const & bit  = _DefaultCordBase.bit,
               uint64_t const & mask = _DefaultCordBase.maskx) const
{
    return (cord >> bit) & mask; 
}

inline uint64_t 
Cord::getCordY(uint64_t const & cord, 
               uint64_t const & mask = _DefaultCordBase.mask) const 
{
    return cord & mask;
}

inline uint64_t 
Cord::createCord(uint64_t const & x, 
                 uint64_t const & y, 
                 uint64_t const & strand,
                 unsigned const & bit = _DefaultCordBase.bit, 
                 unsigned const & bit2 = _DefaultCordBase.flag_bit) const
{
    return (x << bit) + y + (strand << bit2);
}

inline uint64_t 
Cord::hit2Cord(uint64_t const & hit, 
               unsigned const & bit = _DefaultCordBase.bit, 
               uint64_t const & mask = _DefaultCordBase.mask,
               uint64_t const & mask2 = _DefaultCordBase.valueMask
              ) const
{
    return (hit + ((hit & mask) << bit)) & mask2;
}

inline uint64_t 
Cord::hit2Cord_dstr(uint64_t const & hit, 
               unsigned const & bit = _DefaultCordBase.bit, 
               uint64_t const & mask = _DefaultCordBase.mask,
               uint64_t const & mask2 = _DefaultCordBase.valueMask_dstr
              ) const
{
    return (hit + ((hit & mask) << bit)) & mask2;
}

inline uint64_t Cord::cord2Cell(uint64_t const & cord, 
                unsigned const & bit = _DefaultCordBase.cell_bit) const
{
    return cord >> bit;
}

inline uint64_t Cord::cell2Cord(uint64_t const & cell, 
                unsigned const & bit = _DefaultCordBase.cell_bit) const
{
    return cell << bit;
}

inline void Cord::setCordEnd(uint64_t & cord,
            typename CordBase::Flag const & strand = _DefaultCordBase.flag_strand,
            typename CordBase::Flag const & end = _DefaultCordBase.flag_end)
{
    cord |= strand | end;
}

inline typename CordBase::Flag 
Cord::getCordStrand(uint64_t const & cord,
            unsigned const & strand = _DefaultCordBase.flag_bit) const
{
    return (cord >> strand) & 1ULL;
}

inline typename CordBase::Flag 
Cord::isCordEnd(uint64_t const & cord,
                typename CordBase::Flag const & end = _DefaultCordBase.flag_end) const
{
    return cord & end;
}

inline void Cord::setMaxLen(String<uint64_t> & cord, uint64_t const & len, uint64_t const & mask)
{
    if (len > (cord[0] & mask))
        cord[0] = len + ((cord[0]) & (~mask));
}

inline uint64_t Cord::getMaxLen(String<uint64_t> const & cord, uint64_t const & mask)
{
    if (empty(cord))
        return 0;
    return cord[0] & mask;
}

inline uint64_t Cord::shift(uint64_t const & val, int64_t x, int64_t y, unsigned const & bit) //add x and y
{
    if (x < 0)
        return val - ((-x) << bit) + y;
    else
        return val + (x << bit) + y;
}

inline bool Cord::isCordsOverlap(uint64_t & val1, uint64_t & val2, int64_t thd)
{
    int64_t dx = _DefaultCord.getCordX(val2 - val1);
    int64_t dy = _DefaultCord.getCordY(val2 - val1);
    //std::cout << "[]::isCordsOverlap " << dx << " " << dy << " " << _DefaultCord.getCordX(val2) << " " << _DefaultCord.getCordX(val1) << " " << thd << "\n";
    return (dx >= 0) && (dx < thd) && (dy >= 0) && (dy < thd);
}

inline bool Cord::isBlockEnd(uint64_t & val, uint64_t const & flag)
{
    return val & flag;
}

inline uint64_t get_cord_x (uint64_t val) {return _getSA_i2(_DefaultCord.getCordX(val));}
inline uint64_t get_cord_y (uint64_t val) {return _DefaultCord.getCordY(val);}
inline uint64_t get_cord_strand (uint64_t val) {return _DefaultCord.getCordStrand(val);}

inline void cmpRevCord(uint64_t val1, 
                    uint64_t val2,
                    uint64_t & cr_val1,
                    uint64_t & cr_val2,
                    uint64_t read_len)
{
    cr_val1 = (val1 - get_cord_y(val1) + read_len - get_cord_y(val2)) ^ _DefaultCordBase.flag_strand;
    cr_val2 = (val2 - get_cord_y(val1) + read_len - get_cord_y(val2)) ^ _DefaultCordBase.flag_strand;
}
inline uint64_t set_cord_xy (uint64_t val, uint64_t x, uint64_t y)
{
    return (val & (~_DefaultCordBase.valueMask)) + (x << _DefaultCordBase.bit) + y;
}


//======HIndex getIndexMatch()

inline int _scriptDist(int const & s1, int const & s2)
{
    int res = std::abs((s1 & scriptMask)
            - (s2 & scriptMask)) 
            + std::abs(((s1 >> scriptWindow) & scriptMask) 
            - ((s2 >> scriptWindow) & scriptMask)) 
            + std::abs((s1>>scriptWindow2) - (s2>>scriptWindow2));
    return res;
}

inline void createFeatures(TIter5 const & itBegin, TIter5 const & itEnd, String<short> & f)
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

inline uint64_t parallelParm_Static(uint64_t range, 
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
inline void createFeatures(TIter5 const & itBegin, TIter5 const & itEnd, String<short> & f, unsigned threads)
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
inline void createFeatures(StringSet<String<Dna5> > & seq, StringSet<String<short> > & f, unsigned threads)
{
    resize(f, length(seq));
    for (unsigned k = 0; k < length(seq); k++)
        createFeatures(begin(seq[k]), end(seq[k]), f[k], threads);
}

inline void createFeatures(StringSet<String<Dna5> > & seq, StringSet<String<short> > & f)
{
    resize(f, length(seq));
    for (unsigned k = 0; k < length(seq); k++)
        createFeatures(begin(seq[k]), end(seq[k]), f[k]);
}

inline unsigned _windowDist(Iterator<String<short> >::Type const & it1, 
                            Iterator<String<short> >::Type const & it2)
{
    return _scriptDist(*it1, *it2) 
         + _scriptDist(*(it1 + 2), *(it2 + 2)) 
         + _scriptDist(*(it1+4),*(it2+4)) 
         + _scriptDist(*(it1 + 6), *(it2 + 6)) 
         + _scriptDist(*(it1+8), *(it2+8)) 
         + _scriptDist(*(it1 + 10), *(it2 + 10));
}

inline bool nextCord(String<uint64_t> & hit, unsigned & currentIt, String<uint64_t> & cord)
{
    uint64_t cordLY = _DefaultCord.getCordY(back(cord));
    while (++currentIt < length(hit)) 
    {
        uint64_t tmpCord = _DefaultCord.hit2Cord(hit[currentIt]);
        if(_DefaultCord.getCordY(tmpCord) > cordLY + window_delta)
        {
            appendValue(cord, tmpCord);
            return true;
        }
    }
    return false;
}

inline bool initCord(String<uint64_t> & hit, unsigned & currentIt, String<uint64_t> & cord)
{
    currentIt = 0;
    if (empty(hit))
        return false;
    else
        appendValue(cord, _DefaultCord.hit2Cord(hit[0]));
    return true;
}

inline bool previousWindow(String<short> & f1, 
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
            appendValue(cord, _DefaultCord.createCord(_createSANode(genomeId, _DefaultCord.cell2Cord(x_suf - med)),  _DefaultCord.cell2Cord(x_suf - x_min - med + y), strand));
        }
        else
            appendValue(cord, _DefaultCord.createCord(_createSANode(genomeId, _DefaultCord.cell2Cord(x_min)), _DefaultCord.cell2Cord(y), strand));
    }
    return true;
}

inline bool previousWindow(String<short> & f1, 
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
            appendValue(cord, _DefaultCord.createCord(_createSANode(genomeId, _DefaultCord.cell2Cord(x_suf - med)),  _DefaultCord.cell2Cord(x_suf - x_min - med + y), strand));
        }
        else
            appendValue(cord, _DefaultCord.createCord(_createSANode(genomeId, _DefaultCord.cell2Cord(x_min)), _DefaultCord.cell2Cord(y), strand));
    }
    score += min;
    return true;
}

inline uint64_t previousWindow(String<short> & f1, 
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
            return _DefaultCord.createCord(_createSANode(genomeId, _DefaultCord.cell2Cord(x_suf - med)),  _DefaultCord.cell2Cord(x_suf - x_min - med + y), strand);
        }
        else
        {
            return _DefaultCord.createCord(_createSANode(genomeId, _DefaultCord.cell2Cord(x_min)), _DefaultCord.cell2Cord(y), strand);
        } 
    }
    return 0;
}

inline uint64_t previousWindow(String<short> & f1, String<short> & f2, uint64_t cord, unsigned window_threshold)
{
    return previousWindow(f1, 
                          f2, 
                          _DefaultCord.getCordX(cord),
                          _DefaultCord.getCordY(cord), 
                          _DefaultCord.getCordStrand(cord), window_threshold);
}

inline bool nextWindow(String<short> &f1, 
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
            appendValue(cord, _DefaultCord.createCord(_createSANode(genomeId, _DefaultCord.cell2Cord(x_pre + med)),  _DefaultCord.cell2Cord(x_pre + med - x_min + y), strand));
        }
        else
        {
            appendValue(cord, _DefaultCord.createCord(_createSANode(genomeId, _DefaultCord.cell2Cord(x_min)), _DefaultCord.cell2Cord(y), strand));
        }
    return true;
}

inline bool nextWindow(String<short> &f1, 
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
            appendValue(cord, _DefaultCord.createCord(_createSANode(genomeId, _DefaultCord.cell2Cord(x_pre + med)),  _DefaultCord.cell2Cord(x_pre + med - x_min + y), strand));
        }
        else
        {
            appendValue(cord, _DefaultCord.createCord(_createSANode(genomeId, _DefaultCord.cell2Cord(x_min)), _DefaultCord.cell2Cord(y), strand));
        }
    score += min;
    return true;
}

inline uint64_t nextWindow(String<short> & f1, 
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
            return _DefaultCord.createCord(_createSANode(genomeId, _DefaultCord.cell2Cord(x_pre + med)),  _DefaultCord.cell2Cord(x_pre + med - x_min + y), strand);
        }
        else
        {
            return _DefaultCord.createCord(_createSANode(genomeId, _DefaultCord.cell2Cord(x_min)), _DefaultCord.cell2Cord(y), strand);
        }
    return 0;
}

inline uint64_t nextWindow(String<short> & f1, String<short> & f2, uint64_t cord, unsigned window_threshold)
{
    return nextWindow(f1, f2, _DefaultCord.getCordX(cord), _DefaultCord.getCordY(cord), _DefaultCord.getCordStrand(cord), window_threshold);
}

inline bool extendWindow(String<short> &f1, String<short> & f2, String<uint64_t> & cord, uint64_t strand)
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

inline bool path(String<Dna5> & read, String<uint64_t> hit, StringSet<String<short> > & f2, String<uint64_t> & cords)
{
    String<short> f1;
    unsigned currentIt = 0;
    if(!initCord(hit, currentIt, cords))
        return false;
    createFeatures(begin(read), end(read), f1);
    unsigned genomeId = _getSA_i1(_DefaultCord.getCordX(cords[0]));
    while (_DefaultCord.getCordY(back(cords)) < length(read) - window_size)
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
            if (_DefaultCord.getCordY(back(*it)) + window_size * 2 < length(read))
            {
                count++;
            }
        it++;
    }
}

/**================================================================
 *  The following part implements different method of mapping 
 */
inline void Hit::setBlockStart(uint64_t & val, uint64_t const & flag)
{
    val |= flag;
}

inline void Hit::setBlockBody(uint64_t & val, uint64_t const & flag)
{
    val &= (~flag);
}

inline bool Hit::isBlockStart(uint64_t & val, uint64_t const & flag)
{
    return val & flag;
}

inline void Hit::setBlockEnd(uint64_t & val, uint64_t const & flag)
{
    val |= flag;
}

inline void Hit::unsetBlockEnd(uint64_t & val, uint64_t const & flag)
{
    val &= ~flag;
}

inline void Hit::setBlockStrand(uint64_t & val, uint64_t const & strand, uint64_t const & flag)
{
    if (strand)
        val |= flag;
    else
        val &= ~flag;
}

inline bool Hit::isBlockEnd(uint64_t & val, uint64_t const & flag)
{
    return val & flag;
}

inline unsigned Hit::getStrand(uint64_t const & val, uint64_t const & flag)
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
              << _DefaultCord.getCordY(hit[k]) << "\n";
        if (_DefaultHit.isBlockEnd(hit[k]))
        {
            std::cout << "[P]::_printHit() end\n";
        }
    }
}

//===!Note:Need to put this parameterin the mapper threshold
/*
template <typename TDna, typename TSpec>
inline unsigned getIndexMatchAll(typename LIndex & index,
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
inline unsigned getIndexMatchAll(typename LIndex & index,
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
 * search double strand in one round
 */
inline unsigned getIndexMatchAll(LIndex & index,
                               String<Dna5> const & read,
                              String<uint64_t> & set,
                              MapParm & mapParm)
{   
    typedef LIndex::TShape PShape;
    unsigned dt = 0;
    PShape shape;
    uint64_t xpre = 0;
    hashInit(shape, begin(read));
    for (unsigned k = 0; k < length(read); k++)
    {
        hashNexth(shape, begin(read) + k);
        uint64_t pre = ~0;
        if (++dt == mapParm.alpha)
        {
            if(hashNextX(shape, begin(read) + k) ^ xpre)
            {
                xpre = shape.XValue;
                uint64_t pos = getXDir(index, shape.XValue, shape.YValue);
//!Note: The contition is different from single strand \
         which will slightly changes the senstivity.\
         Specifically, in the single strand index, if the size of the block is > mapParm.delta, \
         then the block of the ysa will not be used. \
         While in double strand index, the length of ysa of fixed value includes kmer \
         of f both + and - strands.

        // if (_DefaultHs.getHsBodyY(index.ysa[std::min(pos + mapParm.delta, length(index.ysa) - 1)]) ^ shape.YValue)
            //{
                //while (_DefaultHs.isBodyYEqual(index.ysa[pos], shape.YValue))
                uint64_t ptr = _DefaultHs.getHeadPtr(index.ysa[pos-1]);
                if (ptr < mapParm.delta)
                //if (_DefaultHs.getHeadPtr(index.ysa[pos-1]) < 1000000)
                {
          //          unsigned pr = _DefaultHs.getHeadPtr(index.ysa[pos-1]);
                    //while (_DefaultHs.isBody(index.ysa[pos]))
                    while ((_DefaultHs.getHsBodyY(index.ysa[pos]) == shape.YValue || _DefaultHs.getHsBodyY(index.ysa[pos]) == 0))
                    {
    //!Note: needs change
    //!Note: the sa is in reverse order in hindex. this is different from the generic index
                        if (_DefaultHs.getHsBodyS(pre - index.ysa[pos]) > mapParm.kmerStep)
                        {
                            //[COMT]::condition of complement reverse strand
                            if (((index.ysa[pos] & _DefaultHsBase.bodyCodeFlag) >>_DefaultHsBase.bodyCodeBit) ^ shape.strand)
                            {
                                
                                appendValue(set, (((_DefaultHs.getHsBodyS(index.ysa[pos]) - length(read) + 1 + k) << 20) +  (length(read) - 1 - k)) | _DefaultHitBase.flag2);
                                
                            }
                            //[comt]::condition of normal strand
                            else
                            {     

                                appendValue(set, ((_DefaultHs.getHsBodyS(index.ysa[pos]) - k) << 20) | k);
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
            dt = 0;
        }
    }
    return 0;
}

inline uint64_t getAnchorMatchAll(Anchors & anchors, unsigned const & readLen, MapParm & mapParm, String<uint64_t> & hit)
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

inline uint64_t getAnchorMatchFirst(Anchors & anchors, unsigned const & readLen, MapParm & mapParm, String<uint64_t> & hit)
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

inline uint64_t getAnchorMatchList(Anchors & anchors, unsigned const & readLen, MapParm & mapParm, String<uint64_t> & hit)
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
        return (list[0] >> 40);   
    }
    else
    {
       return 0;
    }
}

inline uint64_t mnMapReadAll( LIndex  & index,
                           String<Dna5> & read,
                          Anchors & anchors,
                          MapParm & mapParm,
                          String<uint64_t> & hit  
                         )
{
    getIndexMatchAll(index, read, anchors.set, mapParm);    
    return getAnchorMatchAll(anchors, length(read), mapParm, hit);
}

inline uint64_t mnMapReadFirst( LIndex  & index,
                           String<Dna5> & read,
                          Anchors & anchors,
                          MapParm & mapParm,
                          String<uint64_t> & hit  
                         )
{
    getIndexMatchAll(index, read, anchors.set,  mapParm);    
    return getAnchorMatchFirst(anchors, length(read), mapParm, hit);
}

inline uint64_t mnMapReadList( LIndex  & index,
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
 * this is initCord for single strand (without strand flag) index
inline bool initCord(typename Iterator<String<uint64_t> >::Type & it, 
                     typename Iterator<String<uint64_t> >::Type & hitEnd,
                     unsigned & preCordStart,
                     String<uint64_t> & cord)
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
        appendValue(cord, _DefaultCord.hit2Cord(*(it)));
        ++it;
        preCordStart = length(cord) - 1;   
    }
    return true;
}
*/

/*
 * this is initCord for double strand index(with flag in cord value)
 */
inline bool initCord(typename Iterator<String<uint64_t> >::Type & it, 
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

inline bool endCord( String<uint64_t> & cord,
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
inline bool endCord( String<uint64_t> & cord,
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
inline bool endCord( String<uint64_t> & cord,
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
inline bool nextCord(typename Iterator<String<uint64_t> >::Type & it, 
                     typename Iterator<String<uint64_t> >::Type const & hitEnd,
                     unsigned & preCordStart,
                     String<uint64_t> & cord)
{
//TODO: add maxlen of anchor to the first node in cord;
    if (it >= hitEnd)
        return false;
    while(!_DefaultHit.isBlockEnd(*(it - 1)))
    {
        if(_DefaultCord.getCordY(*(it))>_DefaultCord.getCordY(back(cord)) +  window_delta) 
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
inline bool nextCord(typename Iterator<String<uint64_t> >::Type & it, 
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
        if(_DefaultCord.getCordY(*(it))>_DefaultCord.getCordY(back(cord)) +  window_delta) 
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
inline bool nextCord(typename Iterator<String<uint64_t> >::Type & it, 
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
        if(_DefaultCord.getCordY(*(it))>_DefaultCord.getCordY(back(cord)) +  window_delta) 
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

inline bool extendWindowAll(String<short> &f1, 
                            String<short> & f2, 
                            String<uint64_t> & cord, 
                            float & score, 
                            uint64_t & strand)
{
    uint64_t preCordY = (_DefaultHit.isBlockEnd(cord[length(cord) - 2]))?
    0:_DefaultCord.getCordY(back(cord)) + window_delta;
    unsigned len = length(cord) - 1;
    while (preCordY<= _DefaultCord.getCordY(back(cord)) && 
        previousWindow(f1, f2, cord, score, strand)){}
    for (unsigned k = len; k < ((length(cord) + len) >> 1); k++) 
    {
        std::swap(cord[k], cord[length(cord) - k + len - 1]);
    }
    while (nextWindow(f1, f2, cord, score, strand)){}
    //std::cout << "[]::extendWindow " << k << "\n";
    return true;
}

inline bool isOverlap (uint64_t cord1, uint64_t cord2, int revscomp_const)
{
    uint64_t strand1 = _DefaultCord.getCordStrand(cord1);
    uint64_t strand2 = _DefaultCord.getCordStrand(cord2);
    int64_t x1 = _DefaultCord.getCordX(cord1);
    //int64_t y1 = revscomp_const * strand1 - _nStrand(strand1) * _DefaultCord.getCordY(cord1);
    int64_t y1 = _DefaultCord.getCordY (cord1);
    int64_t x2 = _DefaultCord.getCordX(cord2);
    //int64_t y2 = revscomp_const * strand2 - _nStrand(strand2) * _DefaultCord.getCordY(cord2);
    int64_t y2 = _DefaultCord.getCordY (cord2);
    (void) revscomp_const;
    return std::abs(x1 - x2) < window_size && std::abs(y1 - y2) < window_size && (!(strand1 ^ strand2));
}
/**
 * cord1 is predecessor of cord2 and they are not overlapped
 */
inline bool isPreGap (uint64_t cord1, uint64_t cord2, int revscomp_const)
{
    //uint64_t strand1 = _DefaultCord.getCordStrand(cord1);
    //uint64_t strand2 = _DefaultCord.getCordStrand(cord2);
    int64_t x1 = _DefaultCord.getCordX(cord1);
    //int64_t y1 = _DefaultCord.getCordY(cord1);
    int64_t x2 = _DefaultCord.getCordX(cord2);
    //int64_t y2 = _DefaultCord.getCordY(cord2);
    //return (x1 + window_size > x2 || y1 + window_size > y2) && (std::abs(x1 - x2) > window_size || std::abs(y1 - y2) > window_size);
    (void) revscomp_const;
    //std::cout << "[]::isPreGap " << x1 + window_size << " " << x2 << "\n";
    return (x1 + window_size < x2);
}
/**
 * cord1 is successor of cord2 and they are not overlapped
 */
inline bool isSucGap (uint64_t cord1, uint64_t cord2, int revscomp_const)
{
    return isPreGap (cord2, cord1, revscomp_const);
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
inline int extendPatch(StringSet<String<short> > & f1, 
                       StringSet<String<short> > & f2, 
                      String<uint64_t> & cords,
                      int k,
                      uint64_t cord1,
                      uint64_t cord2,
                      int revscomp_const
                        )
{
    unsigned window_threshold = 30;
    //std::cout << "[]::extendPatch::coord cord1y " << _DefaultCord.getCordY(cord1) << " cord2y " << _DefaultCord.getCordY(cord2) << " " << _DefaultCord.getCordStrand(cord1 ^ cord2) << "\n";
    if (isOverlap(cord1, cord2, revscomp_const))
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
    uint64_t genomeId1 = _getSA_i1(_DefaultCord.getCordX(pcord));
    uint64_t genomeId2 = _getSA_i1(_DefaultCord.getCordX(scord));
    int len = 0;
    uint64_t cord = pcord;
    String<uint64_t> tmp;
    //std::cout << "[]::extendPatch pcord " << _DefaultCord.getCordY(cord) << " " << _DefaultCord.getCordY(scord) << "\n";
    while (isPreGap(cord, scord, revscomp_const))
    {
        //std::cout << "[]::extendWindow::n cord " << _DefaultCord.getCordY(cord) << "\n";
        cord = nextWindow (f1[strand1], f2[genomeId1], cord, window_threshold);
        if (cord)
        {
            appendValue (tmp, cord);
            //std::cout << "[]::extendWidnow::n append " << _DefaultCord.getCordY(back(tmp)) << " " << _DefaultCord.getCordX(back(tmp)) << " " << strand1 << " " << strand2 << "\n";
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
        insert(cords, k, tmp);
        clear(tmp);
    }
    
    cord = scord;
    while (isSucGap(cord, nw, revscomp_const))
    {
        //std::cout << "[]::extendWindow::p cord " << _DefaultCord.getCordY(cord) << "\n";
        cord = previousWindow(f1[strand2], f2[genomeId2], cord, window_threshold);
        if (cord)
        {
            //TODO do another round previousWindow if cord = 0.
            appendValue (tmp, cord);
            //std::cout << "[]::extendWidnow::p append " << _DefaultCord.getCordY(cord) << " " << _DefaultCord.getCordX(cord) << " " << _DefaultCord.getCordY (nw) << " " << _DefaultCord.getCordX(nw) << "\n";
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
        insert(cords, k + len, tmp);
        len += length(tmp);
    }
    return len;
}

/*
inline bool pathAll(String<Dna5> & read, 
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
            extendWindowAll(f1, f2[_getSA_i1(_DefaultCord.getCordX(back(cords)))], cords, score);
        }
        while (nextCord(it, hitEnd, preBlockPtr, cords));
        return endCord(cords, preBlockPtr);   
    }
    return false;
}

inline bool pathAll(String<Dna5> & read, 
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
            extendWindowAll(f1, f2[_getSA_i1(_DefaultCord.getCordX(back(cords)))], cords, score);
        }
        while (nextCord(it, hitEnd, preBlockPtr, cords, cordLenThr, strand, score));
        return endCord(cords, preBlockPtr, cordLenThr, strand, score);   
    }
    return false;
}
*/

/*
 * path functoin for single strand
inline bool pathAll(String<Dna5> & read, 
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
            extendWindowAll(f1, f2[_getSA_i1(_DefaultCord.getCordX(back(cords)))], cords, score);
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
inline bool path_dst(
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
            extendWindowAll(f1[strand], f2[_getSA_i1(_DefaultCord.getCordX(back(cords)))], cords, score, strand);
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

/*
 * This mapping function use the hash value of double strand 
 * that allow to stream the sequence of double strand for only once
 * pass feature instead of genomes to reduce memory footprint
 *
template <typename TDna, typename TSpec>
int rawMap_dst2_MF(typename LIndex   & index,
            StringSet<String<short> > & f2,
            typename StringSet<String<Dna5> >      & reads,
            MapParm & mapParm,
            StringSet<String<uint64_t> > & cords,
            unsigned & threads,
            StringSet<String<TDna> > & seqs
            )
{
  
    typedef typename String<Dna5> Seq;
    //double time=sysTime();
    float senThr = mapParm.senThr / window_size;
    float cordThr = mapParm.cordThr / window_size;
    MapParm complexParm = mapParm;
    complexParm.alpha = complexParm.alpha2;
    complexParm.listN = complexParm.listN2;
    //double time2 = sysTime();
#pragma omp parallel
{
    unsigned size2 = length(reads) / threads;
    unsigned ChunkSize = size2;
    Seq comStr;
    //Anchors anchors(Const_::_LLTMax, AnchorBase::size);
    Anchors anchors;
    typename String<uint64_t>  crhit;
    StringSet<String<uint64_t> >  cordsTmp;
    StringSet< String<short> > f1;
    unsigned thd_id =  omp_get_thread_num();
    if (thd_id < length(reads) - size2 * threads)
    {
        ChunkSize = size2 + 1;
    }
    resize(cordsTmp, ChunkSize);
    resize(f1, 2);
    unsigned c = 0;
    
    String<uint64_t>  g_hs;
    String<uint64_t>  g_anchor;
    resize (g_hs, 1ULL << 20);
    resize (g_anchor, 1ULL<<20);

    #pragma omp for
    for (unsigned j = 0; j < length(reads); j++)
    {
        if (length(reads[j]) >= mapParm.minReadLen)
        {
            float cordLenThr = length(reads[j]) * cordThr;
            _compltRvseStr(reads[j], comStr);
            createFeatures(begin(reads[j]), end(reads[j]), f1[0]);
            createFeatures(begin(comStr), end(comStr), f1[1]);
            anchors.init(1);
            clear(crhit);
            mnMapReadList<TDna, TSpec>(index, reads[j], anchors, mapParm, crhit);
            path_dst(begin(crhit), end(crhit), f1, f2, cordsTmp[c], cordLenThr);
            if (_DefaultCord.getMaxLen(cordsTmp[c]) < length(reads[j]) * senThr)// && 
               //_DefaultCord.getMaxLen(cordsTmp[c]) > 0)
            {
                clear(cordsTmp[c]);
                anchors.init(1);
                clear(crhit);
                mnMapReadList<TDna, TSpec>(index, reads[j], anchors, complexParm, crhit);
                path_dst(begin(crhit), end(crhit), f1, f2, cordsTmp[c], cordLenThr);
            }   
            //mapGaps ()
            mapGaps(seqs, reads[j], cordsTmp[c], g_hs, g_anchor, 500, 192);
        }
        c += 1;
    } 
    #pragma omp for ordered
    for (unsigned j = 0; j < threads; j++)
        #pragma omp ordered
        {
            append(cords, cordsTmp);
        }
}
    //std::cerr << "    End raw mapping. Time[s]: " << sysTime() - time << std::flush << std::endl;
    return 0;
}
*/
//End all mapper module
//============================================


