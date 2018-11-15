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

#include "base.h"

using namespace seqan;

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
    typedef uint64_t CordType;
    typedef uint64_t CellType;
    typedef unsigned Size;
    
    Bit bit;
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
    
    CordBase():
        bit(20),    
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
    
}_DefaultCordBase;

struct Cord
{
    typedef typename CordBase::CordType CordType;
    typedef typename CordBase::CellType CellType;
    typedef typename CordBase::Flag Flag;
    
    typedef String<CordType> CordString;
    typedef StringSet<CordString> CordSet;
    
    CordType getCordX(CordType const &, typename CordBase::Bit const &, typename CordBase::Mask const &) const;
    CordType getCordY(CordType const &, typename CordBase::Mask const &) const;
    CordType createCord(CordType const &, CordType const &, CordType const &, typename CordBase::Bit const &, typename CordBase::Bit const &) const ;
    CordType hit2Cord(PMRes::HitType const &, typename CordBase::Bit const &, typename CordBase::Mask const &, typename CordBase::Mask const &) const;
    CordType hit2Cord_dstr(PMRes::HitType const &, typename CordBase::Bit const &, typename CordBase::Mask const &, typename CordBase::Mask const &) const;
    CellType cord2Cell(CordType const &, typename CordBase::Bit const &) const;
    CordType cell2Cord(CellType const &, typename CordBase::Bit const &) const;
    void setCordEnd(CordType &, typename CordBase::Flag const &, typename CordBase::Flag const &);
    Flag getCordStrand(CordType const &, CordBase::Bit const &) const;
    Flag AtCordEnd(CordType const &, CordBase::Flag const &)const;
    //void setHead(uint64_t &, uint64_t const &, uint64_t const & = _DefaultCordBase.headFlag);
    void setMaxLen(String<uint64_t> &, uint64_t const &, uint64_t const & = _DefaultCordBase.mask);
    uint64_t getMaxLen(String<uint64_t> const &, uint64_t const & = _DefaultCordBase.mask);
    uint64_t shift(uint64_t & val, int64_t x, int64_t y, unsigned const & = _DefaultCordBase.bit); //add x and y

    
    bool print (CordString const &, std::ostream & = std::cout, CordBase const & = _DefaultCordBase) const;
    bool print (CordSet const &, std::ostream & = std::cout, CordBase const & = _DefaultCordBase) const;
    bool printAlignmentMatrix(CordSet const &,  CordBase const & ) const;
    
}_DefaultCord; 

inline typename Cord::CordType 
Cord::getCordX(typename Cord::CordType const & cord, 
               typename CordBase::Bit const & bit  = _DefaultCordBase.bit,
               typename CordBase::Mask const & mask = _DefaultCordBase.maskx) const
{
    return (cord >> bit) & mask; 
}

inline typename Cord::CordType 
Cord::getCordY(typename Cord::CordType const & cord, 
               typename CordBase::Mask const & mask = _DefaultCordBase.mask) const 
{
    return cord & mask;
}

inline typename Cord::CordType 
Cord::createCord(typename Cord::CordType const & x, 
                 typename Cord::CordType const & y, 
                 typename Cord::CordType const & strand,
                 typename CordBase::Bit const & bit = _DefaultCordBase.bit, 
                 typename CordBase::Bit const & bit2 = _DefaultCordBase.flag_bit) const
{
    return (x << bit) + y + (strand << bit2);
}

inline typename Cord::CordType 
Cord::hit2Cord(typename PMRes::HitType const & hit, 
               typename CordBase::Bit const & bit = _DefaultCordBase.bit, 
               typename CordBase::Mask const & mask = _DefaultCordBase.mask,
               typename CordBase::Mask const & mask2 = _DefaultCordBase.valueMask
              ) const
{
    return (hit + ((hit & mask) << bit)) & mask2;
}

inline typename Cord::CordType 
Cord::hit2Cord_dstr(typename PMRes::HitType const & hit, 
               typename CordBase::Bit const & bit = _DefaultCordBase.bit, 
               typename CordBase::Mask const & mask = _DefaultCordBase.mask,
               typename CordBase::Mask const & mask2 = _DefaultCordBase.valueMask_dstr
              ) const
{
    return (hit + ((hit & mask) << bit)) & mask2;
}

inline typename Cord::CellType 
Cord::cord2Cell(typename Cord::CordType const & cord, 
                typename CordBase::Bit const & bit = _DefaultCordBase.cell_bit) const
{
    return cord >> bit;
}

inline typename Cord::CordType 
Cord::cell2Cord(typename Cord::CellType const & cell, 
                typename CordBase::Bit const & bit = _DefaultCordBase.cell_bit) const
{
    return cell << bit;
}

inline void Cord::setCordEnd(typename Cord::CordType & cord,
            typename CordBase::Flag const & strand = _DefaultCordBase.flag_strand,
            typename CordBase::Flag const & end = _DefaultCordBase.flag_end)
{
    cord |= strand | end;
}

inline typename CordBase::Flag 
Cord::getCordStrand(typename Cord::CordType const & cord,
            typename CordBase::Bit const & strand = _DefaultCordBase.flag_bit) const
{
    return (cord >> strand) & 1ULL;
}

inline typename CordBase::Flag 
Cord::AtCordEnd(typename Cord::CordType const & cord,
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

inline uint64_t Cord::shift(uint64_t & val, int64_t x, int64_t y, unsigned const & bit) //add x and y
{
    return uint64_t((int64_t)val + (x << bit) + y);
}

inline bool Cord::print(typename Cord::CordString const & cords, std::ostream & of, CordBase const & cordBase) const
{
    of << "length of cords " << length(cords) << std::endl;
    for (unsigned j = 1; j < length(cords); j++)
               of << getCordY(cords[j], cordBase.mask) << " " 
                  << _getSA_i1(getCordX(cords[j], cordBase.bit)) << " "
                  << _getSA_i2(getCordX(cords[j], cordBase.bit))  << std::endl;
    of << std::endl;
    return true;
}



inline bool Cord::print(typename Cord::CordSet const & cords, std::ostream & of, CordBase const & cordBase) const
{
    of << "Cord::print() " << std::endl;
    unsigned j = 0;
    for (auto && k : cords)
    {
        of << j++ << std::endl;
        print(k, of, cordBase);
    }
    return true;
}


    

//inline bool Cord::printAlignmentMatrix(typename Cord::CordSet const & cords, unsigned const & readLen, CordBase const & cordBase = _DefaultCordBase) const
//{
//    unsigned it = 0;
//    for (unsigned y = 0; y < readLen / 16; y++)
//    {
//        for (; it < length(cords); it++)
//        {
//            for (unsigned k = 0; k < 100; k++)
//            {
//                for (unsigned j =0; j < 100; j++)
//                    if (cords)
//                    std::cout << "*";
//                std::cout << std::endl;
//            }
//        }
//    }
//    return true;
//}
//


//======HIndex getIndexMatch()

//Note: length of read should be < 1MB;




static const float band_width = 0.25;
//static const unsigned cmask = ((uint64_t)1<<32) - 1;
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
//static const uint64_t scriptMask = (1-3*scriptWindow) -1;
static const int scriptCount[5] = {1, 1<<scriptWindow, 1 <<(scriptWindow * 2), 0, 0};
static const int scriptMask = (1 << scriptWindow) - 1;
static const int scriptMask2 = scriptMask << scriptWindow;
static const int scriptMask3 = scriptMask2 << scriptWindow;

static const uint64_t hmask = (1ULL << 20) - 1;
/**
 * ATTENTION TODO parameter needs tuning: will affect speed, gap extension, clip
 */
static const unsigned windowThreshold = 36; // 36;


template <typename TIter>
inline uint64_t getScript(TIter const & it)
{
    uint64_t script=0;
     
    for (unsigned k=0; k<(1<<scriptWindow); k++)
    {
        script+= scriptCount[ordValue(*(it+k))];
    }
    return script;
}

inline int _scriptDist(int const & s1, int const & s2)
{
    //return std::abs((s1 & scriptMask)- (s2 & scriptMask)) + std::abs(((s1 & scriptMask2) - (s2 & scriptMask2)) >> scriptWindow) + std::abs((s1>>scriptWindow*2) - (s2>>scriptWindow*2));
    int res = std::abs((s1 & scriptMask)- (s2 & scriptMask)) + std::abs(((s1 >> scriptWindow) & scriptMask) - ((s2 >> scriptWindow) & scriptMask)) + std::abs((s1>>scriptWindow2) - (s2>>scriptWindow2));
    //printf("[sc] %d %d %d %d %d\n", res, (s1 & scriptMask), ((s1 >> scriptWindow)&scriptMask), (s1>>scriptWindow2), s1);
    //printf ("[sc] %d %d\n", s1, s2);
    return res;
}

/*
 * calculate features of reverse complement seuqnece
 * !Note: the function is not tested.
inline void _reverseComplementFeature(String<int> & f1, String<int> & f2)
{
    resize(f2, length(f1));
    for (uint64_t k = 0; k < length(f1); k++)
    {
        f2[k] = (uint64_t)window_size - (f1[length(f1) - k - 1] & scriptMask) - 
        (f1[length(f1) - 1 - k] >> scriptWindow2 & scriptMask) - (f1[length(f1) - 
        1 - k] >> scriptWindow3 & scriptMask) + (f1[length(f1) - 1 - k] & 
        scriptMask2 << scriptWindow) + (f1[length(f1) - 1 - k] & scriptMask3 >> 
        scriptWindow);
    }
}
*/

template<typename TIter> 
inline void createFeatures(TIter const & itBegin, TIter const & itEnd, String<short> & f)
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
        //printf("[createFeatures] %d \n", f[next]);
        next++;
    }
}


template <typename T> //unsigned or uint64_t
inline T parallelParm_Static(T range, unsigned threads, unsigned & thd_id, 
                                T & thd_begin, T & thd_end)
{
    T ChunkSize = range / threads;
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
//    printf("[parallelParm_Static] %d, %d\n",thd_id, ChunkSize);
    return ChunkSize;
}

/*
 * parallel
 */
template<typename TIter> 
inline void createFeatures(TIter const & itBegin, TIter const & itEnd, String<short> & f, unsigned threads)
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
    //printf("[createfeature] %d %d %d\n", thd_id, thd_begin, thd_end);
    uint64_t next = thd_begin;
    thd_begin *= scriptStep;
    thd_end *= scriptStep;
    f[next] = 0;
    for (unsigned k = thd_begin; k < thd_begin + window; k++)
    {
        f[next] += scriptCount[ordValue(*(itBegin + k))];
    }
    //printf("%d %d %d\n", next, thd_id, f[next]);
    next++;
    for (unsigned k = thd_begin + scriptStep; k < thd_end ; k+=scriptStep) 
    {
        f[next] = f[next - 1];
        for (unsigned j = k - scriptStep; j < k; j++)
        {
            f[next] += scriptCount[ordValue(*(itBegin + j + window))] - scriptCount[ordValue(*(itBegin + j))];
        }
        //printf("[createFeatures] %d \n", f[next]);
        next++;
    }
    //printf("[createfeature1] %d %d %d\n", thd_id, next);
}
    
}

/*
 * parallel
 */
template<typename TDna> 
inline void createFeatures(StringSet<String<TDna> > & seq, StringSet<String<short> > & f, unsigned threads)
{
    resize(f, length(seq));
    for (unsigned k = 0; k < length(seq); k++)
        createFeatures(begin(seq[k]), end(seq[k]), f[k], threads);
}

template<typename TDna> 
inline void createFeatures(StringSet<String<TDna> > & seq, StringSet<String<short> > & f)
{
    resize(f, length(seq));
    for (unsigned k = 0; k < length(seq); k++)
        createFeatures(begin(seq[k]), end(seq[k]), f[k]);
}

template<typename TIter>
inline unsigned _windowDist(TIter const & it1, TIter const & it2)
{
    return _scriptDist(*it1, *it2)+ _scriptDist(*(it1 + 2), *(it2 + 2)) + _scriptDist(*(it1+4),*(it2+4)) + _scriptDist(*(it1 + 6), *(it2 + 6)) + _scriptDist(*(it1+8), *(it2+8)) + _scriptDist(*(it1 + 10), *(it2 + 10));
}

inline bool nextCord(typename PMRes::HitString & hit, unsigned & currentIt, String<uint64_t> & cord)
{
    //std::cout << "nextCord\n";
    //std::cout << "length(hit) " << length(hit) << std::endl;
    uint64_t cordLY = _DefaultCord.getCordY(back(cord));
    while (++currentIt < length(hit)) 
    //while (--currentIt > 0) 
    {
        uint64_t tmpCord = _DefaultCord.hit2Cord(hit[currentIt]);
        //if(_DefaultCord.getCordY(tmpCord) < cordLY)
        if(_DefaultCord.getCordY(tmpCord) > cordLY + window_delta)
        {
            appendValue(cord, tmpCord);
            return true;
        }
    }
    return false;
}


inline bool initCord(typename PMRes::HitString & hit, unsigned & currentIt, String<uint64_t> & cord)
{
    currentIt = 0;
    //currentIt = length(hit);
    if (empty(hit))
        return false;
    else
        appendValue(cord, _DefaultCord.hit2Cord(hit[0]));
        //appendValue(cord, _DefaultCord.hit2Cord(back(hit)));
    //std::cerr << "init" << (cord[0] >> 20) << " " << (cord[0] & 0xfffff) << std::endl;
    return true;
}

inline bool previousWindow(String<short> & f1, 
                           String<short> & f2, 
                           typename Cord::CordString & cord, 
                           uint64_t & strand)
{
    typedef typename Cord::CordType CordType;
    CordType genomeId = _getSA_i1(_DefaultCord.getCordX(back(cord)));
    CordType x_suf = _DefaultCord.cord2Cell(_getSA_i2(_DefaultCord.getCordX(back(cord))));
    CordType y_suf = _DefaultCord.cord2Cell(_DefaultCord.getCordY(back(cord)));
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
                           typename Cord::CordString & cord, 
                           float & score, 
                           uint64_t & strand)
{
    typedef typename Cord::CordType CordType;
    CordType genomeId = _getSA_i1(_DefaultCord.getCordX(back(cord)));
    CordType x_suf = _DefaultCord.cord2Cell(_getSA_i2(_DefaultCord.getCordX(back(cord))));
    CordType y_suf = _DefaultCord.cord2Cell(_DefaultCord.getCordY(back(cord)));
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

    //printf("[debug]::previousWindow %d\n", min);
    score += min;
    return true;
}

inline uint64_t previousWindow(String<short> & f1, 
                               String<short> & f2, 
                               uint64_t cordx,
                               uint64_t cordy,
                               uint64_t strand,
                               unsigned window_threshold = windowThreshold
                              )
{
    typedef typename Cord::CordType CordType;
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
        //std::cout << "[]::previous score " << min << " " <<_DefaultCord.cell2Cord(x_suf - med) << "\n";
            return _DefaultCord.createCord(_createSANode(genomeId, _DefaultCord.cell2Cord(x_suf - med)),  _DefaultCord.cell2Cord(x_suf - x_min - med + y), strand);
        }
        else
        {
        //std::cout << "[]::previous score " << min << " " << _DefaultCord.cell2Cord(x_min) << "\n";
            return _DefaultCord.createCord(_createSANode(genomeId, _DefaultCord.cell2Cord(x_min)), _DefaultCord.cell2Cord(y), strand);
        } 
    }
    return 0;
}

inline uint64_t previousWindow(String<short> & f1, String<short> & f2, uint64_t cord, unsigned window_threshold)
{
    return previousWindow(f1, f2, _DefaultCord.getCordX(cord), _DefaultCord.getCordY(cord), _DefaultCord.getCordStrand(cord), window_threshold);
}

inline bool nextWindow(String<short> &f1, 
                       String<short> & f2, 
                       typename Cord::CordString & cord, 
                       uint64_t & strand)
{
    typedef typename Cord::CordType CordType;
    CordType genomeId = _getSA_i1(_DefaultCord.getCordX(back(cord)));
    CordType x_pre = _DefaultCord.cord2Cell(_getSA_i2(_DefaultCord.getCordX(back(cord))));
    CordType y_pre = _DefaultCord.cord2Cell(_DefaultCord.getCordY(back(cord)));
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
                       typename Cord::CordString & cord, 
                       float & score, 
                       uint64_t & strand)
{
    typedef typename Cord::CordType CordType;
    CordType genomeId = _getSA_i1(_DefaultCord.getCordX(back(cord)));
    CordType x_pre = _DefaultCord.cord2Cell(_getSA_i2(_DefaultCord.getCordX(back(cord))));
    CordType y_pre = _DefaultCord.cord2Cell(_DefaultCord.getCordY(back(cord)));
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
            //std::cout << "[]::nextWindow min score " << min << " " << _DefaultCord.cell2Cord(x_pre + med) << "\n";
            appendValue(cord, _DefaultCord.createCord(_createSANode(genomeId, _DefaultCord.cell2Cord(x_pre + med)),  _DefaultCord.cell2Cord(x_pre + med - x_min + y), strand));
        }
        else
        {
            //std::cout << "[]::nextWindow min score " << min << " " << _DefaultCord.cell2Cord(x_min) << "\n";
            appendValue(cord, _DefaultCord.createCord(_createSANode(genomeId, _DefaultCord.cell2Cord(x_min)), _DefaultCord.cell2Cord(y), strand));
        }
    score += min;
        //printf("[debug]::nextWindow %d\n", min);

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
    typedef typename Cord::CordType CordType;
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
    //std::cout << "[]::nextWindow score " << min << "\n";
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

inline bool extendWindow(String<short> &f1, String<short> & f2, typename Cord::CordString & cord, uint64_t strand)
{
    Cord::CordType preCord = (length(cord)==1)?0:back(cord);
    unsigned len = length(cord) - 1;
    while (preCord <= back(cord) && previousWindow(f1, f2, cord, strand)){}
    for (unsigned k = len; k < ((length(cord) + len) >> 1); k++) 
    {
        std::swap(cord[k], cord[length(cord) - k + len - 1]);
    }
    while (nextWindow(f1, f2, cord, strand)){}
    return true;
}
/*
inline  bool extendWindow(String<short> & f1, String<short> & f2, String<uint64_t> & cords, uint64_t strand, int k)
{
    //uint64_t preCord = (length(cords) == 1) ? 0 : back(cords);
    ///unsigned len = length(cords) - 1;
    //while (preCord <= back(cords) && previousWindow(f1, f2, cord, strand));
    while (nextWindow(f1, f2, ))
}
*/

inline bool path(String<Dna5> & read, typename PMRes::HitString hit, StringSet<String<short> > & f2, String<uint64_t> & cords)
{
    String<short> f1;
    unsigned currentIt = 0;
    if(!initCord(hit, currentIt, cords))
        return false;
//    std::cerr << "length " << length(read) << std::endl;
    createFeatures(begin(read), end(read), f1);
    unsigned genomeId = _getSA_i1(_DefaultCord.getCordX(cords[0]));
   // uint64_t pre = _getSA_i1(_DefaultCord.getCordX(cords[0]));
   // for (unsigned k = 0; k<length(cords);k++)
   // {
   // 
   // if (pre != _getSA_i1(_DefaultCord.getCordX(cords[k])))
   //     std::cerr << "id error " << pre << " " << _getSA_i1(_DefaultCord.getCordX(cords[k])) << std::endl;
   // pre = _getSA_i1(_DefaultCord.getCordX(cords[k]));
   // }
        

//    std::cerr << "done " << std::endl;
    //std::cout << "[debug] " << _DefaultCord.getCordY(back(cords)) << "\n";
    while (_DefaultCord.getCordY(back(cords)) < length(read) - window_size)
    {
        extendWindow(f1, f2[genomeId], cords, _DefaultCord.getCordStrand(back(cords)));
        //std::cerr << "nextCord" << std::endl;
        if(!nextCord(hit, currentIt, cords))
                return false;
    }
    return true;
}

void path(typename PMRes::HitSet & hits, StringSet<String<Dna5> > & reads, StringSet<String<Dna5> > & genomes, StringSet<String<uint64_t> > & cords)
{
    StringSet<String<short> > f2;
    createFeatures(genomes, f2);
    std::cerr << "raw mapping... " << std::endl;
    for (unsigned k = 0; k < length(reads); k++)
    {
        path(reads[k], hits[k], f2, cords[k]);
    }
//    _DefaultCord.print(cords);
//    _DefaultCord.printAlignmentMatrix(cords);
}

void checkPath(typename Cord::CordSet const & cords, StringSet<String<Dna5> > const & reads)
{
    unsigned count = 0;
    Iterator<typename Cord::CordSet const>::Type it = begin(cords);
    for (auto && read : reads)
    {
        //std::cout << _DefaultCord.getCordY(back(*it)) << " len " << length(reads) << std::endl;
        if(empty(*it))
            count++;
        else
            if (_DefaultCord.getCordY(back(*it)) + window_size * 2 < length(read))
            {
                count++;
            }
        it++;
    }
    //std::cerr << "checkPath " << (float) count / length(reads) << std::endl;
}

/*==================================================
*   this part is for different types of mapping
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
    //val = strand?((1ULL << bit)|val):val;
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
        //printf("[printhit] %d %d %d % " PRIu64 " %" PRIu64 " %" PRIu64 " %" PRIu64 " %" PRIu64 " %d %d\n", j, id1, id2, hit[k], _getSA_i1(_DefaultCord.getCordX(hit[k])), _getSA_i2(_DefaultCord.getCordX(hit[k])), _DefaultCord.getCordY(hit[k]), end, _DefaultCord.getMaxLen(hit), len);
    }
}

void _printHit(String<uint64_t>  & hit)
{
    for (unsigned k = 0; k < length(hit); k++)
    {
    std::cout << "hit " << _getSA_i1(_DefaultCord.getCordX(_DefaultCord.hit2Cord(hit[k]))) << " " << _getSA_i2(_DefaultCord.getCordX(_DefaultCord.hit2Cord(hit[k]))) << " " << _DefaultCord.getCordY(hit[k]) << "\n";
        if (_DefaultHit.isBlockEnd(hit[k]))
            std::cout << "end\n";
    }
}

//===!Note:Need to put this parameterin the mapper threshold
/*
template <typename TDna, typename TSpec>
inline unsigned getIndexMatchAll(typename PMCore<TDna, TSpec>::Index & index,
                              typename PMRecord<TDna>::RecSeq const & read,
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
inline unsigned getIndexMatchAll(typename PMCore<TDna, TSpec>::Index & index,
                              typename PMRecord<TDna>::RecSeq const & read,
                              uint64_t* const set,
                              unsigned & len,
                              MapParm & mapParm
                             )
{   
    typedef typename PMCore<TDna, TSpec>::Index TIndex;
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

/*
 * search double strand in one round
 */
template <typename TDna, typename TSpec>
inline unsigned getIndexMatchAll(typename PMCore<TDna, TSpec>::Index & index,
                              typename PMRecord<TDna>::RecSeq const & read,
                              //uint64_t* const set,
                              String<uint64_t> & set,
                              MapParm & mapParm)
{   
    typedef typename PMCore<TDna, TSpec>::Index TIndex;
    typedef typename TIndex::TShape PShape;
    
    //printf("[debug]\n");
    //unsigned block = (mapParm.blockSize < length(read))?mapParm.blockSize:length(read);
    //unsigned dt = block * (mapParm.alpha / (1 - mapParm.alpha));
    unsigned dt = 0;
    PShape shape;
    uint64_t xpre = 0;
    hashInit(shape, begin(read));
    for (unsigned k = 0; k < length(read); k++)
    {
        //printf("[debug]\n");
        //unsigned count = 0;
        //hashNext(shape, begin(read) + k);
        hashNexth(shape, begin(read) + k);

        uint64_t pre = ~0;
        //uint64_t pos = getXYDir(index, shape.XValue, shape.YValue);
        if (++dt == mapParm.alpha)
        {
            if(hashNextX(shape, begin(read) + k) ^ xpre)
            {
                xpre = shape.XValue;
                uint64_t pos = getXDir(index, shape.XValue, shape.YValue);
//!Note: This contition is different from single strand which will slightly 
//changes the senstivity; In the single strand index, if the size of the block having the same
//is > mapParm.delta, then the block of the ysa will not be used. 
//While in double strand index, the length of ysa of fixed value includes both + 
//and - strands.
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
                        ++pos;
                    }
                }
            }
            dt = 0;
        }
    }
    return 0;
}

inline uint64_t getAnchorMatchAll(Anchors & anchors, unsigned const & readLen, MapParm & mapParm, PMRes::HitString & hit)
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

inline uint64_t getAnchorMatchFirst(Anchors & anchors, unsigned const & readLen, MapParm & mapParm, PMRes::HitString & hit)
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

inline uint64_t getAnchorMatchList(Anchors & anchors, unsigned const & readLen, MapParm & mapParm, PMRes::HitString & hit)
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
    for (unsigned k = 1; k <= anchors.length(); k++)
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

template <typename TDna, typename TSpec>
inline uint64_t mnMapReadAll(typename PMCore<TDna, TSpec>::Index  & index,
                          typename PMRecord<TDna>::RecSeq & read,
                          Anchors & anchors,
                          MapParm & mapParm,
                          PMRes::HitString & hit  
                         )
{
    getIndexMatchAll<TDna, TSpec>(index, read, anchors.set, mapParm);    
    return getAnchorMatchAll(anchors, length(read), mapParm, hit);
}

template <typename TDna, typename TSpec>
inline uint64_t mnMapReadFirst(typename PMCore<TDna, TSpec>::Index  & index,
                          typename PMRecord<TDna>::RecSeq & read,
                          Anchors & anchors,
                          MapParm & mapParm,
                          PMRes::HitString & hit  
                         )
{
    getIndexMatchAll<TDna, TSpec>(index, read, anchors.set,  mapParm);    
    return getAnchorMatchFirst(anchors, length(read), mapParm, hit);
}

template <typename TDna, typename TSpec>
inline uint64_t mnMapReadList(typename PMCore<TDna, TSpec>::Index  & index,
                          typename PMRecord<TDna>::RecSeq & read,
                          Anchors & anchors,
                          MapParm & mapParm,
                          PMRes::HitString & hit)
{
    getIndexMatchAll<TDna, TSpec>(index, read, anchors.set, mapParm);    
    //printf("done getinxmatchall\n");
    return getAnchorMatchList(anchors, length(read), mapParm, hit);
}
/*
 * this is initCord for single strand (without strand flag) index
inline bool initCord(typename Iterator<PMRes::HitString>::Type & it, 
                     typename Iterator<PMRes::HitString>::Type & hitEnd,
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
inline bool initCord(typename Iterator<PMRes::HitString>::Type & it, 
                     typename Iterator<PMRes::HitString>::Type & hitEnd,
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
        //printf("[debug]::nextcord new block %f %d %f\n", (float)score/(length(cord) - preCordStart), length(cord) - preCordStart, cordThr);
    }
    else
    {
        erase(cord, preCordStart, length(cord));
    }
   //printf("[debug]::endcord done\n");
    (void)score;
    return true;
}

/*
inline bool nextCord(typename Iterator<PMRes::HitString>::Type & it, 
                     typename Iterator<PMRes::HitString>::Type const & hitEnd,
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
inline bool nextCord(typename Iterator<PMRes::HitString>::Type & it, 
                     typename Iterator<PMRes::HitString>::Type const & hitEnd, 
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
inline bool nextCord(typename Iterator<PMRes::HitString>::Type & it, 
                     typename Iterator<PMRes::HitString>::Type const & hitEnd, 
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
                            typename Cord::CordString & cord, 
                            float & score, 
                            uint64_t & strand)
{
    Cord::CordType preCordY = (_DefaultHit.isBlockEnd(cord[length(cord) - 2]))?
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
                 typename Iterator<PMRes::HitString>::Type hitBegin, 
                 typename Iterator<PMRes::HitString>::Type hitEnd, 
                 StringSet<String<int> > & f2, 
                 String<uint64_t> & cords
                   )
{
    String<int> f1;
    typename Iterator<PMRes::HitString>::Type it = hitBegin;
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
                 typename Iterator<PMRes::HitString>::Type hitBegin, 
                 typename Iterator<PMRes::HitString>::Type hitEnd, 
                 StringSet<String<int> > & f2, 
                 String<uint64_t> & cords,
                 float const & cordThr, 
                 uint64_t const & strand
                   )
{
    String<int> f1;
    typename Iterator<PMRes::HitString>::Type it = hitBegin;
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
                 typename Iterator<PMRes::HitString>::Type hitBegin, 
                 typename Iterator<PMRes::HitString>::Type hitEnd, 
                 StringSet<String<int> > & f2, 
                 String<uint64_t> & cords,
                 float const & cordThr, 
                 uint64_t const & strand
                   )
{
    //std::cerr << "[pathall]\n";
    String<int> f1;
    typename Iterator<PMRes::HitString>::Type it = hitBegin;
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
                 typename Iterator<PMRes::HitString>::Type hitBegin, 
                 typename Iterator<PMRes::HitString>::Type hitEnd, 
                 StringSet<String<short> > & f1,
                 StringSet<String<short> > & f2, 
                 String<uint64_t> & cords,
                 float const & cordLenThr
                )
{
    typename Iterator<PMRes::HitString>::Type it = hitBegin;
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
template <typename TDna, typename TSpec>
int rawMap_dst(typename PMCore<TDna, TSpec>::Index   & index,
            typename PMRecord<TDna>::RecSeqs      & reads,
            typename PMRecord<TDna>::RecSeqs & genomes,
            MapParm & mapParm,
            StringSet<String<uint64_t> > & cords,
            unsigned & threads)
{
  
    typedef typename PMRecord<TDna>::RecSeq Seq;
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
    typename PMRes::HitString crhit;
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
            mnMapReadList<TDna, TSpec>(index, reads[j], anchors, mapParm, crhit);
            path_dst(begin(crhit), end(crhit), f1, f2, cordsTmp[c], cordLenThr);
            //printf("done1\n");
            if (_DefaultCord.getMaxLen(cordsTmp[c]) < length(reads[j]) * senThr)// && 
               //_DefaultCord.getMaxLen(cordsTmp[c]) > 0)
            {
                clear(cordsTmp[c]);
                anchors.init(1);
                clear(crhit);
                mnMapReadList<TDna, TSpec>(index, reads[j], anchors, complexParm, crhit);
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
int rawMap_dst2_MF(typename PMCore<TDna, TSpec>::Index   & index,
            StringSet<String<short> > & f2,
            typename PMRecord<TDna>::RecSeqs      & reads,
            MapParm & mapParm,
            StringSet<String<uint64_t> > & cords,
            unsigned & threads,
            StringSet<String<TDna> > & seqs
            )
{
  
    typedef typename PMRecord<TDna>::RecSeq Seq;
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
    typename PMRes::HitString crhit;
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


