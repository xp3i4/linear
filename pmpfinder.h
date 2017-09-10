// ==========================================================================
//                          Mappeing SMRT reads
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



template <typename TDna, typename TSpec>
inline unsigned getIndexMatch(typename PMCore<TDna, TSpec>::Index  & index,
                              typename PMRecord<TDna>::RecSeq & read,
//                              Anchors & anchor,
        uint64_t* const set,
        unsigned & len,
                              MapParm & mapParm,
                            double & time
                             )
{    
    double t=sysTime();
    unsigned block = (mapParm.blockSize < length(read))?mapParm.blockSize:length(read);
    unsigned dt = block * (mapParm.alpha / (1 - mapParm.alpha));

    hashInit(index.shape, begin(read));
    for (unsigned h=0; h <= length(read) - block; h += dt)
    //for (unsigned h=0; h <= length(read) - block; h += block)
    {
        hashInit(index.shape, begin(read) + h);
        for (unsigned k = h; k < h + block; k++)
        {
            hashNext(index.shape, begin(read) + k);
            uint64_t dn = getDir(index, index.shape);
            uint64_t pre = ~0;
            if(_getBodyCounth(index.dir[dn+1]) - _getBodyCounth(index.dir[dn]) < mapParm.delta)
            {
                for (uint64_t countn = _getBodyCounth(index.dir[dn]); countn < _getBodyCounth(index.dir[dn + 1]); countn++)
                {
                    if (index.sa[countn] - pre > mapParm.kmerStep)
                    {
//==
//appendValue performance improving
//==
//                      anchor.appendValue(index.sa[countn]- k, k);
                        set[len++] = ((index.sa[countn]-k)<<20)+k;
                        pre = index.sa[countn];
                    }
                }
            }
        }
    }
    time += sysTime() - t;
    return time;
}

inline unsigned getAnchorMatch(Anchors & anchors, MapParm & mapParm, PMRes::HitString & hit, double & time )
{
    double t=sysTime();
    uint64_t ak;
    unsigned maxLen = 0, c_b=mapParm.shapeLen, cbb=0, sb=0, end=0, start = 0;
    anchors[0] = anchors[1];
    ak=anchors[0];
    anchors.sort(anchors.begin(), anchors.end());
    for (unsigned k = 1; k <= anchors.length(); k++)
    {
        if (anchors[k] - ak < AnchorBase::AnchorValue)
            cbb++;
        else
        {
            anchors.sortPos2(anchors.begin() + sb, anchors.begin() + k);
            for (uint64_t m = sb+1; m < k; m++)
            {
                if(anchors.deltaPos2(m, m-1) >  mapParm.shapeLen)
                    c_b += mapParm.shapeLen;
                else
                    c_b += anchors.deltaPos2(m, m-1); 
            }
            if (c_b > maxLen)
            {
                maxLen = c_b;
                start = sb;
                end = k;
            }
            sb = k;
            ak = anchors[k];
            cbb = 1;
            c_b = mapParm.shapeLen;
        }

    }

//    std::cerr << start << " " << end << " " << anchors.len<< std::endl;
    for (unsigned k = start; k < end; k++)
        appendValue(hit, anchors[k]);

    time += sysTime() - t;
    return start + (maxLen << 20) ;
}

template <typename TDna, typename TSpec>
inline unsigned mnMapRead(typename PMCore<TDna, TSpec>::Index  & index,
                          typename PMRecord<TDna>::RecSeq & read,
                          Anchors & anchors,
                          MapParm & mapParm,
                          PMRes::HitString & hit,  
                            double & time,
                            double & time2 
                         )
{
    
    getIndexMatch<TDna, TSpec>(index, read, anchors.set, anchors.len, mapParm, time);    
    return getAnchorMatch(anchors, mapParm, hit, time2);
}

template <typename TDna, typename TSpec>
void mnMap(typename PMCore<TDna, TSpec>::Index   & index,
           typename PMRecord<TDna>::RecSeqs      & reads,
           MapParm                      & mapParm,
            typename PMRes::HitSet & hits)
      //      Anchors & anchors)
           //PMResSet             & rs )
{
//    typedef typename Mapper<TDna, TSpec>::Seq Seq;
    typedef typename PMRecord<TDna>::RecSeq Seq;
    typedef typename PMCore<TDna, TSpec>::Anchors Anchors;
    double time=sysTime();

    std::cerr << "Filtering reads\n"; 
    Anchors anchors(Const_::_LLTMax, AnchorBase::size);
    Seq comStr;
    
    double tm=0;
    double tm2=0;
    unsigned count = 0;
    for (unsigned j = 0; j< length(reads); j++)
    {
//        anchors.init();
//==
// anchor fucntion improving
//==
        anchors.len=1;
        if (mnMapRead<TDna, TSpec>(index, reads[j], anchors, mapParm, hits[j], tm, tm2) 
                    < (mapParm.threshold << 20))
        {
            count++;
            _compltRvseStr(reads[j], comStr);
            mnMapRead<TDna, TSpec>(index, comStr, anchors, mapParm, hits[j], tm, tm2);
        }

//========
//rs
//========
        //rs.appendValue(_getSA_i2(anchor[res & mask1].getPos1()), _getSA_i2(anchor[res & mask1].getPos1()) + length);
         
    }
    std::cerr << "Time[s]: " << sysTime() - time << std::endl;
    std::cerr << "tm [s]: " << tm << std::endl;
    std::cerr << "tm2 [s]: " << tm2 << std::endl;
    std::cerr << "count " << count << std::endl;
}


struct CordBase
{
    //Cord(C): coordinates of the vertext of sliding windows
    //C := genomeCord [40] |readCord [20bits]
    //cell [4] is the minimum length the window is allowed to slide in the alignment matrix.
    //genomeCord(gC or xC): = position in the genome >> cell_bit << cell_bit. the last cell_bit bits maybe set to 0
    //gC:= SA node = Seq num i1 [10] | Base num i2 [30]  
    //readCord(rC or yC): ~= position in the read >> cell_bit << cell_bit. the last cell_bit bits maybe set to 0 during process.
    //rC:= Base num [20]
    
    typedef unsigned Bit;
    typedef uint64_t Mask;
    typedef uint64_t CordType;
    typedef uint64_t CellType;
    typedef unsigned Size;
    
    Bit bit;
    Mask mask;
    Bit cell_bit;
    Size cell_size;
    
    CordBase():
        bit(20),    
        mask(0xfffff),
        cell_bit(4),
        cell_size(16)
        {}
    
}_DefaultCordBase;

struct Cord
{
    typedef typename CordBase::CordType CordType;
    typedef typename CordBase::CellType CellType;
    
    typedef String<CordType> CordString;
    typedef StringSet<CordString> CordSet;
    
    CordType getCordX(CordType const &, typename CordBase::Bit const &) const;
    CordType getCordY(CordType const &, typename CordBase::Mask const &) const;
    CordType createCord(CordType const &, CordType const &, typename CordBase::Bit const &) const ;
    CordType hit2Cord(PMRes::HitType const &, typename CordBase::Bit const &, typename CordBase::Mask const &) const;
    CellType cord2Cell(CordType const &, typename CordBase::Bit const &) const;
    CordType cell2Cord(CellType const &, typename CordBase::Bit const &) const;
    
    bool print (CordString const &, CordBase const & ) const;
    bool print (CordSet const &, CordBase const & ) const;
    
}_DefaultCord; 

inline typename Cord::CordType Cord::getCordX(typename Cord::CordType const & cord, typename CordBase::Bit const & bit  = _DefaultCordBase.bit) const
{
    return cord >> bit; 
}

inline typename Cord::CordType Cord::getCordY(typename Cord::CordType const & cord, typename CordBase::Mask const & mask = _DefaultCordBase.mask) const 
{
    return cord & mask;
}

inline typename Cord::CordType Cord::createCord(typename Cord::CordType const & x, Cord::CordType const & y, typename CordBase::Bit const & bit = _DefaultCordBase.bit) const
{
    return (x << bit) + y;
}

inline typename Cord::CordType Cord::hit2Cord(typename PMRes::HitType const & hit, typename CordBase::Bit const & bit = _DefaultCordBase.bit, typename CordBase::Mask const & mask = _DefaultCordBase.mask) const
{
    return hit + ((hit & mask) << bit);
}

inline typename Cord::CellType Cord::cord2Cell(typename Cord::CordType const & cord, typename CordBase::Bit const & bit = _DefaultCordBase.cell_bit) const
{
    return cord >> bit;
}

inline typename Cord::CordType Cord::cell2Cord(typename Cord::CellType const & cell, typename CordBase::Bit const & bit = _DefaultCordBase.cell_bit) const
{
    return cell << bit;
}

inline bool Cord::print(typename Cord::CordString const & cords, CordBase const & cordBase = _DefaultCordBase) const
{
    for (auto && j : cords)
       std::cout << getCordX(j, cordBase.bit) << " " << getCordY(j, cordBase.mask); 
    std::cout << std::endl;
    return true;
}

inline bool Cord::print(typename Cord::CordSet const & cords, CordBase const & cordBase = _DefaultCordBase) const
{
    for (auto && k : cords)
       print(k, cordBase);
    return true;
}

static const float band_width = 0.25;
//static const unsigned cmask = ((uint64_t)1<<32) - 1;
static const unsigned cmask = ((uint64_t)1<<20) - 1;
static const unsigned cell_size = 16;
static const unsigned cell_num = 12;
static const unsigned window_size = cell_size * cell_num; //16*12
static const unsigned sup = cell_num;
static const unsigned med =ceil((1 - band_width) * cell_num);
static const unsigned inf = ceil((1 - 2 * band_width) * cell_num);

static const unsigned initx = 5; 
static const unsigned inity = 5;

//======================================================
static const unsigned scriptStep=16;
static const unsigned scriptWindow=6; //2^6
//static const uint64_t scriptMask = (1-3*scriptWindow) -1;
static const int scriptCount[5] = {1, 1<<scriptWindow, 1 <<(scriptWindow * 2), 0, 0};
static const int scriptMask = (1 << scriptWindow) - 1;
static const int scriptMask2 = scriptMask << scriptWindow;

static const uint64_t hmask = (1ULL << 20) - 1;

static const unsigned windowThreshold = 30;


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
    return std::abs((s1 & scriptMask)- (s2 & scriptMask)) + std::abs(((s1 & scriptMask2) - (s2 & scriptMask2)) >> scriptWindow) + std::abs((s1>>scriptWindow*2) - (s2>>scriptWindow*2));
}


template<typename TIter> 
inline void createFeatures(TIter const & itBegin, TIter const & itEnd, String<int> & f)
{
    unsigned next = 1;
    unsigned window = 1 << scriptWindow;
    resize (f, (itEnd - itBegin -window) / scriptStep);
    f[0] = 0;
    
    for (unsigned k = 0; k < window; k++)
    {
        f[0] += scriptCount[ordValue(*(itBegin + k))];
    }
    for (unsigned k = scriptStep; k < itEnd - itBegin - window ; k+=scriptStep) 
    {
        f[next] = f[next - 1];
        for (unsigned j = k - scriptStep; j < k; j++)
            f[next] += scriptCount[ordValue(*(itBegin + j + window))] - scriptCount[ordValue(*(itBegin + j))];
        next++;
    }
    

}

template<typename TDna> 
inline void createFeatures(StringSet<String<TDna> > & seq, StringSet<String<int> > & f)
{
    resize(f, length(seq));
    for (unsigned k = 0; k < length(seq); k++)
    createFeatures(begin(seq[k]), end(seq[k]), f[k]);
}

template<typename TIter>
inline unsigned _windowDist(TIter const & it1, TIter const & it2)
{
    return _scriptDist(*it1, *it2)+ _scriptDist(*(it1+4),*(it2+4)) + _scriptDist(*(it1+8), *(it2+8));
}

//inline bool previousWindow(String<int> & f1, String<int> & f2, typename Cord::Cord & cord)
//{
//}

/*
inline uint64_t nextWindow(String<int> &f1, String<int> & f2, String<uint64_t> & cord)
{
    uint64_t x_pre = *(end(cord) - 1) >> 24;
    uint64_t y_pre = *(end(cord) - 1) & cmask;
    uint64_t y = y_pre + med;
    unsigned min = ~0;
    unsigned x_min = 0;
    for (uint64_t x = x_pre + inf; x < x_pre + sup; x += 1) 
    {
        unsigned tmp = _windowDist(begin(f1) + y, begin(f2) + x);
        if (tmp < min)
        {
            min = tmp;
            x_min = x;
        }
    }
    if (min > windowThreshold)
        return 1;
    else 
        if ( x_min - x_pre > med)
            appendValue(cord, ((uint64_t)(x_pre + med) << 24) + x_pre + med - x_min + y);
        else
            appendValue(cord, ((uint64_t)x_min << 24) + y);
    return 0;
}
*/

inline bool nextWindow(String<int> &f1, String<int> & f2, typename Cord::CordString & cord)
{
    typedef typename Cord::CordType CordType;
    CordType x_pre = _DefaultCord.cord2Cell(_DefaultCord.getCordX(back(cord)));
    CordType y_pre = _DefaultCord.cord2Cell(_DefaultCord.getCordY(back(cord)));
    CordType y = y_pre + med;
    
    unsigned min = ~0;
    unsigned x_min = 0;
    //std::cout << "nextWindow() x_pre " << x_pre << std::endl;
    for (CordType x = x_pre + inf; x < x_pre + sup; x += 1) 
    {
        unsigned tmp = _windowDist(begin(f1) + y, begin(f2) + x);
        if (tmp < min)
        {
            min = tmp;
            x_min = x;
        }
    }
    std::cout << "nextWindow " << min << " windowThreshold " << windowThreshold << std::endl;
    if (min > windowThreshold)
        return false;
    else 
        if ( x_min - x_pre > med)
            //appendValue(cord, ((uint64_t)(x_pre + med) << 24) + x_pre + med - x_min + y);
            appendValue(cord, _DefaultCord.createCord(_DefaultCord.cell2Cord(x_pre + med),  _DefaultCord.cell2Cord(x_pre + med - x_min + y)));
        else
            appendValue(cord, _DefaultCord.createCord(_DefaultCord.cell2Cord(x_min), _DefaultCord.cell2Cord(y)));
    return true;
}

//inline bool extendWindow()
//{
//    return previousWindow() && nextWindow();
//}



inline bool nextCord(typename PMRes::HitString & hit, unsigned & currentIt, String<uint64_t> & cord)
{
    std::cout << "nextCord\n";
    std::cout << "length(hit) " << length(hit) << std::endl;
    uint64_t cordLY = _DefaultCord.getCordY(back(cord));
    while (++currentIt < length(hit)) 
    {
        uint64_t tmpCord = _DefaultCord.hit2Cord(hit[currentIt]);
        if(_DefaultCord.getCordY(tmpCord) > cordLY)
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
    if (empty(hit))
        return false;
    else
        appendValue(cord, _DefaultCord.hit2Cord(hit[0]));
    std::cerr << "init" << (cord[0] >> 20) << " " << (cord[0] & 0xfffff) << std::endl;
    return true;
}

//inline bool endWindow (uint64_t cord, unsigned length)
//{
//    if ((cord & cmask) < length - windowSize)
//        return true;
//    else 
//        return false;
//}

inline bool path(String<Dna5> & read, typename PMRes::HitString hit, StringSet<String<int> > & f2, String<uint64_t> & cords)
{
    String<int> f1;
    unsigned currentIt = 0;
    if(!initCord(hit, currentIt, cords))
        return false;
    createFeatures(begin(read), end(read), f1);
    unsigned genomeId = _getSA_i1(_DefaultCord.getCordX(cords[0]));
    std::cerr << "path() genomeId " << genomeId << std::endl;
    while (_DefaultCord.getCordY(back(cords)) < length(read) - window_size)
    {
        
    std::cout << "_DefaultCord.getCordY " << _DefaultCord.getCordY(back(cords)) << " " << length(read) - window_size << std::endl;
        if (!nextWindow(f1, f2[genomeId], cords)) // threshold = 30
        
//==
// need to do previousWindow
//==
//        if (extendWind(f1, f2[genomeId], cords))
            if(!nextCord(hit, currentIt, cords))
                return false;
            
    }
    return true;
}

void path(typename PMRes::HitSet & hits, StringSet<String<Dna5> > & reads, StringSet<String<Dna5> > & genomes, StringSet<String<uint64_t> > & cords)
{
    StringSet<String<int> > f2;
    createFeatures(genomes, f2);
    for (unsigned k = 0; k < length(reads); k++)
    {
        std::cout << "k " << k << std::endl;
        path(reads[k], hits[k], f2, cords[k]);
    }
}




