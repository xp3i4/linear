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
    typedef uint64_t Flag;
    typedef uint64_t CordType;
    typedef uint64_t CellType;
    typedef unsigned Size;
    
    Bit bit;
    Mask mask;
    Mask maskx;
    Flag flag_strand;
    Flag flag_end;
    Bit cell_bit;
    Size cell_size;
    
    CordBase():
        bit(20),    
        mask(0xfffff),
        maskx(0xffffffffff),
        flag_strand(0x2000000000000000),
        flag_end(0x1000000000000000),
        cell_bit(4),
        cell_size(16)
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
    CordType createCord(CordType const &, CordType const &, typename CordBase::Bit const &) const ;
    CordType hit2Cord(PMRes::HitType const &, typename CordBase::Bit const &, typename CordBase::Mask const &) const;
    CellType cord2Cell(CordType const &, typename CordBase::Bit const &) const;
    CordType cell2Cord(CellType const &, typename CordBase::Bit const &) const;
    void setCordEnd(CordType &, typename CordBase::Flag const &, typename CordBase::Flag const &);
    Flag getCordStrand(CordType const &, CordBase::Flag const &) const;
    Flag AtCordEnd(CordType const &, CordBase::Flag const &)const;
    
    
    bool print (CordString const &, std::ostream & = std::cout, CordBase const & = _DefaultCordBase) const;
    bool print (CordSet const &, std::ostream & = std::cout, CordBase const & = _DefaultCordBase) const;
    bool printAlignmentMatrix(CordSet const &,  CordBase const & ) const;
    
}_DefaultCord; 

inline typename Cord::CordType 
Cord::getCordX(typename Cord::CordType const & cord, 
               typename CordBase::Bit const & bit  = _DefaultCordBase.bit,
               typename CordBase::Mask const & mask = _DefaultCordBase.maskx) const
{
    return cord >> bit & mask; 
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
                 typename CordBase::Bit const & bit = _DefaultCordBase.bit) const
{
    return (x << bit) + y;
}

inline typename Cord::CordType 
Cord::hit2Cord(typename PMRes::HitType const & hit, 
               typename CordBase::Bit const & bit = _DefaultCordBase.bit, 
               typename CordBase::Mask const & mask = _DefaultCordBase.mask) const
{
    return hit + ((hit & mask) << bit);
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
            typename CordBase::Flag const & strand = _DefaultCordBase.flag_strand) const
{
    return cord & strand;
}

inline typename CordBase::Flag 
Cord::AtCordEnd(typename Cord::CordType const & cord,
                typename CordBase::Flag const & end = _DefaultCordBase.flag_end) const
{
    return cord & end;
}


inline bool Cord::print(typename Cord::CordString const & cords, std::ostream & of, CordBase const & cordBase) const
{
    of << "length of cords " << length(cords) << std::endl;
    for (auto && j : cords)
       //std::cout << getCordX(j, cordBase.bit) << " , " << getCordY(j, cordBase.mask) << " ";
               of << getCordY(j, cordBase.mask) << " " 
                  << _getSA_i1(getCordX(j, cordBase.bit)) << " "
                  << _getSA_i2(getCordX(j, cordBase.bit))  << std::endl;
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

inline unsigned getAnchorMatch(Anchors & anchors, MapParm & mapParm, float const & thrd, PMRes::HitString & hit, double & time )
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
    return getAnchorMatch(anchors, mapParm, length(read) * 0.25, hit, time2);
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
    typename PMRes::HitString crhit;
    //double time=sysTime();

    std::cerr << "Filtering reads     \r"; 
    
    Anchors anchors(Const_::_LLTMax, AnchorBase::size);
    Seq comStr;
    
    double tm=0;
    double tm2=0;
    unsigned count = 0;
    unsigned countj=0; 
    unsigned ln = length(reads) / 20;
    unsigned jt = 5;
    for (unsigned j = 0; j< length(reads); j++)
    {
//        anchors.init();
//==
// anchor fucntion improving
//==
        anchors.len=1;
//==
//need polishing: only mapping reads > length threshold 
//==
   // if (countj++ == ln)
   // {
   //     std::cerr << "Filtering reads:" << jt << "% \r";
   //     jt += 5; 
   //     countj = 0;
   // }
        if(length(reads[j]) < 1000)
            continue;
//        if (mnMapRead<TDna, TSpec>(index, reads[j], anchors, mapParm, hits[j], tm, tm2) 
//                    < (mapParm.threshold << 20))
if (mnMapRead<TDna, TSpec>(index, reads[j], anchors, mapParm, hits[j], tm, tm2) 
                    < (300 << 20))

        {
            count++;
            _compltRvseStr(reads[j], comStr);
            clear(crhit);
            mnMapRead<TDna, TSpec>(index, comStr, anchors, mapParm, crhit, tm, tm2);
            if (length(crhit) > length(hits[j]))
            {
                clear(hits[j]);
                hits[j] = crhit;
            }

        }

//========
//rs
//========
        //rs.appendValue(_getSA_i2(anchor[res & mask1].getPos1()), _getSA_i2(anchor[res & mask1].getPos1()) + length);
         
    }
    //std::cerr << "Filtering reads Time[s]: " << sysTime() - time << std::endl;
    //std::cerr << "tm [s]: " << tm << std::endl;
    //std::cerr << "tm2 [s]: " << tm2 << std::endl;
    std::cerr << "    rsc strand " << count << std::endl;
}


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
static const unsigned scriptWindow=6; //2^6
//static const uint64_t scriptMask = (1-3*scriptWindow) -1;
static const int scriptCount[5] = {1, 1<<scriptWindow, 1 <<(scriptWindow * 2), 0, 0};
static const int scriptMask = (1 << scriptWindow) - 1;
static const int scriptMask2 = scriptMask << scriptWindow;

static const uint64_t hmask = (1ULL << 20) - 1;

static const unsigned windowThreshold = 25;


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
 //   std::cerr << "createFeatures " << std::endl;
    unsigned next = 1;
    unsigned window = 1 << scriptWindow;
 //   std::cerr << "length() " << itEnd - itBegin << " " << ((itEnd - itBegin -window) >> scriptBit) + 1 << std::endl;
    resize (f, ((itEnd - itBegin -window) >> scriptBit) + 1);
//    std::cerr << length(f) << std::endl;
    f[0] = 0;
    for (unsigned k = 0; k < window; k++)
    {
        f[0] += scriptCount[ordValue(*(itBegin + k))];
    }
    for (unsigned k = scriptStep; k < itEnd - itBegin - window ; k+=scriptStep) 
    {
//    std::cout << k << " " << next << " " << " " << length(f) << std::endl;
        f[next] = f[next - 1];
        for (unsigned j = k - scriptStep; j < k; j++)
        {
            f[next] += scriptCount[ordValue(*(itBegin + j + window))] - scriptCount[ordValue(*(itBegin + j))];
        }
        next++;
    }
//    std::cout << "createFeatures done " << std::endl;

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

inline bool previousWindow(String<int> & f1, String<int> & f2, typename Cord::CordString & cord)
{
    //std::cerr << "previous " << std::endl;
    typedef typename Cord::CordType CordType;
    CordType genomeId = _getSA_i1(_DefaultCord.getCordX(back(cord)));
    CordType x_suf = _DefaultCord.cord2Cell(_getSA_i2(_DefaultCord.getCordX(back(cord))));
    CordType y_suf = _DefaultCord.cord2Cell(_DefaultCord.getCordY(back(cord)));
    //std::cerr << _DefaultCord.getCordY(back(cord)) << " " << y_suf << " " << med << std::endl;
    CordType y;
    if (y_suf < med || x_suf < sup)
        return false;
    else 
       y = y_suf - med;

    unsigned min = ~0;
    unsigned x_min = 0;
    //std::cerr << length(cord) << " " << x_suf << std::endl;
    //std::cerr << "nextWindow() x_pre " << x_pre << std::endl;
    //std::cerr << "pre " << std::endl;
    for (CordType x = x_suf - sup; x < x_suf - inf; x += 1) 
    {
//        std::cerr << "for " << y << " " << length(f1) << " " << length(f2) << " " << x << " " << x_suf - sup << " " << x_suf - inf<< std::endl;
        unsigned tmp = _windowDist(begin(f1) + y, begin(f2) + x);
        
     //   std::cerr << tmp << " ";
        if (tmp < min)
        {
            min = tmp;
            x_min = x;
        }
    }
    //std::cerr << std::endl;
    //std::cerr << "previousWindow " << min << " windowThreshold " << windowThreshold << std::endl;
    if (min > windowThreshold)
        return false;
    else 
    {
        if ( x_suf - x_min > med)
        {
            appendValue(cord, _DefaultCord.createCord(_createSANode(genomeId, _DefaultCord.cell2Cord(x_suf - med)),  _DefaultCord.cell2Cord(x_suf - x_min - med + y)));
     //       std::cerr << "previous " << x_suf - x_min + y << std::endl;
        }
        else
            appendValue(cord, _DefaultCord.createCord(_createSANode(genomeId, _DefaultCord.cell2Cord(x_min)), _DefaultCord.cell2Cord(y)));
      //      std::cerr << "previous y" << y << std::endl;
        //if (_DefaultCord.cord2Cell(_DefaultCord.getCordY(back(cord))) > length(f1))
        //std::cout <<_DefaultCord.cord2Cell(_DefaultCord.getCordY(back(cord))) << std::endl; 
    }
    return true;
}



inline bool nextWindow(String<int> &f1, String<int> & f2, typename Cord::CordString & cord)
{
    typedef typename Cord::CordType CordType;
    CordType genomeId = _getSA_i1(_DefaultCord.getCordX(back(cord)));
    CordType x_pre = _DefaultCord.cord2Cell(_getSA_i2(_DefaultCord.getCordX(back(cord))));
    CordType y_pre = _DefaultCord.cord2Cell(_DefaultCord.getCordY(back(cord)));
    CordType y;
    if (y_pre + sup > length(f1) || x_pre + sup > length(f2))
        return false;
    else 
        y = y_pre + med;
    
    unsigned min = ~0;
    unsigned x_min = 0;
    //std::cout << "nextWindow() x_pre " << x_pre << std::endl;
   // std::cerr << "next " ;
    for (CordType x = x_pre + inf; x < x_pre + sup; x += 1) 
    //for (CordType x = x_pre ; x < x_pre + sup; x += 1) 
    {

//if (y > length(f1)||x > length(f2))
//{
//    std::cerr << "fxy error " << y << " " << length(f1) << " " << x << " " << length(f2) << std::endl;
//}
        unsigned tmp = _windowDist(begin(f1) + y, begin(f2) + x);
//std::cerr << tmp << " ";
        if (tmp < min)
        {
            min = tmp;
            x_min = x;
        }
    }
//    std::cerr << std::endl;
    //std::cout << "nextWindow " << min << " windowThreshold " << _DefaultCord.getCordY(back(cord))<< std::endl;
    if (min > windowThreshold)
       return false;
    else 
        if ( x_min - x_pre > med)
            appendValue(cord, _DefaultCord.createCord(_createSANode(genomeId, _DefaultCord.cell2Cord(x_pre + med)),  _DefaultCord.cell2Cord(x_pre + med - x_min + y)));
        else
            appendValue(cord, _DefaultCord.createCord(_createSANode(genomeId, _DefaultCord.cell2Cord(x_min)), _DefaultCord.cell2Cord(y)));
    return true;
}



inline bool extendWindow(String<int> &f1, String<int> & f2, typename Cord::CordString & cord)
{
    Cord::CordType preCord = (length(cord)==1)?0:back(cord);
    unsigned len = length(cord) - 1;
    //std::cout << preCord << " " << _getSA_i1(_DefaultCord.getCordX(back(cord))) << " " << _getSA_i2(_DefaultCord.getCordX(back(cord))) << std::endl;
    while (preCord < back(cord) && previousWindow(f1, f2, cord)){}
    for (unsigned k = len; k < ((length(cord) + len) >> 1); k++) // when k < length(cord) - 1 - k + len -> swap  to keep cord in assend order
    {
        std::swap(cord[k], cord[length(cord) - k + len - 1]);
       // std::cerr << "swap " << k << " cord[k] " << _getSA_i2(_DefaultCord.getCordY(cord[k])) << " " << length(cord) - k + len - 1 << " " << _getSA_i2(_DefaultCord.getCordY(cord[length(cord) - k + len - 1])) << std::endl;
    }
    while (nextWindow(f1, f2, cord)){}
    return true;
}

inline bool path(String<Dna5> & read, typename PMRes::HitString hit, StringSet<String<int> > & f2, String<uint64_t> & cords)
{
    String<int> f1;
    unsigned currentIt = 0;
    if(!initCord(hit, currentIt, cords))
        return false;
//    std::cerr << "length " << length(read) << std::endl;
    createFeatures(begin(read), end(read), f1);
    unsigned genomeId = _getSA_i1(_DefaultCord.getCordX(cords[0]));
    uint64_t pre = _getSA_i1(_DefaultCord.getCordX(cords[0]));
   // for (unsigned k = 0; k<length(cords);k++)
   // {
   // 
   // if (pre != _getSA_i1(_DefaultCord.getCordX(cords[k])))
   //     std::cerr << "id error " << pre << " " << _getSA_i1(_DefaultCord.getCordX(cords[k])) << std::endl;
   // pre = _getSA_i1(_DefaultCord.getCordX(cords[k]));
   // }
        

//    std::cerr << "done " << std::endl;
    while (_DefaultCord.getCordY(back(cords)) < length(read) - window_size)
    {
        extendWindow(f1, f2[genomeId], cords);
        //std::cerr << "nextCord" << std::endl;
        if(!nextCord(hit, currentIt, cords))
                return false;
    }
    return true;
}

void path(typename PMRes::HitSet & hits, StringSet<String<Dna5> > & reads, StringSet<String<Dna5> > & genomes, StringSet<String<uint64_t> > & cords)
{
    StringSet<String<int> > f2;
    createFeatures(genomes, f2);
    std::cerr << "raw mapping... " << std::endl;
    for (unsigned k = 0; k < length(reads); k++)
    {
        path(reads[k], hits[k], f2, cords[k]);
    }
//    _DefaultCord.print(cords);
//    _DefaultCord.printAlignmentMatrix(cords);
}

template <typename TDna, typename TSpec>
void rawMap(typename PMCore<TDna, TSpec>::Index   & index,
           typename PMRecord<TDna>::RecSeqs      & reads,
            typename PMRecord<TDna>::RecSeqs & genomes,
           MapParm                      & mapParm,
            typename PMRes::HitSet & hits,
            StringSet<String<uint64_t> > & cords)
      //      Anchors & anchors)
           //PMResSet             & rs )
{
//    typedef typename Mapper<TDna, TSpec>::Seq Seq;
    typedef typename PMRecord<TDna>::RecSeq Seq;
    typedef typename PMCore<TDna, TSpec>::Anchors Anchors;
    typename PMRes::HitString crhit;
    String<uint64_t> crcord;
    double time=sysTime();

    std::cerr << "Raw mapping... \r"; 
    
    Anchors anchors(Const_::_LLTMax, AnchorBase::size);
    Seq comStr;
    
    double tm=0;
    double tm2=0;
    unsigned count = 0;
    unsigned countj=0; 
    unsigned ln = length(reads) / 20;
    unsigned jt = 5;

    StringSet<String<int> > f2;
    createFeatures(genomes, f2);

    for (unsigned j = 0; j< length(reads); j++)
    {
//        anchors.init();
//==
// anchor fucntion improving
//==
        anchors.len=1;
//==
//need polishing: only mapping reads > length threshold 
//==
        if(length(reads[j]) < 1000)
            continue;
        mnMapRead<TDna, TSpec>(index, reads[j], anchors, mapParm, hits[j], tm, tm2);
        path(reads[j], hits[j], f2, cords[j]);            
        //_DefaultCord.print(cords[j]);
        if (length(cords[j]) < 20)
        {
            anchors.len=1;
            //anchors.init(Const_::_LLTMax, AnchorBase::size);
            count++;
            _compltRvseStr(reads[j], comStr);
            clear(crhit);
            clear(crcord);
            mnMapRead<TDna, TSpec>(index, comStr, anchors, mapParm, crhit, tm, tm2);
            path(comStr, crhit, f2, crcord);            
           // for (unsigned k = 0; k < length(crhit); k++)
           //     std::cout << _getSA_i2(crhit[k] >> AnchorBase::bit ) << " " << (crhit[k]& AnchorBase::mask) << std::endl;
            if (length(crhit) > length(hits[j]))
            {
                clear(hits[j]);
                clear(cords[j]);
                hits[j] = crhit;
                cords[j] = crcord;
                _DefaultCord.setCordEnd(back(cords[j]));
                //_DefaultCord.print(cords[j]);
            }
            else
            {
                if (!empty(cords[j]))
                    _DefaultCord.setCordEnd(back(cords[j]),0);
            }

        }

//========
//rs
//========
        //rs.appendValue(_getSA_i2(anchor[res & mask1].getPos1()), _getSA_i2(anchor[res & mask1].getPos1()) + length);
         
    }
    std::cerr << "Raw map reads:            " << std::endl;
    std::cerr << "    rsc strand " << count << std::endl;
    std::cerr << "    End raw mapping. Time[s]: " << sysTime() - time << std::endl;
    //std::cerr << "tm [s]: " << tm << std::endl;
    //std::cerr << "tm2 [s]: " << tm2 << std::endl;
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




