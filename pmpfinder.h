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
    Flag flag_strand;
    Flag flag_end;
    Bit cell_bit;
    Size cell_size;
    Mask headFlag;
    
    CordBase():
        bit(20),    
        mask(0xfffff),
        maskx(0xffffffffff),
        valueMask((1ULL<< 60) - 1),
        flag_strand(0x2000000000000000),
        flag_end(0x1000000000000000),
        cell_bit(4),
        cell_size(16),
        headFlag((1ULL<<63))
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
    CordType hit2Cord(PMRes::HitType const &, typename CordBase::Bit const &, typename CordBase::Mask const &, typename CordBase::Mask const &) const;
    CellType cord2Cell(CordType const &, typename CordBase::Bit const &) const;
    CordType cell2Cord(CellType const &, typename CordBase::Bit const &) const;
    void setCordEnd(CordType &, typename CordBase::Flag const &, typename CordBase::Flag const &);
    Flag getCordStrand(CordType const &, CordBase::Flag const &) const;
    Flag AtCordEnd(CordType const &, CordBase::Flag const &)const;
    //void setHead(uint64_t &, uint64_t const &, uint64_t const & = _DefaultCordBase.headFlag);
    void setMaxLen(String<uint64_t> &, uint64_t const &, uint64_t const & = _DefaultCordBase.mask);
    uint64_t getMaxLen(String<uint64_t> const &, uint64_t const & = _DefaultCordBase.mask);

    
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
                 typename CordBase::Bit const & bit = _DefaultCordBase.bit) const
{
    return (x << bit) + y;
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

template <typename TDna, typename TSpec>
inline unsigned getIndexMatch(typename PMCore<TDna, TSpec>::Index & index,
                              typename PMRecord<TDna>::RecSeq const & read,
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
            uint64_t pre = ~0;
            uint64_t pos = getXYDir(index, index.shape.XValue, index.shape.YValue);
            if (_DefaultHs.getHsBodyY(index.ysa[std::min(pos + mapParm.delta, length(index.ysa) - 1)]) ^ index.shape.YValue)
            {
                //while (_DefaultHs.isBodyYEqual(index.ysa[pos], index.shape.YValue))
                while (((index.ysa[pos] >> 41) & ((1ULL << 20) - 1)) == index.shape.YValue && _DefaultHs.isBody(index.ysa[pos]))
                {
//!Note: needs change
                    //std::cout << "[DEBUG] " << (_DefaultHs.getHsBodyS(index.ysa[pos]) >> 30) << " " << (_DefaultHs.getHsBodyS(index.ysa[pos]) & ((1ULL << 30) - 1)) << " " << _DefaultHs.getHsBodyS(index.ysa[pos]) - k << " " << ((_DefaultHs.getHsBodyS(index.ysa[pos]) - k ) & ((1ULL << 30) - 1)) << std::endl;
                    //std::cout << "[DEBUG]" <<  " ref id " << (_DefaultHs.getHsBodyS(index.ysa[pos]) >>30)<< " ref pos " << (_DefaultHs.getHsBodyS(index.ysa[pos]) & ((1ULL << 30) - 1)) << " anchor " << _DefaultHs.getHsBodyS(index.ysa[pos] - k) << " y value "<< _DefaultHs.getHsBodyY(index.ysa[pos])<< " "<< pre - index.ysa[pos] << " pos " << mapParm.kmerStep << "\n";
                    //std::cout << "[DEBUG]" << _DefaultHs.getHsBodyS(pre) << " " << _DefaultHs.getHsBodyY(pre) << " " << _DefaultHs.getHsBodyS(index.ysa[pos]) << " " << _DefaultHs.getHsBodyY(index.ysa[pos])<< " "<< pre - index.ysa[pos] << " pos " << mapParm.kmerStep << "\n";
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
    time += sysTime() - t;
    return time;
}

/*
template <typename TDna, typename TSpec>
inline unsigned getIndexMatch(typename PMCore<TDna, TSpec>::Index & index,
                              typename PMRecord<TDna>::RecSeq const & read,
                              uint64_t* const set,
                              unsigned & len,
                              MapParm & mapParm,
                              double & time
                             )
{    
    double t=sysTime();
    unsigned block = (mapParm.blockSize < length(read))?mapParm.blockSize:length(read);
    std::cerr << " block = " << block << "\n";
    unsigned dt = block * (mapParm.alpha / (1 - mapParm.alpha));
    hashInit(index.shape, begin(read));
    for (unsigned h=0; h <= length(read) - block; h += dt)
    {
        hashInit(index.shape, begin(read) + h);
        for (unsigned k = h; k < h +  block; k++)
        {
            //std::cout << "\n" << k << "\n";
            hashNext(index.shape, begin(read) + k);
            //std::cout << k << " " <<  length(read) << " " << index.shape.hValue << std::endl;
            uint64_t pre = ~0;
            uint64_t pos = getXYDir(index, index.shape.XValue, index.shape.YValue);
                //while (_DefaultHs.isBodyYEqual(index.ysa[pos], index.shape.YValue))
                while (((index.ysa[pos] >> 41) & ((1ULL << 20) - 1)) == index.shape.YValue && _DefaultHs.isBody(index.ysa[pos]))
                {
//!Note: needs change
                    //std::cout << "[DEBUG] " << (_DefaultHs.getHsBodyS(index.ysa[pos]) >> 30) << " " << (_DefaultHs.getHsBodyS(index.ysa[pos]) & ((1ULL << 30) - 1)) << " " << _DefaultHs.getHsBodyS(index.ysa[pos]) - k << " " << ((_DefaultHs.getHsBodyS(index.ysa[pos]) - k ) & ((1ULL << 30) - 1)) << std::endl;
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
    time += sysTime() - t;
    return time;
}
*/

//====End HIndex.getIndexMatch()

//Generic index for getIndexMatch
//to use it please the index type in PMCore to generic qgram index.
/*
 
template <typename TDna, unsigned TSpan,typename TSpec>
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
*/



inline unsigned getAnchorMatch(Anchors & anchors, MapParm & mapParm, float const & thrd, PMRes::HitString & hit, double & time )
{
    double t=sysTime();
    uint64_t ak;
    unsigned maxLen = 0, c_b=mapParm.shapeLen, cbb=0, sb=0, end=0, start = 0;
    unsigned maxlen = 0, start2 = 0;
    unsigned sum = 0;
    anchors[0] = anchors[1];
    ak=anchors[0];
    anchors.sort(anchors.begin(), anchors.end());
    for (unsigned k = 1; k <= anchors.length(); k++)
    {
        //if(thrd == 2366)
        //std::cout << k << " " << (anchors[k]>>20) << " " << (ak>>20)<< " " << start << " " << ((anchors[start])>>20) <<" "<< (anchors[start2] >>20) <<std::endl;
//        if (anchors[k] - ak < AnchorBase::AnchorValue)
        if (anchors[k] - ak < AnchorBase::AnchorValue)
            cbb++;
        else
        {
            anchors.sortPos2(anchors.begin() + sb, anchors.begin() + k);
            //if(thrd == 2366)
            //std::cout << "sb " << sb << " k " << k << std::endl;
            for (uint64_t m = sb+1; m < k; m++)
            {
                if(anchors.deltaPos2(m, m-1) >  mapParm.shapeLen)
                    c_b += mapParm.shapeLen;
                else
                    c_b += anchors.deltaPos2(m, m-1); 
            }
            //if (c_b > 100)
            //std::cout << "[DEBUG] " << (anchors[sb + 1] >> 20) << " " << c_b << " \n";
            if (c_b > maxLen)
            {
                maxlen = maxLen;
                maxLen = c_b;
                start2 = start;
                start = sb;
                end = k;
            }
            if (c_b > 200)
                sum++;
                
            sb = k;
            ak = anchors[k];
           // if (thrd==2366)
           // std::cout << "ak " << (ak >> 20) << " " << k << "\n";
            cbb = 1;
            c_b = mapParm.shapeLen;
        }

    }
    
    //std::cerr << "[DEBUG] len" << maxLen << " "<< start << " " << end << " " << anchors.len<< std::endl;
    for (unsigned k = start; k < end; k++)
        appendValue(hit, anchors[k]);
        uint64_t d =(anchors[start] > anchors[start2])?(anchors[start] - anchors[start2])>>20:(anchors[start2] - anchors[start])>>20;
    //if (thrd == 2366)
    if (sum==0)
        sum=1;
    std::cout << maxlen << " " << sum << " " << d << " " << thrd << " " << (anchors[start2] >> 20) << " " << maxLen << " " << (anchors[start] >> 20) << " "<< thrd * 4 << " " << (float)maxLen /thrd /4<< std::endl;
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
    std::cerr << "length(reads)" << length(reads) << "\n";
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
        std::cout << j << "\n";
        if(length(reads[j]) < 1000)
            continue;
//        if (mnMapRead<TDna, TSpec>(index, reads[j], anchors, mapParm, hits[j], tm, tm2) 
//                    < (mapParm.threshold << 20))
        std::cout << j << "\n";
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
static const unsigned scriptWindow=5; //script_length = 2^scriptWindow
static const unsigned scriptWindow2 = scriptWindow << 1;
//static const uint64_t scriptMask = (1-3*scriptWindow) -1;
static const int scriptCount[5] = {1, 1<<scriptWindow, 1 <<(scriptWindow * 2), 0, 0};
static const int scriptMask = (1 << scriptWindow) - 1;
static const int scriptMask2 = scriptMask << scriptWindow;

static const uint64_t hmask = (1ULL << 20) - 1;

static const unsigned windowThreshold = 35;


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

template<typename TIter> 
inline void createFeatures(TIter const & itBegin, TIter const & itEnd, String<int> & f)
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

inline bool previousWindow(String<int> & f1, String<int> & f2, typename Cord::CordString & cord)
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
            appendValue(cord, _DefaultCord.createCord(_createSANode(genomeId, _DefaultCord.cell2Cord(x_suf - med)),  _DefaultCord.cell2Cord(x_suf - x_min - med + y)));
        }
        else
            appendValue(cord, _DefaultCord.createCord(_createSANode(genomeId, _DefaultCord.cell2Cord(x_min)), _DefaultCord.cell2Cord(y)));
        //if (_DefaultCord.cord2Cell(_DefaultCord.getCordY(back(cord))) > length(f1))
        //std::cout <<_DefaultCord.cord2Cell(_DefaultCord.getCordY(back(cord))) << std::endl; 
    }
    return true;
}

inline bool previousWindow(String<int> & f1, String<int> & f2, typename Cord::CordString & cord, float & score)
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
            appendValue(cord, _DefaultCord.createCord(_createSANode(genomeId, _DefaultCord.cell2Cord(x_suf - med)),  _DefaultCord.cell2Cord(x_suf - x_min - med + y)));
        }
        else
            appendValue(cord, _DefaultCord.createCord(_createSANode(genomeId, _DefaultCord.cell2Cord(x_min)), _DefaultCord.cell2Cord(y)));
        //if (_DefaultCord.cord2Cell(_DefaultCord.getCordY(back(cord))) > length(f1))
        //std::cout <<_DefaultCord.cord2Cell(_DefaultCord.getCordY(back(cord))) << std::endl; 
    }
    score += min;
    return true;
}


inline bool nextWindow(String<int> &f1, String<int> & f2, typename Cord::CordString & cord)
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

inline bool nextWindow(String<int> &f1, String<int> & f2, typename Cord::CordString & cord, float & score)
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
        if (y + 10 > length(f1))
        printf ("[nwe] %d %d %d %d\n", x, y, length(f1), y_pre);
        unsigned tmp = _windowDist(begin(f1) + y, begin(f2) + x);
        if (tmp < min)
        {
            min = tmp;
            x_min = x;
        }
    }
    //printf("[nextCord] %" PRIu64 " %" PRIu64 " %d %d %d\n", x_pre, y_pre, min, x_min, windowThreshold);
    if (min > windowThreshold)
       return false;
    else 
        if ( x_min - x_pre > med)
        {
            //printf("[nextCord] %" PRIu64 " %d %d %d\n", x_pre + med, x_min - x_pre, min, windowThreshold);
            appendValue(cord, _DefaultCord.createCord(_createSANode(genomeId, _DefaultCord.cell2Cord(x_pre + med)),  _DefaultCord.cell2Cord(x_pre + med - x_min + y)));
        }
        else
        {
            //printf("[nextCord] %" PRIu64 " %d %d %d\n", x_pre + med, x_min - x_pre, min, windowThreshold);
            appendValue(cord, _DefaultCord.createCord(_createSANode(genomeId, _DefaultCord.cell2Cord(x_min)), _DefaultCord.cell2Cord(y)));
        }
    score += min;
    //std::cout << "[nxt] " << min << "\n";
    return true;
}

inline bool extendWindow(String<int> &f1, String<int> & f2, typename Cord::CordString & cord)
{
    Cord::CordType preCord = (length(cord)==1)?0:back(cord);
    unsigned len = length(cord) - 1;
    //std::cout << preCord << " " << _getSA_i1(_DefaultCord.getCordX(back(cord))) << " " << _getSA_i2(_DefaultCord.getCordX(back(cord))) << std::endl;
    while (preCord <= back(cord) && previousWindow(f1, f2, cord)){}
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
   // uint64_t pre = _getSA_i1(_DefaultCord.getCordX(cords[0]));
   // for (unsigned k = 0; k<length(cords);k++)
   // {
   // 
   // if (pre != _getSA_i1(_DefaultCord.getCordX(cords[k])))
   //     std::cerr << "id error " << pre << " " << _getSA_i1(_DefaultCord.getCordX(cords[k])) << std::endl;
   // pre = _getSA_i1(_DefaultCord.getCordX(cords[k]));
   // }
        

//    std::cerr << "done " << std::endl;
    std::cout << "[debug] " << _DefaultCord.getCordY(back(cords)) << "\n";
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
    unsigned countl = 0;
    StringSet<String<int> > f2;
std::cerr << "[Debug]done1\n";
    createFeatures(genomes, f2);
std::cerr << "[Debug]done2\n";
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
        //std::cerr << j << "\n";
       // std::cout << "[DEBUG] " << j << "\n";
        if(length(reads[j]) < 1000)
            continue;
        mnMapRead<TDna, TSpec>(index, reads[j], anchors, mapParm, hits[j], tm, tm2);
        path(reads[j], hits[j], f2, cords[j]);            
        //_DefaultCord.print(cords[j]);
        if (anchors.len < 30)
            countl++;
//!Need modify
//!influence speed
        if (length(cords[j]) < (length(reads[j]) / 192 * 0.5))
        //if (length(hits[j]) < 100)
        {
            std::cout << "reverse "<< std::endl;
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
    std::cerr << "[Debug]countl = " << countl << " \n";
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

//========================================
//===thi part is for all map module
//extend the structure Cord;
//NA[3]|head[1]|genome pos[40]|read pos[20]
//NodeType: 1 Head, 0 Body

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
        flag2(1Ull<<bit2),
        mask(flag - 1)
        {}
}_DefaultHitBase;

struct Hit
{
    void setBlockStart(uint64_t &, uint64_t const & = _DefaultHitBase.flag);
    void setBlockBody(uint64_t &, uint64_t const & = _DefaultHitBase.flag);
    bool isBlockStart(uint64_t &, uint64_t const & = _DefaultHitBase.flag);
    void setBlockEnd(uint64_t &, uint64_t const & = _DefaultHitBase.flag);
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
            //std::cout <<"[debug] done\n";
            if (c_b > mapParm.anchorLenThr * readLen)
            {
                anchors.sortPos2(anchors.begin() + sb, anchors.begin() + k);
                //std::cout <<"[debug] done2\n";
                for (unsigned m = sb; m < k; m++)
                {
//to do: replace appendValue? since it's not efficient
                    appendValue(hit, anchors[m]);
                    //std::cout << "[DEBUG]anchors " << _DefaultCord.getCordY(anchors[m]) << "\n";
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

template <typename TDna, typename TSpec>
inline uint64_t mnMapReadAll(typename PMCore<TDna, TSpec>::Index  & index,
                          typename PMRecord<TDna>::RecSeq & read,
                          Anchors & anchors,
                          MapParm & mapParm,
                          PMRes::HitString & hit  
                         )
{
    getIndexMatchAll<TDna, TSpec>(index, read, anchors.set, anchors.len, mapParm);    
    return getAnchorMatchAll(anchors, length(read), mapParm, hit);
}

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

inline bool endCord( String<uint64_t> & cord,
                     unsigned & preCordStart
                   )
{
    _DefaultHit.setBlockEnd(back(cord));
    _DefaultCord.setMaxLen(cord, length(cord) - preCordStart);
    return true;
}

inline bool endCord( String<uint64_t> & cord,
                     unsigned & preCordStart,
                     float const & cordThr,
                     uint64_t const & strand,
                     float & score
                   )
{
    //std::cout << "[endCord] " << length(cord) - preCordStart << "\n";
    if (length(cord) - preCordStart > cordThr)
    {

        _DefaultHit.setBlockEnd(back(cord));
        _DefaultHit.setBlockStrand(cord[preCordStart], strand);
        _DefaultCord.setMaxLen(cord, length(cord) - preCordStart);   
    }
    else
    {
        erase(cord, preCordStart, length(cord));
    }
    //std::cout <<"[nextCord] " << length(cord) - preCordStart << " " << score / (length(cord) - preCordStart) << "\n";
    score = 0;
    return true;
}

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
    //std::cout << "[debug] nextCord() " << length(cord) << "\n";
    if (back(cord)==3390469883040592)
        std::cerr << "[0000000000000000000000000000000000]\n";
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
        //std::cout << "[Cord] " << length(cord) - preCordStart << "\n";
        int64_t b1 = _DefaultCord.getCordX(back(cord)) - _DefaultCord.getCordY(back(cord));
        int64_t b2 = _DefaultCord.getCordX(cord[preCordStart - 1]) - _DefaultCord.getCordY(cord[preCordStart - 1]);
        //std::cout << "[last_delta] " << b1 - b2 << " " << _getSA_i2(_DefaultCord.getCordX(cord[preCordStart]))<<"\n";
            
        if (length(cord) - preCordStart < cordThr )//|| std::abs(b1 - b2) < 10)
        {
            erase(cord, preCordStart, length(cord));
        }
        else
        {
            
         //   for (unsigned k =preCordStart; k<length(cord);k++ )
         //       std::cout << "cord " << score++ << _getSA_i2(_DefaultCord.getCordX(cord[k])) << " " << _DefaultCord.getCordY(cord[k]) << "\n";
         //   std::cout <<"[nextCord] " << length(cord) - preCordStart << " " << score / (length(cord) - preCordStart) << "\n";
            _DefaultCord.setMaxLen(cord, length(cord) - preCordStart);
        }
        preCordStart = length(cord);
        appendValue(cord, _DefaultCord.hit2Cord(*(it)));
        ++it;
        return true;
    }
    else
        return false;
}

/*
//extend the window for the last element in the cord.
inline bool extendWindowAll(String<int> &f1, String<int> & f2, typename Cord::CordString & cord, float & score)
{
    Cord::CordType preCordY = (_DefaultHit.isBlockEnd(cord[length(cord) - 2]))?0:_DefaultCord.getCordY(back(cord)) + window_delta;
    unsigned len = length(cord) - 1;
    while (preCordY<= _DefaultCord.getCordY(back(cord)) && 
        previousWindow(f1, f2, cord, score)){}
    for (unsigned k = len; k < ((length(cord) + len) >> 1); k++) 
    {
        std::swap(cord[k], cord[length(cord) - k + len - 1]);
        // when k < length(cord) - 1 - k + len -> swap, keeping cord in assending order
    }
    //score = _getcord
    while (nextWindow(f1, f2, cord, score)){}
    return true;
}*/
inline bool extendWindowAll(String<int> &f1, String<int> & f2, typename Cord::CordString & cord, float & score)
{
    Cord::CordType preCordY = (_DefaultHit.isBlockEnd(cord[length(cord) - 2]))?0:_DefaultCord.getCordY(back(cord)) + window_delta;
    unsigned len = length(cord) - 1;
   
    while (preCordY<= _DefaultCord.getCordY(back(cord)) && 
        previousWindow(f1, f2, cord, score)){}
    for (unsigned k = len; k < ((length(cord) + len) >> 1); k++) 
    {
        std::swap(cord[k], cord[length(cord) - k + len - 1]);
        // when k < length(cord) - 1 - k + len -> swap, keeping cord in assending order
    }
    //score = _getcord
    while (nextWindow(f1, f2, cord, score)){}
    return true;
}
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
/*

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

template <typename TDna, typename TSpec>
void rawMapAll(typename PMCore<TDna, TSpec>::Index   & index,
            typename PMRecord<TDna>::RecSeqs      & reads,
            typename PMRecord<TDna>::RecSeqs & genomes,
            MapParm & mapParm,
            typename PMRes::HitSet & hits,
            StringSet<String<uint64_t> > & cords)
{
    typedef typename PMRecord<TDna>::RecSeq Seq;
    typedef typename PMCore<TDna, TSpec>::Anchors Anchors;
    typename PMRes::HitString crhit;
    String<uint64_t> crcord;
    double time=sysTime();
    uint64_t count = 0;
    float rcThr = mapParm.rcThr / window_size;
    
    std::cerr << "Raw mapping... \r"; 
    
    Anchors anchors(Const_::_LLTMax, AnchorBase::size);
    Seq comStr;
    
    StringSet<String<int> > f2;
    createFeatures(genomes, f2);
    for (unsigned j = 0; j < length(reads); j++)
    {
        anchors.init(1);
        if (length(reads[j]) < mapParm.minReadLen) // skip reads length < 1000
            continue;
        mnMapReadAll<TDna, TSpec>(index, reads[j], anchors, mapParm, hits[j]);
        pathAll(reads[j], begin(hits[j]), end(hits[j]), f2, cords[j]);
       // for (unsigned k = 0; k < length(hits[j]); k++)
       //     if (_DefaultHit.isBlockStart(hits[j][k]))
       //         std::cerr << "[DEBUG] hit start" <<  k << " " << _DefaultCord.getCordY(hits[j][k]) << "\n";
       //     else 
       //         std::cerr << "[DEBUG] hit " <<  k << " " << _DefaultCord.getCordY(hits[j][k]) << "\n";

            
        //path(reads[j], hits[j], f2, cords[j]);
        
        //std::cerr << "[debug] rc length " << rcThr * length(reads[j]) << "\n";
        if (_DefaultCord.getMaxLen(cords[j]) < rcThr * length(reads[j]))
        {
        
            anchors.init(1);
            count++;
            _compltRvseStr(reads[j], comStr);
            clear(crhit);
            clear(crcord);
            mnMapReadAll<TDna, TSpec>(index, comStr, anchors, mapParm, crhit);
            pathAll(comStr, begin(crhit), end(crhit), f2, crcord);            
            //std::cerr <<"[debug] length " << length(crcord) << " " << _DefaultCord.getMaxLen(crcord)<< std::endl;
            if (_DefaultCord.getMaxLen(crcord) > _DefaultCord.getMaxLen(cords[j]))
            {
                //clear(hits[j]);
                //clear(cords[j]);
                hits[j] = crhit;
                cords[j] = crcord;
                _DefaultCord.setCordEnd(back(cords[j]));
            }
            else
            {
                if (!empty(cords[j]))
                    _DefaultCord.setCordEnd(back(cords[j]),0);
            }
        }
        
    }
    std::cerr << "Raw map reads:            " << std::endl;
    std::cerr << "    rsc strand " << count << std::endl;
    std::cerr << "    End raw mapping. Time[s]: " << sysTime() - time << std::endl;
}

template <typename TDna, typename TSpec>
void rawMapAllComplex(typename PMCore<TDna, TSpec>::Index   & index,
            typename PMRecord<TDna>::RecSeqs      & reads,
            typename PMRecord<TDna>::RecSeqs & genomes,
            MapParm & mapParm,
            typename PMRes::HitSet & hits,
            StringSet<String<uint64_t> > & cords)
{
    typedef typename PMRecord<TDna>::RecSeq Seq;
    typedef typename PMCore<TDna, TSpec>::Anchors Anchors;
    typename PMRes::HitString hit, crhit;
    String<uint64_t> crcord;
    String<unsigned> missId;
    double time=sysTime();
    uint64_t count = 0;
    float rcThr = mapParm.rcThr / window_size;
    float senThr = window_size / 0.85;
    mapParm.alpha = 0.65;
    MapParm complexParm = mapParm;
    complexParm.alpha = 0.5;
    std::cerr << "Raw mapping... \r"; 
    
    Anchors anchors(Const_::_LLTMax, AnchorBase::size);
    Seq comStr;
    
    StringSet<String<int> > f2;
    createFeatures(genomes, f2);
    for (unsigned j = 0; j < length(reads); j++)
    {
        anchors.init(1);
        if (length(reads[j]) < mapParm.minReadLen) // skip reads length < 1000
            continue;
        clear(crhit);
        //mnMapReadAll<TDna, TSpec>(index, reads[j], anchors, mapParm, hits[j]);
        mnMapReadAll<TDna, TSpec>(index, reads[j], anchors, mapParm, crhit);
        pathAll(reads[j], begin(crhit), end(crhit), f2, cords[j]);
        //if (_DefaultCord.getMaxLen(cords[j]) < rcThr * length(reads[j]))
        //{
        
            anchors.init(1);
            count++;
            _compltRvseStr(reads[j], comStr);
            clear(crhit);
            clear(crcord);
            mnMapReadAll<TDna, TSpec>(index, comStr, anchors, mapParm, crhit);
            //pathAll(comStr, begin(crhit), end(crhit), f2, crcord);            
            pathAll(comStr, begin(crhit), end(crhit), f2, cords[j]);            
            //shrinkToFit(cords[j]);
            
            if (_DefaultCord.getMaxLen(crcord) > _DefaultCord.getMaxLen(cords[j]))
            {
                hits[j] = crhit;
                cords[j] = crcord;
                _DefaultCord.setCordEnd(back(cords[j]));
            }
            else
            {
                if (!empty(cords[j]))
                    _DefaultCord.setCordEnd(back(cords[j]),0);
            }   
            
        //}
        //if (_DefaultCord.getMaxLen(crcord) < length(reads[j]) / senThr && 
         //       _DefaultCord.getMaxLen(cords[j]) < length(reads[j]) / senThr) 
         //std::cout << _DefaultCord.getMaxLen(cords[j]) << " " << length(reads[j]) / senThr << std::endl;
         if (_DefaultCord.getMaxLen(cords[j]) < length(reads[j]) / senThr)
            {
                appendValue(missId, j);
                clear(cords[j]);
            }
    }
    std::cout << "welcome 2 complex " << senThr << " "<< length(missId) << " \n";
    for (unsigned k = 0; k < length(missId); k++)
    {
        anchors.init(1);
        if (length(reads[missId[k]]) < complexParm.minReadLen) // skip reads length < 1000
            continue;
        clear(crhit);
        //mnMapReadAll<TDna, TSpec>(index, reads[missId[k]], anchors, complexParm, hits[missId[k]]);
        mnMapReadAll<TDna, TSpec>(index, reads[missId[k]], anchors, complexParm, crhit);
        pathAll(reads[missId[k]], begin(crhit), end(crhit), f2, cords[missId[k]]);

        if (_DefaultCord.getMaxLen(cords[missId[k]]) < rcThr * length(reads[missId[k]]))
        {
        
            anchors.init(1);
            count++;
            _compltRvseStr(reads[missId[k]], comStr);
            clear(crhit);
            clear(crcord);
            mnMapReadAll<TDna, TSpec>(index, comStr, anchors, complexParm, crhit);
            pathAll(comStr, begin(crhit), end(crhit), f2, crcord);            
            if (_DefaultCord.getMaxLen(crcord) > _DefaultCord.getMaxLen(cords[missId[k]]))
            {
                //hits[missId[k]] = crhit;
                cords[missId[k]] = crcord;
                _DefaultCord.setCordEnd(back(cords[missId[k]]));
            }
            else
            {
                if (!empty(cords[missId[k]]))
                    _DefaultCord.setCordEnd(back(cords[missId[k]]),0);
            }   
            
        }
    }
        
    
    
    std::cerr << "Raw map reads:            " << std::endl;
    std::cerr << "    rsc strand " << count << std::endl;
    std::cerr << "    End raw mapping. Time[s]: " << sysTime() - time << std::endl;
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

void _printHit(unsigned j, unsigned id1, unsigned id2, String<uint64_t> & hit, unsigned len)
{
    unsigned end;
    for (unsigned k = 0; k < length(hit); k++)
    {
        if (_DefaultHit.isBlockEnd(hit[k]))
            end = 1;
        else
            end = 0;
        printf("[printhit] %d %d %d % " PRIu64 " %" PRIu64 " %" PRIu64 " %" PRIu64 " %" PRIu64 " %d %d\n", j, id1, id2, hit[k], _getSA_i1(_DefaultCord.getCordX(hit[k])), _getSA_i2(_DefaultCord.getCordX(hit[k])), _DefaultCord.getCordY(hit[k]), end, _DefaultCord.getMaxLen(hit), len);
    }
}


template <typename TDna, typename TSpec>
void rawMapAllComplex2(typename PMCore<TDna, TSpec>::Index   & index,
            typename PMRecord<TDna>::RecSeqs      & reads,
            typename PMRecord<TDna>::RecSeqs & genomes,
            MapParm & mapParm,
            typename PMRes::HitSet & hits,
            StringSet<String<uint64_t> > & cords)
{
    std::cerr << "complex2 \n";
    typedef typename PMRecord<TDna>::RecSeq Seq;
    typedef typename PMCore<TDna, TSpec>::Anchors Anchors;
    typename PMRes::HitString hit, crhit;
    String<uint64_t> crcord;
    String<unsigned> missId;
    double time=sysTime();
    float rcThr = mapParm.rcThr / window_size;
    float senThr = window_size /0.85;
    //float senThr = window_size / 0.85;
    //mapParm.alpha = 0.6;
    MapParm complexParm = mapParm;
    //complexParm.alpha = 0.5;
    complexParm.alpha = complexParm.alpha2;
    std::cerr << "complex 2\n";
    std::cerr << "Raw mapping... \r"; 
    //float cordThr = 0.9;
    
    std::cerr << "mapParm \n";
    mapParm.print();
    std::cerr << "\ncompexParm \n";
    complexParm.print();
    Anchors anchors(Const_::_LLTMax, AnchorBase::size);
    Seq comStr;
    
    StringSet<String<int> > f2;
    createFeatures(genomes, f2);
//    SEQAN_OMP_PRAGMA(parallel for)

    std::cerr << "complex2 \n";
    for (unsigned j = 0; j < length(reads); j++)
    {
        if (length(reads[j]) < mapParm.minReadLen) // skip reads length < 1000
            continue;
        anchors.init(1);
        clear(crhit);
        mnMapReadAll<TDna, TSpec>(index, reads[j], anchors, mapParm, crhit);
        pathAll(reads[j], begin(crhit), end(crhit), f2, cords[j], mapParm.cordThr, 0);
        _printHit(crhit);
        anchors.init(1);
        _compltRvseStr(reads[j], comStr);
        std::cout << "[reverse]\n";
        clear(crhit);
        mnMapReadAll<TDna, TSpec>(index, comStr, anchors, mapParm, crhit);
        pathAll(comStr, begin(crhit), end(crhit), f2, cords[j], mapParm.cordThr, 1);            
        shrinkToFit(cords[j]);
        _printHit(crhit);
        if (_DefaultCord.getMaxLen(cords[j]) < length(reads[j]) / senThr)
        {
            appendValue(missId, j);
            clear(cords[j]);
        }
    }
   
    for (unsigned k = 0; k < length(missId); k++)
    {
        anchors.init(1);
        if (length(reads[missId[k]]) < complexParm.minReadLen) // skip reads length < minReadLen
            continue;
        clear(crhit);
        //mnMapReadAll<TDna, TSpec>(index, reads[missId[k]], anchors, complexParm, hits[missId[k]]);
        mnMapReadAll<TDna, TSpec>(index, reads[missId[k]], anchors, complexParm, crhit);
        pathAll(reads[missId[k]], begin(crhit), end(crhit), f2, cords[missId[k]], mapParm.cordThr, 0);

            anchors.init(1);
            _compltRvseStr(reads[missId[k]], comStr);
            clear(crhit);
            mnMapReadAll<TDna, TSpec>(index, comStr, anchors, complexParm, crhit);
            pathAll(comStr, begin(crhit), end(crhit), f2, cords[missId[k]], mapParm.cordThr, 1);            
            shrinkToFit(cords[missId[k]]);
    }
    
    std::cerr << "Raw map reads:            " << std::endl;
    std::cerr << "    End raw mapping. Time[s]: " << sysTime() - time << std::endl;
}

inline StringSet<String<uint64_t> > & _join(StringSet<String<uint64_t> > & s1,
            StringSet<String<uint64_t> >& s2)
{
        append(s1, s2);
        return s1;
}

inline void _appendMissHit(String<unsigned> & missId, unsigned id, unsigned maxLen, float thr)
{
    if (maxLen < thr)
    {
        appendValue(missId, id);
    }
}

template <typename TDna, typename TSpec>
void rawMapAllComplex2Parallel(typename PMCore<TDna, TSpec>::Index   & index,
            typename PMRecord<TDna>::RecSeqs      & reads,
            typename PMRecord<TDna>::RecSeqs & genomes,
            MapParm & mapParm,
            typename PMRes::HitSet & hits,
            StringSet<String<uint64_t> > & cords,
            unsigned thread)
{
    typedef typename PMRecord<TDna>::RecSeq Seq;
    typedef typename PMCore<TDna, TSpec>::Anchors Anchors;
    typename PMRes::HitString hit;
    String<uint64_t> crcord;
    String<unsigned> missId;
    double time=sysTime();
    float rcThr = mapParm.rcThr / window_size;
    float senThr = window_size /0.85;
    //float senThr = window_size / 0.85;
    //mapParm.alpha = 0.6;
    MapParm complexParm = mapParm;
    //complexParm.alpha = 0.5;
    complexParm.alpha = complexParm.alpha2;
    std::cerr << "Raw Mapping \n"; 
    //float cordThr = 0.9;
    
    std::cerr << "mapParm \n";
    mapParm.print();
    std::cerr << "\ncompexParm \n";
    complexParm.print();
    Anchors anchors(Const_::_LLTMax, AnchorBase::size);
    Seq comStr;
    typename PMRes::HitString crhit;

    StringSet<String<int> > f2;
    createFeatures(genomes, f2);

    omp_set_num_threads(thread);
    
//#pragma omp declare reduction(_joinCords : StringSet<String<uint64_t> > : omp_out = _join(omp_out, omp_in)) 
    
#pragma omp parallel//reduction(_joinCords: cords)
{
    unsigned size2 = length(reads) / omp_get_num_threads();
    unsigned ChunkSize = size2;
    Seq comStr;
    Anchors anchors(Const_::_LLTMax, AnchorBase::size);
    typename PMRes::HitString crhit;
    StringSet<String<uint64_t> >  cordsTmp;
    String<unsigned> missIdTmp;
    if (omp_get_thread_num() < length(reads) - size2 * omp_get_num_threads())
    {
        ChunkSize = size2 + 1;
    }
    resize(cordsTmp, ChunkSize);
    //printf("[threads] %d %d %d\n", omp_get_thread_num(), ChunkSize, length(cordsTmp));
    
    unsigned c = 0;

    #pragma omp for
    for (unsigned j = 0; j < length(reads); j++)
    {
        if (length(reads[j]) >= mapParm.minReadLen )
        {
            anchors.init(1);
            clear(crhit);
            mnMapReadAll<TDna, TSpec>(index, reads[j], anchors, mapParm, crhit);
            pathAll(reads[j], begin(crhit), end(crhit), f2, cordsTmp[c], mapParm.cordThr, 0);
            anchors.init(1);
            _compltRvseStr(reads[j], comStr);
            clear(crhit);
            mnMapReadAll<TDna, TSpec>(index, comStr, anchors, mapParm, crhit);
            pathAll(comStr, begin(crhit), end(crhit), f2, cordsTmp[c], mapParm.cordThr, 1);            
            //shrinkToFit(cords[c]);
            //_appendMissHit(missIdTmp, j, _DefaultCord.getMaxLen(cordsTmp[c]), length(reads[j]) / senThr);
            //_printHit(j, omp_get_thread_num(), c, cordsTmp[c]);
            //printf("[misshit] %d %f\n", _DefaultCord.getMaxLen(cordsTmp[c]), length(reads[j]) / senThr );
            //appendValue(missIdTmp, j);
           if (_DefaultCord.getMaxLen(cordsTmp[c]) < length(reads[j]) / senThr)
           {
               appendValue(missId, j);
               clear(cordsTmp[c]);
           }   
        }
        c += 1;
    } 
    #pragma omp for ordered
    for (unsigned j = 0; j < omp_get_num_threads(); j++)
        #pragma omp ordered
        {
            append(cords, cordsTmp);
            append(missId, missIdTmp);
        }
}

    

    printf("[missId] %d\n", length(missId));

   
    for (unsigned k = 0; k < length(missId); k++)
    {

        anchors.init(1);
        if (length(reads[missId[k]]) < complexParm.minReadLen) // skip reads length < minReadLen
            continue;
        clear(crhit);
        clear(cords[missId[k]]);
        //mnMapReadAll<TDna, TSpec>(index, reads[missId[k]], anchors, complexParm, hits[missId[k]]);
        mnMapReadAll<TDna, TSpec>(index, reads[missId[k]], anchors, complexParm, crhit);
        pathAll(reads[missId[k]], begin(crhit), end(crhit), f2, cords[missId[k]], mapParm.cordThr, 0);

            anchors.init(1);
            _compltRvseStr(reads[missId[k]], comStr);
            clear(crhit);
            mnMapReadAll<TDna, TSpec>(index, comStr, anchors, complexParm, crhit);
            pathAll(comStr, begin(crhit), end(crhit), f2, cords[missId[k]], mapParm.cordThr, 1);            
            shrinkToFit(cords[missId[k]]);
    }
    
    std::cerr << "Raw map reads:            " << length(missId) << std::endl;
    std::cerr << "    End raw mapping. Time[s]: " << sysTime() - time << std::endl;
}



//End all mapper module
//============================================


