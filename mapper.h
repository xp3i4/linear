// ==========================================================================
//                           Mapping SMRT reads 
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
#include "pmpfinder.h"

#ifndef SEQAN_HEADER_PACMAPPER_H
#define SEQAN_HEADER_PACMAPPER_H

template <typename TDna = typename MapperBase<>::DefaultAlphabet, 
    typename TSpec = typename MapperBase<>::DefaultShape>
class Mapper {
    typedef MapperBase<TDna, TSpec> Base;
    typedef typename Base::MRecord   Record;
    typedef typename Base::MParm     Parm;
    typedef typename Base::MIndex    Index;
    typedef typename Base::MAnchors    Anchors;
    typedef typename Base::MRes   Res;
    typedef typename Base::MSeq      Seq;
    typedef typename Base::MSeqs     Seqs;
    typedef typename Cord::CordSet      CordSet;
    typedef typename Cord::CordType CordType;
    typedef typename Res::HitSet    HitSet;
    typedef typename Res::HitType   HitType; 

    Record  record;
    Parm    parm;
    Res     res;
    Index   qIndex;
    CordSet cordSet;
    std::ofstream of;

public:
    Mapper();
    Mapper(Options & options);
    Seqs & reads() {return record.seq1;}             
    Seqs & genomes() {return record.seq2;}             
    Parm & mapParm() {return parm;}
    Res & result() {return res;}
    Index & index() {return qIndex;}
    HitSet & hits() {return res.hits;}             //returns hit set
    Pair<HitType, HitType> getAliX(unsigned) const;
    HitType getHitX (HitType const &) const;     //returns coordinates x,y of the hits 
    HitType getHitY (HitType const &) const;     //type uint64_t
    CordSet & cords() {return cordSet;}            //returns cord set 
    CordType getCordX(CordType const &) const;   // returns coordinates x,y of the vertex of sliding window
    CordType getCordY(CordType const &) const;   // type uint64_t 
    
    void printHits();
    void printBestHitsStart();
    void printResult();
    void printParm();
    void printCords(std::ostream & );
    void printCords();
    int createIndex();
     
    //Mapper(Options const & options)
    //{
    //    loadRecords(options);
    //    setMapParm(options);
    //};
};


template <typename TDna, typename TSpec>
Mapper<TDna, TSpec>::Mapper(Options & options):
    record(options),
    parm(options),
    qIndex(genomes()),
    of(toCString(options.getOutputPath()))
{}

template <typename TDna, typename TSpec>
int Mapper<TDna, TSpec>::createIndex()
{
    std::cerr << "Creating index \n";
    _createQGramIndex(qIndex);
    return 0;
}


template <typename TDna, typename TSpec>
Pair<typename Mapper<TDna,TSpec>::HitType, typename Mapper<TDna, TSpec>::HitType>
Mapper<TDna, TSpec>::getAliX(unsigned k) const
{
    typedef typename Mapper<TDna,TSpec>::HitType HitType;
    typedef Pair<HitType, HitType> Pair;
    if (empty(res.hits[k]))
        return Pair::Pair(~(HitType)0, ~(HitType)0);
    else
        return Pair::Pair(_getSA_i1(res.hits[k][0]), _getSA_i2(res.hits[k][0])) << std::endl;
}

template <typename TDna, typename TSpec>
typename Mapper<TDna,TSpec>::HitType 
Mapper<TDna, TSpec>::getHitX(typename Mapper<TDna, TSpec>::HitType const & hit) const
{
    return _DefaultCord.getCordX(_DefaultCord.hit2Cord(hit));
}

template <typename TDna, typename TSpec>
typename Mapper<TDna,TSpec>::HitType 
Mapper<TDna, TSpec>::getHitY(typename Mapper<TDna, TSpec>::HitType const & hit) const
{
    return _DefaultCord.getCordY(_DefaultCord.hit2Cord(hit));
}

template <typename TDna, typename TSpec>
typename Mapper<TDna,TSpec>::CordType 
Mapper<TDna, TSpec>::getCordX(typename Mapper<TDna, TSpec>::CordType const & cord) const
{
    return _DefaultCord.getCordX(cord);
}

template <typename TDna, typename TSpec>
typename Mapper<TDna,TSpec>::CordType 
Mapper<TDna, TSpec>::getCordY(typename Mapper<TDna, TSpec>::CordType const & cord) const
{
    return _DefaultCord.getCordY(cord);
}

template <typename TDna, typename TSpec>
void Mapper<TDna, TSpec>::printHits()
{
    unsigned j = 0;
    std::cout << "printHits: " << lengthSum(res.hits) << " in sum " << std::endl;
    for (auto && hitStr : res.hits)
    {
        std::cout << j++ << "th hits"<< std::endl;
        for (auto && hit : hitStr)
            std::cout << _getSA_i1(getHitX(hit)) << " " << _getSA_i2(getHitX(hit)) << " " << getHitY(hit) << std::endl;
        std::cout << std::endl;
    }
}

template <typename TDna, typename TSpec>
void Mapper<TDna, TSpec>::printBestHitsStart()
{
    unsigned k = 0;
    for (auto && hitStr : res.hits)
        if (empty(hitStr))
            std::cout << std::endl;
        else
            std::cout << k++ << " " << getHitX(hitStr[0]) - getHitY(hitStr[0]) << std::endl;
}

template <typename TDna, typename TSpec>
void Mapper<TDna, TSpec>::printResult()
{}

template <typename TDna, typename TSpec>
void Mapper<TDna, TSpec>::printParm()
{
    parm.print();
}

template <typename TDna, typename TSpec>
void Mapper<TDna, TSpec>::printCords(std::ostream & of)
{
    double time = sysTime();
    unsigned strand;
    for (unsigned k = 0; k < length(cordSet); k++)
    {
        if (empty(cordSet))
            of << k << "th Strand " << strand << " 2 " << length(reads()[k]) << "\n";
        else
        {
            if (_DefaultCord.getCordStrand(back(cordSet[k]))) 
                strand = 1;
            else 
                strand = 0;
            of << k << " " << strand << " " << length(reads()[k]) << "\n";
            _DefaultCord.print(cordSet[k],of);
        }
    }
    std::cerr << ">Write results to disk          " << std::endl;
    std::cerr << "    End writing results. Time[s]" << sysTime() - time << std::endl;
}

template <typename TDna, typename TSpec>
void Mapper<TDna, TSpec>::printCords()
{
    std::cerr << "Writing results to disk \r";
    double time = sysTime();
    unsigned strand;
    for (unsigned k = 0; k < length(cordSet); k++)
    {
        if (empty(cordSet[k]))
            of << k << "th Strand " << "2 length " << length(reads()[k]) << "\n";
        else
        {
            if (_DefaultCord.getCordStrand(back(cordSet[k]))) 
                strand = 1;
            else 
                strand = 0;
            of << k << "th Strand " << strand << " length " << length(reads()[k]) << "\n";
            _DefaultCord.print(cordSet[k], of);
        }
    }
    of.close();
    std::cerr << ">Write results to disk        " << std::endl;
    std::cerr << "    End writing results. Time[s]" << sysTime() - time << std::endl;
}

template <typename TDna, typename TSpec>
void map(Mapper<TDna, TSpec> & mapper)
{
    //printStatus();
    double time = sysTime();
    mapper.createIndex();
    resize(mapper.hits(), length(mapper.reads()));
    resize(mapper.cords(), length(mapper.reads()));
//  mnMap<TDna, TSpec>(mapper.index(), mapper.reads(), _DefaultMapParm, mapper.hits());
// path(mapper.hits(), mapper.reads(), mapper.genomes(), mapper.cords());
    rawMap<TDna, TSpec>(mapper.index(), mapper.reads(), mapper.genomes(), _DefaultMapParm, mapper.hits(), mapper.cords());
    //mapper.printHits();
    mapper.printCords();
    std::cerr << (length(mapper.index().dir) >>27) << " " << (length(mapper.index().sa)>>27) << std::endl;
    std::cerr << "Time in sum[s] " << sysTime() - time << std::endl;
}


#endif
