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
struct Mapper {
    typedef MapperBase<TDna, TSpec> Base;
    typedef typename Base::MRecord   Record;
    typedef typename Base::MParm     Parm;
    typedef typename Base::MIndex    Index;
    typedef typename Base::MAnchors    Anchors;
    typedef typename Base::MRes   Res;
    typedef typename Base::MSeq      Seq;
    typedef typename Base::MSeqs     Seqs;
    typedef typename Cord::CordSet      CordSet;

    Record  record;
    Parm    parm;
    Res     res;
    Index   qIndex;
    CordSet cordSet;

    //StringSet<String<uint64_t> > cords;

    Mapper();
    Mapper(Options & options);
    Seqs & reads(){return record.seq1;}
    Seqs & genomes(){return record.seq2;}
    Parm & mapParm(){return parm;}
    Res & result(){return res;}
    Index & index(){return qIndex;}
    CordSet & cords(){return cordSet;}
    void printHits();
    void printResult();    
    void printParm();
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
    qIndex(genomes())
{}

template <typename TDna, typename TSpec>
int Mapper<TDna, TSpec>::createIndex()
{
    std::cerr << "Creating index \n";
    _createQGramIndex(qIndex);
    return 0;
}

template <typename TDna, typename TSpec>
void Mapper<TDna, TSpec>::printHits()
{
    std::cout << "Hits: " << lengthSum(res.hits) << " in sum " << std::endl;
    for (auto && hitStr : res.hits)
    {
        for (auto && hit : hitStr)
            std::cout << _DefaultCord.getCordX(_DefaultCord.hit2Cord(hit)) << " " << _DefaultCord.getCordY(_DefaultCord.hit2Cord(hit)) << ", ";
        std::cout << std::endl;
    }
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
void map(Mapper<TDna, TSpec> mapper)
{
    //map.printParm();
    std::cerr << "Encapsulated version " << std::endl;
    double time = sysTime();
    _DefaultMapParm.print();
    mapper.createIndex();
    resize(mapper.res.hits, length(mapper.reads()));
    resize(mapper.cords(), length(mapper.reads()));
    mnMap<TDna, TSpec>(mapper.index(), mapper.reads(), _DefaultMapParm, mapper.res.hits);//, map.result());
    //mapper.printHits();
    //path(mapper.res.hits, mapper.reads(), mapper.genomes(), mapper.cords());
    
  //  _DefaultCord.print(mapper.cords);
    
    std::cerr << "Time of mapping in sum [s] " << sysTime() - time << std::endl;
    mapper.printResult();
}

#endif
