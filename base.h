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

#ifndef SEQAN_HEADER_BASE_H
#define SEQAN_HEADER_BASE_H

#include <seqan/seq_io.h>
#include <seqan/stream.h>
#include <seqan/index.h>
#include <seqan/store.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <seqan/basic.h>
#include <bitset>
#include <climits>
#include <seqan/arg_parse.h>

#include "shape_extend.h"
#include "index_extend.h"

using namespace seqan;

//===================================================================
// variable and type def
//===================================================================
//template <typename TSpec = void>

struct status
{
    double time;
 
};

//struct ConstBase{
//    
//    const unsigned _SHAPELEN;
//    const unsigned _SHAPEWHT;
//    const unsigned _BLOCKSIZE;
//    const unsigned _DELAT; 
//    const unsigned _THRESHOLD; 
//    const float    _ALPHA ;     
//    const unsigned _KMERSTEP;
//    const uint64_t _LLTMax;
//
//    ConstBase()
//        _SHAPELEN(25),
//        _SHAPEWHT(18),
//        _BLOCKSIZE(1000),
//        _DELAT(32),
//        _THRESHOLD(30),
//        _ALPHA(0.8),
//        _KMERSTEP(1000),
//        _LLTMax(~0)
//        {}
//    
//}_DefaultConstBase;

struct Const_{
    
    typedef uint64_t    BIT_INT_;
    typedef CharString  PATH_;
    typedef seqan::Dna5 DEFAULT_ALPHABET_ ;
    typedef CharString  ID_;

    static const unsigned _SHAPELEN;
    static const unsigned _SHAPEWHT;
    static const unsigned _BLOCKSIZE;
    static const unsigned _DELAT; 
    static const unsigned _THRESHOLD; 
    static const float    _ALPHA ;
    static const unsigned _KMERSTEP;
    static const uint64_t _LLTMax;
        
};

const float Const_::_ALPHA = 0.8;
const unsigned Const_::_SHAPELEN = 25;
const unsigned Const_::_SHAPEWHT = 18;
const unsigned Const_::_BLOCKSIZE = 100;
const unsigned Const_::_DELAT = 32; 
const unsigned Const_::_THRESHOLD = 30; 
const unsigned Const_::_KMERSTEP = 1000;
const uint64_t Const_::_LLTMax = ~0;
     
struct Options{

    unsigned    kmerLen;
    unsigned    MiKmLen;
    typename    Const_::PATH_ rPath;
    typename    Const_::PATH_ gPath;
    typename    Const_::PATH_ oPath;
    bool        Sensitive; 
    Options():
        kmerLen(Const_::_SHAPELEN),
        MiKmLen(Const_::_SHAPEWHT),
        rPath(""),
        gPath(""),
        oPath(""),
        Sensitive(false)
        {}
    Const_::PATH_ getGenomePaht() const {return gPath;};
    Const_::PATH_ getReadPaht() const {return rPath;};
    Const_::PATH_ getOutPutPaht() const {return oPath;};
    int print();
}; 

template <typename TDna = Const_::DEFAULT_ALPHABET_>
struct RecordBase
{
    typedef Const_::DEFAULT_ALPHABET_ DefaultAlphabet;
    typedef CharString RecId;
    typedef String<TDna> RecSeq; 
};

template<typename TDna = typename RecordBase<>::DefaultAlphabet>
struct PMRecord
{
    typedef typename RecordBase<TDna>::RecId RecId;
    typedef typename RecordBase<TDna>::RecSeq RecSeq;
    typedef StringSet<RecId> RecIds;
    typedef StringSet<RecSeq> RecSeqs;

    PMRecord(){}
    PMRecord(Options & options);
    RecIds id1, id2;
    RecSeqs seq1, seq2; //seq1=read, seq2=ref

    int loadRecord(Options & options);
};

struct AnchorBase{
    typedef typename Const_::BIT_INT_ AnchorType; 
    static const unsigned size = 131072;
    static const unsigned bit = 20;
    static const uint64_t AnchorValue = 1000ULL << bit;
    static const typename Const_::BIT_INT_ mask = (1ULL<<bit) - 1;
};

struct Anchors{
    typedef typename AnchorBase::AnchorType AnchorType;
    typedef typename AnchorBase::AnchorType * Iter;

    AnchorType set[AnchorBase::size];
    unsigned len;

    Anchors(){len = 1;};
    Anchors(AnchorType val, unsigned range);
    void init(AnchorType val, unsigned k);
    void init();
    void setAnchor(unsigned p, AnchorType pos1,  AnchorType pos2);
    AnchorType getPos1(unsigned p) const;
    AnchorType getPos2(unsigned p) const;
    AnchorType deltaPos1(unsigned p1, unsigned p2);
    AnchorType deltaPos2(unsigned p1, unsigned p2);
    void sort(Iter begin, Iter end);
    void sortPos2(Iter begin, Iter end);
    void appendValue(AnchorType val);
    void appendValue(AnchorType val1, AnchorType val2);
    unsigned size() const {return AnchorBase::size;};
    AnchorType & operator [](unsigned p){return set[p];};
    Iter begin(); 
    Iter end();
    unsigned & length() {return len;};
};


template <typename TDna = Const_::DEFAULT_ALPHABET_, 
        typename CoreMinimizer = Minimizer<Const_::_SHAPELEN> > 
struct CoreBase{
    typedef typename Const_::DEFAULT_ALPHABET_ DefaultAlphabet;
    typedef Minimizer<Const_::_SHAPELEN> DefaultShape;
    typedef typename PMRecord<TDna>::RecSeqs RecSeqs; 
    typedef Shape<TDna, CoreMinimizer> CoreShape;
    typedef Index<RecSeqs, IndexQGram<CoreMinimizer, OpenAddressing> > CoreIndex;
    typedef Anchors AnchorSet;

};

template <typename TDna = typename CoreBase<>::DefaultAlphabet, 
        typename Minimizer = typename CoreBase<>::DefaultShape>
struct PMCore
{
    typedef typename CoreBase<TDna, Minimizer>::CoreIndex Index;
    typedef typename CoreBase<TDna, Minimizer>::RecSeqs Seqs;
    typedef typename CoreBase<TDna, Minimizer>::AnchorSet Anchors;

    //Index index;
    //Anchor anchor;
};


//template <typename TSpec = void>
struct ResBase{
    typedef unsigned SeqLen;
    typedef typename Const_::BIT_INT_ SeqId;
    typedef typename Const_::BIT_INT_ MapPos;
    typedef typename Const_::BIT_INT_ MapScore;
    typedef typename Const_::BIT_INT_ HitType;
    typedef bool MapStrand;

    static const unsigned bit = 32;
    static const Const_::BIT_INT_ mask = (1ULL << bit) - 1;
    static const unsigned hitBit = AnchorBase::bit;
    static const unsigned hitMask = AnchorBase::mask;
    
};

struct PMRes
{
    typedef typename ResBase::SeqId Id;
    typedef typename ResBase::MapPos Pos;
    typedef typename ResBase::MapScore Score;
    typedef typename ResBase::MapStrand Strand;
    
    typedef typename ResBase::HitType  HitType;
    typedef String<HitType> HitString;
    typedef StringSet<HitString> HitSet;

//    String<Id>    id;
//    String<Pos>   pos;
    StringSet<String<Pos> > pos;
    StringSet<String<Score> > score;  
    StringSet<String<Strand> > strand;
    HitSet hits;

    PMRes(){};
    Id getId1(unsigned); 
    Id getId2(unsigned);
    Pos getPosBegin(unsigned);
    Pos getPosEnd(unsigned);
    Strand getStrand();
    Id createId(Id, Id);
    Pos createPos(Pos, Pos);
    void appendValue(unsigned, Pos, Score, Strand);

};

//inline PMRes::Pos PMRes::getPosBegin(unsigned k)
//{
//    return pos[k] >> ResBase::bit;
//};
//
//inline PMRes::Pos PMRes::getPosEnd(unsigned k)
//{
//    return pos[k] & ResBase::mask;
//};
//
//inline PMRes::Id createId(PMRes::Id id1, PMRes::Id id2)
//{
//    return (id1 << ResBase::bit) + id2;
//}
//
//inline PMRes::Pos createPos(PMRes::Pos p1, PMRes::Pos p2)
//{
//    return (p1 << ResBase::bit) + p2;
//}

inline void PMRes::appendValue(unsigned id, PMRes::Pos rpos, PMRes::Score rscore = 0, PMRes::Strand rstrand = false)
{
    //seqan::appendValue(id, rid);
    seqan::appendValue(pos[id], rpos);
    seqan::appendValue(score[id], rscore);
    seqan::appendValue(strand[id], rstrand);
}

struct MapParm{
    unsigned    blockSize;
    unsigned    delta;
    unsigned    threshold;
    unsigned    kmerStep;
    unsigned    shapeLen;
    float       alpha;  
    
    MapParm():
        blockSize(Const_::_BLOCKSIZE),
        delta(Const_::_DELAT),
        threshold(Const_::_THRESHOLD),
        kmerStep(Const_::_KMERSTEP),
        shapeLen(Const_::_SHAPELEN),
        alpha(Const_::_ALPHA)
        {}
// ====
//temp: need modify
    MapParm(Options & options){MapParm();}
    MapParm(MapParm & parm):
        blockSize(parm.blockSize),
        alpha(parm.alpha),
        delta(parm.delta),
        threshold(parm.threshold),
        kmerStep(parm.kmerStep),
        shapeLen(parm.shapeLen)
        {}
    void setMapParm(Options & options);
    void print ();
    
}_DefaultMapParm;

template <typename TDna = typename Const_::DEFAULT_ALPHABET_, 
        typename TSpec = Minimizer<Const_::_SHAPELEN> >
struct MapperBase
{
    typedef Const_::DEFAULT_ALPHABET_ DefaultAlphabet;
    typedef Minimizer<Const_::_SHAPELEN> DefaultShape;
    typedef PMRecord<TDna>  MRecord;
    typedef PMRes           MRes;
    typedef MapParm          MParm;
    typedef PMCore<TDna, TSpec>    MCore;
    typedef typename PMCore<TDna, TSpec>::Index MIndex;
    typedef typename PMCore<TDna, TSpec>::Anchors MAnchors;
    typedef typename PMRecord<TDna>::RecSeq MSeq; 
    typedef typename PMRecord<TDna>::RecSeqs MSeqs;
};
/*
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

    Record  record;
    Parm    parm;
    Res     res;
    Index   qIndex;

    StringSet<String<uint64_t> > cords;

    Mapper();
    Mapper(Options & options);
    Seqs & reads(){return record.seq1;};
    Seqs & genomes(){return record.seq2;};
    Parm & mapParm(){return parm;};
    Res & result(){return res;};
    Index & index(){return qIndex;};
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
*/
int Options::print()
{
    
    std::cerr << "kmerLen " << kmerLen << std::endl
              << "MiKmLen " << MiKmLen << std::endl
              << "reads path " << rPath << std::endl
              << "genomes Path " << gPath << std::endl
              << "output path " << oPath << std::endl
              << "Sensitive " << Sensitive << std::endl;
    return 0;
}


template <typename TDna>
int PMRecord<TDna>::loadRecord(Options & options)
{
    double time = sysTime();
    std::cerr << "loading sequences from files... " << std::endl;
    SeqFileIn rFile(toCString(options.rPath));
    SeqFileIn gFile(toCString(options.gPath));
    //try{
        
        seqan::readRecords(id1, seq1, rFile);
            
        seqan::readRecords(id2, seq2, gFile);
    //   throw 
    //}
    std::cerr << "loading [s] " << sysTime() - time << std::endl;
    return 0;
}

template <typename TDna>
PMRecord<TDna>::PMRecord(Options & options)
{
   loadRecord(options);
}

Anchors::Anchors(Anchors::AnchorType val, unsigned range)
{
    init(val, range);
}

inline void Anchors::init(AnchorType val, unsigned range){
    for (unsigned k = 0; k < range; k++)
        set[k] = val;
    len = 1;
}

inline void Anchors::init()
{
    init(AnchorType(0), size());
}

inline void Anchors::setAnchor(unsigned p, 
    Anchors::AnchorType pos1,  Anchors::AnchorType pos2)
{
    set[p] = (pos1 << AnchorBase::bit) + pos2;
}

inline Anchors::AnchorType Anchors::getPos1(unsigned p) const 
{
    return set[p] >> AnchorBase::bit;
}

inline Anchors::AnchorType Anchors::getPos2(unsigned p) const
{
    return set[p] & AnchorBase::mask;
}

inline Anchors::AnchorType Anchors::deltaPos1(unsigned p1, unsigned p2)
{
    return (set[p1] >> AnchorBase::bit) - (set[p2] >> AnchorBase::bit);
}

inline Anchors::AnchorType Anchors::deltaPos2(unsigned p1, unsigned p2)
{
    return AnchorBase::mask & (set[p1] - set[p2]);
}

inline Anchors::Iter Anchors::begin()
{
    return set;
}

inline Anchors::Iter Anchors::end()
{
    return set + len;
}

inline void Anchors::sort(Anchors::Iter sortBegin, Anchors::Iter sortEnd)
{
    std::sort(sortBegin, sortEnd);
}

inline void Anchors::sortPos2(Anchors::Iter sortBegin, Anchors::Iter sortEnd){
    AnchorBase::AnchorType mask = AnchorBase::mask;
    std::sort(sortBegin, sortEnd,
    [& mask](AnchorBase::AnchorType & a, 
                AnchorBase::AnchorType & b)
    {
        return (a & mask) < (b & mask);
    }) ;
}

inline void Anchors::appendValue(Anchors::AnchorType val)
{
    set[len++]=val;
}

inline void Anchors::appendValue(Anchors::AnchorType val1, Anchors::AnchorType val2)
{
    setAnchor(len++, val1, val2);
    //std::cout << "len " << len << std::endl;
}

void MapParm::print()
{
    std::cerr << "blockSize " << blockSize << std::endl
             << "alpha " << alpha << std::endl
             << "delta " << delta << std::endl
             << "threshold " << threshold << std::endl
             << "kmerStep " << kmerStep << std::endl
             << "shapeLen " << shapeLen << std::endl;
    
}
/*
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
    //std::cout << "Time[s]: " << sysTime() - time << std::endl;
    return 0;
}

template <typename TDna, typename TSpec>
void Mapper<TDna, TSpec>::printHits()
{
    std::cout << "Hits: " << lengthSum(res.hits) << " in sum " << std::endl;
    for (unsigned k = 0; k < length(res.hits); k++)
    {
        std::cout << "reads " << k << " : ";
        for (unsigned j = 0; j < length(res.hits[k]); j++)
            std::cout << (((res.hits[k][j] >> 20) + res.hits[k][j] & 0xfffff))<< " " << (res.hits[k][j] & 0xfffff) << ", ";
        std::cout << std::endl;
    }
        //std::cout << length(res.hits[k]) << std::endl;
}

template <typename TDna, typename TSpec>
void Mapper<TDna, TSpec>::printResult()
{}

template <typename TDna, typename TSpec>
void Mapper<TDna, TSpec>::printParm()
{
    parm.print();
}
*/

static String<Dna5> _complt = "tgcan";
inline void _compltStr(String<Dna5> & str, String<Dna5> & res)
{
    resize(res, length(str));
    for (unsigned k = 0; k < length(str); k++)
     //   res[k]=_complt[str[k] - 'A'];
        res[k] = _complt[(unsigned)ordValue(str[k])];
}

inline void _compltRvseStr(String<Dna5> & str, String<Dna5> & res)
{
    resize(res, length(str));
    for (unsigned k = 0; k < length(str); k++)
     //   res[k]=_complt[str[k] - 'A'];
    {
        res[k] = _complt[(unsigned)ordValue(str[length(str) - k - 1])];
        //std::cout << (unsigned)ordValue(str[length(str) - k - 1]) << std::endl;
    }
}


#endif

