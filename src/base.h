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
#include <seqan/basic.h>
#include <seqan/arg_parse.h>
#include <seqan/vcf_io.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <bitset>
#include <climits>
#include <thread>
#include <chrono>
#include <atomic>   
#include <iomanip>
#include <functional>   // for std::ref()
#include <chrono>
#include <thread>
#include <ctime>
   
#include "shape_extend.h"
#include "index_extend.h"

using namespace seqan;

//===================================================================
// variable and type def
//===================================================================
//template <typename TSpec = void>

/*
struct Status
{
    seqCount;
    String <int> step; 
    String <double> time;
    Status ()
    {
        
    };
}status;
*/


struct Const_{
    
    typedef uint64_t    BIT_INT_;
    typedef CharString  PATH_;
    typedef seqan::Dna5 DEFAULT_ALPHABET_ ;
    typedef CharString  ID_;

    static const unsigned _SHAPELEN;
    static const unsigned _SHAPEWHT;
    static const unsigned _BLOCKSIZE;
    static const unsigned _DELTA; 
    static const unsigned _THRESHOLD; 
    static const float    _ALPHA ;
    static const unsigned _KMERSTEP;
    static const uint64_t _LLTMax;
        
};

const float Const_::_ALPHA = 0.75;
const unsigned Const_::_SHAPELEN = 25;
const unsigned Const_::_SHAPEWHT = 17;
const unsigned Const_::_BLOCKSIZE = 100;
const unsigned Const_::_DELTA = 32; 
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
    unsigned    sensitivity;
    unsigned    thread;
//map tuning
    unsigned listN;
    unsigned listN2;
    float alpha;
    float alpha2;
    float cordThr;
    float senThr;
//
    
    Options():
        kmerLen(Const_::_SHAPELEN),
        MiKmLen(Const_::_SHAPEWHT),
        rPath(""),
        gPath(""),
        oPath("mapper_result.txt"),
        Sensitive(false),
        sensitivity(0),
        thread(4)
        {}
    Const_::PATH_ getGenomePath() const {return gPath;};
    Const_::PATH_ getReadPat() const {return rPath;};
    Const_::PATH_ getOutputPath() const {return oPath;};
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
    
    Const_::PATH_ readPath;
    Const_::PATH_ genomePath; 
    RecIds id1, id2;
    RecSeqs seq1, seq2; //seq1=read, seq2=ref

    int loadRecord(Options & options);
};

struct AnchorBase{
    typedef typename Const_::BIT_INT_ AnchorType; 
    static const unsigned size = 131072 * 1000;
    static const unsigned bit = 20;
    static const uint64_t AnchorValue = 1000ULL << bit;
    static const typename Const_::BIT_INT_ mask = (1ULL<<bit) - 1;
};

struct Anchors{
    
    typedef typename Iterator<String<uint64_t> >::Type Iter;
    typedef AnchorBase::AnchorType AnchorType;
    typedef String<AnchorType> AnchorString;
    
    AnchorString set;

    

    void init(AnchorType val, unsigned k);
    void init(int length);
    void init();
    void setAnchor(unsigned p, AnchorType pos1,  AnchorType pos2);
    AnchorType getPos1(unsigned p) const;
    AnchorType getPos2(unsigned p) const;
    AnchorType deltaPos1(unsigned p1, unsigned p2);
    AnchorType deltaPos2(unsigned p1, unsigned p2);
    void sort(Iter begin, Iter end);
    void sortPos2(Iter begin, Iter end);
    void appendValue(AnchorType val);
    unsigned size() const {return AnchorBase::size;};
    AnchorType & operator [](unsigned p){return set[p];};
    Iter begin(); 
    Iter end();
    unsigned length() {return seqan::length(set);};
};


template <typename TDna = Const_::DEFAULT_ALPHABET_, 
        typename CoreMinimizer = Minimizer<Const_::_SHAPELEN> > 
struct CoreBase{
    typedef typename Const_::DEFAULT_ALPHABET_ DefaultAlphabet;
    typedef Minimizer<Const_::_SHAPELEN> DefaultShape;
    typedef typename PMRecord<TDna>::RecSeqs RecSeqs; 
    typedef Shape<TDna, CoreMinimizer> CoreShape;
    //typedef Index<RecSeqs, IndexQGram<CoreMinimizer, OpenAddressing> > CoreIndex;
//==============
//  change the type of index from generic indx to optimized 25-mer HIndex;
    typedef HIndex<Const_::_SHAPELEN> CoreIndex;
//==============   
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
    //static const uint64_t hitStrandFlag = 0x8ffffffff;
    //static const uint64_t hitEndFlag = 0x4ffffffff;
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
    //void setHitStrand(HitType &);
    //uint64_t getHitStrand(HitType const &) const;
    //void setHitEnd(HitType &);
    //uint64_t getHitEnd(HitType const &) const;


};

/*
 * flip strand from 0, 1 to -1, 1;
 * strand = 0, 1, other values is not allowed
 * return -1 , 1
 */
inline uint64_t _nStrand(uint64_t strand)
{
    return (strand << 1) - 1;
}

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

//inline void PMRes::setHitStrand(typename PMRes::HitType & hit)
//{
//    hit |= ResBase::hitStrandFlag;
//}
//
//inline uint64_t PMRes::getHitStrand(typename PMRes::HitType const & hit) const
//{
//    return hit & ResBase::hitStrandFlag;
//}
//
//inline void PMRes::setHitEnd(typename PMRes::HitType & hit)
//{
//    hit |= ResBase::hitEndFlag;
//}
//
//inline uint64_t PMRes::getHitEnd(typename PMRes::HitType const & hit) const
//{
//    return hit & ResBase::hitEndFlag;
//}

struct MapParm{
    unsigned    blockSize;
    unsigned    delta;
    unsigned    threshold;
    unsigned    kmerStep;
    unsigned    shapeLen;
    unsigned    sensitivity;
    unsigned    anchorDeltaThr;
    unsigned    minReadLen;
    unsigned    listN;
    unsigned    listN2;
    float       alpha;
    float       alpha2;
    float       anchorLenThr;
    float       rcThr;
    float       cordThr;
    float       senThr;
    float       clsThr;
      
    
    MapParm():
        blockSize(Const_::_BLOCKSIZE),
        delta(Const_::_DELTA),
        threshold(Const_::_THRESHOLD),
        kmerStep(Const_::_KMERSTEP),
        shapeLen(Const_::_SHAPELEN),
        sensitivity(0),
        anchorDeltaThr(),
        minReadLen(1000),
        listN(1),
        listN2(1),
        alpha(Const_::_ALPHA),
        alpha2(0.5),
        anchorLenThr(0.02),                  // anchors with lenghth > this parameter is pushed into the queue
        rcThr(0.8),                        // when max anchors in the queue with length < this parameters, reverse complement search will be conducted
        cordThr(0.8),
        senThr(0.8),
        clsThr(0.1)
        {}
        
    MapParm(unsigned bs, unsigned dt, unsigned thr, 
            unsigned ks, unsigned sl, unsigned st,
            unsigned ad, unsigned mr, unsigned listn,
            unsigned listn2,
            float ap, float ap2, float alt, float rt, float ct, float sent, float clst):
        blockSize(bs),
        delta(dt),
        threshold(thr),
        kmerStep(ks),
        shapeLen(sl),
        sensitivity(st),
        anchorDeltaThr(ad),
        minReadLen(mr),
        listN(listn),
        listN2(listn2),
        alpha(ap),
        alpha2(ap2),
        anchorLenThr(alt),                  // anchors with lenghth > this parameter is pushed into the queue
        rcThr(rt),                        // when max anchors in the queue with length < this parameters, reverse complement search will be conducted
        cordThr(ct),
        senThr(sent),
        clsThr(clst)
        {} 


    MapParm(MapParm & parm):
        blockSize(parm.blockSize),
        delta(parm.delta),
        threshold(parm.threshold),
        kmerStep(parm.kmerStep),
        shapeLen(parm.shapeLen),
        sensitivity(parm.sensitivity),
        anchorDeltaThr(),
        minReadLen(parm.minReadLen),
        listN(parm.listN),
        listN2(parm.listN2),
        alpha(parm.alpha),
        alpha2(parm.alpha2),
        anchorLenThr(parm.anchorLenThr),
        rcThr(parm.rcThr),
        cordThr(parm.cordThr),
        senThr(parm.senThr),
        clsThr(parm.clsThr)
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

std::ifstream::pos_type _filesize(const char* filename)
{
    std::ifstream in(filename,  std::ifstream::binary | std::ifstream::ate);
    return in.tellg(); 
}

/*
 *[]::lr
 *
template <typename TDna>
int PMRecord<TDna>::loadRecord(Options & options)
{
    double time = sysTime();
    std::cerr <<"[]::lr_Loading sequences from files \r";
    std::fstream fin (toCString(options.gPath), std::fstream::in);
    std::string buffer;
    float filesize = _filesize (toCString(options.gPath));
    uint64_t sno = 0;
    for (unsigned k = 0; getline(fin, buffer); k++)
    {
        if (buffer.compare(0, 1, ">") == 0)
        {
            sno++;
        }
    }
    SeqFileIn gFile(toCString(options.gPath));
    unsigned k = 0;
    resize (id2, sno);
    resize (seq2, sno);
    unsigned pre = 0;
    float step = 1;
    double rt = sysTime ();
    float len = 0;
    while (!atEnd(gFile))
    {
        //if (k + 1 - pre > step / 100 * sno)
        if (sysTime() - rt > 0.5)
        {
            std::cerr << std::fixed;
            std::cerr << ">>reading genomes " << std::setprecision(1) << (k + 1)/ (float) sno * 100 << "% " << "\r";
            //pre = k + 1;
            rt = sysTime();
        }
        //std::cerr << "reading genomes " << (k + 1) << " of " << sno << "\r";
        readRecord(id2[k], seq2[k], gFile);
        ++k;
    }
    std::cerr << "                                    \r";
    //readRecords(id2, seq2, gFile);
    std::cerr << "--Read genomes 100%                     " << std::endl;
    std::cerr << "  " << length(seq2) << " sequences loaded from " << options.gPath << std::endl;
    std::cerr << "  Elapsed time [s] " << sysTime() - time << std::endl;
    return 0;
}
*/

int readRecords_block (StringSet<CharString> & ids, StringSet<String<Dna5> > & reads, String<int> & lens, SeqFileIn & fin, int blockSize)
{
    int start = length(reads);
    readRecords(ids, reads, fin, blockSize);
    for (unsigned k = 0; k < length(reads) - start; k++)
    {
        appendValue (lens, length(reads[start + k]));
    }
    return 0;
}

/*
 *[]::lr
 */
template <typename TDna>
int PMRecord<TDna>::loadRecord(Options & options)
{
    double time = sysTime();
    SeqFileIn gFile(toCString(options.gPath));
    double fileSize = _filesize (toCString(options.gPath));
    bool flag = false;
    unsigned seqCount = 0;
    double currentFileSize = 0;
    std::fstream fin (toCString(options.gPath), std::fstream::in);
    StringSet<String<char> > dotstatus;
    resize(dotstatus, 3);
    dotstatus[0] = ".   ";
    dotstatus[1] = "..  ";
    dotstatus[2] = "... ";
#pragma omp parallel
{
    #pragma omp sections
    {
        #pragma omp section
        {
            //unsigned preSeqCount = 0;
            String <char> probar;
            float prepercent = 0, percent = 0, showpercent = 0, v = 0.87 ;
            unsigned k = 1;
            while (!flag)
            {
                prepercent = percent;
                percent = currentFileSize / fileSize * 100;
                percent = (percent > 100)?prepercent:percent;
                showpercent += v;
                showpercent = (showpercent > percent)?percent:showpercent;
                
                std::cerr << "                                                            \r";
                if (seqCount > 2)
                {
                    std::cerr << ">>Read genomes" << dotstatus[(k - 1)/10 %3] << "            " << seqCount << "/" << std::setprecision(2) << std::fixed << showpercent << "%\r";
                }
                else
                {
                    std::cerr << ">>Read genomes" << dotstatus[(k - 1)/10 %3] << "\r";
                }
                std::this_thread::sleep_for(std::chrono::milliseconds(100));
                k++;
            }
        }
        #pragma omp section
        {
            PMRecord<TDna>::RecId tmp_id;
            PMRecord<TDna>::RecSeq tmp_seq;
            while (!atEnd(gFile))
            {
                clear (tmp_id);
                clear (tmp_seq);
                readRecord (tmp_id, tmp_seq, gFile);
                currentFileSize += length(tmp_seq);
                appendValue (id2, tmp_id);
                appendValue (seq2, tmp_seq);
                ++seqCount;
            }
            flag = true;
        }
    }
}
    std::cerr << "--Read genomes                "<< length(seq2) <<"/100%                   \n";
    std::cerr << "  File: " << options.gPath ;
    std::cerr << "  Elapsed time [s] " << sysTime() - time << std::endl;
    return 0;
}

template <typename TDna>
PMRecord<TDna>::PMRecord(Options & options)
{
    readPath = options.rPath;
    genomePath = options.gPath;
    loadRecord(options);
}

inline void Anchors::init(AnchorType val, unsigned range)
{
    for (unsigned k = 0; k < range; k++)
        seqan::appendValue(set,val);
}

inline void Anchors::init(int length)
{
    clear(set);
    seqan::appendValue(set, 0);
    (void)length;
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
    return seqan::begin(set);
}

inline Anchors::Iter Anchors::end()
{
    return seqan::end(set);
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
    seqan::appendValue(set,val);
}


void MapParm::print()
{
    std::cerr << "blockSize " << blockSize << std::endl
            << "alpha " << alpha << std::endl
            << "alpha2 " << alpha2 << "\n"
            << "listN " << listN << "\n"
            << "listN2 " << listN2 << "\n"
            << "senThr " << senThr << "\n"
            << "delta " << delta << std::endl
            << "threshold " << threshold << std::endl
            << "kmerStep " << kmerStep << std::endl
            << "shapeLen " << shapeLen << std::endl
            //<<  "sensitivity " << sensitivity << "\n"
            << "anchorDeltaThr " << anchorDeltaThr << "\n"
            << "minReadLen " << minReadLen << "\n"
            << "anchorLenThr" << anchorLenThr << "\n"
            << "rcThr " << rcThr << "\n"
            << "cordThr" << cordThr << "\n";
}

static const String<Dna5> _complt = "tgcan";
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

seqan::ArgumentParser::ParseResult
parseCommandLine(Options & options, int argc, char const ** argv)
{
    // Setup ArgumentParser.
    seqan::ArgumentParser parser("Linear");
    // Set short description, version, and date.
    setShortDescription(parser, "Alignment of SMRT sequencing read");
    setVersion(parser, "1.0");
    setDate(parser, "May 2018");

    // Define usage line and long description.
    addUsageLine(parser,
                    "[\\fIOPTIONS\\fP] \"\\fIread.fa\\fP\" \"\\fIgnome.fa\\fP\"");
    addDescription(parser,
                    "Program for mapping raw SMRT sequencing reads to reference genome.");

    // Argument.
    addArgument(parser, seqan::ArgParseArgument(
        seqan::ArgParseArgument::INPUT_FILE, "read"));
    setHelpText(parser, 0, "Reads file .fa, .fasta");

    addArgument(parser, seqan::ArgParseArgument(
        seqan::ArgParseArgument::INPUT_FILE, "genome", true));
    setHelpText(parser, 1, "Reference file .fa, .fasta");

    addSection(parser, "Mapping Options");
    addOption(parser, seqan::ArgParseOption(
        "o", "output", "choose output file.",
            seqan::ArgParseArgument::STRING, "STR"));
    addOption(parser, seqan::ArgParseOption(
        "s", "sensitivity", "Sensitivity mode. -s 0 normal {DEFAULT} -s 1 fast  -s 2 sensitive",
            seqan::ArgParseArgument::INTEGER, "INT"));
    addOption(parser, seqan::ArgParseOption(
        "t", "thread", "Default -t 4",
            seqan::ArgParseArgument::INTEGER, "INT"));
    
// mapping parameters for tunning 
    addOption(parser, seqan::ArgParseOption(
        "l1", "listn1", "mapping::listn1",
            seqan::ArgParseArgument::INTEGER, "INT"));     
    addOption(parser, seqan::ArgParseOption(
        "l2", "listn2", "mapping::listn2",
            seqan::ArgParseArgument::INTEGER, "INT"));   
    addOption(parser, seqan::ArgParseOption(
        "a1", "alpha1", "mapping::alpha1",
            seqan::ArgParseArgument::DOUBLE, "DOUBLE"));        
    addOption(parser, seqan::ArgParseOption(
        "a2", "alpha2", "mapping::alpha2",
            seqan::ArgParseArgument::DOUBLE, "DOUBLE"));  
    addOption(parser, seqan::ArgParseOption(
        "t1", "cordThr", "mapping::cordThr",
            seqan::ArgParseArgument::DOUBLE, "DOUBLE"));
    addOption(parser, seqan::ArgParseOption(
        "t2", "senThr", "mapping::senThr",
            seqan::ArgParseArgument::DOUBLE, "DOUBLE"));        
        
    // Add Examples Section.
////////////////////////
    addTextSection(parser, "Examples");
    addListItem(parser,
                "\\fBlinear \\fP \\fIreads.fa genomes.fa\\fP",
                "raw map reads.fa to genomes.fa");
    addTextSection(parser, "Examples");
    addListItem(parser,
                "\\fBlinear\\fP \\fP-a \\fIreads.fa genomes.fa\\fP",
                "align reads.fa to genomes.fa");

    // Parse command line.
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

    if (res != seqan::ArgumentParser::PARSE_OK)
        return res;

    getOptionValue(options.oPath, parser, "output");
    getOptionValue(options.sensitivity, parser, "sensitivity");
    getOptionValue(options.thread, parser, "thread");
    
    getOptionValue(options.listN, parser, "listn1");
    getOptionValue(options.listN2, parser, "listn2");
    getOptionValue(options.alpha, parser, "alpha1");
    getOptionValue(options.alpha2, parser, "alpha2");
    getOptionValue(options.cordThr, parser, "cordThr");
    getOptionValue(options.senThr, parser, "senThr");
    
    seqan::getArgumentValue(options.rPath, parser, 0);
    seqan::getArgumentValue(options.gPath, parser, 1);

    return seqan::ArgumentParser::PARSE_OK;

}

#endif

