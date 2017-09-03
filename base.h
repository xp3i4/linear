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

using namespace seqan;

//===================================================================
// variable and type def
//===================================================================
//template <typename TSpec = void>
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
const unsigned Const_::_BLOCKSIZE = 1000;
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
    AnchorType setAnchor(unsigned p, AnchorType pos1,  AnchorType pos2);
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
    typedef bool MapStrand;

    static const unsigned bit = 32;
    static const Const_::BIT_INT_ mask = (1ULL << bit) - 1;

};

struct PMRes
{
    typedef typename ResBase::SeqId Id;
    typedef typename ResBase::MapPos Pos;
    typedef typename ResBase::MapScore Score;
    typedef typename ResBase::MapStrand Strand;

//    String<Id>    id;
//    String<Pos>   pos;
    StringSet<String<Pos> > pos;
    StringSet<String<Score> > score;  
    StringSet<String<Strand> > strand;

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
        alpha(Const_::_ALPHA),
        delta(Const_::_DELAT),
        threshold(Const_::_THRESHOLD),
        kmerStep(Const_::_KMERSTEP),
        shapeLen(Const_::_SHAPELEN)
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

    Mapper();
    Mapper(Options & options);
    Seqs & reads(){return record.seq1;};
    Seqs & genomes(){return record.seq2;};
    Parm & mapParm(){return parm;};
    Res & result(){return res;};
    Index & index(){return qIndex;};
    void printResult();    
    void printParm();
    int createIndex();
     
    //Mapper(Options const & options)
    //{
    //    loadRecords(options);
    //    setMapParm(options);
    //};
};

int Options::print()
{
    
    std::cerr << "kmerLen " << kmerLen << std::endl
              << "MiKmLen " << MiKmLen << std::endl
              << "reads path " << rPath << std::endl
              << "genomes Path " << gPath << std::endl
              << "output path " << oPath << std::endl
              << "Sensitive " << Sensitive << std::endl;
}


template <typename TDna>
int PMRecord<TDna>::loadRecord(Options & options)
{
    std::cerr << "loading sequences from files... " << std::endl;
    SeqFileIn rFile(toCString(options.rPath));
    SeqFileIn gFile(toCString(options.gPath));
    //try{
        
        seqan::readRecords(id1, seq1, rFile);
            
        seqan::readRecords(id2, seq2, gFile);
    //   throw 
    //}
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
    len = 0;
}

inline void Anchors::init()
{
    init(AnchorType(0), size());
}

inline Anchors::AnchorType Anchors::setAnchor(unsigned p, 
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
void Mapper<TDna, TSpec>::printResult()
{}

template <typename TDna, typename TSpec>
void Mapper<TDna, TSpec>::printParm()
{
    parm.print();
}

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
        res[k] = _complt[(unsigned)ordValue(str[length(str) - k])];
}


template <typename TDna, typename TSpec>
inline unsigned getIndexMatch(typename PMCore<TDna, TSpec>::Index  & index,
                              typename PMRecord<TDna>::RecSeq & read,
                              Anchors & anchor,
                              MapParm & mapParm
                             )
{    
    uint64_t dn, pre;
    unsigned block = (mapParm.blockSize < length(read))?mapParm.blockSize:length(read);
    unsigned dt = block * (mapParm.alpha / (1 - mapParm.alpha));
    //std::cout << read << std::endl;

    hashInit(index.shape, begin(read));
    anchor.init();
    for (unsigned h=0; h <= length(read) - block; h += dt)
    //for (unsigned h=0; h <= length(read) - block; h += block)
    {
        hashInit(index.shape, begin(read) + h);
        for (unsigned k = h; k < h + block; k++)
        {
            hashNext(index.shape, begin(read) + k);
        //std::cout << "shape " << index.shape.XValue << std::endl;
            dn = getDir(index, index.shape);
            pre = ~0;
            //std::cout << "dn " << dn << std::endl;
            if(_getBodyCounth(index.dir[dn+1]) - _getBodyCounth(index.dir[dn]) != 0)
    //        std::cout << h << " " << _getBodyCounth(index.dir[dn+1]) - _getBodyCounth(index.dir[dn]) << std::endl;
            if(_getBodyCounth(index.dir[dn+1]) - _getBodyCounth(index.dir[dn]) < mapParm.delta)
            {
                uint64_t countn = _getBodyCounth(index.dir[dn]);
               //std::cout << countn << " countn " << std::endl;
                while ( countn < _getBodyCounth(index.dir[dn + 1]))
                {
      //              std::cout << index.sa[countn] - pre << std::endl;
                    if (index.sa[countn] - pre > mapParm.kmerStep)
                    {
                        //#anchor[x++] = (((index.sa[n])- k) << 20) + k;
    //std::cout << countn << " countn " << index.sa[countn] - k << " " << k << std::endl;
     //                   std::cout << index.sa[countn]- k << " index.sa - k" << std::endl;
                        anchor.appendValue(index.sa[countn]- k, k);
                        pre = index.sa[countn];
                    }
                    countn++;
                }
            }
        }
    }
    return 0;
}

inline unsigned getAnchorMatch(Anchors & anchors, MapParm & mapParm)
{
    uint64_t maxLen = 0, c_b=0, ak, cbb=0, sb=0, start = 0;
    anchors[0] = anchors[1];
    ak=anchors[0];
    anchors.sort(anchors.begin(), anchors.end());
    //std::cout << anchors.length() << " length ";
    for (uint64_t k = 1; k <= anchors.length(); k++)
    {
        //std::cout << anchors[k] << std::endl;
        //std::cout << "getAnchors " << anchors.getPos1(k) << " " << anchors.getPos2(k) << std::endl;
        //std::cout << anchors[k]- ak << " " << std::endl;
        if (anchors[k] - ak < AnchorBase::AnchorValue)
            cbb++;
        else
        {
            anchors.sortPos2(anchors.begin() + sb, anchors.begin() + k);
            //std::cout << sb << " " << k << " sb" << std::endl;
            for (uint64_t m = sb+1; m < k; m++)
            {
                if(anchors.deltaPos2(m, m-1) >  mapParm.shapeLen)
                    c_b += mapParm.shapeLen;
                else
                    c_b += anchors.deltaPos2(m, m-1); 
            //std::cout << anchors.deltaPos2(m, m-1) << " shapelen " << mapParm.shapeLen << std::endl;
            }
            if (c_b > maxLen)
            {
                maxLen = c_b;
                start = sb;
            }
            sb = k;
            ak = anchors[k];
            cbb = 1;
            c_b = mapParm.shapeLen;
        }

    }
    //std::cout << "start maxlen " << anchors.getPos1(start) << " " << anchors.getPos2(start) << " " << maxLen << std::endl;
    return start + (maxLen << 20) ;
}

template <typename TDna, typename TSpec>
inline unsigned mnMapRead(typename PMCore<TDna, TSpec>::Index  & index,
                          typename PMRecord<TDna>::RecSeq & read,
                          Anchors & anchors,
                          MapParm & mapParm
                         )
{
    
    getIndexMatch<TDna, TSpec>(index, read, anchors, mapParm);    
    return getAnchorMatch(anchors, mapParm);
}

template <typename TDna, typename TSpec>
void mnMap(typename PMCore<TDna, TSpec>::Index   & index,
           typename PMRecord<TDna>::RecSeqs      & reads,
           MapParm                      & mapParm)
      //      Anchors & anchors)
           //PMResSet             & rs )
{
    typedef typename Mapper<TDna, TSpec>::Seq Seq;
    typedef typename Mapper<TDna, TSpec>::Anchors Anchors;
    double time=sysTime();
    unsigned mask1 = (1<<20)-1;
    //#unsigned res = 0, tmp = 0;
    //#uint64_t anchor[mask + 1] = {1ULL<<63};

    std::cerr << "Filtering reads\n"; 
    Anchors anchors(Const_::_LLTMax, AnchorBase::size);
    //Anchors anchor;
    Seq comStr;
    
    for (unsigned j = 0; j< length(reads); j++)
    {
        unsigned res = mnMapRead<TDna, TSpec>(index, reads[j], anchors, mapParm);
        //std::cout << (res & mask1) << std::endl;
        if (res < (mapParm.threshold << 20))
        {
            _compltRvseStr(reads[j], comStr);
            unsigned tmp = mnMapRead<TDna, TSpec>(index, comStr, anchors, mapParm);
            res = (tmp > (mapParm.threshold << 20))?tmp:~0;
        }
        if (res != ~0 )
          std::cout << j << " " << anchors.getPos1(res & mask1) << " " << anchors.getPos2(res & mask1) << " " << std::endl;//<< " " << maxLen << std::endl;
        else
            std::cout << j <<" 0 0 " << std::endl;
        //std::cout << (res & mask1) << std::endl;
//        else
//            rs.appendValue(j, anchors[res ], );
        ///rs.appendRes(res);
//            std::cout << j << " " << (res & mask1) << " " << (res >> 20) << std::endl;
        //#assignValueI1(rs[j], _getSA_i2(anchor[res & mask1] >> 20));
        //#assignValueI2(rs[j], _getSA_i2(anchor[res & mask1] >> 20) + length(reads[j]));
//========
//rs
//========
        //rs.appendValue(_getSA_i2(anchor[res & mask1].getPos1()), _getSA_i2(anchor[res & mask1].getPos1()) + length);
    }
    std::cerr << "Time[s]: " << sysTime() - time << std::endl;
}

template <typename TDna, typename TSpec>
void map(Mapper<TDna, TSpec> mapper)
{
//    map.printParm();
    _DefaultMapParm.print();
    double time = sysTime();
    std::cerr << "map" << std::endl;
    mapper.createIndex();
    mnMap<TDna, TSpec>(mapper.index(), mapper.reads(), _DefaultMapParm);//, mapper.result());
    std::cerr << "Total time [s] " << sysTime() - time << std::endl;
    mapper.printResult();
    
}

#endif
