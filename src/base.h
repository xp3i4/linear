#ifndef LINEAR_HEADER_F_IO_H
#define LINEAR_HEADER_F_IO_H

using namespace seqan;

const float base_alpha_ = 0.75;
const unsigned base_shape_len_ = 25;
const unsigned base_shape_weight_ = 17;
const unsigned base_block_size_ = 100;
const unsigned base_delta_ = 32; 
const unsigned base_threshold_= 30; 
const unsigned base_kmer_step_ = 1000;
const uint64_t base_llt_max_ = ~0;
     
struct Options{
    unsigned    kmerLen;
    unsigned    MiKmLen;
    typename    std::string rPath;
    typename    std::string gPath;
    typename    std::string oPath;
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
    int p1;

    Options();
    std::string getGenomePath(); 
    std::string getReadPath();
    std::string getOutputPath();
    int print();
}; 
struct RecordBase
{
    typedef Dna5 DefaultAlphabet;
    typedef CharString RecId;
    typedef String<TDna> RecSeq; 
};

struct PMRecord
{
    typedef CharString RecId;
    typedef String<Dna5> RecSeq;
    typedef StringSet<RecId> RecIds;
    typedef StringSet<RecSeq> RecSeqs;

    PMRecord(){}
    PMRecord(Options & options);
    
    CharString readPath;
    CharString genomePath; 
    RecIds id1, id2;
    RecSeqs seq1, seq2; //seq1=read, seq2=ref

    int loadRecord(Options & options);
};


struct ResBase{
    typedef unsigned SeqLen;
    typedef typename uint64_t SeqId;
    typedef typename uint64_t MapPos;
    typedef typename uint64_t MapScore;
    typedef typename uint64_t HitType;
    typedef bool MapStrand;

    static const unsigned bit = 32;
    static const uint64_t mask = (1ULL << bit) - 1;
    static const unsigned hitBit = AnchorBase::bit;
    static const unsigned hitMask = AnchorBase::mask;
    //static const uint64_t hitStrandFlag = 0x8ffffffff;
    //static const uint64_t hitEndFlag = 0x4ffffffff;
};

struct AnchorBase{
    typedef uint64_t AnchorType;
    const unsigned bit = 20;
    const uint64_t mask = (1ULL<<bit) - 1;
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
    unsigned size() const;
    AnchorType & operator [](unsigned p);
    Iter begin(); 
    Iter end();
    unsigned length();
};

struct PMCore
{
    typedef HIndex<base_shape_len_> Index;
    typedef StringSet<StringDna5> Seqs;
    typedef Anchors Anchors;
};

struct PMRes
{
    CharString Id;
    typedef String<uint64_t> HitString;
    typedef StringSet<HitString> HitSet;

    StringSet<String<uint64_t> > pos;
    StringSet<String<uint64_t> > score;  
    StringSet<String<uint64_t> > hits;
};

inline uint64_t _nStrand(uint64_t strand);
inline uint64_t _flipCoord (uint64_t coord, uint64_t len, uint64_t strand)

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
      
    MapParm();
    MapParm(unsigned bs, unsigned dt, unsigned thr, 
            unsigned ks, unsigned sl, unsigned st,
            unsigned ad, unsigned mr, unsigned listn,
            unsigned listn2,
            float ap, float ap2, float alt, float rt, float ct, float sent, float clst);
    MapParm(MapParm & parm);
    void setMapParm(Options & options);
    void print ();
}_DefaultMapParm;

struct MapperBase
{
    typedef Dna5 DefaultAlphabet;
    typedef Minimizer<base_shape_len_> DefaultShape;
    typedef PMRecord MRecord;
    typedef PMRes    MRes;
    typedef MapParm   MParm;
    typedef PMCore    MCore;
    typedef typename PMCore::Index MIndex;
    typedef typename PMCore::Anchors MAnchors;
    typedef typename PMRecord::RecSeq MSeq; 
    typedef typename PMRecord::RecSeqs MSeqs;
};

static const String<Dna5> _complt = "tgcan";
inline void _compltStr(String<Dna5> & str, String<Dna5> & res)

seqan::ArgumentParser::ParseResult
parseCommandLine(Options & options, int argc, char const ** argv);
#endif