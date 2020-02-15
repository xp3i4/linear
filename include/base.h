#ifndef LINEAR_HEADER_BASE_H
#define LINEAR_HEADER_BASE_H

#include <seqan/seq_io.h>
#include <seqan/stream.h>
#include <seqan/basic.h>
#include <iostream>

using namespace seqan;

extern const float    base_alpha_;
extern const unsigned base_shape_len_;
extern const unsigned base_block_size_;
extern const unsigned base_delta_; 
extern const unsigned base_threshold_; 
extern const unsigned base_kmer_step_;
extern const uint64_t base_llt_max_;

extern uint64_t const ULLMAX;
extern int64_t const LLMAX;
extern int64_t const LLMIN;
extern unsigned const UMAX;

typedef Iterator <String <Dna5> >::Type TIter5;
typedef unsigned uint;
typedef uint64_t uintll;

struct Options{
    typedef std::string PathType;
    typedef StringSet<PathType> PathsType;
    PathType oPath;
    PathsType r_paths;
    PathsType g_paths;
    
    uint gap_len; //0 to turn off gap mapping module, set > 0 to map gaps whose length > this value
    uint apx_chain_flag;
    uint aln_flag;
    uint sam_flag;
    uint bal_flag;

    unsigned    sensitivity;
    unsigned    thread;
    int         index_t;
    int         feature_t;

//map tuning
    unsigned    listN;
    unsigned    listN2;
    float       alpha;
    float       alpha2;
    float       cordThr;
    float       senThr;
    int         p1;

//global options status
    std::string versions;
    std::string date; 

    Options();
    std::string getOutputPath() const;
}; 

std::pair<uint, uint> 
loadRecords(StringSet<String<Dna5> > & seqs, 
            StringSet<CharString> & ids, 
            Options::PathType path);
int loadRecords(StringSet<String<Dna5> > & seqs, 
                StringSet<CharString> & ids, 
                Options::PathsType & paths);

struct RecordBase
{
    typedef Dna5 DefaultAlphabet;
    typedef CharString RecId;
    typedef String<Dna5> RecSeq; 
};

struct PMRecord
{
    typedef CharString RecId;
    typedef String<Dna5> RecSeq;
    typedef StringSet<RecId> RecIds;
    typedef StringSet<RecSeq> RecSeqs;

    PMRecord(){}
    PMRecord(Options & options);
    
    RecIds id1, id2;
    RecSeqs seq1, seq2; //seq1=read, seq2=ref

    //int loadRecord(Options & options);
};

struct AnchorBase{
    typedef uint64_t AnchorType;
    static const unsigned bit = 20;
    static const uint64_t mask = (1ULL<<bit) - 1;
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
    void sort(Iter begin, Iter end);
    void sortPos2(Iter begin, Iter end);
    void appendValue(AnchorType val);
    AnchorType & operator [](unsigned p);
    Iter begin(); 
    Iter end();
    unsigned length();
};

struct ResBase{
    typedef unsigned SeqLen;
    typedef uint64_t SeqId;
    typedef uint64_t MapPos;
    typedef uint64_t MapScore;
    typedef uint64_t HitType;
    typedef bool MapStrand;

    static const unsigned bit = 32;
    static const uint64_t mask = (1ULL << bit) - 1;
    static const unsigned hitBit = 20;
    static const unsigned hitMask = (1ULL << hitBit) - 1;
    //static const uint64_t hitStrandFlag = 0x8ffffffff;
    //static const uint64_t hitEndFlag = 0x4ffffffff;
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

 uint64_t _nStrand(uint64_t strand);
 uint64_t _flipCoord (uint64_t coord, uint64_t len, uint64_t strand);

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
};

int readRecords_block (StringSet<CharString> & ids, StringSet<String<Dna5> > & reads, String<size_t> & lens, SeqFileIn & fin, int blockSize);
void _compltStr(String<Dna5> & str, String<Dna5> & res);
void _compltRvseStr(String<Dna5> & str, String<Dna5> & res);

struct Dout //debug cout utility
{
    std::ostringstream buffer;
    Dout & operator << (int);
    Dout & operator << (unsigned);
    Dout & operator << (int64_t);
    Dout & operator << (uint64_t);
    Dout & operator << (CharString);
    Dout & operator << (String<int64_t> &);
    Dout & operator << (double);
};
extern Dout dout;

class ostreamWapper
{
    CharString contents; 
public:
    void print_message(std::string strs, size_t start, int end_type, std::ostream & os);
    void print_message(double data, size_t start, int end_type, std::ostream & os);
    void print_message(unsigned data, size_t start, int end_type, std::ostream & os);
};
extern ostreamWapper serr;
/*
class status
{
    String<float> percents;
    String<double> times;
public :
    float getPercents();v
    double getTimes();
    int registPercents(); //return numeric iterator;
    int registTime();
}
*/

std::string _2stdstring (CharString str);

struct CmpInt64
{
    int64_t * p_rslt;
    CmpInt64 & init (int64_t & rslt, int64_t init_val);
    CmpInt64 & min (int64_t & rslt, int64_t val = (~(1LL << 63)));
    CmpInt64 & max (int64_t & rslt, int64_t val = ((1LL << 63) + 1));

    CmpInt64 & operator << (int64_t); //get min of all the suffix 
    CmpInt64 & operator >> (int64_t); //get max of ...
};

void sort_ska(Iterator<String<uint64_t> >::Type it_str, Iterator<String<uint64_t> >::Type it_end);

int erf_1k(int64_t val);

int print_seq(String<Dna5> & seq, uint64_t str, uint64_t end, std::string header);



#endif