#ifndef LINEAR_HEADER_PMP_FINDER_H
#define LINEAR_HEADER_PMP_FINDER_H
#include <seqan/sequence.h>
#include "index_util.h"
#include "base.h"
using namespace seqan;

//NOTE:Length of read < 1M;
typedef std::array<int, 3> int96;
typedef int96 FeatureType;
typedef std::pair<uint64_t, uint64_t> UPair;

extern unsigned window_size; //_apx_parm_base.window_size

extern int const typeFeatures1_16;
extern int const typeFeatures1_32;
extern int const typeFeatures2_48;

struct ApxMapParmBase
{
    //common window parm
    float band_width;
    unsigned cell_size;
    unsigned cell_num;
    unsigned windowThreshold;
    unsigned windowThresholdReject;
    unsigned windowSize;  //cell_size * cell_num
    unsigned windowDelta; //
    unsigned sup;
    unsigned med;
    unsigned inf; 
    //--script common
    unsigned scpt_step;
    unsigned scpt_bit;
    unsigned scpt_size;
    unsigned scpt_num;
    unsigned scpt_int_step;
    unsigned abort_score;
    ApxMapParmBase (float, 
                    unsigned, unsigned, unsigned,
                    unsigned, unsigned, unsigned,
                    unsigned, unsigned);
};

struct ApxMapParm1_16 : ApxMapParmBase
{
    unsigned scpt_len;
    unsigned scpt_len2;
    int scriptMask;
    int scriptMask2;
    int scptCount[5];
    ApxMapParm1_16 ();
};

struct ApxMapParm1_32 : ApxMapParmBase //32mer script
{
  //script parm
  //Do not modify variable independently
    unsigned scpt_len;
    unsigned scpt_len2;
    short scriptMask;
    int scriptMask2;
    int scptCount[5];
    ApxMapParm1_32 ();
};

struct ApxMapParm2_48 : ApxMapParmBase
{
    ApxMapParm2_48();
};

struct FeaturesDynamic
{
  int fs_type; //features type
  String<short> fs1_16;
  String<short> fs1_32;
  String<int96> fs2_48;
  ApxMapParm1_16 * apx_parm1_16;
  ApxMapParm1_32 * apx_parm1_32;
  ApxMapParm2_48 * apx_parm2_48;

  int isFs1_16();
  int isFs1_32();
  int isFs2_48();
  void setFs1_16();
  void setFs1_32();
  void setFs2_48();
  void setFeatureType(int);
  int init(int type);
  unsigned length();
  FeaturesDynamic();
  FeaturesDynamic(int type);
};

unsigned getWindowThreshold(FeaturesDynamic &);
unsigned getWindowThreshold(StringSet<FeaturesDynamic> &);
unsigned getWindowThresholdReject(FeaturesDynamic &);
unsigned getWindowThresholdReject(StringSet<FeaturesDynamic> &);
unsigned getFeatureWindowSize(FeaturesDynamic & fs);
unsigned getFeatureWindowSize(StringSet<FeaturesDynamic> & fss);

int printScript(FeatureType & val, CharString);

//A wrapper that is(only) used in the gap.cpp
//Do not call this function frequently since the condition branch will drain the performance.
unsigned _windowDist(FeaturesDynamic & f1,
                     FeaturesDynamic & f2,
                     uint64_t x1, uint64_t x2);

bool path_dst(typename Iterator<String<uint64_t> >::Type, 
              typename Iterator<String<uint64_t> >::Type, 
              StringSet<FeaturesDynamic> &,
              StringSet<FeaturesDynamic> &, 
              String<uint64_t> &,
              float const & );

int extendPatch(StringSet<FeaturesDynamic> & f1, 
                StringSet<FeaturesDynamic> & f2, 
                String<uint64_t> & cords,
                int k,
                uint64_t cord1,
                uint64_t cord2,
                int revscomp_const,
                int overlap_size = window_size,
                int gap_size = window_size);

void printInt96(int96 val, CharString header);

int createFeatures(TIter5, TIter5, FeaturesDynamic & ); //serial
int createFeatures(TIter5, TIter5, FeaturesDynamic &, unsigned); //parallel
//@int feature_type, @unsigned threads
int createFeatures(StringSet<String<Dna5> > &, 
                   StringSet<FeaturesDynamic > &, int, unsigned); //parallel
int createFeatures(StringSet<String<Dna5> > &, 
                   StringSet<FeaturesDynamic > &, int); //serial
/*
uint64_t mnMapReadList(LIndex  & index,
                       String<Dna5> & read,
                       Anchors & anchors,
                       MapParm & mapParm,
                       String<uint64_t> & hit);
uint64_t mnMapReadList(IndexDynamic & index,
                       String<Dna5> & read,
                       Anchors & anchors,
                       MapParm & mapParm,
                       String<uint64_t> & hit);
                       */

//Chainning Score metric wrapper: including a score function with corresponding parms.
struct ChainScoreMetric 
{
    int thd_abort_score; //lower bound of average chain score per anchor

    int (*getScore) (uint64_t const &, uint64_t const &);
    int (*getScore2)(uint64_t const &, uint64_t const &, uint64_t const &, uint64_t const &, uint64_t const & read_len);
    int getAbortScore();
    ChainScoreMetric();
    ChainScoreMetric(int abort_score, int (*scoreFunc) (uint64_t const &, uint64_t const &));
    ChainScoreMetric(int abort_score, int (*scoreFunc)(uint64_t const &, uint64_t const &, uint64_t const &, uint64_t const &, uint64_t const & read_len));
};

struct ChainsRecord
{
    int score; //chain score 
    int score2; //copy of score 
    int len;
    int p2anchor;
    int root_ptr;
    int f_leaf;

    int isLeaf();
};

int getBestChains(String<uint64_t> & anchor, String<ChainsRecord> & chains,
                  int (*getScore) (uint64_t const &, uint64_t const &));

int chainAnchorsBase(String<uint64_t> &, StringSet<String<uint64_t> > &, String<int> &, uint, uint, uint, uint64_t, ChainScoreMetric &, uint64_t (*get_anchor_x)(uint64_t));
int getForwarChainDxDy(uint64_t const & cord11, uint64_t const & cord12, uint64_t const & cord21, uint64_t const & cord22, uint64_t const & read_len, int64_t & dx, int64_t & dy);
int getApxChainScore3(uint64_t const & cord11, uint64_t const & cord12, uint64_t const & cord21, uint64_t const & cord22, uint64_t const & read_len);
int chainBlocksCords(String<uint64_t> & cords, String<UPair> & str_ends_p, ChainScoreMetric & chn_score,
   uint64_t read_len, uint thd_init_cord_score, uint64_t thd_major_limit, 
   void (*unsetEndFunc)(uint64_t &), void (*setEndFunc)(uint64_t &), int f_header);
 
uint64_t apxMap (IndexDynamic & index,
                 String<Dna5> & read,
                 Anchors & anchors,
                 MapParm & mapParm,
                 String<uint64_t> & hit, 
                 StringSet<FeaturesDynamic> & f1,
                 StringSet<FeaturesDynamic> & f2,
                 String<UPair> & apx_gaps,
                 String<uint64_t> & cords, 
                 float cordLenThr,
                 int f_chain);
uint64_t filterGenomes (IndexDynamic & index,
                 String<Dna5> & read,
                 Anchors & anchors,
                 MapParm & mapParm,
                 String<uint64_t> & hit, 
                 StringSet<FeaturesDynamic> & f1,
                 StringSet<FeaturesDynamic> & f2,
                 String<UPair> & apx_gaps,
                 String<uint64_t> & cords, 
                 float cordLenThr,
                 int f_chain);
int gather_blocks_ (String<uint64_t> & cords, 
                    String<UPair> & str_ends,   //result [] closed 
                    String<UPair> & str_ends_p, //result pointer [,) right open
                    uint64_t str_,
                    uint64_t end_,
                    uint64_t readLen,
                    uint64_t thd_large_gap,
                    uint64_t thd_cord_size,
                    int f_set_end,  //will set block end sgn at end of evey detected new noncontinous block
                    uint64_t (*isEndFunc)(uint64_t),
                    void (*setEndFunc)(uint64_t &)); 
int gather_gaps_y_ (String<uint64_t> & cords, 
                    String<UPair> & str_ends,
                    String<UPair> & gaps,
                    uint64_t readLen,
                    uint64_t thd_gap_size);

int preFilterChains2(String<uint64_t> & hits,  String<UPair> & str_ends_p, void (*setEndFunc)(uint64_t &));
UPair getUPForwardy(UPair str_end, uint64_t read_len);

#endif
