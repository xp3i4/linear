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
                    unsigned);
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
    int thd_abort_score;

    int (*getScore) (uint64_t const &, uint64_t const &);
    int getAbortScore();
    ChainScoreMetric();
    ChainScoreMetric(int abort_score, int (*scoreFunc) (uint64_t const &, uint64_t const &));
};

struct ChainsRecord
{
    int score; //chain score = sum of score2 of all records in the score
    int score2; // this record score
    int len;
    int p2anchor;
};
int getBestChains(String<uint64_t> & anchor,
                  String<ChainsRecord> & chains,
                  int anchor_end,
                  int (*getScore) (uint64_t const &, uint64_t const &));

int createChainsFromAnchors(StringSet<String<uint64_t> > &, String<ChainsRecord> &, String<uint64_t> &, int, ChainScoreMetric &);
 
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
                    uint64_t readLen,
                    uint64_t thd_large_gap,
                    uint64_t thd_cord_size,
                    int f_set_end); //will set block end sgn at end of evey detected new noncontinous block
int gather_gaps_y_ (String<uint64_t> & cords, 
                    String<UPair> & str_ends,
                    String<UPair> & gaps,
                    uint64_t readLen,
                    uint64_t thd_gap_size);
#endif
