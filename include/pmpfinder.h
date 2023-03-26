#ifndef LINEAR_HEADER_PMP_FINDER_H
#define LINEAR_HEADER_PMP_FINDER_H
#include <seqan/sequence.h>
#include "cords.h"
#include "base.h"
#include "index_util.h"
using namespace seqan;
/*----------  Parms  ----------*/
struct GetDHitListParms : public Parms
{
    int thd_list_n;
    int thd_best_n;
    String<int> thd_list_ns;
    String<int> thd_best_ns;

    GetDHitListParms();
    void toggle (int i);
};


struct GetIndexMatchAllParms : public Parms
{
    int thd_alpha;
    int thd_delta;
    String<int> thd_alphas;

    GetIndexMatchAllParms();
    void toggle (int i);
};

struct GetHIndexMatchAllParms : public GetIndexMatchAllParms {};
struct GetDIndexMatchAllParms : public GetIndexMatchAllParms {};
struct GetSIndexMatchAllParms : public GetIndexMatchAllParms {};

struct ChainAnchorsHitsParms : public Parms
{
    int thd_best_n;
    int thd_drop_score; 
    int thd_min_chain_len;
    uint64_t thd_chain_depth;
    uint64_t thd_chain_dx_depth;
    float thd_stop_chain_len_ratio;

    ChainAnchorsHitsParms();
};

struct ApxParms : public Parms
{
    float thd_sen;

    ApxParms();
};

struct PMPParms : public Parms
{
    GetDHitListParms       pm_gdl;
    GetHIndexMatchAllParms pm_ghima;
    GetDIndexMatchAllParms pm_gdima;
    GetSIndexMatchAllParms pm_gsima;
    ChainAnchorsHitsParms  pm_cah;
    ApxParms               pm_apx;

    CharString             read_id; 
    PMPParms();
    void toggle (int i);
};

//NOTE:Length of read < 1M;
typedef std::array<int, 3> int96;
typedef int96 FeatureType;

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
  unsigned getAbortScore();
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
                int overlap_size,
                int gap_size,
                unsigned thd_accept_score);

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
                       MapParms & mapParm,
                       String<uint64_t> & hit);
uint64_t mnMapReadList(IndexDynamic & index,
                       String<Dna5> & read,
                       Anchors & anchors,
                       MapParms & mapParm,
                       String<uint64_t> & hit);
                       */
 
uint64_t apxMap (IndexDynamic & index,
                 String<Dna5> & read,
                 CharString & read_id,
                 Anchors & anchors,
                 String<uint64_t> & hit, 
                 StringSet<FeaturesDynamic> & f1,
                 StringSet<FeaturesDynamic> & f2,
                 String<UPair> & apx_gaps,
                 String<uint64_t> & cords_str, 
                 String<uint64_t> & cords_end, 
                 String<CordInfo> & cords_info,
                 int f_chain,
                 GlobalParms & pm_g,
                 PMPParms & pm_pmp);
/*
uint64_t filterGenomes (IndexDynamic & index,
                 String<Dna5> & read,
                 Anchors & anchors,
                 MapParms & mapParm,
                 String<uint64_t> & hit, 
                 StringSet<FeaturesDynamic> & f1,
                 StringSet<FeaturesDynamic> & f2,
                 String<UPair> & apx_gaps,
                 String<uint64_t> & cords, 
                 int f_chain);
                 */
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
/*
uint64_t filterAnchorsList(
    String<uint64_t> & anchors, 
    String<std::pair<unsigned, unsigned> > & anchors_list, 
    uint64_t shape_len, 
    uint64_t thd_anchor_accept_density, 
    uint64_t thd_anchor_accept_min, 
    unsigned thd_anchor_err_bit);
    */
#endif
