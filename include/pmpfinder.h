#ifndef LINEAR_HEADER_PMP_FINDER_H
#define LINEAR_HEADER_PMP_FINDER_H
#include <seqan/sequence.h>
#include "index_util.h"
#include "base.h"
using namespace seqan;

//NOTE:Length of read < 1M;
typedef std::array<int, 3> int96;
typedef int96 FeatureType;

extern const unsigned window_size; //16*12

extern int const typeFeatures1_32;
extern int const typeFeatures2_48;
struct FeaturesDynamic
{
  int fs_type; //features type
  String<short> fs1_32;
  String<int96> fs2_48;

  int isFs1_32();
  int isFs2_48();
  void setFs1_32();
  void setFs2_48();
  void setFeatureType(int);
  FeaturesDynamic(int type = typeFeatures2_48);
};

unsigned get_windowThreshold(FeaturesDynamic &);
unsigned get_windowThreshold(StringSet<FeaturesDynamic> &);

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
uint64_t apxMap (IndexDynamic & index,
                 String<Dna5> & read,
                 Anchors & anchors,
                 MapParm & mapParm,
                 String<uint64_t> & hit, 
                 StringSet<FeaturesDynamic> & f1,
                 StringSet<FeaturesDynamic> & f2,
                 String<uint64_t> & cords, 
                 float cordLenThr);
#endif
