#ifndef SEQAN_HEADER_CHAIN_MAP_H
#define SEQAN_HEADER_CHAIN_MAP_H
#include "index_util.h"
using namespace seqan;
uint64_t mnMapReadList( LIndex  & index,
                        String<Dna5> & read,
                        Anchors & anchors,
                        MapParm & mapParm,
                        String<uint64_t> & hit);
uint64_t mnMapReadList( IndexDynamic  & index,
                        String<Dna5> & read,
                        Anchors & anchors,
                        MapParm & mapParm,
                        String<uint64_t> & hit);
#endif