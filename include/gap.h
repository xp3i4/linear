#ifndef LINEAR_HEADER_GAP_H
#define LINEAR_HEADER_GAP_H
#include <seqan/sequence.h>
#include "pmpfinder.h"

using namespace seqan;
//ClipRecords : operator struct cord

extern uint64_t EmptyClipConst;

uint64_t getClipStr(String<uint64_t> &, int); //seq1
uint64_t getClipEnd(String<uint64_t> &, int);
int getClipsLen(String<uint64_t> &);
void insertClipStr(String<uint64_t> &, uint64_t);
void insertClipEnd(String<uint64_t> &, uint64_t);
bool isClipEmpty(uint64_t);



/**
 * Re-map gaps in cords.
 * Gaps at the beginning and the end of the read are also included.
 */
int mapGaps(StringSet<String<Dna5> > & seqs, 
            String<Dna5> & read, 
            String<Dna5> & comstr,
            String<uint64_t> & cords_str, 
            String<uint64_t> & cords_end, 
            String<uint64_t> & g_hs,
            String<uint64_t> & g_anchor,
            String<uint64_t> & clips,
            StringSet<FeaturesDynamic> & f1,
            StringSet<FeaturesDynamic>& f2,
            int thd_gap, 
            int thd_tileSize,
            float thd_err_rate
           );

int print_clips_gvf_(StringSet<String<uint64_t> > & clips, 
              StringSet<CharString> & readsId, 
              StringSet<CharString> & genomesId,
              std::ofstream & of);

#endif