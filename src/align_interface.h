#ifndef LINEAR_HEADER_ALIGN_INT_H
#define LINEAR_HEADER_ALIGN_INT_H
#include <seqan/align.h>
#include "pmpfinder.h"
#include "f_io.h"

using namespace seqan;
class GapRecords
{
    typedef Row<Align<String<Dna5>, ArrayGaps> >::Type TRow;
    typedef std::pair<uint64_t, uint64_t> TCPair;
    typedef std::pair<TRow, TRow> TRPair;
    typedef String<std::pair<uint64_t, uint64_t> > TCPairs;
    typedef String<std::pair<TRow, TRow> > TRPairs;

public:
    TCPairs c_pairs;  //pair of start and end coordinate of gaps 
    TRPairs r_pairs;  //pair of alingments of blocks contains the start and end of gaps 
    String<int> bam_segs_id;
    uint64_t dx, dy;  //shift from start and end of gaps to the start and end of of the alignment

    TCPair & get_c_pair(int i);  
    TRPair & get_r1_pair(int i); //get rows of first(start) cord of i_th gaps
    TRPair & get_r2_pair(int i); //get rows of second(end) cord of i_th gaps
    int getBamSegIdHead (int i);
    int getBamSegIdTail (int i);
    uint64_t getJointHeadCord(int i); 
    uint64_t getJointTailCord(int i);
    int clear_();
};

struct GapParm
{
    int thd_clip_score;
    int thd_reject_score;
    int thd_accept_score;
    float thd_accept_density;
    Score<int, Simple> thd_clip_scheme;
    GapParm ();
} _gap_parm;

int align_cords (StringSet<String<Dna5> >& genomes,
                 String<Dna5> & read, 
                 String<Dna5> & comrevRead,
                 String<uint64_t> & cords,
                 String<BamAlignmentRecordLink> & bam_records,
                 int p,
                 int block_size = window_size,
                 int band = window_size / 2
                );

#endif