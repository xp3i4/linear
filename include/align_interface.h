#ifndef LINEAR_HEADER_ALIGN_INT_H
#define LINEAR_HEADER_ALIGN_INT_H
#include <seqan/align.h>
#include "f_io.h"
using namespace seqan;
/*
 * Struct of gap recods containing start, end coordinates and alignment of gaps
 * Each record of c_pairs has two rows in the a_rows;
 * Gaps are defined as pairs of its start and end cords.
 * NOTE::structure for align_cords() function, different from gaps defined in 'gaps.h';
 */
class GapRecords
{
    
public:
    typedef Row<Align<String<Dna5>, ArrayGaps> >::Type TRow;
    typedef std::pair<uint64_t, uint64_t> TCPair;
    typedef std::pair<TRow, TRow> TRPair;
    typedef String<std::pair<uint64_t, uint64_t> > TCPairs;
    typedef String<std::pair<TRow, TRow> > TRPairs;

    TCPairs c_pairs;  //pair of start and end coordinate of gaps 
    TRPairs r_pairs;  //pair of alingments of blocks contains the start and end of gaps 
    String<int> bam_segs_id;
    String<int> clip_flags;
    int unset_bit;
    uint64_t dx, dy;  //shift from start and end of gaps to the start and end of of the alignment

    TCPair & get_c_pair(int i);  
    TRPair & get_r1_pair(int i); //get rows of first(start) cord of i_th gaps
    TRPair & get_r2_pair(int i); //get rows of second(end) cord of i_th gaps
    int getBamSegIdHead (int i);
    int getBamSegIdTail (int i);
    int get_clip_flag(int i);
    uint64_t getJointHeadCord(int i); 
    uint64_t getJointTailCord(int i);
    int set_clip_flag(int clip_flag, int pos);
    int clear_();
};

struct GapParm
{
    int thd_clip_score;
    int thd_reject_score;
    int thd_accept_score;
    int thd_min_interval;
    float thd_accept_density;
    Score<int, Simple> thd_clip_scheme;
    GapParm ();
};
extern GapParm _gap_parm;
extern unsigned _default_block_size_;
int align_cords (StringSet<String<Dna5> >& genomes,
                 String<Dna5> & read, 
                 String<Dna5> & comrevRead,
                 String<uint64_t> & cords,
                 String<BamAlignmentRecordLink> & bam_records,
                 int block_size = _default_block_size_,
                 int band = _default_block_size_ / 2
                );
void printCigar(String<CigarElement< > > &, std::string);


#endif