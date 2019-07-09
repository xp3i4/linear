#ifndef LINEAR_HEADER_ALIGN_UTIL_H
#define LINEAR_HEADER_ALIGN_UTIL_H

#include <seqan/bam_io.h>

using namespace seqan;

class BamAlignmentRecordLink : public BamAlignmentRecord 
{ 
public:
    int next_id; //next records id

    BamAlignmentRecordLink();
    void addNext(int id);
    int isEnd() const;
    int next() const;
};

void align2cigar(String<CigarElement< > > &cigar,
                 Row<Align<String<Dna5>,ArrayGaps> >::Type &gaps1,
                 Row<Align<String<Dna5>,ArrayGaps> >::Type &gaps2
                );
int insertCigar(String<CigarElement< > > &cigar1, 
                int pos,
                String<CigarElement< > > &cigar2
         );

/*
 * insert cigar to the original cigar 
 */
int insertBamRecordCigar (BamAlignmentRecord & bam_record,
                    Row<Align<String<Dna5>, ArrayGaps> >::Type & row1,
                    Row<Align<String<Dna5>, ArrayGaps> >::Type & row2, 
                    int pos = -1
                   );

/*
 * Cigar of row1 and row2 are inserted to the cigar from the bam_record
 * beginPos are always updated by g_beginPos 
 * soft/Hard clip cigar are updated only if insert at the front(pos == 0)
 */
int  insertBamRecord (BamAlignmentRecord & bam_record,
                      Row<Align<String<Dna5>, ArrayGaps> >::Type & row1,
                      Row<Align<String<Dna5>, ArrayGaps> >::Type & row2, 
                      int g_id,
                      int g_beginPos,
                      int r_beginPos,
                      int pos = -1,
                      int f_soft = 1
                      );

/**
 * @g_beignPos, @r_beginPos
 * 1-based leftmost exact start coordinates 
 */
int insertNewBamRecord (String<BamAlignmentRecordLink> & bam_records,
                        int g_id, 
                        int g_beginPos,
                        int r_beginPos,
                        int strand,
                        int insert_pos = -1,
                        int f_soft = 1 /*hard and soft clip flag*/
                        );

int insertNewBamRecord (String<BamAlignmentRecordLink> & bam_records,
                        Row<Align<String<Dna5>, ArrayGaps> >::Type & row1,
                        Row<Align<String<Dna5>, ArrayGaps> >::Type & row2, 
                        int g_id,
                        int g_beginPos,
                        int r_beginPos,
                        int strand,
                        int insert_pos = -1,
                        int f_soft = 1 
                        );


#endif
