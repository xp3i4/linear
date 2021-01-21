#ifndef LINEAR_HEADER_ALIGN_UTIL_H
#define LINEAR_HEADER_ALIGN_UTIL_H

#include <seqan/bam_io.h>

using namespace seqan;

extern uint16_t bam_flag_rvcmp;
extern uint16_t bam_flag_rvcmp_nxt;
extern uint16_t bam_flag_suppl;
typedef Row<Align<String<Dna5>,ArrayGaps> >::Type TRow5A;  
struct SAZTag
{
    String<unsigned> bam_records_i;
    String<bool> bam_records_is_chimeric;
    String<bool> bam_records_is_block_end;
};
class BamAlignmentRecordLink : public BamAlignmentRecord 
{ 
public:
    //the heads refers to the first record of each line 
    //line refers to each line in.sam
    //line is compsed of several records in the String<bamAlig...>
    int next_id; //next records id
    String<CigarElement<> > saz_cigar;
    String<unsigned> heads_table; //table pointing to first record of each line in .sam 
    CharString genome_id;
    int nm_i; //nm:i tag in .sam, the heads holds the whole(sum of) value of the line
    BamAlignmentRecordLink();
    void addNext(int id);
    int isEnd() const;
    int next() const;
};
/*
 * The set of functions of manipulating String<BamAlignmentRecordLink> 
 * The functions are supposed to be wrappered with the String<BamAlignmentRecordLink> to
   constitute a new class. However the String<BamAlignmentRecordLink> has been widely used
   by many functions. To avoid modifications of interface this struct is thus declared.
 */
struct BamLinkStringOperator
{
    int updateHeadsTable(String<BamAlignmentRecordLink> & bam_records);
    int getHeadNum(String<BamAlignmentRecordLink> & bam_records);
    int getHead(String<BamAlignmentRecordLink> & bam_records, int i);
    int getLineStrand(String<BamAlignmentRecordLink> & bam_records, int i);
    int writeBamRecordLinkCigar(
            std::ofstream target,
            String<BamAlignmentRecordLink> & bam_records,
            int it);
    int createSAZTagCigarOneChimeric(
            String<BamAlignmentRecordLink> & bam_records,
            String<CigarElement<> > & cigar,
            int it,
            bool f_force);
    int createSAZTagOneChimeric(
            String<BamAlignmentRecordLink> & bam_records,
            CharString & saz_tag, 
            int it);
    int createSAZTagOneLine(
            String<BamAlignmentRecordLink> & bam_records,
            int it);
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
int insertBamRecordCigar (BamAlignmentRecordLink & bam_record,
                    Row<Align<String<Dna5>, ArrayGaps> >::Type & row1,
                    Row<Align<String<Dna5>, ArrayGaps> >::Type & row2, 
                    int pos = -1
                   );

/*
 * Cigar of row1 and row2 are inserted to the cigar from the bam_record
 * beginPos are always updated by g_beginPos 
 * soft/Hard clip cigar are updated only if insert at the front(pos == 0)
 */
int  insertBamRecord (BamAlignmentRecordLink & bam_record,
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
                        int f_soft = 1, /*hard and soft clip flag*/
                        uint16_t flag = 0
                        );

int insertNewBamRecord (String<BamAlignmentRecordLink> & bam_records,
                        Row<Align<String<Dna5>, ArrayGaps> >::Type & row1,
                        Row<Align<String<Dna5>, ArrayGaps> >::Type & row2, 
                        int g_id,
                        int g_beginPos,
                        int r_beginPos,
                        int strand,
                        int insert_pos = -1,
                        int f_soft = 1,
                        uint16_t flag = 0
                        );
/*----------  align pos cache  ----------*
 * Tempory struct to buffer coordinates of rows
   to reduce times of translation of coordinates of rows.
 */
struct AlignCache_
{
    int64_t _uc_view;
    int64_t _g_src_x;
    int64_t _g_src_y;
    AlignCache_(int64_t uc_view, int64_t g_src_x, int64_t g_src_y, int64_t cord_0);
};
class AlignCache
{
    String<AlignCache_> cache;
public:
    void appendValue(int64_t uc_view, int64_t g_src_x, int64_t g_src_y, int64_t cord_0);
    int64_t getUCView(int i);
    int64_t getGSrcX(int i);
    int64_t getGSrcY(int i);
    unsigned length();
    bool empty();
};
/*----------  Viewer of pos of rows ----------*/
class RowPosViewer
{
  public:
    TRow5A *rp; 
    int64_t uc_view;
    int64_t src;
    int64_t array_i; //points to _array[array_i] 
    int f_array; //0 array_i point to gap bucket otherwise1 to sequence bucket
    int64_t view_bucket_sum;
    int64_t init(TRow5A & row);
    //uint64_t nextArrayBucket();
    int64_t nextView();
    int64_t findSrc(int64_t src);
    int64_t findView(int64_t view);
    int64_t getView();
    int64_t getSrc();
    int64_t getUCView();
    RowPosViewer(TRow5A & row);
};
int mergeAlign2_(Row<Align<String<Dna5>,ArrayGaps> >::Type & row11,
                 Row<Align<String<Dna5>,ArrayGaps> >::Type & row12,
                 Row<Align<String<Dna5>,ArrayGaps> >::Type & row21,
                 Row<Align<String<Dna5>,ArrayGaps> >::Type & row22,
                 AlignCache & align1,
                 AlignCache & align2,
                 String<Dna5> & ref,
                 String<Dna5> & read,
                 String<Dna5> & comrev_read,
                 uint64_t & cord1,
                 uint64_t & cord2);
std::pair<int, int> cigar2SeqLen(CigarElement<> & cigar);
std::pair<int, int> cigars2SeqsLen(String<CigarElement<> > & cigars,
                                   unsigned c_str, unsigned c_end);
#endif
