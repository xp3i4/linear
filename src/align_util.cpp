#include "cords.h"
#include "align_util.h"

uint16_t bam_flag_rvcmp = 16;
uint16_t bam_flag_rvcmp_nxt = 32;
uint16_t bam_flag_suppl = 2048;

BamAlignmentRecordLink::BamAlignmentRecordLink()
{
    next_id = -1;
    BamAlignmentRecord();
}
void BamAlignmentRecordLink::addNext(int id)
{
    next_id = id;
}
void addNextBamLink(String<BamAlignmentRecordLink> & bam_records,
                    int id, int next_id)
{
    bam_records[id].addNext(next_id);

    String<CigarElement<> > & cigar = bam_records[next_id].cigar;
    if (!empty(cigar))
    {
        if (cigar[0].operation == 'S' || cigar[0].operation == 'H')
        {
            //eraseFront(cigar);
            erase(cigar,0);
        }    
    }
    //todo::merge 1=|2= to 3=
}
int BamAlignmentRecordLink::isEnd() const 
{
    return next_id < 0;
}
int BamAlignmentRecordLink::next() const 
{
    return next_id;
}

void align2cigar_(String<CigarElement< > > &cigar,
        Row<Align<String<Dna5>,ArrayGaps> >::Type &gaps1,
        Row<Align<String<Dna5>,ArrayGaps> >::Type &gaps2,
        unsigned splicedGapThresh
        )
{
    typedef Row<Align<String<Dna5>,ArrayGaps> >::Type TRow;
    typename Iterator<TRow>::Type it1 = begin(gaps1);
    typename Iterator<TRow>::Type it2 = begin(gaps2);
    char op = '?', lastOp = ' ';
    char last_op;
    unsigned numOps = 0;
    unsigned last_count;
    if (!empty(cigar))
    {
        last_op = back(cigar).operation;
        last_count = back(cigar).count;
    }
    int flag = 0;
    for (; !atEnd(it1) && !atEnd(it2); goNext(it1), goNext(it2))
    {
        if (isGap(it1))
        {
            if (isGap(it2))
                op = 'P';
            else if (isClipped(it2))
                op = '?';
            else
                op = 'I';
        }
        else if (isClipped(it1))
        {
            op = '?';
        }
        else
        {
            if (isGap(it2))
                op = 'D';
            else if (isClipped(it2))
                op = 'S';
            else
                op = (*it1 == *it2)? '=': 'X';
//                op = 'M';
        }
        if (lastOp != op)
        {
            //if (lastOp == 'D' && numOps >= (unsigned)splicedGapThresh)
            //    lastOp = 'N';
            if (numOps > 0)
            {
                if (!flag)
                {
                    if (last_op == lastOp)
                    {
                        back(cigar).count += numOps;
                    } 
                    else
                    {
                        appendValue(cigar, CigarElement<>(lastOp, numOps));
                    }
                    flag = 1;
                }
                else
                {
                    appendValue(cigar, CigarElement<>(lastOp, numOps));
                } 
            }
            numOps = 0;
            lastOp = op;
        }
        ++numOps;
    }
    SEQAN_ASSERT_EQ(atEnd(it1), atEnd(it2));
    //if (lastOp == 'D' && numOps >= splicedGapThresh)
    //    lastOp = 'N';
    if (numOps > 0)
        appendValue(cigar, CigarElement<>(op, numOps));
}
void align2cigar(String<CigarElement< > > &cigar,
                 Row<Align<String<Dna5>,ArrayGaps> >::Type &gaps1,
                 Row<Align<String<Dna5>,ArrayGaps> >::Type &gaps2
                )
{
    align2cigar_(cigar, gaps1, gaps2, 1000);
}

int clip_cigar (String<CigarElement<> > & cigar)
{
    int x = 0, y = 0;
    int score = 0; 
    for (int i = 0; i < length(cigar); i++)
    {

        switch(cigar[i].operation)
        {
            case 'D':
                x += cigar[i].count;
                break;
            case 'I':
                y += cigar[i].count;
                break;
            case '=':
                x += cigar[i].count;
                y += cigar[i].count;
                break;
            case 'X':
                x += cigar[i].count;
                y += cigar[i].count;  
                break;
            case 'M':
                return 1;        //'M' is not allowed in the function
            case 'S': 
                break;
            //case 'N':
            //    break;
            case 'P':
                break;
            default:
                return 2;
        }
        std::cout << "[]::clip_cigar "  << " " << x << " " << y << " " << cigar[i].count << cigar[i].operation << "\n";
    }
    return 0;
}
/*
 * Insert @cigr2 to @cigar1 at @pos 
 */
int insertCigar(String<CigarElement< > > &cigar1, 
                int pos,
                String<CigarElement< > > &cigar2)
{
    int p = pos;
    if (empty(cigar1))
    {
        cigar1 = cigar2;
        return 0;
    }
    if (pos < 0)
    {
        return 1;
    }
    else if (pos > length(cigar1))  
    {
        p = length(cigar1);
    }
    if (p == 0) //insert at front 
    {
        if (cigar1[0].operation == back(cigar2).operation)
        {
            std::cout << "insertCigar p = 0 " << length(cigar1) << " " << length(cigar2) << "\n";
            cigar1[0].count += back(cigar2).count;
            eraseBack(cigar2);
        }
        insert(cigar1, p, cigar2);
        return 0;
    }
    if (p == length(cigar1)) //insert at end
    {
        if (back(cigar1).operation == cigar2[0].operation)
        {
            cigar2[0].count += back(cigar1).count;
            eraseBack(cigar1);
        }
        append(cigar1, cigar2);
        return 0;
    }
    if (cigar1[p - 1].operation == cigar2[0].operation)  //insert in middle
    {
        cigar1[p - 1].count += cigar2[0].count;
        erase(cigar2, 0);
    }
    if (cigar1[p].operation == back(cigar2).operation) 
    {
        cigar1[p].count += back(cigar2).count;
        eraseBack(cigar2);
    }
    insert(cigar1, p, cigar2);
    return 0;
}

/*
 * insert cigar to the original cigar 
 */
int insertBamRecordCigar (BamAlignmentRecord & bam_record,
                    Row<Align<String<Dna5>, ArrayGaps> >::Type & row1,
                    Row<Align<String<Dna5>, ArrayGaps> >::Type & row2, 
                    int pos
                   )
{
    if (pos < 0)
    {
        align2cigar(bam_record.cigar, row1, row2);
    }
    else
    {
        if (pos > length(bam_record.cigar ) - 1)
        {
            return 1;
        }
        String<CigarElement< > > tmp;
        align2cigar(tmp, row1, row2);
        insertCigar(bam_record.cigar, pos, tmp);
    }
    return 0;
}
/**
 * @g_beignPos, @r_beginPos
 * 1-based leftmost exact start coordinates 
 */
int insertNewBamRecord (String<BamAlignmentRecordLink> & bam_records,
                        int g_id, 
                        int g_beginPos,
                        int r_beginPos,
                        int strand,
                        int insert_pos,
                        int f_soft/*hard and soft clip flag*/
                        )
{
    BamAlignmentRecordLink bam_record;
    if (g_id >= 0)
    {
        bam_record.rID = g_id;
    }
    if (g_beginPos >= 0)
    {
        bam_record.beginPos = g_beginPos; 
    }
    if (strand)
    {
        bam_record.flag |= bam_flag_rvcmp | bam_flag_suppl;
    }
    if (r_beginPos != 0)
    {
        char op = 'S';
        if (!f_soft)
        {
            op = 'H';
        }
        appendValue(bam_record.cigar, CigarElement<>(op, r_beginPos));
    }
    if (insert_pos < 0)
    {
        appendValue(bam_records, bam_record);
    }
    else 
    {
        insert(bam_records, insert_pos, bam_record);
    }
    return 0;
}
int insertNewBamRecord (String<BamAlignmentRecordLink> & bam_records,
                        Row<Align<String<Dna5>, ArrayGaps> >::Type & row1,
                        Row<Align<String<Dna5>, ArrayGaps> >::Type & row2, 
                        int g_id,
                        int g_beginPos,
                        int r_beginPos,
                        int strand,
                        int insert_pos,
                        int f_soft 
                        )
{
    BamAlignmentRecordLink bam_record;
    insertNewBamRecord(bam_records, g_id, g_beginPos, r_beginPos, strand, insert_pos, f_soft);
    int i = insert_pos;
    if (insert_pos < 0)
    {
        i = length(bam_records) - 1;
    }
    insertBamRecordCigar(bam_records[i], row1, row2);
    return 0;
}

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
                      int pos,
                      int f_soft
                      )
{
    String<CigarElement<> > & cigar = bam_record.cigar;
    char op = 'S';
    if (!f_soft)
    {
        op = 'H';
    }

    if (empty (bam_record.cigar))
    {
        insertBamRecordCigar(bam_record, row1, row2, 0);
        if (r_beginPos != 0)
        {
            insertValue (cigar, 0, CigarElement<>(op, r_beginPos));
        }
    }
    else
    {
        if (pos == 0 && (cigar[0].operation == 'S' || cigar[0].operation == 'H'))
        {
            insertBamRecordCigar(bam_record, row1, row2, 1);
            cigar[0].count = r_beginPos;
            cigar[0].operation = op;
        }
        else
        {
            insertBamRecordCigar(bam_record, row1, row2, pos);
            if (cigar[0].operation == 'S' || cigar[0].operation == 'H')
            {
                cigar[0].count = r_beginPos;
                cigar[0].operation = op;
            }
            else
            {
                insertValue (cigar, 0, CigarElement<>(op, r_beginPos));
            }
        }
    }
    if (g_id >= 0)
    {
        bam_record.rID = g_id;
    }
    else
    {
        return 1;
    }
    if (g_beginPos >= 0)
    {
        bam_record.beginPos = g_beginPos; 
    }
    else
    {
        return 1;
    }
    (void)r_beginPos;
    return 0;
}