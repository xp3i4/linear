#include "base.h"
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
                        int f_soft,/*hard and soft clip flag*/
                        uint16_t flag
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
    bam_record.flag |=  flag;
    if (strand)
    {
        bam_record.flag |= bam_flag_rvcmp ;
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
/*------------------- BamLinkStringOperator and SAZ tag --------------------*/
/*
 * shrink the @cigar to @cigar_simple keeping start and end identical to @cigar.
 * Format of @cigar_simple is x1Sx2Mx3I(D)x4S : 
   x1=S of cigar x2 = sum(M) of cigar;
   x2=sum(I)-sum(D) of cigar, if x2 > 0 then I if x2=0 none else D
   x3=S of cigar
 * @f_remove : remove 0 operation of cigar elment in @cigar_simple;
 * Hard clip ('H') is not supported yet
 * 'M' is not suppoerted yet either
 * @nm_i is to create nm:i tag. 
    It's to save the times of iterating cigar that the tag is created in this function 
    simultaneously rather than creating antoher indepent function;
    The 'X' and '=' is reuqired to create this tag, hence the 'M' is not supported.
 */
int createSAZTagCigar(String<CigarElement<> > & cigar,
                      String<CigarElement<> > & cigar_simple,
                      int & nm_i,
                      int f_remove) 
{
    clear(cigar_simple);
    appendValue(cigar_simple, CigarElement<>('S', 0));
    appendValue(cigar_simple, CigarElement<>('M', 0));
    appendValue(cigar_simple, CigarElement<>('I', 0));
    appendValue(cigar_simple, CigarElement<>('S', 0));
    unsigned cm = 0, ci = 0;
    nm_i = 0;
    for (unsigned i = 0; i < length(cigar); i++)
    {
        if (i == 0 && cigar[i].operation == 'S')
        {
            cigar_simple[0].count = cigar[i].count;
        }
        else if (cigar[i].operation == '=')
        {
            cm += cigar[i].count;
        }
        else if (cigar[i].operation == 'X')
        {
            cm += cigar[i].count;
            nm_i += cigar[i].count;
        }
        else if (cigar[i].operation == 'I')
        {
            ci -= cigar[i].count;
            nm_i += cigar[i].count;
        }
        else if (cigar[i].operation == 'D')
        {
            ci += cigar[i].count;
            nm_i += cigar[i].count;
        }
        else if (i == length(cigar[i]) - 1 && cigar[i].operation == 'S')
        {
            cigar_simple[3].count = cigar[i].count;
        }
    }
    cigar_simple[1] = CigarElement<>('M', cm);
    cigar_simple[2] = ci < 0 ? CigarElement<>('I', ci) : CigarElement<>('D', ci);
    unsigned it = 0;
    //remove 0operation in cigar
    if (f_remove)
    {
        for (unsigned i = 0; i < length(cigar_simple); i++)
        {
            if (cigar_simple[i].count == 0)
            {
                it++;
            }
            else
            {
                cigar_simple[i - it] = cigar_simple[i];
            }
        }
        resize (cigar_simple, length(cigar_simple) - it);
    }

    return 0;
}
/*
 * merge saz_cigar with cigar_result
   saz_cigar and cigar_result are in the format of cigar in the function 
   createSAZCigar
 */
int mergeSAZTagCigar(String<CigarElement<> > & cigar_result,
                  String<CigarElement<> > & saz_cigar)
{
    for (unsigned i = 0; i < length(saz_cigar); i++)
    {
        if (i >= length(cigar_result))
        {
            appendValue(cigar_result, saz_cigar[i]);
        }
        else
        {
            cigar_result[i].count += saz_cigar[i].count;
        }
    }
    return 0;
}
/*
 *  Return length of head table
 */
int BamLinkStringOperator::getHeadNum(String<BamAlignmentRecordLink> & bam_records)
{
    return empty(bam_records) ? 0 : length(bam_records[0].heads_table);
}
/*
 *  Get pointer of ith head from the head table
 */
int BamLinkStringOperator::getHead(String<BamAlignmentRecordLink> & bam_records, int i)
{
    return empty(bam_records) ? -1 : bam_records[0].heads_table[i];
}
/*
 *  Update heads table of @bam_records.
    The heads table is store in @bam_records[0].
 */
int BamLinkStringOperator::updateHeadsTable(String<BamAlignmentRecordLink> & bam_records)
{
    dout << "<<<<uht" << "\n";
    if (empty(bam_records))
    {
        return 0;
    }
    String<bool> visit_tmp;
    resize(visit_tmp, length(bam_records), false);
    BamAlignmentRecordLink & info_record = bam_records[0];
    clear(info_record.heads_table);
    for (unsigned i = 0; i < length(bam_records); i++)
    {
        if (!visit_tmp[i])
        {
            unsigned it = i;
            while (1)
            {
                visit_tmp[it] = true;
                if (bam_records[it].isEnd())
                {
                    break;
                }
                else
                {
                    it = bam_records[it].next();
                }
            }
            appendValue(info_record.heads_table, i);
        }
    }
    return 0;
}
/*
 * Write the cigars of the link of bam records (one line in .sam) 
   whose first record is specified by it
 */
int BamLinkStringOperator::writeBamRecordLinkCigar(
        std::ofstream target,
        String<BamAlignmentRecordLink> & records,
        int it)
{
    while (1)
    {
        for (unsigned i = 0; i < length(records[it].cigar); ++i)
        {
            appendNumber(target, records[it].cigar[i].count);
            writeValue(target, records[it].cigar[i].operation);
        }
        if (records[it].isEnd())
        {
            break;
        }
        else 
        {
            it = records[it].next();
        }
    }
    return 0;
}
/*
 * Write cigar used in saz for one line(record) in .sam
 * set @f_force = ture to update the saz_cigar anyway.
 * This function creates nm:i tag at the same time to save times of iterating cigars.
 */
int BamLinkStringOperator::createSAZTagCigarOneChimeric(
        String<BamAlignmentRecordLink> & bam_records,
        String<CigarElement<> > & cigar,
        int it,
        bool f_force)
{
    clear(cigar);
    int neg_infi = -999; //random value (<0 required) used to initie nm_i_sum;
    int nm_i_sum = neg_infi, nm_i = -1;
    BamAlignmentRecordLink & record = bam_records[it];
    while (1)
    {
        if (f_force || empty(bam_records[it].saz_cigar))
        {
            createSAZTagCigar(bam_records[it].cigar, bam_records[it].saz_cigar, nm_i, 0);
            nm_i_sum += nm_i;
        }
        mergeSAZTagCigar(cigar, bam_records[it].saz_cigar);
        if (bam_records[it].isEnd())
        {
            break;
        }
        else 
        {
            it = bam_records[it].next();
        }
    }
    if (nm_i_sum != neg_infi)
    {
        record.nm_i = nm_i_sum - neg_infi;
    }
    return it;
}
/*
 * Create sa:z tag for one chimeric alignment of sa:z
 */ 
int BamLinkStringOperator::createSAZTagOneChimeric(
        String<BamAlignmentRecordLink> & bam_records,
        CharString & saz_tag,
        int it)
{
    BamAlignmentRecordLink & record = bam_records[it];
    CharString cigar_string;
    String<CigarElement<> > cigar;
    createSAZTagCigarOneChimeric(bam_records, cigar, it, false);
    for (unsigned i = 0; i < length(cigar); i++)
    {
        //convert cigar to CharString
        append(cigar_string, std::to_string(cigar[i].count));
        appendValue(cigar_string, cigar[i].operation);
    }
    append(saz_tag, record.genome_id);
    appendValue(saz_tag, ',');
    append(saz_tag, std::to_string(record.beginPos + 1)); //+1 to 1-based
    appendValue(saz_tag, ',');
    char strand = int(record.flag) & 16 ? '-' : '+';
    append(saz_tag, strand);   
    appendValue(saz_tag, ',' );
    append(saz_tag, cigar_string);
    appendValue(saz_tag, ',' );
    append(saz_tag, std::to_string(record.mapQ));
    appendValue(saz_tag, ',');
    append(saz_tag, std::to_string(record.nm_i));
    appendValue(saz_tag, ';' );

    return 0;
} 
/*
 *  create  sa:z tag for one line, whcih here refers to line in .sam format;
 *  The function simple makes all supplement alignments as chimeric alignments 
    and insert them to the sa:z tag.
 *  @it is supposed to be pointing to head.
 */
int BamLinkStringOperator::createSAZTagOneLine(
        String<BamAlignmentRecordLink> & bam_records,
        int it)
{
    if (empty(bam_records))
    {
        return 0;
    }
    CharString saz_tag("SA:Z:");
    updateHeadsTable(bam_records);
    for (int i = 0; i < getHeadNum(bam_records); i++)
    {
        int j = getHead(bam_records, i);
        if (it != j)
        {
            createSAZTagOneChimeric(bam_records, saz_tag, j);
        }
    }
    append(bam_records[it].tags, saz_tag);
    return 0;
}
