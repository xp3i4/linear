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
        //std::cout << "[]::clip_cigar "  << " " << x << " " << y << " " << cigar[i].count << cigar[i].operation << "\n";
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
            //std::cout << "insertCigar p = 0 " << length(cigar1) << " " << length(cigar2) << "\n";
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
int insertBamRecordCigar (BamAlignmentRecordLink & bam_record,
                          Row<Align<String<Dna5>, ArrayGaps> >::Type & row1,
                          Row<Align<String<Dna5>, ArrayGaps> >::Type & row2, 
                          int pos)
{
    if (pos < 0)
    {
        align2cigar(bam_record.cigar, row1, row2);
    }
    else
    {
        if (pos > length(bam_record.cigar ))
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
                        uint16_t flag)
{
    //dout << "ib3" << g_beginPos << r_beginPos << "\n";
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
                        int f_soft,
                        uint16_t flag)
{
    BamAlignmentRecordLink bam_record;
    insertNewBamRecord(bam_records, g_id, g_beginPos, r_beginPos, strand, insert_pos, f_soft, flag);
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
int insertBamRecord (BamAlignmentRecordLink & bam_record,
                     Row<Align<String<Dna5>, ArrayGaps> >::Type & row1,
                     Row<Align<String<Dna5>, ArrayGaps> >::Type & row2, 
                     int g_id,
                     int g_beginPos,
                     int r_beginPos,
                     int pos,
                     int f_soft)
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
            dout << "is1" << g_beginPos << length(bam_record.cigar) << "\n";
            insertBamRecordCigar(bam_record, row1, row2, 1);
            cigar[0].count = r_beginPos;
            cigar[0].operation = op;
            dout << "is1" << g_beginPos << length(bam_record.cigar) << "\n";
            std::cout << "is1" << row1 << "\n";
            std::cout << "is1" << row2 << "\n";
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
 * Get the clippeed start position of y(read) of the line whose head is the ith one
 */
int BamLinkStringOperator::getLineStrand(String<BamAlignmentRecordLink> & bam_records, int i)
{
    return bam_records[getHead(bam_records, i)].flag & 16;
}
/*
 *  Update heads table of @bam_records.
    The heads table is stored in @bam_records[0].
 */
int BamLinkStringOperator::updateHeadsTable(String<BamAlignmentRecordLink> & bam_records)
{
    //dout << "<<<<uht" << "\n";
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
    else
    {
        record.nm_i = 0;
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
 *  @it is supposed to be pointer to head.
 */
int BamLinkStringOperator::createSAZTagOneLine(
        String<BamAlignmentRecordLink> & bam_records,
        int it)
{
    if (empty(bam_records))
    {
        return 0;
    }
    updateHeadsTable(bam_records);
    CharString saz_tag("SA:Z:");
    CharString empty_tag = saz_tag;
    for (int i = 0; i < getHeadNum(bam_records); i++)
    {
        int j = getHead(bam_records, i);
        if (it != j)
        {
            createSAZTagOneChimeric(bam_records, saz_tag, j);
        }
    }
    if (saz_tag != empty_tag)
    {
        append(bam_records[it].tags, saz_tag);
    }
    return 0;
}
/*----------  row viewer customize  ----------*/
int64_t RowPosViewer::init(TRow5A & row) //cordxy:coordinates of x or y
{
    rp = &row ;
    src = 0;
    array_i = 0;
    f_array = 0;
    view_bucket_sum = empty(rp->_array) ? 0 : rp->_array[0];
    uc_view = (view_bucket_sum == 0) ? 0 : view_bucket_sum - 1;
    return 0;
}
RowPosViewer::RowPosViewer(TRow5A & row)
{
    init (row);
}
int64_t RowPosViewer::nextView()
{
    ++uc_view;
    if (uc_view > view_bucket_sum) 
    {
        ++array_i;
        view_bucket_sum += rp->_array[array_i];
        f_array = 1 - f_array;
    }
    if (f_array)
    {
        ++src;
    }
    //dout << "nextView" << uc_view << array_i << length(rp->_array) << src << endPosition(*rp) << "\n";
    return 0;
}
int64_t RowPosViewer::findSrc(int64_t find_src)
{
    init(*rp);
    //rp->_clippingBeginPos = 25;
    for (; array_i < length(rp->_array); )
    {
        if (f_array && src >= find_src)
        {
            uc_view -= src - find_src;
            src = find_src;
            break;
        }
        if (++array_i < length(rp->_array))
        {
            view_bucket_sum += rp->_array[array_i];
            uc_view = view_bucket_sum - 1;
            f_array = 1 - f_array;
            src = f_array ? src + rp->_array[array_i] - 1 : src + 1;
        }
    }
    return 1;
}
int64_t RowPosViewer::findView(int64_t find_view)
{
    init(*rp);
    for (;array_i < length(rp->_array);)
    {
        
        if (uc_view - rp->_clippingBeginPos >= find_view)
        {
            if (f_array)
            {
                src -= uc_view - find_view - rp->_clippingBeginPos;
            }
            uc_view = find_view + rp->_clippingBeginPos;
            break;
        }
        if (++array_i < length(rp->_array))
        {
            view_bucket_sum += rp->_array[array_i];
            uc_view = view_bucket_sum - 1;
            f_array = 1- f_array;
            src = f_array ? src + rp->_array[array_i] - 1: src + 1;
        }
    }
    return 0;
}
int64_t RowPosViewer::getView()
{
    return uc_view - rp->_clippingBeginPos;
}
int64_t RowPosViewer::getSrc()
{
    return src;
}
int64_t RowPosViewer::getUCView()
{
    return uc_view;
}
/*----------  align pos cache  ----------*/
AlignCache_::AlignCache_(int64_t uc_view, int64_t g_src_x, 
    int64_t g_src_y, int64_t cord_0)
{
    _uc_view = uc_view;
    _g_src_x = g_src_x + get_cord_x(cord_0);
    _g_src_y = g_src_y + get_cord_y(cord_0);
}
void AlignCache::appendValue(int64_t uc_view, int64_t g_src_x, int64_t g_src_y, 
    int64_t cord_0)
{
    seqan::appendValue(cache, AlignCache_(uc_view, g_src_x, g_src_y, cord_0));
}
int64_t AlignCache::getUCView(int i)
{
    return cache[i]._uc_view;
}
int64_t AlignCache::getGSrcX(int i)
{
    return cache[i]._g_src_x;
}
int64_t AlignCache::getGSrcY(int i)
{
    return cache[i]._g_src_y;
}
unsigned AlignCache::length()
{
    return seqan::length(cache);
}
bool AlignCache::empty()
{
    return seqan::empty(cache);
}

/*----------  mergeAlign2_  ----------*/
/* WARN:!
 * Remove row (._source and ._array) from left to @view_pos. 
 * This function depends on the definition of seqan::ArrayGaps of seqan 2
   Hence keep the function consitent with ArrayGaps in seqan
 *
   seqan::Gaps._array definition: 
   even bucket is count of chars, odd bucket is count of gaps.
   The first and the last bucket are for gaps.
 * 
   The function keeps clipped view pos unchanged,
   namley they point to the original bucket.
 */
int64_t seqanRemoveHead_(Row<Align<String<Dna5>,ArrayGaps> >::Type & row)
{
    uint64_t view_pos = clippedBeginPosition(row);
    unsigned pos = 0; 
    unsigned src_count = row._sourceBeginPos;
    /*
    RowPosViewer rp(row);
    rp.findView(view_pos);
    row._array[rp.array_i] = rp.view_bucket_sum - view_pos;
    if (rp.array_i % 2 == 0)
    {
        erase(row._array, 0, rp.array_i);
    }
    else
    {
        if (row._array[rp.array_i] == 0)
        {
            erase(row._array, 0, rp.array_i + 1);
        }
        else
        {
            erase(row._array, 0, rp.array_i - 1);
            row._array[0] = 0;
        }
    }
    if (length(row._array) < 3)
    {
        resize (row._array, 3);
        row._array[0] = 0;
        row._array[1] = 0;
        row._array[2] = 0;
    }
    setValue(row._source, infix(value(row._source), row._sourceBeginPos, length(value(row._source))));
    src_count = row._sourceBeginPos;
    row._sourceEndPos -= row._sourceBeginPos;
    row._sourceBeginPos = 0;  //0
    row._clippingEndPos -= row._clippingBeginPos;
    row._clippingBeginPos -= 0; //0
    dout << "srh1" << row._sourceBeginPos << " " << row._sourceEndPos << " " << row._clippingBeginPos << " " << row._clippingEndPos << " " << src_count << " " << view_pos << length(row._array) << "\n"; 
    return src_count;
    */
    
    for (int i = 0 ; i < length(row._array); i++) 
    {
        pos += row._array[i];
        //std::cerr << pos << " "<< view_pos << " " << i << " " << length(row._array) << "\n";
        int f_g = (i + 1) % 2; //if i is the gap bucket
        if (pos > view_pos)
        {
            row._array[i] = pos - view_pos;
            if (!f_g)
            {
                //src_count -= row._array[i];
                erase(row._array, 0, i - 1);
                row._array[0] = 0;
            }
            else
            {
                erase(row._array, 0, i);
            }
            setValue(row._source, infix(value(row._source), row._sourceBeginPos, length(value(row._source))));
            row._sourceEndPos -= row._sourceBeginPos;
            row._sourceBeginPos = 0;  //0
            row._clippingEndPos -= row._clippingBeginPos;
            row._clippingBeginPos = 0; //0
      //      dout << "srh1" << row._sourceBeginPos << " " << row._sourceEndPos << " " << row._clippingBeginPos << " " << row._clippingEndPos << " " << src_count << " " << view_pos << "\n";
            break;
        }
    }
    //std::cout << row << "\n";
     //       dout << "srh11" << row._sourceBeginPos << " " << row._sourceEndPos << " " << row._clippingBeginPos << " " << row._clippingEndPos << " " << src_count << " " << view_pos << length(row._array) << "\n";
    return src_count;
    
    /*
    int d = 0;
    for (int i = 0; i < length(row._array); i++)
    {
        if (i%2 ==1)
        d += row._array[i];
        std::cerr << "row " << row._array[i] << "\n" ;
    }
    std::cerr << "d" << d << " " << row._sourceEndPos << " " << row._sourceBeginPos << " " << length(value(row._source)) << "\n";
    std::cerr << value(row._source) << "\n";
    return src_count;
    */
}
/*
 * Same as the seqanRemoveHead_
 * Extend leftward from clippedBeginPosition(@row) with character in @seq by @extend_len
 * @extend_len > beginSource,
 * The function keeps clipped view pos unchanged.
 */
int64_t seqanExtendHead_(Row<Align<String<Dna5>,ArrayGaps> >::Type & row,
                     String<Dna5> & seq, 
                     uint64_t extend_len,
                     uint64_t src_cord_0)
{
    if (empty(row._array) || row._sourceBeginPos > extend_len)
    {
        return 0;
    }
    //for (int i = 0; i < length(row._array); i++)
    //{
    //    dout << "seq1" << row._array[i] << "\n";
    //}
    uint64_t seq_pos = src_cord_0 + row._sourceBeginPos - extend_len;
    //std::cout << "seq1" << value(row._source) << "\n";
    setValue(row._source, infix(seq, seq_pos, seq_pos + length(value(row._source)) + extend_len - row._sourceBeginPos));
    //std::cout << "seq2" << value(row._source) << "\n";
    if (row._array[0] == 0)
    {
        row._array[1] += extend_len - row._sourceBeginPos; 
    }
    else
    {
        typename seqan::Gaps<String<Dna5>, ArrayGaps>::TArray_ new_array;
        appendValue(new_array, 0);
        appendValue(new_array, extend_len - row._sourceBeginPos);
        insert(row._array, 0, new_array);
    } 
    row._sourceEndPos += extend_len - row._sourceBeginPos;
    row._clippingEndPos += extend_len - row._sourceBeginPos;
    /*
    for (int i = 0; i < length(row._array); i++)
    {
        dout << "seq2" << row._array[i] << "\n";
    }
    */
    return -extend_len + row._sourceBeginPos;
}
/*
 * The function replace row (._source, ._array) from left to @view_pos with continuous gaps 
   in @row and continuous seqs in row2
 * This transformation keeps the orignial view pos of @row1 @row2 same.
 * @row1 @row2 are supposed to be the rows of the alignment
   namely, their unclipped coordinates are of the same size,
   and the clippedBeginPos and clippedEndPos of @row1 @row2 are required to be equivalent
 * @src_cord_0 is the coordinates in @seq of the value(@row2._source[0])
    just used to transform the unclipped coordinates to the coodinates in @seq
 */
int replaceHead_(Row<Align<String<Dna5>,ArrayGaps> >::Type & row1,
                 Row<Align<String<Dna5>,ArrayGaps> >::Type & row2,
                 String<Dna5> & seq,
                 int64_t gaps_len,
                 uint64_t src_cord_0,
                 int64_t & src_shift1,
                 int64_t & src_shift2)
{
    if (gaps_len > 0)
    {
        src_shift1 = 0;
        src_shift2 = 0;
        ///insertGaps at head of row1
        src_shift1 = seqanRemoveHead_(row1);
            //dout << "rh11" << length(row1._array) << gaps_len << "\n";
        if (!empty(row1._array))
        {

            seqan::insertGaps(row1, 0, gaps_len);
            //row1._array[0] += gaps_len;
            //row1._clippingEndPos += gaps_len;
        }
        else
        {
        //    dout << "rh11" << length(row1._array) << "\n";
        } 
        ///insert source at head of row2;
        src_shift2 = seqanRemoveHead_(row2);
        src_cord_0 += src_shift2;
        src_shift2 += seqanExtendHead_(row2, seq, gaps_len, src_cord_0); 
        setClippedBeginPosition(row2, 0);
        //std::cout << "rh1 " << src_shift1 << " " << src_shift2 << "\n";
        //std::cout << "\nrh1 " << gaps_len << " " << row1 << "\n";
        //std::cout << "\nrh2 " << gaps_len << " " << row2 << "\n";
        return 0;
    }
    else
    {
        return 1;
    }
}
/*
 * Simple score
   > 50 (ins or del)
   todo::benchmark and regress the score 
 */
int64_t mergeScore1_(int64_t gaps_sum, int64_t dx, int64_t dy)
{
    int64_t score = 0; 
    int64_t precision = 10000;  //precision = 1/10000 float to int
    if (dx < 50 && dy < 50)
    {
        score += gaps_sum + dx + dy;
    }
    else
    {
        score += gaps_sum;
    }
    score *= precision;
    return score;
}
int64_t mergeScore_(int64_t gaps_sum, int64_t dx, int64_t dy)
{
    return mergeScore1_(gaps_sum, dx, dy);
}
int findBestMerge_(AlignCache & align1,
                  AlignCache & align2,
                  int64_t & min_clip1,
                  int64_t & min_clip2,
                  int64_t & min_gaps_len,
                  int64_t & f_xy_min)
{
    int bit = 20, bit2 = 40;
    int64_t mask = (1ULL << bit) - 1;
    int64_t gaps_sum = 0; //sum of gap in row11 and row12
    int64_t merge_score = 0;
    unsigned l = align2.length() - 1;
    int64_t gaps_sum_const = 
        - align1.getUCView(0) * 2 + align1.getGSrcX(0) + align1.getGSrcY(0)
        + align2.getUCView(l) * 2 - align2.getGSrcX(l) - align2.getGSrcY(l);
    int64_t min_gap_sum = 1LL << 60; //just a large number
    int64_t min_merge_score = 1LL << 60;
    int64_t f_min = 0;
    int64_t f_xy;
    int64_t gaps_len = 0;
    unsigned j = 0;
    //dout << "ma4<<<<<<<<" << "\n";
    for (unsigned i = 0; i < align1.length() - 1; i++)    
    {
        int64_t x1 = align1.getGSrcX(i);
        int64_t y1 = align1.getGSrcY(i);
        for (; j < align2.length(); j++)  
        {
            int64_t x2 = align2.getGSrcX(j);
            int64_t y2 = align2.getGSrcY(j);
            int64_t dx = x2 - x1;
            int64_t dy = y2 - y1;
            //dout << "ma4" << min_gap_sum << i << j << dx << dy << x1 << y1 << x2 << y2 << "\n";
            if ((dx == 0 || dy == 0) && dx >= 0 && dy >= 0)
            {
                int clip1 = align1.getUCView(i);
                int clip2 = align2.getUCView(j);
                gaps_sum =  gaps_sum_const + clip1 * 2  - clip2 * 2 + dx + dy;
                //dout << "ma3" << gaps_sum << gaps_sum_const << clip1 << clip2 << x1 << x2 << y1 << y2 << "\n";
                merge_score = mergeScore_(gaps_sum, dx, dy);
                if (dx == 0)
                {
                    //min_gaps_len = dy;
                    gaps_len = dy;
                    f_xy = 1;
                }
                else if (dy == 0)
                {
                    //min_gaps_len = dx;
                    gaps_len = dx;
                    f_xy = 2;
                }
                if (min_merge_score > merge_score)
                {
                   min_merge_score = merge_score; 
                   min_clip1 = clip1;
                   min_clip2 = clip2;
                   min_gaps_len = gaps_len;
                   f_xy_min = f_xy;
                   f_min = 1;
                   //dout << "ma333" << min_clip2 << min_gaps_len << x1 << x2 << "\n";
                }
            }
            else if (dx > 0 && dy > 0)
            {
                break;
            }
        }
    }
    //dout << "fm" << min_gaps_len << min_clip2 << f_min << "\n";
    return f_min ? 0 : (16 | 1);
}
int createMergedRows_(Row<Align<String<Dna5>,ArrayGaps> >::Type & row11,
                      Row<Align<String<Dna5>,ArrayGaps> >::Type & row12,
                      Row<Align<String<Dna5>,ArrayGaps> >::Type & row21,
                      Row<Align<String<Dna5>,ArrayGaps> >::Type & row22,
                      String<Dna5> & ref,
                      String<Dna5> & read,
                      String<Dna5> & comrev_read,
                      int64_t min_clip1,
                      int64_t min_clip2,
                      int64_t min_gaps_len,
                      int64_t f_xy_min,
                      uint64_t & cord1,
                      uint64_t & cord2)
{
    //print_cord(cord2, "rh1c");
/*
    dout << "min1" << min_clip1 << min_clip2 << min_gaps_len << "\n";
                if (min_gaps_len > 100)
            {
                dout << "shiftm5" << min_clip1 << min_clip2 << min_gaps_len << get_cord_x(cord1) + endPosition(row11) << get_cord_x(cord2) + beginPosition(row21)  << "\n";
            print_cord(cord1, "shiftm5");
            print_cord(cord2, "shiftm5");
            std::cout << "shiftm5" << row11 << "\n";
            std::cout << "shiftm5" << row12 << "\n";
            std::cout << "shiftm5" << row21 << "\n";
            std::cout << "shiftm5" << row22 << "\n";
            }
    print_cord(cord1, "min1");
    print_cord(cord2, "min1");
            */
    setClippedEndPosition(row11, min_clip1);
    setClippedEndPosition(row12, min_clip1); 
    setClippedBeginPosition(row21, min_clip2);
    setClippedBeginPosition(row22, min_clip2);
    int64_t src_shift1 = 0, src_shift2 = 0;
    //<<debug
    String<Dna5> & seqr = !get_cord_strand(cord2) ? read : comrev_read;
    //std::cout << "malign21 " << value(row11._source)[0] << ref[get_cord_x(cord1)] << " " << value(row21._source)[0] << ref[get_cord_x(cord2)] << " " << value(row12._source)[0]<< seqr[get_cord_y(cord1)] << " " << value(row22._source)[0] << seqr[get_cord_y(cord2)] << "\n";
    //>>debug
    if (f_xy_min == 1) // dx == 1 ins
    {
        String<Dna5> & seq = !get_cord_strand(cord2) ? read : comrev_read;
        if (!replaceHead_(row21, row22, seq, min_gaps_len, get_cord_y(cord2), src_shift1, src_shift2))
        
        {
            //<<debug
            //dout << "shiftm1" << min_gaps_len << "\n";
            /*
            if (min_gaps_len > 100)
            {
            print_cord(cord1, "shiftm1");
            print_cord(cord2, "shiftm1");
            std::cout << "shiftm1" << row11 << "\n";
            std::cout << "shiftm1" << row12 << "\n";
            std::cout << "shiftm1" << row21 << "\n";
            std::cout << "shiftm1" << row22 << "\n";
            }
            //>>debug 
            */
            cord2 = shift_cord(cord2, src_shift1, src_shift2);
        }
    }
    else // dy == 1 del
    {
        String<Dna5> & seq = ref;
        //<<debug
            //std::cout << "malign22 " << value(row11._source)[0] << ref[get_cord_x(cord1)] << " " << value(row21._source)[0] << ref[get_cord_x(cord2)]<< "\n";
        //>>debug
        /*
       if (min_gaps_len > 100)
            {
                dout << "shiftm4" << min_clip1 << min_clip2 << min_gaps_len << get_cord_x(cord1) + endPosition(row11) << get_cord_x(cord2) + beginPosition(row21) << "\n";
            print_cord(cord1, "shiftm4");
            print_cord(cord2, "shiftm4");
            std::cout << "shiftm4" << row11 << "\n";
            std::cout << "shiftm4" << row12 << "\n";
            std::cout << "shiftm4" << row21 << "\n";
            std::cout << "shiftm4" << row22 << "\n";
            } 
*/
        if(!replaceHead_(row22, row21, seq, min_gaps_len, get_cord_x(cord2), src_shift2, src_shift1))
        {
            /*
                        //<<debug
            //dout << "shiftm1" << min_gaps_len << "\n";
            if (min_gaps_len > 100)
            {
                dout << "shiftm2" << min_clip1 << min_clip2 << min_gaps_len << get_cord_x(cord1) + endPosition(row11) << get_cord_x(cord2) + beginPosition(row21)  << "\n";
            print_cord(cord1, "shiftm2");
            print_cord(cord2, "shiftm2");
            std::cout << "shiftm2" << row11 << "\n";
            std::cout << "shiftm2" << row12 << "\n";
            std::cout << "shiftm2" << row21 << "\n";
            std::cout << "shiftm2" << row22 << "\n";
            }
            //>>debug 
            */
            cord2 = shift_cord(cord2, src_shift1, src_shift2);
            //dout << "shiftm22" << length(seq) << src_shift1 << src_shift2 << get_cord_x(cord1) + endPosition(row11) << "\n";
        }
    }   
    //print_cord(cord2, "rh2c");
    //dout << "fmin1 " << min_gaps_len << " " << beginPosition(row21) << endPosition(row21) << clippedBeginPosition(row21) << clippedEndPosition(row21) << "\n";
    //std::cout << "fmin11" << row21 << "\n";
    //std::cout << "fmin11" << row22 << "\n";
    return 0;
}
/*
 * @cord2 is coordinates of value(@row21._source)[0] and value(@row22._source)[0] 
   in sequences.
   if (only if) mergeAling2_ succeed (return 0) then the
   @cord2 is no longer the coordinates of value(...), since new sequence of length min_gaps_len has been inserted into @row21 and @row22 
 */
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
                 uint64_t & cord2)
{
    int64_t min_clip1;
    int64_t min_clip2;
    int64_t min_gaps_len;
    int64_t f_xy_min;
    int f_flag = findBestMerge_(align1, align2, min_clip1, min_clip2, min_gaps_len, f_xy_min);
    //dout << "ma2x" << get_cord_y(cord2) << f_flag << min_clip1 << min_clip2 << min_gaps_len << "\n";
    if (!f_flag)
    {
        f_flag |= createMergedRows_(row11, row12, row21, row22,  
            ref, read, comrev_read, min_clip1, min_clip2, min_gaps_len, f_xy_min,
            cord1, cord2);
    }
        //std::cout << " merged2 " << get_cord_x(cord1) << " " << f_flag <<" "<< row11 << "\n";
    //dout << "ma2_f" << f_flag << "\n";
    return f_flag;
}
std::pair<int, int> cigar2SeqLen(CigarElement<> & cigar)
{
    std::pair<int, int> seq_len(0, 0);
    if (cigar.operation == 'D')
    {
        seq_len.first += cigar.count;
    }
    else if (cigar.operation == 'I')
    {
        seq_len.second += cigar.count;
    }
    else if (cigar.operation == 'X' || 
             cigar.operation == 'M' || 
             cigar.operation == '=')
    {
        seq_len.first += cigar.count;
        seq_len.second += cigar.count;
    }  
    return seq_len;
}
/*
 * Calculate ref and read len of cigars in [@c_str, @c_end)
 */
std::pair<int, int> cigars2SeqsLen(String<CigarElement<> > & cigars,
                                   unsigned c_str, unsigned c_end)
{
    std::pair<int, int> seqs_len(0, 0);
    for (unsigned i = c_str; i < c_end; i++)
    {
        std::pair<int, int> tmp = cigar2SeqLen(cigars[i]);
        seqs_len.first += tmp.first;
        seqs_len.second += tmp.second;
    }
    return seqs_len;
}