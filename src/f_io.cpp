#include <iostream>
#include <fstream>
#include <seqan/basic.h>
#include "cords.h"
#include "pmpfinder.h"
#include "align_util.h"
#include "f_io.h"

using namespace seqan;
using std::cout;
/*============================================================
=               print Approximate mapping records            =
============================================================*/
void print_cords_apf(CordsSetType & cords, 
                     StringSet<String<Dna5> > & genomes,
                     StringSet<String<Dna5> > & reads,
                     StringSet<CharString> & genomesId,
                     StringSet<CharString> & readsId,
                     std::ofstream & of)
{
    unsigned cordCount = 0;
    uint64_t readCordEnd;
    uint64_t seqsCordEnd;
    char main_icon_strand = '+', icon_strand = '+';
    int fflag = 0;
    
    for (unsigned k = 0; k < length(cords); k++)
    {
        if (!empty(cords[k]))
        {
            for (unsigned j = 1; j < length(cords[k]); j++)
            {
                if (_DefaultHit.isBlockEnd(cords[k][j-1]))
                {
                    unsigned m = j; 
                    int main_strand_count = 0;
                    int block_len = 0;
                    ///>determine the main strand
                    while (!_DefaultHit.isBlockEnd(cords[k][m]))
                    {
                        if (_DefaultCord.getCordStrand(cords[k][m]))
                        {
                            main_strand_count++;
                        }
                        block_len++;
                        m++;
                    }
                    if (main_strand_count > (block_len >> 1))
                    {
                        main_icon_strand = '-';
                    }
                    else
                    {
                        main_icon_strand = '+';
                    }
                    ///>print the header
                    for (unsigned i = j; ; i++)
                    {
                        if (_DefaultHit.isBlockEnd(cords[k][i]) || i == length(cords[k]) - 1)
                        {
                            readCordEnd = get_cord_y(cords[k][i]) + window_size;
                            seqsCordEnd = get_cord_x(cords[k][i]) + window_size;
                            break;
                        }
                    }
                    if (k > 0)
                    {
                        of << "\n";
                    } 
                    of << "@> "
                       << readsId[k] << " " 
                       << length(reads[k]) << " "
                       << get_cord_y(cords[k][j]) << " " 
                       << std::min(readCordEnd, (uint64_t)length(reads[k])) << " " 
                       << main_icon_strand<< " "
                       << genomesId[get_cord_id(cords[k][j])] << " " 
                       << length(genomes[get_cord_id(cords[k][j])]) << " "
                       << get_cord_x(cords[k][j]) << " " 
                       << seqsCordEnd << "\n";
                    cordCount = 0;
                    fflag = 1;
                }
                ///>print the coordinates
                icon_strand = (get_cord_strand(cords[k][j]))?'-':'+';
                CharString mark = "| ";
                if (icon_strand != main_icon_strand)
                {
                    mark = (icon_strand == '+') ? "|**+++++++++++ " :"|**----------- ";
                }
                int64_t d1 = 0;//_DefaultCord.getCordY(cords[k][1]);
                int64_t d2 = 0;
                if (!fflag)
                {
                    d1 = int64_t(get_cord_x(cords[k][j]) - get_cord_x(cords[k][j - 1]));
                    d2 = int64_t(get_cord_y(cords[k][j]) - get_cord_y(cords[k][j - 1]));
                }
                of << mark  
                   << get_cord_y(cords[k][j]) << " " 
                   << get_cord_x(cords[k][j]) << " " 
                   << d2 << " " 
                   << d1 << " " 
                   << j << " \n";
                cordCount++;
                fflag = 0;
            }
        }
    }
}

/*=================================================
=                print SAM records                =
===================================================*/
std::string getFileName(std::string s, std::string delimiter, uint count) 
{
    size_t pos = 0;
    std::string token;
    uint i = 1;
    while ((pos = s.find(delimiter)) != std::string::npos) {
        token = s.substr(0, pos);
        s.erase(0, pos + delimiter.length());
        if (i > count)
        {
            token;
            return token;
        }
        i++;
    }
    return s;
}


//Lightweight sam function of Seqan::write(bamAlignmentRecord)
void writeSam(std::ofstream & target,
              BamAlignmentRecord const & record,
              CharString genome_id,
              CharString genome_id_next
            )
{
    write(target, record.qName);
    writeValue(target, '\t');

    appendNumber(target, record.flag);
    writeValue(target, '\t');

    write(target, genome_id);

    writeValue(target, '\t');

    SEQAN_ASSERT_EQ((__int32)BamAlignmentRecord::INVALID_POS + 1, (__int32)0);
    appendNumber(target, record.beginPos + 1);

    writeValue(target, '\t');

    appendNumber(target, static_cast<__uint16>(record.mapQ));
    writeValue(target, '\t');

    if (empty(record.cigar))
        writeValue(target, '*');
    else
        for (unsigned i = 0; i < length(record.cigar); ++i)
        {
            appendNumber(target, record.cigar[i].count);
            writeValue(target, record.cigar[i].operation);
        }

    writeValue(target, '\t');

    if (record.rNextId == BamAlignmentRecord::INVALID_REFID)
        writeValue(target, '*');
    else if (record.rID == record.rNextId)
        writeValue(target, '=');
    else
        write(target, genome_id_next);

    writeValue(target, '\t');

    appendNumber(target, record.pNext + 1);

    writeValue(target, '\t');

    if (record.tLen == BamAlignmentRecord::INVALID_LEN)
        writeValue(target, '0');
    else
        appendNumber(target, record.tLen);

    writeValue(target, '\t');

    if (empty(record.seq))
        writeValue(target, '*');  // Case of empty seq string / "*".
    else
        write(target, record.seq);

    writeValue(target, '\t');


    if (empty(record.qual))  // Case of empty quality string / "*".
        writeValue(target, '*');
    else
        write(target, record.qual);

    if (!empty(record.tags))
    {
        writeValue(target, '\t');
        appendTagsBamToSam(target, record.tags);
    }

    writeValue(target, '\n');
}

//Lightweight sam function of Seqan::write(bamAlignmentRecord)
int writeSam(std::ofstream & target,
             String<BamAlignmentRecordLink> const & records,
             int & it,
             CharString genome_id,
             CharString genome_id_next
            )
{
    int it_count = -1;
    BamAlignmentRecordLink const & record = records[it];
    write(target, record.qName);
    writeValue(target, '\t');

    appendNumber(target, record.flag);
    writeValue(target, '\t');

    write(target, genome_id);

    writeValue(target, '\t');

    SEQAN_ASSERT_EQ((__int32)BamAlignmentRecord::INVALID_POS + 1, (__int32)0);
    appendNumber(target, record.beginPos + 1);

    writeValue(target, '\t');

    appendNumber(target, static_cast<__uint16>(record.mapQ));
    writeValue(target, '\t');

    if (empty(record.cigar))
        writeValue(target, '*');
    else
    {
        int end = 0;
        while (1)
        {
            for (unsigned i = 0; i < length(records[it].cigar); ++i)
            {
                switch (records[it].cigar[i].operation)
                {
                    case 'D':
                        end += records[it].cigar[i].count;
                        break;
                    case 'M':
                        end += records[it].cigar[i].count;
                        break;
                }
                /*
                writeValue(target, ' ');
                appendNumber(target, end);
                writeValue(target, ' ');
                */

                appendNumber(target, records[it].cigar[i].count);
                writeValue(target, records[it].cigar[i].operation);
            }
            it_count++;
            if (records[it].isEnd())
            {
                break;
            }
            else 
            {
                it = records[it].next();
            }
        }
        //std::cout << "[]::print_sam " << end << "\n";
    }
    writeValue(target, '\t');

    if (record.rNextId == BamAlignmentRecord::INVALID_REFID)
        writeValue(target, '*');
    else if (record.rID == record.rNextId)
        writeValue(target, '=');
    else
        write(target, genome_id_next);

    writeValue(target, '\t');

    appendNumber(target, record.pNext + 1);

    writeValue(target, '\t');

    if (record.tLen == BamAlignmentRecord::INVALID_LEN)
        writeValue(target, '0');
    else
        appendNumber(target, record.tLen);

    writeValue(target, '\t');

    if (empty(record.seq))
        writeValue(target, '*');  // Case of empty seq string / "*".
    else
        write(target, record.seq);

    writeValue(target, '\t');


    if (empty(record.qual))  // Case of empty quality string / "*".
        writeValue(target, '*');
    else
        write(target, record.qual);

    if (!empty(record.tags))
    {
        writeValue(target, '\t');
        appendTagsBamToSam(target, record.tags);
    }

    writeValue(target, '\n');
    return it_count;
}

std::pair<int, int> countCigar(String<CigarElement<> > & cigar)
{
    int len1 = 0, len2 = 0;
    for (int i = 0; i < length(cigar); i++)
    {
        //std::cerr << cigar[i].operation << cigar[i].count << " ";
        switch (cigar[i].operation)
        {
            case 'D':
                len1 += cigar[i].count;
                break;
            case 'I':
                len2 += cigar[i].count;
                break;
            case '=':
                len1 += cigar[i].count;
                len2 += cigar[i].count;
                break;
            case 'X':
                len1 += cigar[i].count;
                len2 += cigar[i].count;
                break;
            default:
                break;  
        }
    }
    //std::cerr << "\n";
    return std::pair<int,int>(len1, len2);
}

void printRows(Row<Align<String<Dna5>,ArrayGaps> >::Type & row1,
               Row<Align<String<Dna5>,ArrayGaps> >::Type & row2,
               CharString header)
{
    std::cout << "printRows::" << header << "\n";
    int len = std::min (clippedEndPosition(row1) - clippedBeginPosition(row1),
                        clippedEndPosition(row2) - clippedBeginPosition(row2));
    std::string line0, line00, line1, line2, line3, line4;
    int css1 = 0, css2 = 0;
    for (int i = 0; i < len; i++)
    {
        if (row1[i] != '-')
        {
            css1++;
        }
        if (row2[i] != '-')
        {
            css2++;
        }
        if (i % 10 == 9)
        {
            append (line1, ":");
        }
        else if (i % 5 == 4)
        {
            append (line1, ".");
        }
        else
        {
            append (line1, " ");
        }
        appendValue(line2, row1[i]);
        if (row1[i] == row2[i])
        {
            append (line3, "|");
        }
        else
        {
            append (line3, " ");
        }
        appendValue(line4, row2[i]);
        if (i % 50 == 49 || i == len - 1)
        {
            line1 += "  " + std::to_string(i + 1) + " " + std::to_string(css1) + " " + std::to_string(css2);
            std::cout << line1 << "\n" << line2 << "\n" << line3 << "\n" << line4 << "\n\n";
            line1.clear();
            line2.clear();
            line3.clear();
            line4.clear();
        }
    }
    clear(line0);
    clear(line00);
    std::cout << "\n";
}

int print_align_sam_header_ (StringSet<CharString> & genomesId,
                             StringSet<String<Dna5> > & genomes,
                             std::ofstream & of
                            )
{
    of << "@HD\tVN:1.6\n";
    for (int k = 0; k < length(genomesId); k++)
    {
        of << "@SQ\tSN:" << genomesId[k] << "\tLN:" << length(genomes[k]) << "\n";
    }
    of << "@PG\tPN:" << "Linear\n";
}

int print_align_sam_record_(StringSet<String<BamAlignmentRecord > > & records, 
                            StringSet<CharString> & genomesId,
                            StringSet<CharString> & readsId, 
                            std::ofstream & of
                            )
{
    for (int i = 0; i < length(records); i++)
    {
        for (int j = 0; j < length(records[i]); j++)
        {
            records[i][j].qName = readsId[i];
            CharString g_id = genomesId[records[i][j].rID];
            writeSam(of, records[i][j], g_id);
        }
    }
}
int print_align_sam_record_(StringSet<String<BamAlignmentRecordLink> > & records, 
                            StringSet<CharString> & genomesId,
                            StringSet<CharString> & readsId, 
                            std::ofstream & of
                            )
{
    for (int i = 0; i < length(records); i++)
    {
        for (int j = 0; j < length(records[i]); j++)
        {
            records[i][j].qName = readsId[i];
            CharString g_id = genomesId[records[i][j].rID];
            if (length(records[i][j].cigar) == 0 ||
                ((length(records[i][j].cigar) == 1) && (records[i][j].cigar[0].operation == 'S' || records[i][j].cigar[0].operation == 'H')))
            {
                continue;
            }
            int dt = writeSam(of, records[i], j, g_id);
        }
    }
}

int print_align_sam (StringSet<String<Dna5> > & genms,
                     StringSet<CharString> & genmsId,
                     StringSet<CharString> & readsId,
                     StringSet<String<BamAlignmentRecordLink> > & bam_records,
                     std::ofstream & of
                     )
{
    print_align_sam_header_(genmsId, 
                            genms,
                            of);
    print_align_sam_record_(bam_records,
                            genmsId,
                            readsId,
                            of); 
    return 0;
}

/*----------  Convert Cords to Bam  ----------*/

/*
 *shortcut to append cigar
 */
void appendCigar(String<CigarElement< > > & cigars, char ops, int opn)
{
    appendValue(cigars, CigarElement<>(ops, opn));
}
void appendCigar(String<CigarElement< > > & cigars, CigarElement<> cigar)
{
    appendCigar(cigars, cigar.operation, cigar.count);
}
/*
 * If need to create new bam record given the @cords
 */
int ifCreateNew_(uint64_t cord1_str, uint64_t cord1_end, uint64_t cord2_str, uint64_t cord2_end)
{
    uint64_t x11 = get_cord_x(cord1_str);
    uint64_t y11 = get_cord_y(cord1_str);
    uint64_t x12 = get_cord_x(cord1_end);
    uint64_t y12 = get_cord_y(cord1_end);
    uint64_t x21 = get_cord_x(cord2_str);
    uint64_t y21 = get_cord_y(cord2_str);
    (void) cord2_end;
    int flag = is_cord_block_end (cord1_str) ||  
                             (x11 > x21) || 
                             (y11 > y21) ||
                (x12 > x21 && y12 < y21) ||
                (x12 < x21 && y12 > x21) ||
                get_cord_strand (cord1_str ^ cord2_str);
    return flag;
}
/*
 * create a cigar pair start from cord1 ends at cord2 as *=*D or *=*I
 * @f_m = 0 use '=', else use 'X' cigar
 */
void createRectangleCigarPair (uint64_t cord1, uint64_t cord2, 
                               CigarElement<> & cigar1,
                               CigarElement<> & cigar2,
                               int f_m)
{
    uint64_t dx = get_cord_x (cord2 - cord1);
    uint64_t dy = get_cord_y (cord2 - cord1);
    cigar1.operation = (!f_m) ? '=' : 'X'; 
    if (dx >= dy) 
    {
        cigar2.operation = 'D';
        cigar1.count = dy;
        cigar2.count = dx - dy;
    }
    else 
    {
        cigar2.operation = 'I';
        cigar1.count = dx;
        cigar2.count = dy - dx;
    }
}

/*
 * NOTE::@cords are required to meet the conditons declared at function of ifCreateNew_()
 */
uint64_t cord2cigar_ (uint64_t cigar_str, //coordinates where the first cigar starts 
                      uint64_t cord1_str, 
                      uint64_t cord1_end,
                      uint64_t cord2_str, 
                      String<CigarElement<> > & cigar)
{
    uint64_t x0 = get_cord_x (cigar_str);
    uint64_t y0 = get_cord_y (cigar_str);
    uint64_t x11 = get_cord_x(cord1_str);
    uint64_t y11 = get_cord_y(cord1_str);
    uint64_t x12 = get_cord_x(cord1_end);
    uint64_t y12 = get_cord_y(cord1_end);
    uint64_t x21 = get_cord_x(cord2_str);
    uint64_t y21 = get_cord_y(cord2_str);

    CigarElement<> cigar1, cigar2;
    uint64_t next_cigar_str;
    int opn;
    char ops;
    if (x0 - y0 != x11 - y11) 
    {
        return ~0; //return error
    }

    uint64_t mstrx = get_cord_x(cigar_str); //match ('M'='X' + '=') start
    uint64_t mstry = get_cord_y(cigar_str); 
    uint64_t dx = x21 - mstrx;
    uint64_t dy = y21 - mstry;
    int e_upper = std::min(x12 - mstrx, y12 - mstry); //'=' len upper bound
    uint64_t m_len = std::min(dx, dy); //'=' + 'X' len
    if (m_len <= e_upper)
    {
        createRectangleCigarPair(cord1_str, cord2_str, cigar1, cigar2, 0); //'='
        appendCigar (cigar, cigar1);
        appendCigar (cigar, cigar2);
        next_cigar_str = cord2_str;
    }
    else
    {
        createRectangleCigarPair(cord1_str, cord1_end, cigar1, cigar2, 0); //'='
        appendCigar (cigar, cigar1);
        appendCigar (cigar, cigar2);
        createRectangleCigarPair(cord1_end, cord2_str, cigar1, cigar2, 1); //'X'
        appendCigar (cigar, cigar1);
        appendCigar (cigar, cigar2);
        next_cigar_str = cord2_str;

    }
    //dout << "next_cigar_str" << get_cord_y(next_cigar_str) << get_cord_y(cord2_str) << "\n";
    return next_cigar_str;
}

/*
 *  Function to convert cords to bam
 *  WARN::The @cords_str[0] and back(@cords_str) are required to have block end sign
    Otherwise will cause seg fault. 
 *  NOTE::addjacent cords, cord1 and cord2, will be break into different bams if cord1y > cord2y || cord1x >
    cord2x
 */
void cords2BamLink(String<uint64_t> & cords_str, 
                   String<uint64_t> & cords_end,
                   String<BamAlignmentRecordLink> & bam_link_records)
{
    uint64_t cigar_str;
    uint64_t cord1_str;
    uint64_t cord2_str;
    uint64_t cord1_end;
    int f_soft = 1; //soft clip in cigar;
    int f_new = 1;
    for (int i = 1; i < length(cords_str); i++)
    {
        if (f_new) //initiate a record for new block 
        {
            f_new = 0;
            int g_id = get_cord_id(cords_str[i]);
            int g_beginPos = get_cord_x(cords_str[i]);
            int r_beginPos = get_cord_y(cords_str[i]);
            int strand = get_cord_strand (cords_str[i]);
            insertNewBamRecord (bam_link_records, g_id, g_beginPos, r_beginPos, strand);
            cigar_str = cords_str[i];
        }
        if (i == length(cords_str) - 1 ||
            ifCreateNew_ (cords_str[i], cords_end[i], 
                          cords_str[i + 1], cords_end[i + 1])) // last cord of current block
        {
            cord1_str = cords_str[i];
            cord1_end = cords_end[i];
            cord2_str = cords_end[i];
            f_new = 1; //next cord[i + 1] will start a new recordd
        }
        else
        {
            cord1_str = cords_str[i];
            cord1_end = cords_end[i];
            cord2_str = cords_str[i + 1];
        }
        cigar_str = cord2cigar_ (cigar_str, 
                                 cord1_str, cord1_end, cord2_str, 
                                 back(bam_link_records).cigar);
        if (cigar_str == ~0) //error
        {
            break;
        }
    }
}

void cords2BamLink(StringSet<String<uint64_t> > & cords_str, 
                   StringSet<String<uint64_t> > & cords_end,
                   StringSet<String<BamAlignmentRecordLink> > & bam_link_records,
                   int thd_cord_size)
{
    if (empty(bam_link_records))
    {
        resize (bam_link_records, length(cords_str));
    }
    for (int i = 0; i < length(cords_str); i++)
    {
        if (empty(cords_end) || empty(cords_end[i]))
        {
            String<uint64_t> tmp_end; 
            for (int j = 0; j < length(cords_str[i]); j++)
            {
                appendValue(tmp_end, shift_cord(cords_str[i][j], thd_cord_size, thd_cord_size));
            }
            cords2BamLink (cords_str[i], tmp_end, bam_link_records[i]);
        }
        else
        {
            cords2BamLink (cords_str[i], cords_end[i], bam_link_records[i]);
        }
    }
}

void print_cords_sam
    (StringSet<String<uint64_t> > & cordset_str,    
     StringSet<String<uint64_t> > & cordset_end,    
     StringSet<String<BamAlignmentRecordLink> > & bam_records,
     StringSet<CharString> & genmsId, 
     StringSet<CharString> & readsId,
     StringSet<String<Dna5> > & genms,
     int thd_cord_size,
     std::ofstream & of)
{
    cords2BamLink (cordset_str, cordset_end, bam_records, thd_cord_size);
    print_align_sam (genms, genmsId, readsId, bam_records, of);
}   