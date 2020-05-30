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
FIOParms::FIOParms()
{
    thd_rcb_xy = 15;
    f_reform_ccs = 0;
    f_print_seq = 0;
}

void print_cords_paf(CordsSetType & cords, 
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
                    if (main_strand_count > (block_len / 2))
                    {
                        main_icon_strand = '-';
                    }
                    else if (main_strand_count == (block_len / 2))
                    {
                        main_icon_strand = get_cord_strand(cords[k][j]) ? '-' : '+';
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

                fflag = 0;
            }
        }
    }
}

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
    
    std::stringstream stream; // #include <sstream> for this
    for (unsigned k = 0; k < length(cords); k++)
    {
        if (!empty(cords[k]))
        {
            for (unsigned j = 1; j < length(cords[k]); j++)
            {
                if (_DefaultHit.isBlockEnd(cords[k][j-1]))
                {
                    //dout << "ifoend\n";
                    unsigned m = j; 
                    int main_strand_count = 0;
                    int block_len = 0;
                    ///>determine the main strand
                    while (m < length(cords[k]) && !_DefaultHit.isBlockEnd(cords[k][m]))
                    {
                        if (_DefaultCord.getCordStrand(cords[k][m]))
                        {
                            main_strand_count++;
                        }
                        block_len++;
                        m++;
                    }
                    if (main_strand_count > (block_len / 2))
                    {
                        main_icon_strand = '-';
                    }
                    else if (main_strand_count == (block_len / 2))
                    {
                        main_icon_strand = get_cord_strand(cords[k][j]) ? '-' : '+';
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
                        stream << "\n";
                    } 
                    stream << "@> "
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
                std::string mark = "| ";
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
                stream << mark  
                       << get_cord_y(cords[k][j]) << " " 
                    << get_cord_x(cords[k][j]) << " " 
                    << d2 << " " 
                    << d1 << " " 
                    << get_cord_x (cords[k][j]) - get_cord_y(cords[k][j]) << "\n";
                //of << stream.str();
                //of << buffer;
                cordCount++;
                fflag = 0;
            }
        }
    } 
    of << stream.str();
}

/*=================================================
=                print SAM records                =
===================================================*/
//infer seq bases (10th column) of SAM according to the given cigar element
void cigar2SamSeq(CigarElement<> & cigar, IupacString & result, 
    Iterator<String<Dna5> >::Type & it1, Iterator<String<Dna5> >::Type & it2)
{
    if (cigar.operation == 'D')
    {
        it1 += cigar.count;
    }
    else if (cigar.operation == 'I')
    {
        for (int i = 0; i < cigar.count; i++)
        {
            //appendValue(result, 'N');
            appendValue(result, *it2);
            it2++;
        }
    }
    else if (cigar.operation == 'M')
    {
        for (int i = 0; i < cigar.count; i++)
        {
            appendValue(result, *it1);
            it1++;
            it2++;
        }
    }
    else if (cigar.operation == '=')
    {
        for (int i = 0; i < cigar.count; i++)
        {
            appendValue(result, *it1);
            it1++;
            it2++;
        }
    }
    else if (cigar.operation == 'X')
    {
        for (int i = 0; i < cigar.count; i++)
        {
            appendValue(result, 'N');
            it1++;
            it2++;
        }
    }
    else if (cigar.operation == 'S')
    {
        for (int i = 0; i < cigar.count; i++)
        {
            appendValue(result, *it2);
            it2++;
        }
    }
    else if (cigar.operation == 'H')
    {
        it2 += cigar.count;
    }
}

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
              BamAlignmentRecord & record,
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
             String<BamAlignmentRecordLink>  & records,
             String<Dna5> & genome,
             String<Dna5> & read,
             int & it,
             CharString genome_id,
             CharString genome_id_next,
             FIOParms & fio_parms)
{
    int it_count = -1;
    BamAlignmentRecordLink & record = records[it];
    if (fio_parms.f_print_seq)
    {
        reserve(record.seq, length(read));
    }
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
    //initiate it1 it2 to print seq 
    Iterator<String<Dna5> >::Type it1 = begin(genome) + record.beginPos + 1; 
    Iterator<String<Dna5> >::Type it2 = begin(read); 
    String<Dna5> comp_reverse;
    if (int(record.flag) & int(16))
    {
        _compltRvseStr(read, comp_reverse);
        it2 = begin(comp_reverse);
    }
    if (empty(record.cigar))
        writeValue(target, '*');
    else
    {
        while (1)
        {
            for (unsigned i = 0; i < length(records[it].cigar); ++i)
            {
                if (fio_parms.f_print_seq)
                {
                    cigar2SamSeq(records[it].cigar[i], record.seq, it1, it2);
                }
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
        writeValue(target, '*');  
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
                             std::ofstream & of,
                             FIOParms & fio_parms
                            )
{
    of << "@HD\tVN:1.6\n";
    for (int k = 0; k < length(genomesId); k++)
    {
        of << "@SQ\tSN:" << genomesId[k] << "\tLN:" << length(genomes[k]) << "\n";
    }
    of << "@PG\tPN:" << "Linear\n";
    if (fio_parms.read_group != "" && fio_parms.sample_name != "")
    {
        std::string sample_name = fio_parms.sample_name == "" ? "sample1" : fio_parms.sample_name;
        of << "@RG\t ID:" << fio_parms.read_group << "\tSM:" << sample_name << "\n";
    }
    return 0;
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
            writeSam(of, records[i][j], g_id, "*");
        }
    }
    return 0;
}
int print_align_sam_record_(StringSet<String<BamAlignmentRecordLink> > & records, 
                            StringSet<String<Dna5> > & genomes,
                            StringSet<String<Dna5> > & reads,
                            StringSet<CharString> & genomesId,
                            StringSet<CharString> & readsId, 
                            std::ofstream & of,
                            FIOParms & fio_parms
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
            CharString gnext = "*";
            int dt = writeSam(of, records[i], genomes[records[i][j].rID], reads[i], j, g_id, gnext, fio_parms);
        }
    }
    return 0;
}

int print_align_sam (StringSet<String<Dna5> > & genms,
                     StringSet<String<Dna5> > & reads,
                     StringSet<CharString> & genmsId,
                     StringSet<CharString> & readsId,
                     StringSet<String<BamAlignmentRecordLink> > & bam_records,
                     std::ofstream & of,
                     int f_header,
                     FIOParms & fio_parms)
{
    if (f_header)
    {
        print_align_sam_header_(genmsId, genms, of, fio_parms);
    }
    print_align_sam_record_(bam_records, genms, reads, genmsId, readsId, of, fio_parms); 
    return 0;
}

/*----------  Convert Cords to Bam  ----------*/

/*
 *shortcut to append cigar
 */
void appendCigar(String<CigarElement< > > & cigars, CigarElement<> & cigar)
{
    appendValue(cigars, CigarElement<>(cigar.operation, cigar.count));
}
void appendCigarShrink(String<CigarElement< > > & cigars, CigarElement<> & cigar)
{
    if (!empty(cigars) && back(cigars).operation == cigar.operation)
    {
        back(cigars).count += cigar.count;
    }
    else 
    {
        appendCigar(cigars, cigar);
    }
}

/*
 * If necessary to create new bam record given the @cords
 */
int ifCreateNew_(uint64_t cord1_str, uint64_t cord1_end, uint64_t cord2_str, uint64_t cord2_end, uint64_t thd_large_X)
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
               //     (x12 > x21 && y12 < y21) ||
               //     (x12 < x21 && y12 > x21) ||
                (int64_t(x21 - x12) > int64_t(thd_large_X) && int64_t(y21 - y12) > int64_t(thd_large_X)) ||
                    get_cord_strand (cord1_str ^ cord2_str);
    //dout << "ifnew" << x11 << y11 << x12 << y12 << x21 << y21 << flag << "\n";
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
 * Find cigar from cord1_str to cord2_str
 * NOTE::@cords are required to satisfy the conditons declared at function of ifCreateNew_()
 * NOTE::Don't change the return value, since it's (internal) called by others.
 * Checking cord2_str > cord1_str before calling is required
 */
uint64_t cord2cigar_ (uint64_t cigar_str, //coordinates where the first cigar starts 
                      uint64_t cord1_str, 
                      uint64_t cord1_end,
                      uint64_t cord2_str, 
                      String<CigarElement<> > & cigar)
{
    CigarElement<> cigar1, cigar2;
    uint64_t next_cigar_str;
    uint64_t x0 = get_cord_x (cigar_str);
    uint64_t y0 = get_cord_y (cigar_str);
    uint64_t x11 = get_cord_x(cord1_str);
    uint64_t y11 = get_cord_y(cord1_str);
    uint64_t x12 = get_cord_x(cord1_end);
    uint64_t y12 = get_cord_y(cord1_end);
    uint64_t x21 = get_cord_x(cord2_str);
    uint64_t y21 = get_cord_y(cord2_str);

    if (x0 - y0 != x11 - y11) 
    {
        return ~0; //error
    }
    uint64_t dx = x21 - x0;
    uint64_t dy = y21 - y0;
    //uint64_t len1 = std::min(x12 - x0, y12 - y0); //'=' len upper bound
    //uint64_t len2 = std::min(x21 - x0, y21 - y0); //'=' + 'X' len
    //if (len2 <= len1)
    if (x12 >= x21  && y12 >= y21)
    {
        createRectangleCigarPair(cord1_str, cord2_str, cigar1, cigar2, 0); //'='
        if (cigar1.count) appendCigarShrink(cigar, cigar1);
        if (cigar2.count) appendCigarShrink(cigar, cigar2);
        //appendCigar (cigar, cigar1);
        //appendCigar (cigar, cigar2);
        next_cigar_str = cord2_str;
    }
    else if (x12 < x21 && y12 < y21)
    {
        createRectangleCigarPair(cord1_str, cord1_end, cigar1, cigar2, 0); //'='
        if (cigar1.count) appendCigarShrink(cigar, cigar1);
        if (cigar2.count) appendCigarShrink(cigar, cigar2);
        //appendCigar (cigar, cigar1);
        //appendCigar (cigar, cigar2);
        createRectangleCigarPair(cord1_end, cord2_str, cigar1, cigar2, 1); //'X'
        if (cigar1.count) appendCigarShrink(cigar, cigar1);
        if (cigar2.count) appendCigarShrink(cigar, cigar2);
        //appendCigar (cigar, cigar1);
        //appendCigar (cigar, cigar2);
        next_cigar_str = cord2_str;
    }
    else if (x12 >= x21 && y12 < y21)
    {
        createRectangleCigarPair(cord1_str, cord2_str, cigar1, cigar2, 0);
        if (cigar1.count) appendCigarShrink(cigar, cigar1);
        if (cigar2.count) appendCigarShrink(cigar, cigar2);
        next_cigar_str = cord2_str;
    }
    else if (x12 < x21 && y12 >= y21)
    {
        createRectangleCigarPair(cord1_str, cord2_str, cigar1, cigar2, 0);
        if (cigar1.count) appendCigarShrink(cigar, cigar1);
        if (cigar2.count) appendCigarShrink(cigar, cigar2);
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
                   String<BamAlignmentRecordLink> & bam_link_records,
                   String<Dna5> & read,
                   uint64_t thd_large_X)
{
    uint64_t cigar_str;
    uint64_t cord1_str;
    uint64_t cord2_str;
    uint64_t cord1_end;
    int g_id, g_beginPos, r_beginPos, strand;
    int f_soft = 1; //soft clip in cigar;
    int f_new = 1;
    for (int i = 1; i < length(cords_str); i++)
    {
        if (f_new) //initiate a record for new block 
        {
            f_new = 0;
            g_id = get_cord_id(cords_str[i]);
            g_beginPos = get_cord_x(cords_str[i]);
            r_beginPos = get_cord_y(cords_str[i]);
            strand = get_cord_strand (cords_str[i]);
            insertNewBamRecord (bam_link_records, g_id, g_beginPos, r_beginPos, strand, -1, 1);
            cigar_str = cords_str[i];
        }
        if (i == length(cords_str) - 1 ||
            ifCreateNew_ (cords_str[i], cords_end[i], 
                          cords_str[i + 1], cords_end[i + 1],
                          thd_large_X)) // last cord of current block
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
//serial without omp
void cords2BamLink(StringSet<String<uint64_t> > & cords_str, 
                   StringSet<String<uint64_t> > & cords_end,
                   StringSet<String<BamAlignmentRecordLink> > & bam_link_records,
                   StringSet<String<Dna5> > & reads,
                   int thd_cord_size,
                   uint64_t thd_large_X)
{
    String<BamAlignmentRecordLink> emptyRecords;
    int ii = 0; //parallel 0-based i for each thread
    for (int i = 0; i < length(cords_str); i++)
    {
        appendValue(bam_link_records, emptyRecords);
        if (empty(cords_end) || empty(cords_end[i]))
        {
            String<uint64_t> tmp_end; 
            for (int j = 0; j < length(cords_str[i]); j++)
            {
                appendValue(tmp_end, shift_cord(cords_str[i][j], thd_cord_size, thd_cord_size));
            }
            cords2BamLink (cords_str[i], tmp_end, bam_link_records[ii], reads[i], thd_large_X);
        }
        else
        {
            cords2BamLink (cords_str[i], cords_end[i], bam_link_records[ii], reads[i], thd_large_X);
        }
        ii++;
    }
}

void cords2BamLink(StringSet<String<uint64_t> > & cords_str, 
                   StringSet<String<uint64_t> > & cords_end,
                   StringSet<String<BamAlignmentRecordLink> > & bam_link_records,
                   StringSet<String<Dna5> > & reads,
                   int thd_cord_size,
                   uint64_t thd_large_X,
                   unsigned threads)
{
    #pragma omp parallel
    {
        StringSet<String<BamAlignmentRecordLink> >  bam_link_records_tmp;
        String<BamAlignmentRecordLink> emptyRecords;
        int ii = 0; //parallel 0-based i for each thread
        #pragma omp for
        for (int i = 0; i < length(cords_str); i++)
        {
            appendValue(bam_link_records_tmp, emptyRecords);
            if (empty(cords_end) || empty(cords_end[i]))
            {
                String<uint64_t> tmp_end; 
                for (int j = 0; j < length(cords_str[i]); j++)
                {
                    appendValue(tmp_end, shift_cord(cords_str[i][j], thd_cord_size, thd_cord_size));
                }
                cords2BamLink (cords_str[i], tmp_end, bam_link_records_tmp[ii], reads[i], thd_large_X);
            }
            else
            {
                cords2BamLink (cords_str[i], cords_end[i], bam_link_records_tmp[ii], reads[i], thd_large_X);
            }
            ii++;
        }
        #pragma omp for ordered
        for (unsigned i = 0; i < threads; i++)
        #pragma omp ordered
        {
            append(bam_link_records, bam_link_records_tmp);
        }
    }
}

void shrink_cords_cigar(String<BamAlignmentRecordLink>  & bam_records)
{
    for (int i = 0; i < length(bam_records); i++)
    {
        int dj = 0;
        for (int j = 1; j < length(bam_records[i].cigar); j++)
        {
            if (bam_records[i].cigar[j - dj - 1].operation == bam_records[i].cigar[j].operation)
            {
                bam_records[i].cigar[j - dj - 1].count += bam_records[i].cigar[j].count;
                dj++;
            }
            else if (bam_records[i].cigar[j].count == 0)
            {
                std::cout << "shrinkdj" << dj << "\n";
                dj++;
            }
            else
            {
                bam_records[i].cigar[j - dj] = bam_records[i].cigar[j];
            } 
        }
        resize (bam_records[i].cigar, length(bam_records[i].cigar) - dj);
    }
}

void shrink_cords_cigar(StringSet<String<BamAlignmentRecordLink> > & bam_records)
{
    for (auto & records : bam_records){shrink_cords_cigar(records);}
}
/*
int trimCordsPair(String<uint64_t> & cords_str, String<uint64_t> & cords_end, FIOParms & fio_parms)
{
    if (length(cords_str) < 3)
    {
        return 0;
    }
    uint64_t pre_x1 = get_cord_x(cords_str[1]);
    uint64_t pre_y1 = get_cord_y(cords_str[1]);
    for (int i = 3; i < length(cords_str); i++)
    {
        uint64_t x1 = get_cord_x(cords_str[i]);
        uint64_t y1 = get_cord_y(cords_str[i]);
        int64_t dist = std::min(int64_t(x1 - pre_x1), int64_t(y1 - pre_y1));
        if (dist < 0)
        {

        }
        pre_x1 = get_cord_x(cords_str[i]);
        pre_y1 = get_cord_y(cords_str[i]);
    }
    return 0;
}
*/

int reformCCSBams(String<BamAlignmentRecordLink> & bam_records, 
                 FIOParms &fio_parms)
{
    unsigned it;
    int f_compress = 0;
    for (unsigned i = 0; i < length(bam_records); i = it + 1)
    {
        it = i;
        int end = 0;
        int x1 = 0, y1 = 0; //original  coordinates
        int x2 = 0, y2 = 0; //revised coordinates
        int xy = 0;
        int f_compress = 0;
        int compress_count = 0;
        char compress_operation;
        //dout << "rcs2" << length(bam_records) << "\n";
        while (1)
        {
        //dout << "rcs1" << it << "\n";
            int j2 = 0;
            for (unsigned j1 = 0; j1 < length(bam_records[it].cigar); ++j1)
            {
                int new_count = bam_records[it].cigar[j1].count;
                int compress_count = new_count;
                char new_operation = bam_records[it].cigar[j1].operation;
                char compress_operation = new_operation;
                std::cout << "rcs3" << j1 << " " << new_count << " " << new_operation << "\n";
                if (new_operation == 'I')
                {
                    if (std::abs(xy + new_count) < fio_parms.thd_rcb_xy)
                    {
                        xy += new_count;
                        f_compress = 1;
                        compress_operation = '=';
                        compress_count = new_count;
                    }
                }
                else if (new_operation == 'D')
                {//deletion is ommited
                    if (std::abs(xy - new_count) < fio_parms.thd_rcb_xy)
                    {//allow to compress (replace D with =)
                        xy -= new_count;
                        f_compress = 1;
                        compress_operation = '=';
                        compress_count = 0; 
                    }
                }
                //std::cout << "rcs4" << j1 << " " << j2 << " " << bam_records[it].cigar[j2 - 1].count<< bam_records[it].cigar[j2 - 1].operation << " " << new_count << new_operation << " " << compress_count<< compress_operation << " " << xy << "\n";
                if (j2 > 0 && bam_records[it].cigar[j2 - 1].operation == compress_operation)
                {//merge to last operation
                    bam_records[it].cigar[j2 - 1].count += compress_count;
                }
                else if (compress_count != 0)
                {
                    bam_records[it].cigar[j2].operation = compress_operation;
                    bam_records[it].cigar[j2].count = compress_count;
                    j2++;
                }
                f_compress = 0;
            }
            resize(bam_records[it].cigar, j2);
            if (bam_records[it].isEnd())
            {
                break;
            }
            else
            {
                it = bam_records[it].next();
            }
        }
    }
    return 0;
}

int reformCCSBams(StringSet<String<BamAlignmentRecordLink> > & bam_records, 
                  FIOParms & fio_parms)
{
    for (int i = 0; i < length(bam_records); i++)
    {
        reformCCSBams(bam_records[i], fio_parms);
    }
    return 0;
}
/*
 * Record containing operation 'X' > @thd_large_X is clipped as two records
 */
void print_cords_sam
    (StringSet<String<uint64_t> > & cordset_str,    
     StringSet<String<uint64_t> > & cordset_end,    
     StringSet<String<BamAlignmentRecordLink> > & bam_records,
     StringSet<CharString> & genmsId, 
     StringSet<CharString> & readsId,
     StringSet<String<Dna5> > & genms,
     StringSet<String<Dna5> > & reads,
     int thd_cord_size,
     std::ofstream & of,
     uint64_t thd_large_X,
     unsigned threads,
     int f_header,
     FIOParms & fio_parms,
     int f_parallel)
{
    //thd_large_X = 500;
    if (f_parallel)
    {
        cords2BamLink (cordset_str, cordset_end, bam_records, reads, thd_cord_size, thd_large_X, threads);
    }
    else
    {
        cords2BamLink (cordset_str, cordset_end, bam_records, reads, thd_cord_size, thd_large_X);
    }
    //shrink_cords_cigar(bam_records);
    //fio_parms.f_print_seq = 1;
    if (fio_parms.f_reform_ccs)
    {
        reformCCSBams(bam_records, fio_parms);
    }
    print_align_sam (genms, reads, genmsId, readsId, bam_records, of, f_header, fio_parms);
}   
