#include <iostream>
#include <fstream>
#include <seqan/basic.h>
#include "cords.h"
#include "pmpfinder.h"
#include "f_io.h"

using namespace seqan;
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
            std::cout << "movesH" << cigar[0].count << " "<< cigar[0].operation << " " << cigar[1].count << " "<< cigar[1].operation <<"\n";
            //eraseFront(cigar);
            erase(cigar,0);
            std::cout << "movesH1 " << cigar[0].count << " "<< cigar[0].operation << "\n";
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
void align2cigar_(String<CigarElement< > > &cigar,
        Row<Align<String<Dna5>,ArrayGaps> >::Type &gaps1,
        Row<Align<String<Dna5>,ArrayGaps> >::Type &gaps2,
        unsigned splicedGapThresh
        )
{
    typedef Row<Align<String<Dna5>,ArrayGaps> >::Type TRow;
    typename Iterator<TRow>::Type it1 = begin(gaps1);
    typename Iterator<TRow>::Type it2 = begin(gaps2);
    std::cout << "align2cigar_ " << clippedEndPosition(gaps1) - clippedBeginPosition(gaps1) << "  " << clippedEndPosition(gaps2) - clippedBeginPosition(gaps2) << "\n";
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
            if (lastOp == 'D' && numOps >= (unsigned)splicedGapThresh)
                lastOp = 'N';
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
    if (lastOp == 'D' && numOps >= splicedGapThresh)
        lastOp = 'N';
    if (numOps > 0)
        appendValue(cigar, CigarElement<>(op, numOps));
    std::cout << "align2cigar_end " << length(cigar) << "\n";
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
            case 'N':
                break;
            case 'P':
                break;
            default:
                return 2;
        }
        std::cout << "[]::clip_cigar "  << " " << x << " " << y << " " << cigar[i].count << cigar[i].operation << "\n";
    }
    return 0;
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
                break;
            else 
                it = records[it].next();
        }
        std::cout << "[]::print_sam " << end << "\n";
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

int insertCigar(String<CigarElement< > > &cigar1, 
                int pos,
                String<CigarElement< > > &cigar2
         )
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
                            StringSet<String<uint64_t> > & cordSet,
                            StringSet<CharString> & readsId, 
                            StringSet<CharString> & genomesId,
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
                            StringSet<String<uint64_t> > & cordSet,
                            StringSet<CharString> & readsId, 
                            StringSet<CharString> & genomesId,
                            std::ofstream & of
                            )
{
    for (int i = 0; i < length(records); i++)
    {
        for (int j = 0; j < length(records[i]); j++)
        {
            records[i][j].qName = readsId[i];
            CharString g_id = genomesId[records[i][j].rID];
            int dt = writeSam(of, records[i], j, g_id);
        }
    }
}

int print_align_sam (StringSet<String<Dna5> > & genms,
                     StringSet<CharString> & readsId,
                     StringSet<CharString> & genmsId,
                     StringSet<String<BamAlignmentRecordLink> > & bam_records,
                     StringSet<String<uint64_t> > & cordset,
                     std::ofstream & of
                     )
{
    print_align_sam_header_(genmsId, 
                            genms,
                            of);
    print_align_sam_record_(bam_records,
                            cordset,
                            readsId,
                            genmsId,
                            of); 
    return 0;
}

void cords2cigar(StringSet<String<uint64_t> > & cordset,     
                 StringSet<String<BamAlignmentRecordLink> > & bam_records
                 //int thd_window_size
                 )
{/*
    for (uint i = 0; i < length(cords); i++)
    {
        for (int j = 0; j < length(cords[i]); j++) 
        {
            int dx = get_cord_x(cords[j + 1]) - get_cord_x(cords[j]);
            int dy = get_cord_y(cords[j + 1]) - get_cord_y(cords[j]);
            if (dx < 0)
            {

            }
        }
    }
*/
    dout << "cords2cigar\n";
}

void print_cords_sam
    (StringSet<String<uint64_t> > & cordset,    
     StringSet<String<BamAlignmentRecordLink> > & bam_records,
     StringSet<CharString> & genmsId, 
     StringSet<CharString> & readsId,
     StringSet<String<Dna5> > & genms,
     std::ofstream & of
     )
{
    cords2cigar(cordset, bam_records);
    print_align_sam (genms, readsId, genmsId, bam_records, cordset, of);
    of << "Hello sam file\n";
    //print
}