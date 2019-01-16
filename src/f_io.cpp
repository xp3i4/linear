#include <seqan/seq_io.h>
#include <seqan/stream.h>
#include <seqan/index.h>
#include <seqan/store.h>
#include <seqan/basic.h>
#include <seqan/arg_parse.h>
#include <seqan/vcf_io.h>
#include <iostream>
#include <fstream>
#include "f_io.h"
using namespace seqan;

std::string getFileName(const std::string& s, char sep, int flag) {

    if (flag == 1)
    {
        size_t i = s.rfind(sep, s.length());
        if (i != std::string::npos) 
        {
            return(s.substr(i+1, s.length() - i));
        }   
    }
    else
    {
        size_t i = s.rfind(sep, s.length());
        if (i != std::string::npos) 
        {
            return(s.substr(0, i));
        }   
    }
    return s;
}

std::string & operator<< (std::string & s, int i)
{
    s += std::to_string(i);
    return s;
}

std::string & operator<< (std::string & s, char s2)
{
    s += s2;
    return s;
}

std::string & operator<< (std::string & s, std::string s2)
{
    s += s2;
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

    char op = '?', lastOp = ' ';
    unsigned numOps = 0;
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
//                op = (*it1 == *it2)? '=': 'X';
                op = 'M';
        }
        if (lastOp != op)
        {
            if (lastOp == 'D' && numOps >= (unsigned)splicedGapThresh)
                lastOp = 'N';
            if (numOps > 0)
                appendValue(cigar, CigarElement<>(lastOp, numOps));
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
}
void align2cigar(String<CigarElement< > > &cigar,
                 Row<Align<String<Dna5>,ArrayGaps> >::Type &gaps1,
                 Row<Align<String<Dna5>,ArrayGaps> >::Type &gaps2
                )
{
    align2cigar_(cigar, gaps1, gaps2, 1000);
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

