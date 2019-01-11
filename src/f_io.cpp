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


int align2cigar_(Align<String<Dna5>,ArrayGaps> & align, 
                 std::string & cigar, 
                 std::string & mutations)
                 //int source_start = 0; 
                 //int source_end = toSourcePosition(length(row(align, 0))))
{

    typedef Align<String<Dna5>, ArrayGaps> TAlign;
    typedef typename Source<TAlign>::Type TSource;
    typedef typename Iterator<TSource, Rooted>::Type TStringIterator;

    typedef typename Row<TAlign>::Type TRow;
    typedef typename Iterator<TRow, Rooted>::Type TAlignIterator;

    TAlignIterator ali_it0_stop = iter(row(align, 0), endPosition(cols(align)));
    TAlignIterator ali_it1_stop = iter(row(align, 1), endPosition(cols(align)));
    TAlignIterator ali_it0 = iter(row(align, 0), beginPosition(cols(align)));
    TAlignIterator ali_it1 = iter(row(align, 1), beginPosition(cols(align)));
    TStringIterator readBase = begin(source(row(align, 0)));
    
    int readPos = 0; 
    bool first = true;
    while (ali_it0 != ali_it0_stop && ali_it1 != ali_it1_stop)
    {    
        int matched = 0; 
        int inserted = 0; 
        int deleted = 0; 
        while (ali_it0 != ali_it0_stop && ali_it1 != ali_it1_stop && !isGap(ali_it0) && !isGap(ali_it1))
        {
            ++readPos;
            if (*ali_it1 != *ali_it0)
            {
                if (first)
                    first = false;
                else
                    mutations << ",";
                mutations << readPos << (char)*readBase;
            }
            ++readBase;
            ++ali_it0;
            ++ali_it1;
            ++matched;
        }
        if (matched > 0)
        {std::to_string(matched);
            cigar << matched << "M";
        }
        while (ali_it0 != ali_it0_stop && isGap(ali_it0))
        {
            ++ali_it0;
            ++ali_it1;
            ++deleted;
        }
        if (deleted > 0)
            cigar << deleted << "D";
        while (isGap(ali_it1) && ali_it1 != ali_it1_stop)
        {
            ++ali_it0;
            ++ali_it1;
            ++readPos;
            if (first)
                first = false;
            else
                mutations << ",";
            mutations << readPos << (char)*readBase;
            ++readBase;
            ++inserted;
        }
        if (inserted > 0)
            cigar << inserted << "I";
    }
    return 0;
}
int align2cigar(Align<String<Dna5>,ArrayGaps> & align, 
                 std::string & cigar, 
                 std::string & mutations)
{
    align2cigar_(align, cigar, mutations);
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
//                op = ((TVal1)*it1 == (TVal2)*it2)? '=': 'X';
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