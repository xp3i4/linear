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

///print one gff record
/*
int print_gff_(GffFileOut & out, 
               CharString ref,
               CharString source,
               CharString type,
               uint64_t beginPos,
               uint64_t endPos,
               uint64_t strand,
               int score,
               StringSet<CharString> tagNames, 
               StringSet<CharString> tagValues
              )
{
    GffRecord record;
    record.ref = ref;
    record.source = source;
    record.type = type;
    record.beginPos = beginPos;
    record.endPos = endPos;
    record.strand = strand;
    record.score = score;
    append(record.tagNames, tagNames);
    append(record.tagValues, tagValues);
    writeRecord(out, record);
    return 0;
}
*/


/**
 *concat two blocks specified by 'start' and 'end' to cigar according to align1 and align2
 * NOTE:j 
 *
int aligns2cigar_(Align<String<Dna5>,ArrayGaps> & align1, 
                 Align<String<Dna5>,ArrayGaps> & align2,
                 uint64_t start1,
                 uint64_t start2,
                 std::string & cigar, 
                 std::string & mutations)
{
    
    typedef Align<String<Dna5>, ArrayGaps> TAlign;
    typedef Row<TAlign>::Type TRow; 
    TRow & row11 = row(align1, 0);
    TRow & row12 = row(align1, 1);
    
    int x1 = start2;
    int x2 = start1 + toSourcePosition(align1, length(row11) - 1) + 1 // + 1 to get region [,)
    int y1 = y2;
    int y2 = start1 + toSourcePosition(align1, length(row12) - 1) + 1;
    
    if (x1 > x2 || y1 > y2)
    {
        return 1;
    }
    
    

    return 0;
}
*/



