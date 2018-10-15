#ifndef LINEAR_HEADER_F_IO_H
#define LINEAR_HEADER_F_IO_H

using namespace seqan;

std::string & operator<< (std::string & s, int i);
std::string & operator<< (std::string & s, char s2);
std::string & operator<< (std::string & s, std::string s2);

int align2cigar_(Align<String<Dna5>,ArrayGaps> & align, 
                 std::string & cigar, 
                 std::string & mutations);

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



/*
//print one gff record
int print_gff_(GffFileOut & out, 
               String ref,
               String source,
               String type,
               uint64_t beginPos,
               uint64_t endPos,
               uint64_t stand,
               int score,
               String tagName, 
               String tagValue
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
    appendValue(record.tagNames, tagName);
    appendValue(record.tagValues, tagValue);
    writeRecord(out, record);
    return 0;
}
*/
#endif
