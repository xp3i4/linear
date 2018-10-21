#ifndef LINEAR_HEADER_F_IO_H
#define LINEAR_HEADER_F_IO_H

using namespace seqan;

std::string & operator<< (std::string & s, int i);
std::string & operator<< (std::string & s, char s2);
std::string & operator<< (std::string & s, std::string s2);

int align2cigar_(Align<String<Dna5>,ArrayGaps> & align, 
                 std::string & cigar, 
                 std::string & mutations);
#endif
