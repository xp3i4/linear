#ifndef LINEAR_HEADER_F_IO_H
#define LINEAR_HEADER_F_IO_H

using namespace seqan;

class BamAlignmentRecordLink : public BamAlignmentRecord 
{ // Used as String<BamAlignmentRecordLink>, in which cigars of different records
  // can be concated without modifyfing memory.
public:
    int next_id; //next records id

    BamAlignmentRecordLink();
    void addNext(int id);
    int isEnd() const;
    int next() const;
};


std::string & operator<< (std::string & s, int i);
std::string & operator<< (std::string & s, char s2);
std::string & operator<< (std::string & s, std::string s2);

std::string getFileName(const std::string&, char sep = '/', int flag = 1);
void align2cigar(String<CigarElement< > > &cigar,
                 Row<Align<String<Dna5>,ArrayGaps> >::Type &gaps1,
                 Row<Align<String<Dna5>,ArrayGaps> >::Type &gaps2
                );
void writeSam(std::ofstream & target,
              BamAlignmentRecord const & record,
              CharString genome_id,
              CharString genome_id_next = "*"
             );
int writeSam(std::ofstream & target,
              String<BamAlignmentRecordLink> const & record,
              int & it,
              CharString genome_id,
              CharString genome_id_next = "*"
             );
int clip_cigar (String<CigarElement<> > & cigar);

int insertCigar(String<CigarElement< > > &cigar1, 
                int pos,
                String<CigarElement< > > &cigar2
         ); //insert cirgar2 to cigar1 at pos
std::pair<int, int> countCigar(String<CigarElement<> > & cigar);
void printRows(Row<Align<String<Dna5>,ArrayGaps> >::Type & row1,
               Row<Align<String<Dna5>,ArrayGaps> >::Type & row2,
               int i = -1
               );

#endif
