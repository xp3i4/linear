#ifndef LINEAR_HEADER_F_IO_H
#define LINEAR_HEADER_F_IO_H

#include <seqan/bam_io.h>
#include "cords.h"
#include "align_util.h"
using namespace seqan;
using std::ofstream;

void print_cords_apf(CordsSetType & cords, 
                     StringSet<String<Dna5> > & genomes,
                     StringSet<String<Dna5> > & reads,
                     StringSet<CharString> & genomesId,
                     StringSet<CharString> & readsId,
                     std::ofstream & of);



std::string & operator<< (std::string & s, int i);
std::string & operator<< (std::string & s, char s2);
std::string & operator<< (std::string & s, std::string s2);

std::string getFileName(std::string, std::string sep = "/", uint flag = ~0);

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

std::pair<int, int> countCigar(String<CigarElement<> > & cigar);
void printRows(Row<Align<String<Dna5>,ArrayGaps> >::Type & row1,
               Row<Align<String<Dna5>,ArrayGaps> >::Type & row2,
               CharString header = ""
               );

int print_align_sam_header_ (StringSet<CharString> & genomesId,
                             StringSet<String<Dna5> > & genomes,
                             std::ofstream & of
                            );

int print_align_sam_record_(StringSet<String<BamAlignmentRecord > > & records, 
                            StringSet<String<uint64_t> > & cordSet,
                            StringSet<CharString> & readsId, 
                            StringSet<CharString> & genomesId,
                            std::ofstream & of
                            );

int print_align_sam_record_(StringSet<String<BamAlignmentRecordLink> > & records, 
                            StringSet<String<uint64_t> > & cordSet,
                            StringSet<CharString> & readsId, 
                            StringSet<CharString> & genomesId,
                            std::ofstream & of
                            );
int print_align_sam (StringSet<String<Dna5> > & genms,
                     StringSet<CharString> & readsId,
                     StringSet<CharString> & genmsId,
                     StringSet<String<BamAlignmentRecordLink> > & bam_records,
                     StringSet<String<uint64_t> > & cordset,
                     std::ofstream & of
                     );
void print_cords_sam
    (StringSet<String<uint64_t> > & cordset,    
     StringSet<String<BamAlignmentRecordLink> > & bam_records,
     StringSet<CharString> & genmsId, 
     StringSet<CharString> & readsId,
     StringSet<String<Dna5> > & genms,
     std::ofstream & of);
void addNextBamLink(String<BamAlignmentRecordLink> & bam_records,
                    int id, int next_id);
void cords2cigar(String<uint64_t> & cords, String<CigarElement< > > &cigar);

#endif
