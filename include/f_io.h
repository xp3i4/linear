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
                            StringSet<CharString> & readsId, 
                            StringSet<CharString> & genomesId,
                            std::ofstream & of
                            );

int print_align_sam_record_(StringSet<String<BamAlignmentRecordLink> > & records, 
                            StringSet<CharString> & readsId, 
                            StringSet<CharString> & genomesId,
                            std::ofstream & of
                            );
int print_align_sam (StringSet<String<Dna5> > & genms,
                     StringSet<CharString> & readsId,
                     StringSet<CharString> & genmsId,
                     StringSet<String<BamAlignmentRecordLink> > & bam_records,
                     std::ofstream & of
                     );

void addNextBamLink(String<BamAlignmentRecordLink> & bam_records,
                    int id, int next_id);

uint64_t cord2cigar_ (uint64_t cigar_str, //coordinates where the first cigar starts 
                      uint64_t cord1_str, 
                      uint64_t cord1_end,
                      uint64_t cord2_str, 
                      String<CigarElement<> > & cigar,
                      int f_first,
                      int f_soft); //flag of first cords

/*
 *  Function to convert cords to bam
 *  WARN::The @cords_str[0] and back(@cords_str) are required to have block end sign
    Otherwise will cause seg fault. 
 *  NOTE::addjacent cords, cord1 and cord2, will be break into different bams if cord1y > cord2y || cord1x >
    cord2x
 */
void cords2BamLink(String<uint64_t> & cords_str, 
                   String<uint64_t> & cords_end,
                   String<BamAlignmentRecordLink> & bam_link_records);

void cords2BamLink(StringSet<String<uint64_t> > & cords_str, 
                   StringSet<String<uint64_t> > & cords_end,
                   StringSet<String<BamAlignmentRecordLink> > & bam_link_records);

void print_cords_sam (StringSet<String<uint64_t> > & cordset_str,    
                      StringSet<String<uint64_t> > & cordset_end,    
                      StringSet<String<BamAlignmentRecordLink> > & bam_records,
                      StringSet<CharString> & genmsId, 
                      StringSet<CharString> & readsId,
                      StringSet<String<Dna5> > & genms,
                      int thd_cord_size,
                      std::ofstream & of
                      );
#endif
