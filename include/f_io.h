#ifndef LINEAR_HEADER_F_IO_H
#define LINEAR_HEADER_F_IO_H

#include <seqan/bam_io.h>
#include "base.h"
#include "cords.h"
#include "align_util.h"
using namespace seqan;
using std::ofstream;

struct FIOParms
{
    uint64_t thd_large_X; 
    uint64_t thd_trim_bound;
    //reformCCSBams
    int thd_rcb_xy;
    int f_reform_ccs;
    int f_print_seq;
    int f_sequence_sam;
    int f_is_align;
    unsigned f_output_type; //call base.h::fp_hanler_ to set unset and resolve the type

    std::string read_group;
    std::string sample_name;

    BamHeader bam_header;
    BamHeader bam_header2; //supposed to be used by nonstand header like those used in pbsv

    FIOParms();
};

void print_cords_paf(CordsSetType & cords, 
                     StringSet<String<Dna5> > & genomes,
                     StringSet<String<Dna5> > & reads,
                     StringSet<CharString> & genomesId,
                     StringSet<CharString> & readsId,
                     std::ofstream & of);

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

std::pair<int, int> countCigar(String<CigarElement<> > & cigar);
void printRows(Row<Align<String<Dna5>,ArrayGaps> >::Type & row1,
               Row<Align<String<Dna5>,ArrayGaps> >::Type & row2,
               CharString header = ""
               );

int printAlignSamBam (StringSet<String<Dna5> > & genms,
                      StringSet<String<Dna5> > & reads,
                      StringSet<CharString> & genmsId,
                      StringSet<CharString> & readsId,
                      StringSet<String<BamAlignmentRecordLink> > & bam_records,
                      std::ofstream & of,
                      int f_header, //if print header
                      FIOParms & fio_parms);

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
void cords2BamLink(StringSet<String<uint64_t> > & cords_str,
                   StringSet<String<uint64_t> > & cords_end,
                   StringSet<String<CordInfo> > & cords_info,
                   StringSet<String<BamAlignmentRecordLink> > & bam_link_records,
                   StringSet<String<Dna5> > & reads,
                   int thd_cord_size,
                   uint64_t thd_large_X);

void cords2BamLink(StringSet<String<uint64_t> > & cords_str, 
                   StringSet<String<uint64_t> > & cords_end,
                   StringSet<String<CordInfo> > & cords_info,
                   StringSet<String<BamAlignmentRecordLink> > & bam_link_records,
                   StringSet<String<Dna5> > & reads,
                   int thd_cord_size,
                   uint64_t thd_large_X,
                   unsigned threads,
                   int f_parallel);

/*
void print_cords_sam (StringSet<String<uint64_t> > & cordset_str,    
                      StringSet<String<uint64_t> > & cordset_end,    
                      StringSet<String<BamAlignmentRecordLink> > & bam_records,
                      StringSet<CharString> & genmsId, 
                      StringSet<CharString> & readsId,
                      StringSet<String<Dna5> > & genms,
                      StringSet<String<Dna5> >& reads,
                      int thd_cord_size,
                      std::ofstream & of,
                      uint64_t thd_large_X,
                      unsigned threads,
                      int f_header,
                      FIOParms & fio_parms,
                      int f_parallel = 1);
                      */
#endif
