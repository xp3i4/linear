#ifndef LINEAR_HEADER_TEST_UNIT_H
#define LINEAR_HEADER_TEST_UNIT_H

#include <seqan/sequence.h>

using seqan::String;

int check_cigar(StringSet<String<Dna5> > & genomes,
                String<Dna5> & read, 
                String<Dna5> & comrevRead,
                String<uint64_t> & cords, //raw cords
                String<BamAlignmentRecordLink> & bam_records);

int check_index2(StringSet<String<Dna5> > & seqs, IndexDynamic & index, unsigned threads);

#endif