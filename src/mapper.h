#ifndef SEQAN_HEADER_PACMAPPER_H
#define SEQAN_HEADER_PACMAPPER_H

#include "base.h"
#include "index_util.h"
#include "f_io.h"
#include "mapparm.h"

extern int const typeDIx;
extern int const typeHIx;

class Mapper {
    PMRecord    record;
    MapParm     parm;
    DIndex      dIndex;
    LIndex      qIndex;
    StringSet<String<uint64_t> >  cordSet;
    std::ofstream of;
    unsigned _thread;
    String<int> rlens;
    StringSet<String<uint64_t> > clip_set;
    StringSet<String<BamAlignmentRecordLink> > bam_records;
    std::string outputPrefix;
    int typeIx;

public:
    Mapper();
    Mapper(Options & options);
    StringSet<String<Dna5> > & reads() {return record.seq1;}             
    StringSet<String<Dna5> > & genomes() {return record.seq2;}             
    MapParm & mapParm() {return parm;}
    DIndex & index_d(){return dIndex;}
    LIndex & index() {return qIndex;}
    StringSet<String<uint64_t> > & cords() {return cordSet;}            //returns cord set 
    
    void printCords(std::ostream & );
    void printCordsRaw();
    void printCordsRaw2();
    int print_vcf();
    int createIndex(bool = false);
    unsigned sens(){return parm.sensitivity;}
    unsigned & thread(){return _thread;}
    CharString & readPath(){return record.readPath;}
    CharString & genomePath(){return record.genomePath;}
    StringSet<CharString> & readsId(){return record.id1;}
    StringSet<CharString> & genomesId(){return record.id2;}
    StringSet<String<uint64_t> > & getClips(){return clip_set;}
    String<int> & readLens(){return rlens;}
    std::ofstream & getOf() {return of;}
    std::string & getOutputPrefix(){return outputPrefix;}
    StringSet<String<BamAlignmentRecordLink> > & getBamRecords() {return bam_records;}
};

#endif
