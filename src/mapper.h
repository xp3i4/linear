#ifndef SEQAN_HEADER_PACMAPPER_H
#define SEQAN_HEADER_PACMAPPER_H

#include "base.h"
#include "index_util.h"
#include "f_io.h"
#include "mapparm.h"


class Mapper {
    PMRecord    record;
    MapParm     parm;
    IndexDynamic index_dynamic;
    StringSet<String<uint64_t> >  cordSet;
    std::ofstream of;
    unsigned _thread;
    String<int> rlens;
    StringSet<String<uint64_t> > clip_set;
    StringSet<String<BamAlignmentRecordLink> > bam_records;
    std::string outputPrefix;

public:
    Mapper();
    Mapper(Options & options);
    StringSet<String<Dna5> > & reads() {return record.seq1;}             
    StringSet<String<Dna5> > & genomes() {return record.seq2;}             
    MapParm & mapParm() {return parm;}
    IndexDynamic & index() {return index_dynamic;}
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
