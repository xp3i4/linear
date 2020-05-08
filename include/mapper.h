#ifndef SEQAN_HEADER_PACMAPPER_H
#define SEQAN_HEADER_PACMAPPER_H

#include "base.h"
#include "index_util.h"
#include "f_io.h"
#include "gap.h"
#include "parallel_io.h"
//#include "mapparm.h"


class Mapper : public P_Mapper 
{
    PMRecord    record;
    MapParms     map_parms;
    GapParms    gap_parms;
    IndexDynamic index_dynamic;
    StringSet<String<uint64_t> >  cordSet;
    StringSet<String<uint64_t> >  cordSet2;
    std::ofstream of;
    unsigned _thread;
    String<size_t> rlens;
    StringSet<String<uint64_t> > clip_set;
    StringSet<String<BamAlignmentRecordLink> > bam_records;
    StringSet<FeaturesDynamic > f2; //features of genomes
    std::string outputPrefix;
    int feature_type;
    int of_type;
    const int OF_APP = 0; //set of std::ios::app
    const int OF_NEW = 1;
    Options::PathsType r_paths;
    Options::PathsType g_paths;
    uint f_map;     //map control flag 
    uint f_print;   //print control flag
    uint64_t gap_len_min; //process gaps of length > this value
    int cord_size; //default cord size 

    //=== pipeline2 of parallel buffer 
    //P_ReadsBuffer reads_buffer;
    //P_ReadsIdsBuffer reads_ids_buffer; 
    //P_CordsBuffer cords_buffer;
    //P_BamLinkBuffer bam_link_buffer;

public:
    Mapper();
    Mapper(Options & options);
    void loadOptions(Options & options);
    StringSet<String<Dna5> > & getReads() {return record.seq1;}             
    StringSet<String<Dna5> > & getGenomes() {return record.seq2;}             
    MapParms & getMapParms() {return map_parms;}
    IndexDynamic & getIndex() {return index_dynamic;}
    StringSet<String<uint64_t> > & getCords() {return cordSet;}             
    StringSet<String<uint64_t> > & getCords2() {return cordSet2;}            
    
    void printCords(std::ostream & );
    void printCordsRaw();
    void printCordsRaw2();
    int  print_vcf();
    int  createIndex(unsigned, unsigned, bool = false);
    uint sens(){return map_parms.sensitivity;}
    uint & getThreads(){return _thread;}
    StringSet<CharString> & getReadsId(){return record.id1;}
    StringSet<CharString> & getGenomesId(){return record.id2;}
    StringSet<String<uint64_t> > & getClips(){return clip_set;}
    String<size_t> & getReadsLen(){return rlens;}
    std::ofstream & getOf() {return of;}
    std::string & getOutputPrefix(){return outputPrefix;}
    StringSet<String<BamAlignmentRecordLink> > & getBamRecords() {return bam_records;}
    StringSet<FeaturesDynamic > & getGenomesFeatures(){return f2;} //features of genomes
    int  getFeatureType();
    void setOfApp();
    void setOfNew();
    bool isOfNew();
    bool isOfApp();
    Options::PathsType & getRPaths();
    Options::PathsType & getGPaths();
    uint & getMapFlag(){return f_map;}
    uint & getPrintFlag(){return f_print;}
    uint getGapLenMin () {return gap_len_min;}
    int  getCordSize() {return cord_size;}
    void loadGenomes();
    void clearIndex();
    GapParms & getGapParms(){return gap_parms;}

    //=== pipeline2 of parallel buffer 
    void initBuffers(int, int, P_Parms & parms);
    int p_calRecords(int, int, int);
    int p_printResults(int, int, int);
};

int map(Mapper & mapper, 
        StringSet<String<short> > & buckets, 
        String<Position<SeqFileIn>::Type> & fin_pos,
        int gid, 
        int f_buckets_enabled,
        int p1,
        bool f_io_append = false);
/*
int map(Mapper & mapper, 
        StringSet<FeaturesDynamic> & f2, 
        StringSet<String<short> > & buckets, 
        String<Position<SeqFileIn>::Type> & fin_pos,
        P_Tasks & p_tasks,
        P_Parms & p_parms,
        int gid, 
        int f_buckets_enabled,
        int p1,
        bool f_io_append = false);
*/
int filter(Mapper & mapper, StringSet<String<short> > & buckets, String<Position<SeqFileIn>::Type> & fin_pos, int p1);
int print_mapper_results(Mapper & mapper, 
    int f_p_mapper = 0, int p_in_id = 0, int p_out_id = 0); //parms to enable P_Mapper

#endif
