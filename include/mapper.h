#ifndef SEQAN_HEADER_PACMAPPER_H
#define SEQAN_HEADER_PACMAPPER_H

#include "base.h"
#include "index_util.h"
#include "f_io.h"
#include "gap.h"
#include "parallel_io.h"
//#include "mapparm.h"
struct MapParms{
    unsigned    blockSize;
    unsigned    delta;
    unsigned    threshold;
    unsigned    kmerStep;
    unsigned    shapeLen;
    unsigned    sensitivity;
    unsigned    anchorDeltaThr;
    unsigned    minReadLen;
    unsigned    listN;
    unsigned    listN2;
    float       alpha;
    float       alpha2;
    float       anchorLenThr;
    float       rcThr;
    float       cordThr;
    float       senThr;
    float       clsThr;

    GlobalParms pm_g;
    PMPParms pm_pmp; 
      
    MapParms();
    MapParms(unsigned bs, unsigned dt, unsigned thr, 
            unsigned ks, unsigned sl, unsigned st,
            unsigned ad, unsigned mr, unsigned listn,
            unsigned listn2,
            float ap, float ap2, float alt, float rt, float ct, float sent, float clst);
    //MapParms(MapParms & parm);
    void setMapParm(Options & options);
    void print ();
};

class Mapper : public P_Mapper 
{
    PMRecord    record;
    MapParms     map_parms;
    String<MapParms> map_parms_set;
    String<GapParms> gap_parms_set; //each thread keep a copy of gap parms
    IndexDynamic index_dynamic;
    StringSet<String<uint64_t> > cordSet;
    StringSet<String<uint64_t> > cordSet2;
    StringSet<String<CordInfo> > cords_info;
    std::ofstream of;
    unsigned _thread;
    String<size_t> rlens;
    StringSet<String<uint64_t> > clip_set;
    StringSet<String<BamAlignmentRecordLink> > bam_records;
    StringSet<FeaturesDynamic > f2; //features of genomes
    std::string outputPrefix;
    int feature_type;
    int of_type;
    int f_new_file;
    int f_output_set;
    const int OF_APP = 0; //set of std::ios::app
    const int OF_NEW = 1;
    Options::PathsType r_paths;
    Options::PathsType g_paths;
    uint f_map;     //map control flag 
    uint f_print;   //print control flag
    int cord_size; //default cord size 
    FIOParms fio_parms;

    //=== pipeline2 of parallel buffer 
    //P_ReadsBuffer reads_buffer;
    //P_ReadsIdsBuffer reads_ids_buffer; 
    //P_CordsBuffer cords_buffer;
    //P_BamLinkBuffer bam_link_buffer;

public:
    Mapper();
    Mapper(Options & options);
    int loadOptions(Options & options);
    int setMapperBamHeaders(Options & options);
    StringSet<String<Dna5> > & getReads() {return record.seq1;}             
    StringSet<String<Dna5> > & getGenomes() {return record.seq2;}             
    MapParms & getMapParms() {return map_parms;}
    IndexDynamic & getIndex() {return index_dynamic;}
    StringSet<String<uint64_t> > & getCords() {return cordSet;}             
    StringSet<String<uint64_t> > & getCords2() {return cordSet2;}            
    StringSet<String<CordInfo> > & getCordsInfo() {return cords_info;}
    
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
    int  getCordSize() {return cord_size;}
    FIOParms & getFIOParms(){return fio_parms;}
    void loadGenomes();
    void clearIndex();
    String<GapParms> & getGapParms(){return gap_parms_set;}

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
