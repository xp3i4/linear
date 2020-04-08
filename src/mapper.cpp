#include <iostream>
#include <fstream>
#include <ctime>
#include "cords.h"
#include "pmpfinder.h"
//#include "gap.h"
#include "align_interface.h"
#include "mapper.h"
//#include "test_units.h"

using namespace seqan; 
using std::cerr;

//efficient 
MapParms parm1 ( 
        base_block_size_,     //blockSize,
        //Const_::_DELTA,          //delta(Const_::_DELTA),
        64,                          //delta
        base_threshold_,     //threshold(Const_::_THRESHOLD),
        base_kmer_step_,       //kmerStep(Const_::_KMERSTEP),
        base_shape_len_,      //shapeLen(Const_::_SHAPELEN),
        1,                      //senstivity(0),
        0,                      //anchorDeltaThr(),
        1000,                   //minReadLen(1000),
        10,                      //listN
        20,                      //listN2
        15,                     //alpha(Const_::_ALPHA),
        5,                      //alpha2 for complex mapping 
        0.02,                    //anchorLenThr(0.02),    anchors with lenghth > this parameter is pushed into the queue
        0.5,                     //rcThr(0.75)
        0.7,                     //cordThr length of cord < cordThr are abandone
        0.7,                     //senthr: perfrom next filter on cords of length < senthr 
        0.1                      //clsthr: thread of cluster
); 

//normal
MapParms parm0 ( 
        base_block_size_,     //blockSize,
        base_delta_,          //delta(Const_::_DELTA),
        base_threshold_,     //threshold(Const_::_THRESHOLD),
        base_kmer_step_,       //kmerStep(Const_::_KMERSTEP),
        base_shape_len_,      //shapeLen(Const_::_SHAPELEN),
        0,                      //senstivity(0),
        0,                      //anchorDeltaThr(),
        1000,                   //minReadLen(1000),
        10,                      //listN
        20,                      //listN2
        15,                     //alpha(Const_::_ALPHA),
        5,                      //alpha2 for complex mapping
        0.02,                    //anchorLenThr(0.02),    anchors with lenghth > this parameter is pushed into the queue
        0.5,                     //rcThr(0.8)
        0.2,                     //cordThr length of cord < cordThr are abandoned
        0.2,                     //senthr: length of cord < senthr are erased duing path
        0.1                      //clsthr: thread of cluster

); 

//sensitive
MapParms parm2 ( 
        base_block_size_,     //blockSize,
        //Const_::_DELTA,          //delta(Const_::_DELTA),
        64,                      //delta
        base_threshold_,     //threshold(Const_::_THRESHOLD),
        base_kmer_step_,       //kmerStep(Const_::_KMERSTEP),
        base_shape_len_,      //shapeLen(Const_::_SHAPELEN),
        2,                      //senstivity(0),
        0,                      //anchorDeltaThr(),
        1000,                   //minReadLen(1000),
        4,                      //listN
        50,                      //listN2
        0.65,                     //alpha(Const_::_ALPHA),
        0.5,                      //alpha2 for complex mapping
        0.02,                    //anchorLenThr(0.02),    anchors with lenghth > this parameter is pushed into the queue
        0.8,                     //rcThr(0.75)
        0.8,                      //cordThr length of cord < cordThr are abandoned
        0.8,                     //senthr: length of cord < senthr are erased duing path
        0.1                      //clsthr: thread of cluster

        
); 


MapParms parmt ( 
        base_block_size_,     //blockSize,
        base_delta_,          //delta(Const_::_DELTA),
        base_threshold_,     //threshold(Const_::_THRESHOLD),
        base_kmer_step_,       //kmerStep(Const_::_KMERSTEP),
        base_shape_len_,      //shapeLen(Const_::_SHAPELEN),
        0,                      //senstivity(0),
        0,                      //anchorDeltaThr(),
        1000,                   //minReadLen(1000),
        2,                         //listN
        2,                      //listN2
        0.75,                     //alpha(Const_::_ALPHA),
        0.65,                      //alpha2 for complex mapping
        0.02,                    //anchorLenThr(0.02),    anchors with lenghth > this parameter is pushed into the queue
        0.8,                     //rcThr(0.8)
        0.8,                       //cordThr length of cord < cordThr are abandoned
        0.8,                     //senthr: length of cord < senthr are erased duing path
        0.1                      //clsthr: thread of cluster
        
); 


/**
 * flags controlling print func;
 */ 
struct F_Print_
{
    void setPrintSam(uint & f){f |= 2;}
    int isPrintSam(uint f){return f & 2;}
}fp_handler_;
/**
 * flags controlling map func;
 */
struct F_Map_
{
    void setApxChainOFF(uint & f){f &= ~8;}
    void setApxChainON(uint & f){f |= 8;}
    void setMapGapOFF(uint & f){f &= ~2;}
    void setMapGapON(uint & f){f |= 2;}
    void setAlignOFF(uint & f){f &= ~4;}
    void setAlignON(uint & f){f |= 4;}
    uint isApxChain(uint & f){return f & 8;}
    uint isMapGap (uint & f){return f & 2;}
    uint isAlign (uint & f){return f & 4;}
}fm_handler_;

Mapper::Mapper(Options & options):
               record(options),
               gap_parms(0.1),
               index_dynamic(getGenomes())
{
    loadOptions(options);
    loadGenomes();
}

void Mapper::loadOptions(Options & options)
{
    uint64_t thd_gap_default = 50; //minium length of gap > 50 bases.
    uint64_t thd_gap_lower = 10; 
    r_paths = options.r_paths;
    g_paths = options.g_paths; 
    cord_size = window_size;
    //dout << "loadoptions" << options.bal_flag << "\n";
    switch (options.sensitivity)
    {
        case 0: 
        {
            map_parms = parm0; //normal
            break;
        }
        case 1:
        {
            map_parms =  parm1; //fast
            break;
        }
        case 2:
        {
            map_parms = parm2; //sensitive
            break;
        }
    }
    _thread = options.thread;
    //-------index---------
    index_dynamic.setIndexType (options.index_t);
    //-------feature-------
    feature_type = options.feature_t;
    of_type = OF_NEW;
    f_map = 0;
    gap_len_min = options.gap_len;
    f_print = 0;
    if (options.gap_len == 0)
    {
        fm_handler_.setMapGapOFF(f_map);
    }
    else
    {
        if (options.gap_len == 1)                  //set to default
        {
            gap_len_min = thd_gap_default;       //set default 50
        }
        else if (options.gap_len == 2)      //just another default option
        {
            gap_len_min = thd_gap_lower;
        }
        else if (options.gap_len < 10)
        {
            gap_len_min = thd_gap_lower; 
        }
        else
        {
            gap_len_min = options.gap_len;
        }
        fm_handler_.setMapGapON(f_map);
    }
    dout << "gap_len"<< gap_len_min << options.gap_len << "\n";
    if (options.apx_chain_flag == 0)
    {
        fm_handler_.setApxChainOFF(f_map);
    }
    else
    {
        fm_handler_.setApxChainON(f_map);
    }
    if (options.aln_flag == 0)
    {
        fm_handler_.setAlignOFF(f_map);
    }
    else
    {
        fm_handler_.setAlignON(f_map);
        fp_handler_.setPrintSam(f_print);
    }
    dout << "sam_flag " << options.sam_flag << f_print << options.apx_chain_flag << (f_map & 8) << "\n";
    if (options.sam_flag)
    {
        fp_handler_.setPrintSam(f_print);
    }
}
int Mapper::createIndex(unsigned gstr, unsigned gend, bool efficient)
{
    createIndexDynamic(getGenomes(), index_dynamic, gstr, gend, _thread, efficient);
    return 0;
}

int Mapper::getFeatureType()
{
    return feature_type;
}

void Mapper::setOfNew ()
{
    of_type = OF_NEW;
}

void Mapper::setOfApp ()
{
    of_type = OF_APP;
}

bool Mapper::isOfNew()
{
    return of_type == OF_NEW;
}

bool Mapper::isOfApp()
{
    return of_type == OF_APP;
}

Options::PathsType & Mapper::getRPaths()
{
    return r_paths;
}

Options::PathsType & Mapper::getGPaths()
{
    return g_paths;
}

void Mapper::loadGenomes()
{
    loadRecords(getGenomes(), getGenomesId(), getGPaths());
}
void Mapper::clearIndex()
{
    index_dynamic.clearIndex();
}

//=== pipeline2 of parallel buffer 

int Mapper::p_calRecords(int in_id, int out_id) 
{
    StringSet<String<Dna5> > & reads = reads_buffer[in_id];
    StringSet<CharString> & reads_id = reads_ids_buffer[in_id]; 
    StringSet<String<uint64_t> > & cords_str = cords1_buffer[out_id];
    StringSet<String<uint64_t> > & cords_end = cords2_buffer[out_id];
    StringSet<String<BamAlignmentRecordLink> > & bam_records = bam_link_buffer[out_id];
    clear(cords_str);
    clear(cords_end);
    clear(bam_records);
    resize(cords_str, length(reads));
    resize(cords_end, length(reads));
    resize(bam_records, length(reads));

    Anchors anchors;
    String<uint64_t> crhit;
    String<Dna5> comStr; //complement revers of read seq
    String<UPair> apx_gaps; 
    StringSet<String<uint64_t> > clips;
    StringSet<FeaturesDynamic> f1;
    StringSet<FeaturesDynamic> & f2 = this->getGenomesFeatures();
    unsigned feature_window_size = getFeatureWindowSize(f2);
    float thd_err_rate = 0.2;    //todo::wrapper parms here
    uint thd_min_read_len = 200; //todo::wrapper parms here

    resize(f1, 2);
    f1[0].init(f2[0].fs_type);
    f1[1].init(f2[0].fs_type);
    resize (clips, length(cords_str));
    int f_chain = fm_handler_.isApxChain(f_map) ? 1 : 0; 
    for (unsigned j = 0; j < length(reads); j++)
    {
        if (length(reads[j]) > thd_min_read_len)
        {
            _compltRvseStr(reads[j], comStr);
            createFeatures(begin(reads[j]), end(reads[j]), f1[0]);
            createFeatures(begin(comStr), end(comStr), f1[1]);
            apxMap(this->getIndex(), reads[j], anchors, this->getMapParms(), crhit, f1, f2, apx_gaps, cords_str[j], f_chain);
            if (fm_handler_.isMapGap(f_map))
            {
                //<<debug
                gap_parms.read_id = reads_id[j];
                //>>debug
                mapGaps(this->getGenomes(), reads[j], comStr, cords_str[j], cords_end[j], clips[j], apx_gaps, f1, f2, gap_len_min, feature_window_size, thd_err_rate, this->getGapParms());
            }
            if (fm_handler_.isAlign(f_map))
            {
                align_cords(this->getGenomes(), reads[j], comStr, cords_str[j], bam_records[j]);
                //check_cigar (seqs, reads[j], comStr, cordsTmp[c], bam_records_tmp[c]);
            }
        }
    } 
    return 0;
}

int Mapper::p_printResults(int it1, int it2)
{
    std::string path = this->getPReadsPathsBuffer()[it1];
    std::string prefix = getFileName(path, "/", ~0);
    this->outputPrefix = getFileName(prefix, ".", 0);
    //prefix = path;
    std::cout << "reads_path"  << this->outputPrefix << "\n"; // (this->outputPrefix) << "\n";
    print_mapper_results(*this, 1, it1, it2);
    return 0;
}

/*----------  Wrapper of file I/O   ----------*/
int print_cords_apf(Mapper & mapper, int f_p_mapper, int p_in_id, int p_out_id)
{
    if (f_p_mapper) //P buffer enabled
    {
        print_cords_apf(mapper.getPCords1Buffer()[p_out_id],
                        mapper.getGenomes(),
                        mapper.getPReadsBuffer()[p_in_id],
                        mapper.getGenomesId(),
                        mapper.getPReadsIdBuffer()[p_in_id],
                        mapper.getOf());
    }
    else
    {
        print_cords_apf(mapper.getCords(), 
                        mapper.getGenomes(),
                        mapper.getReads(),
                        mapper.getGenomesId(),
                        mapper.getReadsId(),
                        mapper.getOf());
    }
    return 0;
}
int print_align_sam (Mapper & mapper, int f_p_mapper, int p_in_id, int p_out_id)
{
    if (f_p_mapper)
    {
        print_align_sam(mapper.getGenomes(),
                        mapper.getGenomesId(),
                        mapper.getPReadsIdBuffer()[p_in_id],
                        mapper.getPBamLinksBuffer()[p_out_id],
                        mapper.getOf(),
                        mapper.isOfNew());
    }
    else
    {
        print_align_sam (mapper.getGenomes(),
                         mapper.getGenomesId(),
                         mapper.getReadsId(),
                         mapper.getBamRecords(),
                         mapper.getOf(),
                         mapper.isOfNew());  
    }
    return 0;
}
int print_clips_gvf(Mapper & mapper)
{
    print_clips_gvf_(mapper.getClips(), 
                     mapper.getReadsId(), 
                     mapper.getGenomesId(), 
                     mapper.getOf());
    return 0;
}
int print_cords_sam(Mapper & mapper, int f_p_mapper, int p_in_id, int p_out_id)
{
//    uint64_t thd_large_X = 80; //cigar containing X > this will be clipped into 2 records
    uint64_t thd_large_X = 8000; //cigar containing X > this will be clipped into 2 records
    //todo:: fix this parm, too large
    int f_parallel = 1;
    if (f_p_mapper)
    {
        clear(mapper.getPBamLinksBuffer()[p_out_id]);
        f_parallel = 0;
        dout << "pcss1" << length(mapper.getPCords1Buffer()[p_out_id])
                        << length(mapper.getPCords2Buffer()[p_out_id])
                        << length(mapper.getPBamLinksBuffer()[p_out_id])
                        << length(mapper.getPReadsIdBuffer()[p_in_id])
                        << length(mapper.getPReadsIdBuffer()[p_in_id])
                        << "\n";
        print_cords_sam(mapper.getPCords1Buffer()[p_out_id],
                        mapper.getPCords2Buffer()[p_out_id],
                        mapper.getPBamLinksBuffer()[p_out_id],
                        mapper.getGenomesId(),
                        mapper.getPReadsIdBuffer()[p_in_id],
                        mapper.getGenomes(),
                        mapper.getPReadsBuffer()[p_in_id],
                        mapper.getCordSize(),
                        mapper.getOf(),
                        thd_large_X,
                        mapper.getThreads(),
                        mapper.isOfNew(),
                        f_parallel);
    }
    else
    {
        dout << "pcssd1" << length(mapper.getCords()) << length(mapper.getBamRecords()) << "\n";
        print_cords_sam(mapper.getCords(),
                        mapper.getCords2(),
                        mapper.getBamRecords(),
                        mapper.getGenomesId(),
                        mapper.getReadsId(),
                        mapper.getGenomes(),
                        mapper.getReads(),
                        mapper.getCordSize(),
                        mapper.getOf(),
                        thd_large_X,
                        mapper.getThreads(),
                        mapper.isOfNew(),
                        f_parallel);
    }
    return 0;
}
/**
 * Open new or append to original 
 */
void open_mapper_of(Mapper & mapper, std::string file_path)
{
    if (mapper.isOfNew())
    {
        //std::cerr << "new............\n" << file_path << "\n\n";
        mapper.getOf().open(toCString(file_path));
    }
    else if (mapper.isOfApp())
    {
        //std::cerr << "n..............\n" << file_path << "\n\n";
        mapper.getOf().open(toCString(file_path), std::ios::app);
    }
}
void close_mapper_of (Mapper & mapper)
{
    close (mapper.getOf());
}
/** 
 * Print main apf sam and gvf
 */
int print_mapper_results(Mapper & mapper, 
    int f_p_mapper, int p_in_id, int p_out_id) //parms to enable P_Mapper
{
    ///.apf
    std::string file1 = mapper.getOutputPrefix() + ".apf";
    std::cout << "file ====== "  << "\n";
    open_mapper_of (mapper, file1);
    print_cords_apf(mapper, f_p_mapper, p_in_id, p_out_id);
    close_mapper_of(mapper);
    ///.gvf
    /*
    std::string file2 = mapper.getOutputPrefix() + ".gvf";
    open_mapper_of (mapper, file2);
    print_clips_gvf(mapper);
    close_mapper_of(mapper);
    */

    ///.sam
    std::string file3 = mapper.getOutputPrefix() + ".sam";
    open_mapper_of (mapper, file3);
    if (fp_handler_.isPrintSam(mapper.getPrintFlag()))
    {
        if (fm_handler_.isAlign(mapper.getMapFlag()))
        {
            print_align_sam(mapper, f_p_mapper, p_in_id, p_out_id);
        }
        else
        {
            print_cords_sam(mapper, f_p_mapper, p_in_id, p_out_id);
        }
    }
    close_mapper_of(mapper);

    mapper.setOfApp(); //set of_type to std::ios::app;
    return 0;
}
//read records from fin_pos and buckckets
int readRecords4FinPosbuckets(StringSet<CharString> & ids, StringSet<String<Dna5> > & reads, 
                      StringSet<String<short> >& buckets, 
                      String<Position<SeqFileIn>::Type> & fin_pos,
                      SeqFileIn & fin, uint rstr, uint rend, uint bucketId)
{
    CharString tmp_id;
    String<Dna5> tmp_read;
    for (int i = rstr; i < rend & !atEnd(fin); i++)
    {
        readRecord  (tmp_id, tmp_read, fin);
        if (buckets[i][bucketId]) //ith read and bucketid genome 
        {
            //setPosition (fin, fin_pos[i]);
            appendValue (ids, tmp_id);
            appendValue (reads, tmp_read);
        }
    }     
    return 0;
}

int readRecords2FinPosBuckets(StringSet<CharString> & ids, StringSet<String<Dna5> > & reads, 
                      StringSet<String<short> >& buckets, 
                      String<Position<SeqFileIn>::Type> & fin_pos, SeqFileIn & fin, 
                      uint rstr, uint rend)
{
    CharString tmp_id;
    String<Dna5> tmp_read;
    if (empty(fin_pos)) 
    {
        appendValue(fin_pos, Position<SeqFileIn>::Type(0));
    }
    for (int i = rstr; i < rend && !atEnd(fin); i++)
    {
        readRecord (tmp_id, tmp_read, fin);
        appendValue(ids, tmp_id);
        appendValue(reads, tmp_read);
        appendValue(fin_pos, seqan::position(fin));
    }
    return 0;
}

/*----------  Map main funcion  ----------*/

int map_(IndexDynamic & index,
         StringSet<FeaturesDynamic > & f2,
         StringSet<String<Dna5> > & reads,
         StringSet<CharString> & readsId,
         MapParms & mapParm,
         StringSet<String<uint64_t> > & cords_str,
         StringSet<String<uint64_t> > & cords_end,
         StringSet<String<uint64_t> > & clips,
         StringSet<String<Dna5> > & seqs,
         StringSet<String<BamAlignmentRecordLink> >& bam_records,
         uint f_map,   //control flags
         uint gap_len_min,
         uint threads,
         int p1)
{
    unsigned feature_window_size = getFeatureWindowSize(f2);
    dout << "fe" << feature_window_size << "\n";
    float senThr = mapParm.senThr / feature_window_size;  //map for 2 roun if cords cover len <
    float cordThr = 0.3 / feature_window_size; //cords cover length < are aborted
    uint thd_min_read_len = 200;
    //todo::tune the cordThr try to merge cords of blocks 
    MapParms complexParm = mapParm;
    complexParm.alpha = complexParm.alpha2;
    complexParm.listN = complexParm.listN2;
    String<uint64_t> gap_len;
    String<uint64_t> red_len;
    resize (gap_len, threads, 0);
    resize (red_len, threads, 0);
    float thd_err_rate = 0.2;
    int f_chain = 1; 
    if (fm_handler_.isApxChain(f_map))
    {
        f_chain = 1;
    }
    else
    {
        f_chain = 0;
    }
#pragma omp parallel
{
    //GapParms gap_parms(0.85);
    GapParms gap_parms(0.1);
    unsigned size2 = length(reads) / threads;
    unsigned ChunkSize = size2;
    String<Dna5> comStr;
    Anchors anchors;
    String<uint64_t> crhit;
    StringSet<String<uint64_t> > cordsTmp;  //cords_str
    StringSet<String<uint64_t> > cordsTmp2; //cords_end
    StringSet<FeaturesDynamic> f1;
    StringSet<String<uint64_t> > clipsTmp;
    StringSet<String<BamAlignmentRecordLink> > bam_records_tmp;
    unsigned thd_id = omp_get_thread_num();
    if (thd_id < length(reads) - size2 * threads)
    {
        ChunkSize = size2 + 1;
    }
    resize(cordsTmp, ChunkSize);
    resize(cordsTmp2, ChunkSize);
    resize(clipsTmp, ChunkSize);
    resize(bam_records_tmp, ChunkSize);
    resize(f1, 2);
    f1[0].init(f2[0].fs_type);
    f1[1].init(f2[0].fs_type);
    unsigned c = 0;
    
    String<UPair> apx_gaps; 
    #pragma omp for
    for (unsigned j = 0; j < length(reads); j++)
    {
        std::cout << "readid " << j << readsId[j] << "\n\n";
        double t1 = sysTime ();
        red_len[thd_id] += length(reads[j]);
        
        if (length(reads[j]) > thd_min_read_len)
        {
            _compltRvseStr(reads[j], comStr);
            createFeatures(begin(reads[j]), end(reads[j]), f1[0]);
            createFeatures(begin(comStr), end(comStr), f1[1]);
            apxMap(index, reads[j], anchors, mapParm, crhit, f1, f2, apx_gaps, cordsTmp[c], f_chain);
            if (fm_handler_.isMapGap(f_map))
                {
                //<<debug
                gap_parms.read_id = readsId[j];
                //>>debug
                mapGaps(seqs, reads[j], comStr, cordsTmp[c], cordsTmp2[c], clipsTmp[c], apx_gaps, f1, f2, gap_len_min, feature_window_size, thd_err_rate, gap_parms);
                }
            if (fm_handler_.isAlign(f_map))
            {
                align_cords(seqs, reads[j], comStr, cordsTmp[c], bam_records_tmp[c]);
                //check_cigar (seqs, reads[j], comStr, cordsTmp[c], bam_records_tmp[c]);
            }
        }
        
        c += 1;
    } 
    #pragma omp for ordered
    for (unsigned j = 0; j < threads; j++)
    #pragma omp ordered
    {
        append(cords_str, cordsTmp);
        append(cords_end, cordsTmp2);
        if (fm_handler_.isMapGap(f_map))
        {
            append(clips, clipsTmp);
        }
        if (fm_handler_.isAlign(f_map))
        {
            append(bam_records, bam_records_tmp);
        }
    }
}
    return 0;
}

/**
 * Map main 
 * Stream all reads records in reads files specified by @path 
 * !!CRITICAL::this function is required to stream the read records of each read file (@path) 
   in the same way the function filter() did. Otherwise the @buckets can't work properly.
 */
int map(Mapper & mapper, 
        StringSet<String<short> > & buckets, 
        String<Position<SeqFileIn>::Type> & fin_pos,
        int gid, 
        int f_buckets_enabled,
        int p1,
        bool f_io_append)
{
    //std::cout << "mapf " << mapper.getFeatureType() << "\n";
    unsigned blockSize = 50000;
    uint rstr = 0;
    SeqFileIn rFile;
    StringSet<std::string> file1s;
    StringSet<std::string> file2s;
    StringSet<std::string> file3s;
    for (auto path : mapper.getRPaths())
    {
        if(!open(rFile, toCString(path)))
        {
            serr.print_message("\033[1;31mError:\033[0m can't open read file ", 2, 0, std::cerr);
            serr.print_message(toCString(path), 0, 1, std::cerr);
            continue; 
        }
        std::string outputPrefix = getFileName(path, "/", ~0);
        outputPrefix = getFileName(outputPrefix, ".", 0);
        mapper.getOutputPrefix() = outputPrefix;
        if (f_io_append){
            mapper.setOfApp();
        }
        else {mapper.setOfNew();}
        unsigned k = 1;
        while (!atEnd(rFile))
        {
            double time1 = sysTime();
            serr.print_message("=>Map::file_I/O", 0, 0, std::cerr);
            serr.print_message(k, 0, 2, std::cerr);
            try
            {
                if (f_buckets_enabled)
                {
                    readRecords4FinPosbuckets(mapper.getReadsId(), mapper.getReads(), buckets, fin_pos, rFile, 
                        rstr, rstr + blockSize, gid);
                    rstr += blockSize;
                }
                else
                {
                    readRecords(mapper.getReadsId(), mapper.getReads(), rFile, blockSize);
                }
            }
            catch (Exception const & e)
            {

            }
            //serr.print_message("", 50, 2, std::cerr); 
            serr.print_message("=>Map::mapping ", 0, 0, std::cerr);
            serr.print_message(path, 0, 0, std::cerr);
            serr.print_message (" block ", 0, 0, std::cerr);
            serr.print_message(k, 0, 0, std::cerr);
            serr.print_message(" Size ", 0, 0, std::cerr);
            serr.print_message(unsigned(length(mapper.getReads())), 0, 2, std::cerr);
            time1 = sysTime() - time1;
            double time2 = sysTime();
            map_(mapper.getIndex(), 
                 mapper.getGenomesFeatures(),
                 mapper.getReads(), 
                 mapper.getReadsId(),
                 mapper.getMapParms(), 
                 mapper.getCords(),  //cords_str 
                 mapper.getCords2(), //cords_end
                 mapper.getClips(),
                 mapper.getGenomes(),
                 mapper.getBamRecords(),
                 mapper.getMapFlag(),
                 mapper.getGapLenMin(),
                 mapper.getThreads(), 
                 p1);
            time2 = sysTime() - time2;
            double time3 = sysTime();
            serr.print_message("=>Write results to disk", 0, 2, std::cerr);

            print_mapper_results(mapper);

            clear (mapper.getCords());
            clear (mapper.getCords2());
            clear (mapper.getClips());
            clear (mapper.getBamRecords());
            clear (mapper.getReads());
            clear (mapper.getReadsId());
            clear (mapper.getClips());
            clear (mapper.getBamRecords());
            time3 = sysTime() - time3;
            std::cerr <<  "--Map::file " << path << " block "<< k << " Size " << length(mapper.getReads()) << " Elapsed Time[s]: file_I/O " << time1 + time3 << " map "<< time2 << "\n";
            k++;
        }      
        std::string file1 = mapper.getOutputPrefix() + ".apf";
        std::string file2 = mapper.getOutputPrefix() + ".gvf";
        std::string file3 = mapper.getOutputPrefix() + ".sam";
        appendValue (file1s, file1);
        appendValue (file2s, file2);
        appendValue (file3s, file3);
        close(rFile);
    }
    serr.print_message("--Write results to disk 100%", 0, 1, cerr);
    for (uint i = 0; i < length(file1s); i++)
    {
        serr.print_message("Result files: \033[1;31m" + file1s[i] + "\033[0m ", 2, 0, cerr);
        serr.print_message("\033[1;31m" + file2s[i] + "\033[0m ", 0, 0, cerr);
        serr.print_message("\033[1;31m" + file3s[i] + "\033[0m ", 16, 1, cerr); 
    }
    //!!TODO::clear index;
    //mapper.index().clear(); 
    //std::cerr << ">Write results to disk        \r";
    //print_mapper_results(mapper);
    return 0;
}
/*
 * parallel io with dynamic balancing task scheduling
 *
int map(Mapper & mapper, 
        StringSet<String<short> > & buckets, 
        String<Position<SeqFileIn>::Type> & fin_pos,
        P_Tasks & p_tasks,
        P_Parms & p_parms,
        int gid, 
        int f_buckets_enabled,
        int p1,
        bool f_io_append)
{
    dout << "p_taskend" << p_tasks.isTaskEnd(0) << p_tasks.tasks[0].task_type << "\n";
#pragma omp parallel
{
    std::cout << "parallel\n";
    p_threadprocess(p_tasks, p_parms, omp_get_thread_num());
}
    return 0;
}
*/
//Shortcut called within filter function to marked the genome requiring map
//@bucket[i][j] == 1:the the ith read should map to the jth genome; otherwise not
void append_genome_bucket(StringSet<String<short> > & buckets, 
                          StringSet<String<uint64_t> > & cords, 
                          unsigned gnmu)
{
    String<short> new_bucket;
    resize (new_bucket, gnmu, 0);
    for (uint i = 0; i < length(cords); i++)
    {
        appendValue(buckets, new_bucket);
        for (uint j = 1; j < length(cords[i]); j++)
        {
            back(buckets)[get_cord_id(cords[i][j])] = 1;
        }
    }
}

int filter_(IndexDynamic & index,
            StringSet<FeaturesDynamic > & f2,
            StringSet<String<Dna5> > & reads,
            MapParms & mapParm,
            StringSet<String<uint64_t> > & cords_str,
            StringSet<String<uint64_t> > & cords_end,
            StringSet<String<uint64_t> > & clips,
            StringSet<String<Dna5> > & seqs,
            StringSet<String<BamAlignmentRecordLink> > & bam_records,
            uint f_map,   //control flags
            uint gap_len_min,
            uint threads,
            int p1)
{
    float senThr = mapParm.senThr / window_size;  //map for 2 roun if cords cover len <
    float cordThr = 0.3 / window_size; //cords cover length < are aborted
    uint thd_min_read_len = 200;
    //todo::tune the cordThr try to merge cords of blocks 
    MapParms complexParm = mapParm;
    String<uint64_t> gap_len;
    String<uint64_t> red_len;
    resize (gap_len, threads, 0);
    resize (red_len, threads, 0);
    float thd_err_rate = 0.2;
    int f_chain = 1; 
#pragma omp parallel
{
    unsigned size2 = length(reads) / threads;
    unsigned ChunkSize = size2;
    String<Dna5> comStr;
    Anchors anchors;
    String<uint64_t> crhit;
    StringSet<String<uint64_t> > cordsTmp;  //cords_str
    StringSet<String<uint64_t> > cordsTmp2; //cords_end
    StringSet<FeaturesDynamic> f1;
    unsigned thd_id = omp_get_thread_num();
    if (thd_id < length(reads) - size2 * threads)
    {
        ChunkSize = size2 + 1;
    }
    resize(cordsTmp, ChunkSize);
    resize(cordsTmp2, ChunkSize);
    resize(f1, 2);
    f1[0].init(f2[0].fs_type);
    f1[1].init(f2[0].fs_type);
    unsigned c = 0;

    String<UPair> apx_gaps; 
    #pragma omp for
    for (unsigned j = 0; j < length(reads); j++)
    {
        double t1 = sysTime ();
        red_len[thd_id] += length(reads[j]);
        if (length(reads[j]) > thd_min_read_len)
        {
            _compltRvseStr(reads[j], comStr);
            createFeatures(begin(reads[j]), end(reads[j]), f1[0]);
            createFeatures(begin(comStr), end(comStr), f1[1]);
            apxMap(index, reads[j], anchors, mapParm, crhit, f1, f2, apx_gaps, cordsTmp[c], f_chain);
            //filterGenomes(index, reads[j], anchors, mapParm, crhit, f1, f2, apx_gaps, cordsTmp[c], cordLenThr, f_chain);
        }
        c += 1;
    } 
    #pragma omp for ordered
    for (unsigned j = 0; j < threads; j++)
    #pragma omp ordered
    {
        append(cords_str, cordsTmp);
        append(cords_end, cordsTmp2);
    }
}
    return 0;
}

/*
 * Filter control interface
 */
int filter(Mapper & mapper, 
          StringSet<String<short> > & buckets, 
          String<Position<SeqFileIn>::Type> & fin_pos, int p1)
{
    unsigned blockSize = 50000;
    SeqFileIn rFile;
    uint rstr = 0;
    //for (auto path : mapper.getRPaths())
    for (int i = 0; i < length(mapper.getRPaths()); i++)
    {
        auto & path = mapper.getRPaths()[i];
        if(!open(rFile, toCString(path)))
        {
            serr.print_message("\033[1;31mError:\033[0m can't open read file ", 2, 0, std::cerr);
            serr.print_message(toCString(path), 0, 1, std::cerr);
            continue; 
        }
        std::string outputPrefix = getFileName(path, "/", ~0);
        outputPrefix = getFileName(outputPrefix, ".", 0);
        mapper.getOutputPrefix() = outputPrefix;
        mapper.setOfNew();
        unsigned k = 1;
        while (!atEnd(rFile))
        {
            double time1 = sysTime();
            serr.print_message("=>Filter::I/O::reading records", 0, 0, std::cerr);
            serr.print_message(k, 0, 2, std::cerr);
            try
            {
                //readRecords_buckets(mapper.getReadsId(), mapper.getReads(), rFile, blockSize);
                readRecords2FinPosBuckets(mapper.getReadsId(), mapper.getReads(), buckets, fin_pos, rFile, 
                    rstr, rstr + blockSize);

            }
            catch (Exception const & e)
            {
                std::cerr << "\033[1;31mError:\033[0m can't read records " << e.what () << "\n";
                //none;
            }
            //serr.print_message("", 50, 2, std::cerr); 
            serr.print_message ("=>Filiter::filtering genomes ", 0, 0, std::cerr);
            serr.print_message (path, 0, 0, std::cerr);
            serr.print_message (" block ", 0, 0, std::cerr);
            serr.print_message (k, 0, 0, std::cerr);
            serr.print_message (" Size ", 0, 0, std::cerr);
            serr.print_message (unsigned(length(mapper.getReads())), 0, 2, std::cerr);
            time1 = sysTime() - time1;
            double time2 = sysTime();
            filter_(mapper.getIndex(), 
                mapper.getGenomesFeatures(),
                 mapper.getReads(), 
                 mapper.getMapParms(), 
                 mapper.getCords(),  //cords_str 
                 mapper.getCords2(), //cords_end
                 mapper.getClips(),
                 mapper.getGenomes(),
                 mapper.getBamRecords(),
                 mapper.getMapFlag(),
                 mapper.getGapLenMin(),
                 mapper.getThreads(), 
                 p1);
            time2 = sysTime() - time2;
            double time3 = sysTime();
            serr.print_message("=>Recording geonme buckets", 0, 2, std::cerr);
            //std::cerr << std::flush;
            append_genome_bucket(buckets, mapper.getCords(), length(mapper.getGenomes()));
            std::cerr <<  "--Filter::genomes " << path << " block "<< k << " Size " << length(mapper.getReads()) << " Elapsed Time[s]: file_I/O " << time1 << " filter "<< time2 << " append buckets " << sysTime() - time3  << "\n";
            clear (mapper.getCords());
            clear (mapper.getCords2());
            clear (mapper.getClips());
            clear (mapper.getBamRecords());
            clear (mapper.getReads());
            clear (mapper.getReadsId());
            clear (mapper.getClips());
            clear (mapper.getBamRecords());
            k++;
        }      
        close(rFile);
    }
    return 0;
}

