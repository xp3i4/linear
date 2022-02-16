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
               index_dynamic(getGenomes())
{
    loadOptions(options);
    loadGenomes();
    setMapperBamHeaders(options);
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
    f_print = 0;
    //-------gap_parms-----
    GapParms gap_parms_template(0.2); //initial error rate 
    if (options.gap_len == 0)
    {
        fm_handler_.setMapGapOFF(f_map);
    }
    else
    {
        if (options.gap_len == 1)                  //set to default
        {
            gap_parms_template.thd_gap_len_min = thd_gap_default; //set default 50
        }
        else if (options.gap_len == 2)      //just another default option
        {
            gap_parms_template.thd_gap_len_min = thd_gap_lower;
        }
        else if (options.gap_len < 10)
        {
            gap_parms_template.thd_gap_len_min = thd_gap_lower; 
        }
        else
        {
            gap_parms_template.thd_gap_len_min = options.gap_len;
        }
        fm_handler_.setMapGapON(f_map);
    }
    for (unsigned i = 0; i < _thread; i++)
    { 
        appendValue(gap_parms_set, gap_parms_template); //a copy for each thread
    }
    //gap_parms_template.printParms("gap_parms");
    //dout << "gap_len"<< gap_len_min << options.gap_len << "\n";
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
    //dout << "sam_flag " << options.sam_flag << f_print << options.apx_chain_flag << (f_map & 8) << "\n";
    //-------f_io parms-----
    /*
    if (options.sam_flag)
    {
        fp_handler_.setPrintSam(f_print);
    }
    else 
    {
        fp_handler_.unsetPrintSam(f_print);
    }
    if (options.apf_flag)
    {
        fp_handler_.setPrintApf(f_print);
    }
    else
    {
        fp_handler_.unsetPrintApf(f_print);
    }
    */
    f_print = options.f_output_type;
    fio_parms.f_reform_ccs = options.reform_ccs_cigar_flag;
    fio_parms.read_group = options.read_group;
    fio_parms.sample_name = options.sample_name;
    fio_parms.f_print_seq = options.sequence_sam_flag;
    fio_parms.f_is_align = options.aln_flag;
    fio_parms.f_output_type = options.f_output_type;
    //fio_parms.bam_header()
    //dout << "par" << options.sequence_sam_flag << fio_parms.f_sequence_sam << "\n";
}
int Mapper::setMapperBamHeaders(Options & options)
{
    if (options.isOutputSam() || options.isOutputBamPbsv() || options.isOutputBamStandard())
    {
        BamHeaderRecord tmp_bam_header;
        for (unsigned i = 0; i < length(getGenomesId()); i++)
        {
            clear(tmp_bam_header);
            tmp_bam_header.type = seqan::BAM_HEADER_REFERENCE;
            setTagValue("SN", getGenomesId()[i], tmp_bam_header);
            setTagValue("LN", std::to_string(length(getGenomes()[i])), tmp_bam_header);
            appendValue(fio_parms.bam_header, tmp_bam_header);
            appendValue(fio_parms.bam_header2, tmp_bam_header);
        }
        clear(tmp_bam_header);
        tmp_bam_header.type = seqan::BAM_HEADER_READ_GROUP;
        setTagValue("ID", options.read_group, tmp_bam_header);
        setTagValue("SM", options.sample_name, tmp_bam_header);
        appendValue(fio_parms.bam_header, tmp_bam_header);

        clear(tmp_bam_header);
        tmp_bam_header.type = seqan::BAM_HEADER_READ_GROUP;
        setTagValue(" ID", options.read_group, tmp_bam_header);
        setTagValue("SM", options.sample_name, tmp_bam_header);
        appendValue(fio_parms.bam_header2, tmp_bam_header);

        clear(tmp_bam_header);
        tmp_bam_header.type = seqan::BAM_HEADER_PROGRAM;
        setTagValue("ID", "M1-3", tmp_bam_header);
        setTagValue("PN", "Linear", tmp_bam_header);
        setTagValue("CL", options.cmd_line, tmp_bam_header);
        appendValue(fio_parms.bam_header, tmp_bam_header);
        appendValue(fio_parms.bam_header2, tmp_bam_header);
    }

    return 0;
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
void Mapper::initBuffers(int buffer_size1, int buffer_size2, P_Parms & p_parms)
{
    P_Mapper::initBuffers(buffer_size1, buffer_size2, getGPaths(), 
        getRPaths(), getThreads(), p_parms);
}
int Mapper::p_calRecords(int in_id, int out_id, int thread_id) 
{
    Counters & counters = this->getPTask(thread_id).counters;
    StringSet<String<Dna5> > & reads = reads_buffer[in_id];
    StringSet<CharString> & reads_id = reads_ids_buffer[in_id]; 
    StringSet<String<uint64_t> > & cords_str = cords1_buffer[out_id];
    StringSet<String<uint64_t> > & cords_end = cords2_buffer[out_id];
    StringSet<String<CordInfo> > & cords_info = cords_info_buffer[out_id];
    StringSet<String<BamAlignmentRecordLink> > & bam_records = bam_link_buffer[out_id];
    clear(cords_str);
    clear(cords_end);
    clear(cords_info);
    clear(bam_records);
    resize(cords_str, length(reads));
    resize(cords_end, length(reads));
    resize(cords_info, length(reads));
    resize(bam_records, length(reads));

    Anchors anchors;
    String<uint64_t> crhit;
    String<Dna5> comStr; //complement revers of read seq
    String<UPair> apx_gaps; 
    StringSet<String<uint64_t> > clips;
    StringSet<FeaturesDynamic> f1;
    StringSet<FeaturesDynamic> & f2 = this->getGenomesFeatures();
    uint thd_min_read_len = 200; //todo::wrapper parms here
    resize(f1, 2);
    f1[0].init(f2[0].fs_type);
    f1[1].init(f2[0].fs_type);
    resize (clips, length(cords_str));
    int f_chain = fm_handler_.isApxChain(f_map) ? 1 : 0; 
    CordsParms cords_parms;
    dout << "fdone3" << "\n";
    for (unsigned j = 0; j < length(reads); j++)
    {
        if (length(reads[j]) > thd_min_read_len)
        {
            _compltRvseStr(reads[j], comStr);
            //clear(f1[0].fs2_48);
            //clear(f1[1].fs2_48);
            createFeatures(begin(reads[j]), end(reads[j]), f1[0]);
            createFeatures(begin(comStr), end(comStr), f1[1]);
            apxMap(getIndex(), reads[j], anchors, crhit, f1, f2, apx_gaps, cords_str[j], cords_end[j], cords_info[j], f_chain, getMapParms().pm_g, getMapParms().pm_pmp);
            if (fm_handler_.isMapGap(f_map))
            {
                gap_parms_set[thread_id].read_id = reads_id[j];
                mapGaps(this->getGenomes(), reads[j], comStr, cords_str[j], cords_end[j], clips[j], apx_gaps, f1, f2, gap_parms_set[thread_id]);
                print_cords(cords_str[j], "pcal1");
                reformCords(cords_str[j], cords_end[j], &reformCordsDxDy1, cords_parms);
                print_cords(cords_str[j], "pcal2");
            }
            if (fm_handler_.isAlign(f_map))
            {
                double t = sysTime();
                alignCords(this->getGenomes(), reads[j], comStr, cords_str[j], cords_end[j], bam_records[j]);
                t = sysTime() - t;
                //check_cigar (seqs, reads[j], comStr, cordsTmp[c], bam_records_tmp[c]);
            }
        }
    } 
    //dout << "fdone2" << "\n";
    if (!fm_handler_.isAlign(f_map)) 
    {
        uint64_t thd_large_X = 8000; 
        clear(bam_records);
        //print_cords(cords_str[0], "pca1");
        cords2BamLink (cords_str, cords_end, cords_info, bam_records, reads, this->getCordSize(), thd_large_X);
    }
    counters.setCalCounter(counters.getCalCounter() + length(reads));
    //dout << "fdone1" << "\n";

    return 0;
}
//it1->reads buffer
//it2->cords buffer
int Mapper::p_printResults(int it1, int it2, int thread_id)
{
    Counters & counters = this->getPTask(thread_id).counters;
    std::string path = this->getPReadsPathsBuffer()[it1];
    std::string prefix = getFileName(path, "/", ~0);
    std::string out_prefix  = getFileName(prefix, ".", 0);
    if (this->outputPrefix != out_prefix)
    {
        this->setOfNew();
        this->outputPrefix = out_prefix;
    }
    else
    {
        this->setOfApp();
    }
    print_mapper_results(*this, 1, it1, it2);
    counters.setOutCounter(counters.getOutCounter() + length(this->getPReadsBuffer()[it1]));

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
int print_clips_gvf(Mapper & mapper)
{
    print_clips_gvf_(mapper.getClips(),
                     mapper.getReadsId(),
                     mapper.getGenomesId(),
                     mapper.getOf());
    return 0;
}
int print_align_sam_bam (Mapper & mapper, int f_p_mapper, int p_in_id, int p_out_id)
{
    if (f_p_mapper)
    {
        printAlignSamBam(mapper.getGenomes(),
                         mapper.getPReadsBuffer()[p_in_id],
                         mapper.getGenomesId(),
                         mapper.getPReadsIdBuffer()[p_in_id],
                         mapper.getPBamLinksBuffer()[p_out_id],
                         mapper.getOf(),
                         mapper.isOfNew(),
                         mapper.getFIOParms());
    }
    else
    {
        printAlignSamBam(mapper.getGenomes(),
                         mapper.getReads(),
                         mapper.getGenomesId(),
                         mapper.getReadsId(),
                         mapper.getBamRecords(),
                         mapper.getOf(),
                         mapper.isOfNew(),
                         mapper.getFIOParms());
    }
    return 0;
}
int print_align_sam (Mapper & mapper, int f_p_mapper, int p_in_id, int p_out_id)
{
    uint f_output_type_original = mapper.getFIOParms().f_output_type;
    fp_handler_.clear(mapper.getFIOParms().f_output_type);
    fp_handler_.setPrintSam(mapper.getFIOParms().f_output_type);
    print_align_sam_bam(mapper,f_p_mapper, p_in_id, p_out_id);
    mapper.getFIOParms().f_output_type = f_output_type_original;
    return 0;
}
int print_align_bam_std (Mapper & mapper, int f_p_mapper, int p_in_id, int p_out_id)
{
    uint f_output_type_original = mapper.getFIOParms().f_output_type;
    fp_handler_.clear(mapper.getFIOParms().f_output_type);
    fp_handler_.setPrintBamStd(mapper.getFIOParms().f_output_type);
    print_align_sam_bam(mapper, f_p_mapper, p_in_id, p_out_id);
    mapper.getFIOParms().f_output_type = f_output_type_original;
    return 0;
}
int print_align_bam_pbsv (Mapper & mapper, int f_p_mapper, int p_in_id, int p_out_id)
{
    uint f_output_type_original = mapper.getFIOParms().f_output_type;
    fp_handler_.clear(mapper.getFIOParms().f_output_type);
    fp_handler_.setPrintBamPbsv(mapper.getFIOParms().f_output_type);
    print_align_sam_bam(mapper, f_p_mapper, p_in_id, p_out_id);
    mapper.getFIOParms().f_output_type = f_output_type_original;
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
    //std::cout << "file ====== "  << "\n";
    if (fp_handler_.isPrintApf(mapper.getPrintFlag()))
    {
        open_mapper_of (mapper, file1);
        print_cords_apf(mapper, f_p_mapper, p_in_id, p_out_id);
        close_mapper_of(mapper);
    }
    ///.gvf
    /*
    std::string file2 = mapper.getOutputPrefix() + ".gvf";
    open_mapper_of (mapper, file2);
    print_clips_gvf(mapper);
    close_mapper_of(mapper);
    */
    ///.sam
    if (fp_handler_.isPrintSam(mapper.getPrintFlag()))
    {
        if (!fm_handler_.isAlign(mapper.getMapFlag()))
        {
            //formatCordsBam4Print(mapper, f_p_mapper);
        }
        std::string file3 = mapper.getOutputPrefix() + ".sam";
        open_mapper_of (mapper, file3);
        print_align_sam(mapper, f_p_mapper, p_in_id, p_out_id);
        close_mapper_of(mapper);
    }

    ///.bam
    if (fp_handler_.isPrintBamStd(mapper.getPrintFlag()))
    {
        if (!fm_handler_.isAlign(mapper.getMapFlag()))
        {
            //formatCordsBam4Print(mapper, f_p_mapper);
        }
        std::string file4 = mapper.getOutputPrefix() + ".bam";
        open_mapper_of (mapper, file4);
        print_align_bam_std(mapper, f_p_mapper, p_in_id, p_out_id);
        close_mapper_of(mapper);
    }
    //non-standard.bam for pbsv
    if (fp_handler_.isPrintBamPbsv(mapper.getPrintFlag()))
    {
        if (!fm_handler_.isAlign(mapper.getMapFlag()))
        {
            //formatCordsBam4Print(mapper, f_p_mapper);
        }
        std::string file5 = mapper.getOutputPrefix() + "_pbsv.bam";
        open_mapper_of (mapper, file5);
        print_align_bam_pbsv(mapper, f_p_mapper, p_in_id, p_out_id);
        close_mapper_of(mapper);
    }
    //add new file here

    mapper.setOfApp(); //set of_type to std::ios::app;
    return 0;
}
//read records from fin_pos and buckckets
int readRecords4FinPosbuckets(StringSet<CharString> & ids, 
        StringSet<String<Dna5> > & reads, 
        StringSet<String<short> >& buckets, 
        String<Position<SeqFileIn>::Type> & fin_pos,
        SeqFileIn & fin, uint rstr, uint rend, uint bucketId)
{
    unused(fin_pos);
    CharString tmp_id;
    String<Dna5> tmp_read;
    for (unsigned i = rstr; i < rend && !atEnd(fin); i++)
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

int readRecords2FinPosBuckets(StringSet<CharString> & ids, 
    StringSet<String<Dna5> > & reads, 
    StringSet<String<short> >& buckets, 
    String<Position<SeqFileIn>::Type> & fin_pos, SeqFileIn & fin, 
    uint rstr, uint rend)
{
    unused(buckets);
    CharString tmp_id;
    String<Dna5> tmp_read;
    if (empty(fin_pos)) 
    {
        appendValue(fin_pos, Position<SeqFileIn>::Type(0));
    }
    for (unsigned i = rstr; i < rend && !atEnd(fin); i++)
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
         StringSet<CharString> & reads_id,
         MapParms & mapParm,
         StringSet<String<uint64_t> > & cords_str,
         StringSet<String<uint64_t> > & cords_end,
         StringSet<String<CordInfo> > & cords_info,
         StringSet<String<uint64_t> > & clips,
         StringSet<String<Dna5> > & seqs,
         StringSet<String<BamAlignmentRecordLink> >& bam_records,
         String<GapParms> & gap_parms,
         uint f_map,   //control flags
         uint threads,
         int cord_size,
         int p1)
{
    unused(reads_id);
    unused(p1);
    uint thd_min_read_len = 200;
    //todo::tune the cordThr try to merge cords of blocks 
    MapParms complexParm = mapParm;
    complexParm.alpha = complexParm.alpha2;
    complexParm.listN = complexParm.listN2;
    String<uint64_t> gap_len;
    String<uint64_t> red_len;
    resize (gap_len, threads, 0);
    resize (red_len, threads, 0);
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
    unsigned size2 = length(reads) / threads;
    unsigned ChunkSize = size2;
    String<Dna5> comStr;
    Anchors anchors;
    String<uint64_t> crhit;
    StringSet<String<uint64_t> > cordsTmp;  //cords_str
    StringSet<String<uint64_t> > cordsTmp2; //cords_end
    StringSet<String<CordInfo> > cords_info_tmp;
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
    resize(cords_info_tmp, ChunkSize);
    resize(clipsTmp, ChunkSize);
    resize(bam_records_tmp, ChunkSize);
    resize(f1, 2);
    f1[0].init(f2[0].fs_type);
    f1[1].init(f2[0].fs_type);
    unsigned c = 0;
    CordsParms cords_parms;    
    String<UPair> apx_gaps; 
    #pragma omp for
    for (unsigned j = 0; j < length(reads); j++)
    {
        //double t1 = sysTime ();
        red_len[thd_id] += length(reads[j]);
        
        if (length(reads[j]) > thd_min_read_len)
        {
            _compltRvseStr(reads[j], comStr);
            createFeatures(begin(reads[j]), end(reads[j]), f1[0]);
            createFeatures(begin(comStr), end(comStr), f1[1]);
            apxMap(index, reads[j], anchors, crhit, f1, f2, apx_gaps, cordsTmp[c], cordsTmp2[c], cords_info_tmp[c], f_chain, mapParm.pm_g, mapParm.pm_pmp);
            if (fm_handler_.isMapGap(f_map))
            {
                mapGaps(seqs, reads[j], comStr, cordsTmp[c], cordsTmp2[c], clipsTmp[c], apx_gaps, f1, f2, gap_parms[thd_id]);
                reformCords(cordsTmp[c], cordsTmp2[c], &reformCordsDxDy1, cords_parms);
            }
            if (fm_handler_.isAlign(f_map))
            {
                alignCords(seqs, reads[j], comStr, cordsTmp[c], cordsTmp2[c], bam_records_tmp[c]);
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
    if (!fm_handler_.isAlign(f_map))
    {
        uint64_t thd_large_X = 8000;
        int f_parallel = 1;
        clear(bam_records);
        cords2BamLink(cords_str, cords_end, cords_info, bam_records, reads, cord_size, thd_large_X, threads, f_parallel);
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
    StringSet<std::string> files_prefix;
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
                 mapper.getCordsInfo(),
                 mapper.getClips(),
                 mapper.getGenomes(),
                 mapper.getBamRecords(),
                 mapper.getGapParms(),
                 mapper.getMapFlag(),
                 mapper.getThreads(), 
                 mapper.getCordSize(),
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
        appendValue (files_prefix, mapper.getOutputPrefix());
        close(rFile);
    }
    serr.print_message("--Write results to disk 100%", 0, 1, cerr);
    for (uint i = 0; i < length(files_prefix); i++)
    {
        serr.print_message("Result files: ", 2, 0, cerr);
        if(fp_handler_.isPrintApf(mapper.getPrintFlag()))
        {
            serr.print_message("\033[1;31m" + files_prefix[i] + ".apf\033[0m ", 0, 0, cerr);
        }
        //serr.print_message("\033[1;31m" + file2s[i] + "\033[0m ", 0, 0, cerr);
        if (fp_handler_.isPrintSam(mapper.getPrintFlag()))
        {
            serr.print_message("\033[1;31m" + files_prefix[i] + ".sam\033[0m ", 0, 0, cerr);
        }
        if (fp_handler_.isPrintBamStd(mapper.getPrintFlag()))
        {
            serr.print_message("\033[1;31m" + files_prefix[i] + ".bam\033[0m ", 0, 0, cerr);
        }
        if (fp_handler_.isPrintBamPbsv(mapper.getPrintFlag()))
        {
            serr.print_message("\033[1;31m" + files_prefix[i] + "_pbsv.bam\033[0m ", 0, 0, cerr);
        }
        serr.print_message(" ", 16, 1, cerr);  
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
            StringSet<String<CordInfo> > & cords_info,
            StringSet<String<uint64_t> > & clips,
            StringSet<String<Dna5> > & seqs,
            StringSet<String<BamAlignmentRecordLink> > & bam_records,
            uint f_map,   //control flags
            uint threads,
            int p1)
{
    unused(clips);
    unused(seqs);
    unused(bam_records);
    unused(f_map);
    unused(p1);
    uint thd_min_read_len = 200;
    //todo::tune the cordThr try to merge cords of blocks 
    //MapParms complexParm = mapParm;
    String<uint64_t> gap_len;
    String<uint64_t> red_len;
    resize (gap_len, threads, 0);
    resize (red_len, threads, 0);
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
        //double t1 = sysTime ();
        red_len[thd_id] += length(reads[j]);
        if (length(reads[j]) > thd_min_read_len)
        {
            _compltRvseStr(reads[j], comStr);
            createFeatures(begin(reads[j]), end(reads[j]), f1[0]);
            createFeatures(begin(comStr), end(comStr), f1[1]);
            apxMap(index, reads[j], anchors, crhit, f1, f2, apx_gaps, cordsTmp[c], cordsTmp2[c], cords_info[c], f_chain, mapParm.pm_g, mapParm.pm_pmp);
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
    for (unsigned i = 0; i < length(mapper.getRPaths()); i++)
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
                 mapper.getCordsInfo(),
                 mapper.getClips(),
                 mapper.getGenomes(),
                 mapper.getBamRecords(),
                 mapper.getMapFlag(),
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

MapParms::MapParms():
        blockSize(base_block_size_),
        delta(base_delta_),
        threshold(base_threshold_),
        kmerStep(base_kmer_step_),
        shapeLen(base_shape_len_),
        sensitivity(0),
        anchorDeltaThr(),
        minReadLen(1000),
        listN(1),
        listN2(1),
        alpha(base_alpha_),
        alpha2(0.5),
        anchorLenThr(0.02),                  // anchors with lenghth > this parameter is pushed into the queue
        rcThr(0.8),                        // when max anchors in the queue with length < this parameters, reverse complement search will be conducted
        cordThr(0.8),
        senThr(0.8),
        clsThr(0.1)
    {
        GlobalParms();
        PMPParms();
    }

MapParms::MapParms(unsigned bs, unsigned dt, unsigned thr, 
            unsigned ks, unsigned sl, unsigned st,
            unsigned ad, unsigned mr, unsigned listn,
            unsigned listn2,
            float ap, float ap2, float alt, float rt, float ct, float sent, float clst):
        blockSize(bs),
        delta(dt),
        threshold(thr),
        kmerStep(ks),
        shapeLen(sl),
        sensitivity(st),
        anchorDeltaThr(ad),
        minReadLen(mr),
        listN(listn),
        listN2(listn2),
        alpha(ap),
        alpha2(ap2),
        anchorLenThr(alt),                  // anchors with lenghth > this parameter is pushed into the queue
        rcThr(rt),                        // when max anchors in the queue with length < this parameters, reverse complement search will be conducted
        cordThr(ct),
        senThr(sent),
        clsThr(clst)
        {
            GlobalParms();
            PMPParms();
        }

MapParms::MapParms(MapParms & parm):
        blockSize(parm.blockSize),
        delta(parm.delta),
        threshold(parm.threshold),
        kmerStep(parm.kmerStep),
        shapeLen(parm.shapeLen),
        sensitivity(parm.sensitivity),
        anchorDeltaThr(),
        minReadLen(parm.minReadLen),
        listN(parm.listN),
        listN2(parm.listN2),
        alpha(parm.alpha),
        alpha2(parm.alpha2),
        anchorLenThr(parm.anchorLenThr),
        rcThr(parm.rcThr),
        cordThr(parm.cordThr),
        senThr(parm.senThr),
        clsThr(parm.clsThr)
        {
            GlobalParms();
            PMPParms();
        }

void MapParms::print()
{
    std::cerr << "blockSize " << blockSize << std::endl
              << "alpha " << alpha << std::endl
              << "alpha2 " << alpha2 << "\n"
              << "listN " << listN << "\n"
              << "listN2 " << listN2 << "\n"
              << "senThr " << senThr << "\n"
              << "delta " << delta << std::endl
              << "threshold " << threshold << std::endl
              << "kmerStep " << kmerStep << std::endl
              << "shapeLen " << shapeLen << std::endl
              //<<  "sensitivity " << sensitivity << "\n"
              << "anchorDeltaThr " << anchorDeltaThr << "\n"
              << "minReadLen " << minReadLen << "\n"
              << "anchorLenThr" << anchorLenThr << "\n"
              << "rcThr " << rcThr << "\n"
              << "cordThr" << cordThr << "\n";
}
