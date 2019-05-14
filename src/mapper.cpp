#include <seqan/arg_parse.h>
#include <iostream>
#include <fstream>
#include <ctime>
#include "mapper.h"
#include "args_parser.h"
#include "cords.h"
#include "pmpfinder.h"
#include "chain_map.h"
#include "gap.h"
#include "align_interface.h"

using namespace seqan; 
using std::cerr;

Mapper::Mapper(Options & options):
               record(options),
               index_dynamic(getGenomes()), 
               of(toCString(options.getOutputPath()))
{
    r_paths = options.r_paths;
    g_paths = options.g_paths; 
    loadRecords(getGenomes(), getGenomesId(), g_paths);
    //for (auto g : getGenomes())
    //std::cout << g << "\n";
    switch (options.sensitivity)
    {
        case 0: 
        {
            parm = parm0; //normal
            break;
        }
        case 1:
        {
            parm =  parm1; //fast
            break;
        }
        case 2:
        {
            parm = parm2; //sensitive
            break;
        }
    }
    _thread = options.thread;
    if (options.index_t == 1)
    {
        index_dynamic.setHIndex();
    }
    else if (options.index_t == 2)
    {
        index_dynamic.setDIndex();
    }
    if (options.feature_t == 1)
    {
        feature_type = typeFeatures1_32; 
    }
    else if (options.feature_t == 2)
    {
        feature_type = typeFeatures2_48;
    }
    of_type = OF_NEW;
}
int Mapper::createIndex(bool efficient)
{
    createIndexDynamic(getGenomes(), index_dynamic, _thread, efficient);
//  createDIndex(genomes(), dIndex, thd_min_step, thd_max_step, _thread);
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

int print_align_sam (Mapper & mapper)
{
    std::string filePath = mapper.getOutputPrefix() + ".sam";
    mapper.getOf().open(toCString(filePath));
    print_align_sam_header_(mapper.getGenomesId(), 
                            mapper.getGenomes(),
                            mapper.getOf()
                           );
    print_align_sam_record_(mapper.getBamRecords(),
                            mapper.getCords(),
                            mapper.getReadsId(),
                            mapper.getGenomesId(),
                            mapper.getOf()
                           ); 
    mapper.getOf().close();
}

int print_clips_gvf(Mapper & mapper)
{
    print_clips_gvf_(mapper.getClips(), mapper.getReadsId(), mapper.getGenomesId(), mapper.getOf(), mapper.getOutputPrefix());
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
 * Print apf sam and gvf
 * Append to the end of the file  
 */
int print_mapper_results(Mapper & mapper, 
                         std::string outputPrefix,
                         int f_prints = 0 /*print options control flags*/)
{
    double time = sysTime();
    if (mapper.getOutputPrefix() != outputPrefix)
    {
        mapper.getOutputPrefix() = outputPrefix;
        mapper.setOfNew();
    }
    else
    {
        mapper.setOfApp();
    }
    std::string file1 = mapper.getOutputPrefix() + ".apf";
    open_mapper_of (mapper, file1);
    print_cords_apf(mapper.getCords(), 
                    mapper.getGenomes(),
                    mapper.getReads(),
                    mapper.getGenomesId(),
                    mapper.getReadsId(),
                    mapper.getOf()
                );
    close_mapper_of(mapper);

    std::string file2 = mapper.getOutputPrefix() + ".gvf";
    open_mapper_of (mapper, file2);
    print_clips_gvf(mapper);
    close_mapper_of(mapper);

    std::string file3 = mapper.getOutputPrefix() + ".sam";
    open_mapper_of (mapper, file3);
    print_align_sam(mapper);
    close_mapper_of(mapper);

    mapper.setOfApp(); //set of_type to std::ios::app;
    return 0;
}

int map_(IndexDynamic & index,
         StringSet<FeaturesDynamic > & f2,
         StringSet<String<Dna5> > & reads,
         MapParm & mapParm,
         StringSet<String<uint64_t> > & cords,
         StringSet<String<uint64_t> > & clips,
         StringSet<String<Dna5> > & seqs,
         StringSet<String<BamAlignmentRecordLink> >& bam_records,
         unsigned & threads,
         int p1
        )
{
    typedef String<Dna5> Seq;
    float senThr = mapParm.senThr / window_size;  //map for 2 roun if cords cover len <
    float cordThr = mapParm.cordThr / window_size; //cords cover length < are aborted
    MapParm complexParm = mapParm;
    complexParm.alpha = complexParm.alpha2;
    complexParm.listN = complexParm.listN2;
    String<uint64_t> gap_len;
    String<uint64_t> red_len;
    resize (gap_len, threads, 0);
    resize (red_len, threads, 0);
#pragma omp parallel
{
    unsigned size2 = length(reads) / threads;
    unsigned ChunkSize = size2;
    Seq comStr;
    Anchors anchors;
    String<uint64_t> crhit;
    StringSet<String<uint64_t> >  cordsTmp;
    StringSet<FeaturesDynamic> f1;
    StringSet<String<uint64_t> > clipsTmp;
    StringSet<String<BamAlignmentRecordLink> > bam_records_tmp;
    unsigned thd_id =  omp_get_thread_num();
    if (thd_id < length(reads) - size2 * threads)
    {
        ChunkSize = size2 + 1;
    }
    resize(cordsTmp, ChunkSize);
    resize(clipsTmp, ChunkSize);
    resize(bam_records_tmp, ChunkSize);
    resize(f1, 2);
    f1[0].fs_type = f2[0].fs_type;
    f1[1].fs_type = f2[0].fs_type;
    unsigned c = 0;
    
    String<uint64_t> g_hs;
    String<uint64_t> g_anchor;
    resize (g_hs, 1ULL << 20);
    resize (g_anchor, 1ULL<<20);
    #pragma omp for
    for (unsigned j = 0; j < length(reads); j++)
    {
        if (length(reads[j]) >= mapParm.minReadLen)
        {
            red_len[thd_id] += length(reads[j]);
            //std::cout << "[]::rawmap::j " << j <<"\n";
            float cordLenThr = length(reads[j]) * cordThr;
            _compltRvseStr(reads[j], comStr);
            createFeatures(begin(reads[j]), end(reads[j]), f1[0]);
            createFeatures(begin(comStr), end(comStr), f1[1]);
            anchors.init(1);
            clear(crhit);
            mnMapReadList(index, reads[j], anchors, mapParm, crhit);
            path_dst(begin(crhit), end(crhit), f1, f2, cordsTmp[c], cordLenThr);
            if (_DefaultCord.getMaxLen(cordsTmp[c]) < length(reads[j]) * senThr)
            {
                clear(cordsTmp[c]);
                anchors.init(1);
                clear(crhit);
                mnMapReadList(index, reads[j], anchors, complexParm, crhit);
                path_dst(begin(crhit), end(crhit), f1, f2, cordsTmp[c], cordLenThr);
            }   
            gap_len[thd_id] += mapGaps(seqs, reads[j], comStr, cordsTmp[c], g_hs, g_anchor, clipsTmp[c], f1, f2, p1, window_size);
            //align_cords(seqs, reads[j], comStr, cordsTmp[c], bam_records_tmp[c], p1);
        }   
        c += 1;
    } 
    #pragma omp for ordered
    for (unsigned j = 0; j < threads; j++)
        #pragma omp ordered
        {
            append(cords, cordsTmp);
            append(clips, clipsTmp);
            append(bam_records, bam_records_tmp);
        }
    
}
    return 0;
}

/**
 * Map control flow
 */
int map(Mapper & mapper, int p1)
{
    int thd_buffer_size = 10000; // every thd_buffer_size reads triggers print results.
    //serr.print_message("=>Create index", 0, 1, std::cerr);
    mapper.createIndex(false); // true: destruct genomes string to reduce memory footprint
    StringSet<FeaturesDynamic> f2;
    createFeatures(mapper.getGenomes(), f2, mapper.getFeatureType(), mapper.thread());
    std::cout << "mapf " << mapper.getFeatureType() << "\n";
    unsigned blockSize = 50000;
    StringSet<String<char> > dotstatus;
    resize(dotstatus, 3);
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
        std::string outputPrefix = getFileName(path);
        //outputPrefix = getFileName(outputPrefix, ".", 1);
        unsigned k = 1;
        while (!atEnd(rFile))
        {
            double time1 = sysTime();
            serr.print_message("=>Map::file_I/O", 0, 0, std::cerr);
            serr.print_message(k, 0, 2, std::cerr);
            try
            {
                readRecords(mapper.getReadsId(), mapper.getReads(), rFile, blockSize);
            }
            catch (Exception const & e)
            {

            }
            //serr.print_message("", 50, 2, std::cerr); 
            serr.print_message("=>Map::mapping  block", 0, 0, std::cerr);
            serr.print_message(k, 0, 0, std::cerr);
            serr.print_message("Size ", 0, 0, std::cerr);
            serr.print_message(unsigned(length(mapper.getReads())), 0, 2, std::cerr);
            time1 = sysTime() - time1;
            double time2 = sysTime();
            map_(mapper.index(), 
                 f2, 
                 mapper.getReads(), 
                 mapper.mapParm(), 
                 mapper.getCords(), 
                 mapper.getClips(),
                 mapper.getGenomes(),
                 mapper.getBamRecords(),
                 mapper.thread(), 
                 p1);
            time2 = sysTime() - time2;
            std::cerr <<  "--Map::file_I/O+Map block "<< k << " Size " << length(mapper.getReads()) << " Elapsed Time[s]: file_I/O " << time1 << " map "<< time2 << "\n";
            //if (buffer_size > thd_buffer_size)
            //{
            serr.print_message("=>Write results to disk", 0, 2, std::cerr);
            print_mapper_results(mapper, outputPrefix);
            //}
            clear (mapper.getCords());
            clear (mapper.getClips());
            clear (mapper.getBamRecords());
            clear (mapper.getReads());
            clear (mapper.getReadsId());
            clear (mapper.getClips());
            clear (mapper.getBamRecords());
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

int main(int argc, char const ** argv)
{
    double time = sysTime();
    //(void)argc;
    Options options;
    options.versions = "1.8.2";
    std::cerr << "["<< options.versions
              << "]\nEncapsulated: Mapping reads efficiently" << std::endl;
    seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;
    Mapper mapper(options);
    omp_set_num_threads(mapper.thread());
    map(mapper, options.p1);

    std::cerr << "Time in sum[s] " << sysTime() - time << std::endl;
    return 0;
}
