#include <seqan/arg_parse.h>
#include "args_parser.h"
#include "pmpfinder.h"
#include "parallel_io.h"
#include "mapper.h"
using namespace seqan; 

int process1 (Mapper & mapper, Options & options, int p1)
{
    StringSet<String<short> > empty_buckets; 
    String<Position<SeqFileIn>::Type>   empty_fin_pos; 

    omp_set_num_threads(mapper.getThreads());
    createFeatures(mapper.getGenomes(), mapper.getGenomesFeatures(), 
        mapper.getFeatureType(), mapper.getThreads());
    mapper.createIndex(0, length(mapper.getGenomes()), false); 
    return map (mapper, empty_buckets, empty_fin_pos, 0, 0, p1);
}
/**
 * Supposed to be called when lengthSum(mapper.getGenomes ()) > 300M
 */
int process2(Mapper & mapper, Options & options, int p1)
{
    //init filter mapper
    StringSet<String<short> > buckets;
    String<Position<SeqFileIn>::Type> fin_pos; 
    Options filter_options = options;
    filter_options.aln_flag = 0;
    filter_options.sam_flag = 0;
    filter_options.apx_chain_flag = 0;
    filter_options.index_t = typeDIx;
    filter_options.feature_t = 1;
    mapper.loadOptions (filter_options);

    //filter chr bucket 
    omp_set_num_threads(mapper.getThreads());
    mapper.getIndex().setMIndex(); //enable filter index && disable mapper index (DIndex for filter by default)
    mapper.createIndex(0, length(mapper.getGenomes()), false); 
    createFeatures(mapper.getGenomes(), mapper.getGenomesFeatures(), mapper.getFeatureType(), mapper.getThreads());
    filter (mapper, buckets, fin_pos, p1);
    
    //clear filter index
    //mapper.clearIndex();

    //map for each chr bucket
    mapper.loadOptions(options);
    mapper.getIndex().setMHIndex(); //disable filter index && enable mapper index
    bool f_io_append = false;       //flag if result append to or open as new 
    for (unsigned i = 0; i < length(mapper.getGenomes()); i++)
    {
        serr.print_message("Map::\033[1;35mref:", 0, 0, std::cerr);
        serr.print_message(i + 1, 0, 0, std::cerr);
        serr.print_message(" " + _2stdstring(mapper.getGenomesId()[i]), 0, 0, std::cerr);
        serr.print_message("\033[0m        ", 0, 1, std::cerr);
        //mapper.
        mapper.createIndex(i, i + 1, false); 
        map (mapper, buckets, fin_pos, i, 1, p1, f_io_append);
        f_io_append = true;
        std::cerr << "\n";
    }
    return 0;
}
/**
 * ----------  Process 3  ----------
 * Dynamic balancing with parallel I/O 
 */
int process3(Mapper & mapper, Options & options, int p1)
{
    int num_buffers = std::max(30, (int)(mapper.getThreads() * 1.5));
    P_Parms p_parms(1, 1, 1, 100);
    mapper.initBuffers(num_buffers, num_buffers, p_parms);
    omp_set_num_threads(mapper.getThreads());
    omp_set_num_threads(mapper.getThreads());
    createFeatures(mapper.getGenomes(), mapper.getGenomesFeatures(), mapper.getFeatureType(), mapper.getThreads());
    mapper.createIndex(0, length(mapper.getGenomes()), false); 
    mapper.getPTasks().startRunTasks();
    #pragma omp parallel
    {
        p_ThreadProcess(mapper, p_parms, omp_get_thread_num());
    }
    //mapper.getTasks().printInfos();
    return 0;
}

int main(int argc, char const ** argv)
{
    double time = sysTime();
    Options options;
    options.printRunInfo();
    std::cerr << std::fixed << std::setprecision(2);
    seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;
    Mapper mapper(options);
    if (options.bal_flag)
    {
        process3 (mapper, options, 0);
    }
    else
    {
        process1 (mapper, options, 0);
    }
    std::cerr << "Time in sum[s] " << sysTime() - time << "      \n";
    return 0;
}
