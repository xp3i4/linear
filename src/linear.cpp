#include <seqan/arg_parse.h>
#include "args_parser.h"
#include "pmpfinder.h"
#include "mapper.h"
using namespace seqan; 


int process1 (Options & options, int p1)
{
    StringSet<FeaturesDynamic> f2;
    StringSet<String<short> > emptyBuckets; 

    omp_set_num_threads(mapper.thread());
    createFeatures(mapper.getGenomes(), f2, mapper.getFeatureType(), mapper.thread());
    mapper.createIndex(0, length(mapper.getGenomes()), false); 
    return map (mapper, f2, emptyBuckets, 0, 0, p1);
    //return map (mapper, 0, length(mapper.getGenomes()), f2, p1);
}

int process2(Mapper & mapper, Options & options, int p1)
{
    //init filter mapper
    StringSet<FeaturesDynamic> f2;
    StringSet<String<short> > buckets;
    Options filter_options = options;
    filter_options.aln_flag = 0;
    filter_options.sam_flag = 0;
    filter_options.apx_chain_flag = 0;
    filter_options.index_t = 2;
    filter_options.feature_t = 1;
    mapper.loadOptions (filter_options);

    //filter chr bucket 
    omp_set_num_threads(mapper.thread());
    mapper.createIndex(0, length(mapper.getGenomes()), false); 
    createFeatures(mapper.getGenomes(), f2, mapper.getFeatureType(), mapper.thread());
    filter (mapper, f2, buckets, p1);

    //clear filer index
    mapper.clearIndex();

    //map for each chr bucket
    mapper.loadOptions(options);
    for (unsigned i = 0; i < length(mapper.getGenomes()); i++)
    {
        dout << "mgs" << i << "\n";
        serr.print_message("Map::\033[1;35mref:", 0, 0, std::cerr);
        serr.print_message(i + 1, 0, 0, std::cerr);
        serr.print_message(" " + _2stdstring(mapper.getGenomesId()[i]), 0, 0, std::cerr);
        serr.print_message("\033[0m        ", 0, 1, std::cerr);
        mapper.createIndex(i, i + 1, false); 
        map (mapper, f2, buckets, i, 1, p1);
        mapper.clearIndex();
        std::cerr << "\n";
    }
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
    uint thd_g_size = 300 << 20; //300M bases 
    if (lengthSum(mapper.getGenomes ()) > thd_g_size)
    {
        process2 (mapper, options, p1);
    }
    else
    {
        process1 (mapper, p1);
    }
    //process2 (options, options.p1);
    std::cerr << "Time in sum[s] " << sysTime() - time << "      \n";
    return 0;
}
