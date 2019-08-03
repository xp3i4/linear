#include <seqan/arg_parse.h>
#include "args_parser.h"
#include "pmpfinder.h"
#include "mapper.h"
using namespace seqan; 


int process1 (Mapper & mapper, int p1)
{
    StringSet<FeaturesDynamic> f2;
    createFeatures(mapper.getGenomes(), f2, mapper.getFeatureType(), mapper.thread());
    omp_set_num_threads(mapper.thread());
    mapper.createIndex(0, length(mapper.getGenomes()), false); // true: destruct genomes string to reduce memory footprint
    return map (mapper, 0, length(mapper.getGenomes()), f2, p1);
}

int process2(Mapper & mapper, Options & options, int p1)
{
    omp_set_num_threads(mapper.thread());
    StringSet<FeaturesDynamic> f2;
    createFeatures(mapper.getGenomes(), f2, mapper.getFeatureType(), mapper.thread());
    Options filter_options = options;
    filter_options.aln_flag = 0;
    filter_options.sam_flag = 0;
    filter_options.apx_chain_flag = 0;
    filter_options.index_t = 2;
    filter_options.feature_t = 1;
    mapper.loadOptions (filter_options);
    mapper.createIndex(0, length(mapper.getGenomes()), false); 
    filter (mapper, f2, p1);
    //mapper.clearIndex();
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
    process1 (mapper, options.p1);
    std::cerr << "Time in sum[s] " << sysTime() - time << std::endl;
    return 0;
}
