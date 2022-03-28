#include <seqan/arg_parse.h>
#include <iostream>
#include <fstream>
#include <ctime>
#include "args_parser.h"
#include "cords.h"
#include "pmpfinder.h"
#include "gap.h"
#include "align_interface.h"
#include "mapper.h"
#include "test_units.h"

using namespace seqan; 
using std::cerr;

int main(int argc, const char ** argv)
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
    omp_set_num_threads(mapper.getThreads());
    std::cerr << "running test1\n";
    //check_index1 (mapper.getGenomes(), mapper.index(), mapper.getThreads());
    //check_index2 (mapper.getGenomes(), mapper.index(), mapper.getThreads());
    check_apx_feature(mapper.getGenomes(), mapper.getThreads());


    std::cerr << "Time in sum[s] " << sysTime() - time << std::endl;
    return 0;
}
