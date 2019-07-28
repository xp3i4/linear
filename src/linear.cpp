#include <seqan/arg_parse.h>
#include "args_parser.h"
#include "mapper.h"
using namespace seqan; 


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
    map(mapper, options.p1);
    std::cerr << "Time in sum[s] " << sysTime() - time << std::endl;
    return 0;
}
