//#include <gtest/gtest.h>
#include <iostream>
#include <fstream>
#include <ctime>
#include <seqan/basic.h>
#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>
#include "base.h"
#include "shape_extend.h"
//#include "args_parser.h"
#include "cords.h"
#include "index_util.h"
#include "pmpfinder.h"
//#include "chain_map.h"
#include "gap.h"
#include "align_interface.h"
#include "f_io.h"
//#include "mapper.h"

using namespace seqan;


int main(int argc, const char ** argv)
{
    std::cerr << "Linear unit test[]" << std::endl;
    Options options;
    //LShape shape(20);
    /*
    seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;
        */

    SeqFileIn fin (toCString("1.fa"));
    LShape shape (20);
    String<Dna5> ss = "ACGT";
    hashInit (shape, begin(ss));
    hashNext (shape, begin(ss) + 1);

    //std::cout << "done\n";
    return 0;
}
