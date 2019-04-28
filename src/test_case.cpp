//#include <gtest/gtest.h>
#include <seqan/arg_parse.h>
#include <iostream>
#include <fstream>
#include <ctime>
#include "mapper.h"
#include "args_parser.h"
#include "pmpfinder.h"
#include "chain_map.h"
#include "gap.h"
#include "align_interface.h"

using namespace seqan;

int testCreateDIndx(DIndex & index, StringSet<String<Dna5> > & seqs)
{
    createDIndex (seqs, index, 0, 0);
}
int testQueryDIndx (DIndex & index, StringSet<String<Dna5> > & seqs)
{
    LShape & shape = index.getShape();
    int64_t pre_xval = ~0;
    bool find_f = false;
    int count_f = 0;
    createDIndex (seqs, index, 0, 0);
    std::cerr << "[testQueryDIndx]::create index done\n";
    for (int j = 0 ; j < length(seqs); j++)
    {
        hashInit(shape, begin(seqs[j]));
        for (int i = 0; i < length(seqs[j]) - shape.span; i++) 
        {
            hashNext(shape, begin(seqs[j]) + i);
            int64_t str = queryHsStr(index, shape.XValue);
            int64_t end = queryHsEnd(index, shape.XValue);

            for (int64_t j = str; j < end; j++)
            {
                int64_t val = index.getHs()[j];
                if (get_cord_x(val) == i)
                {
                    find_f = true;
                    break;
                }
            }
            if (!find_f)
            {
                ++count_f;
            }
            find_f = false;
        }
    }
    std::cout << "count_f " << count_f << "\n";
}
int testMinimizerDIndx (DIndex & index, StringSet<String<Dna5> > & seqs)
{
    LShape & shape = index.getShape();
    int64_t pre_xval = ~0;
    bool find_f = false;
    int count_f = 0;
    createDIndex (seqs, index, 0, 0);
    std::cerr << "[testQueryDIndx]::create index done\n";
    String<int> counts;
    resize(counts, 1000000, 0);
    for (int i = 0; i < length(index.getDir()) - 1; i++)
    {
        int j = index.getDir()[i + 1] - index.getDir()[i];
        ++counts[j];

    }
    for (int i = 0; i < length(counts); i++)
    {
        if (counts[i] != 0)
        {
            std::cout << "size buckets " << i << " " << (float) counts[i] / length(index.getDir()) << "\n";
        }
    }
}
int main(int argc, const char ** argv)
{
    std::cerr << "Linear unit test[]" << std::endl;
    Options options;
    seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;
    Mapper mapper(options);
    DIndex index_d(22);
    //testCreateDIndx(index_d, mapper.genomes());
    testQueryDIndx(index_d, mapper.genomes());
    //testMinimizerDIndx(index_d, mapper.genomes());
    std::cout << "Lshape.len " << index_d.getShape().weight << "\n";


    //testing::InitGoogleTest(&argc, argv);
    //return RUN_ALL_TESTS();
    return 0;
}
