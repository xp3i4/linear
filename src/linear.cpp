#include <seqan/arg_parse.h>
#include "args_parser.h"
#include "pmpfinder.h"
#include "mapper.h"
using namespace seqan; 


int process1 (Mapper & mapper, Options & options, int p1)
{
    StringSet<FeaturesDynamic> f2;
    StringSet<String<short> > empty_buckets; 
    String<Position<SeqFileIn>::Type>   empty_fin_pos; 

    omp_set_num_threads(mapper.getThreads());
    createFeatures(mapper.getGenomes(), f2, mapper.getFeatureType(), mapper.getThreads());
    mapper.createIndex(0, length(mapper.getGenomes()), false); 
    return map (mapper, f2, empty_buckets, empty_fin_pos, 0, 0, p1);
}

int process2(Mapper & mapper, Options & options, int p1)
{
    //init filter mapper
    StringSet<FeaturesDynamic> f2;
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
    createFeatures(mapper.getGenomes(), f2, mapper.getFeatureType(), mapper.getThreads());
    filter (mapper, f2, buckets, fin_pos, p1);
    //<<ddebug
    for (int i = 0; i < length(buckets); i++)
    {

        std::cout << "\nbuckets " << " ";
        for (int j = 0; j < length(buckets[i]); j++)
        {
            std::cout << buckets[i][j] << " " ;
        }
    }
    dout << "indextype" << mapper.getIndex().typeIx << "\n";
    //>>debug

    //clear filter index
    //mapper.clearIndex();

    //map for each chr bucket
    
    mapper.loadOptions(options);
    mapper.getIndex().setMHIndex(); //disable filter index && enable mapper index
    bool f_io_append = false;       //flag if result append to or open as new 
    for (unsigned i = 0; i < length(mapper.getGenomes()); i++)
    {
        dout << "mgs" << i << "\n";
        serr.print_message("Map::\033[1;35mref:", 0, 0, std::cerr);
        serr.print_message(i + 1, 0, 0, std::cerr);
        serr.print_message(" " + _2stdstring(mapper.getGenomesId()[i]), 0, 0, std::cerr);
        serr.print_message("\033[0m        ", 0, 1, std::cerr);
        //mapper.
        mapper.createIndex(i, i + 1, false); 
        map (mapper, f2, buckets, fin_pos, i, 1, p1, f_io_append);
        f_io_append = true;
        std::cerr << "\n";
    }
    return 0;
}

int process3(Mapper & mapper, Options & options, int p1)
{

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
        int p1 = 0;
        process2 (mapper, options, p1);
    }
    else
    {
        int p1 = 0; //temp var for test config
        process1 (mapper, options, p1);
    }
    std::cerr << "Time in sum[s] " << sysTime() - time << "      \n";
    return 0;
}
