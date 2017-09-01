#include <seqan/seq_io.h>
#include <seqan/stream.h>
#include <seqan/index.h>
#include <seqan/store.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <seqan/basic.h>
#include <bitset>
#include <climits>
#include <seqan/arg_parse.h>

#include "shape_extend.h"
#include "index_extend.h"
//#include "pmpfinder.h"
#include "pacmapper.h"

using namespace seqan; 


    seqan::ArgumentParser::ParseResult
    parseCommandLine(Options & options, int argc, char const ** argv)
    {
        // Setup ArgumentParser.
        seqan::ArgumentParser parser("pacMapper");
        // Set short description, version, and date.
        setShortDescription(parser, "Alignment of SMRT sequencing read");
        setVersion(parser, "1.0");
        setDate(parser, "May 2017");

        // Define usage line and long description.
        addUsageLine(parser,
                     "[\\fIOPTIONS\\fP] \"\\fIread.fa\\fP\" \"\\fIgnome.fa\\fP\"");
        addDescription(parser,
                       "Program for mapping raw SMRT sequencing reads to reference genome.");

        // Argument.
        addArgument(parser, seqan::ArgParseArgument(
            seqan::ArgParseArgument::INPUT_FILE, "read"));

        addArgument(parser, seqan::ArgParseArgument(
            seqan::ArgParseArgument::INPUT_FILE, "genome"));

        addSection(parser, "Mapping Options");
        addOption(parser, seqan::ArgParseOption(
            "o", "output", "choose output file.",
             seqan::ArgParseArgument::STRING, "STR"));
        
        // Add Examples Section.
        addTextSection(parser, "Examples");
        addListItem(parser,
                    "\\fBpacMapper\\fP \\fB-U\\fP \\fIchr1.fa reads.fa\\fP",
                    "Print version of \"rd\"");
        addListItem(parser,
                    "\\fBpacMapper\\fP \\fB-L\\fP \\fB-i\\fP \\fI3\\fP "
                    "\\fIchr1.fa reads.fa\\fP",
                    "Print \"\" with every third character "
                    "converted to upper case.");

        // Parse command line.
        seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

        if (res != seqan::ArgumentParser::PARSE_OK)
            return res;

        getOptionValue(options.oPath, parser, "output");

        seqan::getArgumentValue(options.rPath, parser, 0);
        seqan::getArgumentValue(options.gPath, parser, 1);

        //std::cout << options.rFile << " " << options.gFile << " " << options.oFile<< std::endl;

        return seqan::ArgumentParser::PARSE_OK;
    }

/*
void pmpfinder(StringSet<String<Dna5> > & gnome, StringSet<String<Dna5> > & reads)// uint64_t & gbegin, uint64_t &gend)
{
    Shape<Dna5, UngappedShape> > shape;
    String<uint64_t> readFtr;
    uint64_t sum = 0;
    double time = sysTime(); 
    for (unsigned j = 0; j < length(reads); j++) 
    {   
        hashInit(shape, begin(reads[j]));
        for(uint64_t k = 0; k < length(reads[j]); k++)     
        {
            sum ^= hashNext(shape, begin(reads[j]) + k);
        }
    }   
    std::cout << sum << " " <<sysTime()-time << std::endl;
}
*/

int main(int argc, char const ** argv)
{
//    ArgumentParser parser;
//    Options options; 
//    setupArgumentParser(parser, options);
//
//    ArgumentParser::ParseResult res = parseCommandLine(options, parser, argc, argv);

    //time = sysTime();
    (void)argc;
    SeqFileIn gFile(toCString(argv[1]));
    SeqFileIn rFile(toCString(argv[2]));
    StringSet<CharString> ids_r;
    StringSet<CharString> ids_g;
    //StringSet<DnaString> reads;
    //StringSet<DnaString> genome;
    StringSet<String<Dna5> > reads5;
    StringSet<String<Dna5> > genome5;
    //MapParm mapParm;

    //readRecords(ids_r, reads, rFile);
    //readRecords(ids_g, genome, gFile);
    readRecords(ids_r, reads5, rFile);
    readRecords(ids_g, genome5, gFile);
    
    
//    std::fstream myfile(toCString(argv[3]), std::ios_base::in);  
//    String<int> s;
//    uint64_t a[10000];
//    int d,f;
//    String<uint64_t> as;
//    unsigned k =0;
//    while(k < 10000)
//    {
//        myfile >> d >> f;
//        
//        a[k] = ((uint64_t)(d / (float)16 + 0.5)) + ((uint64_t)(f/(float)16+0.5) << 32);
//        k++;
//    }
//    while(k < 10000)
//    {
//        myfile >> d >> f;
//        
//        appendValue(as, ((uint64_t)(d /16)<< 32) + (uint64_t)(f/16));
//        //appendValue(as, (((uint64_t)(d /(float)16+0.5))<< 32) + (uint64_t)(f/(float)16+0.5));
////        std::cout << a[k] << std::endl;
//        k++;
//    }

    //opKmerAnalysis5(genome, reads);
    //mnKmerMap(genome, reads);
    //setupArgumentParser(argc, argv);
    //mnKmerMap5(genome5, reads5);
    //mnMap(genome5, reads5, mapParm);
    map(genome5, reads5);
    //pmpfinder(genome5, reads5);
        //if (k % 100 == 0)
        //    std::cout << std::endl;
    //feTest(reads5, begin(genome5[0]) + 33688421, begin(genome5[0]) + 33695356);
    //feTest(reads5, begin(genome5[0]) + 121275915, begin(genome5[0]) + 121284556);
    //feTest(reads5, begin(genome5[0]) + 194279648, begin(genome5[0]) + 194286590);
//    pTest(reads5, begin(genome5[0]) + 33688421, begin(genome5[0]) + 33695356);
    //pTest(reads5, begin(genome5[0]) + 33688544, begin(genome5[0]) + 33695349);
    //feTest(reads5[0], begin(genome5[0]), end(genome5[0]));
//    pTest1(genome5[0]);
    //pTest2(reads5, begin(genome5[0]) + 33688544, begin(genome5[0]) + 33695349);
    //pTest2(reads5, begin(genome5[0]) + 33688544 - 127, begin(genome5[0]) + 33695349);
    
    //pTest3(reads5, genome5, begin(genome5[0]) + 33688544 - 300, begin(genome5[0]) + 33695349);
    //pTest4(reads5, genome5, ids_r, a);
    //dTest(reads5, genome5, as);

          // Parse the command line.
/*
    Options options;
    seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;
    std::cout << options.rFile << std::endl;
    std::cerr << options.gFile << std::endl;
    
*/
 
    return 0;
}
