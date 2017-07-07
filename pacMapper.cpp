#include <seqan/seq_io.h>
#include <seqan/stream.h>
#include <seqan/index.h>
#include <seqan/store.h>
#include <iostream>
#include <math.h>
#include <seqan/basic.h>
#include <bitset>
#include <climits>
#include <seqan/arg_parse.h>

#include "shape_pm.h"
#include "index_pm.h"
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

        getOptionValue(options.oFile, parser, "output");

        seqan::getArgumentValue(options.rFile, parser, 0);
        seqan::getArgumentValue(options.gFile, parser, 1);

        //std::cout << options.rFile << " " << options.gFile << " " << options.oFile<< std::endl;

        return seqan::ArgumentParser::PARSE_OK;
    }

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

    //opKmerAnalysis5(genome, reads);
    //mnKmerMap(genome, reads);
    //setupArgumentParser(argc, argv);
    //mnKmerMap5(genome5, reads5);
    //mnMap(genome5, reads5, mapParm);
    map(genome5, reads5);
    
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
