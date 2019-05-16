#include <seqan/arg_parse.h>
#include "base.h"
std::string CARTESIAN = "x";
seqan::ArgumentParser::ParseResult
parseCommandLine(Options & options, int argc, char const ** argv)
{
    // Setup ArgumentParser.
    seqan::ArgumentParser parser("Linear");
    // Set short description, version, and date.
    setShortDescription(parser, "Alignment of SMRT sequencing read");
    setVersion(parser, options.versions);
    setDate(parser, options.date);

    // Define usage line and long description.
    addUsageLine(parser,
                    "[\\fIOPTIONS\\fP] \"\\fIread.fa(stq)\\fP\" \"\\fIgnome.fa\\fP\"");
    addDescription(parser,
                    "Program for mapping raw SMRT sequencing reads to reference genome.");

    // Argument.
    addArgument(parser, seqan::ArgParseArgument(
        seqan::ArgParseArgument::INPUT_FILE, "read",true));
    setHelpText(parser, 0, "Reads file .fa(.gz), .fasta(.gz), .fq(.gz), .fastq(.gz) ");

/*
    addArgument(parser, seqan::ArgParseArgument(
        seqan::ArgParseArgument::INPUT_FILE, "genome"));
    setHelpText(parser, 1, "Reference file .fa(.gz), .fasta(.gz)");
*/
    addSection(parser, "Mapping Options");
    addOption(parser, seqan::ArgParseOption(
        "o", "output", "choose output file.",
            seqan::ArgParseArgument::STRING, "STR"));
    addOption(parser, seqan::ArgParseOption(
        "s", "sensitivity", "Sensitivity mode. -s 0 normal {DEFAULT} -s 1 fast  -s 2 sensitive",
            seqan::ArgParseArgument::INTEGER, "INT"));
    addOption(parser, seqan::ArgParseOption(
        "t", "thread", "Default -t 4",
            seqan::ArgParseArgument::INTEGER, "INT"));
    addOption(parser, seqan::ArgParseOption(
        "i", "index_type", "Default -idx 1",
            seqan::ArgParseArgument::INTEGER, "INT"
        ));
    addOption(parser, seqan::ArgParseOption(
        "f", "feature_type", "{1,2}. Default -f 2 (2-mer, 48bases)",
            seqan::ArgParseArgument::INTEGER, "INT"
        )); 
    addOption(parser, seqan::ArgParseOption(
        "g", "gap_len", "0 to turn off gap mapping module, set > 0 to map gaps whose length > this value",
            seqan::ArgParseArgument::INTEGER, "INT"
        )); 
    addOption(parser, seqan::ArgParseOption(
        "a", "aln_flag", "0 to turn of alignment module",
            seqan::ArgParseArgument::INTEGER, "INT"
        )); 



// mapping parameters for tunning 
    addOption(parser, seqan::ArgParseOption(
        "l1", "listn1", "mapping::listn1",
            seqan::ArgParseArgument::INTEGER, "INT"));     
    addOption(parser, seqan::ArgParseOption(
        "l2", "listn2", "mapping::listn2",
            seqan::ArgParseArgument::INTEGER, "INT"));   
    addOption(parser, seqan::ArgParseOption(
        "a1", "alpha1", "mapping::alpha1",
            seqan::ArgParseArgument::DOUBLE, "DOUBLE"));        
    addOption(parser, seqan::ArgParseOption(
        "a2", "alpha2", "mapping::alpha2",
            seqan::ArgParseArgument::DOUBLE, "DOUBLE"));  
    addOption(parser, seqan::ArgParseOption(
        "t1", "cordThr", "mapping::cordThr",
            seqan::ArgParseArgument::DOUBLE, "DOUBLE"));
    addOption(parser, seqan::ArgParseOption(
        "t2", "senThr", "mapping::senThr",
            seqan::ArgParseArgument::DOUBLE, "DOUBLE"));        
    addOption(parser, seqan::ArgParseOption(
        "p1", "par1", "options::p1",
            seqan::ArgParseArgument::INTEGER, "INT")); 
        
    // Add Examples Section.
////////////////////////
    addTextSection(parser, "Examples");
    addListItem(parser,
                "\\fBlinear \\fP \\fIreads.fa genomes.fa\\fP",
                "raw map reads.fa to genomes.fa");
    addTextSection(parser, "Examples");
    addListItem(parser,
                "\\fBlinear\\fP \\fP-a \\fIreads.fa genomes.fa\\fP",
                "align reads.fa to genomes.fa");

    // Parse command line.
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

    if (res != seqan::ArgumentParser::PARSE_OK)
        return res;

    getOptionValue(options.oPath, parser, "output");
    getOptionValue(options.sensitivity, parser, "sensitivity");
    getOptionValue(options.thread, parser, "thread");
    getOptionValue(options.index_t, parser, "index_type");
    getOptionValue(options.feature_t, parser, "feature_type");
    getOptionValue(options.gap_len, parser, "gap_len");
    getOptionValue(options.aln_flag, parser, "aln_flag");
    std::vector<std::string> args;
    args = getArgumentValues(parser, 0);
    if (length(args) < 2)
    {
        std::cerr << "[Err]::At least two argments required to specify the reads and genomes \n";
        return seqan::ArgumentParser::PARSE_ERROR;
    }
    else if (length(args) == 2)
    {
        //appendValue(options.rPath, args[0]);
        appendValue (options.r_paths, args[0]);
        appendValue (options.g_paths, args[1]);
    }
    else
    {
        Options::PathsType * pp = & options.r_paths;
        for (size_t i = 0; i < length(args); i++)
        {
            if (args[i] == CARTESIAN)
            {
                pp = & options.g_paths;
            }
            else
            {
                appendValue (*pp, args[i]);
            }
        }
    }
    getOptionValue(options.listN, parser, "listn1");
    getOptionValue(options.listN2, parser, "listn2");
    getOptionValue(options.alpha, parser, "alpha1");
    getOptionValue(options.alpha2, parser, "alpha2");
    getOptionValue(options.cordThr, parser, "cordThr");
    getOptionValue(options.senThr, parser, "senThr");
    getOptionValue(options.p1, parser, "p1");

    return seqan::ArgumentParser::PARSE_OK;

}