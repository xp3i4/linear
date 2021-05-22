#include <seqan/arg_parse.h>
#include "base.h"

bool is_number(const std::string& s)
{
    return !s.empty() && std::find_if(s.begin(), 
        s.end(), [](char c) { return !std::isdigit(c); }) == s.end();
}

std::string CARTESIAN = "x";
seqan::ArgumentParser::ParseResult
parseCommandLine(Options & options, int argc, char const ** argv)
{
    char const * str_g = "-g";
    //return 0;
    std::vector<const char*> new_args;
    char* new_vals;
    for (int i = 0; i < argc; i++)
    {
        new_args.push_back(argv[i]);
        if (std::string(argv[i]) == "-a" && (i + 1 >= argc || !is_number (argv[i + 1])))
        {
            new_vals = "1";
            new_args.push_back ((const char*)new_vals);
        }
        if (std::string(argv[i]) == "-g" && (i + 1 >= argc || !is_number (argv[i + 1])))
        {
            new_vals = "1";
            new_args.push_back ((const char*)new_vals);
        }
        if (std::string(argv[i]) == "-s" && (i + 1 >= argc || !is_number (argv[i + 1])))
        {
            new_vals = "1";
            new_args.push_back ((const char*)new_vals);
        }
        if (std::string(argv[i]) == "-r" && (i + 1 >= argc || !is_number (argv[i + 1])))
        {
            new_vals = "1";
            new_args.push_back ((const char*)new_vals);
        }
        if (std::string(argv[i]) == "-ss" && (i + 1 >= argc || !is_number (argv[i + 1])))
        {
            new_vals = "1";
            new_args.push_back ((const char*)new_vals);
        }
    }
    argc = length(new_args); 
    // Setup ArgumentParser.
    seqan::ArgumentParser parser("linear");
    // Set short description, version, and date.
    setShortDescription(parser, "Options & arguments ");
    setVersion(parser, options.version);
    setDate(parser, options.date);
    addUsageLine(parser, "[\\fIOPTIONS\\fP] \\fIread.fa/fastq(.gz)\\fP \\fIgenome.fa(.gz)\\fP ");
    //addDescription(parser,
                    //"Extensible Framework of extensifor noisy reads");
    // Argument.
    addArgument(parser, seqan::ArgParseArgument(
        seqan::ArgParseArgument::INPUT_FILE, "read",true));
    setHelpText(parser, 0, "Reads file .fa(.gz), .fasta(.gz), .fq(.gz), .fastq(.gz) ");

    // Add Examples Section.
    //addTextSection(parser, "Examples");
    addTextSection(parser, "Examples");
    addListItem(parser,
                "\\fPlinear \\fIreads_dir/*.fa.gz grch37/chr1.fa\\fP",
                "Map the set of reads to the reference chr1.fa");
    addListItem(parser,
                "\\fPlinear \\fIreads_dir/*.fa.gz grch37/chr1.fa -g 0\\fP",
                "Map the set of reads to the reference chr1.fa with the mapping of gaps disabled.\
                In such case, Linear generates approximate range of the reference where the reads are\
                 supposed to be mapped to");
    addListItem(parser,
                "\\fPlinear \\fIreads_dir/*.fa.gz x grch37/*.fa\\fP",
                "Use the argumnet \\fBx \\fPto map the set of reads to the set of genomes");
    addListItem(parser,
                "\\fPlinear \\fIreads.fa genome.fa -g -a\\fP",
                "Use the option \\fB-a \\fPto enable the alignment"
        );

/*
    addArgument(parser, seqan::ArgParseArgument(
        seqan::ArgParseArgument::INPUT_FILE, "genome"));
    setHelpText(parser, 1, "Reference file .fa(.gz), .fasta(.gz)");
*/
    addSection(parser, "Basic Options");
    addOption(parser, seqan::ArgParseOption(
        "o", "output", "Set the path of output.",
            seqan::ArgParseArgument::STRING, "STR"));
    addOption(parser, seqan::ArgParseOption(
        "p", "preset", "Set preset of parms. -s 0 normal {DEFAULT} -s 1 efficient  -s 2 additional",
            seqan::ArgParseArgument::INTEGER, "INT"));
    addOption(parser, seqan::ArgParseOption(
        "t", "thread", "Set threads to run -t 4 {DEFAULT}",
            seqan::ArgParseArgument::INTEGER, "INT"));
    addOption(parser, seqan::ArgParseOption(
        "g", "gap_len", "Minimal length of gaps for mapping. -g 50 {DEFAULT}. -g 0 to turn off mapping of gaps.",
            seqan::ArgParseArgument::INTEGER, "INT"
        )); 
    addOption(parser, seqan::ArgParseOption(
        "a", "aln_flag", "0 to turn off alignment module",
            seqan::ArgParseArgument::INTEGER, "INT"
        )); 
    addOption (parser, seqan::ArgParseOption(
        "s", "sam_flag", "Set to Enable/Disbale the output in the format of SAM/BAM. -s or -s 1(Enable) {DEFAULT} ",
        seqan::ArgParseArgument::INTEGER, "INT"
        ));
    addOption (parser, seqan::ArgParseOption(
        "rg", "read_group", "Set the tag of name of read group specified in the SAM header",
        seqan::ArgParseArgument::STRING, "STR"
        ));
    addOption (parser, seqan::ArgParseOption(
        "sn", "sample_name", "Set the name of sample specified in the SAM header",
        seqan::ArgParseArgument::STRING, "STR"
        ));
    addOption (parser, seqan::ArgParseOption(
        "ss", "sequence_sam", "Set to Enable/Disable printing sequence segment of reads in the SAM/BAM format. -ss 0(Disabled) {DEFAULT}",
        seqan::ArgParseArgument::INTEGER, "INT"
        ));
    addOption(parser, seqan::ArgParseOption(
        "b", "bal_flag", "Set to Enable/Disable dynamic balancing tasks schedule. -b 1(Enable) {DEFAULT}",
            seqan::ArgParseArgument::INTEGER, "INT" 
        ));   
//    addDefaultValue(parser, "gap_len", "1");

// Advanced parms for mapping
    addSection(parser, "Advanced optoins(for tweak & debug)");
    addOption(parser, seqan::ArgParseOption(
        "i", "index_type", "Choose the type of indices{1, 2}. -i 1 {DEFAULT}",
            seqan::ArgParseArgument::INTEGER, "INT"
        ));
    addOption(parser, seqan::ArgParseOption(
        "c", "apx_chain_flag", "0 to turn off chaining in apx mapping",
            seqan::ArgParseArgument::INTEGER, "INT"
        )); 
    addOption(parser, seqan::ArgParseOption(
        "f", "feature_type", "Set types of features {1,2}. -f 2 (2-mer, 48bases){DEFAULT}",
            seqan::ArgParseArgument::INTEGER, "INT"
        )); 
    addOption (parser, seqan::ArgParseOption(
        "r", "reform_ccs_cigar_flag", "Enable/Disable compressing the cigar string for Pacbio CCS reads. -r 0{Disale} {DEFAULT}",
        seqan::ArgParseArgument::INTEGER, "INT"
        ));
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



    // Parse command line.
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, &new_args[0]);

    if (res != seqan::ArgumentParser::PARSE_OK)
        return res;
    getOptionValue(options.oPath, parser, "output");
    getOptionValue(options.sensitivity, parser, "preset");
    getOptionValue(options.thread, parser, "thread");
    getOptionValue(options.index_t, parser, "index_type");
    getOptionValue(options.feature_t, parser, "feature_type");
    getOptionValue(options.gap_len, parser, "gap_len");
    getOptionValue(options.apx_chain_flag, parser, "apx_chain_flag");
    getOptionValue(options.aln_flag, parser, "aln_flag");
    getOptionValue(options.sam_flag, parser, "sam_flag");
    getOptionValue(options.reform_ccs_cigar_flag, parser, "reform_ccs_cigar_flag");
    getOptionValue(options.read_group, parser, "read_group");
    getOptionValue(options.sample_name, parser, "sample_name");
    getOptionValue(options.sequence_sam_flag, parser, "sequence_sam");
    getOptionValue(options.bal_flag, parser, "bal_flag");

    //std::cerr << "xxxxx " << options.gap_len << isSet(parser, "gap_len") << "\n";
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