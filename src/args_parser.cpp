#include <seqan/arg_parse.h>
#include "base.h"

bool is_number(const std::string& s)
{
    return !s.empty() && std::find_if(s.begin(), 
        s.end(), [](char c) { return !std::isdigit(c); }) == s.end();
}

std::string CARTESIAN = "x";
seqan::ArgumentParser::ParseResult
parseCommandLineEmpty(Options & options, std::vector<const char*> new_args);
seqan::ArgumentParser::ParseResult 
parseCommandLineFilter(Options & options, std::vector<const char*> & new_args);
 
seqan::ArgumentParser::ParseResult
parseCommandLine(Options & options)
{
    //add submodules here
    int submodules = 0;
    //char const * str_g = "-g";
    //return 0;
    int argc = options.op_argc;
    char const ** & argv = options.op_argv;
    std::vector<const char*> new_args;
    char* new_vals;
    for (int i = 0; i < argc; i++)
    {
        if (i == 1)
        {
            if (strcmp(argv[1], "filter") == 0) 
            {
                submodules = 1;
                continue;
            }
            else 
            {
                submodules = 0;
            }   
        }
        new_args.push_back(argv[i]);
        if (std::string(argv[i]) == "-a" && (i + 1 >= argc || !is_number (argv[i + 1])))
        {
            new_vals = (char*)"1";
            new_args.push_back ((const char*)new_vals);
        }
        if (std::string(argv[i]) == "-g" && (i + 1 >= argc || !is_number (argv[i + 1])))
        {
            new_vals = (char*)"1";
            new_args.push_back ((const char*)new_vals);
        }
        if (std::string(argv[i]) == "-os" && (i + 1 >= argc || !is_number (argv[i + 1])))
        {
            new_vals = (char*)"1";
            new_args.push_back ((const char*)new_vals);
        }
        if (std::string(argv[i]) == "-oa" && (i + 1 >= argc || !is_number (argv[i + 1])))
        {
            new_vals = (char*)"1";
            new_args.push_back ((const char*)new_vals);
        }
        if (std::string(argv[i]) == "-r" && (i + 1 >= argc || !is_number (argv[i + 1])))
        {
            new_vals = (char*)"1";
            new_args.push_back ((const char*)new_vals);
        }
        if (std::string(argv[i]) == "-ss" && (i + 1 >= argc || !is_number (argv[i + 1])))
        {
            new_vals = (char*)"1";
            new_args.push_back ((const char*)new_vals);
        }
    }
    /*dout << "new_args len" << length(new_args) << "\n";*/
    if (length(new_args) < 3)
    {
        new_vals = (char*)"-h";
        new_args.push_back((const char*)new_vals);
    }
    else if (length(new_args) >= 3)
    {
        submodules = 1;
    }
    if (submodules == 0)
    {
        return parseCommandLineEmpty(options, new_args);
    }
    else if (submodules == 1)
    {   
        return parseCommandLineFilter(options, new_args);
    }
    return  seqan::ArgumentParser::PARSE_OK;
;
}

seqan::ArgumentParser::ParseResult
parseCommandLineEmpty(Options & options, std::vector<const char*> new_args)
{
    seqan::ArgumentParser parser("Linear");

    // Set short description, version, and date.
    setShortDescription(parser, "options and arguments. 🐼");
    setVersion(parser, options.version);
    setDate(parser, options.date);
    addUsageLine(parser, "<submodules> -h for help");
    seqan::addTextSection(parser,
        "Available submodules:");
    addListItem(parser,
                "filter: The submodule is to detect SVs signals hidden in long reads.",
                "It takes input as long reads and outputs SAM/BAM. Type \"linear filter -h\" for more info.");
    addArgument(parser, seqan::ArgParseArgument(
        seqan::ArgParseArgument::STRING, "empty",true)); 
    seqan::parse(parser, length(new_args), &new_args[0]);
    return seqan::ArgumentParser::PARSE_ERROR;
    //return res;
}


seqan::ArgumentParser::ParseResult parseCommandLineFilter(Options & options, 
        std::vector<const char*> & new_args)
{
    // Setup ArgumentParser.
    int argc = length(new_args); 
    seqan::ArgumentParser parser("Linear filter");

    // Set short description, version, and date.
    setShortDescription(parser, "options and arguments. 🐼");
    setVersion(parser, options.version);
    setDate(parser, options.date);
    addUsageLine(parser, "[\\fIOPTIONS\\fP] \\fIread.fa/fastq(.gz)\\fP \\fIgenome.fa(.gz)\\fP ");
    // Argument.
    addArgument(parser, seqan::ArgParseArgument(
        seqan::ArgParseArgument::INPUT_FILE, "read",true));
    //setHelpText(getArgument(parser, 1), "filter1 ");

    //Add descriptions 
    addDescription (parser, "The submodule filter is an ultra-fast SVs filter for detecting SVs in long reads. Unlike conventional assembly- or alignment-based pipelines, it uses generative models to detect SVs hidden in long reads. The filter does not compute assembly or alignment. However it takes as input long reads and outputs SAM/BAM*, which is compatible with alignment-based software including but not limited to samtools, SVs callers. The module is orders of magnitude faster than conventional pipelines for long-read SVs detection to compute and mainly aims for population-scale applications.");

    // Add Examples Section.
    addTextSection(parser, "Examples");
    addListItem(parser,
                "\\fPlinear filter \\fIreads_dir/example.fa.gz ref.fa\\fP",
                "Filter SVs for the reads example.fa using the reference genome ref.fa");                
    addListItem(parser,
                "\\fPlinear filter \\fIreads_dir/*.fa.gz x grch38/*.fa\\fP",
                "Set the additional argumnet \\fBx \\fPto filter SVs for more than one files of reads *.fa.gz using reference genomes grch38/*.fa");
    /*
    addArgument(parser, seqan::ArgParseArgument(
        seqan::ArgParseArgument::INPUT_FILE, "genome"));
    setHelpText(parser, 1, "Reference file .fa(.gz), .fasta(.gz)");
    */
    addSection(parser, "Basic options");
    addOption(parser, seqan::ArgParseOption(
        "o", "output", "Set the prefix of output. \
        The filter will use the prefix of the filename of reads as the prefix of output if the option isn't set",
            seqan::ArgParseArgument::STRING, "STR"));
    /*
    addOption(parser, seqan::ArgParseOption(
        "ot", "output_type", "Set the type of output. 1 to enable .apf; 2 to enable .sam; 4 to enable \
        standard bam; 8 to enable .bam for pbsv; Adding values, such as 1+2=3{DEFAULT} to enable .apf, .sam \
        \033[1;33mNote\033[m: The sv caller pbsv uses a non-standard format of sam/bam, \
        in which delimiter of header is space rather than tab. The non-standard sam/bam can be called \
        only by pbsv, while it's incompatible to tools following strictly to the sam/bam format, such as\
        samtools. Hence enable the option with value 8 only when the bam is supposed to be called by pbsv",
            seqan::ArgParseArgument::INTEGER, "INT"));
    */
    addOption(parser, seqan::ArgParseOption(
        "ot", "output_type", "Set the output format: 1 to enable .APF, an approximate map file for non-standard application; \
         2 to enable .SAM {DEFAULT}; 4 to enable .BAM; Example: Set -ot 3 (=1+2) to enable .apf and .sam",
            seqan::ArgParseArgument::INTEGER, "INT"));
    addOption(parser, seqan::ArgParseOption(
        "t", "thread", "Set the number of threads to run. -t 4 {DEFAULT}",
            seqan::ArgParseArgument::INTEGER, "INT"));
    addOption(parser, seqan::ArgParseOption(
        "g", "gap_len", "Set the minimal length of gaps. -g 50 {DEFAULT}. -g 0 to turn off mapping gaps.",
            seqan::ArgParseArgument::INTEGER, "INT"
        )); 
    /*
    addOption (parser, seqan::ArgParseOption(
        "os", "output_sam", "Set to Enable/Disbale the output in the format of SAM. -os or -s 1(Enable) {DEFAULT} ",
        seqan::ArgParseArgument::INTEGER, "INT"
        ));
    addOption (parser, seqan::ArgParseOption(
        "oa", "output_apf", "Set to Enable/Disbale the output in the format of APF. -oa or -oa 1(Enable) {DEFAULT} ",
        seqan::ArgParseArgument::INTEGER, "INT"
        ));
    */
    addOption (parser, seqan::ArgParseOption(
        "rg", "read_group", "Set the name of the read group specified in the SAM header",
        seqan::ArgParseArgument::STRING, "STR"
        ));
    addOption (parser, seqan::ArgParseOption(
        "sn", "sample_name", "Set the name of the sample specified in the SAM header",
        seqan::ArgParseArgument::STRING, "STR"
        ));
    addOption (parser, seqan::ArgParseOption(
        "ss", "sequence_sam", "Set to Enable/Disable printing sequence segment of reads in the SAM/BAM format. -ss 0(Disable) {DEFAULT}",
        seqan::ArgParseArgument::INTEGER, "INT"
        ));
    addSection(parser, "More options (tweak)");
    addOption (parser, seqan::ArgParseOption(
        "dup", "duplication", "Redetect duplications for signals of insertions. Enabling (-dup 1) this option will treat many insertions as duplications. This option is off (-dup 0) {DEFAULT}",
        seqan::ArgParseArgument::INTEGER, "INT"
        )); 
    addOption(parser, seqan::ArgParseOption(
        "b", "bal_flag", "Set to Enable/Disable dynamic balancing tasks schedule. -b 1(Enable) {DEFAULT}",
            seqan::ArgParseArgument::INTEGER, "INT" 
        ));   
    //addDefaultValue(parser, "gap_len", "1");

    //Advanced parms for mapping
    addOption(parser, seqan::ArgParseOption(
        "p", "preset", "Set predefined sets of parameters. -p 0 {DEFAULT} -p 1 Suggest to use for HiFi reads with cuteSV -p 2 Suggest to use for HiFi reads for cuteSV () or SVIM",
            seqan::ArgParseArgument::INTEGER, "INT"));
    /*
    addOption(parser, seqan::ArgParseOption(
        "a", "aln_flag", "Set to Enable/Disable alignment. -a 0(Disable) {DEFAULT}",
            seqan::ArgParseArgument::INTEGER, "INT"
        )); 
        */
    addOption(parser, seqan::ArgParseOption(
        "i", "index_type", "Choose the type of indices{1, 2}. -i 1 {DEFAULT}",
            seqan::ArgParseArgument::INTEGER, "INT"
        ));
    addOption(parser, seqan::ArgParseOption(
        "c", "apx_c_flag", "0 to turn off apx map",
            seqan::ArgParseArgument::INTEGER, "INT"
        )); 
    addOption(parser, seqan::ArgParseOption(
        "f", "feature_type", "Set types of features {1,2}. -f 2 (2-mer, 48bases){DEFAULT}",
            seqan::ArgParseArgument::INTEGER, "INT"
        )); 
    addOption (parser, seqan::ArgParseOption(
        "r", "reform_ccs_cigar_flag", "Enable/Disable compressing the cigar string for Pacbio CCS reads. -r 0(Disable) {DEFAULT}",
        seqan::ArgParseArgument::INTEGER, "INT"
        ));
    /*
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
        "p1", "par1", "mapping::p1",
            seqan::ArgParseArgument::INTEGER, "INT")); 
    */
    // Parse command line.
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, &new_args[0]);
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res;
    getOptionValue(options.oPath, parser, "output");
    getOptionValue(options.f_output_type, parser, "output_type");
    getOptionValue(options.sensitivity, parser, "preset");
    getOptionValue(options.thread, parser, "thread");
    getOptionValue(options.index_t, parser, "index_type");
    getOptionValue(options.feature_t, parser, "feature_type");
    getOptionValue(options.gap_len, parser, "gap_len");
    getOptionValue(options.apx_chain_flag, parser, "apx_c_flag");
    //getOptionValue(options.aln_flag, parser, "aln_flag");
    /*
    getOptionValue(options.sam_flag, parser, "output_sam");
    getOptionValue(options.apf_flag, parser, "output_apf");
    */
    getOptionValue(options.reform_ccs_cigar_flag, parser, "reform_ccs_cigar_flag");
    getOptionValue(options.read_group, parser, "read_group");
    getOptionValue(options.sample_name, parser, "sample_name");
    getOptionValue(options.sequence_sam_flag, parser, "sequence_sam");
    getOptionValue(options.bal_flag, parser, "bal_flag");
    getOptionValue(options.f_dup, parser, "duplication"); 

    //std::cerr << "xxxxx " << options.gap_len << isSet(parser, "gap_len") << "\n";
    String<std::string> args;
    args = getArgumentValues(parser, 0);
    if (length(args) < 2)
    {
        serr.print_message("\033[1;31mE[01]:\033[0m: Please specify the files of reads and genomes", 0, 1, std::cerr);
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
        bool f_CARTESIAN = false;
        for (size_t i = 0; i < length(args); i++)
        {
            if (args[i] == CARTESIAN)
            {
                pp = & options.g_paths;
                f_CARTESIAN = true;
            }
            else
            {
                appendValue (*pp, args[i]);
            }
        }
        if (!f_CARTESIAN)
        {
            options.op_status = 2;
            serr.print_message("\033[1;31mE[02]:\033[0mPlease add '\033[1;31mx\033[0m' between files of reads and genomes.", 0, 1, std::cerr);
            res = seqan::ArgumentParser::PARSE_ERROR;
            return res;
        }
    }
     
    /*
    getOptionValue(options.listN, parser, "listn1");
    getOptionValue(options.listN2, parser, "listn2");
    getOptionValue(options.alpha, parser, "alpha1");
    getOptionValue(options.alpha2, parser, "alpha2");
    getOptionValue(options.cordThr, parser, "cordThr");
    getOptionValue(options.senThr, parser, "senThr");
    getOptionValue(options.p1, parser, "p1");
    */
    /*
    for (int i = 0; i < length(options.r_paths); i++)
    {
        dout << "rpaths" << options.r_paths[i] << "\n";
    }
    for (int i = 0; i < length(options.g_paths); i++)
    {
        dout << "gpaths" << options.g_paths[i] << "\n";
    }
    */
    return seqan::ArgumentParser::PARSE_OK;

}
