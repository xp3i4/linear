#include <seqan/arg_parse.h>
#include <iostream>
#include <fstream>
#include <ctime>
#include "mapper.h"
#include "pmpfinder.h"
#include "chain_map.h"
#include "gap.h"
#include "align_interface.h"

using namespace seqan; 

typedef StringSet<String<uint64_t> > CordsType;

seqan::ArgumentParser::ParseResult
parseCommandLine(Options & options, int argc, char const ** argv)
{
    // Setup ArgumentParser.
    seqan::ArgumentParser parser("Linear");
    // Set short description, version, and date.
    setShortDescription(parser, "Alignment of SMRT sequencing read");
    setVersion(parser, "1.0");
    setDate(parser, "May 2018");

    // Define usage line and long description.
    addUsageLine(parser,
                    "[\\fIOPTIONS\\fP] \"\\fIread.fa\\fP\" \"\\fIgnome.fa\\fP\"");
    addDescription(parser,
                    "Program for mapping raw SMRT sequencing reads to reference genome.");

    // Argument.
    addArgument(parser, seqan::ArgParseArgument(
        seqan::ArgParseArgument::INPUT_FILE, "read"));
    setHelpText(parser, 0, "Reads file .fa, .fasta");

    addArgument(parser, seqan::ArgParseArgument(
        seqan::ArgParseArgument::INPUT_FILE, "genome", true));
    setHelpText(parser, 1, "Reference file .fa, .fasta");

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
    
    getOptionValue(options.listN, parser, "listn1");
    getOptionValue(options.listN2, parser, "listn2");
    getOptionValue(options.alpha, parser, "alpha1");
    getOptionValue(options.alpha2, parser, "alpha2");
    getOptionValue(options.cordThr, parser, "cordThr");
    getOptionValue(options.senThr, parser, "senThr");
    getOptionValue(options.p1, parser, "p1");
    seqan::getArgumentValue(options.rPath, parser, 0);
    seqan::getArgumentValue(options.gPath, parser, 1);

    return seqan::ArgumentParser::PARSE_OK;

}

Mapper::Mapper(Options & options):
    record(options),
    qIndex(genomes()), 
    of(toCString(options.getOutputPath()))
{
    outputPrefix = getFileName(getFileName(options.getReadPath()), '.', 0);
    switch (options.sensitivity)
    {
        case 0: 
        {
            parm = parm0; //normal
            break;
        }
        case 1:
        {
            parm =  parm1; //fast
            break;
        }
        case 2:
        {
            parm = parm2; //sensitive
            break;
        }
    }
    _thread = options.thread;
}

int Mapper::createIndex(bool efficient)
{
    std::cerr << ">>Create index \r";
    createHIndex(genomes(), qIndex, _thread, efficient);
    return 0;
}

int print_align_sam_header_ (StringSet<CharString> & genomesId,
                             StringSet<String<Dna5> > & genomes,
                             std::ofstream & of
                            )
{
    of << "@HD\tVN:1.6\n";
    for (int k = 0; k < length(genomesId); k++)
    {
        of << "@SQ\tSN:" << genomesId[k] << "\tLN:" << length(genomes[k]) << "\n";
    }
    of << "@PG\tPN:" << "Linear\n";
}
int print_align_sam_record_(StringSet<String< BamAlignmentRecord > > & records, 
                     StringSet<String<uint64_t> > & cordSet,
                     StringSet<CharString> & readsId, 
                     StringSet<CharString> & genomesId,
                     std::ofstream & of
                    )
{
    for (int i = 0; i < length(records); i++)
    {
        for (int j = 0; j < length(records[i]); j++)
        {
            records[i][j].qName = readsId[i];
            CharString g_id = genomesId[records[i][j].rID];
            writeSam(of, records[i][j], g_id);
        }
    }
}
int print_align_sam_record_(StringSet<String< BamAlignmentRecordLink> > & records, 
                     StringSet<String<uint64_t> > & cordSet,
                     StringSet<CharString> & readsId, 
                     StringSet<CharString> & genomesId,
                     std::ofstream & of
                    )
{
    for (int i = 0; i < length(records); i++)
    {
        for (int j = 0; j < length(records[i]); j++)
        {
            records[i][j].qName = readsId[i];
            CharString g_id = genomesId[records[i][j].rID];
            int dt = writeSam(of, records[i], j, g_id);
            //<<<debug 
        }

        for (int j = 0; j < length(records[i]); j++)
        {
            std::pair<int,int> lens;
            lens = countCigar(records[i][j].cigar);
            if (records[i][j].isEnd())
                std::cout << "\n";
            std::cout << "cigarLen " << lens.first << " " << lens.second << "\n";
            //>>>debug
        }
    }
}
int print_align_sam (Mapper & mapper)
{
    std::string filePath = mapper.getOutputPrefix() + ".sam";
    mapper.getOf().open(toCString(filePath));
    print_align_sam_header_(mapper.genomesId(), 
                            mapper.genomes(),
                            mapper.getOf()
                           );
    print_align_sam_record_(mapper.getBamRecords(),
                            mapper.cords(),
                            mapper.readsId(),
                            mapper.genomesId(),
                            mapper.getOf()
                           ); 
    mapper.getOf().close();
}

void Mapper::printCordsRaw()
{
    double time = sysTime();
    //unsigned strand;
    unsigned cordCount = 0;

    cordCount = 0;
    for (unsigned k = 0; k < length(cordSet); k++)
    {
        //unsigned recordCount = 0;
        if (!empty(cordSet[k]))
        {
            for (unsigned j = 1; j < length(cordSet[k]); j++)
            {
                if (_DefaultHit.isBlockEnd(cordSet[k][j-1]) )//&& ++recordCount < 10)
                {
                    of <<"@S1_"<< k+1 << " " << length(reads()[k]) << " "
                    << _DefaultCord.getCordY(cordSet[k][j]) << " " << length(cordSet[k]) << " x " 
                    << _getSA_i1(_DefaultCord.getCordX(cordSet[k][j])) << " " << cordCount << " "
                    << _getSA_i2(_DefaultCord.getCordX(cordSet[k][j]))  << " " 
                    << "\n";   
                    cordCount = 0;
                }
                cordCount++;
            }
        }
    }
    std::cerr << ">Write results to disk          " << std::endl;
    std::cerr << "    End writing results. Time[s]" << sysTime() - time << std::endl; 
}
/**
 * print all cords with cordinates
 */
void print_cords_txt(Mapper & mapper)
{
    std::cerr << ">>Write results to disk        \r";
    double time = sysTime();
    unsigned cordCount = 0;
    uint64_t readCordEnd;
    uint64_t seqsCordEnd;
    char main_icon_strand = '+', icon_strand = '+';
    int fflag = 0;
    CordsType & cords = mapper.cords();
    for (unsigned k = 0; k < length(cords); k++)
    {
        if (!empty(cords[k]))
        {
            for (unsigned j = 1; j < length(cords[k]); j++)
            {
                if (_DefaultHit.isBlockEnd(cords[k][j-1]))
                {
                    unsigned m = j; 
                    int main_strand_count = 0;
                    int block_len = 0;
                    ///>determine the main strand
                    while (!_DefaultHit.isBlockEnd(cords[k][m]))
                    {
                        if (_DefaultCord.getCordStrand(cords[k][m]))
                        {
                            main_strand_count++;
                        }
                        block_len++;
                        m++;
                    }
                    if (main_strand_count > (block_len >> 1))
                    {
                        main_icon_strand = '-';
                    }
                    else
                    {
                        main_icon_strand = '+';
                    }
                    ///>print the header
                    for (unsigned i = j; ; i++)
                    {
                        if (_DefaultHit.isBlockEnd(cords[k][i]) || i == length(cords[k]) - 1)
                        {
                            readCordEnd = _DefaultCord.getCordY(cords[k][i]) + window_size;
                            seqsCordEnd = _getSA_i2(_DefaultCord.getCordX(cords[k][i])) + window_size;
                            break;
                        }
                    }
                    if (k > 0)
                    {
                        mapper.getOf() << "\n";
                    } 
                    mapper.getOf() << "@> "
                        << mapper.readsId()[k] << " " 
                        << mapper.readLens()[k] << " "
                        << _DefaultCord.getCordY(cords[k][j]) << " " 
                        << std::min(readCordEnd, (uint64_t)mapper.readLens()[k]) << " " 
                        << main_icon_strand<< " "
                        << mapper.genomesId()[get_cord_id(cords[k][j])] << " " 
                        << length(mapper.genomes()[get_cord_id(cords[k][j])]) << " "
                        << get_cord_x(cords[k][j]) << " " 
                        << seqsCordEnd << "\n";
                    cordCount = 0;
                    fflag = 1;
                }
                ///>print the coordinates
                icon_strand = (get_cord_strand(cords[k][j]))?'-':'+';
                CharString mark = "| ";
                if (icon_strand != main_icon_strand)
                {
                    mark = (icon_strand == '+') ? "|**+++++++++++ " :"|**----------- ";
                }
                int64_t d1 = 0;//_DefaultCord.getCordY(cords[k][1]);
                int64_t d2 = 0;
                if (!fflag)
                {
                    d1 = int64_t(get_cord_x(cords[k][j]) - get_cord_x(cords[k][j - 1]));
                    d2 = int64_t(get_cord_y(cords[k][j]) - get_cord_y(cords[k][j - 1]));
                }
                mapper.getOf() << mark  
                               << _DefaultCord.getCordY(cords[k][j]) << " " 
                               << get_cord_x(cords[k][j]) << " " 
                               << d2 << " " 
                               << d1 << " " 
                               << j << " \n";

                cordCount++;
                fflag = 0;
            }
        }

    }
    close(mapper.getOf());
    std::cerr << "--Write results to disk       100% Elapsed Time[s] " 
              << sysTime() - time << "\n";
}

int print_clips_gff_(StringSet<String<uint64_t> > & clips, 
              StringSet<CharString> & readsId, 
              std::ofstream & of, 
              std::string & outputPrefix)
{
    std::string file_path = outputPrefix + ".gff";
    std::cerr << "[]::filepath " << file_path << "\n";
    of.open(toCString(file_path));
    for (unsigned i = 0; i < length(clips); i++)
    {
        if (!empty(clips[i]))
        {
            of << i << " " << readsId[i] << " ";
            for (unsigned j = 0; j < length(clips[i]); j++)
            {
                uint64_t cord_x = _DefaultCord.getCordX(clips[i][j]);
                uint64_t cord_y = _DefaultCord.getCordY(clips[i][j]);
                of << cord_x << " ";   
            }
            of << '\n';
        }
    }
    of.close();
    return 0;
}

int print_clips_gvf_(StringSet<String<uint64_t> > & clips, 
              StringSet<CharString> & readsId, 
              StringSet<CharString> & genomesId,
              std::ofstream & of, 
              std::string outputPrefix)
{
    std::string file_path = outputPrefix + ".gvf";
    //std::cerr << "[]::filepath " << file_path << "\n";
    of.open(toCString(file_path));
    of << "##gvf-version 1.10\n ";
    std::string source = ".";
    std::string type = ".";
    for (unsigned i = 0; i < length(clips); i++)
    {
        if (!empty(clips[i]))
        {
            for (unsigned j = 0; j < length(clips[i]); j++)
            {
                uint64_t cord_x = _getSA_i2(_DefaultCord.getCordX(clips[i][j]));
                //uint64_t cord_y = _DefaultCord.getCordY(clips[i][j]);
                CharString genomeId = genomesId[get_cord_id(clips[i][j])];
                if ((j >> 1) << 1 == j)
                {
                    of  << genomeId << "\t" 
                        << source << "\t" 
                        << type << "\t" 
                        << cord_x << "\t";   
                    if (j == length(clips[i]) - 1)
                    {
                        of << " . readId=" << readsId[i] << ";" << i << "\n";
                    }
                }
                else
                {
                    of << cord_x << " readId=" << readsId[i] << ";" << i << "\n";
                }
            }
        }
    }
    of.close();
    return 0;
}
int print_clips_gff(Mapper & mapper)
{
    print_clips_gff_(mapper.getClips(), mapper.readsId(), mapper.getOf(), mapper.getOutputPrefix());
    return 0;
}

int print_clips_gvf(Mapper & mapper)
{
    print_clips_gvf_(mapper.getClips(), mapper.readsId(), mapper.genomesId(), mapper.getOf(), mapper.getOutputPrefix());
    return 0;
}

int rawMap_dst2_MF(LIndex & index,
                   StringSet<String<short> > & f2,
                   StringSet<String<Dna5> > & reads,
                   MapParm & mapParm,
                   StringSet<String<uint64_t> > & cords,
                   StringSet<String<uint64_t> > & clips,
                   StringSet<String<Dna5> > & seqs,
                   StringSet<String<BamAlignmentRecordLink> >& bam_records,
                   unsigned & threads,
                   int p1
                  )
{
  
    typedef String<Dna5> Seq;
    //double time=sysTime();
    float senThr = mapParm.senThr / window_size;
    float cordThr = mapParm.cordThr / window_size;
    MapParm complexParm = mapParm;
    complexParm.alpha = complexParm.alpha2;
    complexParm.listN = complexParm.listN2;
    String<uint64_t> gap_len;
    String<uint64_t> red_len;
    resize (gap_len, threads);
    resize (red_len, threads);
    for (int i = 0; i < length(gap_len); i++)
{
    gap_len[i]  = 0;
    red_len[i] = 0;
}
    //double time2 = sysTime();
int su = 0;
int64_t len = 0;
#pragma omp parallel
{
    unsigned size2 = length(reads) / threads;
    unsigned ChunkSize = size2;
    Seq comStr;
    //Anchors anchors(Const_::_LLTMax, AnchorBase::size);
    Anchors anchors;
    String<uint64_t> crhit;
    StringSet<String<uint64_t> >  cordsTmp;
    StringSet< String<short> > f1;
    StringSet<String<uint64_t> > clipsTmp;
    StringSet<String<BamAlignmentRecordLink> > bam_records_tmp;
    unsigned thd_id =  omp_get_thread_num();
    if (thd_id < length(reads) - size2 * threads)
    {
        ChunkSize = size2 + 1;
    }
    resize(cordsTmp, ChunkSize);
    resize(clipsTmp, ChunkSize);
    resize(bam_records_tmp, ChunkSize);
    resize(f1, 2);
    unsigned c = 0;
    
    String<uint64_t> g_hs;
    String<uint64_t> g_anchor;
    String<uint64_t> bands;
    resize (g_hs, 1ULL << 20);
    resize (g_anchor, 1ULL<<20);
    #pragma omp for
    for (unsigned j = 0; j < length(reads); j++)
    {
        if (length(reads[j]) >= mapParm.minReadLen)
        {
            red_len[thd_id] += length(reads[j]);
            std::cout << "[]::rawmap::j " << j <<"\n";
            float cordLenThr = length(reads[j]) * cordThr;
            _compltRvseStr(reads[j], comStr);
            createFeatures(begin(reads[j]), end(reads[j]), f1[0]);
            createFeatures(begin(comStr), end(comStr), f1[1]);
            anchors.init(1);
            clear(crhit);
            mnMapReadList(index, reads[j], anchors, mapParm, crhit);
            path_dst(begin(crhit), end(crhit), f1, f2, cordsTmp[c], cordLenThr);
            if (_DefaultCord.getMaxLen(cordsTmp[c]) < length(reads[j]) * senThr)
            {
                clear(cordsTmp[c]);
                anchors.init(1);
                clear(crhit);
                mnMapReadList(index, reads[j], anchors, complexParm, crhit);
                path_dst(begin(crhit), end(crhit), f1, f2, cordsTmp[c], cordLenThr);
            }   
            gap_len[thd_id] += mapGaps(seqs, reads[j], comStr, cordsTmp[c], g_hs, g_anchor, clipsTmp[c], f1, f2, p1, 192);
            //align_cords(seqs, reads[j], comStr, cordsTmp[c], bam_records_tmp[c], p1);
        }   
        c += 1;
    } 
    #pragma omp for ordered
    for (unsigned j = 0; j < threads; j++)
        #pragma omp ordered
        {
            append(cords, cordsTmp);
            append(clips, clipsTmp);
            append(bam_records, bam_records_tmp);
        }
    
}
    //std::cerr << "    End raw mapping. Time[s]: " << sysTime() - time << std::flush << std::endl;
    return 0;
}

/*
 *[]::map
 */
int map(Mapper & mapper, int p1)
{
    //printStatus();
    StringSet<String<short> > f2;
    mapper.createIndex(false); // true: destruct genomes string to reduce memory footprint
    createFeatures(mapper.genomes(), f2, mapper.thread());
    SeqFileIn rFile(toCString(mapper.readPath()));
    unsigned k = 1, j = 0;
    unsigned blockSize = 50000;
    StringSet<String<char> > dotstatus;
    resize(dotstatus, 3);
    dotstatus[0] = ".   ";
    dotstatus[1] = "..  ";
    dotstatus[2] = "... ";
    while (!atEnd(rFile))
    {
        double time1 = sysTime();
        clear (mapper.reads());
        std::cerr <<  ">>Map::file_I/O  block " << k << dotstatus[j++ % length(dotstatus)] << "\r";
        readRecords_block(mapper.readsId(), mapper.reads(), mapper.readLens(), rFile, blockSize);
        std::cerr << "                                    \r";
        std::cerr <<  ">>Map::mapping  block "<< k << " Size " << length(mapper.reads()) << " " << dotstatus[j++ % length(dotstatus)] << "\r";
        time1 = sysTime() - time1;
        double time2 = sysTime();
        rawMap_dst2_MF(mapper.index(), 
                       f2, 
                       mapper.reads(), 
                       mapper.mapParm(), 
                       mapper.cords(), 
                       mapper.getClips(),
                       mapper.genomes(),
                       mapper.getBamRecords(),
                       mapper.thread(), 
                       p1);
        time2 = sysTime() - time2;
        std::cerr <<  "--Map::file_I/O+Map block "<< k << " Size " << length(mapper.reads()) << " Elapsed Time[s]: file_I/O " << time1 << " map "<< time2 << "\n";
        k++;
    }
    mapper.index().clear(); 
    print_cords_txt(mapper);
    print_align_sam(mapper);
    print_clips_gvf(mapper);
    return 0;
}
int main(int argc, char const ** argv)
{
    double time = sysTime();
    std::cerr << "[]\n";
    (void)argc;
    Options options;
    seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;
    std::cerr << "Encapsulated version: Mapping reads efficiently" << std::endl;
    Mapper mapper(options);
    omp_set_num_threads(mapper.thread());
    map(mapper, options.p1);

    //mapper.print_vcf();
    std::cerr << "  Result Files: \033[1;31m" << options.oPath << "\033[0m" << std::endl;
    std::cerr << "                \033[1;31m" << (mapper.getOutputPrefix() + ".gvf") << "\033[0m" << std::endl;
    std::cerr << "                \033[1;31m" << (mapper.getOutputPrefix() + ".sam") << "\033[0m" << std::endl;
    std::cerr << "Time in sum[s] " << sysTime() - time << std::endl;
    return 0;
}
