// ==========================================================================
//                           Mapping SMRT reads 
// ==========================================================================
// Copyright (c) 2006-2016, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: cxpan <chenxu.pan@fu-berlin.de>
// ==========================================================================

#include <csignal>
#include "mapper.h"

using namespace seqan; 

template <typename TDna, typename TSpec>
Mapper<TDna, TSpec>::Mapper(Options & options):
    record(options),
    parm(options),
    qIndex(genomes()),
    of(toCString(options.getOutputPath()))
{
        std::cerr << "parm init\n " << options.sensitivity << " \n" ;
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
        // parmt for test 
        parm = parmt;
}

template <typename TDna, typename TSpec>
int Mapper<TDna, TSpec>::createIndex()
{
    std::cerr << "Creating index \n";
    _createQGramIndex(qIndex, genomes());
    return 0;
}

template <typename TDna, typename TSpec>
Pair<typename Mapper<TDna,TSpec>::HitType, typename Mapper<TDna, TSpec>::HitType>
Mapper<TDna, TSpec>::getAliX(unsigned k) const
{
    typedef typename Mapper<TDna,TSpec>::HitType HitType;
    typedef Pair<HitType, HitType> Pair;
    if (empty(res.hits[k]))
        return Pair::Pair(~(HitType)0, ~(HitType)0);
    else
        return Pair::Pair(_getSA_i1(res.hits[k][0]), _getSA_i2(res.hits[k][0])) << std::endl;
}

template <typename TDna, typename TSpec>
typename Mapper<TDna,TSpec>::HitType 
Mapper<TDna, TSpec>::getHitX(typename Mapper<TDna, TSpec>::HitType const & hit) const
{
    return _DefaultCord.getCordX(_DefaultCord.hit2Cord(hit));
}

template <typename TDna, typename TSpec>
typename Mapper<TDna,TSpec>::HitType 
Mapper<TDna, TSpec>::getHitY(typename Mapper<TDna, TSpec>::HitType const & hit) const
{
    return _DefaultCord.getCordY(_DefaultCord.hit2Cord(hit));
}

template <typename TDna, typename TSpec>
typename Mapper<TDna,TSpec>::CordType 
Mapper<TDna, TSpec>::getCordX(typename Mapper<TDna, TSpec>::CordType const & cord) const
{
    return _DefaultCord.getCordX(cord);
}

template <typename TDna, typename TSpec>
typename Mapper<TDna,TSpec>::CordType 
Mapper<TDna, TSpec>::getCordY(typename Mapper<TDna, TSpec>::CordType const & cord) const
{
    return _DefaultCord.getCordY(cord);
}

template <typename TDna, typename TSpec>
void Mapper<TDna, TSpec>::printHits()
{
    unsigned j = 0;
    std::cout << "printHits: " << lengthSum(res.hits) << " in sum " << std::endl;
    for (auto && hitStr : res.hits)
    {
        std::cout << j++ << "th hits"<< std::endl;
        for (auto && hit : hitStr)
            std::cout << _getSA_i1(getHitX(hit)) << " " << _getSA_i2(getHitX(hit)) << " " << getHitY(hit) << std::endl;
        std::cout << std::endl;
    }
}

template <typename TDna, typename TSpec>
void Mapper<TDna, TSpec>::printBestHitsStart()
{
    unsigned k = 0;
    for (auto && hitStr : res.hits)
        if (empty(hitStr))
            std::cout << std::endl;
        else
            std::cout << k++ << " " << getHitX(hitStr[0]) - getHitY(hitStr[0]) << std::endl;
}

template <typename TDna, typename TSpec>
void Mapper<TDna, TSpec>::printResult()
{}

template <typename TDna, typename TSpec>
void Mapper<TDna, TSpec>::printParm()
{
    parm.print();
}

template <typename TDna, typename TSpec>
void Mapper<TDna, TSpec>::printCords(std::ostream & of)
{
    double time = sysTime();
    unsigned strand;
    for (unsigned k = 0; k < length(cordSet); k++)
    {
        if (empty(cordSet))
            of << k << " th Strand " << strand << " 2 " << length(reads()[k]) << "\n";
        else
        {
            if (_DefaultCord.getCordStrand(back(cordSet[k]))) 
                strand = 1;
            else 
                strand = 0;
            of << k << " th Strand " << strand << " " << length(reads()[k]) << "\n";
            _DefaultCord.print(cordSet[k],of);
        }
    }
    std::cerr << ">Write results to disk          " << std::endl;
    std::cerr << "    End writing results. Time[s]" << sysTime() - time << std::endl;
}

template <typename TDna, typename TSpec>
void Mapper<TDna, TSpec>::printCordsAll()
{
    double time = sysTime();
    unsigned strand;
    for (unsigned k = 0; k < length(cordSet); k++)
    {
        if (empty(cordSet[k]))
            of << k << " th Strand " << " 2 length " << length(reads()[k]) << "\nlength of cords -\n\n";
        else
        {
            if (_DefaultCord.getCordStrand(back(cordSet[k]))) 
                strand = 1;
            else 
                strand = 0;
            of << k << " th Strand " << strand << " length " 
            << length(reads()[k]) << "\nlength of cords " << "\n";
            unsigned cordCount = 0;
            unsigned first = 0;
            unsigned cover = 0;
            for (unsigned j = 1; j < length(cordSet[k]); j++)
            {
                of << _DefaultCord.getCordY(cordSet[k][j]) << " " 
                    << _getSA_i1(_DefaultCord.getCordX(cordSet[k][j])) << " "
                    << _getSA_i2(_DefaultCord.getCordX(cordSet[k][j]))  << std::endl;
                if (_DefaultCord.getCordY(cordSet[k][j]) - first < 192)
                    cover += _DefaultCord.getCordY(cordSet[k][j]) - first;
                else
                    cover += 192;
                first = _DefaultCord.getCordY(cordSet[k][j]);
                cordCount++;
                if (_DefaultHit.isBlockEnd(cordSet[k][j]))
                {
                   
                    of << "coverage " << (float)cover / (length(reads()[k])) << "\n";
                    if (j < length(cordSet[k]) - 1)
                    {
                        of << "\n" << k << " th Strand " << strand << " length " 
                        << length(reads()[k]) << "\nlength of cords " << "\n";
                    }   
                     cordCount =0;
                    first = 0;
                    cover = 0;
                }
                
 
            }
            of << "\n";
            //_DefaultCord.print(cordSet[k],of);
        }
    }
    std::cerr << ">Write results to disk          " << std::endl;
    std::cerr << "    End writing results. Time[s]" << sysTime() - time << std::endl;
}

template <typename TDna, typename TSpec>
void Mapper<TDna, TSpec>::printCords()
{
    std::cerr << "Writing results to disk \r";
    double time = sysTime();
    unsigned strand;
    unsigned count = 0;
    for (unsigned k = 0; k < length(cordSet); k++)
    {
        count++;
        if (empty(cordSet[k]))
            of << k << " th Strand " << "2 length " << length(reads()[k]) << "\n";
        else
        {
            if (_DefaultCord.getCordStrand(back(cordSet[k]))) 
                strand = 1;
            else 
                strand = 0;
            of << k << " th Strand " << strand << " length " << length(reads()[k]) << "\n";
            _DefaultCord.print(cordSet[k], of);
        }
    }
    of.close();
    std::cerr << ">Write results to disk        " << count << std::endl;
    std::cerr << "    End writing results. Time[s]" << sysTime() - time << std::endl;
}

template <typename TDna, typename TSpec>
unsigned Mapper<TDna, TSpec>::sens()
{
    return parm.sensitivity;
}


template <typename TDna, typename TSpec>
void map(Mapper<TDna, TSpec> & mapper)
{
    //printStatus();
    double time = sysTime();
    mapper.createIndex();
    resize(mapper.hits(), length(mapper.reads()));
    resize(mapper.cords(), length(mapper.reads()));
    //rawMap<TDna, TSpec>(mapper.index(), mapper.reads(), mapper.genomes(),
    //                     mapper.mapParm(), mapper.hits(), mapper.cords());
    rawMapAllComplex<TDna, TSpec>(mapper.index(), mapper.reads(), mapper.genomes(),
                         mapper.mapParm(), mapper.hits(), mapper.cords());
        //mapper.printHits();
    mapper.printCordsAll();
    std::cerr << length(mapper.cords()) << " " << length(mapper.reads()) << " \n";
    std::cerr << "Time in sum[s] " << sysTime() - time << std::endl;
}


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
    getOptionValue(options.sensitivity, parser, "sensitivity");

    seqan::getArgumentValue(options.rPath, parser, 0);
    seqan::getArgumentValue(options.gPath, parser, 1);


    return seqan::ArgumentParser::PARSE_OK;

}

int main(int argc, char const ** argv)
{
    std::cerr << "Encapsulated version: Mapping reads efficiently" << std::endl;
    (void)argc;
    // Parse the command line.
    Options options;
    seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;
    Mapper<> mapper(options);
    mapper.printParm();
    map(mapper);

    return 0;
}
