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


#include "base.h"
#include "pmpfinder.h"

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
        setHelpText(parser, 0, "Reads file .fa, .fasta");

        addArgument(parser, seqan::ArgParseArgument(
            seqan::ArgParseArgument::INPUT_FILE, "genome", true));
        setHelpText(parser, 1, "Reference file .fa, .fasta");

        addSection(parser, "Mapping Options");
        addOption(parser, seqan::ArgParseOption(
            "o", "output", "choose output file.",
             seqan::ArgParseArgument::STRING, "STR"));
        addOption(parser, seqan::ArgParseOption(
            "s", "sensitive", "Sensitive mode. Default closed",
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


        return seqan::ArgumentParser::PARSE_OK;

    }

/*
    typedef Iterator<String<Dna5> >::Type TIter;
    typedef Shape<Dna5, Minimizer<30> > TShape;
    typedef Shape<Dna5, UngappedShape<30> > TShape_u;


typedef Index<StringSet<Dna5String>, IndexQGram<Minimizer<30>, OpenAddressing > > TIndex;
typedef Index<StringSet<Dna5String>, IndexQGram<UngappedShape<30>, OpenAddressing > > TIndex_u;

int uTest3(StringSet<Dna5String> & reads, StringSet<Dna5String> & genome)
{
    TShape_u t_shape;
    TIndex_u index(genome);
    unsigned kmerLength = t_shape.span;
    uint64_t sum=0, count=0, p = 0;
    double time = sysTime();
    std::cout << "uTest3():\n";
    std::cout << "    fullDirLength " << _fullDirLength(index) << std::endl; 

    indexCreate(index, FibreSADir());
    std::cout << "    getSA start sysTime(): " << sysTime() - time<< std::endl;
    for(unsigned k = 0; k < length(reads); k++)
    {
        TIter it = begin(reads[k]);
        hashInit(t_shape, it);
        for (uint64_t j = 0; j < length(reads[k]) - kmerLength + 1; j++)
        //for (uint64_t j = 0; j < 100; j++)
        {
            hashNext(t_shape, it + j);
            p = getBucket(index.bucketMap, t_shape.hValue);
            for (uint64_t k = index.dir[p]; k < index.dir[p+1]; k++)
                sum ^= index.sa[k].i2;
        }
    }
    std::cout << "    sum = " << sum << " count = " << count << std::endl;
    std::cout << "    getSA end sysTime(): "<< sysTime() - time<< std::endl;
    std::cout << "    End uTest()" << std::endl;
    return 0;
}


int mTest3(StringSet<Dna5String> & reads, StringSet<Dna5String> & genome)
{
    TShape shape, shape1;
    TIndex index(genome);
    uint64_t sum = 0, p = 0;
    double time = sysTime();
    std::cout << "mTest3(): " << std::endl;
    _createQGramIndex(index);
    std::cout << "    getDir start sysTime(): " << sysTime() - time << std::endl;
    std::cout << "    length Dir = " << length(index.dir) << std::endl;
    std::cout << "    length Text = " << lengthSum(indexText(index)) << std::endl;
    std::cout << "    length SA = " << length(index.sa) << std::endl;
    for(uint64_t k = 0; k < length(reads); k++)
    {
        TIter it = begin(reads[k]);
        hashInit(shape, it);
        for (uint64_t j = 0; j < length(reads[k]) - shape.span + 1; j++)
        {
            hashNext(shape, it + j);
            p = getDir(index, shape);
            for (uint64_t n = _getBodyCounth(index.dir[p]); n < _getBodyCounth(index.dir[p + 1]); n++)
            {
                sum ^= _getSA_i2(index.sa[n]);
            }
        }
    }
    std::cout << "    sum = " << sum << std::endl;
    std::cout << "    getDir end sysTime(): " << sysTime() - time << std::endl;
    std::cout << "    End mTest1()" << std::endl;

    return 0;
}
*/

int main(int argc, char const ** argv)
{
    (void)argc;
    // Parse the command line.
    Options options;
    seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;
    double t=sysTime();
    Mapper<> mapper(options);
    //options.print(); 
    map(mapper);
    std::cerr << "total time " << sysTime() - t << std::endl;
    
    //mTest3(mapper.reads(), mapper.genomes());   
    //uTest3(mapper.reads(), mapper.genomes());   
    return 0;
}
