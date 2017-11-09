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

<<<<<<< Updated upstream:mapper.cpp
=======
template <typename TDna>   
bool test_0_1(StringSet<String<TDna> > & seqs, StringSet<String<TDna> > & seqs2, bool flag1 = true, bool flag2 = true, bool flag3 = true)
{
    std::cerr << "test_0_1 count dir test\n";
    
    const unsigned shapelength = 25;
    double time;
   /* 
    if (flag1)
    {
        uint64_t sum = 0, counth = 0, county2, preX = 0;
        bool vflag = false;
        HIndex<shapelength> index;
        typename HIndexBase<shapelength>::TShape shape;
        //=====HIndex
        createHIndex(seqs2, index);
        time = sysTime();
        for(uint64_t j = 0; j < length(seqs); j++)
        {
            hashInit(shape, begin(seqs[j]));
            for (uint64_t k =0; k < length(seqs[j]) - shape.span + 1; k++)
            {
                //std::cout << " k = " << k << "\n";
                if(ordValue(*(begin(seqs[j]) + k + shape.span - 1)) == 4)
                {
                    k += hashInit(shape, begin(seqs[j]) + k);
                }
                hashNext(shape, begin(seqs[j]) + k);
                uint64_t pos = getXDir(index.xstr, shape.XValue, shape.YValue, vflag);
                uint64_t m = pos;
                do
                {
                    if(_DefaultHs.getHsBodyY(index.ysa[m]) == shape.YValue)
                    {
                        counth++;
                        continue;
                    }
                    if (_DefaultHs.getHsBodyY(index.ysa[m]) < shape.YValue)
                        break;
                }
                while (_DefaultHs.isBody(index.ysa[++m]));
            }
        }
        std::cerr << "      countm = " << counth << " [s]" << sysTime() - time << "\n";
    }
*/
    if(flag2)
    {
        //=====MIndex
        uint64_t countm = 0;
        uint64_t mask_msa = (1ULL << 40) - 1;
        typedef Index<StringSet<String<Dna5> >, IndexQGram<Minimizer<shapelength>, OpenAddressing > > TIndex_m;
        TIndex_m index_m(seqs2);
        Shape<Dna5, Minimizer<shapelength> > shape_m;
        time = sysTime();
        
        std::cerr << "  MIndex \n";
        _createQGramIndex(index_m);
        std::cerr << "      createing index time [s] " << sysTime() - time << "\n";
        
        time = sysTime();
        for(uint64_t j = 0; j < length(seqs); j++)
        {
            hashInit(shape_m, begin(seqs[j]));
            for (uint64_t k =0; k < length(seqs[j]) - shape_m.span + 1; k++)
            {
                //std::cout << " k = " << k << "\n";
                if(ordValue(*(begin(seqs[j]) + k + shape_m.span - 1)) == 4)
                {
                    k += hashInit(shape_m, begin(seqs[j]) + k);
                }
                hashNext(shape_m, begin(seqs[j]) + k);
                uint64_t pos = getDir(index_m, shape_m);
                if (_getBodyCounth(index_m.dir[pos + 1]) - _getBodyCounth(index_m.dir[pos]) == 0)
                    std::cout << "      MIndex:0 count " << j << " " << k << "\n";
                countm += _getBodyCounth(index_m.dir[pos + 1]) - _getBodyCounth(index_m.dir[pos]);
            }
        }
        std::cerr << "      countm = " << countm << " [s]" << sysTime() - time << "\n";  
    }

    if (flag3)
    {
        //=====OpenAddressing index
        typedef Shape<Dna5, UngappedShape<shapelength> > TShape_u;
        typedef Index<StringSet<String<Dna5> >, IndexQGram<UngappedShape<shapelength>, OpenAddressing > > TIndex_u;
        typedef Iterator<String<Dna5> >::Type TIter;

        TShape_u shape2;
        TIndex_u index2(seqs2);
        uint64_t sum2 = 0, p = 0, counto;
        indexCreate(index2, FibreDir());
        
        unsigned occ = 0;
        for(uint64_t k = 0; k < length(seqs); k++)
        {
            TIter it = begin(seqs[k]);
            hashInit(shape2, it);
            for (uint64_t j = 0; j < length(seqs[k]) - shape2.span + 1; j++)
            {
                hashNext(shape2, it + j);
                p = getBucket(index2.bucketMap, shape2.hValue);
                counto += index2.dir[p + 1] - index2.dir[p];
            }
        }   
        std::cout << "    counto = " << counto << "\n";
    }
    
    std::cerr << "    End test_0_1()\n\n";
    

    return true;
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
    double t=sysTime();
    Mapper<> mapper(options);
    map(mapper);

    return 0;
}
