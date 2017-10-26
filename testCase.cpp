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

template <typename TDna>   
void _testIndex(StringSet<String<TDna> > & seqs, typename Mapper<>::Index & index)   
{
    double time = sysTime();
    uint64_t sum = 0;
    for(uint64_t j = 0; j < length(seqs); j++)
    {
        hashInit(index.shape, begin(seqs[j]));
        for (uint64_t k =0; k < length(seqs[j]) - index.shape.span + 1; k++)
        {
           // if(ordValue(*(begin(seqs[j]) + k + shape.span - 1)) == 4)
           // {
           //     k += hashInit(shape, begin(seqs[j]) + k);
           //     if(k > length(seqs[j]) - shape.span + 1)
           //         break;
           // }
            
            hashNext(index.shape, begin(seqs[j]) + k);
            sum += getDir(index, index.shape);
        }
    }
    std::cout << sysTime() - time << " " << sum << std::endl;
}



template <typename TDna>   
bool _testIndexGetDir(StringSet<String<TDna> > & seqs, StringSet<String<TDna> > & seqs2)
{
    std::cerr << "testIndexGetDir1 \n";
    double time = sysTime();
    uint64_t sum = 0, county, county2, sum2=0;
    
    
    HIndex<25> index;
    createHIndex(seqs2, index);
    std::cerr << "time of creating index " << sysTime() - time << "\n";
    for(uint64_t j = 0; j < length(seqs); j++)
    {
        hashInit(index.shape, begin(seqs[j]));
        for (uint64_t k =0; k < length(seqs[j]) - index.shape.span + 1; k++)
        {
            if(ordValue(*(begin(seqs[j]) + k + index.shape.span - 1)) == 4)
            {
                k += hashInit(index.shape, begin(seqs[j]) + k);
            }
            
            hashNext(index.shape, begin(seqs[j]) + k);
//            uint64_t pos = getXDir(index.xstr, index.ysa, index.shape.XValue, index.shape.YValue) & ((1ULL << 32) - 1);
            uint64_t pos = getXDir(index.xstr, index.ysa, index.shape.XValue, index.shape.YValue);
            //if (pos != _DefaultXNodeBase._Empty_Dir_)
            //std::cout << index.shape.XValue << " " << index.shape.YValue << "\n";
            //if (pos <= length(index.ysa))
            //{
                
                //while(_DefaultHs.isBody(index.ysa[pos]))
                while(_DefaultHs.getHsBodyY(index.ysa[pos]) == index.shape.YValue)
                {
                    //if (_DefaultHs.getHsBodyY(index.ysa[pos]) == index.shape.YValue)
                    //{
                        //sum ^= index.ysa[pos];
                        //pos++;
                        //continue;
                        //break;
                    //}
                    //if (_DefaultHs.getHsBodyY(index.ysa[pos]) < index.shape.YValue)
                     //   break;
                    sum++;
                    pos++;
                }
            //}
            //else 
            //{
            //    std::cerr << k << " " << pos << " error " << length(hs) << "\n";
            //    return false;
            //} 
            //std::cerr << k << " " << length(seqs[j]) << std::endl;
        }
    }
    std::cout << "End testIndexGetDir time:" << sysTime() - time << " sum = " << sum << " " << sum / (float) lengthSum(seqs)<< " " << sum2 / (float) lengthSum(seqs) << std::endl;
   
    
    
    /*
    const unsigned shapelength = 25;
    typedef typename Iterator<String<TDna> >::Type TIter;
    Index<StringSet<String<TDna> >, IndexQGram<UngappedShape<shapelength>, OpenAddressing > > index2(seqs);
    Shape<TDna, UngappedShape<shapelength> > shape2;
    time = sysTime();
    //indexCreate(index2, FibreDir());
    createQGramIndex(index2);
    std::cerr << "Creating OpenAddressing index2 time " << sysTime() - time << "\n";
    
    time = sysTime();
    sum2 = 0;
    uint64_t p = 0;
    for(unsigned k = 0; k < length(seqs); k++)
    {
        TIter it = begin(seqs[k]);
        hashInit(shape2, it);
        for (uint64_t j = 0; j < length(seqs[k]) - shapelength + 1; j++)
        {
            hashNext(shape2, it + j);
            p = getBucket(index2.bucketMap, shape2.hValue);
            std::cerr << j << " " << p << " \n";

            for (uint64_t n = index2.dir[p]; n < index2.dir[p+1]; n++)
                sum2 ^= index2.sa[n].i2;
        }
    }
    std::cout << "    sum2 = " << sum2 << std::endl;
    std::cout << "    getSA end sysTime(): "<< sysTime() - time<< std::endl;
    */
    return true;
}

void mTest3(StringSet<String<Dna5> > & reads, StringSet<String<Dna5> > & genome)
{
    const unsigned shapelength = 25;
    typedef Shape<Dna5, Minimizer<shapelength> > TShape;
    typedef Index<StringSet<String<Dna5> >, IndexQGram<Minimizer<shapelength>, OpenAddressing > > TIndex;
    typedef Iterator<String<Dna5> >::Type TIter;

    TShape shape, shape1;
    TIndex index(reads);
    uint64_t sum = 0, p = 0;
    double time = sysTime();
    std::cout << "mTest3(): " << std::endl;
    _createQGramIndex(index);
    std::cout << "    getDir start sysTime(): " << sysTime() - time << std::endl;
    std::cout << "    length Dir = " << length(index.dir) << std::endl;
    std::cout << "    length Text = " << lengthSum(indexText(index)) << std::endl;
//    std::cout << "    length SA = " << length(index.sa) << std::endl << _getSortPara_i1(shape.weight) << " " << _getSortPara_i2(shape.weight);
     
    for(uint64_t k = 0; k < length(genome); k++)
    {
        TIter it = begin(genome[k]);
        hashInit(shape, it);
        for (uint64_t j = 0; j < length(genome[k]) - shape.span + 1; j++)
        {
            hashNext(shape, it + j);
            p = getDir(index, shape);
            sum += _getBodyCounth(index.dir[p + 1]) - _getBodyCounth(index.dir[p]);
            //for (uint64_t n = _getBodyCounth(index.dir[p]); n < _getBodyCounth(index.dir[p + 1]); n++)
            //{
            //    sum ^= _getSA_i2(index.sa[n]);
            //}
        }
    }
    std::cout << "    sum = " << sum << std::endl;
    std::cout << "    getDir end sysTime(): " << sysTime() - time << std::endl;
    std::cout << "    End mTest1()" << std::endl;
    
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
    //_createIndexDir(mapper.genomes(), mapper.index().dir, mapper.index().shape);
    //_createQGramIndexDirSA(mapper.genomes(), mapper.index().dir, mapper.index().bucketMap, mapper.index().sa, mapper.index().shape, true);
    //mapper.createIndex();
    //_testIndex(mapper.reads(), mapper.index());
    String<uint64_t> dir, bucketMap, sa, hs;
    XString xstr;
    //_createQGramIndexDirSA(mapper.genomes(), xstr, hs, mapper.index().shape, true);
    
   // for (unsigned k = 0 ; k < length(hs); k++)
   //     if (_DefaultHs.isBody(hs[k]))
   //         //std::cerr << "hs " << _DefaultHs.getHsBodyY(hs[k]) << std::endl;
   //         std::cerr << "hs " << ((hs[k] >> 41) & ((1ULL << 20) - 1)) << std::endl;
   //     else
   //         std::cerr << "0\n";
    //map(mapper);
    //_testIndexGetDir(mapper.reads(), mapper.index(), xstr, hs);
    _testIndexGetDir(mapper.reads(), mapper.genomes());
    mTest3(mapper.genomes(), mapper.reads());
    std::cerr << mapper.index().shape.weight << std::endl;
    std::cerr << "results saved to " << options.getOutputPath() << "\n";
    
    return 0;
}
