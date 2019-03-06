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

#ifndef SEQAN_HEADER_BASE_H
#define SEQAN_HEADER_BASE_H

#include <seqan/seq_io.h>
#include <seqan/stream.h>
#include <seqan/index.h>
#include <seqan/store.h>
#include <seqan/basic.h>
#include <seqan/arg_parse.h>
#include <seqan/vcf_io.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <bitset>
#include <climits>
#include <thread>
#include <chrono>
#include <atomic>   
#include <iomanip>
#include <functional>   // for std::ref()
#include <chrono>
#include <thread>
#include <ctime>
   
#include "shape_extend.h"
#include "index_extend.h"

using namespace seqan;

//===================================================================
// variable and type def
//===================================================================

Options::Options():
        kmerLen(base_shape_len_),
        MiKmLen(base_shape_weight_),
        rPath(""),
        gPath(""),
        oPath("mapper_result.txt"),
        Sensitive(false),
        sensitivity(1),
        thread(16){}

std::string Options::getGenomePath() const {return gPath;};
std::string Options::getReadPath() const {return rPath;};
std::string Options::getOutputPath() const {return oPath;};
int Options::print()
{
    
    std::cerr << "kmerLen " << kmerLen << std::endl
              << "MiKmLen " << MiKmLen << std::endl
              << "reads path " << rPath << std::endl
              << "genomes Path " << gPath << std::endl
              << "output path " << oPath << std::endl
              << "Sensitive " << Sensitive << std::endl;
    return 0;
}

/*
 * flip strand from 0, 1 to -1, 1;
 * strand = 0, 1, other values is not allowed
 * return -1 , 1
 */
inline uint64_t _nStrand(uint64_t strand)
{
    return (strand << 1) - 1;
}
/*
 * flip the coordinates if strand = -1
 * do nothing if strand = 0
 * len is the length of the sequences;
 */
inline uint64_t _flipCoord (uint64_t coord, uint64_t len, uint64_t strand)
{
    return len * strand - _nStrand(strand) * coord;
}

MapParm::MapParm():
        blockSize(base_block_size_),
        delta(base_delta_),
        threshold(base_threshold_),
        kmerStep(base_kmer_step_),
        shapeLen(base_shape_len_),
        sensitivity(0),
        anchorDeltaThr(),
        minReadLen(1000),
        listN(1),
        listN2(1),
        alpha(base_alpha_),
        alpha2(0.5),
        anchorLenThr(0.02),                  // anchors with lenghth > this parameter is pushed into the queue
        rcThr(0.8),                        // when max anchors in the queue with length < this parameters, reverse complement search will be conducted
        cordThr(0.8),
        senThr(0.8),
        clsThr(0.1)
    {}
MapParm::MapParm(unsigned bs, unsigned dt, unsigned thr, 
            unsigned ks, unsigned sl, unsigned st,
            unsigned ad, unsigned mr, unsigned listn,
            unsigned listn2,
            float ap, float ap2, float alt, float rt, float ct, float sent, float clst):
        blockSize(bs),
        delta(dt),
        threshold(thr),
        kmerStep(ks),
        shapeLen(sl),
        sensitivity(st),
        anchorDeltaThr(ad),
        minReadLen(mr),
        listN(listn),
        listN2(listn2),
        alpha(ap),
        alpha2(ap2),
        anchorLenThr(alt),                  // anchors with lenghth > this parameter is pushed into the queue
        rcThr(rt),                        // when max anchors in the queue with length < this parameters, reverse complement search will be conducted
        cordThr(ct),
        senThr(sent),
        clsThr(clst)
        {} 
MapParm::MapParm(MapParm & parm):
        blockSize(parm.blockSize),
        delta(parm.delta),
        threshold(parm.threshold),
        kmerStep(parm.kmerStep),
        shapeLen(parm.shapeLen),
        sensitivity(parm.sensitivity),
        anchorDeltaThr(),
        minReadLen(parm.minReadLen),
        listN(parm.listN),
        listN2(parm.listN2),
        alpha(parm.alpha),
        alpha2(parm.alpha2),
        anchorLenThr(parm.anchorLenThr),
        rcThr(parm.rcThr),
        cordThr(parm.cordThr),
        senThr(parm.senThr),
        clsThr(parm.clsThr)
        {}

std::ifstream::pos_type _filesize(const char* filename)
{
    std::ifstream in(filename,  std::ifstream::binary | std::ifstream::ate);
    return in.tellg(); 
}

int readRecords_block (StringSet<CharString> & ids, StringSet<String<Dna5> > & reads, String<int> & lens, SeqFileIn & fin, int blockSize)
{
    int start = length(reads);
    readRecords(ids, reads, fin, blockSize);
    for (unsigned k = 0; k < length(reads) - start; k++)
    {
        appendValue (lens, length(reads[start + k]));
    }
    return 0;
}

/*
 *[]::lr
 */
int PMRecord::loadRecord(Options & options)
{
    double time = sysTime();
    SeqFileIn gFile(toCString(options.gPath));
    double fileSize = _filesize (toCString(options.gPath));
    bool flag = false;
    unsigned seqCount = 0;
    double currentFileSize = 0;
    std::fstream fin (toCString(options.gPath), std::fstream::in);
    StringSet<String<char> > dotstatus;
    resize(dotstatus, 3);
    dotstatus[0] = ".   ";
    dotstatus[1] = "..  ";
    dotstatus[2] = "... ";
#pragma omp parallel
{
    #pragma omp sections
    {
        #pragma omp section
        {
            //unsigned preSeqCount = 0;
            String <char> probar;
            float prepercent = 0, percent = 0, showpercent = 0, v = 0.87 ;
            unsigned k = 1;
            while (!flag)
            {
                prepercent = percent;
                percent = currentFileSize / fileSize * 100;
                percent = (percent > 100)?prepercent:percent;
                showpercent += v;
                showpercent = (showpercent > percent)?percent:showpercent;
                
                std::cerr << "                                                            \r";
                if (seqCount > 2)
                {
                    std::cerr << ">>Read genomes" << dotstatus[(k - 1)/10 %3] << "            " << seqCount << "/" << std::setprecision(2) << std::fixed << showpercent << "%\r";
                }
                else
                {
                    std::cerr << ">>Read genomes" << dotstatus[(k - 1)/10 %3] << "\r";
                }
                std::this_thread::sleep_for(std::chrono::milliseconds(100));
                k++;
            }
        }
        #pragma omp section
        {
            PMRecord<TDna>::RecId tmp_id;
            PMRecord<TDna>::RecSeq tmp_seq;
            while (!atEnd(gFile))
            {
                clear (tmp_id);
                clear (tmp_seq);
                readRecord (tmp_id, tmp_seq, gFile);
                currentFileSize += length(tmp_seq);
                appendValue (id2, tmp_id);
                appendValue (seq2, tmp_seq);
                ++seqCount;
            }
            flag = true;
        }
    }
}
    std::cerr << "--Read genomes                "<< length(seq2) <<"/100%                   \n";
    std::cerr << "  File: " << options.gPath ;
    std::cerr << "  Elapsed time [s] " << sysTime() - time << std::endl;
    return 0;
}

PMRecord::PMRecord(Options & options)
{
    readPath = options.rPath;
    genomePath = options.gPath;
    loadRecord(options);
}

inline void Anchors::init(AnchorType val, unsigned range)
{
    for (unsigned k = 0; k < range; k++)
        seqan::appendValue(set,val);
}

inline void Anchors::init(int length)
{
    clear(set);
    seqan::appendValue(set, 0);
    (void)length;
}

inline void Anchors::init()
{
    init(AnchorType(0), size());
}

inline void Anchors::setAnchor(unsigned p, 
    Anchors::AnchorType pos1,  Anchors::AnchorType pos2)
{
    set[p] = (pos1 << AnchorBase::bit) + pos2;
}

inline Anchors::AnchorType Anchors::getPos1(unsigned p) const 
{
    return set[p] >> AnchorBase::bit;
}

inline Anchors::AnchorType Anchors::getPos2(unsigned p) const
{
    return set[p] & AnchorBase::mask;
}

inline Anchors::AnchorType Anchors::deltaPos1(unsigned p1, unsigned p2)
{
    return (set[p1] >> AnchorBase::bit) - (set[p2] >> AnchorBase::bit);
}

inline Anchors::AnchorType Anchors::deltaPos2(unsigned p1, unsigned p2)
{
    return AnchorBase::mask & (set[p1] - set[p2]);
}
inline Anchors::Iter Anchors::begin()
{
    return seqan::begin(set);
}
inline Anchors::Iter Anchors::end()
{
    return seqan::end(set);
}
inline void Anchors::sort(Anchors::Iter sortBegin, Anchors::Iter sortEnd)
{
    std::sort(sortBegin, sortEnd);
}
inline void Anchors::sortPos2(Anchors::Iter sortBegin, Anchors::Iter sortEnd){
    AnchorBase::AnchorType mask = AnchorBase::mask;
    std::sort(sortBegin, sortEnd,
    [& mask](AnchorBase::AnchorType & a, 
                AnchorBase::AnchorType & b)
    {
        return (a & mask) < (b & mask);
    }) ;
}
inline void Anchors::appendValue(Anchors::AnchorType val)
{
    seqan::appendValue(set,val);
}
unsigned Anchors::size() const 
{
    return AnchorBase::size;
};
AnchorType & Anchors::operator [](unsigned p)
{
    return set[p];
};
unsigned Anchors::length() 
{
    return seqan::length(set);
};

void MapParm::print()
{
    std::cerr << "blockSize " << blockSize << std::endl
            << "alpha " << alpha << std::endl
            << "alpha2 " << alpha2 << "\n"
            << "listN " << listN << "\n"
            << "listN2 " << listN2 << "\n"
            << "senThr " << senThr << "\n"
            << "delta " << delta << std::endl
            << "threshold " << threshold << std::endl
            << "kmerStep " << kmerStep << std::endl
            << "shapeLen " << shapeLen << std::endl
            //<<  "sensitivity " << sensitivity << "\n"
            << "anchorDeltaThr " << anchorDeltaThr << "\n"
            << "minReadLen " << minReadLen << "\n"
            << "anchorLenThr" << anchorLenThr << "\n"
            << "rcThr " << rcThr << "\n"
            << "cordThr" << cordThr << "\n";
}

static const String<Dna5> _complt = "tgcan";
inline void _compltStr(String<Dna5> & str, String<Dna5> & res)
{
    resize(res, length(str));
    for (unsigned k = 0; k < length(str); k++)
     //   res[k]=_complt[str[k] - 'A'];
        res[k] = _complt[(unsigned)ordValue(str[k])];
}

inline void _compltRvseStr(String<Dna5> & str, String<Dna5> & res)
{
    resize(res, length(str));
    for (unsigned k = 0; k < length(str); k++)
     //   res[k]=_complt[str[k] - 'A'];
    {
        res[k] = _complt[(unsigned)ordValue(str[length(str) - k - 1])];
        //std::cout << (unsigned)ordValue(str[length(str) - k - 1]) << std::endl;
    }
}

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

#endif