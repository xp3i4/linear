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

template <typename TDna, typename TSpec>
Mapper<TDna, TSpec>::Mapper(Options & options):
    record(options),
    qIndex(genomes()),
    of(toCString(options.getOutputPath()))
{
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
        std::cerr << "[mapper thread] " << _thread << "\n";
        /*
//map tunning
        parm.listN = options.listN;
        parm.listN2 = options.listN2;
        parm.alpha = options.alpha;
        parm.alpha2 = options.alpha2;
        parm.cordThr = options.cordThr;
        parm.senThr = options.senThr;
        */
}

template <typename TDna, typename TSpec>
int Mapper<TDna, TSpec>::createIndex(bool efficient)
{
    std::cerr << ">[Creating index] \n";
    createHIndex(genomes(), qIndex, _thread, efficient);
    return 0;
}

/*
 * print all cords with cordinates
 */
template <typename TDna, typename TSpec>
void Mapper<TDna, TSpec>::printCordsRaw2()
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
                    of << "\n" << record.id1[k] << " " << length(cordSet[k]) << " "
                    << _DefaultCord.getCordY(cordSet[k][j]) << " " << length(reads()[k]) << " x " 
                    << _getSA_i1(_DefaultCord.getCordX(cordSet[k][j])) << " " << cordCount << " "
                    << _getSA_i2(_DefaultCord.getCordX(cordSet[k][j]))  << " " 
                    ;  
                    //flag = false;
                    cordCount = 0;
                }
                of << _DefaultCord.getCordY(cordSet[k][j]) <<" " << 
                _getSA_i2(_DefaultCord.getCordX(cordSet[k][j])) << " | ";
                cordCount++;
            }
            //_DefaultCord.print(cordSet[k],of);
        }
    }
    std::cerr << ">Write results to disk          " << std::endl;
    std::cerr << "    End writing results. Time[s]" << sysTime() - time << std::endl; 
}

template <typename TDna, typename TSpec>
int map(Mapper<TDna, TSpec> & mapper)
{
    //printStatus();
    omp_set_num_threads(mapper.getThreads());
    StringSet<String<int> > f2;
    createFeatures(mapper.genomes(), f2, mapper.getThreads());
    
    mapper.createIndex(false); // this function will destroy genomes string during the creation to reduce memory footprint
    SeqFileIn rFile(toCString(mapper.readPath()));
    double time = sysTime();
    std::cerr <<  ">reading reads from " << mapper.readPath() << "\r";
    readRecords(mapper.readsId(), mapper.reads(), rFile);//, blockSize);
    std::cerr << ">end reading " <<sysTime() - time << "[s]" << std::endl;
    std::cerr << ">mapping " << length(mapper.reads()) << " reads to reference genomes"<< std::endl;
    rawMap_dst2_MF<TDna, TSpec>(mapper.index(), f2, mapper.reads(), mapper.mapParm(), mapper.cords(), mapper.getThreads());
    mapper.printCordsRaw2();
    return 0;
    //std::cerr << length(mapper.cords()) << " " << length(mapper.reads()) << " \n";
}

int g_test1 (String<Dna5> & seq)
{
    double time = sysTime ();
    GIndex g_index;
    for (unsigned k = 0; k < 100000; k++)
    {
        g_createDir(seq, k, k+1000, g_index);
        hashInit(g_index.shape, begin(seq) + k);
        unsigned count = 0;
        for (unsigned j = k; j < k + 1000; j++)
        {
            if (++count == 10)
            {
                hashNext(g_index.shape, begin(seq) + j);
                if (_defaultGNode.getXValue(g_index.g_hs[g_index.g_dir[g_index.shape.XValue]]) != g_index.shape.XValue)
                {
                    std::cerr << "[]::g_test::error " << k << " " << j << "\n";
                    return 1;
                }
                count = 0;
            }
            else
            {
                hashNexth(g_index.shape, begin(seq) + j);
            }
        }
    }
    std::cerr << "[]::g_test " << sysTime () - time << "\n";
    return 0;
}

int g_test_2 (String<Dna5> & seq)
{
    std::cerr << "[]::g_test \n";
    double time = sysTime ();
    Gap gap (100, 700, 200, 500);
    String <uint64_t> tile;
    for (uint64_t k = 0; k <1; k++)
    {
        std::cout << "[]:: " << k << "\n";
        gap.setValue(k*10, k*10 + 1000, k*10, k*10 + 1000);
        mapGap(seq, seq, gap, tile, 192);
    }
    std::cerr << "[]::g_test " << length(seq) << " " << sysTime () - time << "\n";
    return 0;
}
/*
 * Generate gaps manually and test mapGaps
 */
int g_test3(StringSet<String<Dna5> > & seqs)
{
    double time = sysTime();
    std::cerr << "[]::g_test \n";
    String<uint64_t> cords;
    String<uint64_t> tile;
    unsigned thd_gap = 200;
    unsigned window_size = 192;
    std::cerr << "[] " << length(cords) << "\n";
    for (unsigned k = 0; k < 45; k++)
    {
        appendValue(cords, _DefaultCord.createCord(k * window_size, k * window_size));
    }
    for (unsigned k = 50; k < 53; k++)
    {
        appendValue(cords, _DefaultCord.createCord(k * window_size, k * window_size));
    }
    for (unsigned k =  55; k < 100; k++)
    {
        appendValue(cords, _DefaultCord.createCord(k * window_size, k * window_size));
    }
    mapGaps (seqs, seqs[0], cords, thd_gap, 192);
    std::cerr << "[]::g_test " << sysTime () - time << "\n";
    return 0;
}

/*
 * read genomes and read from files
 *  insert 'deletions' to read[0]  
 * map read
 * map gaps
 */
int g_test(Mapper & mapper)
{
    double time = sysTime();
    std::cerr << "[]::g_test \n";
    String<uint64_t> cords;
    String<uint64_t> tile;
    StringSet<String<TDna> > reads;
    unsigned thd_gap = 200;
    unsigned window_size = 192;
    unsigned g_start = 4000; // gap start and end
    unsigned g_end = 5000;
    omp_set_num_threads(mapper.getThreads());
    StringSet<String<int> > f2;
    createFeatures(mapper.genomes(), f2, mapper.getThreads());
    mapper.createIndex(); 
    SeqFileIn rfile (toCString(mapper.readPath()));
    readRecords(mapper.readsId(), mapper.reads(), rfile);
    resize (reads, 1);
    for (unsigned k = 0; k < g_start; k++)
    {
        appendValue(reads[0], mapper.reads()[0][k]);
    }
    for (unsigned k = g_end; k < length(mapper.reads()[0]); k++)
    {
        appendValue(reads[0], mapper.reads()[0][k]);
    }
    rawMap_dst2_MF<TDna, TSpec>(mapper.index(), f2, reads, mapper.mapParm(), mapper.cords(), mapper.getThreads());
    mapGaps (mapper.genomes(), reads[0], mapper.cords()[0], window_size, window_size);
    mapper.printCordsRaw2();
    std::cerr << "[]::g_test " << sysTime () - time << "\n";
    return 0;
}

int main(int argc, char const ** argv)
{
    std::cerr << "Encapsulated version: Mapping reads efficiently" << std::endl;
    std::cerr << "Test units [] \n";
    (void)argc;
    // Parse the command line.
    Options options;
    seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;
    double t=sysTime();
    Mapper<> mapper(options);
    
    //g_test(mapper);

    std::cerr << "results saved to " << options.getOutputPath() << "\n";
    
    return 0;
}
