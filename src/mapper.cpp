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
    std::ofstream of2("mapper_result2.txt");
    //unsigned strand;
    for (unsigned k = 0; k < length(cordSet); k++)
    {
        if (empty(cordSet[k]))
            of2 << k << " th Strand " << " 2 length " << length(reads()[k]) << "\nlength of cords -\n\n";
        else
        {
            //if (_DefaultCord.getCordStrand(back(cordSet[k]))) 
            //    strand = 1;
            //else 
            //    strand = 0;
            unsigned cordCount = 0;
            unsigned first = 0;
            unsigned cover = 0;
            unsigned firstCord = 1;
            for (unsigned j = 1; j < length(cordSet[k]); j++)
            {
                if (firstCord)
                {
                    firstCord = 0;
                    if (j != 1)
                        of2 << "\n";
                    of2 << k << " th Strand " << _DefaultHit.getStrand(cordSet[k][j]) << " length " 
                        << length(reads()[k]) << "\nlength of cords " << "\n";
                }
                of2 << j << " " << _DefaultCord.getCordY(cordSet[k][j]) << " " 
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
                   
                    of2 << "coverage " << (float)cover / (length(reads()[k])) << "\n";
                    if (j < length(cordSet[k]) - 1)
                    {
                        
                    }   
                    cordCount =0;
                    first = 0;
                    cover = 0;
                    firstCord = 1;
                }
                
 
            }
            of2 << "\n";
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

/*
 * print the cords without detailes 
 */
template <typename TDna, typename TSpec>
void Mapper<TDna, TSpec>::printCordsRaw()
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
                    //of << record.id1[k] << " " << length(cordSet[k]) << " "
                    of <<"@S1_"<< k+1 << " " << length(reads()[k]) << " "
                    << _DefaultCord.getCordY(cordSet[k][j]) << " " << length(cordSet[k]) << " x " 
                    << _getSA_i1(_DefaultCord.getCordX(cordSet[k][j])) << " " << cordCount << " "
                    << _getSA_i2(_DefaultCord.getCordX(cordSet[k][j]))  << " " 
                    //<< cordSet[k][j] 
                    << "\n";   
                    //flag = false;
                    cordCount = 0;
                }
                cordCount++;
            }
            //_DefaultCord.print(cordSet[k],of);
        }
    }
    std::cerr << ">Write results to disk          " << std::endl;
    std::cerr << "    End writing results. Time[s]" << sysTime() - time << std::endl; 
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

/*
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
    rawMapAllComplex2Parallel<TDna, TSpec>(mapper.index(), mapper.reads(), mapper.genomes(),
                         mapper.mapParm(), mapper.hits(), mapper.cords());
    //mapper.printHits();
    mapper.printCordsAll();
    mapper.printCordsRaw();
    std::cerr << length(mapper.cords()) << " " << length(mapper.reads()) << " \n";
    std::cerr << "Time in sum[s] " << sysTime() - time << std::endl;
}
*/

/*
 *[]::map
 */
template <typename TDna, typename TSpec>
int map(Mapper<TDna, TSpec> & mapper)
{
    std::cerr << "[]::map_1\n";
    //printStatus();
    omp_set_num_threads(mapper.thread());
    StringSet<String<int> > f2;
    mapper.createIndex(false); // true: destroy genomes string to reduce memory footprint
    createFeatures(mapper.genomes(), f2, mapper.thread());
    //uint64_t blockSize = 10000;
    //uint64_t lenSum;
    SeqFileIn rFile(toCString(mapper.readPath()));
    //while (!atEnd(rFile))
    //{
        double time = sysTime();
        std::cerr <<  ">reading reads from " << mapper.readPath() << "\r";
        readRecords(mapper.readsId(), mapper.reads(), rFile);//, blockSize);
        std::cerr << ">end reading " <<sysTime() - time << "[s]" << std::endl;
        std::cerr << ">mapping " << length(mapper.reads()) << " reads to reference genomes"<< std::endl;
   //     if (mapper.thread() < 2)
   //     {
   //         rawMap_dst<TDna, TSpec>(mapper.index(), mapper.reads(), mapper.genomes(), mapper.mapParm(), mapper.cords());
   //     }
   //     else
   //     {
            rawMap_dst2_MF<TDna, TSpec>(mapper.index(), f2, mapper.reads(), mapper.mapParm(), mapper.cords(), mapper.thread());
        //clear (mapper.reads());   
    //    }
        
        //clear (mapper.readsId());
    //}
    //mapper.printCordsAll();
    mapper.printCordsRaw2();
    return 0;
    //std::cerr << length(mapper.cords()) << " " << length(mapper.reads()) << " \n";
}

/*
 * Mapping with extra processing of gaps
 */
template <typename TDna, typename TSpec>
int map2(Mapper<TDna, TSpec> & mapper)
{
    omp_set_num_threads(mapper.thread());
    StringSet<String<int> > f2;
    createFeatures(mapper.genomes(), f2, mapper.thread());
    mapper.createIndex(); 
    SeqFileIn rFile(toCString(mapper.readPath()));
    double time = sysTime();
    
    std::cerr <<  ">reading reads from " << mapper.readPath() << "\r";
    readRecords(mapper.readsId(), mapper.reads(), rFile);//, blockSize);
    std::cerr << ">end reading " <<sysTime() - time << "[s]" << std::endl;
    
    std::cerr << ">mapping " << length(mapper.reads()) << " reads to reference genomes"<< std::endl;
    rawMap_dst2_MF<TDna, TSpec>(mapper.index(), f2, mapper.reads(), mapper.mapParm(), mapper.cords(), mapper.thread());
    mapGaps (mapper.genomes(), mapper.reads(), mapper.coords());
    
    mapper.printCordsRaw2();
    return 0;
}

int main(int argc, char const ** argv)
{
    double time = sysTime();

    std::cerr << "Encapsulated version: Mapping reads efficiently" << std::endl;
    (void)argc;
    Options options;
    seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;
    Mapper<> mapper(options);
    //mapper.printParm();
    map(mapper);
    std::cerr << "Time in sum[s] " << sysTime() - time << std::endl;

    return 0;
}
