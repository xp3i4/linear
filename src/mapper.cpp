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
    std::cerr << ">>Create index \r";
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
    std::cerr << ">>Write results to disk        \r";
    double time = sysTime();
    unsigned cordCount = 0;
    CharString first_line = "";
    cordCount = 0;
    uint64_t readCordEnd;
    uint64_t seqsCordEnd;
    //std::cerr << "[]::printCordsRaw2 " << length(rlens) << " " << length(cordSet) << "\n";
    char main_icon_strand = '+', icon_strand = '+';
    int fflag = 0;
    for (unsigned k = 0; k < length(cordSet); k++)
    {
        /*
        if (k - prek > k_delta)
        {
            std::cerr << ">>Write results to disk        " << std::setprecision(2) << (float)k / length(cordSet) * 100 << "% \r";
            prek = k;
        }
        */
        if (!empty(cordSet[k]))
        {
            for (unsigned j = 1; j < length(cordSet[k]); j++)
            {
                if (_DefaultHit.isBlockEnd(cordSet[k][j-1]))
                {
                    for (unsigned i = j; ; i++)
                    {
                        if (_DefaultHit.isBlockEnd(cordSet[k][i]) || i == length(cordSet[k]) - 1)
                        {
                            readCordEnd = _DefaultCord.getCordY(cordSet[k][i]) + window_size;
                            seqsCordEnd = _getSA_i2(_DefaultCord.getCordX(cordSet[k][i])) + window_size;
                            break;
                        }
                    }
                    main_icon_strand = (_DefaultCord.getCordStrand(cordSet[k][j]))?'-':'+';
                    of  << first_line;
                    of  << record.id1[k] << " " 
                        //<< length(record.seqs1[k]) << " "
                        << rlens[k] << " "
                        << _DefaultCord.getCordY(cordSet[k][j]) << " " 
                        << std::min(readCordEnd, (uint64_t)rlens[k]) << " " 
                        << main_icon_strand<< " "
                        << record.id2[_getSA_i1(_DefaultCord.getCordX(cordSet[k][j]))] << " " 
                        << length(record.seq2[_getSA_i1(_DefaultCord.getCordX(cordSet[k][j]))]) << " "
                        << _getSA_i2(_DefaultCord.getCordX(cordSet[k][j]))  << " " 
                        << seqsCordEnd << "\n";
                    first_line = "\n";
                    cordCount = 0;
                    fflag = 1;
                }
                icon_strand = (_DefaultCord.getCordStrand(cordSet[k][j]))?'-':'+';
                CharString mark = "| ";
                if (icon_strand != main_icon_strand)
                    mark = "********** ";
                int64_t d = 0;//_DefaultCord.getCordY(cordSet[k][1]);
                int64_t d2 = 0;
                if (!fflag)
                {
                    d = (int) _DefaultCord.getCordX(cordSet[k][j] - cordSet[k][j - 1]);
                    d2 = (int) _DefaultCord.getCordY(cordSet[k][j] - cordSet[k][j - 1]);
                    
                }
                
                of  << mark  << _DefaultCord.getCordY(cordSet[k][j]) << " " 
                    << _getSA_i2(_DefaultCord.getCordX(cordSet[k][j])) << " " << d2 << " " << d << " \n";
                cordCount++;
                fflag = 0;
            }
        }
    }
    close(of);
    std::cerr << "--Write results to disk       100% Elapsed Time[s] " << sysTime() - time << std::endl;
}

template <typename TDna, typename TSpec>
int Mapper<TDna, TSpec>::print_vcf()
{
    return 0;
}

template <typename TDna, typename TSpec>
int Mapper<TDna, TSpec>::print_gff()
{
 //   for (int k = 0; k < length(cordSet); k++)
 //   {
 //       print_gff_()
 //   }
    return 0;
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

template <typename TDna, typename TSpec>
int rawMap_dst2_MF(typename PMCore<TDna, TSpec>::Index   & index,
            StringSet<String<short> > & f2,
            typename PMRecord<TDna>::RecSeqs      & reads,
            MapParm & mapParm,
            StringSet<String<uint64_t> > & cords,
            unsigned & threads,
            StringSet<String<TDna> > & seqs
            )
{
  
    typedef typename PMRecord<TDna>::RecSeq Seq;
    //double time=sysTime();
    float senThr = mapParm.senThr / window_size;
    float cordThr = mapParm.cordThr / window_size;
    MapParm complexParm = mapParm;
    complexParm.alpha = complexParm.alpha2;
    complexParm.listN = complexParm.listN2;
    //double time2 = sysTime();
#pragma omp parallel
{
    unsigned size2 = length(reads) / threads;
    unsigned ChunkSize = size2;
    Seq comStr;
    //Anchors anchors(Const_::_LLTMax, AnchorBase::size);
    Anchors anchors;
    typename PMRes::HitString crhit;
    StringSet<String<uint64_t> >  cordsTmp;
    StringSet< String<short> > f1;
    unsigned thd_id =  omp_get_thread_num();
    if (thd_id < length(reads) - size2 * threads)
    {
        ChunkSize = size2 + 1;
    }
    resize(cordsTmp, ChunkSize);
    resize(f1, 2);
    unsigned c = 0;
    
    String<uint64_t>  g_hs;
    String<uint64_t>  g_anchor;
    resize (g_hs, 1ULL << 20);
    resize (g_anchor, 1ULL<<20);

    #pragma omp for
    for (unsigned j = 0; j < length(reads); j++)
    {
        if (length(reads[j]) >= mapParm.minReadLen)
        //if (1)
        {
            std::cout << "[]::rawmap::j " << j <<"\n";
            float cordLenThr = length(reads[j]) * cordThr;
            _compltRvseStr(reads[j], comStr);
            createFeatures(begin(reads[j]), end(reads[j]), f1[0]);
            createFeatures(begin(comStr), end(comStr), f1[1]);
            anchors.init(1);
            clear(crhit);
            mnMapReadList<TDna, TSpec>(index, reads[j], anchors, mapParm, crhit);
            path_dst(begin(crhit), end(crhit), f1, f2, cordsTmp[c], cordLenThr);
            if (_DefaultCord.getMaxLen(cordsTmp[c]) < length(reads[j]) * senThr)// && 
            //_DefaultCord.getMaxLen(cordsTmp[c]) > 0)
            {
                clear(cordsTmp[c]);
                anchors.init(1);
                clear(crhit);
                mnMapReadList<TDna, TSpec>(index, reads[j], anchors, complexParm, crhit);
                path_dst(begin(crhit), end(crhit), f1, f2, cordsTmp[c], cordLenThr);
            }   
            mapGaps(seqs, reads[j], comStr, cordsTmp[c], g_hs, g_anchor, f1, f2, 300, 192);
            ///align (seqs, reads[j], comStr, cordsTmp[c]);
        }   
        
        c += 1;
    } 
    #pragma omp for ordered
    for (unsigned j = 0; j < threads; j++)
        #pragma omp ordered
        {
            append(cords, cordsTmp);
        }
}
    //std::cerr << "    End raw mapping. Time[s]: " << sysTime() - time << std::flush << std::endl;
    return 0;
}


/*
 *[]::map
 */
template <typename TDna, typename TSpec>
int map(Mapper<TDna, TSpec> & mapper)
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
        rawMap_dst2_MF<TDna, TSpec>(mapper.index(), 
                                    f2, 
                                    mapper.reads(), 
                                    mapper.mapParm(), 
                                    mapper.cords(), 
                                    mapper.thread(), 
                                    mapper.genomes());
        time2 = sysTime() - time2;
        std::cerr <<  "--Map::file_I/O+Map block "<< k << " Size " << length(mapper.reads()) << " Elapsed Time[s]: file_I/O " << time1 << " map "<< time2 << "\n";
        k++;
    }
    //clear (mapper.genomes());
    mapper.index().clear();
    mapper.printCordsRaw2();
    //mapper.print_gff();
    
    return 0;
}
/*
int map (Mapper<> mapper)
{
    if (mapper.raw)
    {
        return mapRaw(mapper);
    }
    else
    {
        return mapAlign(mapper);
    }
}
*/
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
    Mapper<> mapper(options);
    omp_set_num_threads(mapper.thread());
    map(mapper);
    //mapper.print_vcf();
//    align(mapper.genomes(), mapper.reads(), mapper.cords());
    std::cerr << "  Result File: \033[1;31m" << options.oPath << "\033[0m" << std::endl;
    std::cerr << "Time in sum[s] " << sysTime() - time << std::endl;
    return 0;
}
