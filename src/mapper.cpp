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

template <typename TDna, typename TSpec>
int Mapper<TDna, TSpec>::createIndex(bool efficient)
{
    std::cerr << ">>Create index \r";
    createHIndex(genomes(), qIndex, _thread, efficient);
    return 0;
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
template <typename TDna, typename TSpec>
int print_align_sam (Mapper<TDna, TSpec> & mapper)
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
    char main_icon_strand = '+', icon_strand = '+';
    int fflag = 0;
    for (unsigned k = 0; k < length(cordSet); k++)
    {
        if (!empty(cordSet[k]))
        {
            for (unsigned j = 1; j < length(cordSet[k]); j++)
            {
                if (_DefaultHit.isBlockEnd(cordSet[k][j-1]))
                {
                    unsigned m = j; 
                    int main_strand_count = 0;
                    int block_len = 0;
                    ///>determine the main strand
                    while (!_DefaultHit.isBlockEnd(cordSet[k][m]))
                    {
                        if (_DefaultCord.getCordStrand(cordSet[k][m]))
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
                        if (_DefaultHit.isBlockEnd(cordSet[k][i]) || i == length(cordSet[k]) - 1)
                        {
                            readCordEnd = _DefaultCord.getCordY(cordSet[k][i]) + window_size;
                            seqsCordEnd = _getSA_i2(_DefaultCord.getCordX(cordSet[k][i])) + window_size;
                            break;
                        }
                    }
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
                ///>print the coordinates
                icon_strand = (_DefaultCord.getCordStrand(cordSet[k][j]))?'-':'+';
                CharString mark = "| ";
                if (icon_strand != main_icon_strand)
                    mark = "********** ";
                int64_t d = 0;//_DefaultCord.getCordY(cordSet[k][1]);
                int64_t d2 = 0;
                if (!fflag)
                {
                    d = (int64_t)_DefaultCord.getCordX(cordSet[k][j]) - (int64_t)_DefaultCord.getCordX(cordSet[k][j - 1]);
                    d2 = (int64_t)_DefaultCord.getCordY(cordSet[k][j]) - (int64_t)_DefaultCord.getCordY(cordSet[k][j - 1]);
                    
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

int print_clip_gff_(StringSet<String<uint64_t> > & clips, 
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

int print_clip_gvf_(StringSet<String<uint64_t> > & clips, 
              StringSet<CharString> & readsId, 
              StringSet<CharString> & genomesId,
              std::ofstream & of, 
              std::string outputPrefix)
{
    std::string file_path = outputPrefix + ".gvf";
    //std::cerr << "[]::filepath " << file_path << "\n";
    of.open(toCString(file_path));
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
                CharString genomeId = genomesId[_getSA_i1(_DefaultCord.getCordX(clips[i][j]))];
                if ((j >> 1) << 1 == j)
                {
                    of  << genomeId << " " << source << " " << type << " " << cord_x << " ";   
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

template <typename TDna, typename TSpec>
int print_clip_gff(Mapper<TDna, TSpec> & mapper)
{
    print_clip_gff_(mapper.getClips(), mapper.readsId(), mapper.getOf(), mapper.getOutputPrefix());
    return 0;
}

template <typename TDna, typename TSpec>
int print_clip_gvf(Mapper<TDna, TSpec> & mapper)
{
    print_clip_gvf_(mapper.getClips(), mapper.readsId(), mapper.genomesId(), mapper.getOf(), mapper.getOutputPrefix());
    return 0;
}

template <typename TDna, typename TSpec>
int rawMap_dst2_MF(typename PMCore<TDna, TSpec>::Index & index,
                   StringSet<String<short> > & f2,
                   typename PMRecord<TDna>::RecSeqs & reads,
                   MapParm & mapParm,
                   StringSet<String<uint64_t> > & cords,
                   StringSet<String<uint64_t> > & clips,
                   StringSet<String<TDna> > & seqs,
                   StringSet<String<BamAlignmentRecord> >& bam_records,
                   unsigned & threads,
                   int p1
                  )
{
  
    typedef typename PMRecord<TDna>::RecSeq Seq;
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
    StringSet<String<uint64_t> > clipsTmp;
    StringSet<String<BamAlignmentRecord> > bam_records_tmp;
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
            mnMapReadList<TDna, TSpec>(index, reads[j], anchors, mapParm, crhit);
            path_dst(begin(crhit), end(crhit), f1, f2, cordsTmp[c], cordLenThr);
            if (_DefaultCord.getMaxLen(cordsTmp[c]) < length(reads[j]) * senThr)
            {
                clear(cordsTmp[c]);
                anchors.init(1);
                clear(crhit);
                mnMapReadList<TDna, TSpec>(index, reads[j], anchors, complexParm, crhit);
                path_dst(begin(crhit), end(crhit), f1, f2, cordsTmp[c], cordLenThr);
            }   
            gap_len[thd_id] += mapGaps(seqs, reads[j], comStr, cordsTmp[c], g_hs, g_anchor, clipsTmp[c], f1, f2, 300, 192, p1);
            align_cords(seqs, reads[j], comStr, cordsTmp[c], bam_records_tmp[c]);
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
for (int k = 0; k < length(gap_len); k++)
{
    std::cout << "[] gap len " << gap_len[k] << " " << red_len[k] << " " << float(gap_len[k]) / red_len[k] << "\n";
}
    //std::cerr << "    End raw mapping. Time[s]: " << sysTime() - time << std::flush << std::endl;
    return 0;
}

/*
 *[]::map
 */
template <typename TDna, typename TSpec>
int map(Mapper<TDna, TSpec> & mapper, int p1)
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
                                    mapper.getClips(),
                                    mapper.genomes(),
                                    mapper.getBamRecords(),
                                    mapper.thread(), 
                                    p1
                                   );
        time2 = sysTime() - time2;
        std::cerr <<  "--Map::file_I/O+Map block "<< k << " Size " << length(mapper.reads()) << " Elapsed Time[s]: file_I/O " << time1 << " map "<< time2 << "\n";
        k++;
    }
    mapper.index().clear(); 
    mapper.printCordsRaw2();
    print_align_sam(mapper);
    print_clip_gvf(mapper);
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
    Mapper<> mapper(options);
    omp_set_num_threads(mapper.thread());
    map(mapper, options.p1);

    //mapper.print_vcf();
    std::cerr << "  Result Files: \033[1;31m" << options.oPath << "\033[0m" << std::endl;
    std::cerr << "                \033[1;31m" << (mapper.getOutputPrefix() + ".gvf") << "\033[0m" << std::endl;
    std::cerr << "                \033[1;31m" << (mapper.getOutputPrefix() + ".sam") << "\033[0m" << std::endl;
    std::cerr << "Time in sum[s] " << sysTime() - time << std::endl;
    return 0;
}
