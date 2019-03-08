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
#ifndef SEQAN_HEADER_PACMAPPER_H
#define SEQAN_HEADER_PACMAPPER_H

#include "base.h"
#include "index_util.h"
#include "f_io.h"
#include "mapparm.h"

class Mapper {
    PMRecord      record;
    MapParm     parm;
    LIndex      qIndex;
    StringSet<String<uint64_t> > cordSet;
    std::ofstream of;
    unsigned _thread;
    String<int> rlens;
    StringSet<String<uint64_t> > clip_set;
    StringSet<String<BamAlignmentRecordLink> > bam_records;
    std::string outputPrefix;

public:
    Mapper();
    Mapper(Options & options);
    StringSet<String<Dna5> > & reads() {return record.seq1;}             
    StringSet<String<Dna5> > & genomes() {return record.seq2;}             
    MapParm & mapParm() {return parm;}
    LIndex & index() {return qIndex;}
    StringSet<String<uint64_t> > & cords() {return cordSet;}            //returns cord set 
    
    void printCords(std::ostream & );
    void printCordsRaw();
    void printCordsRaw2();
    int print_vcf();
    int createIndex(bool = false);
    unsigned sens(){return parm.sensitivity;}
    unsigned & thread(){return _thread;}
    CharString & readPath(){return record.readPath;}
    CharString & genomePath(){return record.genomePath;}
    StringSet<CharString> & readsId(){return record.id1;}
    StringSet<CharString> & genomesId(){return record.id2;}
    StringSet<String<uint64_t> > & getClips(){return clip_set;}
    String<int> & readLens(){return rlens;}
    std::ofstream & getOf() {return of;}
    std::string & getOutputPrefix(){return outputPrefix;}
    StringSet<String<BamAlignmentRecordLink> > & getBamRecords() {return bam_records;}
};

#endif
