// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2018, Knut Reinert, FU Berlin
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
// Author:      Enrico Seiler <enrico.seiler@fu-berlin.de>
// ==========================================================================

#include <algorithm>
#include <cassert>
#include <stdio.h>
#include <vector>

#include <seqan/seq_io.h>
#include <seqan/kmer.h>

static const uint32_t filterMetadataSize = 256;
static const uint8_t INT_WIDTH = 0x40;

using namespace seqan;

int main()
{
    //SeqAnBloomFilter<> filter (10, 1, 12, 16777216);
    uint64_t threads{3};
    uint64_t noBins{10};
    uint64_t kmerSize{3};
    uint64_t hashFunc{3};
    uint64_t bits{11829};


    // ==========================================================================
    // Test constructors
    // ==========================================================================
    std::cout << "Testing ctors" << '\n';

    typedef InterleavedBloomFilter TSpec;
    // typedef DirectAddressing       TSpec;

    // Empty default constructor
    KmerFilter<Dna, TSpec> ctor_empty;
    // Default constructor
    // KmerFilter<Dna, TSpec> ctor_default (noBins, hashFunc, kmerSize, bits);

    KmerFilter<Dna, TSpec> ctor_default (noBins, hashFunc, kmerSize, bits);
    KmerFilter<Dna, TSpec> ctor_default_helper1 (noBins, hashFunc, kmerSize, bits);
    KmerFilter<Dna, TSpec> ctor_default_helper2 (noBins, hashFunc, kmerSize, bits);

    // KmerFilter<Dna, TSpec> ctor_default (noBins, kmerSize);
    // KmerFilter<Dna, TSpec> ctor_default_helper1 (noBins, kmerSize);
    // KmerFilter<Dna, TSpec> ctor_default_helper2 (noBins, kmerSize);

    // Copy constructor
    KmerFilter<Dna, TSpec> ctor_copy (ctor_default);
    // Copy assignment
    KmerFilter<Dna, TSpec> assignment_copy;
    assignment_copy = ctor_default;
    // Move constructor
    KmerFilter<Dna, TSpec> ctor_move(std::move(ctor_default_helper1));
    // Move assignment
    KmerFilter<Dna, TSpec> assignment_move;
    assignment_move = std::move(ctor_default_helper2);


    // ==========================================================================
    // Test addKmer()
    // ==========================================================================
    std::cout << "Testing addKmer" << '\n';

        CharString fasta("../test.fasta.gz");
        addFastaFile(ctor_default, toCString(fasta), 1);
        addFastaFile(ctor_default, toCString(fasta), 5);
        addFastaFile(ctor_default, toCString(fasta), 8);
        addFastaFile(ctor_copy, toCString(fasta), 2);
        addFastaFile(assignment_copy, toCString(fasta), 3);
        addFastaFile(ctor_move, toCString(fasta), 0);
        addFastaFile(assignment_move, toCString(fasta), 9);


    // ==========================================================================
    // Test storing and retrieving
    // ==========================================================================
    std::cout << "Testing store/retrieve" << '\n';

    CharString store1("file.dat");
    store(ctor_default, toCString(store1));
    retrieve(ctor_empty, toCString(store1));
    std::remove(toCString(store1));


    // ==========================================================================
    // Test whichBins()
    // ==========================================================================
    std::cout << "Testing whichBins" << '\n';

    std::vector<uint64_t> ctor_default_set;

    std::vector<bool> which = whichBins(ctor_default, DnaString("AAA"), 1);
    (void) whichBins(ctor_default, DnaString("AAA"));
    for (uint64_t i = 0; i < which.size(); ++i)
    {
        if (i == 1 || i == 5 || i == 8)
        {
            assert(which[i]);
            ctor_default_set.push_back(i);
        }
        else
            assert(!which[i]);
    }


    // ==========================================================================
    // Test clearBins()
    // ==========================================================================
    std::cout << "Testing clearBins" << '\n';

    // Check if any elements are set in the filters.
    bool ctor_default_any = false;
    for (uint64_t i = 0; i < ctor_default.noOfBlocks * ctor_default.noOfBins; ++i)
    {
        if (ctor_default.filterVector[i])
        {
            ctor_default_any = true;
            break;
        }
    }
    assert(ctor_default_any == true);

    // Reset the filter vectors.
    clearBins(ctor_default, ctor_default_set, threads);

    // Check if filter Vectors are empty.
    ctor_default_any = false;
    for (uint64_t i = 0; i < ctor_default.noOfBlocks * ctor_default.noOfBins; ++i)
    {
        if (ctor_default.filterVector[i])
        {
            ctor_default_any = true;
            break;
        }
    }
    assert(ctor_default_any == false);

    return 0;
}
