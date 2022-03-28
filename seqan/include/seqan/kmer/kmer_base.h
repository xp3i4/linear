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
// Author:  Temesgen H. Dadi <temesgen.dadi@fu-berlin.de>
//          Enrico Seiler <enrico.seiler@fu-berlin.de>
// ==========================================================================

#ifndef INCLUDE_SEQAN_KMER_KMER_BASE_H_
#define INCLUDE_SEQAN_KMER_KMER_BASE_H_

#include <sdsl/bit_vectors.hpp>
#include <seqan/seq_io.h>
#include <valarray>
#include <algorithm>
#include <future>

namespace seqan {

// ==========================================================================
// Forwards
// ==========================================================================

// ==========================================================================
// Tags, Classes, Enums
// ==========================================================================

// --------------------------------------------------------------------------
// Tag k-mer Filter Tags
// --------------------------------------------------------------------------

//!\brief A tag for the IBF.
struct InterleavedBloomFilter_;
typedef Tag<InterleavedBloomFilter_> InterleavedBloomFilter;

//!\brief A tag for direct addressing.
struct DirectAddressing_;
typedef Tag<DirectAddressing_> DirectAddressing;

// --------------------------------------------------------------------------
// Class KmerFilter
// --------------------------------------------------------------------------

//!\brief The KmerFilter class.
template<typename TValue = Dna, typename TSpec = DirectAddressing>
class KmerFilter;

// ==========================================================================
// Metafunctions
// ==========================================================================

//!\brief Type definition for variables.
template<typename TValue, typename TSpec>
struct Value<KmerFilter<TValue, TSpec> >
{
    typedef uint64_t Type;
};

// --------------------------------------------------------------------------
// Metafunction MetafunctionName
// --------------------------------------------------------------------------

// ==========================================================================
// Functions
// ==========================================================================

// --------------------------------------------------------------------------
// Function functionName()
// --------------------------------------------------------------------------

/*!
 * \brief Adds a k-mer to a bin in a given filter.
 * \param me The KmerFilter instance.
 * \param kmer The k-mer to be added.
 * \param binNo The bin to add the k-mer to.
 */
template<typename TValue, typename TSpec, typename TString, typename TInt>
inline void addKmer(KmerFilter<TValue, TSpec> & me, TString const & kmer, TInt && binNo)
{
    me.addKmer(kmer, binNo);
}

/*!
 * \brief Sets the vectors for given bins to 0.
 * \param me The KmerFilter instance.
 * \param bins A vector containing the bin numbers.
 * \param threads The number of threads to use.
 */
template<typename TValue, typename TSpec, typename TInt1, typename TInt2>
inline void clearBins(KmerFilter<TValue, TSpec> &  me, std::vector<TInt1> & bins, TInt2&& threads)
{
    me.clearBins(bins, static_cast<uint64_t>(threads));
}

/*!
 * \brief Adds all k-mers from a fasta file to a bin of a given KmerFilter.
 * \param me The KmerFilter instance.
 * \param kmer The fasta file to process.
 * \param binNo The bin to add the k-mers to.
 */
template<typename TValue, typename TSpec, typename TInt>
inline void addFastaFile(KmerFilter<TValue, TSpec> &  me, const char * fastaFile, TInt && binNo)
{
    CharString id;
    String<TValue> seq;

    SeqFileIn seqFileIn;
    if (!open(seqFileIn, fastaFile))
    {
        CharString msg = "Unable to open contigs file: ";
        append(msg, CharString(fastaFile));
        throw toCString(msg);
    }
    while(!atEnd(seqFileIn))
    {
        readRecord(id, seq, seqFileIn);
        if(length(seq) < me.kmerSize)
            continue;
        addKmer(me, seq, binNo);
    }
    close(seqFileIn);
}

/*!
 * \brief Calculates the k-mer counts of a given text.
 * \param me The KmerFilter instance.
 * \param counts Vector of length binNo to save counts to.
 * \param text A single text to count all contained k-mers for.
 */
template<typename TValue, typename TSpec>
inline void whichBins(KmerFilter<TValue, TSpec> &  me, std::vector<uint64_t> & counts, String<TValue> const & text)
{
    me.whichBins(counts, text);
}

/*!
 * \brief Returns the k-mer counts of a given text.
 * \param me The KmerFilter instance.
 * \param text A single text to count all contained k-mers for.
 * \returns std::vector<uint64_t> of size binNo containing counts.
 */
template<typename TValue, typename TSpec>
inline std::vector<uint64_t> whichBins(KmerFilter<TValue, TSpec> &  me, String<TValue> const & text)
{
    std::vector<uint64_t> counts(me.noOfBins, 0);
    whichBins(me, counts, text);
    return counts;
}

/*!
 * \brief Checks for which bins the counts of all k-mers in a text exceed a threshold.
 * \param me The KmerFilter instance.
 * \param selected Vector of length binNo to save true/false to.
 * \param text A single text to count all contained k-mers for.
 * \param threshold The minimal number of occurences to return true for the bin.
 */
template<typename TValue, typename TSpec, typename TInt>
inline void whichBins(KmerFilter<TValue, TSpec> &  me, std::vector<bool> & selected, String<TValue> const & text, TInt threshold)
{
    me.whichBins(selected, text, static_cast<uint64_t>(threshold));
}

/*!
 * \brief Returns for which bins the counts of all k-mers in a text exceed a threshold.
 * \param me The KmerFilter instance.
 * \param text A single text to count all contained k-mers for.
 * \param threshold The minimal number of occurences to return true for the bin.
 * \returns std::vector<bool> of size binNo indicating whether the text is in the bin.
 */
template<typename TValue, typename TSpec, typename TInt>
inline std::vector<bool> whichBins(KmerFilter<TValue, TSpec> &  me, String<TValue> const & text, TInt && threshold)
{
    std::vector<bool> selected(me.noOfBins, false);
    whichBins(me, selected, text, threshold);
    return selected;
}

/*!
 * \brief Returns the number of bins.
 * \param me The KmerFilter instance.
 * \returns Value<KmerFilter<TValue, TSpec> >::Type Number of bins.
 */
template<typename TValue, typename TSpec>
inline typename Value<KmerFilter<TValue, TSpec> >::Type getNumberOfBins(KmerFilter<TValue, TSpec> &  me)
{
    return me.noBins;
}

/*!
 * \brief Returns the k-mer size.
 * \param me The KmerFilter instance.
 * \returns Value<KmerFilter<TValue, TSpec> >::Type k-mer size.
 */
template<typename TValue, typename TSpec>
inline typename Value<KmerFilter<TValue, TSpec> >::Type getKmerSize(KmerFilter<TValue, TSpec> &  me)
{
    return me.kmerSize;
}

/*!
 * \brief Reads the metadata.
 * \param me The KmerFilter instance.
 */
template<typename TValue, typename TSpec>
inline void getMetadata(KmerFilter<TValue, TSpec> &  me)
{
    typedef typename Value<KmerFilter<TValue, TSpec> >::Type THValue;

    //-------------------------------------------------------------------
    //| kmer_size | n_hash_func | n_bins |              bf              |
    //-------------------------------------------------------------------
    me.noOfBits = me.filterVector.bit_size();

    THValue metadataStart = me.noOfBits - me.filterMetadataSize;
    me.noOfBins = me.filterVector.get_int(metadataStart);
    me.noOfHashFunc = me.filterVector.get_int(metadataStart+64);
    me.kmerSize = me.filterVector.get_int(metadataStart+128);
}

/*!
 * \brief Writes the metadata.
 * \param me The KmerFilter instance.
 */
template<typename TValue, typename TSpec>
inline void setMetadata(KmerFilter<TValue, TSpec> &  me)
{
    typedef typename Value<KmerFilter<TValue, TSpec> >::Type THValue;

    //-------------------------------------------------------------------
    //| kmer_size | n_hash_func | n_bins |              bf              |
    //-------------------------------------------------------------------

    THValue metadataStart = me.noOfBits - me.filterMetadataSize;

    // TODO also store TValue (alphabet)
    me.filterVector.set_int(metadataStart, me.noOfBins);
    me.filterVector.set_int(metadataStart + 64, me.noOfHashFunc);
    me.filterVector.set_int(metadataStart+128, me.kmerSize);
}

/*!
 * \brief Returns the of the filter vector in MB.
 * \param me The KmerFilter instance.
 * \returns double filter vector size in MB.
 */
template<typename TValue, typename TSpec, typename THValue>
inline double size(KmerFilter<TValue, TSpec> &  me)
{
    return sdsl::size_in_mega_bytes(me.filterVector);
}

/*!
 * \brief Writes the filter vector to a file.
 * \param me The KmerFilter instance.
 * \param fileName Name of the file to write to.
 * \returns bool Indicates if the operation was successful.
 */
template<typename TValue, typename TSpec>
inline bool store(KmerFilter<TValue, TSpec> &  me, const char * fileName)
{
    setMetadata(me);
    return sdsl::store_to_file(me.filterVector, fileName);
}

/*!
 * \brief Loads the filter vector from a file.
 * \param me The KmerFilter instance.
 * \param fileName Name of the file to read from.
 * \returns bool Indicates if the operation was successful.
 */
template<typename TValue, typename TSpec>
inline bool retrieve(KmerFilter<TValue, TSpec> &  me, const char * fileName)
{
    if (!sdsl::load_from_file(me.filterVector, fileName))
    {
        std::cerr << "File \"" << fileName << "\" could not be read." << std::endl;
        exit(1);
    }
    getMetadata(me);
    me.init();
    return true;
}

/*!
 * \brief Indicates whether a bit is set in an integer.
 * \param me The KmerFilter instance.
 * \param num Integer to check.
 * \param bit bit position to check.
 * \returns bool Indicates if the bit is set.
 */
template<typename TValue, typename TSpec, typename TInt1, typename TInt2>
inline bool isBitSet(KmerFilter<TValue, TSpec> &  me, TInt1&& num, TInt2&& bit)
{
    return 1 == ( (num >> bit) & 1);
}

}  // namespace seqan

#endif  // INCLUDE_SEQAN_KMER_KMER_BASE_H_
