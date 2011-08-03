// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2010, Knut Reinert, FU Berlin
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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================

#ifndef EXTRAS_INCLUDE_SEQAN_BAM_IO_BAM_INDEX_BAI_H_
#define EXTRAS_INCLUDE_SEQAN_BAM_IO_BAM_INDEX_BAI_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

struct Bai_;
typedef Tag<Bai_> Bai;

struct BaiBamIndexBinData_
{
    String<Pair<__uint64, __uint64> > chunkBegEnds;
};

/**
.Spec.BAI BamIndex
..cat:BAM I/O
..general:Class.BamIndex
..summary:Access to BAI (samtools-style) Indices.
..signature:BamIndex<Bai>
..include:seqan/bam_io.h
*/

template <>
class BamIndex<Bai>
{
public:
    typedef std::map<__uint32, BaiBamIndexBinData_> TBinIndex_;
    typedef String<__uint64> TLinearIndex_;

    __uint64 _unalignedCount;
    
    String<TBinIndex_> _binIndices;
    String<TLinearIndex_> _linearIndices;
    
    BamIndex() :
            _unalignedCount(0)
    {}
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function jumpToPos()
// ----------------------------------------------------------------------------

/**
.Function.jumpToPos
..cat:BAM I/O
..signature:jumpToPos(bgzfStream, hasAlignments, bamIOContext, refId, pos, bamIndex)
..summary:Seek in BAM BGZF stream using an index.
..param.bgzfStream:The BGZF Stream to seek in.
...type:Spec.BGZF Stream
..param.refId:Reference ID to seek to.
...type:nolink:$__int32$
..param.hasAlignments:Set to $true$ iff there are alignments at this position.
...type:nolink:$bool$
..param.bamIOContext:Context to use for loading alignments.
...type:Class.BamIOContext
..param.pos:Zero-based position in the reference.
...type:nolink:$__int32$
..param.bamIndex:The index to use.
...type:Class.BamIndex
..returns:$bool$ indicating success.
..remarks:This function may fail if the refId/pos is invalid.
..remarks:This function jumps to a "close position left of $pos$".
..include:seqan/bam_io.h
*/

static inline void
_baiReg2bins(String<__uint16> & list, __uint32 beg, __uint32 end)
{
	unsigned k;
	if (beg >= end) return;
	if (end >= 1u<<29) end = 1u<<29;
	--end;
	appendValue(list, 0);
	for (k =    1 + (beg>>26); k <=    1 + (end>>26); ++k) appendValue(list, k);
	for (k =    9 + (beg>>23); k <=    9 + (end>>23); ++k) appendValue(list, k);
	for (k =   73 + (beg>>20); k <=   73 + (end>>20); ++k) appendValue(list, k);
	for (k =  585 + (beg>>17); k <=  585 + (end>>17); ++k) appendValue(list, k);
    for (k = 4681 + (beg>>14); k <= 4681 + (end>>14); ++k) appendValue(list, k);
}

template <typename TNameStore, typename TNameStoreCache>
inline bool
jumpToPos(Stream<Bgzf> & stream, bool & hasAlignments, BamIOContext<TNameStore, TNameStoreCache> /*const*/ & bamIOContext, __int32 refId, __int32 pos, BamIndex<Bai> const & index)
{
    hasAlignments = false;
    if (refId < 0)
        return false;  // Cannot seek to invalid reference.
    if (static_cast<unsigned>(refId) >= length(index._binIndices))
        return false;  // Cannot seek to invalid reference.
    
    // ------------------------------------------------------------------------
    // Compute offset in BGZF file.
    // ------------------------------------------------------------------------
    __uint64 offset = MaxValue<__uint64>::VALUE;

    // Retrieve the candidate bin identifiers for [pos, pos+1).
    String<__uint16> candidateBins;
    _baiReg2bins(candidateBins, pos, pos + 1);

    // Retrieve the smallest required offset from the linear index.
    unsigned windowIdx = pos >> 14;  // Linear index consists of 16kb windows.
    __uint64 linearMinOffset = 0;
    if (windowIdx > length(index._linearIndices[refId]))
        linearMinOffset = back(index._linearIndices[refId]);
    else
        linearMinOffset = index._linearIndices[refId][windowIdx];

    // Combine candidate bins and smallest required offset from linear index into candidate offset.
    typedef std::set<__uint64> TOffsetCandidates;
    TOffsetCandidates offsetCandidates;
    typedef typename Iterator<String<__uint16>, Rooted>::Type TCandidateIter;
    for (TCandidateIter it = begin(candidateBins, Rooted()); !atEnd(it); goNext(it))
    {
        typedef typename std::map<__uint32, BaiBamIndexBinData_>::const_iterator TMapIter;
        TMapIter mIt = index._binIndices[refId].find(*it);
        if (mIt == index._binIndices[refId].end())
            continue;  // Candidate is not in index!

        typedef typename Iterator<String<Pair<__uint64, __uint64> > const, Rooted>::Type TBegEndIter;
        for (TBegEndIter it2 = begin(mIt->second.chunkBegEnds, Rooted()); !atEnd(it2); goNext(it2))
            if (it2->i2 >= linearMinOffset)
                offsetCandidates.insert(it2->i1);
    }

    // Search through candidate offsets, find smallest with a fitting alignment.
    typedef typename TOffsetCandidates::const_iterator TOffsetCandidateIter;
    BamAlignmentRecord record;
    for (TOffsetCandidateIter candIt = offsetCandidates.begin(); candIt != offsetCandidates.end(); ++candIt)
    {
        int res = streamSeek(stream, *candIt, SEEK_SET);
        if (res != 0)
            return false;  // Error while seeking.
        res = readRecord(record, bamIOContext, stream, Bam());
        if (res != 0)
            return false;  // Error while reading.
        __int32 endPos = record.pos + getAlignmentLengthInRef(record);
        if ((record.rId == refId && endPos >= pos) || (record.rId > refId))  // Found!
        {
            hasAlignments = true;
            if (candIt != offsetCandidates.begin())
                --candIt;  // Against overlaps, logic taken from BamTools.
            offset = *candIt;
            break;
        }
    }

    if (offset != MaxValue<__uint64>::VALUE)
    {
        int res = streamSeek(stream, offset, SEEK_SET);
        if (res != 0)
            return false;  // Error while seeking.
    }
    // Finding no overlapping alignment is not an error, hasAlignments is false.
    return true;
}

// ----------------------------------------------------------------------------
// Function getUnalignedCount()
// ----------------------------------------------------------------------------

/**
.Function.getUnalignedCount
..cat:BAM I/O
..signature:load(index, filename)
..summary:Query index for number of unaligned reads.
..param.index:Index to query.
...type:Class.BamIndex
..returns:$__uint64$ with number of unaligned reads.
..include:seqan/bam_io.h
*/

inline __uint64
getUnalignedCount(BamIndex<Bai> const & index)
{
    return index._unalignedCount;
}

// ----------------------------------------------------------------------------
// Function load()
// ----------------------------------------------------------------------------

/**
.Function.load
..cat:BAM I/O
..signature:load(index, filename)
..summary:Load a BAM index from a given file name.
..param.index:Target data structure.
...type:Class.BamIndex
..param.filename:Path to file to load.
...type:nolink:$char const *$
..returns:$bool$ indicating success.
..include:seqan/bam_io.h
 */

inline bool
load(BamIndex<Bai> & index, char const * filename)
{
    std::fstream fin(filename, std::ios::binary | std::ios::in);
    if (!fin.good())
        return false;  // Could not open file.

    // Read magic number.
    CharString buffer;
    resize(buffer, 4);
    fin.read(&buffer[0], 4);
    if (!fin.good())
        return false;
    if (buffer != "BAI\1")
        return false;  // Magic number is wrong.

    __int32 nRef = 0;
    fin.read(reinterpret_cast<char *>(&nRef), 4);
    if (!fin.good())
        return false;

    resize(index._linearIndices, nRef);
    resize(index._binIndices, nRef);
    
    for (int i = 0; i < nRef; ++i)  // For each reference.
    {
        // Read bin index.
        __int32 nBin = 0;
        fin.read(reinterpret_cast<char *>(&nBin), 4);
        if (!fin.good())
            return false;
        index._binIndices[i].clear();
        BaiBamIndexBinData_ data;
        for (int j = 0; j < nBin; ++j)  // For each bin.
        {
            clear(data.chunkBegEnds);

            __uint32 bin = 0;
            fin.read(reinterpret_cast<char *>(&bin), 4);
            if (!fin.good())
                return false;

            __int32 nChunk = 0;
            fin.read(reinterpret_cast<char *>(&nChunk), 4);
            if (!fin.good())
                return false;
            reserve(data.chunkBegEnds, nChunk);
            for (int k = 0; k < nChunk; ++k)  // For each chunk;
            {
                __uint64 chunkBeg = 0;
                __uint64 chunkEnd = 0;
                fin.read(reinterpret_cast<char *>(&chunkBeg), 8);
                fin.read(reinterpret_cast<char *>(&chunkEnd), 8);
                if (!fin.good())
                    return false;
                appendValue(data.chunkBegEnds, Pair<__uint64>(chunkBeg, chunkEnd));
            }

            // Copy bin data into index.
            index._binIndices[i][bin] = data;
        }

        // Read linear index.
        __int32 nIntv = 0;
        fin.read(reinterpret_cast<char *>(&nIntv), 4);
        if (!fin.good())
            return false;
        clear(index._linearIndices[i]);
        reserve(index._linearIndices[i], nIntv);
        for (int j = 0; j < nIntv; ++j)
        {
            __uint64 ioffset = 0;
            fin.read(reinterpret_cast<char *>(&ioffset), 8);
            if (!fin.good())
                return false;
            appendValue(index._linearIndices[i], ioffset);
        }
    }

    if (!fin.good())
        return false;

    // Read (optional) number of alignments without coordinate.
    __uint64 nNoCoord = 0;
    fin.read(reinterpret_cast<char *>(&nNoCoord), 8);
    if (!fin.good())
    {
        fin.clear();
        nNoCoord = 0;
    }
    index._unalignedCount = nNoCoord;

    return true;
}

}  // namespace seqan

#endif  // #ifndef EXTRAS_INCLUDE_SEQAN_BAM_IO_BAM_INDEX_BAI_H_
