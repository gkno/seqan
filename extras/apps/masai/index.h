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
// Author: Enrico Siragusa <enrico.siragusa@fu-berlin.de>
// ==========================================================================
// This file contains index type definitions.
// ==========================================================================

#ifndef SEQAN_EXTRAS_MASAI_INDEX_H_
#define SEQAN_EXTRAS_MASAI_INDEX_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/file.h>
#include <seqan/index.h>
#include <seqan/index_fm.h>

#include "index/index_qgram_stretched.h"
#include "index/find_backtracking_stretched.h"

using namespace seqan;


// ============================================================================
// QGram Shape Definitions
// ============================================================================

typedef UngappedShape<10>                               TShape;


// ============================================================================
// Genome Index Type Definitions
// ============================================================================

// ----------------------------------------------------------------------------
// Genome Suffix Array Value Type Definition
// ----------------------------------------------------------------------------

template <>
struct SAValue<TGenome>
{
    typedef Pair<unsigned char, unsigned int, Pack> Type;
};

// ----------------------------------------------------------------------------
// Genome Enhanced Suffix Array Type Definitions
// ----------------------------------------------------------------------------

typedef Index<TGenome, IndexEsa<> >                     TGenomeEsa;

#if defined(SEQAN_EXTRAS_MASAI_DISABLE_MMAP) || defined(SEQAN_EXTRAS_MASAI_INDEXER_H_)
typedef DefaultIndexStringSpec<TGenomeEsa>::Type        TGenomeEsaStringSpec;
#else
typedef MMap<>                                          TGenomeEsaStringSpec;
#endif

namespace seqan {
template <>
struct Fibre<TGenomeEsa, FibreLcp>
{
    typedef String<unsigned int, TGenomeEsaStringSpec>   Type;
};

template <>
struct Fibre<TGenomeEsa, FibreChildtab>
{
    typedef String<unsigned int, TGenomeEsaStringSpec>   Type;
};
}

// ----------------------------------------------------------------------------
// Genome Suffix Array Type Definitions
// ----------------------------------------------------------------------------

typedef Index<TGenome, IndexSa<> >                      TGenomeSa;

#if defined(SEQAN_EXTRAS_MASAI_DISABLE_MMAP) || defined(SEQAN_EXTRAS_MASAI_INDEXER_H_)
typedef DefaultIndexStringSpec<TGenomeSa>::Type         TGenomeSaStringSpec;
#else
typedef MMap<>                                          TGenomeSaStringSpec;
#endif

// ----------------------------------------------------------------------------
// Genome QGram Index with Bucket Refinement Type Definitions
// ----------------------------------------------------------------------------

typedef IndexQGram<TShape>                              TQGramBaseIndex;
typedef IndexQGram<TShape, BucketRefinement>            TQGramBucketRefinementIndex;

typedef Index<TGenome, TQGramBaseIndex>                 TGenomeBaseQGram;
typedef Index<TGenome, IndexSa<InfixSegment> >          TGenomeInfixSa;
typedef Index<TGenome, TQGramBucketRefinementIndex>     TGenomeQGram;

#if defined(SEQAN_EXTRAS_MASAI_DISABLE_MMAP) || defined(SEQAN_EXTRAS_MASAI_INDEXER_H_)
typedef DefaultIndexStringSpec<TGenomeBaseQGram>::Type  TGenomeQGramStringSpec;
#else
typedef MMap<>                                          TGenomeQGramStringSpec;
#endif

namespace seqan {
template <>
struct Fibre<TGenomeBaseQGram, FibreDir>
{
    typedef String<unsigned int, TGenomeQGramStringSpec>   Type;
};

template <>
struct Fibre<TGenomeInfixSa, FibreSA>
{
    typedef Segment<Fibre<TGenomeBaseQGram, FibreSA>::Type const, InfixSegment>     Type;
};
}

// ----------------------------------------------------------------------------
// Genome FM Index Type Definitions
// ----------------------------------------------------------------------------

typedef Index<TGenome, FMIndex<WT<FmiDollarSubstituted<MultiDollar<void> > >, CompressText> > TGenomeFM;


// ============================================================================
// Reads Index Type Definitions
// ============================================================================

// ----------------------------------------------------------------------------
// Reads Suffix Array Value Type Definition
// ----------------------------------------------------------------------------

template <>
struct SAValue<TReadSeqStore>
{
    typedef Pair<unsigned int, unsigned short, Pack> Type;
};

// ----------------------------------------------------------------------------
// Reads Wotd Type Definitions
// ----------------------------------------------------------------------------

typedef Index<TReadSeqStore, IndexWotd<> >              TReadsWotd;

//namespace seqan {
//	template <>
//	struct Fibre<TReadsWotd, FibreDir>
//	{
//		typedef String< unsigned int, DefaultIndexStringSpec<TReadsWotd>::Type >   Type;
//	};
//}

// ----------------------------------------------------------------------------
// Reads QGram Index with Bucket Refinement Type Definitions
// ----------------------------------------------------------------------------

//typedef Index<TReadSeqStore, TQGramBaseIndex>                 TReadsBaseQGram;
//typedef Index<TReadSeqStore, IndexSa<InfixSegment> >          TReadsInfixSa;
//typedef Index<TReadSeqStore, TQGramBucketRefinementIndex>     TReadsQGram;
//
//namespace seqan
//{
//    template <>
//    struct Fibre<TReadsBaseQGram, FibreDir>
//    {
//        typedef String<unsigned int, DefaultIndexStringSpec<TReadsBaseQGram>::Type >   Type;
//    };
//
//    template <>
//    struct Fibre<TReadsInfixSa, FibreSA>
//    {
//        typedef Segment<Fibre<TReadsBaseQGram, FibreSA>::Type const, InfixSegment>           Type;
//    };
//}

// ----------------------------------------------------------------------------
// Reads QGram Index Type Definitions
// ----------------------------------------------------------------------------

typedef Index<TReadSeqStore, TQGramBaseIndex>               TReadsQGram;

namespace seqan {
template <>
struct Fibre<TReadsQGram, FibreDir>
{
    typedef String<unsigned int, DefaultIndexStringSpec<TReadsQGram>::Type>    Type;
};
}

// ----------------------------------------------------------------------------
// Reads QGram Index with Wotd Subtree Type Definitions
// ----------------------------------------------------------------------------

//typedef Index<TReadSeqStore, IndexWotd< Subtree<> > >   TReadsWotdSubtree;
//
//namespace seqan
//{
//    template <>
//    struct Fibre<TReadsWotdSubtree, FibreSA>
//    {
//        typedef Segment< Fibre<TReadsQGram, FibreSA>::Type, InfixSegment>  Type;
//    };
//}


#endif  // #ifndef SEQAN_EXTRAS_MASAI_INDEX_H_
