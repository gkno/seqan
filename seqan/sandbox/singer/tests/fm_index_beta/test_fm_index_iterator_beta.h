// ==========================================================================
//                               fm_index_beta
// ==========================================================================
// Copyright (c) 2006-2011, Knut Reinert, FU Berlin
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
// Author: Your Name <your.email@example.net>
// ==========================================================================

#ifndef TEST_FM_INDEX_ITERATOR_BETA_H_
#define TEST_FM_INDEX_ITERATOR_BETA_H_

#include <seqan/basic.h>
#include <seqan/index.h>
#include <seqan/sequence.h>
#include <seqan/fm_index.h>

#include <seqan/random.h>

using namespace seqan;

template <typename TIter>
void fmIndexIteratorConstuctor(TIter & /*tag*/)
{
	typedef typename Container<TIter>::Type TIndex;
    typedef typename Value<TIndex>::Type TText;

	TText text;
	generateText(text);

	Index<TText> esa(text);
	typename Iterator<Index<TText>, TopDown<> >::Type esaIt(esa);
	goDown(esaIt);

	TIndex fmIndex(text);

	TIter it(fmIndex);

	SEQAN_ASSERT_EQ(isRoot(it), true);
	SEQAN_ASSERT_EQ(repLength(it), 0u);
	SEQAN_ASSERT_EQ(goDown(it), true);
	SEQAN_ASSERT_EQ(isRoot(it), false);
	SEQAN_ASSERT_EQ(repLength(it), 1u);
	SEQAN_ASSERT_EQ(representative(it), 'A');
	SEQAN_ASSERT_EQ(goRight(it), true);
	SEQAN_ASSERT_EQ(representative(it), 'C');
	SEQAN_ASSERT_EQ(repLength(it), 1u);
	SEQAN_ASSERT_EQ(goRight(it), true);
	SEQAN_ASSERT_EQ(representative(it), 'G');
	SEQAN_ASSERT_EQ(repLength(it), 1u);
	SEQAN_ASSERT_EQ(goRight(it), true);
	SEQAN_ASSERT_EQ(representative(it), 'T');
	SEQAN_ASSERT_EQ(repLength(it), 1u);
	SEQAN_ASSERT_EQ(goRight(it), false);
	SEQAN_ASSERT_EQ(representative(it), 'T');
	SEQAN_ASSERT_EQ(repLength(it), 1u);
	SEQAN_ASSERT_EQ(goDown(it), true);
	SEQAN_ASSERT_EQ(goDown(it), true);
	SEQAN_ASSERT_EQ(goUp(it), true);
	SEQAN_ASSERT_EQ(representative(it), "AT");
	SEQAN_ASSERT_EQ(goUp(it), true);
	SEQAN_ASSERT_EQ(representative(it), "T");
	
}

template <typename TIter>
void fmIndexIteratorGoDown(TIter & /*tag*/)
{
    typedef typename Container<TIter>::Type TIndex;
    typedef typename Value<TIndex>::Type TText;

	TText text = "ACGACGACG";
	TIndex fmIndex(text);
	TIter it(fmIndex);

	SEQAN_ASSERT_EQ(goDown(it), true);
	SEQAN_ASSERT_EQ(representative(it), 'A');

	SEQAN_ASSERT_EQ(goDown(it, 'G'), true);
	SEQAN_ASSERT_EQ(representative(it), "AG");
	
	SEQAN_ASSERT_EQ(goDown(it, "ACGAC"), true);
	SEQAN_ASSERT_EQ(representative(it), "AGCAGCA");

	SEQAN_ASSERT_EQ(goDown(it), false);
	SEQAN_ASSERT_EQ(representative(it), "AGCAGCA");

	SEQAN_ASSERT_EQ(goDown(it, 'G'), false);
	SEQAN_ASSERT_EQ(representative(it), "AGCAGCA");
	
	SEQAN_ASSERT_EQ(goDown(it, "ACGAC"), false);
	SEQAN_ASSERT_EQ(representative(it), "AGCAGCA");
}

SEQAN_DEFINE_TEST(fm_index_iterator_constuctor)
{
    using namespace seqan;

    typedef FmIndex<WaveletTreeBased<FmiDollarSubstituted<SingleDollar<void> > >, void> TDefaultIndex;
    typedef FmIndex<WaveletTreeBased<FmiDollarSubstituted<SingleDollar<void> > >, Compressed> TCompressedIndex;
    typedef TopDown<> TIterSpec;
    typedef TopDown<ParentLinks<> > TParentLinksIterSpec;

    DnaString genome = "AAA";
    
//     {
//         Index<DnaString,TDefaultIndex> index(genome);
//         Iterator<Index<DnaString,TDefaultIndex>, TIterSpec>::Type dnaTag(index);
//         fmIndexIteratorConstuctor(dnaTag);
//     }
//     {
//         Index<DnaString,TCompressedIndex> index(genome);
//         Iterator<Index<DnaString,TCompressedIndex>, TIterSpec>::Type dnaTag(index);
//         fmIndexIteratorConstuctor(dnaTag);
//     }
    {
        Index<DnaString,TDefaultIndex> index(genome);
        Iterator<Index<DnaString,TDefaultIndex>, TParentLinksIterSpec>::Type dnaTag(index);
        fmIndexIteratorConstuctor(dnaTag);
    }
    {
        Index<DnaString,TCompressedIndex> index(genome);
        Iterator<Index<DnaString,TCompressedIndex>, TParentLinksIterSpec>::Type dnaTag(index);
        fmIndexIteratorConstuctor(dnaTag);
    }
}

#endif // TEST_FM_INDEX_ITERATOR_BETA_H_

