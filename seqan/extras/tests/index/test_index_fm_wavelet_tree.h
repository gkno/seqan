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

#ifndef TESTS_WAVELT_TREE_BETA_H_
#define TESTS_WAVELT_TREE_BETA_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/fm_sequence.h>
#include "test_fm_index_beta.h"

using namespace seqan;

template <typename TWaveletTree>
void waveletTreeClear(TWaveletTree & /*tag*/)
{
	typedef typename Fibre<TWaveletTree, FibreWaveletTreeStructure>::Type TWaveletTreeStructure;
	typedef typename Fibre<TWaveletTreeStructure, FibreTreeVertieces>::Type TWaveletTreeVertieces;
	typedef typename Value<TWaveletTreeVertieces>::Type TWaveletTreeVertex;
	typedef typename Value<TWaveletTree>::Type TChar;
	typedef typename Value<TWaveletTreeVertex, 2>::Type TPos;
	
	String<TChar> text = "ACGTNACGTNACGTN";
	TWaveletTree waveletTree(text);

    clear(waveletTree);

	SEQAN_ASSERT_EQ(empty(waveletTree), true);
}

template <typename TWaveletTree>
void waveletTreeConstructor(TWaveletTree & /*tag*/)
{
	typedef typename Fibre<TWaveletTree, FibreWaveletTreeStructure>::Type TWaveletTreeStructure;
	typedef typename Fibre<TWaveletTreeStructure, FibreTreeVertieces>::Type TWaveletTreeVertieces;
	typedef typename Value<TWaveletTreeVertieces>::Type TWaveletTreeVertex;
	typedef typename Value<TWaveletTree>::Type TChar;
	typedef typename Value<TWaveletTreeVertex, 2>::Type TPos;

	{
		TWaveletTree waveletTree;
	}
	{
		String<TChar> text;        
        generateText(text);
		TWaveletTree waveletTree(text);

		for (unsigned i = 0; i < length(text); ++i)
        {
		    SEQAN_ASSERT_EQ(getCharacter(waveletTree, i), (TChar)text[i]);
        }
    }
}

template <typename TWaveletTree>
void waveletTreeDollarPosition(TWaveletTree & /*tag*/)
{
	typedef typename Fibre<TWaveletTree, FibreWaveletTreeStructure>::Type TWaveletTreeStructure;
	typedef typename Fibre<TWaveletTreeStructure, FibreTreeVertieces>::Type TWaveletTreeVertieces;
	typedef typename Value<TWaveletTreeVertieces>::Type TWaveletTreeVertex;
	typedef typename Value<TWaveletTree>::Type TChar;
	typedef typename Value<TWaveletTreeVertex, 2>::Type TPos;
	typedef typename Fibre<TWaveletTree, FibreDollarPosition>::Type TDollarPos;

	TWaveletTree waveletTree;

	TDollarPos dollarPos = 0u;

	setDollarPosition(waveletTree, dollarPos);
	SEQAN_ASSERT_EQ(getDollarPosition(waveletTree), dollarPos);
}

template <typename TWaveletTree>
void waveletTreeDollarSubstitute(TWaveletTree & /*tag*/)
{
	typedef typename Fibre<TWaveletTree, FibreWaveletTreeStructure>::Type TWaveletTreeStructure;
	typedef typename Fibre<TWaveletTreeStructure, FibreTreeVertieces>::Type TWaveletTreeVertieces;
	typedef typename Value<TWaveletTreeVertieces>::Type TWaveletTreeVertex;
	typedef typename Value<TWaveletTree>::Type TChar;
	typedef typename Value<TWaveletTreeVertex, 2>::Type TPos;

	TWaveletTree waveletTree;

	setDollarSubstitute(waveletTree, 'C');
	SEQAN_ASSERT_EQ(getDollarSubstitute(waveletTree), 'C');
}

template <typename TWaveletTree>
void waveletTreeEmpty(TWaveletTree & /*tag*/)
{
	typedef typename Fibre<TWaveletTree, FibreWaveletTreeStructure>::Type TWaveletTreeStructure;
	typedef typename Fibre<TWaveletTreeStructure, FibreTreeVertieces>::Type TWaveletTreeVertieces;
	typedef typename Value<TWaveletTreeVertieces>::Type TWaveletTreeVertex;
	typedef typename Value<TWaveletTree>::Type TChar;
	typedef typename Value<TWaveletTreeVertex, 2>::Type TPos;
	
	String<TChar> text = "ACGTNACGTNACGTN";
	TWaveletTree waveletTree(text);
	
	SEQAN_ASSERT_EQ(empty(waveletTree), false);

    clear(waveletTree);

	SEQAN_ASSERT_EQ(length(getFibre(waveletTree, FibreBitStrings())), 0u);
	SEQAN_ASSERT_EQ(length(getFibre(waveletTree, FibreWaveletTreeStructure())), 0u);
}

template <typename TWaveletTree>
void waveletTreeGetCharacter(TWaveletTree & /*tag*/)
{
	typedef typename Fibre<TWaveletTree, FibreWaveletTreeStructure>::Type TWaveletTreeStructure;
	typedef typename Fibre<TWaveletTreeStructure, FibreTreeVertieces>::Type TWaveletTreeVertieces;
	typedef typename Value<TWaveletTreeVertieces>::Type TWaveletTreeVertex;
	typedef typename Value<TWaveletTree>::Type TChar;
	typedef typename Value<TWaveletTreeVertex, 2>::Type TPos;

	{
		String<TChar> text;
 		generateText(text);

		TWaveletTree waveletTree(text); 

        for (unsigned i = 0; i < length(text); ++i)
        {
		    SEQAN_ASSERT_EQ(getCharacter(waveletTree, i), text[i]);
        }
    }
}

template <typename TWaveletTree>
void waveletTreeGetFibre(TWaveletTree & /*tag*/)
{
	typedef typename Fibre<TWaveletTree, FibreWaveletTreeStructure>::Type TWaveletTreeStructure;
	typedef typename Fibre<TWaveletTreeStructure, FibreTreeVertieces>::Type TWaveletTreeVertieces;
	typedef typename Value<TWaveletTreeVertieces>::Type TWaveletTreeVertex;
	typedef typename Value<TWaveletTree>::Type TChar;
	typedef typename Value<TWaveletTreeVertex, 2>::Type TPos;
	
	String<TChar> text = "ACGTNACGTNACGTN";
	TWaveletTree waveletTree(text);

    typename Fibre<TWaveletTree, FibreBitStrings>::Type & tempBitStrings = getFibre(waveletTree, FibreBitStrings());
    typename Fibre<TWaveletTree, FibreWaveletTreeStructure>::Type & tempWaveletTreeStructure = getFibre(waveletTree, FibreWaveletTreeStructure());

    resize(tempBitStrings, 110);
    resize(tempWaveletTreeStructure, 100);

	SEQAN_ASSERT_EQ(length(getFibre(waveletTree, FibreBitStrings())), 110u);
	SEQAN_ASSERT_EQ(length(getFibre(waveletTree, FibreWaveletTreeStructure())), 100u);
}



template <typename TWaveletTree>
void waveletTreeGetOcc(TWaveletTree & /*tag*/)
{
	typedef typename Fibre<TWaveletTree, FibreWaveletTreeStructure>::Type TWaveletTreeStructure;
	typedef typename Fibre<TWaveletTreeStructure, FibreTreeVertieces>::Type TWaveletTreeVertieces;
	typedef typename Value<TWaveletTreeVertieces>::Type TWaveletTreeVertex;
	typedef typename Value<TWaveletTree>::Type TChar;
	typedef typename Value<TWaveletTreeVertex, 2>::Type TPos;

	{
		String<TChar> text;
 		generateText(text);
 		resize(text, 1000);

		TWaveletTree waveletTree(text);

        for (int i = MinValue<TChar>::VALUE; i <= MaxValue<TChar>::VALUE; ++i)
        {
            unsigned counter = 0;
            for (unsigned j = 0; j < length(text); ++j)
            {
                if(text[j] == (TChar)i)
                    ++counter;
		        SEQAN_ASSERT_EQ(getOccurrences(waveletTree, (TChar)i, j), counter);
            }
        }
    }
}
template <typename TString, typename TSpec>
void waveletTreeGetOcc(WaveletTree<TString, FmiDollarSubstituted<TSpec> > & /*tag*/)
{
    typedef WaveletTree<TString, FmiDollarSubstituted<TSpec> > TWaveletTree;
	typedef typename Fibre<TWaveletTree, FibreWaveletTreeStructure>::Type TWaveletTreeStructure;
	typedef typename Fibre<TWaveletTreeStructure, FibreTreeVertieces>::Type TWaveletTreeVertieces;
	typedef typename Value<TWaveletTreeVertieces>::Type TWaveletTreeVertex;
	typedef typename Value<TWaveletTree>::Type TChar;
	typedef typename Value<TWaveletTreeVertex, 2>::Type TPos;

	{
		String<TChar> text;
 		generateText(text);
 		resize(text, 1000);

		TWaveletTree waveletTree(text);
		setDollarSubstitute(waveletTree, getCharacter(waveletTree, 0u));
		setDollarPosition(waveletTree, 0u);

        for (int i = MinValue<TChar>::VALUE; i <= MaxValue<TChar>::VALUE; ++i)
        {
            unsigned counter = 0;
            for (unsigned j = 0; j < length(text); ++j)
            {
                if(text[j] == (TChar)i)
                {
                    ++counter;
                    if(text[j] == getDollarSubstitute(waveletTree) && j == getDollarPosition(waveletTree))
                        --counter;
                }
		        SEQAN_ASSERT_EQ(getOccurrences(waveletTree, (TChar)i, j), counter);
            }
        }
    }
}


template <typename TWaveletTree>
void waveletTreeNumVertieces(TWaveletTree & /*tag*/)
{
	typedef typename Fibre<TWaveletTree, FibreWaveletTreeStructure>::Type TWaveletTreeStructure;
	typedef typename Fibre<TWaveletTreeStructure, FibreTreeVertieces>::Type TWaveletTreeVertieces;
	typedef typename Value<TWaveletTreeVertieces>::Type TWaveletTreeVertex;
	typedef typename Value<TWaveletTree>::Type TChar;
	typedef typename Value<TWaveletTreeVertex, 2>::Type TPos;

    String<TChar> text;
    generateText(text);

    TWaveletTree waveletTree(text);

    SEQAN_ASSERT_EQ(numVertieces(waveletTree), (unsigned)ValueSize<TChar>::VALUE - 1);    
}

template <typename TWaveletTree>
void waveletTreeFillWaveletTree_(TWaveletTree & /*tag*/)
{
	typedef typename Fibre<TWaveletTree, FibreWaveletTreeStructure>::Type TWaveletTreeStructure;
	typedef typename Fibre<TWaveletTreeStructure, FibreTreeVertieces>::Type TWaveletTreeVertieces;
	typedef typename Value<TWaveletTreeVertieces>::Type TWaveletTreeVertex;
	typedef typename Value<TWaveletTree>::Type TChar;
	typedef typename Value<TWaveletTreeVertex, 2>::Type TPos;

String<TChar> text;
    generateText(text);
    resize(text, 1000);

    TWaveletTree waveletTree(text);

    clear(getFibre(waveletTree, FibreBitStrings()));

    fillWaveletTree_(waveletTree, text);

    for (unsigned i = 0; i < length(text); ++i)
    {
        SEQAN_ASSERT_EQ(getCharacter(waveletTree, i), text[i]);
    }
}

template <typename TWaveletTree>
void waveletTreeOpenSave(TWaveletTree & /*tag*/)
{
	typedef typename Fibre<TWaveletTree, FibreWaveletTreeStructure>::Type TWaveletTreeStructure;
	typedef typename Fibre<TWaveletTreeStructure, FibreTreeVertieces>::Type TWaveletTreeVertieces;
	typedef typename Value<TWaveletTreeVertieces>::Type TWaveletTreeVertex;
	typedef typename Value<TWaveletTree>::Type TChar;
	typedef typename Value<TWaveletTreeVertex, 2>::Type TPos;

    String<TChar> text;
    generateText(text);
    resize(text, 1000);

    TWaveletTree waveletTree(text);
    save(waveletTree, "test-waveletTree");

    TWaveletTree waveletTreeOpen;
    open(waveletTreeOpen, "test-waveletTree");
    SEQAN_ASSERT(waveletTree == waveletTreeOpen);
}

SEQAN_DEFINE_TEST(test_wavelet_tree_clear)
{
    using namespace seqan;

    WaveletTree<String<Dna>, void> dnaTag;
    WaveletTree<String<Dna5>, void> dna5Tag;
    WaveletTree<String<AminoAcid>, void> asTag;
    WaveletTree<String<char>, void> charTag;
    WaveletTree<String<unsigned char>, void> uCharTag;
    waveletTreeClear(dnaTag);
    waveletTreeClear(dna5Tag);
    waveletTreeClear(asTag);
    waveletTreeClear(charTag);
    waveletTreeClear(uCharTag);
}


SEQAN_DEFINE_TEST(test_wavelet_tree_constructor)
{
    using namespace seqan;

    WaveletTree<String<Dna>, void> dnaTag;
    WaveletTree<String<Dna5>, void> dna5Tag;
    WaveletTree<String<AminoAcid>, void> asTag;
    WaveletTree<String<char>, void> charTag;
    WaveletTree<String<unsigned char>, void> uCharTag;
    waveletTreeConstructor(dnaTag);
    waveletTreeConstructor(dna5Tag);
    waveletTreeConstructor(asTag);
    waveletTreeConstructor(charTag);
    waveletTreeConstructor(uCharTag);
}

SEQAN_DEFINE_TEST(test_wavelet_tree_dollar_position)
{
    using namespace seqan;

    WaveletTree<String<Dna>, FmiDollarSubstituted<SingleDollar<void> > > dnaTag;
    WaveletTree<String<Dna5>, FmiDollarSubstituted<SingleDollar<void> > > dna5Tag;
    WaveletTree<String<AminoAcid>, FmiDollarSubstituted<SingleDollar<void> > > asTag;
    WaveletTree<String<char>, FmiDollarSubstituted<SingleDollar<void> > > charTag;
    WaveletTree<String<unsigned char>, FmiDollarSubstituted<SingleDollar<void> > > uCharTag;
    waveletTreeDollarPosition(dnaTag);
    waveletTreeDollarPosition(dna5Tag);
    waveletTreeDollarPosition(asTag);
    waveletTreeDollarPosition(charTag);
    waveletTreeDollarPosition(uCharTag);
}

SEQAN_DEFINE_TEST(test_wavelet_tree_dollar_substitute)
{
    using namespace seqan;

    WaveletTree<String<Dna>, FmiDollarSubstituted<SingleDollar<void> > > dnaTag;
    WaveletTree<String<Dna5>, FmiDollarSubstituted<SingleDollar<void> > > dna5Tag;
    WaveletTree<String<AminoAcid>, FmiDollarSubstituted<SingleDollar<void> > > asTag;
    WaveletTree<String<char>, FmiDollarSubstituted<SingleDollar<void> > > charTag;
    WaveletTree<String<unsigned char>, FmiDollarSubstituted<SingleDollar<void> > > uCharTag;
    waveletTreeDollarSubstitute(dnaTag);
    waveletTreeDollarSubstitute(dna5Tag);
    waveletTreeDollarSubstitute(asTag);
    waveletTreeDollarSubstitute(charTag);
    waveletTreeDollarSubstitute(uCharTag);
}


SEQAN_DEFINE_TEST(test_wavelet_tree_empty)
{
    using namespace seqan;

    WaveletTree<String<Dna>, void> dnaTag;
    WaveletTree<String<Dna5>, void> dna5Tag;
    WaveletTree<String<AminoAcid>, void> asTag;
    WaveletTree<String<char>, void> charTag;
    WaveletTree<String<unsigned char>, void> uCharTag;
    waveletTreeEmpty(dnaTag);
    waveletTreeEmpty(dna5Tag);
    waveletTreeEmpty(asTag);
    waveletTreeEmpty(charTag);
    waveletTreeEmpty(uCharTag);
}




SEQAN_DEFINE_TEST(test_wavelet_tree_get_character)
{
    using namespace seqan;
    {   
        WaveletTree<String<Dna>, void> dnaTag;
        WaveletTree<String<Dna5>, void> dna5Tag;
        WaveletTree<String<AminoAcid>, void> asTag;
        WaveletTree<String<signed char>, void> charTag;
        WaveletTree<String<unsigned char>, void> uCharTag;

        waveletTreeGetCharacter(dnaTag);
        waveletTreeGetCharacter(dna5Tag);
        waveletTreeGetCharacter(asTag);
        waveletTreeGetCharacter(charTag);
        waveletTreeGetCharacter(uCharTag);
    }
    {
        WaveletTree<String<Dna>, FmiDollarSubstituted<SingleDollar<void> > > dnaTag;
        WaveletTree<String<Dna5>, FmiDollarSubstituted<SingleDollar<void> > > dna5Tag;
        WaveletTree<String<AminoAcid>, FmiDollarSubstituted<SingleDollar<void> > > asTag;
        WaveletTree<String<signed char>, FmiDollarSubstituted<SingleDollar<void> > > charTag;
        WaveletTree<String<unsigned char>, FmiDollarSubstituted<SingleDollar<void> > > uCharTag;

        waveletTreeGetCharacter(dnaTag);
        waveletTreeGetCharacter(dna5Tag);
        waveletTreeGetCharacter(asTag);
        waveletTreeGetCharacter(charTag);
        waveletTreeGetCharacter(uCharTag);
    }
}

SEQAN_DEFINE_TEST(test_wavelet_tree_get_fibre)
{
    using namespace seqan;

    WaveletTree<String<Dna>, void> dnaTag;
    WaveletTree<String<Dna5>, void> dna5Tag;
    WaveletTree<String<AminoAcid>, void> asTag;
    WaveletTree<String<char>, void> charTag;
    WaveletTree<String<unsigned char>, void> uCharTag;
    waveletTreeGetFibre(dnaTag);
    waveletTreeGetFibre(dna5Tag);
    waveletTreeGetFibre(asTag);
    waveletTreeGetFibre(charTag);
    waveletTreeGetFibre(uCharTag);
}

SEQAN_DEFINE_TEST(test_wavelet_tree_get_occ)
{
    using namespace seqan;

    {
        WaveletTree<String<Dna>, void> dnaTag;
        WaveletTree<String<Dna5>, void> dna5Tag;
        WaveletTree<String<AminoAcid>, void> asTag;
        WaveletTree<String<signed char>, void> charTag;
        WaveletTree<String<unsigned char>, void> uCharTag;

//         waveletTreeGetOcc(dnaTag);
//         waveletTreeGetOcc(dna5Tag);
//         waveletTreeGetOcc(asTag);
        waveletTreeGetOcc(charTag);
//         waveletTreeGetOcc(uCharTag);
    }
    {
        WaveletTree<String<Dna>, FmiDollarSubstituted<SingleDollar<void> > > dnaTag;
        WaveletTree<String<Dna5>, FmiDollarSubstituted<SingleDollar<void> > > dna5Tag;
        WaveletTree<String<AminoAcid>, FmiDollarSubstituted<SingleDollar<void> > > asTag;
        WaveletTree<String<signed char>, FmiDollarSubstituted<SingleDollar<void> > > charTag;
        WaveletTree<String<unsigned char>, FmiDollarSubstituted<SingleDollar<void> > > uCharTag;

//         waveletTreeGetOcc(dnaTag);
//         waveletTreeGetOcc(dna5Tag);
//         waveletTreeGetOcc(asTag);
        waveletTreeGetOcc(charTag);
//         waveletTreeGetOcc(uCharTag);
    }
}

SEQAN_DEFINE_TEST(test_wavelet_tree_num_vertieces)
{
    using namespace seqan;

    WaveletTree<String<Dna>, void> dnaTag;
    WaveletTree<String<Dna5>, void> dna5Tag;
    WaveletTree<String<AminoAcid>, void> asTag;
    WaveletTree<String<signed char>, void> charTag;
    WaveletTree<String<unsigned char>, void> uCharTag;

    waveletTreeNumVertieces(dnaTag);
    waveletTreeNumVertieces(dna5Tag);
    waveletTreeNumVertieces(asTag);
    waveletTreeNumVertieces(charTag);
    waveletTreeNumVertieces(uCharTag);
}

SEQAN_DEFINE_TEST(test_wavelet_tree_fill_wavelet_tree_)
{
    using namespace seqan;

    WaveletTree<String<Dna>, void> dnaTag;
    WaveletTree<String<Dna5>, void> dna5Tag;
    WaveletTree<String<AminoAcid>, void> asTag;
    WaveletTree<String<signed char>, void> charTag;
    WaveletTree<String<unsigned char>, void> uCharTag;

    waveletTreeFillWaveletTree_(dnaTag);
    waveletTreeFillWaveletTree_(dna5Tag);
    waveletTreeFillWaveletTree_(asTag);
    waveletTreeFillWaveletTree_(charTag);
    waveletTreeFillWaveletTree_(uCharTag);
}

SEQAN_DEFINE_TEST(test_wavelet_tree_open_save)
{
    using namespace seqan;

    WaveletTree<String<Dna>, void> dnaTag;
    WaveletTree<String<Dna5>, void> dna5Tag;
    WaveletTree<String<AminoAcid>, void> asTag;
    WaveletTree<String<signed char>, void> charTag;
    WaveletTree<String<unsigned char>, void> uCharTag;

    waveletTreeOpenSave(dnaTag);
    waveletTreeOpenSave(dna5Tag);
    waveletTreeOpenSave(asTag);
    waveletTreeOpenSave(charTag);
    waveletTreeOpenSave(uCharTag);
}

#endif  // TESTS_WAVELT_TREE_STRUCTURE_BETA_H_
