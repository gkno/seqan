// ==========================================================================
//                               score_blocks.h
//                           breakpoint_calculator
// ==========================================================================
// Copyright (C) 2012 by Birte Kehr
//
// This program is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option)
// any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
// more details.
//
// You should have received a copy of the GNU General Public License along
// with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// ==========================================================================
// Author: Birte Kehr <birte.kehr@fu-berlin.de>
// ==========================================================================

#ifndef SANDBOX_BKEHR_APPS_BREAKPOINT_CALCULATOR_SCORE_BLOCKS_H_
#define SANDBOX_BKEHR_APPS_BREAKPOINT_CALCULATOR_SCORE_BLOCKS_H_

#include <seqan/align.h>

using namespace seqan;

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// Computes the pairwise alignment score of row1 and row2 using specified scoring scheme.
// Igonores columns that contain gaps in both rows.
template<typename TSequence, typename TScoreValue, typename TScoreSpec>
TScoreValue
scoreRows(Gaps<TSequence, ArrayGaps> & row1,
		  Gaps<TSequence, ArrayGaps> & row2,
		  Score<TScoreValue, TScoreSpec> & sc,
		  bool out)
{
	typedef typename Position<TSequence>::Type TPosition;

	TScoreValue pairwiseScore = 0;

	TPosition pos1 = 0;
	TPosition pos2 = 0;
	TPosition len1 = length(row1) + beginPosition(row1);
	TPosition len2 = length(row2) + beginPosition(row2);

	char gap = 0; // 0 - no gap, 1 - gap row1, 2 - gap row2

	while (pos1 < len1 || pos2 < len2)
	{
		if (isGap(row1, pos1))
		{
			if (!isGap(row2, pos2))
			{
				// score gap in row1
				pairwiseScore += (gap == 1) ? scoreGapOpen(sc) : scoreGapExtend(sc);
				gap = 1;
			}
		}
		else if (isGap(row2, pos2))
		{
			// score gap in row2
			pairwiseScore += (gap == 2) ? scoreGapOpen(sc) : scoreGapExtend(sc);
			gap = 2;
		}
		else
		{
			// score match/mismatch
			pairwiseScore += score(sc, toSourcePosition(row1, pos1), toSourcePosition(row2, pos2), source(row1), source(row2));
		}

		++pos1;
		++pos2;
	}
	if (out){
		Align<TSequence, ArrayGaps> align;
		appendValue(rows(align), row1);
		appendValue(rows(align), row2);
		std::cout << align;
		std::cout << "Pairwise Score: " << pairwiseScore << std::endl;
	}
	return pairwiseScore;
}


// Computes the sum of all sum-of-pairs scores for the alignments in aligns using specified scoring scheme.
template<typename TSequence, typename TScoreValue, typename TScoreSpec>
int
sopScore(String<Align<TSequence, ArrayGaps> > & aligns, Score<TScoreValue, TScoreSpec> & score)
{
	typedef Align<TSequence, ArrayGaps> TAlign;
	typedef typename Iterator<typename Rows<TAlign>::Type >::Type TRowsIterator;

	TScoreValue totalScore = 0;
	
	// sum up score for all alignment blocks in aligns
	for (typename Iterator<String<TAlign> >::Type it = begin(aligns); it < end(aligns); ++it)
	{
		TScoreValue sopScore = 0;

		// sum up score for all pairs of alignment rows
		TRowsIterator endIt = end(rows(*it));
		for (TRowsIterator rowIt1 = begin(rows(*it)); rowIt1 < endIt - 1; ++rowIt1)
			for (TRowsIterator rowIt2 = rowIt1 + 1; rowIt2 < endIt; ++rowIt2)
				sopScore += scoreRows(*rowIt1, *rowIt2, score, false);

		totalScore += sopScore;
		//std::cout << *it << "Score: " << sopScore << std::endl;
	}

	std::cout << "Total sum-of-pairs score: " << totalScore << std::endl;

	return 0;
}

template<typename TSize, typename TSequence>
bool
countNucleotidesAndGaps(String<TSize> & numNucleotide,
						TSize & numGapOpen,
						TSize & numGapExtend,
						Gaps<TSequence, ArrayGaps> & row)
{
	typedef typename Position<TSequence>::Type TPosition;
	typedef typename Value<TSequence>::Type TAlph;

	TPosition pos = 0;
	TPosition len = length(row) + beginPosition(row);

	bool returnVal = false;

	bool gap = false;

	while (pos < len)
	{
		if (pos == 2981 && isGap(row, pos)) returnVal = true;
		if (isGap(row, pos))
		{
			if (gap) ++numGapExtend;
			else ++numGapOpen;
			gap = true;
		}
		else
		{
			TAlph nucl = source(row)[toSourcePosition(row, pos)];
			++numNucleotide[ordValue(nucl)];
			gap = false;
		}
		++pos;
	}
	return returnVal;
	// TODO: end gaps
}

template<typename TSize, typename TSequence>
void
countSubstitutions(String<String<TSize> > & numSubst,
				   Gaps<TSequence, ArrayGaps> & row1,
				   Gaps<TSequence, ArrayGaps> & row2)
{
	typedef typename Position<TSequence>::Type TPosition;
	typedef typename Value<TSequence>::Type TAlph;

	TPosition pos1 = 0;
	TPosition pos2 = 0;
	TPosition len1 = length(row1) + beginPosition(row1);
	TPosition len2 = length(row2) + beginPosition(row2);

	while (pos1 < len1 || pos2 < len2)
	{
		if (!isGap(row1, pos1) && !isGap(row2, pos2))
		{
				TAlph nucl1 = source(row1)[toSourcePosition(row1, pos1)];
				TAlph nucl2 = source(row2)[toSourcePosition(row2, pos2)];
				++numSubst[ordValue(nucl1)][ordValue(nucl2)];
		}
		++pos1;
		++pos2;
	}
}

template<typename TValue, typename TAlph, typename TSize>
void
computeScoreMatrix(Score<TValue, ScoreMatrix<TAlph> > & score,
				   String<String<TSize> > & numSubst)
{
	TSize alphSize = ValueSize<TAlph>::VALUE;

	// count total number of substitutions (totalSubst) and count nucleotide occurences (numNucl)
	double totalSubst = 0;
	String<TSize> numNucl;
	resize(numNucl, alphSize, 0);

	for (TSize i = 0; i < alphSize; ++i)
	{
		for (TSize j = 0; j < alphSize; ++j)
		{
			totalSubst += numSubst[i][j];
			numNucl[i] += numSubst[i][j] + numSubst[j][i];
		}
	}

	// compute log(q_ij/p_j) and set score
	// q_ij = numSubst / totalSubst
	// p_j = numNucl / totalNucl
	// totalNucl = 2 * totalSubst
	for (TSize i = 0; i < alphSize; ++i)
	{
		for (TSize j = 0; j < alphSize; ++j)
		{
			double probSubst = (numSubst[i][j] + numSubst[j][i]) / totalSubst;
			double probNucl = numNucl[j] / (2 * totalSubst);
			double sc = log(probSubst / probNucl);
			setScore(score, TAlph(i), TAlph(j), sc);
			//std::cout << sc << "  ";
		}
		//std::cout << std::endl;
	}
}

template<typename TValue, typename TAlphabet, typename TStream>
void
writeScoreMatrix(Score<TValue, ScoreMatrix<TAlphabet> > & matrix, TStream & file)
{
	file << "# gap open: " << scoreGapOpen(matrix) << "\n";
	file << "# gap extend: " << scoreGapExtend(matrix) << "\n";
	for (int i = 0; i < ValueSize<TAlphabet>::VALUE; ++i)
		file << "\t" << TAlphabet(i);
	file << "\n";

	for (int i = 0; i < ValueSize<TAlphabet>::VALUE; ++i)
	{
		file << TAlphabet(i);
		for (int j = 0; j < ValueSize<TAlphabet>::VALUE; ++j)
			file << "\t" << score(matrix, TAlphabet(i), TAlphabet(j));
		file << "\n";
	}
}

// Computes the sum of all sum-of-pairs scores for the alignments in aligns computing its own scoring matrix.
template<typename TSequence, typename TScoreValue, typename TAlphabet>
int
sopScore(String<Align<TSequence, ArrayGaps> > & aligns,
		 Score<TScoreValue, ScoreMatrix<TAlphabet> > & scoreMatrix,
		 CharString & matrixFile)
{
	typedef typename Size<TSequence>::Type TSize;

	typedef Align<TSequence, ArrayGaps> TAlign;
	typedef typename Iterator<typename Rows<TAlign>::Type >::Type TRowsIter;

	TSize alphSize = ValueSize<TAlphabet>::VALUE;

	// count number of substitions and gaps
	TSize numGapOpen = 0, numGapExtend = 0;
	String<TSize> numNucleotide;
	String<String<TSize> > numSubst;

	resize(numNucleotide, alphSize, 0);
	resize(numSubst, alphSize);

	typedef typename Iterator<String<String<TSize> > >::Type TSubstIter;
	for (TSubstIter it = begin(numSubst); it < end(numSubst); ++it)
		resize(*it, alphSize, 0);

	// sum up for all alignment blocks in aligns
	for (typename Iterator<String<TAlign> >::Type it = begin(aligns); it < end(aligns); ++it)
	{
		// sum up for all pairs of alignment rows
		TRowsIter endIt = end(rows(*it));
		for (TRowsIter rowIt1 = begin(rows(*it)); rowIt1 < endIt - 1; ++rowIt1)
		{
			countNucleotidesAndGaps(numNucleotide, numGapOpen, numGapExtend, *rowIt1);
			for (TRowsIter rowIt2 = rowIt1 + 1; rowIt2 < endIt; ++rowIt2)
				countSubstitutions(numSubst, *rowIt1, *rowIt2);
		}
	}

	// compute the scoring matrix
	computeScoreMatrix(scoreMatrix, numSubst);

	// compute the score from frequencies
	TScoreValue totalScore = numGapOpen * scoreGapOpen(scoreMatrix) + numGapExtend * scoreGapExtend(scoreMatrix);
	for (TSize i = 0; i < alphSize; ++i)
	{
		if (numSubst[i][i] > 0)
			totalScore += numSubst[i][i] * score(scoreMatrix, TAlphabet(i), TAlphabet(i));
		for (TSize j = i + 1; j < alphSize; ++j)
		{
			if (numSubst[i][j] > 0 || numSubst[j][i] > 0)
				totalScore += (numSubst[i][j] + numSubst[j][i]) * score(scoreMatrix, TAlphabet(i), TAlphabet(j));
		}
	}

	std::cout << "Total sum-of-pairs score: " << totalScore << std::endl;

	// Open matrix file and write matrix
	std::fstream stream(toCString(matrixFile), std::ios::out | std::ios::binary);
	if (!stream.good())
	{
		std::cerr << "WARNING: Could not open " << matrixFile << "!\n";
		return 1;
	}
	else writeScoreMatrix(scoreMatrix, stream);

	return 0;
}

#endif  // #ifndef SANDBOX_BKEHR_APPS_BREAKPOINT_CALCULATOR_SCORE_BLOCKS_H_
