 /*==========================================================================
                SeqAn - The Library for Sequence Analysis
                          http://www.seqan.de 
 ============================================================================
  Copyright (C) 2007

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  Lesser General Public License for more details.

 ============================================================================
  $Id$
 ==========================================================================*/

#ifndef SEQAN_HEADER_SHAPE_THRESHOLD_H
#define SEQAN_HEADER_SHAPE_THRESHOLD_H

namespace SEQAN_NAMESPACE_MAIN
{

struct ThreshQGramLemma_;
struct ThreshExact_;
struct ThreshHeuristic_;

typedef Tag<ThreshQGramLemma_> const	ThreshQGramLemma;
typedef Tag<ThreshExact_> const			ThreshExact;
typedef Tag<ThreshHeuristic_> const		ThreshHeuristic;

template <typename TShape, typename TPatternLength, typename TErrors, typename TDistance>
inline int qgramThreshold(TShape const & shape, TPatternLength patternLength, TErrors errors, TDistance const, ThreshQGramLemma const)
{
	return patternLength - length(shape) + 1 - errors * weight(shape);
}

template <typename TShape, typename TPatternSize, typename TErrors, typename TDistance>
int qgramThreshold(TShape const & shape, TPatternSize patternLength, TErrors errors, TDistance const, ThreshHeuristic const)
{
	String<unsigned char> coverage;
	String<bool> preserved;
	String<unsigned> ones;
	CharString bitmap;

	// initialize coverage map and bitmap of preserved q-grams
	fill(preserved, patternLength - length(shape) + 1, true);
	fill(coverage, patternLength, 0);

	shapeToString(bitmap, shape);
	for (unsigned i = 0; i < length(bitmap); ++i)
		if (bitmap[i] == '1')
		{
			appendValue(ones, i);
			for (unsigned j = 0; j < length(preserved); ++j)
				++coverage[i + j];
		}

	// greedily destroy a maximum number of q-grams
	for (; errors > 0; --errors)
	{
		// find position that destroys a maximum number of q-grams
		unsigned maxCoverage = 0;
		unsigned maxCoveragePos = 0;
		for (unsigned i = 0; i < length(coverage); ++i)
			if (maxCoverage < coverage[i])
			{
				maxCoverage = coverage[i];
				maxCoveragePos = i;
			}

		// destroy q-grams
		for (unsigned k = 0; k < length(ones); ++k)
			if (ones[k] <= maxCoveragePos)
			{
				unsigned startPos = maxCoveragePos - ones[k];
				if (startPos < length(preserved) && preserved[startPos])
				{
					preserved[startPos] = false;
					for (unsigned l = 0; l < length(ones); ++l)
						--coverage[startPos + l];
				}
			}
	}

	unsigned thresh = 0;
	for (unsigned i = 0; i < length(preserved); ++i)
		if (preserved[i])
			++thresh;

	return thresh;
}	

}	// namespace seqan

#endif
