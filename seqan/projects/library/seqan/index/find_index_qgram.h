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

#ifndef SEQAN_HEADER_INDEX_QGRAM_FIND_H
#define SEQAN_HEADER_INDEX_QGRAM_FIND_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// QGram finders

	struct FinderQGramLookup_; //Finder that simply looks up the q-gram in the hash table

/**
.Tag.Index Find Algorithm
..tag.QGramFindLookup:q-gram search.
Finds q-grams in a @Spec.IndexQGram@ index using the hash table.
..include:seqan/index.h
*/

	typedef Tag<FinderQGramLookup_> const QGramFindLookup;

//____________________________________________________________________________


	template < typename TText, typename TShapeSpec, typename TSpec >
	struct DefaultFinder<Index<TText, IndexQGram<TShapeSpec, TSpec> > > {
        typedef QGramFindLookup Type;
    };


//////////////////////////////////////////////////////////////////////////////
// _findFirstIndex implementation

	template < typename TText, typename TSpec, typename TSpecFinder, typename TPattern >
	inline void _findFirstIndex(
		Finder< Index<TText, TSpec>, TSpecFinder > &finder,
		TPattern const &pattern,
		QGramFindLookup const)
	{
		typedef Index<TText, TSpec>									TIndex;
		typedef typename Fibre<TIndex, QGram_SA>::Type				TSA;
		typedef typename Fibre<TIndex, QGramShape>::Type			TShape;
		typedef typename Fibre<TIndex, QGramDir>::Type				TDir;
		typedef typename Iterator<TSA const, Standard>::Type		TSAIterator;
		typedef typename Iterator<TPattern const, Standard>::Type	TPatternIterator;

		TIndex &index = haystack(finder);
		indexRequire(index, QGramSADir());

		TSAIterator saIt = begin(indexSA(index), Standard());
		TPatternIterator pIt = begin(pattern, Standard());
		TDir const &dir = indexDir(index);
		TShape &shape = indexShape(index);

		finder.range.i1 = saIt + dir[getBucket(index.bucketMap, hash(shape, pIt, length(pattern)))];
		finder.range.i2 = saIt + dir[getBucket(index.bucketMap, hashUpper(shape, pIt, length(pattern)))];
	}


//////////////////////////////////////////////////////////////////////////////

} //namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_
