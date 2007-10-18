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

	struct _Finder_QGramLookup; //Finder that simply looks up the q-gram in the hash table

/**
.Tag.QGram_FIND_Lookup:
..summary:Finding q-grams in index using a hash table.
..general:Class.Finder
..cat:Index
..signature:Finder<TIndex>
..signature:Finder<TIndex, QGram_FIND_Lookup>
..param.TIndex:The index type.
...type:Spec.Index_QGram
*/

	typedef Tag<_Finder_QGramLookup> const QGram_FIND_Lookup;

//____________________________________________________________________________


	template < typename TText, typename TSpec >
	struct DefaultFinder<Index<TText, Index_QGram<TSpec> > > {
        typedef QGram_FIND_Lookup Type;
    };


//////////////////////////////////////////////////////////////////////////////
// _findFirstIndex implementation

	template < typename TText, typename TSpec, typename TSpecFinder, typename TPattern >
	inline void _findFirstIndex(
		Finder< Index<TText, TSpec>, TSpecFinder > &finder,
		TPattern const &pattern,
		QGram_FIND_Lookup const)
	{
		typedef Index<TText, TSpec>									TIndex;
		typedef typename Fibre<TIndex, QGram_SA>::Type				TSA;
		typedef typename Fibre<TIndex, QGram_Shape>::Type			TShape;
		typedef typename Fibre<TIndex, QGram_Dir>::Type				TDir;
		typedef typename Iterator<TSA const, Standard>::Type		TSAIterator;
		typedef typename Iterator<TPattern const, Standard>::Type	TPatternIterator;

		TIndex &index = haystack(finder);
		indexRequire(index, QGram_SADir());

		TSAIterator saIt = begin(indexSA(index), Standard());
		TPatternIterator pIt = begin(pattern, Standard());
		TDir const &dir = indexDir(index);
		TShape &shape = indexShape(index);

		finder.range.i1 = saIt + dir[hash(shape, pIt, length(pattern))];
		finder.range.i2 = saIt + dir[hashUpper(shape, pIt, length(pattern))];
	}


//////////////////////////////////////////////////////////////////////////////

} //namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_
