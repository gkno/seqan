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

#ifndef SEQAN_HEADER_STORE_LIBRARY_H
#define SEQAN_HEADER_STORE_LIBRARY_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Library Store
//////////////////////////////////////////////////////////////////////////////

/**
.Class.LibraryStoreElement
..summary:Represents a fragment library.
..cat:Fragment Store
..signature:LibraryStoreElement<>
..signature:LibraryStoreElement<TMean[, TStd[, TSpec]]>
..param.TMean:The type to represent the library size mean.
...default:$double$
..param.TStd:The type to represent the library size standard deviation.
...default:$double$
..param.TSpec:The specialization type.
...default:$void$
..remarks:A fragment library is a set of mate-pairs having a certain distribution of insert sizes.
..remarks:Value type of the @Memvar.FragmentStore#libraryStore@ string.

.Memfunc.LibraryStoreElement#LibraryStoreElement
..summary:Constructor
..signature:LibraryStoreElement<>()
..signature:LibraryStoreElement<TMean[, TStd[, TSpec]]> ()
..remarks:Sets $mean$ and $std$ to $0$.
..class:Class.LibraryStoreElement
.Memvar.LibraryStoreElement#mean
..summary:The library size mean.
..class:Class.LibraryStoreElement
.Memvar.LibraryStoreElement#std
..summary:The library size standard deviation.
..class:Class.LibraryStoreElement
*/

template <typename TMean = double, typename TStd = double, typename TSpec = void>
struct LibraryStoreElement
{
	TMean		mean;		// mean library size in bps
	TStd		std;	// library size variance

	LibraryStoreElement() : mean(0), std(0) {}
};

//////////////////////////////////////////////////////////////////////////////


}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
