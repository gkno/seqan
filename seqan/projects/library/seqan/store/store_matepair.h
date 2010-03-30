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

#ifndef SEQAN_HEADER_STORE_MATEPAIR_H
#define SEQAN_HEADER_STORE_MATEPAIR_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Mate Store
//////////////////////////////////////////////////////////////////////////////

/**
.Class.MatePairStoreElement
..summary:Represents a mate-pair.
..cat:Fragment Store
..signature:MatePairStoreElement<>
..signature:MatePairStoreElement<TSpec>
..param.TSpec:The specialization type.
...default:$void$
..remarks:A mate-pair consists of two reads sequenced from opposite ends and strands of the same fragment.
The insert size of a mate-pair is the size of the fragment.
..remarks:Value type of the @Memvar.FragmentStore#matePairStore@ string.

.Memfunc.MatePairStoreElement#MatePairStoreElement
..summary:Constructor
..signature:MatePairStoreElement<> ()
..signature:MatePairStoreElement<TSpec> ()
..remarks:Sets $readId[0]$, $readId[1]$ and $libId$ to $INVALID_ID$.
..class:Class.MatePairStoreElement
.Memvar.MatePairStoreElement#readId[2]
..summary:Refers to two paired reads in the @Memvar.FragmentStore#readStore@ or contains $INVALID_ID$ values.
..type:Metafunction.Id
..class:Class.MatePairStoreElement
.Memvar.MatePairStoreElement#libId
..summary:Refers to a library in the @Memvar.FragmentStore#libraryStore@ or is $INVALID_ID$ if the mate-pair has no library.
..type:Metafunction.Id
..class:Class.MatePairStoreElement
.Memvar.MatePairStoreElement#INVALID_ID
..summary:Constant to represent an invalid id.
..type:Metafunction.Id
..class:Class.MatePairStoreElement
*/

template <typename TSpec = void>
struct MatePairStoreElement
{
	typedef typename Id<MatePairStoreElement>::Type TId;

	static const TId INVALID_ID;
	
	TId		readId[2];	// refers to the two reads of a mate-pair, INVALID_ID if this is a singleton fragment (e.g. in afg: reads refer to fragments (mate pairs) and these refer to libraries, singletons refer to an empty fragment)
	TId		libId;

	MatePairStoreElement() : libId(INVALID_ID) 
	{
		readId[0] = INVALID_ID;
		readId[1] = INVALID_ID;
	}
};

//////////////////////////////////////////////////////////////////////////////

template <typename TSpec>
const typename Id<MatePairStoreElement<TSpec> >::Type
MatePairStoreElement<TSpec>::INVALID_ID = SupremumValue<typename Id<MatePairStoreElement<TSpec> >::Type>::VALUE;

//////////////////////////////////////////////////////////////////////////////

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
