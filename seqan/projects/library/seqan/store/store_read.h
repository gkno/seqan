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

#ifndef SEQAN_HEADER_STORE_READ_H
#define SEQAN_HEADER_STORE_READ_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Read Store
//////////////////////////////////////////////////////////////////////////////

/**
.Class.ReadStoreElement
..summary:Represents a single read (without sequence).
..cat:Fragment Store
..signature:ReadStoreElement<>
..signature:ReadStoreElement<TSpec>
..param.TSpec:The specialization type.
...default:$void$
..remarks:Value type of the @Memvar.FragmentStore#readStore@ string.

.Memfunc.ReadStoreElement#ReadStoreElement
..summary:Constructor
..signature:ReadStoreElement<>()
..signature:ReadStoreElement<TSpec> ()
..remarks:Sets $matePairId$ to $INVALID_ID$.
..class:Class.ReadStoreElement
.Memvar.ReadStoreElement#matePairId
..summary:Refers to a mate-pair in the @Memvar.FragmentStore#matePairStore@ or is $INVALID_ID$ if the read is not paired.
..type:Metafunction.Id
..class:Class.ReadStoreElement
.Memvar.ReadStoreElement#INVALID_ID
..summary:Constant to represent an invalid id.
..type:Metafunction.Id
..class:Class.ReadStoreElement
*/

template <typename TSpec = void>
struct ReadStoreElement
{
	typedef typename Id<ReadStoreElement>::Type TId;
	
	static const TId INVALID_ID = SupremumValue<typename Id<ReadStoreElement<TSpec> >::Type>::VALUE;

	TId matePairId;				// refers to the mate-pair, INVALID_ID if not part of a mate-pair

	ReadStoreElement() : matePairId(INVALID_ID) {}
};

//////////////////////////////////////////////////////////////////////////////


}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
