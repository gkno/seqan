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

#ifndef SEQAN_HEADER_STORE_ALIGN_H
#define SEQAN_HEADER_STORE_ALIGN_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Aligned Read Store
//////////////////////////////////////////////////////////////////////////////

template <typename TPos, typename TGapAnchor, typename TSpec = void>
struct AlignedReadStoreElement
{
	typedef typename Id<AlignedReadStoreElement>::Type TId;
	
	TId					readId;
	TId					contigId;
	TId					pairMatchId;	// unique id. for multiple mate-pair matches
	TPos				beginPos;		// begin position of the gapped sequence in gapped contig sequence
	TPos				endPos;			// end position of ..., for reverse aligned reads holds end < begin
	String<TGapAnchor>	gaps;
};

//////////////////////////////////////////////////////////////////////////////


}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
