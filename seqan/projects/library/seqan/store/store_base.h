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

#ifndef SEQAN_HEADER_STORE_BASE_H
#define SEQAN_HEADER_STORE_BASE_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Base structs
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

// We store gap anchors only for the first text character behind a gap or a clipped sequence character

template <typename TPos>
struct GapAnchor {
	TPos	seqPos;			// sequence character position in the ungapped sequence
	TPos	gapPos;			// sequence character position in the gapped sequence

	GapAnchor() : seqPos(0), gapPos(0) {}
	GapAnchor(TPos sP, TPos gP) : seqPos(sP), gapPos(gP) {}
};

//////////////////////////////////////////////////////////////////////////////


}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
