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

template <typename TPos>
struct GapAnchor {
	TPos	seqPos;			// sequence character position in the ungapped sequence
	TPos	gapPos;			// sequence character position in the gapped sequence
};

// We store gap anchors
//  * of the first text character
//  * of the last text character before a gap
//  * of the last alignment character

//////////////////////////////////////////////////////////////////////////////


}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
