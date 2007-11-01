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
  $Id: shape_gapped.h 996 2007-08-06 16:38:01Z weese@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

#ifndef SEQAN_HEADER_SHAPE_PREDEFINED_H
#define SEQAN_HEADER_SHAPE_PREDEFINED_H

//SEQAN_NO_GENERATED_FORWARDS: no forwards are generated for this file

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// some predefined gapped shapes

	// seed shape of PatternHunter 111010010100110111
	typedef FixedGappedShape< 
		HardwiredShape< 1, 1, 2, 3, 2, 3, 1, 2, 1, 1 > 
	> PatternHunterShape;


}	// namespace seqan

#endif
