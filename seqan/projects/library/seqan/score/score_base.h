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

#ifndef SEQAN_HEADER_SCORE_BASE_H
#define SEQAN_HEADER_SCORE_BASE_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////

struct Simple;

//////////////////////////////////////////////////////////////////////////////


/**
.Class.Score:
..cat:Miscellaneous
..summary:A scoring scheme.
..signature:Score<TValue, TSpec>
..param.TValue:The value type.
...default:int
..param.TSpec:The specializing type.
...default:@Spec.Simple Score@
*/
template <typename TValue = int, typename TSpec = Simple>
class Score;

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSpec>
struct Value< Score<TValue, TSpec> >
{
	typedef TValue Type;
};

//////////////////////////////////////////////////////////////////////////////

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
