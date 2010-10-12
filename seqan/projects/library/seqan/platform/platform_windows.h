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

//SEQAN_NO_GENERATED_FORWARDS: no forwards are generated for this file

#define PLATFORM_WINDOWS
#define PLATFORM_WINDOWS_VS

// Disable warning "'function' : resolved overload was found by
// argument-dependent lookup".  Visual Studio warns because Koenig
// lookup was introduced in later version and behaviour has changed at some
// point.
#pragma warning( disable : 4675 )
// Disable warning for identifer name truncation.
#pragma warning( disable : 4503 )

#define finline __forceinline

typedef unsigned __int64 __uint64;

// Workaround for missing round() from C99 in Visual Studio.
template <typename T>
inline T round(T const & x)
{
	return floor(x + 0.5);
}

//define SEQAN_SWITCH_USE_FORWARDS to use generated forwards 
//#define SEQAN_SWITCH_USE_FORWARDS
