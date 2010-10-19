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

// Disabling warning 4267 assigning variables with different size on 32 and 64 bit.
#pragma warning( disable : 4267 )
// Disabling warning 4244, loss of data when values with different domain sizes.
#pragma warning( disable : 4244 )

#define finline __forceinline

typedef unsigned __int64 __uint64;

// Workaround for missing round() from C99 in Visual Studio.
template <typename T>
inline T round(T const & x)
{
	return floor(x + 0.5);
}

// Rename some underscore-functions in Windows.
#ifndef snprintf
#define snprintf _snprintf
#endif  // #ifndef snprintf

//define SEQAN_SWITCH_USE_FORWARDS to use generated forwards 
//#define SEQAN_SWITCH_USE_FORWARDS

// Define warning disabling macros as empty.
#ifndef SEQAN_PUSH_WARNING_DISABLE
// C4675: Disable warning "'function' : resolved overload was found by
// argument-dependent lookup".  Visual Studio warns because Koenig
// lookup was introduced in later version and behaviour has changed at
// some point.
//
// C4503: Disable warning for identifer name truncation.
//
// C4267: Disabling warning 4267 assigning variables with different
// size on 32 and 64 bit.  Need to re-enable this later.
//
// C4244: Disabling warning 4244, loss of data when values with
// different domain sizes.
//
// C4996: Do not warn against deprecated functions.
// TODO(holtgrew): Disable here, not on command line!
#define SEQAN_PUSH_WARNING_DISABLE \
    _Pragma("warning ( push )") \
    _Pragma("warning ( disable : 4675 )") \
    _Pragma("warning ( disable : 4503 )") \
    _Pragma("warning ( disable : 4267 )") \
    _Pragma("warning ( disable : 4244 )") \
    _Pragma("warning ( disable : 4996 )")
#endif  // #ifndef SEQAN_PUSH_WARNING_DISABLE

#ifndef SEQAN_POP_WARNING_DISABLE
#define SEQAN_POP_WARNING_DISABLE \
    _Pragma("warning ( pop )")
#endif  // #ifndef SEQAN_POP_WARNING_DISABLE
