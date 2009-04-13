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

#define PLATFORM_WINDOWS
#define PLATFORM_WINDOWS_MINGW

#include <inttypes.h>

#define finline __inline__

typedef unsigned __int64 __uint64;

//define SEQAN_SWITCH_USE_FORWARDS to use generated forwards 
#define SEQAN_SWITCH_USE_FORWARDS
