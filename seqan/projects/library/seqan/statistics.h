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
  $Id: sequence.h 2596 2008-08-25 16:17:13Z rausch@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

#ifndef SEQAN_HEADER_STATISTICS_H
#define SEQAN_HEADER_STATISTICS_H

//____________________________________________________________________________
// prerequisites

#include <cmath>
#include <seqan/align.h>
#include <seqan/index.h>
//____________________________________________________________________________

SEQAN_PUSH_WARNING_DISABLE  // Disable warnings from here.

#ifdef SEQAN_SWITCH_USE_FORWARDS
#include <seqan/statistics/statistics_generated_forwards.h>
#endif

#include <seqan/statistics/statistics_markov_model.h>
#include <seqan/statistics/statistics_base.h>

//____________________________________________________________________________

SEQAN_POP_WARNING_DISABLE  // Enable warnings again.

#endif //#ifndef SEQAN_HEADER_...
