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

#ifndef SEQAN_HEADER_INDEX_MANUAL_FORWARDS_H 
#define SEQAN_HEADER_INDEX_MANUAL_FORWARDS_H 

//SEQAN_NO_GENERATED_FORWARDS: no forwards are generated for this file

//////////////////////////////////////////////////////////////////////////////
// CLASSES
//////////////////////////////////////////////////////////////////////////////

namespace SEQAN_NAMESPACE_MAIN {

	struct FibreText_;		// Original text. Can be a String or a StringSet
	struct FibreRawText_;	// Concatenation of the strings above
	struct FibreSA_;		// suffix array (of raw text with virtual $-delimiters) with Pair entries
	struct FibreRawSA_;	// suffix array with integer entries
	struct FibreSae_;		// suffix array reordered in a b-tree
	struct FibreLcp_;		// lcp table of raw text
	struct FibreLcpe_;		// lcp interval tree
	struct FibreChildtab_;	// childtab (Kurtz et al.) of raw text
	struct FibreBwt_;		// burrows wheeler table of raw text

	typedef Tag<FibreText_> const		FibreText;
	typedef Tag<FibreRawText_> const	FibreRawText;
	typedef Tag<FibreSA_> const		FibreSA;
	typedef Tag<FibreRawSA_> const		FibreRawSA;
	typedef Tag<FibreSae_> const		FibreSae;
	typedef Tag<FibreLcp_> const		FibreLcp;
	typedef Tag<FibreLcpe_> const		FibreLcpe;
	typedef Tag<FibreChildtab_> const	FibreChildtab;
	typedef Tag<FibreBwt_> const		FibreBwt;

}

#endif

