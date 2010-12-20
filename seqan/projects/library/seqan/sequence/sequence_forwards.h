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

#ifndef SEQAN_HEADER_SEQUENCE_FORWARDS_H 
#define SEQAN_HEADER_SEQUENCE_FORWARDS_H 

//SEQAN_NO_GENERATED_FORWARDS: no forwards are generated for this file

//////////////////////////////////////////////////////////////////////////////

namespace SEQAN_NAMESPACE_MAIN 
{

// Workaround (copied from generated forwards) for VS 2003.
#if defined(_MSC_VER) && (_MSC_VER < 1400)
template <unsigned int SPACE > struct Block;       	// "projects\library\seqan/sequence\string_stack.h"(48)
template <typename THostspec > struct Packed;       	// "projects\library\seqan/sequence\string_packed.h"(33)
template <typename TValue, typename TSpec > class String;       	// "projects\library\seqan/sequence\string_base.h"(54)
template <typename TString, typename TSpec > class StringSet;       	// "projects\library\seqan/sequence\sequence_multiple.h"(98)

template <typename TValue, typename THostspec, typename TTag> inline typename Iterator<String<TValue, Packed<THostspec> >, Tag<TTag> const>::Type end(String<TValue, Packed<THostspec> > & me, Tag<TTag> const tag_);       	// "projects\library\seqan/sequence\string_packed.h"(470)
template <typename TValue, typename THostspec, typename TTag> inline typename Iterator<String<TValue, Packed<THostspec> > const, Tag<TTag> const>::Type end(String<TValue, Packed<THostspec> > const & me, Tag<TTag> const tag_);       	// "projects\library\seqan/sequence\string_packed.h"(478)
template <typename TValue, unsigned int SPACE, typename TSpec> inline typename Iterator<String<TValue, Block<SPACE> >, Tag<TSpec> const >::Type end(String<TValue, Block<SPACE> > &me, Tag<TSpec> const);       	// "projects\library\seqan/sequence\string_stack.h"(209)
template <typename TValue, unsigned int SPACE, typename TSpec> inline typename Iterator<String<TValue, Block<SPACE> > const, Tag<TSpec> const>::Type end(String<TValue, Block<SPACE> > const &me, Tag<TSpec> const);       	// "projects\library\seqan/sequence\string_stack.h"(217)
template <typename TString, typename TSpec, typename TTag> inline typename Iterator< StringSet< TString, TSpec >, Tag<TTag> const>::Type end(StringSet< TString, TSpec > & me, Tag<TTag> const tag);       	// "projects\library\seqan/sequence\sequence_multiple.h"(1398)
template <typename TString, typename TSpec, typename TTag> inline typename Iterator< StringSet< TString, TSpec > const, Tag<TTag> const>::Type end(StringSet< TString, TSpec > const & me, Tag<TTag> const tag);       	// "projects\library\seqan/sequence\sequence_multiple.h"(1405)
#endif  // defined(_MSC_VER) && (_MSC_VER < 1400)

//////////////////////////////////////////////////////////////////////////////

template <typename TFile, typename TData> 
void read(TFile & file, TData & data);       	// "projects/library/seqan/file/file_format_raw.h"(307)

//////////////////////////////////////////////////////////////////////////////

template <typename TFile, typename TData> 
void write(TFile & file, TData & data);       	// "projects/library/seqan/file/file_format_raw.h"(327)

template <typename TFile, typename TData> 
void write(TFile & file, TData const & data);   // "projects/library/seqan/file/file_format_raw.h"(335)


//////////////////////////////////////////////////////////////////////////////

} //namespace SEQAN_NAMESPACE_MAIN

//////////////////////////////////////////////////////////////////////////////

#endif

