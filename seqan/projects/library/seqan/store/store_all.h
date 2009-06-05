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

#ifndef SEQAN_HEADER_STORE_ALL_H
#define SEQAN_HEADER_STORE_ALL_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Contig Store Configuration
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template <typename TSpec = void>
struct FragmentStoreConfig 
{
	typedef String<Dna5Q>	TReadSeq;
	typedef String<Dna5Q>	TContigSeq;
	
	typedef double			TMean;
	typedef double			TStd;
		
	typedef void			TReadStoreElementSpec;
	typedef void			TMatePairStoreElementSpec;
	typedef void			TLibraryStoreElementSpec;
	typedef void			TContigStoreElementSpec;
	typedef void			TAlignedReadStoreElementSpec;
};

//////////////////////////////////////////////////////////////////////////////
// Contig Store
//////////////////////////////////////////////////////////////////////////////

template <typename TSpec = void, typename TConfig = FragmentStoreConfig<TSpec> >
struct _FragmentStore
{
	typedef typename TConfig::TMean					TMean;
	typedef typename TConfig::TStd					TStd;
	
	typedef typename TConfig::TReadSeq				TReadSeq;
	typedef typename TConfig::TContigSeq			TContigSeq;
	typedef typename Position<TReadSeq>::Type		TReadPos;
	typedef typename Position<TContigSeq>::Type		TContigPos;
	
	typedef GapAnchor<TReadPos>						TReadGapAnchor;
	typedef GapAnchor<TContigPos>					TContigGapAnchor;

	typedef typename TConfig::TReadStoreElementSpec			TReadStoreElementSpec;
	typedef typename TConfig::TMatePairStoreElementSpec		TMatePairStoreElementSpec;
	typedef typename TConfig::TLibraryStoreElementSpec		TLibraryStoreElementSpec;
	typedef typename TConfig::TContigStoreElementSpec		TContigStoreElementSpec;
	typedef typename TConfig::TAlignedReadStoreElementSpec	TAlignedReadStoreElementSpec;
	
	typedef String< ReadStoreElement< TReadSeq, TReadPos, TReadStoreElementSpec > >							TReadStore;
	typedef String< MatePairStoreElement< TMatePairStoreElementSpec > >										TMatePairStore;
	typedef String< LibraryStoreElement< TMean, TStd, TLibraryStoreElementSpec > >							TLibraryStore;
	typedef String< ContigStoreElement< TContigSeq, TContigGapAnchor, TContigStoreElementSpec > >			TContigStore;
	typedef String< AlignedReadStoreElement< TContigPos, TReadGapAnchor, TAlignedReadStoreElementSpec > >	TAlignedReadStore;
	
	TReadStore			readStore;
	TMatePairStore		matePairStore;
	TLibraryStore		libraryStore;
	TContigStore		contigStore;
	TAlignedReadStore	alignedReadStore;
	
	StringSet<CharString>	readNameStore;
	StringSet<CharString>	matePairNameStore;
	StringSet<CharString>	libraryNameStore;
	StringSet<CharString>	contigNameStore;
};

//////////////////////////////////////////////////////////////////////////////


}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
