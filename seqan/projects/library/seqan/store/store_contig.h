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

#ifndef SEQAN_HEADER_STORE_CONTIG_H
#define SEQAN_HEADER_STORE_CONTIG_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Contig Store
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template <typename _TContigSeq, typename _TGapAnchor, typename _TSpec = void>
struct ContigStoreElement
{
	typedef typename Id<ContigStoreElement>::Type	TId;
	
	typedef _TContigSeq			TContigSeq;
	typedef _TGapAnchor			TGapAnchor;
	typedef _TSpec				TSpec;
	typedef __int64				TPos;
	typedef String<TGapAnchor>	TGapAnchors;

	static const TId INVALID_ID;

	TContigSeq	seq;
	TGapAnchors	gaps;
	
// dynamic loading and disposing of contigs
	unsigned	usage;			// number of threads,... using this contig
	TId			fileId;
	TPos		fileBeginPos;
	TPos		fileEndPos;
};

template <typename _TSpec = void>
struct ContigFile
{
	typedef typename Id<ContigFile>::Type	TId;

	static const TId INVALID_ID;

	CharString		fileName;
	AutoSeqFormat	format;
	TId				firstContigId;	// first sequence of the file corresponds to this contigId
};

//////////////////////////////////////////////////////////////////////////////

template <typename _TContigSeq, typename _TGapAnchor, typename _TSpec>
const typename Id<ContigStoreElement<_TContigSeq, _TGapAnchor, _TSpec> >::Type
ContigStoreElement<_TContigSeq, _TGapAnchor, _TSpec>::INVALID_ID = SupremumValue<typename Id<ContigStoreElement<_TContigSeq, _TGapAnchor, _TSpec> >::Type>::VALUE;

template <typename _TSpec>
const typename Id<ContigFile<_TSpec> >::Type
ContigFile<_TSpec>::INVALID_ID = SupremumValue<typename Id<ContigFile<_TSpec> >::Type>::VALUE;

//////////////////////////////////////////////////////////////////////////////

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
