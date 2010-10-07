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

/**
.Class.ContigStoreElement
..summary:Represents a single contig.
..cat:Fragment Store
..signature:ContigStoreElement<>
..signature:ContigStoreElement<TContigSeq[, TGapAnchor[, TSpec]]>
..param.TContigSeq:Type to store the contig sequence.
..param.TGapAnchor:Type of a contig gap anchor.
...type:Class.GapAnchor
..param.TSpec:The specialization type.
...default:$void$
..remarks:Value type of the @Memvar.FragmentStore#contigStore@ string.

.Typedef.ContigStoreElement#TContigSeq
..summary:Type of the $seq$ member.
..class:Class.ContigStoreElement
.Typedef.ContigStoreElement#TGapAnchors
..summary:Type of the $gaps$ member.
..class:Class.ContigStoreElement
.Typedef.ContigStoreElement#TPos
..summary:Type of the $fileBeginPos$ and $fileEndPos$ members.
..class:Class.ContigStoreElement
.Typedef.ContigStoreElement#TSpec
..summary:The specialization type.
..class:Class.ContigStoreElement


.Memfunc.ContigStoreElement#ContigStoreElement
..summary:Constructor
..signature:ContigStoreElement<> ()
..signature:ContigStoreElement<TContigSeq[, TGapAnchor[, TSpec]]> ()
..remarks:Sets $fileId$ to $INVALID_ID$ and $usage$, $fileBeginPos$ and $fileEndPos$ to $0$.
..class:Class.ContigStoreElement
.Memvar.ContigStoreElement#seq
..summary:Contig sequence.
..type:Typedef.ContigStoreElement#TContigSeq
..class:Class.ContigStoreElement
.Memvar.ContigStoreElement#gaps
..summary:String of contig gap anchors. Can be used to create a $Spec.AnchorGaps$ alignment row.
..type:Typedef.ContigStoreElement#TGapAnchors
..class:Class.ContigStoreElement
.Memvar.ContigStoreElement#usage
..summary:Counts the number of locks, see @Function.lockContigs@.
..class:Class.ContigStoreElement
.Memvar.ContigStoreElement#fileId
..summary:Refers to a file in the @Memvar.FragmentStore#contigFileStore@ or is $INVALID_ID$ if the contig has no file association.
..type:Metafunction.Id
..class:Class.ContigStoreElement
.Memvar.ContigStoreElement#fileBeginPos
..summary:Begin position of the contig sequence fragment in the file.
..type:Typedef.ContigStoreElement#TPos
..class:Class.ContigStoreElement
.Memvar.ContigStoreElement#fileEndPos
..summary:End position of the contig sequence fragment in the file.
..type:Typedef.ContigStoreElement#TPos
..class:Class.ContigStoreElement
.Memvar.ContigStoreElement#INVALID_ID
..summary:Constant to represent an invalid id.
..type:Metafunction.Id
..class:Class.ContigStoreElement
..include:seqan/store.h
*/

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

	ContigStoreElement() : usage(0), fileId(INVALID_ID), fileBeginPos(0), fileEndPos(0) {}
};

//////////////////////////////////////////////////////////////////////////////

template <typename _TContigSeq, typename _TGapAnchor, typename _TSpec> 
const typename Id<ContigStoreElement<_TContigSeq, _TGapAnchor, _TSpec> >::Type 
ContigStoreElement<_TContigSeq, _TGapAnchor, _TSpec>::INVALID_ID = SupremumValue<typename Id<ContigStoreElement<_TContigSeq, _TGapAnchor, _TSpec> >::Type>::VALUE; 

//////////////////////////////////////////////////////////////////////////////

/**
.Class.ContigFile
..summary:Represents a file containing contigs.
..cat:Fragment Store
..signature:ContigFile<>
..signature:ContigFile<TSpec>
..param.TSpec:The specialization type.
...default:$void$
..remarks:Value type of the @Memvar.FragmentStore#contigFileStore@ string.

.Memvar.ContigFile#fileName
..summary:Contig file name.
..type:Shortcut.CharString
..class:Class.ContigFile
.Memvar.ContigFile#format
..summary:Stores the contig file format, auto-detected in $Function.loadContigs$.
..type:Class.AutoSeqFormat
..class:Class.ContigFile
.Memvar.ContigFile#firstContigId
..summary:The $contigId$ of the first sequence in the file. Subsequent contig sequences have an increasing $contigId$.
..type:Metafunction.Id
..class:Class.ContigFile
..include:seqan/store.h
*/

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

template <typename _TSpec> 
const typename Id<ContigFile<_TSpec> >::Type 
ContigFile<_TSpec>::INVALID_ID = SupremumValue<typename Id<ContigFile<_TSpec> >::Type>::VALUE; 

//////////////////////////////////////////////////////////////////////////////

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
