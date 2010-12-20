// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2010, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: David Weese <david.weese@fu-berlin.de>
// ==========================================================================

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
ContigStoreElement<_TContigSeq, _TGapAnchor, _TSpec>::INVALID_ID = MaxValue<typename Id<ContigStoreElement<_TContigSeq, _TGapAnchor, _TSpec> >::Type>::VALUE; 

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
ContigFile<_TSpec>::INVALID_ID = MaxValue<typename Id<ContigFile<_TSpec> >::Type>::VALUE; 

//////////////////////////////////////////////////////////////////////////////

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
