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

#ifndef SEQAN_HEADER_STORE_IO_UCSC_H
#define SEQAN_HEADER_STORE_IO_UCSC_H

namespace SEQAN_NAMESPACE_MAIN
{

/**
.Tag.File Format.tag.UCSC:
	UCSC Genome Browser annotation file (a.k.a. knownGene format).
..include:seqan/store.h
*/
struct TagUCSC_;
typedef Tag<TagUCSC_> const UCSC;

	
//////////////////////////////////////////////////////////////////////////////
// _parse_readUCSCIdentifier
    
    template<typename TFile, typename TString, typename TChar>
    inline void
    _parse_readUCSCIdentifier(TFile & file, TString & str, TChar& c)
    {
        if (c == ' ' || c == '\t' || c == '\n') return;
        append(str, c);
        while (!_streamEOF(file)) 
		{
            c = _streamGet(file);
            if (c == ' ' || c == '\t' || c == '\n') return;
            append(str, c);
        }
    }

	template<typename TFile, typename TChar>
	inline void 
	_parse_skipWhiteComma(TFile& file, TChar& c)
	{
		if (c != ',' && c != ' ') return;
		while (!_streamEOF(file)) {
			c = _streamGet(file);
			if (c != ',' && c != ' ') break;
		}
	}
	
//////////////////////////////////////////////////////////////////////////////
// Read UCSC
//////////////////////////////////////////////////////////////////////////////

template <typename TFragmentStore, typename TSpec = void>
struct _IOContextUCSC
{
	typedef typename TFragmentStore::TAnnotationStore   TAnnotationStore;
	typedef typename Value<TAnnotationStore>::Type      TAnnotation;
	typedef typename TAnnotation::TId                   TId;

	CharString		transName;
	CharString		contigName;
	__int64			cdsBegin;
	__int64			cdsEnd;
	String<__int64>	exonBegin;
	String<__int64>	exonEnd;
	CharString		proteinName;
	
	enum { KNOWN_GENE, KNOWN_ISOFORMS } format;
	TAnnotation annotation;
};

template <typename TFragmentStore, typename TSpec>
inline void clear(_IOContextUCSC<TFragmentStore, TSpec> &ctx)
{
	typedef typename TFragmentStore::TAnnotationStore   TAnnotationStore;
	typedef typename Value<TAnnotationStore>::Type      TAnnotation;

	clear(ctx.transName);
	clear(ctx.contigName);
	clear(ctx.exonBegin);
	clear(ctx.exonEnd);
	clear(ctx.proteinName);
	clear(ctx.annotation.values);
}

//////////////////////////////////////////////////////////////////////////////
// _readOneAnnotation
//
// reads in one annotation line from a GFF file

template <typename TFile, typename TChar, typename TFragmentStore, typename TSpec>
inline bool 
_readOneAnnotation (
	TFile & file,
	TChar & c,
	_IOContextUCSC<TFragmentStore, TSpec> & ctx)
{
	typedef typename TFragmentStore::TContigPos         TContigPos;	
	typedef typename TFragmentStore::TAnnotationStore   TAnnotationStore;
	typedef typename Value<TAnnotationStore>::Type      TAnnotation;
	typedef typename TAnnotation::TId                   TId;
		
	clear(ctx);

	// read fields of alignments line        
	_parse_skipWhitespace(file, c);
	
	// read column 1: transcript name
	// The letters until the first whitespace will be read.
	// Then, we skip until we hit the first tab character.
	_parse_readUCSCIdentifier(file, ctx.transName, c);
	if (!empty(ctx.transName) && ctx.transName[0] == '#')
	{
		_parse_skipLine(file, c);
		return false;
	}
	_parse_skipUntilChar(file, '\t', c);
	c = _streamGet(file);
	
	// read column 2: contig name
	_parse_readUCSCIdentifier(file, ctx.contigName, c);
	_parse_skipSpace(file, c);
	
	// read column 3: orientation
	if (c != '+' && c != '-')
	{
		ctx.format = ctx.KNOWN_ISOFORMS;
		assign(prefix(ctx.transName, 0), "GENE");
		_parse_skipLine(file, c);
		return true;
	}
	ctx.format = ctx.KNOWN_GENE;
	char orientation = c;
	c = _streamGet(file);
	_parse_skipWhitespace(file, c);

	// read column 4: transcript begin position
	if (_parse_isDigit(c))
		ctx.annotation.beginPos = _parse_readNumber(file, c);
	else
	{
		ctx.annotation.beginPos = TAnnotation::INVALID_POS;
		_parse_skipEntryUntilWhitespace(file, c);
	}
	_parse_skipWhitespace(file, c);

	// read column 5: transcript end position
	if (_parse_isDigit(c))
		ctx.annotation.endPos = _parse_readNumber(file, c);
	else 
	{
		ctx.annotation.endPos = TAnnotation::INVALID_POS;
		_parse_skipEntryUntilWhitespace(file, c);
	}
	_parse_skipWhitespace(file, c);	

	// read column 6: CDS begin position
	if (_parse_isDigit(c))
		ctx.cdsBegin = _parse_readNumber(file, c);
	else
	{
		ctx.cdsBegin = TAnnotation::INVALID_POS;
		_parse_skipEntryUntilWhitespace(file, c);
	}
	_parse_skipWhitespace(file, c);
	
	// read column 7: CDS end position
	if (_parse_isDigit(c))
		ctx.cdsEnd = _parse_readNumber(file, c);
	else 
	{
		ctx.cdsEnd = TAnnotation::INVALID_POS;
		_parse_skipEntryUntilWhitespace(file, c);
	}
	_parse_skipWhitespace(file, c);	
	
	// read column 8: exon count
	int exons = -1;
	if (_parse_isDigit(c))
		exons = _parse_readNumber(file, c);
	_parse_skipWhitespace(file, c);

	// read column 9: exon begin positions
	for (int i = 0; i < exons; ++i)
	{
		appendValue(ctx.exonBegin, _parse_readNumber(file, c), Generous());
		_parse_skipWhiteComma(file, c);
	}
	_parse_skipUntilChar(file, '\t', c);
	c = _streamGet(file);
	
	// read column 10: exon end positions
	for (int i = 0; i < exons; ++i)
	{
		appendValue(ctx.exonEnd, _parse_readNumber(file, c));
		_parse_skipWhiteComma(file, c);
	}
	_parse_skipUntilChar(file, '\t', c);
	c = _streamGet(file);
	
	// read column 10: protein name
	_parse_readUCSCIdentifier(file, ctx.proteinName, c);
	_parse_skipUntilChar(file, '\t', c);
	c = _streamGet(file);

	// skip column 11
	_parse_skipEntryUntilWhitespace(file, c);
	_parse_skipWhitespace(file, c);

	// adapt positions
	if (orientation == '-')
	{
		TContigPos tmp = ctx.annotation.beginPos;
		ctx.annotation.beginPos = ctx.annotation.endPos;
		ctx.annotation.endPos = tmp;
		tmp = ctx.cdsBegin;
		ctx.cdsBegin = ctx.cdsEnd;
		ctx.cdsEnd = tmp;
		for (int i = 0; i < exons; ++i)
		{
			tmp = ctx.exonBegin[i];
			ctx.exonBegin[i] = ctx.exonEnd[i];
			ctx.exonEnd[i] = tmp;
		}
	}

	return true;
}

template <typename TFragmentStore, typename TSpec>
inline void 
_storeOneAnnotation_KNOWN_GENE (
	TFragmentStore & fragStore,
	_IOContextUCSC<TFragmentStore, TSpec> & ctx)
{
	typedef typename TFragmentStore::TAnnotationStore   TAnnotationStore;
	typedef typename Value<TAnnotationStore>::Type      TAnnotation;
	typedef typename TAnnotation::TId                   TId;
	
	SEQAN_ASSERT_EQ(length(fragStore.annotationStore), length(fragStore.annotationNameStore));

	TId annoStoreLen = length(fragStore.annotationStore);
	TId transId = TAnnotation::INVALID_ID;
	TId cdsId = TAnnotation::INVALID_ID;

	// add transcript and CDS
	_storeAppendAnnotationName(fragStore, transId, ctx.transName);
	cdsId = length(fragStore.annotationNameStore);
	appendName(fragStore.annotationNameStore, ctx.proteinName, fragStore.annotationNameStoreCache);
	
	if (annoStoreLen <= transId)
		annoStoreLen = transId + 1;
	
	if (annoStoreLen <= cdsId)
		annoStoreLen = cdsId + 1;
	
	resize(fragStore.annotationStore, annoStoreLen + length(ctx.exonBegin), Generous());
	resize(fragStore.annotationNameStore, annoStoreLen + length(ctx.exonBegin), Generous());

	// add contig name
	_storeAppendContig(fragStore, ctx.annotation.contigId, ctx.contigName);	

	TAnnotation &transcript = fragStore.annotationStore[transId];
	TId geneId = transcript.parentId;
	if (geneId == TAnnotation::INVALID_ID) geneId = 0;
	transcript = ctx.annotation;
	transcript.parentId = geneId;
	transcript.typeId = TFragmentStore::ANNO_MRNA;
	
	TAnnotation &cds = fragStore.annotationStore[cdsId];
	cds = ctx.annotation;
	cds.parentId = transId;
	cds.typeId = TFragmentStore::ANNO_CDS;
	cds.beginPos = ctx.cdsBegin;
	cds.endPos = ctx.cdsEnd;
	_adjustParent(transcript, cds);
	
	// insert exons
	ctx.annotation.parentId = transId;
	ctx.annotation.typeId = TFragmentStore::ANNO_EXON;
	for (unsigned i = 0; i < length(ctx.exonBegin); ++i)
	{
		ctx.annotation.beginPos = ctx.exonBegin[i];
		ctx.annotation.endPos = ctx.exonEnd[i];
		fragStore.annotationStore[annoStoreLen + i] = ctx.annotation;
		_adjustParent(transcript, ctx.annotation);
	}
	if (geneId != 0)
		_adjustParent(fragStore.annotationStore[geneId], transcript);
}

template <typename TFragmentStore, typename TSpec>
inline void 
_storeOneAnnotation_KNOWN_ISOFORMS (
	TFragmentStore & fragStore,
	_IOContextUCSC<TFragmentStore, TSpec> & ctx)
{
	typedef typename TFragmentStore::TAnnotationStore   TAnnotationStore;
	typedef typename Value<TAnnotationStore>::Type      TAnnotation;
	typedef typename TAnnotation::TId                   TId;
	
	SEQAN_ASSERT_EQ(length(fragStore.annotationStore), length(fragStore.annotationNameStore));
	
	TId annoStoreLen = length(fragStore.annotationStore);
	TId geneId = TAnnotation::INVALID_ID;
	TId transId = TAnnotation::INVALID_ID;
	
	// add transcript and CDS
	_storeAppendAnnotationName(fragStore, geneId, ctx.transName);
	_storeAppendAnnotationName(fragStore, transId, ctx.contigName);
	
	if (annoStoreLen <= geneId)
		annoStoreLen = geneId + 1;
	
	if (annoStoreLen <= transId)
		annoStoreLen = transId + 1;
	
	resize(fragStore.annotationStore, annoStoreLen, Generous());
	resize(fragStore.annotationNameStore, annoStoreLen, Generous());
	
	TAnnotation &locus = fragStore.annotationStore[geneId];
	locus.parentId = 0;
	locus.typeId = TFragmentStore::ANNO_GENE;

	TAnnotation &transcript = fragStore.annotationStore[transId];
	transcript.parentId = geneId;
	transcript.typeId = TFragmentStore::ANNO_MRNA;
	
	_adjustParent(locus, transcript);
}

template <typename TFragmentStore, typename TSpec>
inline void 
_storeOneAnnotation (
	TFragmentStore & fragStore,
	_IOContextUCSC<TFragmentStore, TSpec> & ctx)
{
	if (ctx.format == ctx.KNOWN_GENE)
		_storeOneAnnotation_KNOWN_GENE(fragStore, ctx);
	else
		_storeOneAnnotation_KNOWN_ISOFORMS(fragStore, ctx);
}

template<typename TFile, typename TSpec, typename TConfig>
inline void 
read (
	TFile & file,
	FragmentStore<TSpec, TConfig> & fragStore,
	UCSC)
{
	typedef FragmentStore<TSpec, TConfig> TFragmentStore;
	
	if (_streamEOF(file)) return;

	// get first character from the stream
	char c = _streamGet(file);
	_IOContextUCSC<TFragmentStore> ctx;
	
	refresh(fragStore.contigNameStoreCache);
	refresh(fragStore.annotationNameStoreCache);
	refresh(fragStore.annotationTypeStoreCache);
	
	while (!_streamEOF(file))
	{
		if (_readOneAnnotation(file, c, ctx))
			_storeOneAnnotation(fragStore, ctx);
	}
	_storeCreateAnnoBackLinks(fragStore.annotationStore);
	_storeRemoveTempAnnoNames(fragStore);
}

//////////////////////////////////////////////////////////////////////////////
// Write UCSC
//////////////////////////////////////////////////////////////////////////////

template <typename TFragmentStore, typename TSpec, typename TAnnotation, typename TId>
inline bool 
_retrieveOneAnnotation (
	TFragmentStore & fragStore,
	_IOContextUCSC<TFragmentStore, TSpec> & ctx,
	TAnnotation &annotation,
	TId id)
{	
	if (annotation.typeId == TFragmentStore::ANNO_GENE)
	{
		if (annotation.lastChildId == TAnnotation::INVALID_ID) return false;
		ctx.format = ctx.KNOWN_ISOFORMS;
		ctx.transName = getAnnoUniqueName(fragStore, id);
		ctx.contigName = getAnnoUniqueName(fragStore, annotation.lastChildId);
		return true;
	}
	
	if (annotation.typeId == TFragmentStore::ANNO_MRNA)
	{
		ctx.format = ctx.KNOWN_GENE;
		ctx.transName = getAnnoUniqueName(fragStore, id);
		if (annotation.contigId < length(fragStore.contigNameStore))
			ctx.contigName = fragStore.contigNameStore[annotation.contigId];
		else
			clear(ctx.contigName);
		
		ctx.annotation = annotation;
		clear(ctx.proteinName);
		clear(ctx.exonBegin);
		clear(ctx.exonEnd);
		
		TId lastChildId = annotation.lastChildId;
		TId i = lastChildId;
		do {
			i = fragStore.annotationStore[i].nextSiblingId;
			TAnnotation &anno = fragStore.annotationStore[i];
			if (anno.typeId == TFragmentStore::ANNO_CDS)
			{
				if (i < length(fragStore.annotationNameStore))
					ctx.proteinName = fragStore.annotationNameStore[i];
				ctx.cdsBegin = anno.beginPos;
				ctx.cdsEnd = anno.endPos;
			}
			if (anno.typeId == TFragmentStore::ANNO_EXON)
			{
				appendValue(ctx.exonBegin, anno.beginPos, Generous());
				appendValue(ctx.exonEnd, anno.endPos, Generous());
			}
		} while (i != lastChildId);
		return true;
	}
	return false;
}

template<typename TTargetStream, typename TFragmentStore, typename TSpec>
inline void 
_writeOneAnnotation (
	TTargetStream & file,
	_IOContextUCSC<TFragmentStore, TSpec> & ctx)
{
	typedef typename TFragmentStore::TContigPos         TContigPos;	
	
	unsigned suf = 0;
	if (ctx.format == ctx.KNOWN_ISOFORMS && length(ctx.transName) >= 4 && prefix(ctx.transName, 4) == "GENE")
		suf = 4;
	
	// read column 1: transcript name
	// The letters until the first whitespace will be read.
	// Then, we skip until we hit the first tab character.
	_streamWrite(file, suffix(ctx.transName, suf));
	_streamPut(file, '\t');
	
	// read column 2: contig name
	_streamWrite(file, ctx.contigName);
	if (ctx.format == ctx.KNOWN_ISOFORMS)
	{
		_streamPut(file, '\n');
		return;
	}
	_streamPut(file, '\t');
	
	// read column 3: orientation
	TContigPos transBeginPos, transEndPos;
	TContigPos cdsBeginPos, cdsEndPos;
	if (ctx.annotation.beginPos < ctx.annotation.endPos)
	{
		_streamPut(file, '+');
		transBeginPos = ctx.annotation.beginPos;
		transEndPos = ctx.annotation.endPos;
		cdsBeginPos = ctx.cdsBegin;
		cdsEndPos = ctx.cdsEnd;
	}
	else
	{
		_streamPut(file, '-');
		transEndPos = ctx.annotation.beginPos;
		transBeginPos = ctx.annotation.endPos;
		cdsEndPos = ctx.cdsBegin;
		cdsBeginPos = ctx.cdsEnd;
	}
	_streamPut(file, '\t');

	// read column 4: transcript begin position
	_streamPutInt(file, transBeginPos);
	_streamPut(file, '\t');

	// read column 5: transcript end position
	_streamPutInt(file, transEndPos);
	_streamPut(file, '\t');

	// read column 6: CDS begin position
	_streamPutInt(file, cdsBeginPos);
	_streamPut(file, '\t');
	
	// read column 7: CDS end position
	_streamPutInt(file, cdsEndPos);
	_streamPut(file, '\t');
	
	// read column 8: exon count
	_streamPutInt(file, length(ctx.exonBegin));
	_streamPut(file, '\t');

	// read column 9: exon begin positions
	for (unsigned i = 0; i < length(ctx.exonBegin); ++i)
	{
		_streamPutInt(file, _min(ctx.exonBegin[i], ctx.exonEnd[i]));
		_streamPut(file, ',');
	}
	_streamPut(file, '\t');
	
	// read column 10: exon end positions
	for (unsigned i = 0; i < length(ctx.exonBegin); ++i)
	{
		_streamPutInt(file, _max(ctx.exonBegin[i], ctx.exonEnd[i]));
		_streamPut(file, ',');
	}
	_streamPut(file, '\t');
	
	// read column 10: protein name
	_streamWrite(file, ctx.proteinName);
	_streamPut(file, '\t');

	// skip column 11
	_streamWrite(file, ctx.transName);
	_streamPut(file, '\n');
}

template<typename TTargetStream, typename TSpec, typename TConfig, typename TFormat>
inline void 
_write_KNOWN_GENE (
	TTargetStream & target,
	FragmentStore<TSpec, TConfig> & store,
	TFormat)
{
	typedef FragmentStore<TSpec, TConfig>							TFragmentStore;
	typedef typename TFragmentStore::TAnnotationStore				TAnnotationStore;
	typedef typename Value<TAnnotationStore>::Type					TAnnotation;
	typedef typename Iterator<TAnnotationStore, Standard>::Type		TAnnoIter;
	typedef typename Id<TAnnotation>::Type							TId;

	_IOContextUCSC<TFragmentStore> ctx;

	TAnnoIter it = begin(store.annotationStore, Standard());
	TAnnoIter itEnd = end(store.annotationStore, Standard());
	
	for(TId id = 0; it != itEnd; ++it, ++id)
	{
		if ((*it).typeId == TFragmentStore::ANNO_MRNA && _retrieveOneAnnotation(store, ctx, *it, id))
			_writeOneAnnotation(target, ctx);
	}
}

template<typename TTargetStream, typename TSpec, typename TConfig, typename TFormat>
inline void 
_write_KNOWN_ISOFORMS (
	TTargetStream & target,
	FragmentStore<TSpec, TConfig> & store,
	TFormat)
{
	typedef FragmentStore<TSpec, TConfig>							TFragmentStore;
	typedef typename TFragmentStore::TAnnotationStore				TAnnotationStore;
	typedef typename Value<TAnnotationStore>::Type					TAnnotation;
	typedef typename Iterator<TAnnotationStore, Standard>::Type		TAnnoIter;
	typedef typename Id<TAnnotation>::Type							TId;

	_IOContextUCSC<TFragmentStore> ctx;

	TAnnoIter it = begin(store.annotationStore, Standard());
	TAnnoIter itEnd = end(store.annotationStore, Standard());
	
	for(TId id = 0; it != itEnd; ++it, ++id)
	{
		if ((*it).typeId == TFragmentStore::ANNO_GENE && _retrieveOneAnnotation(store, ctx, *it, id))
			_writeOneAnnotation(target, ctx);
	}
}

template<typename TTargetStream, typename TSpec, typename TConfig>
inline void 
write (
	TTargetStream & target,
	FragmentStore<TSpec, TConfig> & store,
	UCSC format)
{
	_write_KNOWN_GENE(target, store, format);
}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
