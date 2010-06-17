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

#ifndef SEQAN_HEADER_STORE_IO_GFF_H
#define SEQAN_HEADER_STORE_IO_GFF_H

namespace SEQAN_NAMESPACE_MAIN
{

/**
.Tag.File Format.tag.GFF:
	GFF annotation file.
..include:seqan/store.h
*/
struct TagGFF_;
typedef Tag<TagGFF_> const GFF;

/**
.Tag.File Format.tag.GTF:
	GTF annotation file.
..include:seqan/store.h
*/
struct TagGTF_;
typedef Tag<TagGTF_> const GTF;

//////////////////////////////////////////////////////////////////////////////
// _parse_readGffIdentifier
    
    template<typename TFile, typename TString, typename TChar>
    inline void
    _parse_readGffIdentifier(TFile & file, TString & str, TChar& c)
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
	
//////////////////////////////////////////////////////////////////////////////
// skip entry until whitespace
//////////////////////////////////////////////////////////////////////////////

	template<typename TFile, typename TChar>
	inline bool
	_parse_skipEntryUntilWhitespace(TFile& file, TChar& c)
	{
		if (c== ' ' || c== '\t' || c == '\n' || (c == '\r' && _streamPeek(file) != '\n')) return false;
		
		while (!_streamEOF(file)) {
			c = _streamGet(file);
			if (c== ' ' || c== '\t' || c == '\n' || (c == '\r' && _streamPeek(file) != '\n')) break;
		}
		return true; 
	}

    template<typename TFile, typename TKeyString, typename TValueString, typename TChar>
    inline bool
    _parse_readGFFKeyValue(TFile & file, TKeyString & key, TValueString & value, TChar& c)
    {
		if (c == ' ' || c == '\t' || c == '\n' || c == '=') return false;
        append(key, c);
        while (!_streamEOF(file)) 
		{
            c = _streamGet(file);
            if (c == ' ' || c == '\t' || c == '\n' || c == '=') break;
            append(key, c);
        }
		_parse_skipSpace(file, c);
		if (c == '=')
		{
			c = _streamGet(file);
			_parse_skipSpace(file, c);
		}
		
		if (c == '"')
		{
			c = _streamGet(file);
			append(value, c);
			while (!_streamEOF(file)) 
			{
				c = _streamGet(file);
				if (c == '\n') return true;
				if (c == '"')
				{
					if (!_streamEOF(file)) c = _streamGet(file);
					break;
				}
				append(value, c);
			}
			if (c == ';')
			{
				if (!_streamEOF(file)) c = _streamGet(file);
			}
		}
		else
		{
			append(value, c);
			while (!_streamEOF(file)) 
			{
				c = _streamGet(file);
				if (c == ' ' || c == '\t' || c == '\n') return true;
				if (c == ';')
				{
					if (!_streamEOF(file)) c = _streamGet(file);
					return true;
				}
				append(value, c);
			}
		}
		return true;
	}


//////////////////////////////////////////////////////////////////////////////
// Read GFF
//////////////////////////////////////////////////////////////////////////////

template <typename TFragmentStore, typename TSpec = void>
struct _IOContextGFF
{
	typedef typename TFragmentStore::TAnnotationStore   TAnnotationStore;
	typedef typename Value<TAnnotationStore>::Type      TAnnotation;
	typedef typename TAnnotation::TId                   TId;

	CharString contigName;
	CharString typeName;
	CharString annotationName;
	CharString parentKey;
	CharString parentName;
	
	CharString _key;
	CharString _value;
	StringSet<CharString> keys;
	StringSet<CharString> values;
	
	CharString gtfGeneKey;
	CharString gtfGene;

	TId annotationId;
	TAnnotation annotation;
};

template <typename TFragmentStore, typename TSpec>
inline void clear(_IOContextGFF<TFragmentStore, TSpec> &ctx)
{
	typedef typename TFragmentStore::TAnnotationStore   TAnnotationStore;
	typedef typename Value<TAnnotationStore>::Type      TAnnotation;

	clear(ctx.contigName);
	clear(ctx.typeName);
	clear(ctx.annotationName);
	clear(ctx.parentKey);
	clear(ctx.parentName);
	clear(ctx._key);
	clear(ctx._value);
	clear(ctx.gtfGeneKey);
	clear(ctx.gtfGene);
	clear(ctx.keys);
	clear(ctx.values);
	ctx.annotationId = TAnnotation::INVALID_ID;
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
	_IOContextGFF<TFragmentStore, TSpec> & ctx)
{
	typedef typename TFragmentStore::TContigPos         TContigPos;	
	typedef typename TFragmentStore::TAnnotationStore   TAnnotationStore;
	typedef typename Value<TAnnotationStore>::Type      TAnnotation;
	typedef typename TAnnotation::TId                   TId;
		
	clear(ctx);

	// read fields of annotation line        
	_parse_skipWhitespace(file, c);
	
	// read column 1: contig name
	// The letters until the first whitespace will be read.
	// Then, we skip until we hit the first tab character.
	_parse_readGffIdentifier(file, ctx.contigName, c);
	if (!empty(ctx.contigName) && ctx.contigName[0] == '#')
	{
		_parse_skipLine(file, c);
		return false;
	}
	_parse_skipUntilChar(file, '\t', c);
	c = _streamGet(file);
	
	// skip column 2
	_parse_skipEntryUntilWhitespace(file, c);
	_parse_skipWhitespace(file, c);
	
	// read column 3: type
	_parse_readGffIdentifier(file, ctx.typeName, c);
	_parse_skipWhitespace(file, c);
	
	// read column 4: begin position
	if (_parse_isDigit(c))
		ctx.annotation.beginPos = _parse_readNumber(file, c) - 1;
	else
	{
		ctx.annotation.beginPos = TAnnotation::INVALID_POS;
		_parse_skipEntryUntilWhitespace(file, c);
	}
	_parse_skipWhitespace(file, c);

	// read column 5: end position
	if (_parse_isDigit(c))
		ctx.annotation.endPos = _parse_readNumber(file, c);
	else 
	{
		ctx.annotation.endPos = TAnnotation::INVALID_POS;
		_parse_skipEntryUntilWhitespace(file, c);
	}
	_parse_skipWhitespace(file, c);	

	// skip column 6
	_parse_skipEntryUntilWhitespace(file, c);
	_parse_skipWhitespace(file, c);

	// read column 7: orientation
	if (c == '-')
	{
		TContigPos tmp = ctx.annotation.beginPos;
		ctx.annotation.beginPos = ctx.annotation.endPos;
		ctx.annotation.endPos = tmp;
	}
	c = _streamGet(file);
	_parse_skipWhitespace(file, c);

	// skip column 8
	_parse_skipEntryUntilWhitespace(file, c);
	_parse_skipSpace(file, c);
	
	// read column 9: name
	while (!_streamEOF(file) &&	_parse_readGFFKeyValue(file, ctx._key, ctx._value, c))
	{
		if (ctx._key == "ID") 
			ctx.annotationName = ctx._value;
		else
			if (!empty(ctx._key) && !empty(ctx._value))
			{
				appendValue(ctx.keys, ctx._key);
				appendValue(ctx.values, ctx._value);
			}
				
		if (ctx._key == "Parent" || ctx._key == "ParentID" || ctx._key == "transcript_id") 
		{
			ctx.parentKey = ctx._key;
			ctx.parentName = ctx._value;
		}
		else if (ctx._key == "transcript_name")
		{
			if (empty(ctx.parentName)) 
			{
				ctx.parentKey = ctx._key;
				ctx.parentName = ctx._value;
			}
		} 
		else if (ctx._key == "gene_id")
		{
			ctx.gtfGeneKey = ctx._key;
			ctx.gtfGene = ctx._value;
		}
		else if (ctx._key == "gene_name")
		{
			if (empty(ctx.gtfGene))
			{
				ctx.gtfGeneKey = ctx._key;
				ctx.gtfGene = ctx._value;
			}
		} 

		clear(ctx._key);
		clear(ctx._value);
		_parse_skipSpace(file, c);
	}
	return true;
}

template <typename TAnnotation>
inline void 
_adjustParent (
	TAnnotation &parent,
	TAnnotation const &child)
{
	if (child.contigId == TAnnotation::INVALID_ID || child.beginPos == TAnnotation::INVALID_POS || child.endPos == TAnnotation::INVALID_POS)
		return;

	parent.contigId = child.contigId;	
	if ((parent.beginPos == TAnnotation::INVALID_POS) != (parent.endPos == TAnnotation::INVALID_POS))
		return;

	typename TAnnotation::TPos childBegin, childEnd;
	if (child.beginPos < child.endPos)
	{
		childBegin = child.beginPos;
		childEnd = child.endPos;
	} else {
		childBegin = child.endPos;
		childEnd = child.beginPos;
	}

	if (parent.beginPos < parent.endPos)
	{
		if (parent.beginPos == TAnnotation::INVALID_POS || parent.beginPos > childBegin)
			parent.beginPos = childBegin;
		if (parent.endPos == TAnnotation::INVALID_POS || parent.endPos < childEnd)
			parent.endPos = childEnd;
	} else
	{
		if (parent.endPos == TAnnotation::INVALID_POS || parent.endPos > childBegin)
			parent.endPos = childBegin;
		if (parent.beginPos == TAnnotation::INVALID_POS || parent.beginPos < childEnd)
			parent.beginPos = childEnd;
	}
}

template <typename TFragmentStore, typename TSpec>
inline void 
_storeOneAnnotation (
	TFragmentStore & fragStore,
	_IOContextGFF<TFragmentStore, TSpec> & ctx)
{
	typedef typename TFragmentStore::TAnnotationStore   TAnnotationStore;
	typedef typename Value<TAnnotationStore>::Type      TAnnotation;
	typedef typename TAnnotation::TId                   TId;
	
	TId maxId = 0;

	// for lines in GTF format get/add the parent gene first
	TId geneId = TAnnotation::INVALID_ID;
	if (!empty(ctx.gtfGene))
	{
		_storeAppendAnnotationName(fragStore, geneId, ctx.gtfGene);
		if (maxId < geneId)
			maxId = geneId;
	}	

	// if we have a parent transcript, get/add the parent transcript then
	if (!empty(ctx.parentName))
	{
		_storeAppendAnnotationName(fragStore, ctx.annotation.parentId, ctx.parentName);
		if (maxId < ctx.annotation.parentId)
			maxId = ctx.annotation.parentId;
	}
	else
		ctx.annotation.parentId = 0;	// if we have no parent, we are a child of the root

	// add contig and type name
	_storeAppendContig(fragStore, ctx.annotation.contigId, ctx.contigName);	
	_storeAppendType(fragStore, ctx.annotation.typeId, ctx.typeName);

	// add annotation name of the current line
	_storeAppendAnnotationName(fragStore, ctx.annotationId, ctx.annotationName);
	if (maxId < ctx.annotationId)
		maxId = ctx.annotationId;
	
	for (unsigned i = 0; i < length(ctx.keys); ++i)
		if (ctx.keys[i] != ctx.gtfGeneKey && ctx.keys[i] != ctx.parentKey)
			annotationAssignValueByKey(fragStore, ctx.annotation, ctx.keys[i], ctx.values[i]);

	if (length(fragStore.annotationStore) <= maxId)
		resize(fragStore.annotationStore, maxId + 1, Generous());
	fragStore.annotationStore[ctx.annotationId] = ctx.annotation;
	
	if (geneId != TAnnotation::INVALID_ID)
	{
		// link and adjust our gtf ancestors
		TAnnotation &gene = fragStore.annotationStore[geneId];
		TAnnotation &transcript = fragStore.annotationStore[ctx.annotation.parentId];

		gene.parentId = 0;
		gene.typeId = TFragmentStore::ANNO_GENE;
		_adjustParent(gene, ctx.annotation);

		transcript.parentId = geneId;
		transcript.typeId = TFragmentStore::ANNO_MRNA;
		_adjustParent(transcript, ctx.annotation);
	}
}

template<typename TFile, typename TSpec, typename TConfig>
inline void 
read (
	TFile & file,
	FragmentStore<TSpec, TConfig> & fragStore,
	GFF)
{
	typedef FragmentStore<TSpec, TConfig> TFragmentStore;
	
	if (_streamEOF(file)) return;

	// get first character from the stream
	char c = _streamGet(file);
	_IOContextGFF<TFragmentStore> ctx;
	
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

template<typename TFile, typename TSpec, typename TConfig>
inline void 
read (
	TFile & file,
	FragmentStore<TSpec, TConfig> & fragStore,
	GTF)
{
	read (file, fragStore, GFF());
}

//////////////////////////////////////////////////////////////////////////////
// Write GFF
//////////////////////////////////////////////////////////////////////////////

template<typename TTargetStream, typename TSpec, typename TConfig, typename TAnnotation, typename TId>
inline void 
_writeOneAnnotation (
	TTargetStream & target,
	FragmentStore<TSpec, TConfig> & store,
	TAnnotation &annotation,
	TId id,
	GFF)
{
	typedef FragmentStore<TSpec, TConfig>       TFragmentStore;
	typedef typename TFragmentStore::TContigPos TContigPos;
	
	if (id == 0) return;
	
	// write column 1: contig name
	if (annotation.contigId < length(store.contigNameStore))
		_streamWrite(target, store.contigNameStore[annotation.contigId]);
	_streamPut(target, '\t');
	
	// skip column 2: source
	_streamWrite(target, ".\t");
	
	// write column 3: type
	if (annotation.typeId < length(store.annotationTypeStore))
		_streamWrite(target, store.annotationTypeStore[annotation.typeId]);
	_streamPut(target, '\t');
	
	TContigPos beginPos = annotation.beginPos;
	TContigPos endPos = annotation.endPos;
	char orienation = '+';
	if (endPos < beginPos)
	{
		TContigPos tmp = beginPos;
		beginPos = endPos;
		endPos = tmp;
		orienation = '-';
	}
	
	// write column 4: begin position
	if (beginPos != TAnnotation::INVALID_POS)
		_streamPutInt(target, beginPos + 1);
	else
		_streamPut(target, '.');
	_streamPut(target, '\t');

	// write column 5: end position
	if (endPos != TAnnotation::INVALID_POS)
		_streamPutInt(target, endPos);
	else
		_streamPut(target, '.');
	_streamPut(target, '\t');

	// skip column 6: score
	_streamWrite(target, "0\t");

	// write column 7: orientation
	_streamPut(target, orienation);
	_streamPut(target, '\t');

	// skip column 8: frame
	_streamWrite(target, ".\t");
	
	// write column 9: group
	// write column 9.1: annotation id
	bool semicolon = false;
	if (id < length(store.annotationNameStore) && !empty(getAnnoName(store, id)))
	{
		_streamWrite(target, "ID=");
		_streamWrite(target, getAnnoName(store, id));
		semicolon = true;
	} 
	else if (annotation.lastChildId != TAnnotation::INVALID_ID)
	{
		_streamWrite(target, "ID=");
		_streamWrite(target, getAnnoUniqueName(store, id));
		semicolon = true;
	}
	
	// write column 9.2: parent id
	if (store.annotationStore[annotation.parentId].typeId > 1)	// ignore root/deleted nodes
	{
		if (semicolon) _streamPut(target, ';');
		_streamWrite(target, "Parent=");
		if (annotation.parentId < length(store.annotationNameStore) && !empty(getAnnoName(store, annotation.parentId)))
			_streamWrite(target, getAnnoName(store, annotation.parentId));
		else
			_streamWrite(target, getAnnoUniqueName(store, annotation.parentId));
		semicolon = true;
	}
	
	// write column 9.3-...: key, value pairs
	for (unsigned keyId = 0; keyId < length(annotation.values); ++keyId)
		if (!empty(annotation.values[keyId]))
		{
			if (semicolon) _streamPut(target, ';');
			_streamWrite(target, store.annotationKeyStore[keyId]);
			_streamPut(target, '=');
			_streamWrite(target, annotation.values[keyId]);
			semicolon = true;
		}
	
	_streamPut(target, '\n');	
}

template<typename TTargetStream, typename TSpec, typename TConfig, typename TAnnotation, typename TId>
inline void 
_writeOneAnnotation (
	TTargetStream & target,
	FragmentStore<TSpec, TConfig> & store,
	TAnnotation &annotation,
	TId id,
	GTF)
{
	typedef FragmentStore<TSpec, TConfig>				TFragmentStore;
	typedef typename TFragmentStore::TContigPos			TContigPos;
	
	if (annotation.typeId <= TFragmentStore::ANNO_MRNA) return;
	
	// write column 1: contig name
	if (annotation.contigId < length(store.contigNameStore))
		_streamWrite(target, store.contigNameStore[annotation.contigId]);
	_streamPut(target, '\t');
	
	// skip column 2: source
	_streamWrite(target, ".\t");
	
	// write column 3: type
	if (annotation.typeId < length(store.annotationTypeStore))
		_streamWrite(target, store.annotationTypeStore[annotation.typeId]);
	_streamPut(target, '\t');
	
	TContigPos beginPos = annotation.beginPos;
	TContigPos endPos = annotation.endPos;
	char orienation = '+';
	if (endPos < beginPos)
	{
		TContigPos tmp = beginPos;
		beginPos = endPos;
		endPos = tmp;
		orienation = '-';
	}
	
	// write column 4: begin position
	if (beginPos != TAnnotation::INVALID_POS)
		_streamPutInt(target, beginPos + 1);
	else
		_streamPut(target, '.');
	_streamPut(target, '\t');

	// write column 5: end position
	if (endPos != TAnnotation::INVALID_POS)
		_streamPutInt(target, endPos);
	else
		_streamPut(target, '.');
	_streamPut(target, '\t');

	// skip column 6: score
	_streamWrite(target, "0\t");

	// write column 7: orientation
	_streamPut(target, orienation);
	_streamPut(target, '\t');

	// skip column 8: frame
	_streamWrite(target, ".\t");
	
	// write column 9: group
	// write column 9.1: transcript_id
	bool semicolon = false;
	TId transcriptId = annotation.parentId;
	
	// step up until we reach a transcript
	while (transcriptId < length(store.annotationStore) && store.annotationStore[transcriptId].typeId != TFragmentStore::ANNO_MRNA)
		transcriptId = store.annotationStore[transcriptId].parentId;
	
	if (transcriptId < length(store.annotationStore))
	{
		_streamWrite(target, "transcript_id \"");
		if (transcriptId < length(store.annotationNameStore) && !empty(getAnnoName(store, transcriptId)))
			_streamWrite(target, getAnnoName(store, transcriptId));
		else
			_streamWrite(target, getAnnoUniqueName(store, transcriptId));
		_streamPut(target, '"');

		// write column 9.2: gene_id	
		TId geneId = store.annotationStore[transcriptId].parentId;
		if (geneId < length(store.annotationStore))
		{
			_streamWrite(target, "; gene_id \"");
			if (geneId < length(store.annotationNameStore) && !empty(getAnnoName(store, geneId)))
				_streamWrite(target, getAnnoName(store, geneId));
			else
				_streamWrite(target, getAnnoUniqueName(store, geneId));
			_streamPut(target, '"');
		}
		semicolon = true;
	}

	if (id < length(store.annotationNameStore) && !empty(getAnnoName(store, id)))
	{
		if (semicolon) _streamWrite(target, "; ");
		_streamWrite(target, "ID \"");
		_streamWrite(target, getAnnoName(store, id));
		_streamPut(target, '"');
		semicolon = true;
	} 
		
	// write column 9.3-...: key, value pairs
	for (unsigned keyId = 0; keyId < length(annotation.values); ++keyId)
		if (!empty(annotation.values[keyId]))
		{
			if (semicolon) _streamWrite(target, "; ");
			_streamWrite(target, store.annotationKeyStore[keyId]);
			_streamWrite(target, " \"");
			_streamWrite(target, annotation.values[keyId]);
			_streamPut(target, '"');
			semicolon = true;
		}

	if (semicolon) _streamWrite(target, ';');	
	_streamPut(target, '\n');	
}

template<typename TTargetStream, typename TSpec, typename TConfig, typename TFormat>
inline void 
_writeGFFGTF (
	TTargetStream & target,
	FragmentStore<TSpec, TConfig> & store,
	TFormat format)
{
	typedef FragmentStore<TSpec, TConfig>							TFragmentStore;
	typedef typename TFragmentStore::TAnnotationStore				TAnnotationStore;
	typedef typename Value<TAnnotationStore>::Type					TAnnotation;
	typedef typename Iterator<TAnnotationStore, Standard>::Type		TAnnoIter;
	typedef typename Id<TAnnotation>::Type							TId;

	TAnnoIter it = begin(store.annotationStore, Standard());
	TAnnoIter itEnd = end(store.annotationStore, Standard());
	
	for(TId id = 0; it != itEnd; ++it, ++id)
		_writeOneAnnotation(target, store, *it, id, format);
}

template<typename TTargetStream, typename TSpec, typename TConfig>
inline void 
write (
	TTargetStream & target,
	FragmentStore<TSpec, TConfig> & store,
	GFF format)
{
	_writeGFFGTF(target, store, format);
}

template<typename TTargetStream, typename TSpec, typename TConfig>
inline void 
write (
	TTargetStream & target,
	FragmentStore<TSpec, TConfig> & store,
	GTF format)
{
	_writeGFFGTF(target, store, format);
}


}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
