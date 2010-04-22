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
// _appendAnnotation
// 
// adds a new entry to the read store if neccessary. Otherwise it writes the 
// correct Id in the variable using to qname to identify it
    
    template<typename TSpec, typename TConfig, typename TId, typename TName>
    inline void 
	_appendAnnotation (
		FragmentStore<TSpec, TConfig> & fragStore, 
		TId & annotationId, 
		TName & annotationName)
    {
        typedef FragmentStore<TSpec, TConfig> TFragmentStore;
        typedef typename Value<typename TFragmentStore::TContigStore>::Type TContigElement;
        
        if (!getIdByName(fragStore.annotationNameStore, annotationName, annotationId, fragStore.annotationNameStoreCache))
        {
			// if the annotation is not in the store yet
            // set the ID on the last entry after appending
            annotationId = length(fragStore.annotationNameStore);
            // append to annotationName store
			if (!empty(annotationName))
				appendName(fragStore.annotationNameStore, annotationName, fragStore.annotationNameStoreCache);
        }
    }
    
    template<typename TSpec, typename TConfig, typename TId, typename TName>
    inline void 
	_appendType (
		FragmentStore<TSpec, TConfig> & fragStore, 
		TId & typeId, 
		TName & annotationType)
    {
        typedef FragmentStore<TSpec, TConfig> TFragmentStore;
        typedef typename Value<typename TFragmentStore::TContigStore>::Type TContigElement;
        
        if (!getIdByName(fragStore.annotationTypeStore, annotationType, typeId, fragStore.annotationTypeStoreCache))
        {
			// if the annotation is not in the store yet
            // set the ID on the last entry after appending
            typeId = length(fragStore.annotationTypeStore);
            // append to annotationName store
			if (!empty(annotationType))
				appendName(fragStore.annotationTypeStore, annotationType, fragStore.annotationTypeStoreCache);
        }
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
	CharString key;
	CharString value;
	CharString annotationName;
	CharString parentName;

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
	clear(ctx.key);
	clear(ctx.value);
	clear(ctx.annotationName);
	clear(ctx.parentName);
	ctx.annotationId = TAnnotation::INVALID_ID;
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
	
	// read fields of alignments line        
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
		ctx.annotation.beginPos = TAnnotation::INVALID_POS;
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
	while (!_streamEOF(file) &&	_parse_readGFFKeyValue(file, ctx.key, ctx.value, c))
	{
		if (ctx.key == "ID") ctx.annotationName = ctx.value;
		if (ctx.key == "Parent" || ctx.key == "ParentID") ctx.parentName = ctx.value;
		clear(ctx.key);
		clear(ctx.value);
	}
	return true;
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
	
	_appendContig(fragStore, ctx.annotation.contigId, ctx.contigName);	
	_appendType(fragStore, ctx.annotation.typeId, ctx.typeName);
	if (!empty(ctx.annotationName))
		_appendAnnotation(fragStore, ctx.annotationId, ctx.annotationName);
	else
	{
		appendName(fragStore.annotationNameStore, ctx.annotationName, fragStore.annotationNameStoreCache);
		ctx.annotationId = length(fragStore.annotationNameStore) - 1;
	}
	TId maxId = ctx.annotationId;
	if (!empty(ctx.parentName))
	{
		_appendAnnotation(fragStore, ctx.annotation.parentId, ctx.parentName);
		if (maxId < ctx.annotation.parentId)
			maxId = ctx.annotation.parentId;
	}
	if (length(fragStore.annotationStore) <= maxId)
		resize(fragStore.annotationStore, maxId + 1, Generous());
	fragStore.annotationStore[ctx.annotationId] = ctx.annotation;
}

template<typename TFile, typename TSpec, typename TConfig>
inline void 
read (
	TFile & file,
	FragmentStore<TSpec, TConfig> & fragStore,
	GFF)
{
	typedef FragmentStore<TSpec, TConfig> TFragmentStore;
	
	refresh(fragStore.annotationNameStoreCache);
	refresh(fragStore.annotationTypeStoreCache);
	
	//////////////////////////////////////////////////////////////////////////////
	// build annotationNameStore:
	//////////////////////////////////////////////////////////////////////////////
	
	if (_streamEOF(file)) return;

	// get first character from the stream
	char c = _streamGet(file);
	_IOContextGFF<TFragmentStore> ctx;
	
	while (!_streamEOF(file))
	{
		if (_readOneAnnotation(file, c, ctx))
			_storeOneAnnotation(fragStore, ctx);
	}
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
	if (beginPos != TAnnotation::INVALID_POS )
		_streamPutInt(target, beginPos + 1);
	else
		_streamPut(target, '.');
	_streamPut(target, '\t');

	// write column 5: end position
	if (endPos != TAnnotation::INVALID_POS )
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
	if (id < length(store.annotationNameStore))
		if (!empty(store.annotationNameStore[id]))
		{
			_streamWrite(target, "ID=");
			_streamWrite(target, store.annotationNameStore[id]);
			_streamPut(target, ';');
		}
	
	if (annotation.parentId < length(store.annotationNameStore))
	{
		_streamWrite(target, "Parent=");
		_streamWrite(target, store.annotationNameStore[annotation.parentId]);
	}
	_streamPut(target, '\n');	
}

template<typename TTargetStream, typename TSpec, typename TConfig>
inline void 
write (
	TTargetStream & target,
	FragmentStore<TSpec, TConfig> & store,
	GFF)
{
	typedef FragmentStore<TSpec, TConfig>							TFragmentStore;
	typedef typename TFragmentStore::TAnnotationStore				TAnnotationStore;
	typedef typename Value<TAnnotationStore>::Type					TAnnotation;
	typedef typename Iterator<TAnnotationStore, Standard>::Type		TAnnoIter;
	typedef typename Id<TAnnotation>::Type							TId;

	TAnnoIter it = begin(store.annotationStore, Standard());
	TAnnoIter itEnd = end(store.annotationStore, Standard());
	
	for(TId id = 0; it != itEnd; ++it, ++id)
		_writeOneAnnotation(target, store, *it, id, GFF());
}


}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
