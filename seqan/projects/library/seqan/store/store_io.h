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

#ifndef SEQAN_HEADER_STORE_IO_H
#define SEQAN_HEADER_STORE_IO_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// File tags
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

/**
.Tag.File Format.tag.Amos message file:
	Amos message file.
*/
struct TagAmos_;
typedef Tag<TagAmos_> const Amos;


//////////////////////////////////////////////////////////////////////////////
// Read / Write of AMOS message files (*.afg)
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////


template<typename TFile, typename TSpec, typename TConfig>
inline void 
read(TFile & file,
	 _FragmentStore<TSpec, TConfig>& fragStore,
	 Amos) 
{
	SEQAN_CHECKPOINT
	// Basic types
	typedef _FragmentStore<TSpec, TConfig> TFragmentStore;
	typedef typename Id<TFragmentStore>::Type TId;
	typedef typename Size<TFragmentStore>::Type TSize;
	typedef typename Value<TFile>::Type TValue;

	// All fragment store element types
	typedef typename Value<typename TFragmentStore::TContigStore>::Type TContigElement;
	typedef typename Value<typename TFragmentStore::TLibraryStore>::Type TLibraryStoreElement;
	typedef typename Value<typename TFragmentStore::TMatePairStore>::Type TMatePairElement;
	typedef typename Value<typename TFragmentStore::TReadStore>::Type TReadStoreElement;
	typedef typename Value<typename TFragmentStore::TAlignedReadStore>::Type TAlignedElement;

	// All maps to mirror file ids to our ids
	typedef std::map<TId, TId> TIdMap;
	TIdMap libIdMap;
	TIdMap frgIdMap;
	TIdMap readIdMap;

	// Parse the file and convert the internal ids
	TValue c;
	if (_streamEOF(file)) return;
	else c = _streamGet(file);
	while (!_streamEOF(file)) {
		if (_streamEOF(file)) break;

		// New block?
		if (c == '{') {
			c = _streamGet(file);
			String<char> blockIdentifier;
			_parse_readIdentifier(file, blockIdentifier, c);
			_parse_skipLine(file, c);

			// Library block
			if (blockIdentifier == "LIB") {
				TLibraryStoreElement libEl;
				TId id = 0;
				String<char> fieldIdentifier;
				String<char> eid;
				while (c != '}') {
					clear(fieldIdentifier);
					_parse_readIdentifier(file, fieldIdentifier, c);
					if (fieldIdentifier == "iid") {
						c = _streamGet(file);
						id = _parse_readNumber(file, c);
						_parse_skipLine(file, c);
					} else if (fieldIdentifier == "eid") {
						c = _streamGet(file);
						while ((c != '\n') && (c != '\r')) {
							appendValue(eid, c);
							c = _streamGet(file);
						}
						_parse_skipLine(file, c);
					} else if (fieldIdentifier == "mea") {
						c = _streamGet(file);
						libEl.mean = _parse_readNumber(file, c);
						_parse_skipLine(file, c);
					} else if (fieldIdentifier == "std") {
						c = _streamGet(file);
						libEl.std = _parse_readNumber(file, c);
						_parse_skipLine(file, c);
					} else {
						_parse_skipLine(file, c);
					}
				}
				libIdMap.insert(std::make_pair(id, length(fragStore.libraryStore)));
				appendValue(fragStore.libraryStore, libEl);
				appendValue(fragStore.libraryNameStore, eid);
			} else if (blockIdentifier == "FRG") {  // Fragment block
				TMatePairElement matePairEl;
				TId id = 0;
				String<char> fieldIdentifier;
				String<char> eid;
				bool foundRds = false;
				while (c != '}') {
					clear(fieldIdentifier);
					_parse_readIdentifier(file, fieldIdentifier, c);
					if (fieldIdentifier == "iid") {
						c = _streamGet(file);
						id = _parse_readNumber(file, c);
						_parse_skipLine(file, c);
					} else if (fieldIdentifier == "eid") {
						c = _streamGet(file);
						while ((c != '\n') && (c != '\r')) {
							appendValue(eid, c);
							c = _streamGet(file);
						}
						_parse_skipLine(file, c);
					} else if (fieldIdentifier == "lib") {
						c = _streamGet(file);
						matePairEl.libId = _parse_readNumber(file, c);
						_parse_skipLine(file, c);
					} else if (fieldIdentifier == "rds") {
						foundRds = true;
						c = _streamGet(file);
						matePairEl.readId[0] = _parse_readNumber(file, c);
						c = _streamGet(file);
						matePairEl.readId[1] = _parse_readNumber(file, c);
						_parse_skipLine(file, c);
					} else {
						_parse_skipLine(file, c);
					}
				}
				// Only insert valid mate pairs
				if (foundRds) {
					frgIdMap.insert(std::make_pair(id, length(fragStore.matePairStore)));
					appendValue(fragStore.matePairStore, matePairEl);
					appendValue(fragStore.matePairNameStore, eid);
				}
			} else if (blockIdentifier == "RED") {   // Read block
				TReadStoreElement readEl;
				TId id = 0;
				String<char> fieldIdentifier;
				String<char> eid;
				String<char> qual;
				while (c != '}') {
					clear(fieldIdentifier);
					_parse_readIdentifier(file, fieldIdentifier, c);
					if (fieldIdentifier == "iid") {
						c = _streamGet(file);
						id = _parse_readNumber(file, c);
						_parse_skipLine(file, c);
					} else if (fieldIdentifier == "eid") {
						c = _streamGet(file);
						while ((c != '\n') && (c != '\r')) {
							appendValue(eid, c);
							c = _streamGet(file);
						}
						_parse_skipLine(file, c);
					} else if (fieldIdentifier == "frg") {
						c = _streamGet(file);
						readEl.matePairId = _parse_readNumber(file, c);
						_parse_skipLine(file, c);
					} else if (fieldIdentifier == "seq") {
						c = _streamGet(file);
						_parse_skipWhitespace(file, c);
						while (c != '.') {
							_parse_readSequenceData(file,c,readEl.seq);
							_parse_skipLine(file, c);
						}
					} else if (fieldIdentifier == "qlt") {
						clear(qual);
						c = _streamGet(file);
						_parse_skipWhitespace(file, c);
						while (c != '.') {
							if ((c!=' ') && (c != '\t') && (c != '\n') && (c != '\r')) appendValue(qual, c);
							c = _streamGet(file);
						}
					} else {
						_parse_skipLine(file, c);
					}
				}
				// Set quality
				typedef typename Iterator<typename TFragmentStore::TReadSeq>::Type TReadIter;
				typedef typename Iterator<String<char> >::Type TQualIter;
				TReadIter begIt = begin(readEl.seq);
				TQualIter qualIt = begin(qual);
				TQualIter qualItEnd = end(qual);
				for(;qualIt != qualItEnd; goNext(qualIt), goNext(begIt)) assignQualityValue(value(begIt), value(qualIt));

				// Insert the read
				readIdMap.insert(std::make_pair(id, length(fragStore.readStore)));
				appendValue(fragStore.readStore, readEl);
				appendValue(fragStore.readNameStore, eid);
			} else if (blockIdentifier == "CTG") {   // Contig block
				TContigElement contigEl;
				TSize fromAligned = length(fragStore.alignedReadStore);
				TId id = 0;
				String<char> fieldIdentifier;
				String<char> eid;
				String<char> contigSeq;
				String<char> contigQual;
				while (c != '}') {
					// Are we entering a TLE block
					if (c == '{') {
						TAlignedElement alignEl;
						String<char> fdIdentifier;
						typedef typename TFragmentStore::TContigPos TContigPos;
						TContigPos offsetPos = 0;
						TContigPos clr1 = 0;
						TContigPos clr2 = 0;
						String<TContigPos> gaps;
						while (c != '}') {
							clear(fdIdentifier);
							_parse_readIdentifier(file, fdIdentifier, c);
							if (fdIdentifier == "src") {
								c = _streamGet(file);
								alignEl.readId = _parse_readNumber(file, c);
								_parse_skipLine(file, c);
							} else if (fdIdentifier == "off") {
								c = _streamGet(file);
								offsetPos = _parse_readNumber(file, c);
								_parse_skipLine(file, c);
							} else if (fdIdentifier == "clr") {
								c = _streamGet(file);
								clr1 = _parse_readNumber(file, c);
								c = _streamGet(file);
								clr2 = _parse_readNumber(file, c);
								_parse_skipLine(file, c);
							} else if (fdIdentifier == "gap") {
								c = _streamGet(file);
								_parse_skipWhitespace(file, c);
								while (c != '.') {
									if ((c!=' ') && (c != '\t') && (c != '\n') && (c != '\r')) {
										TSize nextGap = _parse_readNumber(file, c);
										appendValue(gaps, nextGap);
									}
									c = _streamGet(file);
								}
							} else {
								_parse_skipLine(file, c);
							}
						}
						_parse_skipLine(file, c);

						// Get the length of the read
						TId readId = (readIdMap.find(alignEl.readId))->second;
						TSize lenRead = length((value(fragStore.readStore, readId)).seq);

						// Create the gap anchors
						typedef typename TFragmentStore::TContigGapAnchor TContigGapAnchor;
						int offset = 0;
						if ((clr1 < clr2) && (clr1>0)) offset = clr1;
						else if ((clr1 > clr2) && (clr1 < lenRead)) offset = lenRead - clr1;
						int diff = -1 * (int) (offset);
						// Clipped begin
						if (offset != 0) appendValue(alignEl.gaps, TContigGapAnchor(offset, 0));
						// Internal gaps
						typedef typename Iterator<String<TContigPos> >::Type TPosIter;
						TPosIter posIt = begin(gaps); TPosIter posItEnd = end(gaps);
						TContigPos lastGap = 0;
						TSize gapLen = 0;
						for(;posIt!=posItEnd; goNext(posIt)) {
							if (gapLen == 0) {
								++gapLen;
								++diff;
								lastGap = value(posIt);
							} 
							else if (lastGap == value(posIt)) {
								++gapLen;
								++diff;
							}
							else {
								appendValue(alignEl.gaps, TContigGapAnchor(offset + lastGap, offset + lastGap + diff));								
								gapLen = 1;
								lastGap = value(posIt);
								++diff;
							}
						}
						if (gapLen > 0) appendValue(alignEl.gaps, TContigGapAnchor(offset + lastGap, offset + lastGap + diff));
						// Clipped end
						if ((clr1 < clr2) && (clr2 < lenRead)) {
							diff -= (lenRead - clr2);				
							appendValue(alignEl.gaps, TContigGapAnchor(lenRead, lenRead + diff));

						} else if ((clr1 > clr2) && (clr2 > 0)) {
							diff -= clr2;
							appendValue(alignEl.gaps, TContigGapAnchor(lenRead, lenRead + diff));
						}
						
						// Set begin and end position
						if (clr1 < clr2) {
							alignEl.beginPos = offsetPos;
							alignEl.endPos = offsetPos + length(gaps) + (clr2 - clr1);
						} else {
							alignEl.beginPos = offsetPos + length(gaps) + (clr1 - clr2);
							alignEl.endPos = offsetPos;
						}

						// Append new align fragment, note: contigId must still be set
						appendValue(fragStore.alignedReadStore, alignEl);
					} else {
						clear(fieldIdentifier);
						_parse_readIdentifier(file, fieldIdentifier, c);
						if (fieldIdentifier == "iid") {
							c = _streamGet(file);
							id = _parse_readNumber(file, c);
							_parse_skipLine(file, c);
						} else if (fieldIdentifier == "eid") {
							c = _streamGet(file);
							while ((c != '\n') && (c != '\r')) {
								appendValue(eid, c);
								c = _streamGet(file);
							}
							_parse_skipLine(file, c);
						} else if (fieldIdentifier == "seq") {
							c = _streamGet(file);
							_parse_skipWhitespace(file, c);
							while (c != '.') {
								do {
									_parse_readSequenceData(file,c,contigSeq);
								} while (c == '-');
								_parse_skipLine(file, c);
							}
						} else if (fieldIdentifier == "qlt") {
							c = _streamGet(file);
							_parse_skipWhitespace(file, c);
							while (c != '.') {
								if ((c!=' ') && (c != '\t') && (c != '\n') && (c != '\r')) {
									appendValue(contigQual, c);
								}
								c = _streamGet(file);
							}
						} else {
							_parse_skipLine(file, c);
						}
					}
				}

				// Create the gap anchors
				char gapChar = gapValue<char>();
				typedef typename Iterator<String<char> >::Type TStringIter;
				TStringIter seqIt = begin(contigSeq);
				TStringIter seqItEnd = end(contigSeq);
				TStringIter qualIt = begin(contigQual);
				typedef typename TFragmentStore::TReadPos TPos;
				typedef typename TFragmentStore::TContigGapAnchor TContigGapAnchor;
				TPos ungappedPos = 0;
				TPos gappedPos = 0;
				bool gapOpen = false;
				for(;seqIt != seqItEnd; goNext(seqIt), goNext(qualIt), ++gappedPos) {
					if (value(seqIt) == gapChar) gapOpen = true;				
					else {
						if (gapOpen) {
							appendValue(contigEl.gaps, TContigGapAnchor(ungappedPos, gappedPos));
							gapOpen = false;
						}
						Dna5Q letter = value(seqIt);
						assignQualityValue(letter, value(qualIt));
						appendValue(contigEl.seq, letter);
						++ungappedPos;
					}
				}
				if (gapOpen) appendValue(contigEl.gaps, TContigGapAnchor(ungappedPos, gappedPos));

				// Set the contigId in all aligned reads
				TSize toAligned = length(fragStore.alignedReadStore);
				TId newContigId = length(fragStore.contigStore);
				for(; fromAligned < toAligned; ++fromAligned) {
					(value(fragStore.alignedReadStore, fromAligned)).contigId = newContigId;
				}

				// Insert the contig
				appendValue(fragStore.contigStore, contigEl);
				appendValue(fragStore.contigNameStore, eid);
			} else {
				_parse_skipLine(file, c);
			}	
		} else {
			_parse_skipLine(file, c);
		}
	}

	// Renumber all ids
	typedef typename TIdMap::const_iterator TIdMapIter;
	typedef typename Iterator<typename TFragmentStore::TMatePairStore>::Type TMateIter;
	TMateIter mateIt = begin(fragStore.matePairStore);
	TMateIter mateItEnd = end(fragStore.matePairStore);
	for(;mateIt != mateItEnd; goNext(mateIt)) {
		if (mateIt->libId != TMatePairElement::INVALID_ID) {
			TIdMapIter libIdPos = libIdMap.find(mateIt->libId);
			if (libIdPos != libIdMap.end()) mateIt->libId = libIdPos->second;
			else mateIt->libId = TMatePairElement::INVALID_ID;
		}
		if (mateIt->readId[0] != TMatePairElement::INVALID_ID) {
			TIdMapIter readIdPos = readIdMap.find(mateIt->readId[0]);
			if (readIdPos != readIdMap.end()) mateIt->readId[0] = readIdPos->second;
			else mateIt->readId[0] = TMatePairElement::INVALID_ID;
		}
		if (mateIt->readId[1]!= TMatePairElement::INVALID_ID) {
			TIdMapIter readIdPos = readIdMap.find(mateIt->readId[1]);
			if (readIdPos != readIdMap.end()) mateIt->readId[1] = readIdPos->second;
			else mateIt->readId[0] = TMatePairElement::INVALID_ID;
		}
	}
	typedef typename Iterator<typename TFragmentStore::TReadStore>::Type TReadIter;
	TReadIter readIt = begin(fragStore.readStore);
	TReadIter readItEnd = end(fragStore.readStore);
	for(;readIt != readItEnd; goNext(readIt)) {
		if (readIt->matePairId != TReadStoreElement::INVALID_ID) {
			TIdMapIter mateIdPos = frgIdMap.find(readIt->matePairId);
			if (mateIdPos != frgIdMap.end()) readIt->matePairId = mateIdPos->second;
			else readIt->matePairId = TReadStoreElement::INVALID_ID;
		}
	}
	TId myPairMatchId = 0;  // Dummy variable to count the matches
	typedef typename Iterator<typename TFragmentStore::TAlignedReadStore>::Type TAlignIter;
	TAlignIter alignIt = begin(fragStore.alignedReadStore);
	TAlignIter alignItEnd = end(fragStore.alignedReadStore);
	for(;alignIt != alignItEnd; goNext(alignIt)) {
		if (alignIt->readId != TAlignedElement::INVALID_ID) {
			TIdMapIter readIdPos = readIdMap.find(alignIt->readId);
			if (readIdPos != readIdMap.end()) alignIt->readId = readIdPos->second;
			else alignIt->readId = TAlignedElement::INVALID_ID;
		}
		alignIt->pairMatchId = myPairMatchId++;
	}
}


//////////////////////////////////////////////////////////////////////////////


template <typename TPos, typename TGapAnchor, typename TSpec>
struct _LessContigIdSort :
	public ::std::unary_function<AlignedReadStoreElement<TPos, TGapAnchor, TSpec>, bool>
{
	inline bool 
	operator() (AlignedReadStoreElement<TPos, TGapAnchor, TSpec> const& a1, AlignedReadStoreElement<TPos, TGapAnchor, TSpec> const& a2) const {
		return a1.contigId < a2.contigId;
	}
};


//////////////////////////////////////////////////////////////////////////////

template<typename TFile, typename TSpec, typename TConfig>
inline void 
write(TFile & target,
	  _FragmentStore<TSpec, TConfig>& fragStore,
	  Amos) 
{
	SEQAN_CHECKPOINT
	// Basic types
	typedef _FragmentStore<TSpec, TConfig> TFragmentStore;
	typedef typename Id<TFragmentStore>::Type TId;
	typedef typename Size<TFragmentStore>::Type TSize;
	typedef typename Value<TFile>::Type TValue;

	// All fragment store element types
	typedef typename Value<typename TFragmentStore::TContigStore>::Type TContigElement;
	typedef typename Value<typename TFragmentStore::TLibraryStore>::Type TLibraryStoreElement;
	typedef typename Value<typename TFragmentStore::TMatePairStore>::Type TMatePairElement;
	typedef typename Value<typename TFragmentStore::TReadStore>::Type TReadStoreElement;
	typedef typename Value<typename TFragmentStore::TAlignedReadStore>::Type TAlignedElement;

	// Write Header
	_streamWrite(target,"{UNV\niid:1\neid:seqan\ncom:\nafg file created with SeqAn\n.\n}\n");
	
	// Write Libraries
	typedef typename Iterator<typename TFragmentStore::TLibraryStore>::Type TLibIter;
	TLibIter libIt = begin(fragStore.libraryStore);
	TLibIter libItEnd = end(fragStore.libraryStore);
	bool noNamesPresent = (length(fragStore.libraryNameStore) == 0);
	for(TSize idCount = 0;libIt != libItEnd; goNext(libIt), ++idCount) {
		_streamWrite(target,"{LIB\n");
		_streamWrite(target,"iid:");
		_streamPutInt(target, idCount + 1);
		_streamPut(target, '\n');
		if (!noNamesPresent) {
			_streamWrite(target,"eid:");
			_streamWrite(target, value(fragStore.libraryNameStore, idCount));
			_streamPut(target, '\n');
		}
		_streamWrite(target,"{DST\n");
		_streamWrite(target,"mea:");
		_streamPutFloat(target, libIt->mean);
		_streamPut(target, '\n');
		_streamWrite(target,"std:");
		_streamPutFloat(target, libIt->std);
		_streamPut(target, '\n');
		_streamWrite(target,"}\n");	
		_streamWrite(target,"}\n");
	}

	// Write Fragments / mate pairs
	typedef typename Iterator<typename TFragmentStore::TMatePairStore>::Type TMateIter;
	TMateIter mateIt = begin(fragStore.matePairStore);
	TMateIter mateItEnd = end(fragStore.matePairStore);
	noNamesPresent = (length(fragStore.matePairNameStore) == 0);
	for(TSize idCount = 0;mateIt != mateItEnd; goNext(mateIt), ++idCount) {
		_streamWrite(target,"{FRG\n");
		_streamWrite(target,"iid:");
		_streamPutInt(target, idCount + 1);
		_streamPut(target, '\n');
		if (!noNamesPresent) {
			_streamWrite(target,"eid:");
			_streamWrite(target, value(fragStore.matePairNameStore, idCount));
			_streamPut(target, '\n');
		}
		_streamWrite(target,"lib:");
		_streamPutInt(target, mateIt->libId + 1);
		_streamPut(target, '\n');
		if ((mateIt->readId[0] != TMatePairElement::INVALID_ID) && (mateIt->readId[1] != TMatePairElement::INVALID_ID)) {
			_streamWrite(target,"rds:");
			_streamPutInt(target, mateIt->readId[0] + 1);
			_streamPut(target, ',');
			_streamPutInt(target, mateIt->readId[1] + 1);
			_streamPut(target, '\n');
		}
		_streamWrite(target,"}\n");
	}

	// Write reads
	typedef typename Iterator<typename TFragmentStore::TReadStore>::Type TReadIter;
	TReadIter readIt = begin(fragStore.readStore);
	TReadIter readItEnd = end(fragStore.readStore);
	noNamesPresent = (length(fragStore.readNameStore) == 0);
	for(TSize idCount = 0;readIt != readItEnd; goNext(readIt), ++idCount) {
		_streamWrite(target,"{RED\n");
		_streamWrite(target,"iid:");
		_streamPutInt(target, idCount + 1);
		_streamPut(target, '\n');
		if (!noNamesPresent) {
			_streamWrite(target,"eid:");
			_streamWrite(target, value(fragStore.readNameStore, idCount));
			_streamPut(target, '\n');
		}
		_streamWrite(target,"seq:\n");
		typedef typename Iterator<typename TFragmentStore::TReadSeq>::Type TSeqIter;
		typedef typename Value<typename TFragmentStore::TReadSeq>::Type TAlphabet;
		TSeqIter seqIt = begin(readIt->seq);
		TSeqIter seqItEnd = end(readIt->seq);
		for(TSize k = 0;seqIt!=seqItEnd;goNext(seqIt), ++k) {
			if ((k % 60 == 0) && (k != 0)) _streamPut(target, '\n');
			_streamPut(target, getValue(seqIt));
		}
		_streamWrite(target, "\n.\n");
		_streamWrite(target,"qlt:\n");
		seqIt = begin(readIt->seq);
		for(TSize k = 0;seqIt!=seqItEnd;goNext(seqIt), ++k) {
			if ((k % 60 == 0) && (k != 0)) _streamPut(target, '\n');
			Ascii c = ' ';
			convertQuality(c, getQualityValue(value(seqIt)));
			_streamPut(target, c);
		}
		_streamWrite(target, "\n.\n");
		if (readIt->matePairId != TReadStoreElement::INVALID_ID) {
			_streamWrite(target,"frg:");
			_streamPutInt(target, readIt->matePairId + 1);
			_streamPut(target, '\n');
		}
		_streamWrite(target,"}\n");
	}

	// Sort aligned reads according to contigId
	std::sort(begin(fragStore.alignedReadStore, Standard() ), end(fragStore.alignedReadStore, Standard() ), _LessContigIdSort<typename TFragmentStore::TContigPos, typename TFragmentStore::TReadGapAnchor, typename TFragmentStore::TAlignedReadStoreElementSpec>() );
	

	// Write Contigs
	typedef typename Iterator<typename TFragmentStore::TContigStore>::Type TContigIter;
	TContigIter contigIt = begin(fragStore.contigStore);
	TContigIter contigItEnd = end(fragStore.contigStore);
	typedef typename Iterator<typename TFragmentStore::TAlignedReadStore>::Type TAlignIter;
	TAlignIter alignIt = begin(fragStore.alignedReadStore);
	TAlignIter alignItEnd = end(fragStore.alignedReadStore);
	noNamesPresent = (length(fragStore.contigNameStore) == 0);
	for(TSize idCount = 0;contigIt != contigItEnd; goNext(contigIt), ++idCount) {
		_streamWrite(target,"{CTG\n");
		_streamWrite(target,"iid:");
		_streamPutInt(target, idCount + 1);
		_streamPut(target, '\n');
		if (!noNamesPresent) {
			_streamWrite(target,"eid:");
			_streamWrite(target, value(fragStore.contigNameStore, idCount));
			_streamPut(target, '\n');
		}
		String<char> qlt;
		_streamWrite(target,"seq:\n");
		typedef typename Iterator<typename TFragmentStore::TContigSeq>::Type TContigIter;
		TContigIter seqContigIt = begin(contigIt->seq);
		TContigIter seqContigItEnd = end(contigIt->seq);
		typedef typename Iterator<String<typename TFragmentStore::TContigGapAnchor> >::Type TGapsIter;
		TGapsIter itGaps = begin(contigIt->gaps);
		TGapsIter itGapsEnd = end(contigIt->gaps);
		int diff = 0;
		char gapChar = gapValue<char>();
		typename TFragmentStore::TContigPos mySeqPos = 0;
		TSize k = 0;
		for(;itGaps != itGapsEnd; goNext(itGaps)) {
			while (mySeqPos < itGaps->seqPos) {
				if ((k % 60 == 0) && (k != 0)) _streamPut(target, '\n');
				++k;
				_streamPut(target, value(seqContigIt));
				Ascii c = ' ';
				convertQuality(c, getQualityValue(value(seqContigIt)));
				appendValue(qlt, c);
				goNext(seqContigIt);++mySeqPos;
			}
			for(int i = 0; i < ((int) itGaps->gapPos - (int) itGaps->seqPos) - diff; ++i) {
				if ((k % 60 == 0) && (k != 0)) _streamPut(target, '\n');
				++k;
				_streamPut(target, gapChar);
				appendValue(qlt, '0');
			}
			diff = (itGaps->gapPos - itGaps->seqPos);
		}
		for(;seqContigIt != seqContigItEnd; goNext(seqContigIt)) {
			if ((k % 60 == 0) && (k != 0)) _streamPut(target, '\n');
			++k;
			_streamPut(target, value(seqContigIt));
			Ascii c = ' ';
			convertQuality(c, getQualityValue(value(seqContigIt)));
			appendValue(qlt, c);
		}
		_streamWrite(target, "\n.\n");
		_streamWrite(target,"qlt:\n");
		for(TSize k = 0;k<length(qlt); k+=60) {
			TSize endK = k + 60;
			if (endK > length(qlt)) endK = length(qlt);
			_streamWrite(target, infix(qlt, k, endK));
			_streamPut(target, '\n');
		}
		_streamWrite(target, ".\n");
		
		while ((alignIt != alignItEnd) && (idCount < alignIt->contigId)) goNext(alignIt);
		for(;(alignIt != alignItEnd) && (idCount == alignIt->contigId); goNext(alignIt)) {
			_streamWrite(target,"{TLE\n");
			_streamWrite(target,"src:");
			_streamPutInt(target, alignIt->readId + 1);
			_streamPut(target, '\n');
			typedef typename Iterator<String<typename TFragmentStore::TReadGapAnchor> >::Type TReadGapsIter;
			TReadGapsIter itGaps = begin(alignIt->gaps);
			TReadGapsIter itGapsEnd = end(alignIt->gaps);

			// Create the gaps string and the clear ranges
			typename TFragmentStore::TReadPos lenRead = length((value(fragStore.readStore, alignIt->readId)).seq);
			TSize clr1 = 0;
			TSize clr2 = lenRead;
			// Create first clear range
			if ((itGaps != itGapsEnd) && (itGaps->gapPos == 0)) clr1 = itGaps->seqPos;
			int diff = clr1;
			String<unsigned int> gaps;
			for(;itGaps != itGapsEnd; goNext(itGaps)) {
				for(int i = 0; i< diff - ((int) itGaps->seqPos - (int) itGaps->gapPos); ++i) {
					appendValue(gaps, itGaps->seqPos - clr1);
				}
				// Clipped sequence
				if (diff - ((int) itGaps->seqPos - (int) itGaps->gapPos) < 0) {
					clr2 = lenRead + diff - ((int) itGaps->seqPos - (int) itGaps->gapPos);
				}
				diff = ((int) itGaps->seqPos - (int) itGaps->gapPos);
			}
			if (alignIt->beginPos > alignIt->endPos) {
				clr1 = lenRead - clr1;
				clr2 = lenRead - clr2;
			}
			_streamWrite(target,"off:");
			if (alignIt->beginPos < alignIt->endPos) _streamPutInt(target, alignIt->beginPos);
			else _streamPutInt(target, alignIt->endPos);
			_streamPut(target, '\n');
			_streamWrite(target,"clr:");
			_streamPutInt(target, clr1);
			_streamPut(target, ',');
			_streamPutInt(target, clr2);
			_streamPut(target, '\n');
			_streamWrite(target,"gap:\n");
			for(TSize z = 0;z<length(gaps); ++z) {
				_streamPutInt(target, value(gaps, z));
				_streamPut(target, '\n');
			}
			_streamWrite(target, ".\n");
			_streamWrite(target,"}\n");
		}
		_streamWrite(target,"}\n");
	}
}


//////////////////////////////////////////////////////////////////////////////
// Read: Simple fasta read file with positions (Read simulator format)
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////


template<typename TFile, typename TSpec, typename TConfig, typename TFilePath>
inline bool 
_convertSimpleReadFile(TFile& file,
					   _FragmentStore<TSpec, TConfig>& fragStore,
					   TFilePath& filePath, 
					   bool moveToFront)
{
	SEQAN_CHECKPOINT
	// Basic types
	typedef _FragmentStore<TSpec, TConfig> TFragmentStore;
	typedef typename Id<TFragmentStore>::Type TId;
	typedef typename Size<TFragmentStore>::Type TSize;
	typedef typename Value<TFile>::Type TValue;
	typedef typename TFragmentStore::TContigPos TPos;

	// All fragment store element types
	typedef typename Value<typename TFragmentStore::TContigStore>::Type TContigElement;
	typedef typename Value<typename TFragmentStore::TLibraryStore>::Type TLibraryStoreElement;
	typedef typename Value<typename TFragmentStore::TMatePairStore>::Type TMatePairElement;
	typedef typename Value<typename TFragmentStore::TReadStore>::Type TReadStoreElement;
	typedef typename Value<typename TFragmentStore::TAlignedReadStore>::Type TAlignedElement;

	// All maps to mirror file ids to our internal ids
	typedef std::map<TId, TId> TIdMap;
	TIdMap libIdMap;
	TIdMap frgIdMap;
	TIdMap readIdMap;


	// Parse the file and convert the internal ids
	TPos maxPos = 0;
	TPos minPos = 0;
	TId count = 0;
	TValue c;
	if (_streamEOF(file)) return false;
	else c = _streamGet(file);
	while (!_streamEOF(file)) {
		if (_streamEOF(file)) break;

		// New read?
		if (c == '>') {
			TReadStoreElement readEl;
			TAlignedElement alignEl;
			TId id = count;
			TId fragId = count;
			TId repeatId = 0;
			
			c = _streamGet(file);
			_parse_skipWhitespace(file, c);

			// Get the layout positions
			alignEl.beginPos = _parse_readNumber(file, c);
			if ((count == 0) || (alignEl.beginPos < minPos)) minPos = alignEl.beginPos;
			if (alignEl.beginPos > maxPos) maxPos = alignEl.beginPos;
			c = _streamGet(file);
			_parse_skipWhitespace(file, c);
			alignEl.endPos = _parse_readNumber(file, c);
			if (alignEl.endPos < minPos) minPos = alignEl.endPos;
			if (alignEl.endPos > maxPos) maxPos = alignEl.endPos;
			
			// Any attributes?
			String<char> eid;
			String<char> qlt;
			if (c == '[') {
				String<char> fdIdentifier;
				while (c != ']') {
					c = _streamGet(file);
					_parse_skipWhitespace(file, c);
					clear(fdIdentifier);
					_parse_readIdentifier(file, fdIdentifier, c);
					if (fdIdentifier == "id") {
						c = _streamGet(file);
						id = _parse_readNumber(file, c);
					} else if (fdIdentifier == "fragId") {
						c = _streamGet(file);
						fragId = _parse_readNumber(file, c);
					} else if (fdIdentifier == "repeatId") {
						c = _streamGet(file);
						repeatId = _parse_readNumber(file, c);
					} else if (fdIdentifier == "eid") {
						c = _streamGet(file);
						while ((c != ',') && (c != ']')) {
							appendValue(eid, c);
							c = _streamGet(file);
						}
					} else if (fdIdentifier == "qlt") {
						c = _streamGet(file);
						while ((c != ',') && (c != ']')) {
							appendValue(qlt, c);
							c = _streamGet(file);
						}
					} else {
						// Jump to next attribute
						while ((c != ',') && (c != ']')) {
							c = _streamGet(file);
						}
					}
				}
			}
			_parse_skipLine(file, c);
			_parse_skipWhitespace(file, c);
			while ((!_streamEOF(file)) && (c != '>')) {
				_parse_readSequenceData(file,c,readEl.seq);
				_parse_skipWhitespace(file, c);
			}

			// Set quality
			typedef typename Iterator<typename TFragmentStore::TReadSeq>::Type TReadIter;
			typedef typename Iterator<String<char> >::Type TQualIter;
			TReadIter begIt = begin(readEl.seq);
			TReadIter begItEnd = begin(readEl.seq);
			if (length(qlt)) {
				TQualIter qualIt = begin(qlt);
				TQualIter qualItEnd = end(qlt);
				for(;qualIt != qualItEnd; goNext(qualIt), goNext(begIt)) assignQualityValue(value(begIt), value(qualIt));
			} else {
				for(;begIt != begItEnd; goNext(begIt)) assignQualityValue(value(begIt), 'D');
			}

			// Set eid if not given
			if (empty(eid)) {
				std::stringstream input;
				input << "R" << id;
				input << "-" << repeatId;
				eid = input.str().c_str();
			}

			// Set mate pair id
			readEl.matePairId = fragId;

			// Insert the read
			readIdMap.insert(std::make_pair(id, length(fragStore.readStore)));
			appendValue(fragStore.readStore, readEl);
			appendValue(fragStore.readNameStore, eid);


			// Insert an aligned read
			alignEl.readId = id;
			alignEl.pairMatchId =  fragId;
			alignEl.contigId = 0;
			appendValue(fragStore.alignedReadStore, alignEl);

			++count;
		} else {
			_parse_skipLine(file, c);
		}
	}

	// Read contig or reference sequence
	TContigElement contigEl;
	std::string fileName = filePath + 'S';
	std::fstream strmRef;
	strmRef.open(fileName.c_str(), ::std::ios_base::in | ::std::ios_base::binary);
	String<char> contigEid = "C0";
	if (!_streamEOF(strmRef)) {
		c = _streamGet(strmRef);
		while (!_streamEOF(strmRef)) {
			if (_streamEOF(strmRef)) break;
			if (c == '>') {
				clear(contigEid);
				c = _streamGet(strmRef);
				while ((c != '\r') && (c != '\n')) {
					appendValue(contigEid, c);
					c = _streamGet(strmRef);
				}
				_parse_skipLine(strmRef, c);
				_parse_skipWhitespace(strmRef, c);
				while ((!_streamEOF(strmRef)) && (c != '>')) {
					_parse_readSequenceData(strmRef,c,contigEl.seq);
					_parse_skipWhitespace(strmRef, c);
				}
			} else {
				_parse_skipLine(strmRef, c);
			}
		}
	}
	strmRef.close();
	if (empty(contigEl.seq)) {
		typedef typename TFragmentStore::TContigGapAnchor TContigGapAnchor;
		if (moveToFront) appendValue(contigEl.gaps, TContigGapAnchor(0, maxPos - minPos));
		else appendValue(contigEl.gaps, TContigGapAnchor(0, maxPos));
	}
	appendValue(fragStore.contigStore, contigEl);
	appendValue(fragStore.contigNameStore, contigEid);


	// Read fragments
	fileName = filePath + 'F';
	std::fstream strmFrag;
	strmFrag.open(fileName.c_str(), ::std::ios_base::in | ::std::ios_base::binary);
	if (!_streamEOF(strmFrag)) {
		c = _streamGet(strmFrag);
		while (!_streamEOF(strmFrag)) {
			if (_streamEOF(strmFrag)) break;
			if (c == '>') {
				TMatePairElement matePairEl;
				c = _streamGet(strmFrag);
				_parse_skipWhitespace(strmFrag, c);

				// Get the fragment id
				TId id = _parse_readNumber(strmFrag, c);
			
				// Any attributes?
				std::stringstream input;
				input << "F" << id;
				String<char> eid(input.str().c_str());
				if (c == '[') {
					String<char> fdIdentifier;
					while (c != ']') {
						c = _streamGet(strmFrag);
						_parse_skipWhitespace(strmFrag, c);
						clear(fdIdentifier);
						_parse_readIdentifier(strmFrag, fdIdentifier, c);
						if (fdIdentifier == "libId") {
							c = _streamGet(strmFrag);
							matePairEl.libId = _parse_readNumber(strmFrag, c);
						} else if (fdIdentifier == "eid") {
							clear(eid);
							c = _streamGet(strmFrag);
							while ((c != ',') && (c != ']')) {
								appendValue(eid, c);
								c = _streamGet(strmFrag);
							}
						} else {
							// Jump to next attribute
							while ((c != ',') && (c != ']')) {
								c = _streamGet(strmFrag);
							}
						}
					}
				}
				_parse_skipLine(strmFrag, c);
				_parse_skipWhitespace(strmFrag, c);

				// Read the two reads belonging to this mate pair
				matePairEl.readId[0] = _parse_readNumber(strmFrag, c);
				c = _streamGet(strmFrag);
				_parse_skipWhitespace(strmFrag, c);
				matePairEl.readId[1] = _parse_readNumber(strmFrag, c);
				_parse_skipLine(strmFrag, c);

				// Insert mate pair
				frgIdMap.insert(std::make_pair(id, length(fragStore.matePairStore)));
				appendValue(fragStore.matePairStore, matePairEl);
				appendValue(fragStore.matePairNameStore, eid);
			} else {
				_parse_skipLine(strmFrag, c);
			}
		}
	}
	strmFrag.close();

	// Read libraries
	fileName = filePath + 'L';
	std::fstream strmLib;
	strmLib.open(fileName.c_str(), ::std::ios_base::in | ::std::ios_base::binary);
	if (!_streamEOF(strmLib)) {
		c = _streamGet(strmLib);
		while (!_streamEOF(strmLib)) {
			if (_streamEOF(strmLib)) break;
			if (c == '>') {

				TLibraryStoreElement libEl;
				c = _streamGet(strmLib);
				_parse_skipWhitespace(strmLib, c);

				// Get the fragment id
				TId id = _parse_readNumber(strmLib, c);
			
				// Any attributes?
				std::stringstream input;
				input << "L" << id;
				String<char> eid(input.str().c_str());
				if (c == '[') {
					String<char> fdIdentifier;
					while (c != ']') {
						c = _streamGet(strmLib);
						_parse_skipWhitespace(strmLib, c);
						clear(fdIdentifier);
						_parse_readIdentifier(strmLib, fdIdentifier, c);
						if (fdIdentifier == "eid") {
							clear(eid);
							c = _streamGet(strmLib);
							while ((c != ',') && (c != ']')) {
								appendValue(eid, c);
								c = _streamGet(strmLib);
							}
						} else {
							// Jump to next attribute
							while ((c != ',') && (c != ']')) {
								c = _streamGet(strmLib);
							}
						}
					}
				}
				_parse_skipLine(strmLib, c);
				_parse_skipWhitespace(strmLib, c);

				// Read the mean and standard deviation
				libEl.mean = _parse_readNumber(strmLib, c);
				c = _streamGet(strmLib);
				_parse_skipWhitespace(strmLib, c);
				libEl.std = _parse_readNumber(strmLib, c);
				_parse_skipLine(strmLib, c);

				// Insert mate pair
				libIdMap.insert(std::make_pair(id, length(fragStore.libraryStore)));
				appendValue(fragStore.libraryStore, libEl);
				appendValue(fragStore.libraryNameStore, eid);
			} else {
				_parse_skipLine(strmLib, c);
			}
		}
	}
	strmLib.close();
	
	// Renumber all ids
	typedef typename TIdMap::const_iterator TIdMapIter;
	typedef typename Iterator<typename TFragmentStore::TMatePairStore>::Type TMateIter;
	TMateIter mateIt = begin(fragStore.matePairStore);
	TMateIter mateItEnd = end(fragStore.matePairStore);
	for(;mateIt != mateItEnd; goNext(mateIt)) {
		if (mateIt->libId != TMatePairElement::INVALID_ID) {
			TIdMapIter libIdPos = libIdMap.find(mateIt->libId);
			if (libIdPos != libIdMap.end()) mateIt->libId = libIdPos->second;
			else mateIt->libId = TMatePairElement::INVALID_ID;
		}
		if (mateIt->readId[0] != TMatePairElement::INVALID_ID) {
			TIdMapIter readIdPos = readIdMap.find(mateIt->readId[0]);
			if (readIdPos != readIdMap.end()) mateIt->readId[0] = readIdPos->second;
			else mateIt->readId[0] = TMatePairElement::INVALID_ID;
		}
		if (mateIt->readId[1]!= TMatePairElement::INVALID_ID) {
			TIdMapIter readIdPos = readIdMap.find(mateIt->readId[1]);
			if (readIdPos != readIdMap.end()) mateIt->readId[1] = readIdPos->second;
			else mateIt->readId[0] = TMatePairElement::INVALID_ID;
		}
	}
	typedef typename Iterator<typename TFragmentStore::TReadStore>::Type TReadIter;
	TReadIter readIt = begin(fragStore.readStore);
	TReadIter readItEnd = end(fragStore.readStore);
	for(;readIt != readItEnd; goNext(readIt)) {
		if (readIt->matePairId != TReadStoreElement::INVALID_ID) {
			TIdMapIter mateIdPos = frgIdMap.find(readIt->matePairId);
			if (mateIdPos != frgIdMap.end()) readIt->matePairId = mateIdPos->second;
			else readIt->matePairId = TReadStoreElement::INVALID_ID;
		}
	}
	typedef typename Iterator<typename TFragmentStore::TAlignedReadStore>::Type TAlignIter;
	TAlignIter alignIt = begin(fragStore.alignedReadStore);
	TAlignIter alignItEnd = end(fragStore.alignedReadStore);
	for(;alignIt != alignItEnd; goNext(alignIt)) {
		if (alignIt->readId != TAlignedElement::INVALID_ID) {
			TIdMapIter readIdPos = readIdMap.find(alignIt->readId);
			if (readIdPos != readIdMap.end()) alignIt->readId = readIdPos->second;
			else alignIt->readId = TAlignedElement::INVALID_ID;
		}
		if (moveToFront) {
			alignIt->beginPos -= minPos;
			alignIt->endPos -= minPos;
		}
	}
	return true;
}





}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
