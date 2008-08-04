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
  $Id: graph_align_tcoffee_io.h 1919 2008-05-02 15:54:46Z rausch@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

#ifndef SEQAN_HEADER_GRAPH_ALIGN_TCOFFEE_IO_H
#define SEQAN_HEADER_GRAPH_ALIGN_TCOFFEE_IO_H

namespace SEQAN_NAMESPACE_MAIN
{

/////////////////////////////////////////////////////////////////////////////
// Input and Output of an alignment graph, Tree Reading im Newick Format
//////////////////////////////////////////////////////////////////////////////



/////////////////////////////////////////////////////////////////////////////
// Tags
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

/**
.Tag.Alignment Graph Format.value.TCoffeeLib:
	T-Coffee library format to read and write an alignment graph.
*/

struct TCoffeeLib_;
typedef Tag<TCoffeeLib_> const TCoffeeLib;

//////////////////////////////////////////////////////////////////////////////

/**
.Tag.Alignment Graph Format.value.BlastLib:
	A blast library for matches for an alignment graph.
*/

struct BlastLib_;
typedef Tag<BlastLib_> const BlastLib;

//////////////////////////////////////////////////////////////////////////////

/**
.Tag.Alignment Graph Format.value.MummerLib:
	A mummer library for matches for an alignment graph.
*/

struct MummerLib_;
typedef Tag<MummerLib_> const MummerLib;


//////////////////////////////////////////////////////////////////////////////

/**
.Tag.Alignment Graph Format.value.NewickFormat:
	NewickFormat format to write a guide tree.
*/

struct NewickFormat_;
typedef Tag<NewickFormat_> const NewickFormat;



/////////////////////////////////////////////////////////////////////////////
// T-Coffee Library Reading / Writing
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template<typename TFile, typename TChar, typename TString>
inline void
_parse_readSequenceData(TFile & file,
						TChar & c,
						TString& str)
{
	SEQAN_CHECKPOINT

	append(str, c);

	// Read word
	while (!_streamEOF(file)) {
		c = _streamGet(file);
		if (!_parse_isLetter(c)) break;
		append(str, c);
	}
}


template<typename TFile, typename TStringSet, typename TCargo, typename TSpec>
inline void
_readLibrary(TFile & file,
			 Graph<Alignment<TStringSet, TCargo, TSpec> >& g)
{
	SEQAN_CHECKPOINT
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename Size<TGraph>::Type TSize;
	typedef typename Value<TFile>::Type TValue;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename Id<TGraph>::Type TId;

	TValue c = '#';
	bool seq1ToN = false;
	if (_streamEOF(file)) return;
	
	typedef std::pair<std::pair<TSize, TSize>, TCargo> TResiduePair;
	typedef std::set<TResiduePair> TResiduePairSet;
	TSize nseq = length(stringSet(g));
	String<TResiduePairSet> resPair;
	resize(resPair, nseq * nseq);	
	TSize seq1 = 0;
	TSize seq2 = 0;
	bool firstPass = true;
	while (!_streamEOF(file)) {
		_parse_skipWhitespace(file,c);
		if (_streamEOF(file)) break;
		if (c == '#') {
			c = _streamGet(file);
			_parse_skipWhitespace(file,c);
			seq1 = _parse_readNumber(file, c);
			_parse_skipWhitespace(file,c);
			seq2 = _parse_readNumber(file, c);
			if (firstPass) {
				firstPass = false;
				if ((seq1 != 0) && (seq2 != 0)) seq1ToN = true;
			}
			if (seq1ToN) {
				--seq1;
				--seq2;
			}
		} else if (c == '!') {
			_parse_skipLine(file, c);
		} else {
			TSize res1 = _parse_readNumber(file, c);
			_parse_skipWhitespace(file,c);
			TSize res2 = _parse_readNumber(file, c);
			_parse_skipWhitespace(file,c);
			TSize weight = _parse_readNumber(file, c);
			_parse_skipLine(file,c);

			if (seq1 < seq2) {
				TSize index = seq1 * nseq + seq2;
				resPair[index].insert(std::make_pair(std::make_pair(--res1,--res2), weight));
			} else {
				TSize index = seq2 * nseq + seq1;
				resPair[index].insert(std::make_pair(std::make_pair(--res2,--res1), weight));
			}		
		}
	}

	typedef Fragment<> TFragment;
	typedef String<TFragment> TFragmentString;
	typedef typename Iterator<TFragmentString>::Type TFragmentStringIter;
	TFragmentString matches;
	String<TCargo> score_values;

	for(TSize i = 0; i<length(resPair); ++i) {
		if (resPair[i].empty()) continue;
		TSize seq1 = i / nseq;
		TSize seq2 = i % nseq;
		//std::cout << "#" << seq1 << ',' << seq2 << std::endl;
		typename TResiduePairSet::const_iterator pos = resPair[i].begin();
		typename TResiduePairSet::const_iterator posEnd = resPair[i].end();
		TSize startMatch1 = pos->first.first;
		TSize startMatch2 = pos->first.second;
		TCargo carg = pos->second;
		TSize len = 1;
		//std::cout << pos->first.first << ',' << pos->first.second << ',' << pos->second << std::endl;
		++pos;
		while(pos != posEnd) {
			if ((startMatch1 + len == pos->first.first) &&
				(startMatch2 + len == pos->first.second) &&
				(carg / len == pos->second)) {
					carg += pos->second;
					++len;
			} else {
				appendValue(matches, TFragment(seq1, startMatch1, seq2, startMatch2, len));
				appendValue(score_values, carg);
				startMatch1 = pos->first.first;
				startMatch2 = pos->first.second;
				carg = pos->second;
				len = 1;
			}
			//std::cout << pos->first.first << ',' << pos->first.second << ',' << pos->second << std::endl;
			++pos;
		}
		appendValue(matches, TFragment(seq1, startMatch1, seq2, startMatch2, len));
		appendValue(score_values, carg);
	}

	// Refine all matches, rescore matches and create multiple alignment
	matchRefinement(matches,stringSet(g),g);

	// Adapt edge weights
	TFragmentStringIter endIt = end(matches);
	TSize positionIt = 0;
	for(TFragmentStringIter it = begin(matches); it != endIt; ++it, ++positionIt) {
		TId id1 = sequenceId(*it,0);
		TId id2 = sequenceId(*it,1);
		TSize pos1 = fragmentBegin(*it, id1);
		TSize pos2 = fragmentBegin(*it, id2);
		TSize origFragLen = fragmentLength(*it, id1);
		TSize end1 = pos1 + origFragLen;
		while(pos1 < end1) {
			TVertexDescriptor p1 = findVertex(g, id1, pos1);
			TVertexDescriptor p2 = findVertex(g, id2, pos2);
			TSize fragLen = fragmentLength(g, p1);
			TEdgeDescriptor e = findEdge(g, p1, p2);
			cargo(e) += (TCargo) ( (double) fragLen / (double) origFragLen * (double) getValue(score_values, positionIt));
			pos1 += fragLen;
			pos2 += fragLen;
		}
	}
}

//////////////////////////////////////////////////////////////////////////////

template<typename TFile, typename TStringSet, typename TCargo, typename TSpec>
void 
read(TFile & file,
	 Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
	 TCoffeeLib) 
{
	SEQAN_CHECKPOINT
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename Size<TGraph>::Type TSize;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename Position<TFile>::Type TPosition;
	typedef typename Value<TFile>::Type TValue;
	typedef typename Id<TGraph>::Type TIdType;
	typedef typename Value<TStringSet>::Type TString;

	TValue c;
	if (_streamEOF(file)) return;
	else c = _streamGet(file);

	// Ignore first line
	_parse_skipLine(file, c);
		
	// Read number of sequences
	TSize nSeq = (TSize) _parse_readNumber(file, c);
	_parse_skipLine(file, c);
	
	// Skip header
	for(TSize i=0; i<nSeq; ++i) _parse_skipLine(file, c);

	// Read library
	clearVertices(g);
	_readLibrary(file,g);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TFile, typename TString, typename TSpec, typename TNames>
void 
read(TFile & file,
	 StringSet<TString, TSpec>& oriStr,
	 TNames& names,
	 TCoffeeLib) 
{
	SEQAN_CHECKPOINT
	typedef typename Size<TNames>::Type TSize;
	typedef typename Position<TFile>::Type TPosition;
	typedef typename Value<TFile>::Type TValue;
	
	TValue c;
	if (_streamEOF(file)) return;
	else c = _streamGet(file);

	// Ignore first line
	_parse_skipLine(file, c);
	
	// Read number of sequences
	TSize nSeq = (TSize) _parse_readNumber(file, c);
	resize(oriStr, nSeq);
	_parse_skipLine(file, c);

	// Read sequences
	for(TSize i=0; i<nSeq; ++i) {
		appendValue(names, _parse_readIdentifier(file, c));
		_parse_skipWhitespace(file, c);
		_parse_readNumber(file, c);
		_parse_skipWhitespace(file, c);
		_parse_readSequenceData(file,c,oriStr[i]);
		_parse_skipLine(file, c);
	}
}

//////////////////////////////////////////////////////////////////////////////

template<typename TFile, typename TStringSet, typename TCargo, typename TSpec, typename TNames>
void write(TFile & file, 
		   Graph<Alignment<TStringSet, TCargo, TSpec> > const& g,
		   TNames& names,
		   TCoffeeLib) 
{
	SEQAN_CHECKPOINT
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef typename Value<TStringSet>::Type TString;
	typedef typename Size<TStringSet>::Type TSize;
	
	typedef std::pair<std::pair<TSize, TSize>, TCargo> TResiduePair;
	typedef std::set<TResiduePair> TResiduePairSet;
	TSize nseq = length(stringSet(g));
	String<TResiduePairSet> resPair;
	resize(resPair, nseq * nseq);
	
	typedef typename Iterator<TGraph, EdgeIterator>::Type TIter;
	TIter it(g);
	for(;!atEnd(it);++it) {
		TVertexDescriptor sV = sourceVertex(it);
		TVertexDescriptor tV = targetVertex(it);
		TSize fragLen = fragmentLength(g,sV);
		TSize fragPos1 = fragmentBegin(g,sV);
		TSize fragPos2 = fragmentBegin(g,tV);
		TSize seq1 = sequenceId(g,sV);
		TSize seq2 = sequenceId(g,tV);
		TCargo my_carg =  getCargo(*it);
		if (my_carg <= 0) my_carg = 1;
		else my_carg = (TCargo) ((double) my_carg / (double) fragLen);
		for(TSize i = 0; i<fragLen; ++i) {
			resPair[seq1 * nseq + seq2].insert(std::make_pair(std::make_pair(fragPos1 + i, fragPos2 + i), my_carg));
		}
	}

	_streamWrite(file, "! TC_LIB_FORMAT_01\n");
	TSize len = length(getStringSet(g));
	_streamPutInt(file, len);
	_streamPut(file, '\n');
	for(TSize i=0;i<len;++i) {
		_streamWrite(file, names[i]);
		_streamPut(file, ' ');
		TString str = value(getStringSet(g), i);
		_streamPutInt(file, length(str));
		_streamPut(file, ' ');
		_streamWrite(file, str);
		_streamPut(file, '\n');
	}

	for(TSize i = 0; i<length(resPair); ++i) {
		if (resPair[i].empty()) continue;
		TSize seq1 = i / nseq;
		TSize seq2 = i % nseq;
		_streamPut(file, '#');
		_streamPutInt(file, seq1 + 1);
		_streamPut(file, ' ');
		_streamPutInt(file, seq2 + 1);
		_streamPut(file, '\n');	
		typename TResiduePairSet::const_iterator pos = resPair[i].begin();
		typename TResiduePairSet::const_iterator posEnd = resPair[i].end();
		while(pos != posEnd) {
			_streamPutInt(file, pos->first.first + 1);
			_streamPut(file, ' ');
			_streamPutInt(file, pos->first.second + 1);
			_streamPut(file, ' ');
			_streamPutInt(file, pos->second);
			_streamPut(file, '\n');	
			++pos;
		}
	}
	_streamWrite(file, "! SEQ_1_TO_N");
	_streamPut(file, '\n');
}


/////////////////////////////////////////////////////////////////////////////
// FastaAlign Reading
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template<typename TFile, typename TString, typename TSpec, typename TNames>
void 
read(TFile & file,
	 StringSet<TString, TSpec>& oriStr,
	 TNames& names,
	 FastaAlign) 
{
	SEQAN_CHECKPOINT
	typedef typename Size<TNames>::Type TSize;
	typedef TSize TWord;
	typedef typename Value<TFile>::Type TValue;
	
	TValue c;
	if (_streamEOF(file)) return;
	else c = _streamGet(file);

	// Read sequences
	TString seq;
	while(!_streamEOF(file)) {
		_parse_skipWhitespace(file, c);
		if (_streamEOF(file)) break;
		if (c == '>') {
			if (length(seq)) {
				appendValue(oriStr, seq);
				clear(seq);
			}
			c = _streamGet(file);
			appendValue(names, _parse_readIdentifier(file, c));
			_parse_skipLine(file, c);
		} else if ((c == '-') || (c == '.') || (c == '\n') || (c == '\r')) {
			c = _streamGet(file);
		} else {
			appendValue(seq, c);
			c = _streamGet(file);
		}
	}
	if (length(seq)) appendValue(oriStr, seq);
}


template<typename TValue, typename TScore, typename TStringSet, typename TCargo, typename TSpec>
void 
_buildAlignmentGraph(String<TValue> const& mat,
					 TScore const& scType,
					 Graph<Alignment<TStringSet, TCargo, TSpec> >& g) 
{
	SEQAN_CHECKPOINT
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename Size<TGraph>::Type TSize;
	TSize nseq = length(stringSet(g));
	TSize len = length(mat) / nseq;
	TValue gapChar = gapValue<TValue>();

	// Create the anchor graph
	typedef Fragment<> TFragment;
	typedef String<TFragment> TFragmentString;
	typedef typename Iterator<TFragmentString>::Type TFragmentStringIter;
	TFragmentString matches;

	typedef std::pair<TSize, TSize> TResiduePair;
	typedef std::set<TResiduePair> TResiduePairSet;
	String<TResiduePairSet> resPair;
	resize(resPair, nseq * nseq);	
	for(TSize seq1 = 0; seq1 < nseq - 1; ++seq1) {
		for(TSize seq2 = seq1 + 1; seq2 < nseq; ++seq2) {
			TSize index = seq1 * nseq + seq2;
			TSize offset1 = 0;
			TSize offset2 = 0;
			for(TSize col = 0; col<len; ++col) {
				if (value(mat, seq1 * len + col) != gapChar) {
					if (value(mat, seq2 * len + col) != gapChar) {
						resPair[index].insert(std::make_pair(offset1, offset2));
						++offset1;
						++offset2;
					} else ++offset1;
				} else if (value(mat, seq2 * len + col) != gapChar) ++offset2;
			}
		}
	}
	for(TSize i = 0; i<length(resPair); ++i) {
		if (resPair[i].empty()) continue;
		TSize seq1 = i / nseq;
		TSize seq2 = i % nseq;
		typename TResiduePairSet::const_iterator pos = resPair[i].begin();
		typename TResiduePairSet::const_iterator posEnd = resPair[i].end();
		TSize startMatch1 = pos->first;
		TSize startMatch2 = pos->second;
		TSize len = 1;
		++pos;
		while(pos != posEnd) {
			if ((startMatch1 + len == pos->first) && (startMatch2 + len == pos->second)) ++len;
			else {
				appendValue(matches, TFragment(seq1, startMatch1, seq2, startMatch2, len));
				startMatch1 = pos->first;
				startMatch2 = pos->second;
				len = 1;
			}
			++pos;
		}
		appendValue(matches, TFragment(seq1, startMatch1, seq2, startMatch2, len));
	}
	clearVertices(g);
	matchRefinement(matches,stringSet(g),const_cast<TScore&>(scType),g);
}


//////////////////////////////////////////////////////////////////////////////

template<typename TFile, typename TStringSet, typename TCargo, typename TSpec, typename TNames, typename TScoreType>
void 
read(TFile & file,
	 Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
	 TNames const& origNames,
	 TScoreType const& scType,
	 FastaAlign) 
{
	SEQAN_CHECKPOINT
	typedef typename Size<TNames>::Type TSize;
	typedef typename Value<TFile>::Type TValue;
	
	TValue c;
	if (_streamEOF(file)) return;
	else c = _streamGet(file);

	// Read sequences
	String<TValue> mat;
	TNames names;
	while(!_streamEOF(file)) {
		_parse_skipWhitespace(file, c);
		if (_streamEOF(file)) break;
		if (c == '>') {
			c = _streamGet(file);
			appendValue(names, _parse_readIdentifier(file, c));
			_parse_skipLine(file, c);
		} else if ((c == '\n') || (c == '\r')) {
			c = _streamGet(file);
		} else {
			appendValue(mat, c);
			c = _streamGet(file);
		}
	}
	// Reorder rows according to names order
	String<TValue> finalMat = mat;
	TSize nseq = length(stringSet(g));
	TSize len = length(mat) / nseq;
	for(TSize i = 0; i<length(names); ++i) {
		if (value(names, i) == value(origNames, i)) continue;
		else {
			for(TSize j = 0; j<length(origNames); ++j) {
				if (value(names, i) != value(origNames, j)) continue;
				for(TSize k = 0; k<len; ++k) {
					value(finalMat, j * len + k) = value(mat, i * len + k);
				}
			}
		}
	}
	clear(mat);
	_buildAlignmentGraph(finalMat, scType, g); 
}

/////////////////////////////////////////////////////////////////////////////
// BLAST Library Reading / Writing
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template<typename TSizeSpec, typename TSpec1, typename TSpec2, typename TSize>
inline void 
__includeFragment(String<Fragment<TSizeSpec, ExactReversableFragment<TSpec1> >, TSpec2>& matches, 
				  TSize seq1Id, 
				  TSize beg1, 
				  TSize seq2Id, 
				  TSize beg2, 
				  TSize len, 
				  bool reversed)
{
	SEQAN_CHECKPOINT
	typedef Fragment<TSizeSpec, ExactReversableFragment<TSpec1> > TFragment;
	appendValue(matches, TFragment(seq1Id, beg1, seq2Id, beg2, len, reversed));
}

//////////////////////////////////////////////////////////////////////////////

template<typename TSizeSpec, typename TSpec1, typename TSpec2, typename TSize>
inline void 
__includeFragment(String<Fragment<TSizeSpec, ExactFragment<TSpec1> >, TSpec2>& matches, 
				  TSize seq1Id, 
				  TSize beg1, 
				  TSize seq2Id, 
				  TSize beg2, 
				  TSize len, 
				  bool)
{
	SEQAN_CHECKPOINT
	typedef Fragment<TSizeSpec, ExactFragment<TSpec1> > TFragment;
	appendValue(matches, TFragment(seq1Id, beg1, seq2Id, beg2, len));
}

//////////////////////////////////////////////////////////////////////////////

template<typename TFile, typename TFragment, typename TSpec1, typename TScoreValue, typename TSpec2, typename TNames>
inline void 
read(TFile & file,
	 String<TFragment, TSpec1>& matches,
	 String<TScoreValue, TSpec2>& scores,
	 TNames& names,
	 BlastLib) 
{
	SEQAN_CHECKPOINT
	typedef typename Size<TNames>::Type TSize;
	typedef typename Value<TFile>::Type TValue;
	typedef typename Value<TNames>::Type TName;
	typedef String<TFragment, TSpec1> TFragmentString;
	typedef typename Iterator<TFragmentString>::Type TFragmentStringIter;
	
	// Map the names to slots
	typedef std::map<TName, TSize> TNameToPosition;
	TNameToPosition namePosMap;
	for(TSize i = 0;i<length(names);++i) namePosMap.insert(std::make_pair(names[i], i));
	
	// Read the Blast file
	TValue c;
	if (_streamEOF(file)) return;
	else c = _streamGet(file);

	TName seq1;
	TName seq2;
	while (!_streamEOF(file)) {
		clear(seq1);
		clear(seq2);
		_parse_skipWhitespace(file, c);
		seq1 = _parse_readIdentifier(file, c);
		_parse_skipWhitespace(file, c);
		seq2 = _parse_readIdentifier(file, c);
		if (seq1 == seq2) {
			_parse_skipLine(file, c);
			continue;
		}
		TSize seq1Id = namePosMap[seq1];
		TSize seq2Id = namePosMap[seq2];
		_parse_skipWhitespace(file, c);
		_parse_readDouble(file, c);
		_parse_skipWhitespace(file, c);
		_parse_readNumber(file, c);
		_parse_skipWhitespace(file, c);
		_parse_readNumber(file, c);
		_parse_skipWhitespace(file, c);
		_parse_readNumber(file, c);
		_parse_skipWhitespace(file, c);
		TSize beg1 = _parse_readNumber(file, c);
		_parse_skipWhitespace(file, c);
		TSize end1 = _parse_readNumber(file, c);
		TSize len = end1 - beg1 + 1;
		_parse_skipWhitespace(file, c);
		TSize beg2 = _parse_readNumber(file, c);
		_parse_skipWhitespace(file, c);
		TSize end2 = _parse_readNumber(file, c);
		_parse_skipWhitespace(file, c);
		_parse_readIdentifier(file, c);
		_parse_skipWhitespace(file, c);
		_parse_readDouble(file, c);

		bool reversed = false;
		if (beg1 > end1) { TSize tmp = beg1; beg1 = end1; end1 = tmp; reversed = !reversed; }
		if (beg2 > end2) { TSize tmp = beg2; beg2 = end2; end2 = tmp; reversed = !reversed; }
		//// Debug code
		//std::cout << seq1Id << ',' << beg1 << ',' << seq2Id << ',' << beg2 << ',' << len << std::endl;
		//std::cout << infix(strSet[seq1Id], beg1, beg1+len) << std::endl;
		//std::cout << infix(strSet[seq2Id], beg2, beg2+len) << std::endl;

		__includeFragment(matches, seq1Id, --beg1, seq2Id, --beg2, len, reversed);
		appendValue(scores, (TScoreValue) len);
		_parse_skipLine(file, c);
	}
}

//////////////////////////////////////////////////////////////////////////////

template<typename TFile, typename TStringSet, typename TCargo, typename TSpec, typename TNames, typename TEdgeMap>
void write(TFile & file, 
		   Graph<Alignment<TStringSet, TCargo, TSpec> > const& g,
		   TNames& names,
		   TEdgeMap& edgeMap,
		   BlastLib) 
{
	SEQAN_CHECKPOINT
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef typename Value<TStringSet>::Type TString;
	typedef typename Size<TStringSet>::Type TSize;
	
	TStringSet& str = stringSet(g);
	
	typedef typename Iterator<TGraph, EdgeIterator>::Type TIter;
	TIter it(g);
	for(;!atEnd(it);++it) {
		TVertexDescriptor sV = sourceVertex(it);
		TVertexDescriptor tV = targetVertex(it);
		if (sequenceId(g,sV) < sequenceId(g,tV)) {
			TVertexDescriptor tmp = sV;
			sV = tV;
			tV = tmp;
		}
		TSize fragLen = fragmentLength(g,sV);
		TSize fragPos1 = fragmentBegin(g,sV);
		TSize fragPos2 = fragmentBegin(g,tV);
		TSize seq1 = idToPosition(str, sequenceId(g,sV));
		TSize seq2 = idToPosition(str, sequenceId(g,tV));
		TCargo my_carg =  getCargo(*it);
		_streamWrite(file, names[seq1]);
		_streamPut(file, '\t');	
		_streamWrite(file, names[seq2]);
		_streamPut(file, '\t');	
		_streamPutInt(file, 0);
		_streamPut(file, '\t');	
		_streamPutInt(file, fragLen);
		_streamPut(file, '\t');	
		_streamPutInt(file, 0);
		_streamPut(file, '\t');	
		_streamPutInt(file, 0);
		_streamPut(file, '\t');	
		_streamPutInt(file, fragPos1+1);
		_streamPut(file, '\t');	
		_streamPutInt(file, fragPos1 + fragLen);
		_streamPut(file, '\t');	
		if (!property(edgeMap, *it)) {
			_streamPutInt(file, fragPos2+1);
			_streamPut(file, '\t');	
			_streamPutInt(file, fragPos2 + fragLen);
			_streamPut(file, '\t');	
		} else {
			_streamPutInt(file, fragPos2 + fragLen);
			_streamPut(file, '\t');		
			_streamPutInt(file, fragPos2+1);
			_streamPut(file, '\t');	
		}
		_streamPutInt(file, 0);
		_streamPut(file, '\t');	
		_streamPutInt(file, my_carg);
		_streamPut(file, '\n');
	}
}

//////////////////////////////////////////////////////////////////////////////

template<typename TFile, typename TStringSet, typename TCargo, typename TSpec, typename TNames>
void write(TFile & file, 
		   Graph<Alignment<TStringSet, TCargo, TSpec> > const& g,
		   TNames& names,
		   BlastLib) 
{
	SEQAN_CHECKPOINT
	String<bool> edgeMap;
	fill(edgeMap, getIdUpperBound(_getEdgeIdManager(g)), false);
	write(file, g, names, edgeMap, BlastLib());
}



/////////////////////////////////////////////////////////////////////////////
// MUMMER Format Reading
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////

template<typename TFile, typename TStringSet, typename TCargo, typename TSpec, typename TNames, typename TEdgeMap>
void 
read(TFile & file,
	 Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
	 TNames& names, 
	 TEdgeMap& edgeMap,
	 MummerLib) 
{
	SEQAN_CHECKPOINT
	typedef typename Size<TNames>::Type TSize;
	typedef typename Value<TFile>::Type TValue;
	typedef typename Value<TNames>::Type TName;
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename Id<TGraph>::Type TId;
	
	// Initialization
	typedef Fragment<TSize, ExactReversableFragment<> > TFragment;
	typedef String<TFragment> TFragmentString;
	typedef typename Iterator<TFragmentString>::Type TFragmentStringIter;
	TFragmentString matches;
	String<TCargo> score_values;

	// Map the names to slots
	typedef std::map<TName, TSize> TNameToPosition;
	TNameToPosition namePosMap;
	for(TSize i = 0;i<length(names);++i) {
		namePosMap.insert(std::make_pair(names[i], i));
	}
	
	// Read the Mummer file
	TValue c;
	if (_streamEOF(file)) return;
	else c = _streamGet(file);
	TStringSet& strSet = stringSet(g);

	TName seq1;
	TName seq2;
	TSize seq1Id = 0;
	TSize seq2Id = 0;
	bool reversed = false;
	while (!_streamEOF(file)) {
		if (c == '>') {
			c = _streamGet(file);
			_parse_skipWhitespace(file, c);
			seq1 = _parse_readIdentifier(file, c);
			 seq1Id = namePosMap[seq1];
			 _parse_skipWhitespace(file, c);
			 if (c == 'R') {
				 reversed = true;
				 _parse_skipLine(file, c);
			 } else reversed = false;
		} else {
			_parse_skipWhitespace(file, c);
			if (_streamEOF(file)) {
				break;
			}
			seq2 = _parse_readIdentifier(file, c);
			seq2Id = namePosMap[seq2];
			_parse_skipWhitespace(file, c);
			TSize beg2 = _parse_readNumber(file, c);
			_parse_skipWhitespace(file, c);
			TSize beg1 = _parse_readNumber(file, c);
			_parse_skipWhitespace(file, c);
			TSize len = _parse_readNumber(file, c);
			_parse_skipLine(file, c);
			if (seq1Id == seq2Id) continue;
			
			if (!reversed) appendValue(matches, TFragment(seq1Id, --beg1, seq2Id,  --beg2, len, reversed));
			else appendValue(matches, TFragment(seq1Id, (length(value(strSet, seq1Id)) - (--beg1 + len)), seq2Id,  --beg2, len, reversed));
			appendValue(score_values, (TCargo) len);
		}
	}

	// Refine all matches, rescore matches and create multiple alignment
	matchRefinement(matches,strSet,g);
	resizeEdgeMap(g, edgeMap);

	// Adapt edge weights
	TFragmentStringIter endIt = end(matches);
	TSize positionIt = 0;
	for(TFragmentStringIter it = begin(matches); it != endIt; ++it, ++positionIt) {
		TId id1 = sequenceId(*it,0);
		TSize pos1 = fragmentBegin(*it, id1);
		TSize origFragLen = fragmentLength(*it, id1);
		TSize end1 = pos1 + origFragLen;
		bool reversed = isReversed(*it);
		while(pos1 < end1) {
			TSize pos2 = 0;
			TId id2 = 0;
			getProjectedPosition(*it, id1, pos1, id2, pos2);
			TVertexDescriptor p1 = findVertex(g, id1, pos1);
			TVertexDescriptor p2 = findVertex(g, id2, pos2);
			TSize fragLen = fragmentLength(g, p1);
			TEdgeDescriptor e = findEdge(g, p1, p2);
			cargo(e) += (TCargo) ( ((double) fragLen / (double) origFragLen) * (double) value(score_values, positionIt));
			property(edgeMap, e) = reversed;
			pos1 += fragLen;
			pos2 += fragLen;
		}
	}
}

//////////////////////////////////////////////////////////////////////////////

template<typename TFile, typename TStringSet, typename TCargo, typename TSpec, typename TNames>
inline void 
read(TFile & file,
	 Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
	 TNames& names, 
	 MummerLib) 
{
	String<bool> edgeMap;
	read(file, g, names, edgeMap, MummerLib() );
}


/////////////////////////////////////////////////////////////////////////////
// Newick Format
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template<typename TFile, typename TCargo, typename TSpec, typename TNames>
void 
read(TFile & file,
	 Graph<Tree<TCargo, TSpec> >& guideTree,
	 TNames& names,
	 NewickFormat) 
{
	SEQAN_CHECKPOINT
	typedef Graph<Tree<TCargo, TSpec> > TGuideTree;
	typedef typename VertexDescriptor<TGuideTree>::Type TVertexDescriptor;
	typedef typename EdgeDescriptor<TGuideTree>::Type TEdgeDescriptor;
	typedef typename Size<TGuideTree>::Type TSize;
	typedef typename Id<TGuideTree>::Type TId;
	typedef typename Position<TFile>::Type TPosition;
	typedef typename Value<TFile>::Type TValue;
	typedef typename Value<TNames>::Type TName;
	TVertexDescriptor nilVertex = getNil<TVertexDescriptor>();


	if (length(names) < 3) {
		TVertexDescriptor v1 = addVertex(guideTree);
		TVertexDescriptor v2 = addVertex(guideTree);
		TVertexDescriptor internalVertex = addVertex(guideTree);
		addEdge(guideTree, internalVertex, v1, 1.0);
		addEdge(guideTree, internalVertex, v2, 1.0);
		assignRoot(guideTree, internalVertex);
		return;
	}


	typedef std::map<TName, TId> TNameToId;
	TNameToId nameToId;
	for(TId i=0; i<length(names);++i) {
		addVertex(guideTree);	// Create the sequence vertices
		nameToId.insert(std::make_pair(names[i], i));
	}

	TValue c;
	if (_streamEOF(file)) return;
	else c = _streamGet(file);

	TVertexDescriptor lastVertex = nilVertex;
	TVertexDescriptor lastChild = nilVertex;
	while (!_streamEOF(file)) {
		if (c=='(') {
			if (lastVertex == nilVertex) {
				lastVertex = addVertex(guideTree);
				assignRoot(guideTree, lastVertex);
			} else {
				TVertexDescriptor ch = addChild(guideTree, lastVertex);
				lastVertex = ch;
			}
			c = _streamGet(file);
			_parse_skipWhitespace(file, c);
		} else if (c==')') {
			if (!isRoot(guideTree, lastVertex)) {
				lastChild = lastVertex;
				lastVertex = parentVertex(guideTree, lastVertex);
			} else {
				lastChild = lastVertex;
				lastVertex = nilVertex;
			}
			c = _streamGet(file);
			_parse_skipWhitespace(file, c);
		} else if (c==',') {
			c = _streamGet(file);
			_parse_skipWhitespace(file, c);
		} else if (c==':') {
			c = _streamGet(file);
			cargo(findEdge(guideTree, lastVertex, lastChild)) = _parse_readDouble(file,c);
		} else if (c==';') {
			c = _streamGet(file);
			_parse_skipWhitespace(file, c);
		} else {
			TName tmp = _parse_readIdentifier(file, c);
			//std::cout << tmp << std::endl;
			if (lastVertex == nilVertex) {
				// Tree is rooted at a leaf
				// Create artificial root node
				lastVertex = length(names);
				assignRoot(guideTree, addVertex(guideTree));
				addEdge(guideTree, getRoot(guideTree), lastVertex);
				addEdge(guideTree, getRoot(guideTree), nameToId[tmp]);
			} else {
				addEdge(guideTree, lastVertex, nameToId[tmp]);
			}
			lastChild = nameToId[tmp];
		}
	}
	
	// Root the tree if necessary 
	if (outDegree(guideTree, lastChild) > 2) {
		TVertexDescriptor myRoot = addVertex(guideTree);
		assignRoot(guideTree, myRoot);
		typedef typename Iterator<TGuideTree, OutEdgeIterator>::Type TOutEdgeIterator;
		TOutEdgeIterator it(guideTree, lastChild);
		TVertexDescriptor tV = targetVertex(it);
		TCargo c = cargo(*it);
		removeEdge(guideTree, lastChild, tV);
		addEdge(guideTree, myRoot, tV, c);
		addEdge(guideTree, myRoot, lastChild, (TCargo) 0);
	}

	//std::fstream strm1; // Alignment graph as dot
	//strm1.open("D:\\matches\\test\\tree.dot", std::ios_base::out | std::ios_base::trunc);
	//write(strm1,guideTree,DotDrawing());
	//strm1.close();


	//std::cout << guideTree << std::endl;
}


}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
