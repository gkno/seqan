// ==========================================================================
//                              parse_alignment
// ==========================================================================
// Copyright (c) 2011, Birte Kehr
//
// ==========================================================================
// Author: bkehr
// ==========================================================================

#ifndef SANDBOX_BKEHR_APPS_BREAKPOINT_CALCULATOR_PARSE_ALIGNMENT_
#define SANDBOX_BKEHR_APPS_BREAKPOINT_CALCULATOR_PARSE_ALIGNMENT_

#include <seqan/align.h>

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

struct Xmfa_;
typedef Tag<Xmfa_> Xmfa;

struct XmfaSwap_;
typedef Tag<XmfaSwap_> XmfaSwap;

struct Maf_;
typedef Tag<Maf_> Maf;


template<typename TPosition, typename TSize>
struct AlignmentBlockRow
{
	TSize rowNum; // row number
	
	TPosition startPos; // in source seq
	TPosition endPos; // in source seq

	bool orientation;
	// startPos is always smaller or equal endPos also if orientation is false

	AlignmentBlockRow() {}

	AlignmentBlockRow(TSize id, TPosition start, TPosition end, bool ori)
	{
		rowNum = id;
		startPos = start;
		endPos = end;
		orientation = ori;
	}
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

template<typename TStream, typename TPassSpec>
int
skipHeader(RecordReader<TStream, SinglePass<TPassSpec> > & recordReader,
		   Maf const &)
{
	CharString buffer;
	int res = 0;

	//// Check first line
	//res = readUntilWhitespace(buffer, recordReader);
	//if (res) return res;
	//if (buffer != "##maf")
	//	return 1;  // FORMAT ERROR, should probably be a constant
	//res = skipLine(recordReader);
	//if (res) return res;

	// Skip comment lines
	while (value(recordReader) == '#')
	{
		res = skipLine(recordReader);
		if (res) return res;
		if (atEnd(recordReader)) return 1;
	}

	return 0;
}

template<typename TStream, typename TPassSpec, typename TTag>
int
skipHeader(RecordReader<TStream, SinglePass<TPassSpec> > & recordReader,
		   TTag const &)
{
	CharString buffer;
	int res = 0;

	//// Check first line
	//res = readUntilWhitespace(buffer, recordReader);
	//if (res) return res;
	//if (buffer != "#FormatVersion")
	//	return 1;  // FORMAT ERROR, should probably be a constant
	//res = skipLine(recordReader);
	//if (res) return res;

	// Skip comment lines
	while (value(recordReader) == '#')
	{
		res = skipLine(recordReader);
		if (res) return res;
		if (atEnd(recordReader)) return 1;
	}

	return 0;
}

template<typename TAlignRow>
bool
readGappedSeq(TAlignRow & gapseq,
			  CharString & buffer)
{
	typedef typename Source<TAlignRow>::Type TSequence;
	typedef typename Value<TSequence>::Type TAlphabet;

	typedef typename Position<TSequence>::Type TPosition;
	typedef typename Size<TSequence>::Type TSize;
	typedef Pair<TPosition, TSize> TGapPosAndLen; // source position and length defining a gap
	typedef String<TGapPosAndLen> TGapsString;

	if (length(buffer) == 0) return 0;

	TSequence seq;
	TGapsString gaps;

	// read sequence, gap positions and gap lengths
	TSize gapLen = 0;
	TPosition sourcePos = 0;

	typename Iterator<CharString>::Type it = begin(buffer);
	typename Iterator<CharString>::Type itEnd = end(buffer);
	for (; it != itEnd; ++it)
	{
		if (*it == '-')
		{
			++gapLen;
		}
		else
		{
			if (gapLen > 0)
			{
				appendValue(gaps, TGapPosAndLen(sourcePos, gapLen));
				gapLen = 0;
			}
			appendValue(seq, static_cast<TAlphabet>(*it));
			++sourcePos;
		}
	}

	// set sequence
	assignSource(gapseq, seq);

	// insert gaps from back to front
	if (length(gaps) > 0)
	{
		typename Iterator<TGapsString>::Type itGaps = end(gaps);
		typename Iterator<TGapsString>::Type itGapsBegin = begin(gaps);
		for (--itGaps; itGaps >= itGapsBegin; --itGaps)
			insertGaps(gapseq, (*itGaps).i1, (*itGaps).i2);
	}

	return 0;
}

// Reads one collinear alignment block from a file in XMFA format
template<typename TSequence, typename TSize, typename TStream, typename TPassSpec>
int
readRecord(Align<TSequence, ArrayGaps> & align,
		   std::map<CharString, AlignmentBlockRow<TSize, TSize> > & idToRowMap,
		   RecordReader<TStream, SinglePass<TPassSpec> > & recordReader,
		   bool swapPos,
		   Xmfa const &)
{
	typedef typename Position<TSequence>::Type TPosition;

	CharString buffer;
	int res = 0;

	TPosition startPos, endPos;
	CharString id;
	bool orientation;

	while (!skipChar(recordReader, '>'))
	{
		typename Row<Align<TSequence, ArrayGaps> >::Type gapseq;

		res = skipWhitespaces(recordReader);
		if (res) return res;

		// read seq id
		clear(buffer);
		res = readDigits(buffer, recordReader);
		if (res) return res;
		id = buffer;

		res = skipChar(recordReader, ':');
		if (res) return res;

		// read start position
		clear(buffer);
		res = readDigits(buffer, recordReader);
		if (res) return res;
		startPos = lexicalCast<TPosition>(buffer);

		res = skipChar(recordReader, '-');
		if (res) return res;

		// read end position
		clear(buffer);
		res = readDigits(buffer, recordReader);
		if (res) return res;
		endPos = lexicalCast<TPosition>(buffer);

		res = skipWhitespaces(recordReader);
		if (res) return res;

		// read orientation
		clear(buffer);
		res = readNChars(buffer, recordReader, 1);
		if (res) return res;

		if (buffer  == "+") orientation = true;
		else if (buffer == "-")	orientation = false;
		else return 1;

		if (swapPos && orientation == false)
		{
			TPosition help = startPos;
			startPos = endPos;
			endPos = help;
		}
		if (startPos != 0) --startPos;

		// skip rest of line
		res = skipLine(recordReader);
		if (res) return res;

		clear(buffer);
		while (value(recordReader) != '>' && value(recordReader) != '=')
		{
			res = readLine(buffer, recordReader);
			if (res) return res;
		}
		res = readGappedSeq(gapseq, buffer);
		if (res) return res;

		//std::cout << endPos << "  " << startPos << std::endl;
		SEQAN_ASSERT_EQ(length(source(gapseq)), endPos-startPos);
		
		if (endPos != startPos) {
			idToRowMap[id] = AlignmentBlockRow<TSize, TSize>(length(rows(align)), startPos, endPos, orientation);
			appendValue(rows(align), gapseq);
		}
	}

	if (value(recordReader) != '=')
		return 1;

	skipUntilLineBeginsWithChar(recordReader, '>');
	
	return 0;
}

template<typename TSequence, typename TSize, typename TStream, typename TPassSpec>
int
readRecord(Align<TSequence, ArrayGaps> & align,
		   std::map<CharString, AlignmentBlockRow<TSize, TSize> > & idToRowMap,
		   RecordReader<TStream, SinglePass<TPassSpec> > & recordReader,
		   Xmfa const &)
{
	return readRecord(align, idToRowMap, recordReader, false, Xmfa());
}

template<typename TSequence, typename TSize, typename TStream, typename TPassSpec>
int
readRecord(Align<TSequence, ArrayGaps> & align,
		   std::map<CharString, AlignmentBlockRow<TSize, TSize> > & idToRowMap,
		   RecordReader<TStream, SinglePass<TPassSpec> > & recordReader,
		   XmfaSwap const &)
{
	return readRecord(align, idToRowMap, recordReader, true, Xmfa());
}

// Reads one collinear alignment block from a file in MAF format
template<typename TSequence, typename TSize, typename TStream, typename TPassSpec>
int
readRecord(Align<TSequence, ArrayGaps> & align,
		   std::map<CharString, AlignmentBlockRow<TSize, TSize> > & idToRowMap,
		   RecordReader<TStream, SinglePass<TPassSpec> > & recordReader,
		   Maf const &)
{
	typedef typename Position<TSequence>::Type TPosition;

	CharString buffer;
	int res = 0;

	TPosition startPos, len, seqLen;
	CharString id;
	bool orientation;

	skipWhitespaces(recordReader);
	res = skipChar(recordReader, 'a');
	if (res) return res;

	res = skipLine(recordReader);
	if (res) return res;

	skipWhitespaces(recordReader);
	while (!skipChar(recordReader, 's'))
	{
		typename Row<Align<TSequence, ArrayGaps> >::Type gapseq;

		res = skipWhitespaces(recordReader);
		if (res) return res;

		// read seq id
		clear(buffer);
		res = readUntilWhitespace(buffer, recordReader);
		if (res) return res;
		id = buffer;

		res = skipWhitespaces(recordReader);
		if (res) return res;

		// read start position
		clear(buffer);
		res = readDigits(buffer, recordReader);
		if (res) return res;
		startPos = lexicalCast<TPosition>(buffer);
		
		res = skipWhitespaces(recordReader);
		if (res) return res;

		// read length
		clear(buffer);
		res = readDigits(buffer, recordReader);
		if (res) return res;
		len = lexicalCast<TPosition>(buffer);
		
		res = skipWhitespaces(recordReader);
		if (res) return res;

		// skip orientation
		clear(buffer);
		res = readNChars(buffer, recordReader, 1);
		if (res) return res;

		if (buffer  == '+')	orientation = true;
		else if (buffer == '-') orientation = false;
		else return 1;

		res = skipWhitespaces(recordReader);
		if (res) return res;

		// skip seq length
		clear(buffer);
		res = readDigits(buffer, recordReader);
		if (res) return res;
		seqLen = lexicalCast<TPosition>(buffer);

		res = skipWhitespaces(recordReader);
		if (res) return res;

		// read gapped seq
		clear(buffer);
		res = readUntilWhitespace(buffer, recordReader);
		if (res) return res;

		res = readGappedSeq(gapseq, buffer);
		if (res) return res;

		SEQAN_ASSERT_EQ(length(source(gapseq)), len);

		if (!orientation)
			startPos = seqLen - (startPos+len);

		if (len > 0) {
			if (idToRowMap.count(id) != 0)
			{
				std::cerr << "ERROR: Sequence " << id << " occurs twice in a block!\n";
				return 1;
			}
			idToRowMap[id] = AlignmentBlockRow<TSize, TSize>(length(rows(align)), startPos, startPos+len, orientation);
			appendValue(rows(align), gapseq);
		}

		skipWhitespaces(recordReader);
	}
	return 0;
}


template<typename TDepStringSet, typename TStringSet, typename TSize, typename TSequence>
int
buildAlignmentGraph(Graph<Alignment<TDepStringSet> > & graph,
					TStringSet & seqs,
					String<std::map<CharString, AlignmentBlockRow<TSize, TSize> > > const & idToRowMaps,
					String<Align<TSequence, ArrayGaps> > const & aligns)
{
	typedef typename Id<TDepStringSet>::Type TId;
	typedef Graph<Alignment<TDepStringSet> > TAlignmentGraph;
	typedef typename VertexDescriptor<TAlignmentGraph>::Type TVertexDescriptor;

	typedef std::map<CharString, Triple<TSize> > TMap;
	typedef typename TMap::const_iterator TMapIterator;

	// count the number of sequences
	std::map<TId, std::set<TSize> > startPositionsById; 
	std::map<TId, std::set<TSize> > endPositionsById;
	std::map<CharString, TId> idToSeqId;
	typename Iterator<String<TMap> >::Type it = begin(idToRowMaps);
	for (; it != end(idToRowMaps); ++it)
	{
		for (TMapIterator mapIt = (*it).begin(); mapIt != (*it).end(); ++mapIt) 
		{
			if (idToSeqId.count((*mapIt).first) == 0)
				idToSeqId[(*mapIt).first] = idToSeqId.size();
			TId seqId = idToSeqId[(*mapIt).first];
			startPositionsById[seqId].insert((*mapIt).second.startPos);
			endPositionsById[seqId].insert((*mapIt).second.endPos);
		}
	}

	// set sequences
	TDepStringSet depSeqs;
	typename std::map<TId, std::set<TSize> >::iterator eIt = endPositionsById.begin();
	for (; eIt != endPositionsById.end(); ++eIt)
	{
		TSequence seq;
		resize(seq, *(--(*eIt).second.end()));
		appendValue(seqs, seq);
		appendValue(depSeqs, seqs[length(seqs)-1]);
	}

	graph = Graph<Alignment<TDepStringSet> >(depSeqs);

	// add vertices
	typename std::map<TId, std::set<TSize> >::iterator sIt = startPositionsById.begin();
	for (eIt = endPositionsById.begin(); eIt != endPositionsById.end(); ++sIt, ++eIt)
	{
		TId seqId = (*sIt).first;
		typename std::set<TSize>::iterator startPosIt = (*sIt).second.begin();
		typename std::set<TSize>::iterator endPosIt = (*eIt).second.begin();
		for (; startPosIt != (*sIt).second.end(); ++startPosIt, ++endPosIt)
		{
			if ((*endPosIt) != (*startPosIt))
				addVertex(graph, seqId, *startPosIt, (*endPosIt) - (*startPosIt));
		}
	}

	// add edges and set seqs
	int i = 0;
	for (it = begin(idToRowMaps); it < end(idToRowMaps); ++it, ++i)
	{
		for (TMapIterator mapIt = (*it).begin(); mapIt != (*it).end(); ++mapIt) 
		{
			TVertexDescriptor v1 = findVertex(graph, idToSeqId[(*mapIt).first], (*mapIt).second.startPos);
			if (v1 == getNil<TVertexDescriptor>()) return 1;
			TMapIterator mapIt2 = mapIt;
			for (++mapIt2; mapIt2 != (*it).end(); ++mapIt2)
			{
				TId seqId = idToSeqId[(*mapIt2).first];
				infix(valueById(seqs, seqId), (*mapIt2).second.startPos (*mapIt2).second.endPos) = source(row(aligns[i], (*mapIt2).second.rowNum));
				TVertexDescriptor v2 = findVertex(graph, seqId, (*mapIt2).second.startPos);
				//int score = ;
				addEdge(graph, v1, v2/*, score*/);
				if (v1 == getNil<TVertexDescriptor>()) return 1;
			}
		}
	}

	return 0;
}

template<typename TSequence, typename TStringSet, typename TFile, typename TTag>
int
parseAlignment(String<Align<TSequence, ArrayGaps> > & aligns,
			   String<std::map<CharString, AlignmentBlockRow<typename Size<TStringSet>::Type, typename Size<TStringSet>::Type> > > & idToRowMaps,
			   TStringSet & /*seqs*/,
			   TFile & file,
			   bool verbose,
			   TTag const tag)
{
	typedef typename Size<TStringSet>::Type TSize;
	typedef std::map<CharString, AlignmentBlockRow<TSize, TSize> > TMap;

	typedef Align<TSequence, ArrayGaps> TAlign;

	RecordReader<TFile, SinglePass<> > recordReader(file);
	int res = 0;

	res = skipHeader(recordReader, tag);
	if (res) return res;

	while (!atEnd(recordReader))
	{
		TAlign align;
		TMap idToRow;
		res = readRecord(align, idToRow, recordReader, tag);
		if (res) return res;

		appendValue(aligns, align);
		appendValue(idToRowMaps, idToRow);

		skipWhitespaces(recordReader);
		while (!atEnd(recordReader) && value(recordReader) == '#')
			res = skipLine(recordReader);
	}
	if (!atEnd(recordReader)) return 1;
	
	if (verbose)
	{
		std::cout << length(aligns) << " alignment block";
		if (length(aligns) != 1) std::cout << "s";
		std::cout << " loaded." << std::endl;
	}

	return 0;
}

template<typename TDepStringSet, typename TStringSet, typename TFile, typename TTag>
int
parseAlignment(Graph<Alignment<TDepStringSet> > & graph,
			   TStringSet & seqs,
			   TFile & file,
			   bool verbose,
			   TTag tag)
{	
	typedef typename Size<TStringSet>::Type TSize;
	typedef std::map<CharString, AlignmentBlockRow<TSize, TSize> > TMap;
	String<TMap> idToRowMaps;

	typedef typename Value<TStringSet>::Type TSequence;
	typedef Align<TSequence, ArrayGaps> TAlign;
	String<TAlign> aligns;

	if (parseAlignment(aligns, idToRowMaps, seqs, file, false, tag))
	{
		std::cerr << "ERROR: Failed reading alignment blocks from input file." << std::endl;
		return 1;
	}

	if (buildAlignmentGraph(graph, seqs, idToRowMaps, aligns))
	{
		std::cerr << "ERROR: Failed building alignment graph from alignment blocks." << std::endl;
		return 1;
	}

	if (verbose)
	{
		std::cout << "Finished building alignment graph from alignment blocks." << std::endl;
	}

	return 0;
}

template<typename TAlignmentGraph, typename TStringSet, typename TFile, typename TTag>
int
parseAlignment(TAlignmentGraph & align, TStringSet & seqs, TFile & file, TTag tag)
{
	return parseAlignment(align, seqs, file, false, tag);
}


} // namespace seqan

#endif  // #ifndef SANDBOX_BKEHR_APPS_BREAKPOINT_CALCULATOR_PARSE_ALIGNMENT_
