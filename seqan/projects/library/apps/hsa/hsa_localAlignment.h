
#ifndef SEQAN_HEADER_HSA_LOCALALIGN_H
#define SEQAN_HEADER_HSA_LOCALALIGN_H

#include <apps/stellar/stellar.h>

using namespace seqan;


template<typename TAlignmentGraph, typename TPosition, typename TSize>
void
_findBeginPositions(TAlignmentGraph & parentGraph, 
					String<TPosition> const & segmentBegins,
					TSize i,
					String<TPosition> & beginPos) {
	typedef typename VertexDescriptor<TAlignmentGraph>::Type TVertexDescriptor;
	typedef typename Iterator<TAlignmentGraph, OutEdgeIterator>::Type TOutEdgeIterator;

	if (beginPos[i] == segmentBegins[i]) return;

	TVertexDescriptor v = findVertex(parentGraph, i, beginPos[i]-1);
	SEQAN_ASSERT_GT(outDegree(parentGraph, v), 0u);
	TOutEdgeIterator outEdgeIt(parentGraph, v);

	while (!atEnd(outEdgeIt)) {
		TVertexDescriptor w = targetVertex(outEdgeIt);
		TSize j = sequenceId(parentGraph, w);
		if (fragmentBegin(parentGraph, w) + fragmentLength(parentGraph, w) > value(beginPos, j)) {
			value(beginPos, j) = fragmentBegin(parentGraph, w) + fragmentLength(parentGraph, w);
			_findBeginPositions(parentGraph, segmentBegins, j, beginPos);
		}
		++outEdgeIt;
	}
}

template<typename TAlignmentGraph, typename TPosition, typename TSize>
void
_findEndPositions(TAlignmentGraph & parentGraph,
				  String<TPosition> const & segmentEnds,
				  TSize i,
				  String<TPosition> & endPos) {
	typedef typename VertexDescriptor<TAlignmentGraph>::Type TVertexDescriptor;
	typedef typename Iterator<TAlignmentGraph, OutEdgeIterator>::Type TOutEdgeIterator;

	if (endPos[i] == segmentEnds[i]) return;

	TVertexDescriptor v = findVertex(parentGraph, i, endPos[i]);
	SEQAN_ASSERT_GT(outDegree(parentGraph, v), 0u);
	TOutEdgeIterator outEdgeIt(parentGraph, v);

	while (!atEnd(outEdgeIt)) {
		TVertexDescriptor w = targetVertex(outEdgeIt);
		TSize j = sequenceId(parentGraph, w);
		if (fragmentBegin(parentGraph, w) < value(endPos, j)) {
			value(endPos, j) = fragmentBegin(parentGraph, w);
			_findEndPositions(parentGraph, segmentEnds, j, endPos);
		}
		++outEdgeIt;
	}
}

template<typename TSequence, typename TSegmentSet>
inline void
truncateMatchBegin(StellarMatch<TSequence, CharString> & match, TSegmentSet const & segments) {
	typedef StellarMatch<TSequence, CharString> TMatch;
	typedef typename TMatch::TPos TPos;

	if (match.begin1 < beginPosition(segments[0])) {
		match.begin1 = beginPosition(segments[0]);
		match.begin2 = toSourcePosition(match.row2, toViewPosition(match.row1, match.begin1));
		setClippedBeginPosition(match.row1, match.begin1);
		setClippedBeginPosition(match.row2, match.begin2);
	}
	if (match.begin2 < beginPosition(segments[1])) {
		match.begin2 = beginPosition(segments[1]);
		match.begin1 = toSourcePosition(match.row1, toViewPosition(match.row2, match.begin2));
		setClippedBeginPosition(match.row1, match.begin1);
		setClippedBeginPosition(match.row2, match.begin2);
	}
}

template<typename TSequence, typename TSegmentSet>
inline void
truncateMatchEnd(StellarMatch<TSequence, CharString> & match, TSegmentSet const & segments) {
	typedef StellarMatch<TSequence, CharString> TMatch;
	typedef typename TMatch::TPos TPos;

	if (match.end1 > endPosition(segments[0])) {
		match.end1 = endPosition(segments[0]);
		match.end2 = toSourcePosition(match.row2, toViewPosition(match.row1, match.end1));
		setClippedEndPosition(match.row1, match.end1);
		setClippedEndPosition(match.row2, match.end2);
	}
	if (match.end2 > endPosition(segments[1])) {
		match.end2 = endPosition(segments[1]);
		match.end1 = toSourcePosition(match.row1, toViewPosition(match.row2, match.end2));
		setClippedEndPosition(match.row1, match.end1);
		setClippedEndPosition(match.row2, match.end2);
	}
}

template<typename TSegmentSet, typename TMatch>
inline void
generateLocalMatches(TSegmentSet & segments, StellarParams params, String<TMatch> & matches) {
	typedef typename Value<TSegmentSet>::Type TInfix;
    typedef Finder<TInfix, Swift<SwiftLocal> > TFinder;
	typedef Index<StringSet<TInfix, Dependent<> >, IndexQGram<SimpleShape, OpenAddressing> > TQGramIndex;

	typedef StellarMatch<typename Host<TInfix>::Type, CharString> TEpsMatch;
	typedef typename Iterator<String<TEpsMatch> >::Type TIterator;

	// swift finder
	TFinder swiftFinder(value(segments, 0), 1000, 1);

	// swift pattern
	StringSet<TInfix, Dependent<> > stringSet;
	appendValue(stringSet, value(segments, 1));
    TQGramIndex qgramIndex(stringSet);
    resize(indexShape(qgramIndex), params.qgram);
	Pattern<TQGramIndex, Swift<SwiftLocal> > swiftPattern(qgramIndex);

	// container for eps-matches
	StringSet<QueryMatches<TEpsMatch> > stellarMatches;
    resize(stellarMatches, 1);

	stellar(swiftFinder, swiftPattern, params.epsilon, params.minLength, params.xdrop, stellarMatches, AllLocal());

	// append stellarMatches to matches
	TIterator it = begin(stellarMatches[0].matches);
	TIterator itEnd = end(stellarMatches[0].matches);
	while (it != itEnd) {
		truncateMatchBegin(*it, segments);
		truncateMatchEnd(*it, segments);

		TMatch m;
		resize(rows(m), 2);
		row(m, 0) = (*it).row1;
		row(m, 1) = (*it).row2;
		//std::cout << (*it).begin1 << ".." << (*it).end1 << " , " << (*it).begin2 << ".." << (*it).end2 << std::endl;
		//std::cout << m;

		appendValue(matches, m);
		++it;
	}
	//std::cout << "# STELLAR matches: " << length(matches) << "   n0 = " << params.minLength << "  eps = " << params.epsilon << std::endl;
}

template<typename TSegmentSet, typename TOptions, typename TLevel, typename TMatch>
inline void
findStellarMatches(TSegmentSet & segs,
                   TOptions & options,
				   TLevel recursionLevel,
                   String<TMatch> & matches) {
	typedef typename Size<TSegmentSet>::Type TSize;

	StellarParams params(options.initialEpsilon + options.deltaEpsilon * recursionLevel,
						 options.initialMinLength - options.deltaMinLength * recursionLevel);

	// stellar on each pair of sequences
    for (TSize i = 0; i < length(segs)-1; ++i) {
        for (TSize j = i+1; j < length(segs); ++j) {
			TSegmentSet segmentPair;
			resize(segmentPair, 2);
			value(segmentPair, 0) = segs[i];
			value(segmentPair, 1) = segs[j];
			generateLocalMatches(segmentPair, params, matches);
        }
    }
}

template<typename TSegmentSet, typename TAlignmentGraph, typename TOptions, typename TLevel, typename TMatch>
inline void
findStellarMatches(TSegmentSet const & segs,
				   TAlignmentGraph & parentGraph,
                   TOptions const & options,
				   TLevel recursionLevel,
                   String<TMatch> & matches) {
	typedef typename Size<TSegmentSet>::Type TSize;
	typedef typename Position<typename Value<TSegmentSet>::Type>::Type TPosition;
	typedef typename VertexDescriptor<TAlignmentGraph>::Type TVertexDescriptor;

	StellarParams params(options.initialEpsilon + options.deltaEpsilon * recursionLevel,
						 options.initialMinLength - options.deltaMinLength * recursionLevel);

	// begin and end position for all sequence segments
	String<TPosition> segmentBegins, segmentEnds;
	resize(segmentBegins, length(segs));
	resize(segmentEnds, length(segs));
	for (TSize j = 0; j < length(segs); ++j) {
		value(segmentBegins, j) = beginPosition(value(segs, j));
		value(segmentEnds, j) = endPosition(value(segs, j));
	}

	// for each sequence segment
	for (TSize i = 0; i < length(segs); ++i) {
		TPosition pos = beginPosition(segs[i]);
		// for each unaligned vertex on sequence segment
		while (pos < endPosition(segs[i])) {
			TVertexDescriptor v = findVertex(parentGraph, i, pos);
			pos += fragmentLength(parentGraph, v);
			if (outDegree(parentGraph, v) != 0) continue;

			// initialize begin and end position for all sequence segments
			String<TPosition> beginPos(segmentBegins);
			String<TPosition> endPos(segmentEnds);

			// set begin and end position for current sequence to vertex begin and end
			value(beginPos, i) = fragmentBegin(parentGraph, v);
			value(endPos, i) = pos; // same as: fragmentBegin(parentGraph, v) + fragmentLength(parentGraph, v);

			// walk through alignment graph to find other sequence segments for pairwise comparison to this segment
			_findBeginPositions(parentGraph, segmentBegins, i, beginPos);
			_findEndPositions(parentGraph, segmentEnds, i, endPos);

			// stellar on each pair of sequence segments
			TSegmentSet segmentPair;
			appendValue(segmentPair, infix(host(segs[i]), beginPos[i], endPos[i]));
			for (TSize j = i+1; j < length(segs); ++j) {
				appendValue(segmentPair, infix(host(segs[j]), beginPos[j], endPos[j]));
				generateLocalMatches(segmentPair, params, matches);
				resize(segmentPair, 1);
			}
		}
	}
}

#endif
