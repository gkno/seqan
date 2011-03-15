
#ifndef SEQAN_HEADER_HSA_LOCALALIGN_H
#define SEQAN_HEADER_HSA_LOCALALIGN_H

#include <apps/stellar/stellar.h>

using namespace seqan;

template<typename TStringSet, typename TId>
typename Position<TStringSet>::Type
stringIdToPosition(TStringSet & stringSet, TId & stringId) {
	typedef typename Iterator<TStringSet>::Type TIterator;
	typedef typename Position<TStringSet>::Type TPosition;

	TIterator itEnd = end(stringSet);
	for (TIterator it = begin(stringSet); it != itEnd; ++it) {
		if (id(*it) == stringId) {
			return (TPosition)(it - begin(stringSet));
		}
	}

	return MaxValue<TPosition>::VALUE;
}

template<typename TStringSet, typename TId>
typename Id<TStringSet>::Type
stringIdToStringSetId(TStringSet & stringSet, TId & stringId) {
	typedef typename Iterator<TStringSet>::Type TIterator;

	TIterator itEnd = end(stringSet);
	for (TIterator it = begin(stringSet); it != itEnd; ++it) {
		if (id(*it) == stringId) {
			return positionToId(stringSet, it - begin(stringSet));
		}
	}

	return MaxValue<typename Id<TStringSet>::Type >::VALUE;
}


template<typename TAlignmentGraph, typename TId, typename TPosition, typename TSize>
void
_findBeginPositions(TAlignmentGraph & parentGraph, 
					std::map<TId, TPosition> & segmentBegins,
					TSize i,
					std::map<TId, TPosition> & beginPos) {
	typedef typename VertexDescriptor<TAlignmentGraph>::Type TVertexDescriptor;
	typedef typename Iterator<TAlignmentGraph, OutEdgeIterator>::Type TOutEdgeIterator;
	
	TId id_i = id(value(stringSet(parentGraph), i));
	if (beginPos[id_i] == segmentBegins[id_i]) return;

	// TODO: Schleife über alle davor liegenden Segmente???

	TVertexDescriptor v = findVertex(parentGraph, i, beginPos[id_i]-1);
	SEQAN_ASSERT_GT(outDegree(parentGraph, v), 0u);
	TOutEdgeIterator outEdgeIt(parentGraph, v);

	while (!atEnd(outEdgeIt)) {
		TVertexDescriptor w = targetVertex(outEdgeIt);
		TSize j = sequenceId(parentGraph, w);
		TId id_j = id(value(stringSet(parentGraph), j));
		if (fragmentBegin(parentGraph, w) + fragmentLength(parentGraph, w) > beginPos[id_j]) {
			beginPos[id_j] = fragmentBegin(parentGraph, w) + fragmentLength(parentGraph, w);
			_findBeginPositions(parentGraph, segmentBegins, j, beginPos);
		}
		++outEdgeIt;
	}
}

template<typename TAlignmentGraph, typename TId, typename TPosition, typename TSize>
void
_findEndPositions(TAlignmentGraph & parentGraph,
				  std::map<TId, TPosition> & segmentEnds,
				  TSize i,
				  std::map<TId, TPosition> & endPos) {
	typedef typename VertexDescriptor<TAlignmentGraph>::Type TVertexDescriptor;
	typedef typename Iterator<TAlignmentGraph, OutEdgeIterator>::Type TOutEdgeIterator;

	TId id_i = id(value(stringSet(parentGraph), i));
	if (endPos[id_i] == segmentEnds[id_i]) return;

	// TODO: Schleife über alle dahinter liegenden Segmente???

	TVertexDescriptor v = findVertex(parentGraph, i, endPos[id_i]);
	SEQAN_ASSERT_GT(outDegree(parentGraph, v), 0u);
	TOutEdgeIterator outEdgeIt(parentGraph, v);

	while (!atEnd(outEdgeIt)) {
		TVertexDescriptor w = targetVertex(outEdgeIt);
		TSize j = sequenceId(parentGraph, w);
		TId id_j = id(value(stringSet(parentGraph), j));
		if (fragmentBegin(parentGraph, w) < endPos[id_j]) {
			endPos[id_j] = fragmentBegin(parentGraph, w);
			_findEndPositions(parentGraph, segmentEnds, j, endPos);
		}
		++outEdgeIt;
	}
}

template<typename TSequence, typename TInfix>
inline void
truncateMatchBegin(StellarMatch<TSequence, CharString> & match, Pair<TInfix> & segments) {
	typedef StellarMatch<TSequence, CharString> TMatch;
	typedef typename TMatch::TPos TPos;

	if (match.begin1 < beginPosition(segments.i1)) {
		match.begin1 = beginPosition(segments.i1);
		match.begin2 = toSourcePosition(match.row2, toViewPosition(match.row1, match.begin1));
		setClippedBeginPosition(match.row1, match.begin1);
		setClippedBeginPosition(match.row2, match.begin2);
	}
	if (match.begin2 < beginPosition(segments.i2)) {
		match.begin2 = beginPosition(segments.i2);
		match.begin1 = toSourcePosition(match.row1, toViewPosition(match.row2, match.begin2));
		setClippedBeginPosition(match.row1, match.begin1);
		setClippedBeginPosition(match.row2, match.begin2);
	}
}

template<typename TSequence, typename TInfix>
inline void
truncateMatchEnd(StellarMatch<TSequence, CharString> & match, Pair<TInfix> & segments) {
	typedef StellarMatch<TSequence, CharString> TMatch;
	typedef typename TMatch::TPos TPos;

	if (match.end1 > endPosition(segments.i1)) {
		match.end1 = endPosition(segments.i1);
		match.end2 = toSourcePosition(match.row2, toViewPosition(match.row1, match.end1));
		setClippedEndPosition(match.row1, match.end1);
		setClippedEndPosition(match.row2, match.end2);
	}
	if (match.end2 > endPosition(segments.i2)) {
		match.end2 = endPosition(segments.i2);
		match.end1 = toSourcePosition(match.row1, toViewPosition(match.row2, match.end2));
		setClippedEndPosition(match.row1, match.end1);
		setClippedEndPosition(match.row2, match.end2);
	}
}

template<typename TInfix, typename TMatch>
inline void
generateLocalMatches(Pair<TInfix> & segments, StellarParams params, String<TMatch> & matches) {
    typedef Finder<TInfix, Swift<SwiftLocal> > TFinder;
	typedef Index<StringSet<TInfix, Dependent<> >, IndexQGram<SimpleShape, OpenAddressing> > TQGramIndex;

	typedef StellarMatch<typename Host<TInfix>::Type, CharString> TEpsMatch;
	typedef typename Iterator<String<TEpsMatch> >::Type TIterator;

	// swift finder
	TFinder swiftFinder(segments.i1, 1000, 1);

	// swift pattern
	StringSet<TInfix, Dependent<> > stringSet;
	appendValue(stringSet, segments.i2);
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
		//std::cout << id(source((*it).row1)) << ":" << (*it).begin1 << ".." << (*it).end1 << " , ";
		//std::cout << id(source((*it).row2)) << ":" << (*it).begin2 << ".." << (*it).end2 << std::endl;
		//std::cout << m;

		appendValue(matches, m);
		++it;
	}
	//std::cout << "# STELLAR matches: " << length(matches) << "   n0 = " << params.minLength << "  eps = " << params.epsilon << std::endl;
}

template<typename TId, typename TInfix, typename TOptions, typename TLevel, typename TMatch>
inline void
findStellarMatches(std::map<TId, TInfix> & segs,
                   TOptions & options,
				   TLevel recursionLevel,
                   String<TMatch> & matches) {
	typedef typename std::map<TId, TInfix>::size_type TSize;
	typedef typename std::map<TId, TInfix>::const_iterator TMapIterator;

	StellarParams params(options.initialEpsilon + options.deltaEpsilon * recursionLevel,
						 options.initialMinLength - options.deltaMinLength * recursionLevel);

	// stellar on each pair of sequences
	TMapIterator end1 = --segs.end();
    for (TMapIterator it1 = segs.begin(); it1 != end1; ++it1) {
		TMapIterator it2 = it1;
		TMapIterator end2 = segs.end();
        for (++it2; it2 != end2; ++it2) {
			Pair<TInfix> segmentPair = Pair<TInfix>();
			segmentPair.i1 = it1->second;  // TODO: aus Schleife herausziehen
			segmentPair.i2 = it2->second;
			generateLocalMatches(segmentPair, params, matches);
        }
    }
}

template<typename TId, typename TInfix, typename TAlignmentGraph, typename TOptions, typename TLevel, typename TMatch>
inline void
findStellarMatches(std::map<TId, TInfix> & segs,
				   TAlignmentGraph & parentGraph,
                   TOptions const & options,
				   TLevel recursionLevel,
                   String<TMatch> & matches) {
	typedef typename std::map<TId, TInfix>::const_iterator TMapIterator;
	typedef typename Position<TInfix>::Type TPosition;
	typedef typename VertexDescriptor<TAlignmentGraph>::Type TVertexDescriptor;

	StellarParams params(options.initialEpsilon + options.deltaEpsilon * recursionLevel,
						 options.initialMinLength - options.deltaMinLength * recursionLevel);

	// begin and end position for all sequence segments
	std::map<TId, TPosition> segmentBegins, segmentEnds;
	TMapIterator mapIt = segs.begin();
	TMapIterator mapEnd = segs.end();
	while (mapIt != mapEnd) {
		TId id = mapIt->first;
		TInfix segment = mapIt->second;
		segmentBegins[id] = beginPosition(segment);
		segmentEnds[id] = endPosition(segment);
		++mapIt;
	}

	// for each sequence segment
	for (mapIt = segs.begin(); mapIt != mapEnd; ++mapIt) {
		TPosition pos = beginPosition(mapIt->second);
		// for each unaligned vertex on sequence segment
		while (pos < endPosition(mapIt->second)) {
			typename Id<TAlignmentGraph>::Type idPosition = stringIdToPosition(stringSet(parentGraph), mapIt->first);
			TVertexDescriptor v = findVertex(parentGraph, idPosition, pos);
			pos += fragmentLength(parentGraph, v);
			if (outDegree(parentGraph, v) != 0) continue;

			// initialize begin and end position for all sequence segments
			std::map<TId, TPosition> beginPos(segmentBegins);
			std::map<TId, TPosition> endPos(segmentEnds);

			// set begin and end position for current sequence to vertex begin and end
			beginPos[mapIt->first] = fragmentBegin(parentGraph, v);
			endPos[mapIt->first] = pos; // same as: fragmentBegin(parentGraph, v) + fragmentLength(parentGraph, v);

			// walk through alignment graph to find other sequence segments for pairwise comparison to this segment
			_findBeginPositions(parentGraph, segmentBegins, idPosition, beginPos);
			_findEndPositions(parentGraph, segmentEnds, idPosition, endPos);

			// stellar on each pair of sequence segments
			TMapIterator mapIt2 = mapIt;
			for (mapIt2++; mapIt2 != mapEnd; ++mapIt2) {
				Pair<TInfix> segmentPair = Pair<TInfix>();
				segmentPair.i1 = infix(host(mapIt->second), beginPos[mapIt->first], endPos[mapIt->first]); // TODO: aus Schleife herausziehen
				segmentPair.i2 = infix(host(mapIt2->second), beginPos[mapIt2->first], endPos[mapIt2->first]);
				generateLocalMatches(segmentPair, params, matches);
			}
		}
	}
}

#endif
