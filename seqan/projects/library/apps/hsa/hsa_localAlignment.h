
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


template<typename TAlignmentGraph, typename TId, typename TPosition>
void
_findBeginPositions(TAlignmentGraph & parentGraph, 
					std::map<TId, TPosition> & segmentBegins,
					typename Id<TAlignmentGraph>::Type i,
					std::map<TId, TPosition> & beginPos) {
	typedef typename Id<TAlignmentGraph>::Type TGraphId;
	typedef typename VertexDescriptor<TAlignmentGraph>::Type TVertexDescriptor;
	typedef typename Iterator<TAlignmentGraph, OutEdgeIterator>::Type TOutEdgeIterator;
	
	TId id_i = id(value(stringSet(parentGraph), idToPosition(stringSet(parentGraph), i)));
	if (beginPos[id_i] == segmentBegins[id_i]) return;

	// iterate over all vertices before beginPos[id_i]
	TPosition pos = beginPos[id_i];
	while (pos > segmentBegins[id_i]) {

		TVertexDescriptor v = findVertex(parentGraph, i, pos-1);
		TOutEdgeIterator outEdgeIt(parentGraph, v);

		// iterate over all neighbors of v and recurse if updated beginPos
		while (!atEnd(outEdgeIt)) {
			TVertexDescriptor w = targetVertex(outEdgeIt);
			TGraphId j = sequenceId(parentGraph, w);
			TId id_j = id(value(stringSet(parentGraph), idToPosition(stringSet(parentGraph), j)));
			if (fragmentBegin(parentGraph, w) + fragmentLength(parentGraph, w) > beginPos[id_j]) {
				beginPos[id_j] = fragmentBegin(parentGraph, w) + fragmentLength(parentGraph, w);
				_findBeginPositions(parentGraph, segmentBegins, j, beginPos);
			}
			++outEdgeIt;
		}
		pos -= fragmentLength(parentGraph, v);
	}
}

template<typename TAlignmentGraph, typename TId, typename TPosition>
void
_findEndPositions(TAlignmentGraph & parentGraph,
				  std::map<TId, TPosition> & segmentEnds,
				  typename Id<TAlignmentGraph>::Type i,
				  std::map<TId, TPosition> & endPos) {
	typedef typename Id<TAlignmentGraph>::Type TGraphId;
	typedef typename VertexDescriptor<TAlignmentGraph>::Type TVertexDescriptor;
	typedef typename Iterator<TAlignmentGraph, OutEdgeIterator>::Type TOutEdgeIterator;

	TId id_i = id(value(stringSet(parentGraph), idToPosition(stringSet(parentGraph), i)));
	if (endPos[id_i] == segmentEnds[id_i]) return;

	// iterate over all vertices behind endPos[id_i]
	TPosition pos = endPos[id_i];
	while (pos < segmentEnds[id_i]) {

		TVertexDescriptor v = findVertex(parentGraph, i, pos);
		TOutEdgeIterator outEdgeIt(parentGraph, v);

		// iterate over all neighbors of v and recurse if updated endPos
		while (!atEnd(outEdgeIt)) {
			TVertexDescriptor w = targetVertex(outEdgeIt);
			TGraphId j = sequenceId(parentGraph, w);
			TId id_j = id(value(stringSet(parentGraph), idToPosition(stringSet(parentGraph), j)));
			if (fragmentBegin(parentGraph, w) < endPos[id_j]) {
				endPos[id_j] = fragmentBegin(parentGraph, w);
				_findEndPositions(parentGraph, segmentEnds, j, endPos);
			}
			++outEdgeIt;
		}
		pos += fragmentLength(parentGraph, v);
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

	// TODO implement an alternative version on connected components?

	// begin and end position for all sequence segments
	std::map<TId, TPosition> segmentBegins, segmentEnds;
	TMapIterator mapEnd = segs.end();
	for (TMapIterator mapIt = segs.begin(); mapIt != mapEnd; ++mapIt) {
		TId id = mapIt->first;
		TInfix segment = mapIt->second;
		segmentBegins[id] = beginPosition(segment);
		segmentEnds[id] = endPosition(segment);
	}

	// for each sequence segment
	for (TMapIterator mapIt = segs.begin(); mapIt != mapEnd; ++mapIt) {
		TId id = mapIt->first;
		TPosition pos = beginPosition(mapIt->second);

		std::map<TId, Pair<TPosition> > previousBeginPos;
		std::map<TId, Pair<TPosition> > previousEndPos;
		TMapIterator mapIt2 = mapIt;
		for (mapIt2++; mapIt2 != mapEnd; ++mapIt2) {
			TId id2 = mapIt2->first;
 			previousBeginPos[id2] = Pair<TPosition>(beginPosition(mapIt2->second), beginPosition(mapIt->second));
			previousEndPos[id2] = Pair<TPosition>(beginPosition(mapIt2->second), beginPosition(mapIt->second));
		}
		// for each unaligned vertex on sequence segment
		while (pos < endPosition(mapIt->second)) {
			typename Id<TAlignmentGraph>::Type graphId = stringIdToStringSetId(stringSet(parentGraph), id);
			TVertexDescriptor v = findVertex(parentGraph, graphId, pos);
			pos += fragmentLength(parentGraph, v);
			if (outDegree(parentGraph, v) != 0) continue;

			// initialize begin and end position for all sequence segments
			std::map<TId, TPosition> beginPos(segmentBegins);
			std::map<TId, TPosition> endPos(segmentEnds);

			// set begin and end position for current sequence to vertex begin and end
			beginPos[id] = fragmentBegin(parentGraph, v);
			endPos[id] = pos;
			SEQAN_ASSERT_EQ(pos, fragmentBegin(parentGraph, v) + fragmentLength(parentGraph, v));

			// walk through alignment graph to find other sequence segments for pairwise comparison to this segment
			_findBeginPositions(parentGraph, segmentBegins, graphId, beginPos);
			_findEndPositions(parentGraph, segmentEnds, graphId, endPos);

			// stellar on each pair of sequence segments
			mapIt2 = mapIt;
			for (mapIt2++; mapIt2 != mapEnd; ++mapIt2) {
				TId id2 = mapIt2->first;
				// generate "closed" matches from previous segment
				if (previousBeginPos[id2].i1 != beginPos[id2]) {
					Pair<TInfix> segmentPair = Pair<TInfix>();
					segmentPair.i1 = infix(host(mapIt->second), previousBeginPos[id2].i2, previousEndPos[id2].i2);
					segmentPair.i2 = infix(host(mapIt2->second), previousBeginPos[id2].i1, previousEndPos[id2].i1);
					generateLocalMatches(segmentPair, params, matches);
					previousBeginPos[id2] = Pair<TPosition>(beginPos[id2], beginPos[id]);
					previousEndPos[id2] = Pair<TPosition>(endPos[id2], endPos[id]);
				} else {
					previousBeginPos[id2] = Pair<TPosition>(beginPos[id2], previousBeginPos[id2].i2);
					previousEndPos[id2] = Pair<TPosition>(endPos[id2], endPos[id]);
				}
			}
		}
		// generate remaining matches
		mapIt2 = mapIt;
		for (mapIt2++; mapIt2 != mapEnd; ++mapIt2) {
			TId id2 = mapIt2->first;
			Pair<TInfix> segmentPair = Pair<TInfix>();
			segmentPair.i1 = infix(host(mapIt->second), previousBeginPos[id2].i2, previousEndPos[id2].i2);
			segmentPair.i2 = infix(host(mapIt2->second), previousBeginPos[id2].i1, previousEndPos[id2].i1);
			generateLocalMatches(segmentPair, params, matches);
		}
	}
}

#endif
