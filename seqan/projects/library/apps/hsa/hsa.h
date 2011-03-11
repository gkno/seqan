
#ifndef SEQAN_HEADER_HSA_H
#define SEQAN_HEADER_HSA_H

#include "hsa_types.h"
#include "hsa_localAlignment.h"

using namespace seqan;

template<typename TSegmentSet, typename TAlignmentGraph, typename TMatch>
void
_addParentMatches(TSegmentSet & segments, TAlignmentGraph & parentGraph, String<TMatch> & matches) {
	typedef typename Size<TSegmentSet>::Type TSize;
	typedef typename Position<typename Value<TSegmentSet>::Type>::Type TPosition;
	typedef typename Iterator<TAlignmentGraph, OutEdgeIterator>::Type TOutEdgeIterator;
	typedef typename VertexDescriptor<TAlignmentGraph>::Type TVertexDescriptor;

	// iterate over sequences
	for (TSize id = 0; id < length(segments)-1; ++id) {
		//iterate over vertices of segment
		TPosition pos = beginPosition(value(segments, id));
		while (pos < endPosition(value(segments, id))) {
			TVertexDescriptor source = findVertex(parentGraph, id, pos);
			SEQAN_ASSERT_LEQ(pos+fragmentLength(parentGraph, source), endPosition(segments[id]));

			// iterate over neighbors of vertex
			TOutEdgeIterator edgeIt(parentGraph, source);
			while (!atEnd(edgeIt)) {
				TVertexDescriptor target = targetVertex(edgeIt);
				TSize targetId = sequenceId(parentGraph, target);
				if (targetId <= id) { // corresponding match was already appended
					++edgeIt;
					continue;
				}

				// initialize a match
				TMatch m;
				resize(rows(m), 2);
				setSource(row(m, 0), value(stringSet(parentGraph), id));
				setSource(row(m, 1), value(stringSet(parentGraph), targetId));

				// determine begin and end position of source and target segments
				TPosition sourceBegin = fragmentBegin(parentGraph, source);
				TPosition sourceEnd = sourceBegin + fragmentLength(parentGraph, source);
				TPosition targetBegin = fragmentBegin(parentGraph, target);
				TPosition targetEnd = targetBegin + fragmentLength(parentGraph, target);

				// set the match begin and end positions
				setClippedBeginPosition(row(m, 0), sourceBegin);
				setClippedBeginPosition(row(m, 1), targetBegin);
				setBeginPosition(row(m, 0), 0);
				setBeginPosition(row(m, 1), 0);
				setClippedEndPosition(row(m, 0), sourceEnd);
				setClippedEndPosition(row(m, 1), targetEnd);
				
				std::cout << sourceBegin << ".." << sourceEnd << " , " << targetBegin << ".." << targetEnd << std::endl;
				//std::cout << m;
				
				// append match
				appendValue(matches, m);
				++edgeIt;
			}
			pos += fragmentLength(parentGraph, source);
		}
	}
}

template<typename TSegmentSet, typename TAlignmentGraph>
void
_fixHigherLevelMatches(TSegmentSet & segments, TAlignmentGraph & parentGraph, TAlignmentGraph & g) {
	typedef typename Size<TSegmentSet>::Type TSize;
	typedef typename Position<typename Value<TSegmentSet>::Type>::Type TPosition;
	typedef typename Iterator<TAlignmentGraph, OutEdgeIterator>::Type TOutEdgeIterator;
	typedef typename VertexDescriptor<TAlignmentGraph>::Type TVertexDescriptor;
	typedef typename EdgeDescriptor<TAlignmentGraph>::Type TEdgeDescriptor;

	// determine the weight for edges to be fixed: the length of the longest sequence segment
	TSize maxWeight = 0;
	for (TSize id = 0; id < length(segments); ++id) {
		if (maxWeight < length(value(segments, id))) {
			maxWeight = length(value(segments, id));
		}
	}

	// iterate over sequences
	for (TSize id = 0; id < length(segments)-1; ++id) {
		//iterate over vertices of segment in parentGraph
		TPosition pos = beginPosition(value(segments, id));
		while (pos < endPosition(value(segments, id))) {
			TVertexDescriptor source = findVertex(parentGraph, id, pos);
			SEQAN_ASSERT_LEQ(pos+fragmentLength(parentGraph, source), endPosition(segments[id]));

			// iterate over neighbors of vertex
			TOutEdgeIterator edgeIt(parentGraph, source);
			while (!atEnd(edgeIt)) {
				TSize targetId = sequenceId(parentGraph, targetVertex(edgeIt));
				TPosition targetPos = fragmentBegin(parentGraph, targetVertex(edgeIt));
				
				// find corresponding edges in g
				TSize len = 0;
				do {
					TVertexDescriptor newSource = findVertex(g, id, pos + len);
					TVertexDescriptor newTarget = findVertex(g, targetId, targetPos + len);
					len += fragmentLength(g, newSource);

					// find edge and assign maximal weight
					TEdgeDescriptor e = findEdge(g, newSource, newTarget);
					if (e == 0) continue; // edge was not found
					assignCargo(e, maxWeight);
				} while (len < fragmentLength(parentGraph, source));
				SEQAN_ASSERT_EQ(len, fragmentLength(parentGraph, source));
				++edgeIt;
			}
			pos += fragmentLength(parentGraph, source);
		}
	}
}

template<typename TSegmentSet, typename TOptions, typename TLevel, typename TTree>
void
_computeGuideTree(TSegmentSet & segments, TOptions & options, TLevel recursionLevel, TTree & guideTree) {
	typedef String<double> TDistanceMatrix;
	typedef typename Iterator<TDistanceMatrix, Standard>::Type TMatrixIterator;
	typedef typename Value<typename Value<TSegmentSet>::Type>::Type TAlphabet;
	
	if (options.globalGuideTree && recursionLevel > 0) {
		guideTree = options.guideTree;
	} else {
		// compute distance matrix
		TDistanceMatrix distanceMatrix;		
		/*getDistanceMatrix(g, distanceMatrix, KmerDistance());*/
		getKmerSimilarityMatrix(segments, distanceMatrix, 3, TAlphabet());
	
		// similarity to distance conversion
		TMatrixIterator matIt = begin(distanceMatrix, Standard());
		TMatrixIterator endMatIt = end(distanceMatrix, Standard());
		for(;matIt != endMatIt;++matIt) 
			*matIt = SEQAN_DISTANCE_UNITY - (*matIt);

		// compute tree from distance matrix
		upgmaTree(distanceMatrix, guideTree);
		std::cout << guideTree;
		if (options.globalGuideTree) {
			options.guideTree = guideTree;
		}
	}
}

template<typename TAlignmentGraph>
void
flushNodes(TAlignmentGraph & g) {
	typedef typename VertexDescriptor<TAlignmentGraph>::Type TVertexDescriptor;
	typedef typename EdgeDescriptor<TAlignmentGraph>::Type TEdgeDescriptor;
	typedef typename Cargo<TAlignmentGraph>::Type TCargo;
	typedef typename Iterator<TAlignmentGraph,OutEdgeIterator>::Type TOutEdgeIterator;

	typedef typename StringSetType<TAlignmentGraph>::Type TStringSet;
	typedef typename Size<TStringSet>::Type TSize;
	typedef typename Position<typename Value<TStringSet>::Type>::Type TPosition;

	typedef Pair<TVertexDescriptor, TVertexDescriptor> TVertexPair;
	typedef Pair<TSize, TPosition> TVertexInfo; // seqId, beginPos
	typedef Triple<TVertexInfo, TVertexInfo, TCargo> TEdgeInfo; // seqId, seqId, cargo

	for (TSize i = 0; i < length(stringSet(g)); ++i) {
		TPosition pos = beginPosition(value(stringSet(g), i));
		TVertexDescriptor nextVertex = findVertex(g, i, pos);
		pos += fragmentLength(g, nextVertex);
		while (pos < endPosition(value(stringSet(g), i))) {
			TVertexDescriptor vertex = nextVertex;
			nextVertex = findVertex(g, i, pos);
			pos += fragmentLength(g, nextVertex);

			if (outDegree(g, vertex) != outDegree(g, nextVertex)) continue;

			String<TVertexPair> vertices;
			String<TEdgeInfo> edges;
			appendValue(vertices, TVertexPair(vertex, nextVertex));

			TOutEdgeIterator outEdgeIt(g, vertex);
			bool sameTargets = true;
			while (!atEnd(outEdgeIt)) {
				// check whether edge exists between nextVertex and next of target vertex
				TVertexDescriptor target = targetVertex(outEdgeIt);
				TPosition targetEndPos = fragmentBegin(g, target) + fragmentLength(g, target);
				TVertexDescriptor targetNext = findVertex(g, sequenceId(g, target), targetEndPos);
				TEdgeDescriptor e = findEdge(g, nextVertex, targetNext);
				if (e == 0) {
					sameTargets = false;
					break;
				}

				// append vertex and edge information to lists
				appendValue(vertices, TVertexPair(target, targetNext));
				appendValue(edges, TEdgeInfo(TVertexInfo(i, fragmentBegin(g, vertex)),
					                         TVertexInfo(sequenceId(g, target), fragmentBegin(g,target)),
											 getCargo(*outEdgeIt) + getCargo(e)));
				++outEdgeIt;
			}
			if (!sameTargets) continue;
			
			// check that all edges between target vertices exist and append edgeInfos to list
			TSize numVertices = length(vertices);
			for (TSize j = 1; j < numVertices; ++j) {
				for (TSize k = j+1; k < numVertices; k++) {
					TEdgeDescriptor edge = findEdge(g, vertices[j].i1, vertices[k].i1);
					TEdgeDescriptor nextEdge = findEdge(g, vertices[j].i2, vertices[k].i2);
					SEQAN_ASSERT(edge != 0);
					SEQAN_ASSERT(nextEdge != 0);
					appendValue(edges,
						       TEdgeInfo(TVertexInfo(sequenceId(g, vertices[j].i1), fragmentBegin(g, vertices[j].i1)),
						                 TVertexInfo(sequenceId(g, vertices[k].i1), fragmentBegin(g, vertices[k].i1)),
										 getCargo(edge) + getCargo(nextEdge)));
				}				
			}

			// remove two old vertices per sequence and add one new vertex
			for (TSize j = 1; j < numVertices; ++j) {
				TSize seqId = sequenceId(g, vertices[j].i1);
				TPosition beginPos = fragmentBegin(g, vertices[j].i1);
				TPosition fragLen = fragmentLength(g, vertices[j].i1) + fragmentLength(g, vertices[j].i2);
				removeVertex(g, vertices[j].i1);
				removeVertex(g, vertices[j].i2);
				TVertexDescriptor v = addVertex(g, seqId, beginPos, fragLen);
				if (i == j) nextVertex = v;
			}

			// add new edges with cargo
			for (TSize j = 0; j < length(edges); ++j) {
				TVertexDescriptor v = findVertex(g, edges[j].i1.i1, edges[j].i1.i2);
				TVertexDescriptor w = findVertex(g, edges[j].i2.i1, edges[j].i2.i2);
				addEdge(g, v, w, edges[j].i3);
			}
		}
	}
}

//template<typename TAlignmentGraph>
//void
//flushNodes(TAlignmentGraph & g) {
//	typedef typename VertexDescriptor<TAlignmentGraph>::Type TVertexDescriptor;
//	typedef typename Iterator<TAlignmentGraph,OutEdgeIterator>::Type TOutEdgeIter;
//    typedef typename StringSetType<TAlignmentGraph>::Type TStringSet;
//    typedef typename Infix<typename Value<TStringSet>::Type>::Type TSegment;
//	typedef typename Size<TSegment>::Type TSize;
//	typedef String<TVertexDescriptor> TVertices;
//
//	for (TSize i = 0; i < length(stringSet(g)); i++) {
//		TSegment seg = value(stringSet(g), i);
//		
//		TVertexDescriptor prev = findVertex(g, i, 0);
//		TSize pos = fragmentLength(g, prev);
//		while (pos < length(seg)) {
//			TVertexDescriptor next = findVertex(g, i, pos);
//			pos += fragmentLength(g, next);
//			
//			if (outDegree(g,prev) != outDegree(g, next)) {
//			} else {
//				TVertices targets;
//
//				TOutEdgeIter prevEdgeIt(g, prev);
//				while (!atEnd(prevEdgeIt)) {
//					TVertexDescriptor target = targetVertex(prevEdgeIt);
//					appendValue(targets, target);
//					TVertexDescriptor targetNext = findVertex(g, sequenceId(g, target),
//						fragmentBegin(g, target) + fragmentLength(g, target));
//					if (findEdge(g, next, targetNext) == 0) break;
//					goNext(prevEdgeIt);
//				}
//				if (atEnd(prevEdgeIt)) {
//					TVertices newTargets;
//					for (TSize t = 0; t < length(targets); t++) {
//						TVertexDescriptor target = value(targets, t);
//						TVertexDescriptor targetNext = findVertex(g, sequenceId(g, target),
//							fragmentBegin(g, target) + fragmentLength(g, target));
//
//						removeVertex(g, target);
//						removeVertex(g, targetNext);
//
//						TVertexDescriptor newTarget = addVertex(g, sequenceId(g, target), fragmentBegin(g, target),
//							fragmentLength(g,target) + fragmentLength(g, targetNext));
//						appendValue(newTargets, newTarget);
//					}
//					removeVertex(g, prev);
//					removeVertex(g, next);
//
//					TVertexDescriptor newVertex = addVertex(g, i, fragmentBegin(g, prev), 
//						fragmentLength(g, prev) + fragmentLength(g, next));
//
//					for (TSize t = 0; t < length(newTargets); t++) {
//						addEdge(g, newVertex, value(newTargets, t));
//						for (TSize t2 = t + 1; t2 < length(newTargets); t2++) {
//							addEdge(g, value(newTargets, t), value(newTargets ,t2));
//						}
//					}
//
//					next = newVertex;
//				}
//			}
//			prev = next;
//		}
//	}
//}

template<typename TSegmentSet, typename TAlignmentGraph>
void
_segmentGraph(TSegmentSet & segments, TAlignmentGraph & g) {
	typedef typename Size<TSegmentSet>::Type TSize;
	typedef typename VertexDescriptor<TAlignmentGraph>::Type TVertexDescriptor;

	for (TSize i = 0; i < length(segments); ++i) {
		typename Value<TSegmentSet>::Type seg = value(segments, i);
		TSize len = length(host(seg));

		// additional vertex for unconsidered leading sequence part
		if (beginPosition(seg) != 0) {
			TVertexDescriptor v = findVertex(g, i, 0);
			SEQAN_ASSERT_EQ(outDegree(g, v), 0u);
			SEQAN_ASSERT_EQ(fragmentBegin(g, v), 0u);
			TSize fragLen = fragmentLength(g, v);
			if (fragLen > beginPosition(seg)) {
				removeVertex(g, v);
				addVertex(g, i, 0, beginPosition(seg));
				addVertex(g, i, beginPosition(seg), fragLen - beginPosition(seg));
			}
		}
		// additional vertex for unconsidered trailing sequence part
		if (endPosition(seg) != len) {
			TVertexDescriptor v = findVertex(g, i, len-1);
			SEQAN_ASSERT_EQ(outDegree(g, v), 0u);
			SEQAN_ASSERT_EQ(fragmentBegin(g, v) + fragmentLength(g, v), len);
			TSize fragBegin = fragmentBegin(g, v);
			if (fragBegin < endPosition(seg)) {
				removeVertex(g, v);
				addVertex(g, i, fragBegin, endPosition(seg) - fragBegin);
				addVertex(g, i, endPosition(seg), len - endPosition(seg));
			}
		}
	}
}

// segment alignment on gOut given a string of matches as segment input
template<typename TSegmentSet, typename TOptions, typename TLevel, typename TAlignmentGraph>
void segmentAlignment(TSegmentSet & segments,
					  TAlignmentGraph & parentGraph,
					  TOptions & options,
					  TLevel & recursionLevel,
					  TAlignmentGraph & gOut) {
	typedef typename Value<TSegmentSet>::Type TSegment;
	typedef Align<typename Host<TSegment>::Type> TMatch;

	typedef typename VertexDescriptor<TAlignmentGraph>::Type TVertexDescriptor;
	typedef typename EdgeDescriptor<TAlignmentGraph>::Type TEdgeDescriptor;

	typedef typename Size<TAlignmentGraph>::Type TSize;

	clearEdges(gOut);
	clearVertices(gOut);

	//for( unsigned i = 0; i < length(segments); ++i) {
		std::cout << beginPosition(segments[0]) << " - " << endPosition(segments[0]) << std::endl;
	//}

	// generate local alignments
	String<TMatch> matches;
	if (options.anchoredPairwiseComparison && recursionLevel != 0) {
		findStellarMatches(segments, parentGraph, options, recursionLevel, matches);
	} else {
		findStellarMatches(segments, options, recursionLevel, matches);
	}

	// copy matches from segments of parent graph to the new set of matches
	if (recursionLevel != 0) {
		_addParentMatches(segments, parentGraph, matches);
	}

    // match refinement
	TAlignmentGraph g(stringSet(gOut));
	matchRefinement(matches, stringSet(gOut), g);
	clear(matches);
	
    // triplet extension
	tripletLibraryExtension(g);

	if (options.fixedHigherLevelMatches && recursionLevel != 0) {
		// set the weight of matches from parentGraph to infinity
		_fixHigherLevelMatches(segments, parentGraph, g);
	}

    // guide tree computation
	typename TOptions::TTree guideTree;
	_computeGuideTree(segments, options, recursionLevel, guideTree);
	
    // progressive alignment
	progressiveAlignment(g, guideTree, gOut);

    // combine neighboring vertices if edges are the same
	flushNodes(gOut);

	// add separate vertices for leading and trailing sequence segments
	_segmentGraph(segments, gOut);
}

template<typename TAlignmentGraph, typename TInfixSet, typename TPosition>
String<TPosition>
findUnalignedSegment(TAlignmentGraph & g, TInfixSet const & segments, String<TPosition> startPos, TInfixSet & unalignedSegments) {
	typedef typename VertexDescriptor<TAlignmentGraph>::Type TVertexDescriptor;
	typedef typename Iterator<TAlignmentGraph, OutEdgeIterator>::Type TOutEdgeIter;

	// find the next cut (in first sequence: seqId = 0)
	TPosition seqId = 0;
	TVertexDescriptor v = findVertex(g, seqId, value(startPos, seqId));
	SEQAN_ASSERT_EQ(startPos[seqId], fragmentBegin(g, v));
	TPosition endPos = value(startPos, seqId) + fragmentLength(g, v);
	while (outDegree(g, v) < length(segments)-1 && endPos < endPosition(value(segments, seqId))) {
		v = findVertex(g, seqId, endPos);
		SEQAN_ASSERT_EQ(endPos, fragmentBegin(g, v));
		endPos += fragmentLength(g, v);
	}

	String<TPosition> nextStartPos;
	resize(nextStartPos, length(startPos));
	resize(unalignedSegments, length(startPos));

	// ends of segments reached
	if (outDegree(g, v) < length(segments)-1) {
		for (TPosition i = 0; i < length(segments); ++i) {
			value(unalignedSegments, i) = infix(host(segments[i]), startPos[i], endPosition(segments[i]));
			value(nextStartPos, i) = endPosition(segments[i]);
		}
		return nextStartPos;
	}

	// determine start and end position of unaligned segment in all sequences
	TOutEdgeIter it(g, v);
	value(unalignedSegments, seqId) = infix(host(segments[seqId]), startPos[seqId], endPos - fragmentLength(g, v));
	value(nextStartPos, seqId) = endPos;
	while (!atEnd(it)) {
		TVertexDescriptor u = targetVertex(it);
		SEQAN_ASSERT_NEQ(u, v);

		TPosition id = sequenceId(g, u);
		TPosition start = value(startPos, id);
		TPosition end = fragmentBegin(g, u);
		
		value(unalignedSegments, id) = infix(host(segments[id]), start, end);
		value(nextStartPos, id) = end + fragmentLength(g, u);

		goNext(it);
	}

	return nextStartPos;
}

template<typename TDepSeqSet, typename TSegmentSet>
void
integrateAlignmentGraph(Graph<Alignment<TDepSeqSet> > & parentGraph,
						Graph<Alignment<TDepSeqSet> > & graph,
						TSegmentSet & segments) {
	typedef Graph<Alignment<TDepSeqSet> > TAlignmentGraph;
	typedef typename VertexDescriptor<TAlignmentGraph>::Type TVertexDescriptor;
	typedef typename Iterator<TAlignmentGraph, VertexIterator>::Type TVertexIterator;
	typedef typename Iterator<TAlignmentGraph, EdgeIterator>::Type TEdgeIterator;

	typedef typename Size<TDepSeqSet>::Type TSize;
	typedef typename Position<TDepSeqSet>::Type TPosition;

	// remove vertices on segments from parentGraph
	for (TSize i = 0; i < length(segments); ++i) {
		TPosition pos = beginPosition(segments[i]);
		while (pos < endPosition(segments[i])) {
			TVertexDescriptor v = findVertex(parentGraph, i, pos);
			pos += fragmentLength(parentGraph, v);
			removeVertex(parentGraph, v);
		}
	}

	// copy vertices from graph to parentGraph
	TVertexIterator vertexIt(graph);
	while (!atEnd(vertexIt)) {
		TPosition seqId = sequenceId(graph, *vertexIt);
		TPosition pos = fragmentBegin(graph, *vertexIt);
		TSize len = fragmentLength(graph, *vertexIt);

		if (pos + len > endPosition(segments[seqId])) {
			SEQAN_ASSERT_EQ(pos, endPosition(segments[seqId]));
			SEQAN_ASSERT_EQ(pos+len, length(host(segments[seqId])));
		} else if (pos < beginPosition(segments[seqId])) {
			SEQAN_ASSERT_EQ(pos, 0u);
			SEQAN_ASSERT_EQ(pos+len, beginPosition(segments[seqId]));
		} else {
			addVertex(parentGraph, seqId, pos, len);
		}
		++vertexIt;
	}

	// copy edges from graph to parentGraph
	TEdgeIterator edgeIt(graph);
	while (!atEnd(edgeIt)) {
		TVertexDescriptor source = sourceVertex(graph, *edgeIt);
		TVertexDescriptor target = targetVertex(graph, *edgeIt);

		TPosition sourceSeqId = sequenceId(graph, source);
		TPosition targetSeqId = sequenceId(graph, target);

		TPosition sourcePos = fragmentBegin(graph, source);
		TPosition targetPos = fragmentBegin(graph, target);

		source = findVertex(parentGraph, sourceSeqId, sourcePos);
		target = findVertex(parentGraph, targetSeqId, targetPos);

		addEdge(parentGraph, source, target, getCargo(*edgeIt));
		++edgeIt;
	}
}

// computes segmentalignment and recurses for decreasing matchMinLengths on unaligned subgraphs
template<typename TSegmentSet, typename TOptions, typename TLevel, typename TAlignmentGraph>
void
recurseSegmentAlignment(TSegmentSet segments,
						TAlignmentGraph & parentGraph,
						TOptions & options,
						TLevel recursionLevel,
						TAlignmentGraph & g) {
    typedef typename Value<TSegmentSet>::Type TInfix;
    typedef typename Position<TInfix>::Type TPosition;

	segmentAlignment(segments, parentGraph, options, recursionLevel, g);
	//if (recursionLevel == 0) std::cout << g;

	if(options.recursions <= recursionLevel + 1) {
		return;
	}

	String<TPosition> pos;
	for (TPosition i = 0; i < length(segments); ++i) {
		appendValue(pos, beginPosition(value(segments, i)));
	}

	while (pos[0] < endPosition(value(segments, 0))) {
		TSegmentSet unalignedSegments;
		pos = findUnalignedSegment(g, segments, pos, unalignedSegments);

		//for( int space = 0; space < recursionLevel; ++space) std::cout << ".";
		//std::cout << beginPosition(unalignedSegments[0]) << " - " << endPosition(unalignedSegments[0]) << std::endl;

		TAlignmentGraph subgraph(stringSet(g));
		recurseSegmentAlignment(unalignedSegments, g, options, recursionLevel+1, subgraph);
		integrateAlignmentGraph(g, subgraph, unalignedSegments);
	}
}

#endif
