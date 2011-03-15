
#ifndef SEQAN_HEADER_HSA_H
#define SEQAN_HEADER_HSA_H

#include "hsa_types.h"
#include "hsa_localAlignment.h"

using namespace seqan;

template<typename TId, typename TInfix, typename TAlignmentGraph, typename TMatch>
void
_addParentMatches(std::map<TId, TInfix> & segments, TAlignmentGraph & parentGraph, String<TMatch> & matches) {
	typedef typename std::map<TId, TInfix>::const_iterator TMapIterator;
	typedef typename Position<TInfix>::Type TPosition;
	typedef typename Id<TAlignmentGraph>::Type TGraphId;
	typedef typename Iterator<TAlignmentGraph, OutEdgeIterator>::Type TOutEdgeIterator;
	typedef typename VertexDescriptor<TAlignmentGraph>::Type TVertexDescriptor;

	// iterate over sequences
	TMapIterator itEnd = segments.end();
	for (TMapIterator it = segments.begin(); it != itEnd; ++it) {
		//iterate over vertices of segment
		TPosition pos = beginPosition(it->second);
		while (pos < endPosition(it->second)) {
			TGraphId i = stringIdToPosition(stringSet(parentGraph), it->first);
			TVertexDescriptor source = findVertex(parentGraph, i, pos);
			SEQAN_ASSERT_LEQ(pos+fragmentLength(parentGraph, source), endPosition(it->second));

			// iterate over neighbors of vertex
			TOutEdgeIterator edgeIt(parentGraph, source);
			while (!atEnd(edgeIt)) {
				TVertexDescriptor target = targetVertex(edgeIt);
				TGraphId targetId = sequenceId(parentGraph, target);
				if (targetId <= i) { // corresponding match was already appended
					++edgeIt;
					continue;
				}

				// initialize a match
				TMatch m;
				resize(rows(m), 2);
				setSource(row(m, 0), value(stringSet(parentGraph), i));
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
				
				//std::cout << sourceBegin << ".." << sourceEnd << " , " << targetBegin << ".." << targetEnd << std::endl;
				//std::cout << m;
				
				// append match
				appendValue(matches, m);
				++edgeIt;
			}
			pos += fragmentLength(parentGraph, source);
		}
	}
}

template<typename TId, typename TInfix, typename TAlignmentGraph>
void
_fixHigherLevelMatches(std::map<TId, TInfix> & segments, TAlignmentGraph & parentGraph, TAlignmentGraph & g) {
	typedef typename std::map<TId, TInfix>::const_iterator TMapIterator;
	typedef typename Size<TInfix>::Type TSize;
	typedef typename Position<TInfix>::Type TPosition;
	typedef typename Id<TAlignmentGraph>::Type TGraphId;
	typedef typename Iterator<TAlignmentGraph, OutEdgeIterator>::Type TOutEdgeIterator;
	typedef typename VertexDescriptor<TAlignmentGraph>::Type TVertexDescriptor;
	typedef typename EdgeDescriptor<TAlignmentGraph>::Type TEdgeDescriptor;

	if (segments.size() == 0) return;

	// determine the weight for edges to be fixed: the length of the longest sequence segment
	TSize maxWeight = 0;
	TMapIterator itEnd = segments.end();
	for (TMapIterator it = segments.begin(); it != itEnd; ++it) {
		if (maxWeight < length(it->second)) {
			maxWeight = length(it->second);
		}
	}

	// iterate over sequences
	--itEnd;
	for (TMapIterator it = segments.begin(); it != itEnd; ++it) {
		//iterate over vertices of segment in parentGraph
		TPosition pos = beginPosition(it->second);
		while (pos < endPosition(it->second)) {
			TGraphId idParent = stringIdToStringSetId(stringSet(parentGraph), it->first);
			TVertexDescriptor source = findVertex(parentGraph, idParent, pos);
			SEQAN_ASSERT_LEQ(pos+fragmentLength(parentGraph, source), endPosition(it->second));

			// iterate over neighbors of vertex
			TOutEdgeIterator edgeIt(parentGraph, source);
			while (!atEnd(edgeIt)) {
				TGraphId targetIdParent = sequenceId(parentGraph, targetVertex(edgeIt));
				TId targetId = id(value(stringSet(parentGraph), idToPosition(stringSet(parentGraph), targetIdParent)));
				TPosition targetPos = fragmentBegin(parentGraph, targetVertex(edgeIt));
				
				// find corresponding edges in g
				TSize len = 0;
				do {
					TVertexDescriptor newSource = findVertex(g, stringIdToStringSetId(stringSet(g), it->first), pos + len);
					TVertexDescriptor newTarget = findVertex(g, stringIdToStringSetId(stringSet(g), targetId), targetPos + len);

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

template<typename TSequenceSet, typename TId, typename TInfix, typename TOptions, typename TLevel, typename TTree>
void
_computeGuideTree(Graph<Alignment<TSequenceSet> > & g,
				  std::map<TId, TInfix> & segments,
				  TOptions & options,
				  TLevel recursionLevel,
				  TTree & guideTree) {
	typedef String<double> TDistanceMatrix;
	typedef typename Iterator<TDistanceMatrix>::Type TMatrixIterator;
	typedef typename Size<TSequenceSet>::Type TSize;

	typedef typename VertexDescriptor<TTree>::Type TVertexDescriptor;
	typedef typename Iterator<TTree, VertexIterator>::Type TVertexIterator;
	
	typedef typename Value<TInfix>::Type TAlphabet;
	
	if (options.globalGuideTree && recursionLevel > 0) {
		guideTree = options.guideTree;

		// reduce global tree to tree on stringSet(g)
		for (TVertexIterator it(guideTree); !atEnd(it); ++it) {
			if (!isLeaf(guideTree, *it) || segments.count(options.idMap[*it]) == 1) continue;
			
			// remove vertex from guideTree
			TVertexDescriptor v = *it;
			TVertexDescriptor parent_v = parentVertex(guideTree, v);
			removeVertex(guideTree, v);

			// connect sibling and grandparent of v by an edge or assign new root
			SEQAN_ASSERT_EQ(outDegree(guideTree, parent_v), 1u);
			typename Iterator<TTree, OutEdgeIterator>::Type edge(guideTree, parent_v);
			TVertexDescriptor sibling_v = childVertex(guideTree, *edge);
			if (!isRoot(guideTree, parent_v)) {
				TVertexDescriptor grandparent_v = parentVertex(guideTree, parent_v);
				addEdge(guideTree, sibling_v, grandparent_v);
			} else {
				assignRoot(guideTree, sibling_v);
			}
			removeVertex(guideTree, parent_v);
		}
	} else {
		// compute distance matrix
		TDistanceMatrix distanceMatrix;		
		/*getDistanceMatrix(g, distanceMatrix, KmerDistance());*/
		StringSet<TInfix> segs;
		resize(segs, segments.size());
		for (TSize i = 0; i < length(stringSet(g)); ++i) {
			value(segs, positionToId(stringSet(g), i)) = segments[id(value(stringSet(g), i))];
		}
		getKmerSimilarityMatrix(segs, distanceMatrix, 3, TAlphabet());
		
		// Similarity to distance conversion
		TMatrixIterator matIt = begin(distanceMatrix, Standard());
		TMatrixIterator endMatIt = end(distanceMatrix, Standard());
		for(;matIt != endMatIt;++matIt) 
			*matIt = SEQAN_DISTANCE_UNITY - (*matIt);

		// compute tree from distance matrix
		upgmaTree(distanceMatrix, guideTree);

		if (options.globalGuideTree) {
			// compute map seqId -> graphSeqIds
			for (TSize i = 0; i < length(stringSet(g)); ++i) {
				options.idMap[positionToId(stringSet(g), i)] = id(value(stringSet(g), i));
			}
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
	typedef typename Id<TAlignmentGraph>::Type TId;
	typedef typename Size<TStringSet>::Type TSize;
	typedef typename Position<typename Value<TStringSet>::Type>::Type TPosition;

	typedef Pair<TVertexDescriptor, TVertexDescriptor> TVertexPair;
	typedef Pair<TId, TPosition> TVertexInfo; // seqId, beginPos
	typedef Pair<TVertexInfo> TEdgeInfo; // source, target

	for (TSize i = 0; i < length(stringSet(g)); ++i) {
		TId id_i = positionToId(stringSet(g), i);
		TPosition pos = beginPosition(value(stringSet(g), i));
		TVertexDescriptor nextVertex = findVertex(g, id_i, pos);
		pos += fragmentLength(g, nextVertex);
		while (pos < endPosition(value(stringSet(g), i))) {
			TVertexDescriptor vertex = nextVertex;
			nextVertex = findVertex(g, id_i, pos);
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
				appendValue(edges, TEdgeInfo(TVertexInfo(id_i, fragmentBegin(g, vertex)),
					                         TVertexInfo(sequenceId(g, target), fragmentBegin(g,target))));
				++outEdgeIt;
			}
			if (!sameTargets) continue;
			
			// check that all edges between target vertices exist and append edgeInfos to list
			TSize numVertices = length(vertices);
			for (TSize j = 1; j < numVertices; ++j) {
				for (TSize k = j+1; k < numVertices; k++) {
					SEQAN_ASSERT(findEdge(g, vertices[j].i1, vertices[k].i1) != 0);
					SEQAN_ASSERT(findEdge(g, vertices[j].i2, vertices[k].i2) != 0);
					appendValue(edges,
						        TEdgeInfo(TVertexInfo(sequenceId(g, vertices[j].i1), fragmentBegin(g, vertices[j].i1)),
						                  TVertexInfo(sequenceId(g, vertices[k].i1), fragmentBegin(g, vertices[k].i1))));
				}				
			}

			// remove two old vertices per sequence and add one new vertex
			for (TSize j = 0; j < numVertices; ++j) {
				TId seqId = sequenceId(g, vertices[j].i1);
				TPosition beginPos = fragmentBegin(g, vertices[j].i1);
				TPosition fragLen = fragmentLength(g, vertices[j].i1) + fragmentLength(g, vertices[j].i2);
				removeVertex(g, vertices[j].i1);
				removeVertex(g, vertices[j].i2);
				TVertexDescriptor v = addVertex(g, seqId, beginPos, fragLen);
				if (i == j) nextVertex = v;
			}

			// add new edges
			for (TSize j = 0; j < length(edges); ++j) {
				TVertexDescriptor v = findVertex(g, edges[j].i1.i1, edges[j].i1.i2);
				TVertexDescriptor w = findVertex(g, edges[j].i2.i1, edges[j].i2.i2);
				addEdge(g, v, w);
			}
		}
	}
}

template<typename TId, typename TInfix, typename TAlignmentGraph>
void
_segmentGraph(std::map<TId, TInfix> & segments, TAlignmentGraph & g) {
	typedef typename Size<TInfix>::Type TSize;
	typedef typename std::map<TId, TInfix>::const_iterator TMapIterator;
	typedef typename VertexDescriptor<TAlignmentGraph>::Type TVertexDescriptor;

	TMapIterator itEnd = segments.end();
	for (TMapIterator it = segments.begin(); it != itEnd; ++it) {
		typename Id<TAlignmentGraph>::Type i = stringIdToPosition(stringSet(g), it->first);
		TInfix seg = it->second;
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
template<typename TId, typename TInfix, typename TOptions, typename TLevel, typename TAlignmentGraph>
void segmentAlignment(std::map<TId, TInfix> & segments,
					  TAlignmentGraph & parentGraph,
					  TOptions & options,
					  TLevel & recursionLevel,
					  TAlignmentGraph & gOut) {
	typedef Align<typename Host<TInfix>::Type> TMatch;

	typedef typename VertexDescriptor<TAlignmentGraph>::Type TVertexDescriptor;
	typedef typename EdgeDescriptor<TAlignmentGraph>::Type TEdgeDescriptor;

	typedef typename Size<TAlignmentGraph>::Type TSize;

	clearEdges(gOut);
	clearVertices(gOut);

	std::cout << beginPosition(segments.begin()->second) << " - " << endPosition(segments.begin()->second) << std::endl;

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
		// set high weight on matches from parentGraph
		_fixHigherLevelMatches(segments, parentGraph, g);
	}

    // guide tree computation
	typename TOptions::TTree guideTree;
	_computeGuideTree(g, segments, options, recursionLevel, guideTree);
	
    // progressive alignment
	progressiveAlignment(g, guideTree, gOut);

    // combine neighboring vertices if edges are the same
	flushNodes(gOut);

	// add separate vertices for leading and trailing sequence segments
	_segmentGraph(segments, gOut);
}

template<typename TSequenceSet, typename TId, typename TInfix, typename TLength, typename TPosition>
std::map<TId, TPosition>
findUnalignedSegment(Graph<Alignment<TSequenceSet> > & g,
					 std::map<TId, TInfix> & segments,
					 TLength & minLength,
					 std::map<TId, TPosition> & startPos,
					 std::map<TId, TInfix> & unalignedSegments) {
	typedef typename Id<Graph<Alignment<TSequenceSet> > >::Type TGraphId;
	typedef typename VertexDescriptor<Graph<Alignment<TSequenceSet> > >::Type TVertexDescriptor;
	typedef typename Iterator<Graph<Alignment<TSequenceSet> >, OutEdgeIterator>::Type TOutEdgeIter;
	typedef typename std::map<TId, TInfix>::const_iterator TMapIterator;

	// find the next cut (in first sequence: seqId = 0)
	TId seqId = segments.begin()->first;
	TGraphId i = stringIdToPosition(stringSet(g), seqId);
	TVertexDescriptor v = findVertex(g, i, startPos[seqId]);
	SEQAN_ASSERT_EQ(startPos[seqId], fragmentBegin(g, v));
	TPosition endPos = startPos[seqId] + fragmentLength(g, v);
	while (outDegree(g, v) < length(segments)-1 && endPos < endPosition(segments[seqId])) {
		v = findVertex(g, i, endPos);
		SEQAN_ASSERT_EQ(endPos, fragmentBegin(g, v));
		endPos += fragmentLength(g, v);
	}

	std::map<TId, TPosition> nextStartPos;

	// ends of segments reached
	if (outDegree(g, v) < length(segments)-1) {
		TMapIterator itEnd = segments.end();
		for (TMapIterator it = segments.begin(); it != itEnd; ++it) {
			if (startPos[it->first] + minLength < endPosition(it->second)) {
				unalignedSegments[it->first] = infix(host(it->second), startPos[it->first], endPosition(it->second));
			}
			nextStartPos[it->first] = endPosition(it->second);
		}
		return nextStartPos;
	}

	// determine start and end position of unaligned segment in all sequences
	TOutEdgeIter it(g, v);
	if (startPos[seqId] + minLength < endPos - fragmentLength(g, v)) {
		unalignedSegments[seqId] = infix(host(segments[seqId]), startPos[seqId], endPos - fragmentLength(g, v));
	}
	nextStartPos[seqId] = endPos;
	while (!atEnd(it)) {
		TVertexDescriptor u = targetVertex(it);
		SEQAN_ASSERT_NEQ(u, v);

		TId seq2Id = id(value(stringSet(g), idToPosition(stringSet(g), sequenceId(g, u))));
		TPosition start = startPos[seq2Id];
		TPosition end = fragmentBegin(g, u);
		
		if (start + minLength < end) {
			unalignedSegments[seq2Id] = infix(host(segments[seq2Id]), start, end);
		}
		nextStartPos[seq2Id] = end + fragmentLength(g, u);

		goNext(it);
	}

	return nextStartPos;
}

template<typename TDepSeqSet, typename TId, typename TInfix>
void
integrateAlignmentGraph(Graph<Alignment<TDepSeqSet> > & parentGraph,
						Graph<Alignment<TDepSeqSet> > & graph,
						std::map<TId, TInfix> & segments) {
	typedef Graph<Alignment<TDepSeqSet> > TAlignmentGraph;
	typedef typename Id<TAlignmentGraph>::Type TGraphId;
	typedef typename VertexDescriptor<TAlignmentGraph>::Type TVertexDescriptor;
	typedef typename Iterator<TAlignmentGraph, VertexIterator>::Type TVertexIterator;
	typedef typename Iterator<TAlignmentGraph, EdgeIterator>::Type TEdgeIterator;

	typedef typename std::map<TId, TInfix>::const_iterator TMapIterator;
	typedef typename Size<TInfix>::Type TSize;
	typedef typename Position<TDepSeqSet>::Type TPosition;

	// remove vertices on segments from parentGraph
	TMapIterator itEnd = segments.end();
	for (TMapIterator it = segments.begin(); it != itEnd; ++it) {
		TPosition pos = beginPosition(it->second);
		while (pos < endPosition(it->second)) {
			TGraphId i = stringIdToPosition(stringSet(parentGraph), it->first);
			TVertexDescriptor v = findVertex(parentGraph, i, pos);
			pos += fragmentLength(parentGraph, v);
			removeVertex(parentGraph, v);
		}
	}

	// copy vertices from graph to parentGraph
	TVertexIterator vertexIt(graph);
	while (!atEnd(vertexIt)) {
		TGraphId i = sequenceId(graph, *vertexIt);
		TId seqId = id(value(stringSet(graph), idToPosition(stringSet(graph), i)));
		TPosition pos = fragmentBegin(graph, *vertexIt);
		TSize len = fragmentLength(graph, *vertexIt);

		if (pos + len > endPosition(segments[seqId])) {
			SEQAN_ASSERT_EQ(pos, endPosition(segments[seqId]));
			SEQAN_ASSERT_EQ(pos+len, length(host(segments[seqId])));
		} else if (pos < beginPosition(segments[seqId])) {
			SEQAN_ASSERT_EQ(pos, 0u);
			SEQAN_ASSERT_EQ(pos+len, beginPosition(segments[seqId]));
		} else {
			addVertex(parentGraph, stringIdToStringSetId(stringSet(parentGraph), seqId), pos, len);
		}
		++vertexIt;
	}

	// copy edges from graph to parentGraph
	TEdgeIterator edgeIt(graph);
	while (!atEnd(edgeIt)) {
		TVertexDescriptor source = sourceVertex(graph, *edgeIt);
		TVertexDescriptor target = targetVertex(graph, *edgeIt);

		TGraphId sourceSeqId = sequenceId(graph, source);
		TId sourceId = id(value(stringSet(graph), idToPosition(stringSet(graph), sourceSeqId)));
		TGraphId targetSeqId = sequenceId(graph, target);
		TId targetId = id(value(stringSet(graph), idToPosition(stringSet(graph), targetSeqId)));

		TPosition sourcePos = fragmentBegin(graph, source);
		TPosition targetPos = fragmentBegin(graph, target);

		source = findVertex(parentGraph, stringIdToStringSetId(stringSet(parentGraph), sourceId), sourcePos);
		target = findVertex(parentGraph, stringIdToStringSetId(stringSet(parentGraph), targetId), targetPos);

		addEdge(parentGraph, source, target);
		++edgeIt;
	}
}

// computes segmentalignment and recurses for decreasing matchMinLengths on unaligned subgraphs
template<typename TId, typename TInfix, typename TSequenceSet, typename TOptions, typename TLevel>
void
recurseSegmentAlignment(std::map<TId, TInfix> & segments,
						Graph<Alignment<TSequenceSet> > & parentGraph,
						TOptions & options,
						TLevel recursionLevel,
						Graph<Alignment<TSequenceSet> > & g) {
	typedef typename Size<TInfix>::Type TSize;
    typedef typename Position<TInfix>::Type TPosition;
	typedef typename std::map<TId, TInfix>::const_iterator TMapIterator;

	segmentAlignment(segments, parentGraph, options, recursionLevel, g);

	if(options.recursions <= recursionLevel + 1) {
		return;
	}

	std::map<TId, TPosition> pos;
	TMapIterator itEnd = segments.end();
	for (TMapIterator it = segments.begin(); it != itEnd; ++it) {
		pos[it->first] = beginPosition(it->second);
	}

	while(pos.begin()->second < endPosition(segments.begin()->second)) {
		std::map<TId, TInfix> unalignedSegments;
		TSize minLength = options.initialMinLength - options.deltaMinLength * recursionLevel;
		pos = findUnalignedSegment(g, segments, minLength, pos, unalignedSegments);
		if (unalignedSegments.size() < 2) continue;

		TSequenceSet seqs;
		TMapIterator uEnd = unalignedSegments.end();
		for (TMapIterator uIt = unalignedSegments.begin(); uIt != uEnd; ++uIt) {
			appendValue(seqs, host(uIt->second));
		}
		Graph<Alignment<TSequenceSet> > subgraph(seqs);
		
		//for( int space = 0; space < recursionLevel; ++space) std::cout << ".";
		//std::cout << beginPosition(unalignedSegments[0]) << " - " << endPosition(unalignedSegments[0]) << std::endl;

		recurseSegmentAlignment(unalignedSegments, g, options, recursionLevel+1, subgraph);
		integrateAlignmentGraph(g, subgraph, unalignedSegments);
	}
}

#endif
