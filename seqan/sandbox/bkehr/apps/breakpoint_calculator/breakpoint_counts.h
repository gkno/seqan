// ==========================================================================
//                            breakpoint_counts.h
//                           breakpoint_calculator
// ==========================================================================
// Copyright (C) 2012 by Birte Kehr
//
// This program is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option)
// any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
// more details.
//
// You should have received a copy of the GNU General Public License along
// with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// ==========================================================================
// Author: Birte Kehr <birte.kehr@fu-berlin.de>
// ==========================================================================

#ifndef SANDBOX_BKEHR_APPS_BREAKPOINT_CALCULATOR_BREAKPOINT_COUNTS_H_
#define SANDBOX_BKEHR_APPS_BREAKPOINT_CALCULATOR_BREAKPOINT_COUNTS_H_

#include <seqan/align.h>

#include <lemon/lgf_writer.h>
#include <lemon/smart_graph.h>
#include <lemon/matching.h>


using namespace seqan;

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

struct BlockInMatchingGraph
{
	typedef lemon::SmartGraph::Node TNode;
	TNode head;
	TNode tail;
	TNode headTelomere;
	TNode tailTelomere;
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

template<typename TSeqId, typename TSize>
int
sequencesOfBlocks(std::map<TSeqId, String<int> > & blockSeqs,
				  String<std::map<CharString, AlignmentBlockRow<TSize, TSize> > > & idToRowMaps)
{
	typedef std::map<CharString, AlignmentBlockRow<TSize, TSize> >  TIdRowMap;

	// map of seq ids to maps of start positions to block numbers
	typedef std::map<CharString, std::map<TSize, int> > TSeqMaps;
	TSeqMaps sequenceMaps;

	for (TSize i = 0; i < length(idToRowMaps); ++i)
	{
		TIdRowMap map = value(idToRowMaps, i);

		typename TIdRowMap::const_iterator it = map.begin();
		while (it != map.end())
		{
			if ((*it).second.orientation)
				sequenceMaps[it->first][(*it).second.startPos] = i;
			else
				sequenceMaps[it->first][(*it).second.startPos] = -static_cast<int>(i);

			++it;
		}
	}

	TSize i = 0;
	typename TSeqMaps::const_iterator mapIt = sequenceMaps.begin();
	while (mapIt != sequenceMaps.end())
	{
		typename std::map<TSize, int>::const_iterator seqIt = (*mapIt).second.begin();
		while (seqIt != (*mapIt).second.end())
		{
			appendValue(blockSeqs[mapIt->first], seqIt->second);
			++seqIt;
		}
		++mapIt;
		++i;
	}

	//// Debug output
	//for (TSize i = 0; i < length(blockSeqs); ++i)
	//{
	//	for (TSize j = 0; j < length(blockSeqs[i]); ++j)
	//		std::cout << blockSeqs[i][j] << ",";
	//	std::cout << std::endl;
	//}
	return 0;
}

Position<CharString>::Type
lastOccOfChar(CharString const & str, char c)
{
	if (length(str) == 0) return 0;

	Iterator<CharString const, Rooted>::Type it = end(str, Rooted());
	Iterator<CharString const>::Type itBegin = begin(str);
	
	for (--it; it >= itBegin; --it)
		if (*it == c) return position(it);

	return length(str);
}

Position<CharString>::Type
firstOccOfChar(CharString const & str, char c)
{
	if (length(str) == 0) return 0;

	Iterator<CharString const, Rooted>::Type it = begin(str, Rooted());
	Iterator<CharString const>::Type itEnd = end(str);
	
	for (; it < itEnd; ++it)
		if (*it == c) return position(it);

	return length(str);
}

template<typename TBlockId>
void
collateChromosomes(std::map<CharString, StringSet<String<TBlockId>, Dependent<> > > & blockSeqSets,
				   std::map<CharString, String<TBlockId> > & blockSeqs)
{
	typedef std::map<CharString, String<TBlockId> > TSeqMap;
	typedef typename TSeqMap::const_iterator TIter;
	typedef typename Position<CharString>::Type TPos;

	if (length(blockSeqs) == 0) return;

	TIter it = blockSeqs.begin();
	TIter itEnd = blockSeqs.end();

	CharString prevGenomeId = prefix(it->first, firstOccOfChar(it->first, '.'));
	StringSet<String<TBlockId>, Dependent<> > chromosomes;
	appendValue(chromosomes, it->second);

	CharString prev = it->first;
	for (++it; it != itEnd; ++it)
	{
		CharString genomeId = prefix(it->first, firstOccOfChar(it->first, '.'));
		if (genomeId != prevGenomeId)
		{
			blockSeqSets[prevGenomeId] = chromosomes;
			clear(chromosomes);
			prevGenomeId = genomeId;
		}
		appendValue(chromosomes, it->second);
	}
	blockSeqSets[prevGenomeId] = chromosomes;
}

template<typename TBlockId, typename TSize>
void
commonBlocks(std::map<TBlockId, TSize> & set,
			 StringSet<String<TBlockId>, Dependent<> > const & seq1,
			 StringSet<String<TBlockId>, Dependent<> > const & seq2)
{
	typedef typename Iterator<ConcatenatorManyToOne<StringSet<String<TBlockId>, Dependent<> > > >::Type TIterator;

	std::set<TBlockId> seq1Set;

	// insert all blocks of seq1 to seq1Set
	TIterator end1 = end(concat(seq1));
	for (TIterator it = begin(concat(seq1)); it != end1; ++it)
		seq1Set.insert(abs(*it));

	// insert all blocks of seq2 that are in seq1Set to set
	TIterator end2 = end(concat(seq2));
	TSize pos = 1;
	for (TIterator it = begin(concat(seq2)); it != end2; ++it)
	{
		if (seq1Set.count(abs(*it)) != 0)
		{
			set[abs(*it)] = pos;
			++pos;
		}
	}
}

template<typename TBlockId, typename TSize>
void
commonBlocks(std::map<TBlockId, TSize> & set,
			 StringSet<String<TBlockId>, Dependent<> > const & seq1,
			 StringSet<String<TBlockId>, Dependent<> > const & seq2,
			 StringSet<String<TBlockId>, Dependent<> > const & seq3)
{
	typedef typename Iterator<ConcatenatorManyToOne<StringSet<String<TBlockId>, Dependent<> > > >::Type TIterator;

	std::set<TBlockId> seq1Set, seq2Set;

	// insert all blocks of seq1 to seq1Set
	TIterator end1 = end(concat(seq1));
	for (TIterator it = begin(concat(seq1)); it != end1; ++it)
		seq1Set.insert(abs(*it));

	// insert all blocks of seq2 that are in seq1Set to set
	TIterator end2 = end(concat(seq2));
	TSize pos = 1;
	for (TIterator it = begin(concat(seq2)); it != end2; ++it)
	{
		if (seq1Set.count(abs(*it)) != 0)
		{
			set[abs(*it)] = pos;
			++pos;
		}
		else seq2Set.insert(abs(*it));
	}

	// insert all blocks of seq3 that are in seq1Set or seq2Set to set
	TIterator end3 = end(concat(seq3));
	for (TIterator it = begin(concat(seq3)); it != end3; ++it)
	{
		if (set.count(abs(*it)) == 0 && (seq1Set.count(abs(*it)) != 0 || seq2Set.count(abs(*it)) != 0))
		{
			set[abs(*it)] = pos;
			++pos;
		}
	}
}

template<typename TBlockId>
void
uniqueBlocks(std::set<TBlockId> & set, String<TBlockId> const & seq1, String<TBlockId> const & seq2)
{
	typedef typename Iterator<String<TBlockId> >::Type TIterator;

	// insert all blocks of seq1 to set
	TIterator end1 = end(seq1);
	for (TIterator it1 = begin(seq1); it1 != end1; ++it1)
		set.insert(abs(*it1));

	// erase all blocks of seq2 if present from set, insert otherewise
	TIterator end2 = end(seq2);
	for (TIterator it2 = begin(seq2); it2 != end2; ++it2)
	{
		if (set.count(abs(*it2)) == 0) set.insert(abs(*it2));
		else set.erase(abs(*it2));
	}
}

template<typename TGraph, typename TSize>
void
createMatchingGraph(TGraph & graph, String<BlockInMatchingGraph> & nodes, TSize numBlocks)
{
	resize(nodes, numBlocks);

	for (TSize i = 0; i < numBlocks; ++i)
	{
		// add nodes for one block
		value(nodes, i).head = graph.addNode();
		value(nodes, i).tail = graph.addNode();
		value(nodes, i).headTelomere = graph.addNode();
		value(nodes, i).tailTelomere = graph.addNode();

		// add telomere edges for block
		graph.addEdge(value(nodes, i).head, value(nodes, i).headTelomere);
		graph.addEdge(value(nodes, i).tail, value(nodes, i).tailTelomere);
	}
}

template<typename TGraph, typename TEdgeMap, typename TBlockId, typename TSize>
TSize
addSequenceToMatchingGraph(TGraph & graph,
						   TEdgeMap & eMap,
						   String<TBlockId> const & seq,
						   String<BlockInMatchingGraph> const & nMap,
						   std::map<TBlockId, TSize> & blocks)
{
	typedef typename Iterator<String<TBlockId> >::Type TIterator;
	typedef typename TGraph::Node TNode;

	if (blocks.size() == 0) return 0;

	TSize weight = 0;
	TIterator it = begin(seq);
	TIterator itEnd = end(seq);

	TNode u, v, uTelo, vTelo;

	while (it != itEnd && blocks.count(abs(*it)) == 0) ++it;
	if (it == itEnd) return 0;

	// add head telomere weight
	if (*it > 0)
	{
		++eMap[lemon::findEdge(graph, value(nMap, blocks[abs(*it)]-1).head, value(nMap, blocks[abs(*it)]-1).headTelomere)];
		++weight;

		u = value(nMap, blocks[abs(*it)]-1).tail;
		uTelo = value(nMap, blocks[abs(*it)]-1).tailTelomere;
	}
	else
	{
		++eMap[lemon::findEdge(graph, value(nMap, blocks[abs(*it)]-1).tail, value(nMap, blocks[abs(*it)]-1).tailTelomere)];
		++weight;

		u = value(nMap, blocks[abs(*it)]-1).head;
		uTelo = value(nMap, blocks[abs(*it)]-1).headTelomere;
	}

	// add adjacency weights
	for (++it; it != itEnd; ++it)
	{
		if (blocks.count(abs(*it)) == 0) continue;

		if (*it > 0) v = value(nMap, blocks[abs(*it)]-1).head;
		else v = value(nMap, blocks[abs(*it)]-1).tail;

		typename TGraph::Edge edge = lemon::findEdge(graph, u, v);
		if (edge == lemon::Invalid())
		{
			edge = graph.addEdge(u, v);

			if (*it > 0) vTelo = value(nMap, blocks[abs(*it)]-1).headTelomere;
			else vTelo = value(nMap, blocks[abs(*it)]-1).tailTelomere;
			graph.addEdge(uTelo, vTelo);
		}
		++eMap[edge];
		++weight;

		if (*it > 0)
		{
			u = value(nMap, blocks[abs(*it)]-1).tail;
			uTelo = value(nMap, blocks[abs(*it)]-1).tailTelomere;
		}
		else
		{
			u = value(nMap, blocks[abs(*it)]-1).head;
			uTelo = value(nMap, blocks[abs(*it)]-1).headTelomere;
		}
	}

	// add tail telomere weight
	--it;
	while (blocks.count(abs(*it)) == 0) --it;
	if (*it > 0) vTelo = value(nMap, blocks[abs(*it)]-1).tailTelomere;
	else vTelo = value(nMap, blocks[abs(*it)]-1).headTelomere;
	++eMap[lemon::findEdge(graph, u, vTelo)];
	++weight;
	
	return weight;
}

template<typename TGraph, typename TEdgeMap, typename TBlockId, typename TSize>
TSize
addSequenceToMatchingGraph(TGraph & graph,
						   TEdgeMap & eMap,
						   StringSet<String<TBlockId>, Dependent<> > const & seq,
						   String<BlockInMatchingGraph> & nMap,
						   std::map<TBlockId, TSize> & blocks)
{
	typedef typename Iterator<StringSet<String<TBlockId>, Dependent<> > const>::Type TIterator;

	TSize weight = 0;
	for (TIterator it = begin(seq); it != end(seq); ++it)
		weight += addSequenceToMatchingGraph(graph, eMap, *it, nMap, blocks);

	return weight;
}

template<typename TBlockId, typename TSize>
typename Size<String<TBlockId> >::Type
pwCount(String<TBlockId> const & seq1, String<TBlockId> const & seq2, std::map<TBlockId, TSize> & blocks)
{
	typedef typename Position<String<TBlockId> >::Type TPosition;
	typedef Pair<TPosition, bool> TPair;
	typename Size<String<TBlockId> >::Type count = 0;

	// store positions of seq1 blocks in map
	std::map<TBlockId, TPair > seq1Map;
	TPosition pos = 0;
	for (TPosition i = 0; i < length(seq1); ++i)
	{
		if (blocks.count(abs(seq1[i])) == 0) continue;
		if (seq1[i] >= 0)
			seq1Map[seq1[i]] = TPair(pos, true);
		else seq1Map[-seq1[i]] = TPair(pos, false);
		++pos;
	}

	// sort seq2 according to seq1
	String<TPair> seq2Sorted;
	for (TPosition i = 0; i < length(seq2); ++i)
	{
		if (blocks.count(abs(seq2[i])) == 0) continue;

		bool orientation = true;
		if ((seq1Map[abs(seq2[i])].i2 && seq2[i] < 0) || (!seq1Map[abs(seq2[i])].i2 && seq2[i] > 0))
			orientation = false;

		appendValue(seq2Sorted, TPair(seq1Map[abs(seq2[i])].i1, orientation));
	}

	for (TPosition i = 1; i < length(seq2Sorted); ++i)
	{
		if (seq2Sorted[i-1].i2)
		{
			if (seq2Sorted[i].i2)
			{
				if (seq2Sorted[i-1].i1 + 1 != seq2Sorted[i].i1) ++count;
			}
			else
			{
				++count;
			}
		}
		else
		{
			if (!seq2Sorted[i].i2)
			{
				if (seq2Sorted[i-1].i1 != seq2Sorted[i].i1 + 1) ++count;
			}
			else ++count;
		}
	}

	return count;
}

template<typename TSeqId, typename TBlockId>
typename Size<String<TBlockId> >::Type
pairwiseCounts(std::map<TSeqId, StringSet<String<TBlockId>, Dependent<> > > & blockSeqs, bool detailed)
{
	typedef typename Size<String<TBlockId> >::Type TSize;
	typedef typename std::map<TSeqId, StringSet<String<TBlockId>, Dependent<> > >::const_iterator TSeqIterator;
	typedef typename lemon::SmartGraph TGraph;
	typedef typename TGraph::EdgeMap<int> TEdgeMap;

	if (blockSeqs.size() < 2) return 0;
	if (detailed)
		std::cout << "# 2-way counts: seq1, seq2, count\n";

	TSize pairwiseBreakpoints = 0;

	// iterate over all pairs of sequences
	TSeqIterator endSeqs = --blockSeqs.end();
	for (TSeqIterator itSeq1 = blockSeqs.begin(); itSeq1 != endSeqs; ++itSeq1)
	{
		++endSeqs;

		TSeqIterator itSeq2 = itSeq1;
		for (++itSeq2; itSeq2 != endSeqs; ++itSeq2)
		{
			std::map<TBlockId, TBlockId> blocks;
			commonBlocks(blocks, itSeq1->second, itSeq2->second);

			TGraph graph;
			TEdgeMap eMap(graph);
			String<BlockInMatchingGraph> nMap;
			TSize graphWeight = 0;

			// setup graph
			createMatchingGraph(graph, nMap, length(blocks));

			// set edge weights
			graphWeight += addSequenceToMatchingGraph(graph, eMap, itSeq1->second, nMap, blocks);
			graphWeight += addSequenceToMatchingGraph(graph, eMap, itSeq2->second, nMap, blocks);

			lemon::MaxWeightedPerfectMatching<TGraph, TEdgeMap> matching(graph, eMap);
			if (!matching.run())
			{
				std::cerr << "ERROR: Pairwise perfect Matching failed!";
			}

			SEQAN_ASSERT_LEQ(matching.matchingWeight(), static_cast<__int64>(graphWeight));
			TSize weightRemoved = graphWeight - matching.matchingWeight();
			if (detailed)
			{
				std::cout << itSeq1->first << "\t" << itSeq2->first << "\t" << weightRemoved;
				std::cout << std::endl;
			}

			pairwiseBreakpoints += weightRemoved;
		}
		--endSeqs;
	}

	return pairwiseBreakpoints;
}

template<typename TBlockId, typename TSize>
TSize
pairwiseWeightRemoved(StringSet<String<TBlockId>, Dependent<> > const & seq1,
					  StringSet<String<TBlockId>, Dependent<> > const & seq2,
					  std::map<TBlockId, TSize> & blocks)
{
	typedef typename lemon::SmartGraph TGraph;
	typedef TGraph::EdgeMap<int> TEdgeMap;

	TGraph graph;
	TEdgeMap eMap(graph);
	String<BlockInMatchingGraph> nMap;
	TSize graphWeight = 0;

	// setup graph
	createMatchingGraph(graph, nMap, length(blocks));

	// set edge weights
	graphWeight += addSequenceToMatchingGraph(graph, eMap, seq1, nMap, blocks);
	graphWeight += addSequenceToMatchingGraph(graph, eMap, seq2, nMap, blocks);

	lemon::MaxWeightedPerfectMatching<TGraph, TEdgeMap> matching(graph, eMap);
	if (!matching.run())
	{
		std::cerr << "ERROR: Perfect matching failed while computing triplet breakpoint counts!\n";
		return maxValue<TSize>();
	}

	SEQAN_ASSERT_LEQ(matching.matchingWeight(), static_cast<__int64>(graphWeight));
	TSize weightRemoved = graphWeight - matching.matchingWeight();

	return weightRemoved;
}

template<typename TBlockId, typename TSize>
TSize
tripletWeightRemoved(StringSet<String<TBlockId>, Dependent<> > const & seq1,
					 StringSet<String<TBlockId>, Dependent<> > const & seq2,
					 StringSet<String<TBlockId>, Dependent<> > const & seq3,
					 std::map<TBlockId, TSize> & blocks)
{
	typedef typename lemon::SmartGraph TGraph;
	typedef TGraph::EdgeMap<int> TEdgeMap;

	TGraph graph;
	TEdgeMap eMap(graph);
	String<BlockInMatchingGraph> nMap;
	TSize graphWeight = 0;

	// setup graph
	createMatchingGraph(graph, nMap, length(blocks));

	// set edge weights
	graphWeight += addSequenceToMatchingGraph(graph, eMap, seq1, nMap, blocks);
	graphWeight += addSequenceToMatchingGraph(graph, eMap, seq2, nMap, blocks);
	graphWeight += addSequenceToMatchingGraph(graph, eMap, seq3, nMap, blocks);

	lemon::MaxWeightedPerfectMatching<TGraph, TEdgeMap> matching(graph, eMap);
	if (!matching.run())
	{
		std::cerr << "ERROR: Perfect matching failed while computing triplet breakpoint counts!\n";
		return maxValue<TSize>();
	}

	SEQAN_ASSERT_LEQ(matching.matchingWeight(), static_cast<__int64>(graphWeight));
	TSize weightRemoved = graphWeight - matching.matchingWeight();

	return weightRemoved;
}

template<typename TSeqId, typename TBlockId>
double
tripletCounts(std::map<TSeqId, StringSet<String<TBlockId>, Dependent<> > > & blockSeqs, bool detailed)
{
	typedef typename std::map<TSeqId, StringSet<String<TBlockId>, Dependent<> > >::const_iterator TSeqIterator;
	typedef typename Size<String<TBlockId> >::Type TSize;

	if (blockSeqs.size() < 3) return 0;	
	if (detailed)
		std::cout << "# 3-way counts: seq1, seq2, seq3, count\n";

	double tripletCounts = 0;

	// iterate over all triplets of sequences
	TSeqIterator endSeqs = ----blockSeqs.end();
	for (TSeqIterator itSeq1 = blockSeqs.begin(); itSeq1 != endSeqs; ++itSeq1)
	{
		++endSeqs;
		TSeqIterator itSeq2 = itSeq1;
		for (++itSeq2; itSeq2 != endSeqs; ++itSeq2)
		{
			++endSeqs;
			TSeqIterator itSeq3 = itSeq2;
			for (++itSeq3; itSeq3 != endSeqs; ++itSeq3)
			{
				std::map<TBlockId, TSize> blocks;
				commonBlocks(blocks, itSeq1->second, itSeq2->second, itSeq3->second);
				
				TSize count12 = pairwiseWeightRemoved(itSeq1->second, itSeq2->second, blocks);
				TSize count13 = pairwiseWeightRemoved(itSeq1->second, itSeq3->second, blocks);
				TSize count23 = pairwiseWeightRemoved(itSeq2->second, itSeq3->second, blocks);
				TSize count123 = tripletWeightRemoved(itSeq1->second, itSeq2->second, itSeq3->second, blocks);

				SEQAN_ASSERT_GEQ(count123, (count12 + count13 + count23) / 2);

				double count = count123 - (count12 + count13 + count23) / 2.0;
				if (detailed)
					std::cout << itSeq1->first << "\t" << itSeq2->first << "\t" << itSeq3->first << "\t" << count << "\n";
				tripletCounts += count;
			}
			--endSeqs;
		}
		--endSeqs;
	}
	
	return tripletCounts;
}

#endif  // #ifndef SANDBOX_BKEHR_APPS_BREAKPOINT_CALCULATOR_BREAKPOINT_COUNTS_H_
