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
  $Id: graph_consensus_base.h 2103 2008-05-23 07:57:13Z rausch@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

#ifndef SEQAN_HEADER_CONSENSUS_BASE_H
#define SEQAN_HEADER_CONSENSUS_BASE_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////

/**
.Tag.File Format.tag.CeleraFrg message file:
	Celera fragment message file.
*/
struct TagCeleraFrg_;
typedef Tag<TagCeleraFrg_> const CeleraFrg;

/**
.Tag.File Format.tag.CeleraCgb message file:
	Celera cgb file.
*/
struct TagCeleraCgb_;
typedef Tag<TagCeleraCgb_> const CeleraCgb;



//////////////////////////////////////////////////////////////////////////////
// Segment Match Generation tag
//////////////////////////////////////////////////////////////////////////////

/**
.Tag.Segment Match Generation.value.Overlap_Library:
	Segment matches from overlap alignments.
*/

struct Overlap_Library_;
typedef Tag<Overlap_Library_> const Overlap_Library;



//////////////////////////////////////////////////////////////////////////////
// Consensus tag
//////////////////////////////////////////////////////////////////////////////

/**
.Tag.Consensus Calling:
..summary:A tag that specifies how to call the consensus.
*/


/**
.Tag.Consensus Calling.value.Majority_Vote:
	A consensus based on the most common character.
*/

struct Majority_Vote_;
typedef Tag<Majority_Vote_> const Majority_Vote;

/**
.Tag.Consensus Calling.value.Bayesian:
	A consensus based on bayesian probability.
*/

struct Bayesian_;
typedef Tag<Bayesian_> const Bayesian;



//////////////////////////////////////////////////////////////////////////////
// Read alignment and Consensus Generation
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////

struct ConsensusOptions {
public:
	// Bandwidth of overlap alignment
	int bandwidth;

	// Number of computed overlaps per read (at the beginning and end of a read)
	int overlaps;

	// Minimum match length of a computed overlap
	int matchlength;

	// Minimum quality (in percent identity) of a computed overlap
	int quality;

	// Window size, only relevant for insert sequencing
	// If window == 0, no insert sequencing is assumed
	int window;
	
	// SNP calling
	// 0: majority
	// 1: bayesian
	int snp;

	// Output
	// 0: seqan style
	// 1: afg output format
	int output;

	// Conversion option
	// 0: No conversion, regular consensus computation
	// 1: Creates an afg file
	// 2: Creates a Celera frg file
	// 3: Creates a Celera cgb file
	int convert;

	// Offset all reads, so the first read starts at position 0
	bool moveToFront;

	// Scoring object for overlap alignments
	Score<int> sc;

	// Various input and output files
	std::string readsfile;				// File of reads in FASTA format
	std::string afgfile;				// AMOS afg file input
	std::string source;					// Reference genome
	std::string outfile;				// Output file name
	
	// Initialization
	ConsensusOptions() : bandwidth(8), overlaps(3), matchlength(15), quality(80), window(0), snp(0), output(0), convert(0), moveToFront(false), outfile("readAlign.txt") 
	{
		sc = Score<int>(2,-6,-4,-9);
	}
};


//////////////////////////////////////////////////////////////////////////////////


template<typename TStringSet, typename TCargo, typename TSpec, typename TSize, typename TConfigOptions>
inline void
consensusAlignment(Graph<Alignment<TStringSet, TCargo, TSpec> >& gOut,
				   String<Pair<TSize, TSize> >& begEndPos,
				   TConfigOptions const& consOpt) 
{
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TOutGraph;
	typedef typename Id<TOutGraph>::Type TId;

	// Initialization
	TStringSet& seqSet = stringSet(gOut);
	TSize nseq = length(seqSet);

	// Select all overlapping reads and record the diagonals of the band
	String<Pair<TId, TId> > pList;
	String<Pair<int, int> > diagList;
	if (consOpt.window == 0) selectPairs(seqSet, begEndPos, consOpt.bandwidth, pList, diagList);
	else selectPairsIndel(seqSet, begEndPos, consOpt.window, pList, diagList);

	// Estimate the number of overlaps we want to compute
#ifdef SEQAN_PROFILE
	if (consOpt.window == 0) ::std::cout << "Matchlength: " << consOpt.matchlength << ", " << "Quality: " << consOpt.quality << ", " << "Bandwidth: " << consOpt.bandwidth << ", " << "Overlaps: " << consOpt.overlaps << ::std::endl;
	else ::std::cout << "Matchlength: " << consOpt.matchlength << ", " << "Quality: " << consOpt.quality << ", " << "Window: " << consOpt.window << ", " << "Overlaps: " << consOpt.overlaps << ::std::endl;
	std::cout << "Number of reads: " << nseq << std::endl;
	if (consOpt.window == 0) {
		TSize covEstim = length(pList) / nseq;
		std::cout << "Estimated coverage: " << covEstim << std::endl;
	}
	std::cout << "Pair selection done: " << SEQAN_PROTIMEUPDATE(__myProfileTime) << " seconds" << std::endl;
#endif

	// Set-up a sparse distance matrix
	Graph<Undirected<double> > pairGraph;
	
	// Containers for segment matches and corresponding scores 
	typedef String<Fragment<> > TFragmentString;
	TFragmentString matches;
	typedef String<int> TScoreValues;
	TScoreValues scores;

	// Compute segment matches from global pairwise alignments
	appendSegmentMatches(seqSet, pList, diagList, begEndPos, consOpt.sc, consOpt.matchlength, consOpt.quality, consOpt.overlaps, matches, scores, pairGraph, Overlap_Library() );
	clear(pList);
	clear(diagList);
#ifdef SEQAN_PROFILE
	std::cout << "Overlap done: " << SEQAN_PROTIMEUPDATE(__myProfileTime) << " seconds" << std::endl;
#endif

	// Re-Score the matches
	scoreMatches(seqSet, consOpt.sc, matches, scores);
#ifdef SEQAN_PROFILE
	std::cout << "Re-scoring done: " << SEQAN_PROTIMEUPDATE(__myProfileTime) << " seconds" << std::endl;
#endif

	// Use these segment matches for the initial alignment graph
	typedef Graph<Alignment<TStringSet, TSize> > TGraph;
	TGraph g(seqSet);
	buildAlignmentGraph(matches, scores, g, FractionalScore() );
	clear(matches);
	clear(scores);
#ifdef SEQAN_PROFILE
	std::cout << "Construction of Alignment Graph done: " << SEQAN_PROTIMEUPDATE(__myProfileTime) << " seconds" << std::endl;
#endif

	// Guide Tree
	Graph<Tree<double> > guideTree;
	upgmaTree(pairGraph, guideTree);
#ifdef SEQAN_PROFILE
	std::cout << "Guide tree done: " << SEQAN_PROTIMEUPDATE(__myProfileTime) << " seconds" << std::endl;
#endif
	clear(pairGraph);

	// Triplet library extension
	if ( ((2 * numEdges(g)) / numVertices(g) ) < 50 ) graphBasedTripletLibraryExtension(g);
	else reducedTripletLibraryExtension(g);
#ifdef SEQAN_PROFILE
	std::cout << "Triplet done: " << SEQAN_PROTIMEUPDATE(__myProfileTime) << " seconds" << std::endl;
#endif

	// Perform a progressive alignment
	progressiveAlignment(g, guideTree, gOut);
	clear(g);
	clear(guideTree);
#ifdef SEQAN_PROFILE
	std::cout << "Progressive alignment done: " << SEQAN_PROTIMEUPDATE(__myProfileTime) << " seconds" << std::endl;
#endif
}

//////////////////////////////////////////////////////////////////////////////////


template<typename TStringSet, typename TCargo, typename TSpec, typename TSize>
inline void
consensusAlignment(Graph<Alignment<TStringSet, TCargo, TSpec> >& gOut,
				   String<Pair<TSize, TSize> >& begEndPos) 
{
	ConsensusOptions consOpt;
	consensusAlignment(gOut, begEndPos, consOpt);
}


//////////////////////////////////////////////////////////////////////////////

template <typename TStringSet, typename TCargo, typename TSpec, typename TValue, typename TBegEndRowPos, typename TSize>
inline void
multireadAlignment(Graph<Alignment<TStringSet, TCargo, TSpec> > const& g,
				   String<TValue>& mat,
				   TBegEndRowPos& readBegEndRowPos,
				   TSize& maxCoverage)
{
	SEQAN_CHECKPOINT
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename Iterator<TGraph, VertexIterator>::Type TVertexIterator;
	typedef typename Id<TGraph>::Type TIdType;
	typedef typename Value<TStringSet>::Type TString;
	typedef typename Value<TString>::Type TAlphabet;
	typedef typename Infix<TString>::Type TInfix;
	typedef typename Iterator<TInfix>::Type TInfixIter;
	typedef std::map<std::pair<TIdType, TIdType>, TVertexDescriptor> TPosToVertexMap;
	typedef std::map<TSize, TSize> TComponentLength;
	
	// Strongly Connected Components, topological sort, and length of each component
	String<TSize> component;
	String<TSize> order;
	TComponentLength compLength;
	if (convertAlignment(g, component, order, compLength)) {
		TSize numOfComponents = length(order);
		TStringSet& strSet = stringSet(g);
		TSize nseq = length(strSet);
		clear(readBegEndRowPos);
		resize(readBegEndRowPos, nseq);
		
		// Assign to each sequence the start and end (in terms of component ranks)
		typedef std::map<TSize, TSize> TComponentToRank;
		TComponentToRank compToRank;
		for(TSize compIndex = 0; compIndex < numOfComponents; ++compIndex) {
			compToRank.insert(std::make_pair(order[compIndex], compIndex));
		}
		typedef Pair<TSize, TSize> TRankPair;
		typedef String<TRankPair> TSequenceToRanks;
		TSequenceToRanks seqToRank;
		fill(seqToRank, nseq, TRankPair(0, 0));	
		TVertexIterator itVertex(g);
		for(;!atEnd(itVertex);++itVertex) {
			TVertexDescriptor vert = value(itVertex);
			TSize seq = idToPosition(strSet, sequenceId(g, vert));
			if (fragmentBegin(g, vert) == 0) {
				(value(seqToRank, seq)).i1 = (compToRank.find(getProperty(component,vert)))->second;
			}
			if (fragmentBegin(g, vert) + fragmentLength(g, vert) == length(value(strSet, seq))) {
				(value(seqToRank, seq)).i2 = (compToRank.find(getProperty(component,vert)))->second;
			}
		}
		compToRank.clear();

		// Assign the sequences to rows
		String<TSize> seqToRow;
		resize(seqToRow, nseq);
		maxCoverage = 0;
		typedef std::set<TSize> TLeftOver;
		TLeftOver leftOver;
		for(TSize i=0;i<nseq; ++i) {
			leftOver.insert(i);
		}
		while(!leftOver.empty()) {
			typedef std::set<std::pair<TSize, TSize> > TSeqToBegin;
			TSeqToBegin seqToBegin;
			for(typename TLeftOver::const_iterator pos = leftOver.begin(); pos != leftOver.end(); ++pos) {
				seqToBegin.insert(std::make_pair((seqToRank[*pos]).i1, *pos));
			}
			TSize endPos = 0;
			for(typename TSeqToBegin::const_iterator s = seqToBegin.begin(); s != seqToBegin.end();++s) {
				if (endPos <= (*s).first) {
					TSize currentSeq = (*s).second;
					seqToRow[currentSeq] = maxCoverage;
					endPos = (seqToRank[currentSeq]).i2 + 2;
					leftOver.erase(currentSeq);
				}	
			}
			++maxCoverage;
		}


		// Create the matrix
		TSize len = 0;
		String<TSize> compOffset;
		fill(compOffset, numOfComponents, 0);
		for(TSize compIndex = 0; compIndex < numOfComponents; ++compIndex) {
			compOffset[order[compIndex]] = len;
			len+=compLength[order[compIndex]];
		}
		TValue gapChar = gapValue<TValue>();
		TValue specialGap = '.';
		fill(mat, len * maxCoverage, gapChar);

		// Fill in the segments
		for(typename TPosToVertexMap::const_iterator it = g.data_pvMap.begin();it != g.data_pvMap.end(); ++it) {
			TInfix str = label(g,it->second);
			TSize c = property(component, it->second);
			TSize strPos = idToPosition(strSet, it->first.first);
			TSize row = seqToRow[strPos];
			//if (row == 0) {
			//	std::cout << sequenceId(g, it->second) << ':' << str << ',' << strSet[sequenceId(g, it->second)] << std::endl;
			//	std::cout << getProperty(component, it->second) << ',' << order[compIndex] << std::endl;
			//	std::cout << (seqToRank[sequenceId(g, it->second)]).i1 << ',' << (seqToRank[sequenceId(g, it->second)]).i2 << std::endl;
			//}
			TInfixIter sIt = begin(str);
			TInfixIter sItEnd = end(str);
			TSize i = compOffset[c];
			for(TSize pCol = i;sIt!=sItEnd;goNext(sIt), ++pCol, ++i) assignValue(mat, row * len + pCol, *sIt);
		}
		String<bool> active;
		for(TSize compIndex = 0; compIndex < numOfComponents; ++compIndex) {
			TSize offset = compOffset[order[compIndex]];
			TSize currentCompLength = compLength[order[compIndex]];

			clear(active);
			fill(active, maxCoverage, false);

			// Find the empty rows
			for(TSize i=0;i<nseq; ++i) {
				if (((seqToRank[i]).i1 <= compIndex) && ((seqToRank[i]).i2 >= compIndex)) active[(seqToRow[i])] = true;
			}
			
			// Substitute false gaps with special gap character
			for(TSize i = 0; i < maxCoverage; ++i) {
				if (!(active[i])) {
					for(TSize pCol = offset;pCol < offset + currentCompLength;++pCol) assignValue(mat, i * len + pCol, specialGap);
				}
			}

		}

		// Get the new begin and end positions
		for(TSize i=0;i<nseq; ++i) {
			TVertexDescriptor lastVertex = findVertex(const_cast<TGraph&>(g), positionToId(strSet, i), length(strSet[i]) - 1);
			TSize readBegin = compOffset[getProperty(component, findVertex(const_cast<TGraph&>(g), positionToId(strSet, i), 0))];
			TSize readEnd = compOffset[getProperty(component, lastVertex)] + fragmentLength(const_cast<TGraph&>(g), lastVertex);
			readBegEndRowPos[i].i1 = readBegin;
			readBegEndRowPos[i].i2 = readEnd;
			readBegEndRowPos[i].i3 = seqToRow[i];
		}
	}
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSpec, typename TCounters, typename TCoverage, typename TSize, typename TAlphabet>
inline void
__countLetters(String<TValue, TSpec> const& mat,
			   TCounters& counterValues,
			   TCoverage& coverage,
			   TSize alignDepth,
			   TAlphabet)
{
	SEQAN_CHECKPOINT
	typedef String<TValue, TSpec> TMatrix;
	typedef typename Iterator<TMatrix>::Type TMatIter;
	typedef typename Iterator<TCoverage>::Type TCovIter;

	// Initialization
	TSize len = length(mat) / alignDepth;
	TValue gapChar = gapValue<TValue>();
	TValue specialGap = '.';
	TSize alphabetSize = ValueSize<TAlphabet>::VALUE;
	clear(coverage);
	fill(coverage, len, 0);


	// Set-up counter values
	typedef typename Value<TCounters>::Type TCounter;
	typedef typename Iterator<TCounters>::Type TCounterIt;
	clear(counterValues);
	resize(counterValues, len);
	for(TSize i=0;i<len; ++i) {
		TCounter counter;
		fill(counter, alphabetSize + 1, 0);
		value(counterValues, i) = counter;
	}

	// Count all 
	TMatIter matIt = begin(mat);
	TMatIter matItEnd = end(mat);
	TCounterIt countIt = begin(counterValues);
	TCovIter covIt = begin(coverage);
	TSize pos = 0;
	for(; matIt != matItEnd; goNext(matIt), goNext(countIt), goNext(covIt), ++pos) {
		if (pos % len == 0) {
			countIt = begin(counterValues);
			covIt = begin(coverage);
		}
		TValue c = value(matIt);
		if (c == specialGap) continue;
		else {
			++value(covIt);
			if (c == gapChar) ++value(value(countIt), alphabetSize);
			else ++value(value(countIt), ordValue(TAlphabet(c)));
		}
	}
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSpec, typename TAlphabet, typename TSpec2, typename TGappedConsensus, typename TCoverage, typename TSize>
inline void
consensusCalling(String<TValue, TSpec> const& mat,
				 String<TAlphabet, TSpec2>& consensus,
				 TGappedConsensus& gappedConsensus,
				 TCoverage& coverage,
				 TSize alignDepth,
				 Bayesian)
{
	SEQAN_CHECKPOINT
	typedef double TProbability;
	typedef String<TProbability> TProbabilityDistribution;
	typedef String<TProbabilityDistribution> TPositionalPrDist;
	typedef typename Iterator<TPositionalPrDist>::Type TPosPrDistIter;

	// Initialization
	TSize len = length(mat) / alignDepth;
	TValue gapChar = gapValue<TValue>();
	TValue specialGap = '.';
	TSize alphabetSize = ValueSize<TAlphabet>::VALUE;
	TProbabilityDistribution backroundDist;
	fill(backroundDist, alphabetSize + 1, (1.0 / (TProbability) (alphabetSize + 1)));
	
	// Set-up the counters
	typedef String<TSize> TCounter;
	typedef String<TCounter> TCounters;
	TCounters counterValues;
	__countLetters(mat, counterValues, coverage, alignDepth, TAlphabet() );

	// Get an initial consensus
	typedef typename Iterator<TCounters>::Type TCounterIt;
	TCounterIt countIt = begin(counterValues);
	TCounterIt countItEnd = end(counterValues);
	TPositionalPrDist posPrDist;
	for(;countIt != countItEnd; goNext(countIt)) {
		TSize max = 0;
		TValue c = TValue();
		typedef typename Iterator<TCounter>::Type TCIt;
		TCIt cIt = begin(value(countIt));
		TCIt cItEnd = end(value(countIt));
		TSize pos = 0;
		for(;cIt != cItEnd; goNext(cIt), ++pos) {
			if (value(cIt) > max) {
				max = value(cIt);
				if (pos == alphabetSize) c = gapChar;
				else c = TAlphabet(pos);
			}
		}
		TProbabilityDistribution prDist;
		fill(prDist, alphabetSize + 1, 0);
		if (c == gapChar) value(prDist, alphabetSize) = 1;
		else value(prDist, ordValue(TAlphabet(c))) = 1;
		appendValue(posPrDist, prDist);
	}

	bool run = false;
	TProbabilityDistribution pI;
	TProbabilityDistribution pIJ;
	TProbabilityDistribution pIOld;
	TProbabilityDistribution pIJOld;
	std::cout << "Bayesian Consensus";
	while((run) || (empty(pIOld))) {
		// Store the values from the last iteration
		pIOld = pI;
		pIJOld = pIJ;

		// Count all letters in the consensus
		TProbabilityDistribution nI;
		fill(nI, alphabetSize + 1, 0);
		TPosPrDistIter itPosPrDist = begin(posPrDist);
		TPosPrDistIter itPosPrDistEnd = end(posPrDist);
		for(;itPosPrDist!=itPosPrDistEnd; ++itPosPrDist) {
			for(TSize i = 0; i<(alphabetSize + 1); ++i) {
				value(nI, i) += value(value(itPosPrDist), i);
			}
		}
	
		// Composition probabilities
		clear(pI);
		resize(pI, alphabetSize + 1);
		for(TSize i = 0; i<length(pI); ++i) {
			value(pI, i) = (TProbability) value(nI, i) / (TProbability) length(posPrDist);
		}

		// Count all letters that agree / disagree with the consensus
		TProbabilityDistribution nIJ;
		fill(nIJ, (alphabetSize + 1) * (alphabetSize + 1), 0);
		typedef String<TValue, TSpec> TMatrix;
		typedef typename Iterator<TMatrix>::Type TMatIter;
		TMatIter matIt = begin(mat);
		TMatIter matItEnd = end(mat);
		itPosPrDist = begin(posPrDist);
		TSize pos = 0;
		for(; matIt != matItEnd; goNext(matIt), goNext(itPosPrDist), ++pos) {
			if (pos % len == 0) {
				itPosPrDist = begin(posPrDist);
			}
			TValue c = value(matIt);
			if (c == specialGap) continue;
			else {
				TSize fragJ = alphabetSize;
				if (c != gapChar) fragJ = ordValue(TAlphabet(c));
				for(TSize consI = 0; consI<(alphabetSize + 1); ++consI) {
					value(nIJ, consI * (alphabetSize + 1) + fragJ) += 1.0 * value(value(itPosPrDist), consI);
				}
			}
		}

		// Sequencing error probabilities
		clear(pIJ);
		resize(pIJ, (alphabetSize + 1) * (alphabetSize + 1));
		TProbability sumIJ = 0;
		for(TSize diag = 0; diag<(alphabetSize + 1); ++diag) sumIJ += value(nIJ, diag * (alphabetSize + 1) + diag);
		for(TSize consI = 0; consI<(alphabetSize + 1); ++consI) {
			for(TSize fragJ = 0; fragJ<(alphabetSize + 1); ++fragJ) {
				value(pIJ, consI * (alphabetSize + 1) + fragJ) = value(nIJ, consI * (alphabetSize + 1) + fragJ) / sumIJ;
			}
		}
	
		//// Debug Code
		//std::cout << "A " << value(pI, 0) << std::endl;
		//std::cout << "C " << value(pI, 1) << std::endl;
		//std::cout << "G " << value(pI, 2) << std::endl;
		//std::cout << "T " << value(pI, 3) << std::endl;
		//std::cout << "- " << value(pI, 4) << std::endl;
		//std::cout << "AA " << value(pIJ, 0) << std::endl;
		//std::cout << "AC " << value(pIJ, 1) << std::endl;
		//std::cout << "AG " << value(pIJ, 2) << std::endl;
		//std::cout << "AT " << value(pIJ, 3) << std::endl;
		//std::cout << "A- " << value(pIJ, 4) << std::endl;
		//std::cout << "CA " << value(pIJ, 5) << std::endl;
		//std::cout << "CC " << value(pIJ, 6) << std::endl;
		//std::cout << "CG " << value(pIJ, 7) << std::endl;
		//std::cout << "CT " << value(pIJ, 8) << std::endl;
		//std::cout << "C- " << value(pIJ, 9) << std::endl;
		//std::cout << "GA " << value(pIJ, 10) << std::endl;
		//std::cout << "GC " << value(pIJ, 11) << std::endl;
		//std::cout << "GG " << value(pIJ, 12) << std::endl;
		//std::cout << "GT " << value(pIJ, 13) << std::endl;
		//std::cout << "G- " << value(pIJ, 14) << std::endl;
		//std::cout << "TA " << value(pIJ, 15) << std::endl;
		//std::cout << "TC " << value(pIJ, 16) << std::endl;
		//std::cout << "TG " << value(pIJ, 17) << std::endl;
		//std::cout << "TT " << value(pIJ, 18) << std::endl;
		//std::cout << "T- " << value(pIJ, 19) << std::endl;
		//std::cout << "-A " << value(pIJ, 20) << std::endl;
		//std::cout << "-C " << value(pIJ, 21) << std::endl;
		//std::cout << "-G " << value(pIJ, 22) << std::endl;
		//std::cout << "-T " << value(pIJ, 23) << std::endl;
		//std::cout << "-- " << value(pIJ, 24) << std::endl;

		// Recompute positional probability distribution
		itPosPrDist = begin(posPrDist);
		TSize col = 0;
		for(;itPosPrDist!=itPosPrDistEnd; goNext(itPosPrDist), ++col) {
			TProbabilityDistribution prDist;
			resize(prDist, alphabetSize + 1);
			for(TSize consI = 0; consI<(alphabetSize + 1); ++consI) {
				TProbability numerator = value(pI, consI);
				TProbability denominator = 0;
				for(TSize allI = 0; allI<(alphabetSize + 1); ++allI) {
					TProbability denominatorSub = value(pI, allI);
					for(TSize row = 0; row < alignDepth; ++row) {
						TValue c = value(mat, row * len + col);
						if (c == specialGap) continue;
						TSize fragJ = alphabetSize;
						if (c != gapChar) fragJ = ordValue(TAlphabet(c));
						if (allI == consI) {
							numerator *= value(pIJ, allI * (alphabetSize + 1) + fragJ); 
						}
						denominatorSub *= value(pIJ, allI * (alphabetSize + 1) + fragJ); 
					}
					denominator += denominatorSub;
				}
				value(prDist, consI) = numerator / denominator;
			}
			value(itPosPrDist) = prDist;
		}	

		// Check termination criterion
		TProbability eps = 0.00001;
		typedef typename Iterator<TProbabilityDistribution>::Type TProbIter;
		TProbIter pIter = begin(pIOld);
		TProbIter pIterCompare = begin(pI);
		TProbIter pIterEnd = end(pIOld);
		run = false;
		for(;pIter != pIterEnd; goNext(pIter), goNext(pIterCompare)) {
			if (value(pIter) > value(pIterCompare)) {
				if (value(pIter) - value(pIterCompare) > eps) {
					run = true;
					break;
				}
			} else {
				if (value(pIterCompare) - value(pIter) > eps) {
					run = true;
					break;
				}
			}
		}
		if (!run) {
			pIter = begin(pIJOld);
			pIterCompare = begin(pIJ);
			pIterEnd = end(pIJOld);
			for(;pIter != pIterEnd; goNext(pIter), goNext(pIterCompare)) {
				if (value(pIter) > value(pIterCompare)) {
					if (value(pIter) - value(pIterCompare) > eps) {
						run = true;
						break;
					}
				} else {
					if (value(pIterCompare) - value(pIter) > eps) {
						run = true;
						break;
					}
				}
			}
		}
		std::cout << '.';
	}
	std::cout << std::endl;
	
	// Compute the most likely consensus
	TPosPrDistIter itPosPrDist = begin(posPrDist);
	TPosPrDistIter itPosPrDistEnd = end(posPrDist);
	clear(consensus);
	clear(gappedConsensus);
	for(;itPosPrDist!=itPosPrDistEnd; goNext(itPosPrDist)) {
		TProbability max = 0;
		TSize ind = 0;
		for(TSize consI = 0; consI<(alphabetSize + 1); ++consI) {
			if (value(value(itPosPrDist), consI) > max) {
				max = value(value(itPosPrDist), consI);
				ind = consI;
			}
		}
		if (ind == alphabetSize) appendValue(gappedConsensus, gapChar);
		else {
			appendValue(consensus, TAlphabet(ind));
			appendValue(gappedConsensus, TAlphabet(ind));
		}
	}
}


//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSpec, typename TStringSet, typename TScore, typename TBegEndRow, typename TSize>
inline TSize
fixDisruptedReads(String<TValue, TSpec>& mat,
				  TStringSet const& str,
				  TScore& scType,
				  TBegEndRow& begER,
				  TSize alignDepth)
{
	typedef typename Value< typename Value<TStringSet>::Type >::Type TAlphabet;

	// The limit expand threshold (the lower the more reads will be realigned)
	TSize limitExpand = 60;
	TValue gapChar = gapValue<TValue>();
	TValue specialGap = '.';
	TSize len = length(mat) / alignDepth;

	// Find all bad reads
	typedef typename Iterator<TBegEndRow>::Type TBEIter;
	typedef typename Iterator<TStringSet const>::Type TStrIter;
	TStrIter strIt = begin(str);
	TStrIter strItEnd = end(str);
	TBEIter itBE = begin(begER);
	TBEIter itBEEnd = end(begER);
	typedef String<TSize> TBadReads;
	TBadReads badReads;
	TSize counter = 0;
	for(;strIt != strItEnd; goNext(strIt), goNext(itBE), ++counter) {
		TSize b = value(itBE).i1;
		TSize e = value(itBE).i2;
		TSize diff = e - b;
		if (b > e) diff = b - e;
		if (length(value(strIt)) + limitExpand <= diff) {
			appendValue(badReads, counter);
		}
	}

	// Call the consensus
	String<unsigned int> coverage;
	String<char> gappedConsensus;
	String<Dna> consensusSequence;
	consensusCalling(mat, consensusSequence, gappedConsensus, coverage, alignDepth, Majority_Vote() );
	
	// Realign to the consensus
	typedef typename Iterator<TBadReads>::Type TBadIter;
	TBadIter badIt = begin(badReads);
	TBadIter badItEnd = end(badReads);
	for(;badIt != badItEnd; goNext(badIt)) {
		TSize b = value(begER, value(badIt)).i1;
		TSize e = value(begER, value(badIt)).i2;
		if (b > e) { TSize tmp = b; b = e; e = tmp; }
		String<char> cons = infix(gappedConsensus, b, e);
		String<TAlphabet> noGapCons;
		typedef typename Iterator<String<char> >::Type TCharIter;
		TCharIter charIt = begin(cons);
		TCharIter charItEnd = end(cons);
		for(;charIt != charItEnd; goNext(charIt)) {
			if (value(charIt) != gapChar) appendValue(noGapCons, value(charIt));
		}
		
		// Make a pairwise string-set
		TStringSet pairSet;
		//std::cout << noGapCons << std::endl;
		//std::cout << value(str, value(badIt)) << std::endl;
		appendValue(pairSet, noGapCons);
		appendValue(pairSet, value(str, value(badIt)));
		
		// Re-align
		// Maybe LCS is better???
		Graph<Alignment<TStringSet, TSize> > tmp(pairSet);
		globalAlignment(tmp, pairSet, scType, AlignConfig<true,true,true,true>(), Gotoh() );
		//std::cout << tmp << std::endl;		

		// Walk through all 3 sequences in parallel
		String<char> localAlign;
		convertAlignment(tmp, localAlign);
		String<char> consAlign = infix(localAlign, 0, length(localAlign) / 2);
		String<char> seqAlign = infix(localAlign, length(localAlign) / 2, length(localAlign));
		
		//std::cout << cons << std::endl;
		//std::cout << consAlign << std::endl;
		//std::cout << seqAlign << std::endl;
		TCharIter consIt = begin(cons);
		TCharIter consItEnd = end(cons);
		TCharIter consAlignIt = begin(consAlign);
		TCharIter seqAlignIt = begin(seqAlign);
		TSize row = value(begER, value(badIt)).i3;
		TSize col = b;
		while(consIt != consItEnd) {
			//std::cout << value(consIt) << ',' << value(consAlignIt) << ',' << value(seqAlignIt) << std::endl;
			if (value(consIt) == value(consAlignIt)) {
				// Both consensi have a gap
				value(mat, row * len + col) = value(seqAlignIt);
				goNext(consIt);
				goNext(consAlignIt);
				goNext(seqAlignIt);
				++col;
			} else {
				if (value(consIt) == gapChar) {
					while ((consIt != consItEnd) && (value(consIt) != value(consAlignIt))) {
						value(mat, row * len + col) = gapChar;
						goNext(consIt);
						++col;
					}
				} else {
					// This kind of read should be deleted
					while (value(consIt) != value(consAlignIt)) {
						goNext(consAlignIt);
					}
				}
			}
		}	
	}

	// Fix begin and end gaps
	itBE = begin(begER);
	itBEEnd = end(begER);
	for(;itBE!=itBEEnd; goNext(itBE)) {
		TSize b = value(itBE).i1;
		TSize e = value(itBE).i2;
		TSize row = value(itBE).i3;
		while (value(mat, row * len + b) == gapChar) {
			value(mat, row * len + b) = specialGap;
			++b;
		}
		while (value(mat, row * len + (e-1)) == gapChar) {
			value(mat, row * len + (e-1)) = specialGap;
			--e;
		}
		value(itBE).i1 = b;
		value(itBE).i2 = e;
	}


	return (length(badReads));
}


//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSpec, typename TAlphabet, typename TSpec2, typename TGappedConsensus, typename TCoverage, typename TSize>
inline void
consensusCalling(String<TValue, TSpec> const& mat,
				 String<TAlphabet, TSpec2>& consensus,
				 TGappedConsensus& gappedConsensus,
				 TCoverage& coverage,
				 TSize alignDepth,
				 Majority_Vote)
{
	// Initialization
	TSize alphabetSize = ValueSize<TAlphabet>::VALUE;
	TValue gapChar = gapValue<TValue>();

	// Set-up the counters
	typedef String<TSize> TCounter;
	typedef String<TCounter> TCounters;
	TCounters counterValues;
	__countLetters(mat, counterValues, coverage, alignDepth, TAlphabet() );

	// Get the consensus
	typedef typename Iterator<TCounters>::Type TCounterIt;
	TCounterIt countIt = begin(counterValues);
	TCounterIt countItEnd = end(counterValues);
	clear(consensus);
	clear(gappedConsensus);
	for(;countIt != countItEnd; goNext(countIt)) {
		TSize max = 0;
		TValue c = TValue();
		typedef typename Iterator<TCounter>::Type TCIt;
		TCIt cIt = begin(value(countIt));
		TCIt cItEnd = end(value(countIt));
		TSize pos = 0;
		for(;cIt != cItEnd; goNext(cIt), ++pos) {
			if (value(cIt) > max) {
				max = value(cIt);
				if (pos == alphabetSize) c = gapChar;
				else c = TAlphabet(pos);
			}
		}
		if (c != gapChar) appendValue(consensus, c);
		appendValue(gappedConsensus, c);
	}
}

//////////////////////////////////////////////////////////////////////////////

template <typename TStringSet, typename TCargo, typename TSpec, typename TBegEndRowPos, typename TNewLibraryGraph>
inline unsigned int
realignLowQualityReads(Graph<Alignment<TStringSet, TCargo, TSpec> > const& gIn,
					   TBegEndRowPos const& readBegEndRowPos,
					   TNewLibraryGraph& gOut)
{
	SEQAN_CHECKPOINT
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename Id<TGraph>::Type TId;
	typedef typename Size<TGraph>::Type TSize;
	typedef typename Iterator<TGraph, EdgeIterator>::Type TEdgeIterator;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename Iterator<TBegEndRowPos>::Type TBegEndIter;

	// Initialization
	TStringSet& str = stringSet(gIn);
	clearVertices(gOut);
	
	// Find disrupted reads
	std::set<TId> unalignedRead;
	TBegEndIter beIt = begin(readBegEndRowPos);
	TBegEndIter beItEnd = end(readBegEndRowPos);
	TSize pos = 0;
	for(;beIt != beItEnd; ++beIt, ++pos) {
		TSize lenStr = length(value(str, pos));
		if (((value(beIt)).i2 - (value(beIt)).i1) > (lenStr + lenStr / 3 + 10)) unalignedRead.insert(positionToId(str, pos));
	}

	// Any disrupted reads
	if (unalignedRead.empty()) return 0;
	else exit(0);

	
	//// String of fragments to combine all pairwise alignments into a multiple alignment
	//typedef Fragment<> TFragment;
	//typedef String<TFragment> TFragmentString;
	//typedef typename Iterator<TFragmentString>::Type TFragmentStringIter;
	//TFragmentString matches;
	//reserve(matches, numEdges(gIn));

	//// Insert all overlaps from the previous alignment
	//TEdgeIterator it_tmp(gIn);
	//for(;!atEnd(it_tmp);++it_tmp) appendValue(matches, TFragment(sequenceId(gIn, sourceVertex(it_tmp)),fragmentBegin(gIn, sourceVertex(it_tmp)), sequenceId(gIn, targetVertex(it_tmp)), fragmentBegin(gIn, targetVertex(it_tmp)), fragmentLength(gIn, sourceVertex(it_tmp))));
	//
	//// Recompute the overlap of interesting pairs
	//Score<int> score_type = Score<int>(2,-1,-4,-6);
	//TPairIter pairIt = begin(pList);
	//TPairIter pairItEnd = end(pList);
	//for(;pairIt != pairItEnd; ++pairIt) {
	//	TId id1 = (value(pairIt)).i1;
	//	TId id2 = (value(pairIt)).i2;

	//	if ((unalignedRead.find(id1) != unalignedRead.end()) || (unalignedRead.find(id2) != unalignedRead.end())) {
	//		// Make a pairwise string-set
	//		TStringSet pairSet;
	//		assignValueById(pairSet, str, id1);
	//		assignValueById(pairSet, str, id2);

	//		// Overlap alignment with a small mismatch score
	//		TSize from = length(matches);
	//		globalAlignment(matches, pairSet, score_type, AlignConfig<true,true,true,true>(), Gotoh() );

	//		// Determine a sequence weight
	//		TSize matchLen = 0;
	//		TSize overlapLen = 0;
	//		TSize alignLen = 0;
	//		getAlignmentStatistics(matches, pairSet, from, matchLen, overlapLen, alignLen);
	//		double quality = (double) matchLen / (double) overlapLen;

	//		// Take all overlaps of good quality
	//		if ((quality < 0.75) || (matchLen < 8)) {
	//			resize(matches, from);
	//		}
	//	}
	//}

	//// Refine all matches
	//matchRefinement(matches,stringSet(gOut),score_type, gOut);

	//return unalignedRead.size();
}






//////////////////////////////////////////////////////////////////////////////
// ToDo: Kill these functions, sooner or later :-)
//////////////////////////////////////////////////////////////////////////////

template<typename TSpec>
inline void 
reverseComplementInPlace(String<Dna5Q, TSpec>& str1)
{
	typedef typename Size<String<Dna5Q, TSpec> >::Type TSize;
	TSize pos1 = 0;
	TSize pos2 = length(str1)-1;
	for(;((pos1 < length(str1)) && (pos1<=pos2)); ++pos1, --pos2) {
		if ((str1[pos1] == 'A') && (str1[pos2] == 'A')) {
			Dna5Q c = 'T';
			assignQualityValue(c, getQualityValue(str1[pos2]));
			str1[pos1] = c;
			c = 'T';
			assignQualityValue(c, getQualityValue(str1[pos1]));
			str1[pos2] = c;
		} else if ((str1[pos1] == 'A') && (str1[pos2] == 'C')) {
			Dna5Q c = 'G';
			assignQualityValue(c, getQualityValue(str1[pos2]));
			str1[pos1] = c;
			c = 'T';
			assignQualityValue(c, getQualityValue(str1[pos1]));
			str1[pos2] = c;
		} else if ((str1[pos1] == 'A') && (str1[pos2] == 'G')) {
			Dna5Q c = 'C';
			assignQualityValue(c, getQualityValue(str1[pos2]));
			str1[pos1] = c;
			c = 'T';
			assignQualityValue(c, getQualityValue(str1[pos1]));
			str1[pos2] = c;
		} else if ((str1[pos1] == 'A') && (str1[pos2] == 'T')) {
			Dna5Q c = 'A';
			assignQualityValue(c, getQualityValue(str1[pos2]));
			str1[pos1] = c;
			c = 'T';
			assignQualityValue(c, getQualityValue(str1[pos1]));
			str1[pos2] = c;
		} else if ((str1[pos1] == 'A') && (str1[pos2] == 'N')) {
			Dna5Q c = 'N';
			assignQualityValue(c, getQualityValue(str1[pos2]));
			str1[pos1] = c;
			c = 'T';
			assignQualityValue(c, getQualityValue(str1[pos1]));
			str1[pos2] = c;
		} else if ((str1[pos1] == 'C') && (str1[pos2] == 'A')) {
			Dna5Q c = 'T';
			assignQualityValue(c, getQualityValue(str1[pos2]));
			str1[pos1] = c;
			c = 'G';
			assignQualityValue(c, getQualityValue(str1[pos1]));
			str1[pos2] = c;
		} else if ((str1[pos1] == 'C') && (str1[pos2] == 'C')) {
			Dna5Q c = 'G';
			assignQualityValue(c, getQualityValue(str1[pos2]));
			str1[pos1] = c;
			c = 'G';
			assignQualityValue(c, getQualityValue(str1[pos1]));
			str1[pos2] = c;
		} else if ((str1[pos1] == 'C') && (str1[pos2] == 'G')) {
			Dna5Q c = 'C';
			assignQualityValue(c, getQualityValue(str1[pos2]));
			str1[pos1] = c;
			c = 'G';
			assignQualityValue(c, getQualityValue(str1[pos1]));
			str1[pos2] = c;
		} else if ((str1[pos1] == 'C') && (str1[pos2] == 'T')) {
			Dna5Q c = 'A';
			assignQualityValue(c, getQualityValue(str1[pos2]));
			str1[pos1] = c;
			c = 'G';
			assignQualityValue(c, getQualityValue(str1[pos1]));
			str1[pos2] = c;
		} else if ((str1[pos1] == 'C') && (str1[pos2] == 'N')) {
			Dna5Q c = 'N';
			assignQualityValue(c, getQualityValue(str1[pos2]));
			str1[pos1] = c;
			c = 'G';
			assignQualityValue(c, getQualityValue(str1[pos1]));
			str1[pos2] = c;
		} else if ((str1[pos1] == 'G') && (str1[pos2] == 'A')) {
			Dna5Q c = 'T';
			assignQualityValue(c, getQualityValue(str1[pos2]));
			str1[pos1] = c;
			c = 'C';
			assignQualityValue(c, getQualityValue(str1[pos1]));
			str1[pos2] = c;
		} else if ((str1[pos1] == 'G') && (str1[pos2] == 'C')) {
			Dna5Q c = 'G';
			assignQualityValue(c, getQualityValue(str1[pos2]));
			str1[pos1] = c;
			c = 'C';
			assignQualityValue(c, getQualityValue(str1[pos1]));
			str1[pos2] = c;
		} else if ((str1[pos1] == 'G') && (str1[pos2] == 'G')) {
			Dna5Q c = 'C';
			assignQualityValue(c, getQualityValue(str1[pos2]));
			str1[pos1] = c;
			c = 'C';
			assignQualityValue(c, getQualityValue(str1[pos1]));
			str1[pos2] = c;
		} else if ((str1[pos1] == 'G') && (str1[pos2] == 'T')) {
			Dna5Q c = 'A';
			assignQualityValue(c, getQualityValue(str1[pos2]));
			str1[pos1] = c;
			c = 'C';
			assignQualityValue(c, getQualityValue(str1[pos1]));
			str1[pos2] = c;
		} else if ((str1[pos1] == 'G') && (str1[pos2] == 'N')) {
			Dna5Q c = 'N';
			assignQualityValue(c, getQualityValue(str1[pos2]));
			str1[pos1] = c;
			c = 'C';
			assignQualityValue(c, getQualityValue(str1[pos1]));
			str1[pos2] = c;
		}  else if ((str1[pos1] == 'T') && (str1[pos2] == 'A')) {
			Dna5Q c = 'T';
			assignQualityValue(c, getQualityValue(str1[pos2]));
			str1[pos1] = c;
			c = 'A';
			assignQualityValue(c, getQualityValue(str1[pos1]));
			str1[pos2] = c;
		} else if ((str1[pos1] == 'T') && (str1[pos2] == 'C')) {
			Dna5Q c = 'G';
			assignQualityValue(c, getQualityValue(str1[pos2]));
			str1[pos1] = c;
			c = 'A';
			assignQualityValue(c, getQualityValue(str1[pos1]));
			str1[pos2] = c;
		} else if ((str1[pos1] == 'T') && (str1[pos2] == 'G')) {
			Dna5Q c = 'C';
			assignQualityValue(c, getQualityValue(str1[pos2]));
			str1[pos1] = c;
			c = 'A';
			assignQualityValue(c, getQualityValue(str1[pos1]));
			str1[pos2] = c;
		} else if ((str1[pos1] == 'T') && (str1[pos2] == 'T')) {
			Dna5Q c = 'A';
			assignQualityValue(c, getQualityValue(str1[pos2]));
			str1[pos1] = c;
			c = 'A';
			assignQualityValue(c, getQualityValue(str1[pos1]));
			str1[pos2] = c;
		} else if ((str1[pos1] == 'T') && (str1[pos2] == 'N')) {
			Dna5Q c = 'N';
			assignQualityValue(c, getQualityValue(str1[pos2]));
			str1[pos1] = c;
			c = 'A';
			assignQualityValue(c, getQualityValue(str1[pos1]));
			str1[pos2] = c;
		} else if ((str1[pos1] == 'N') && (str1[pos2] == 'A')) {
			Dna5Q c = 'T';
			assignQualityValue(c, getQualityValue(str1[pos2]));
			str1[pos1] = c;
			c = 'N';
			assignQualityValue(c, getQualityValue(str1[pos1]));
			str1[pos2] = c;
		} else if ((str1[pos1] == 'N') && (str1[pos2] == 'C')) {
			Dna5Q c = 'G';
			assignQualityValue(c, getQualityValue(str1[pos2]));
			str1[pos1] = c;
			c = 'N';
			assignQualityValue(c, getQualityValue(str1[pos1]));
			str1[pos2] = c;
		} else if ((str1[pos1] == 'N') && (str1[pos2] == 'G')) {
			Dna5Q c = 'C';
			assignQualityValue(c, getQualityValue(str1[pos2]));
			str1[pos1] = c;
			c = 'N';
			assignQualityValue(c, getQualityValue(str1[pos1]));
			str1[pos2] = c;
		} else if ((str1[pos1] == 'N') && (str1[pos2] == 'T')) {
			Dna5Q c = 'A';
			assignQualityValue(c, getQualityValue(str1[pos2]));
			str1[pos1] = c;
			c = 'N';
			assignQualityValue(c, getQualityValue(str1[pos1]));
			str1[pos2] = c;
		} else if ((str1[pos1] == 'N') && (str1[pos2] == 'N')) {
			Dna5Q c = 'N';
			assignQualityValue(c, getQualityValue(str1[pos2]));
			str1[pos1] = c;
			c = 'N';
			assignQualityValue(c, getQualityValue(str1[pos1]));
			str1[pos2] = c;
		}
	}
}


template<typename TSpec>
inline void 
reverseComplementInPlace(String<Dna5Q, TSpec> const& str1)
{
	reverseComplementInPlace(const_cast<String<Dna5Q, TSpec>&>(str1));
}


/*
//////////////////////////////////////////////////////////////////////////////


template<typename TFile, typename TAlphabet, typename TSpec, typename TFragmentStore, typename TLibraryStore, typename TContigStore>
inline void 
write(TFile& target,
	  ReadStore<TAlphabet, TSpec>& readSt,
	  TFragmentStore&,
	  TLibraryStore&,
	  TContigStore&,
	  CeleraFrg) 
{
	SEQAN_CHECKPOINT
	typedef ReadStore<TAlphabet, TSpec> TReadStore;
	typedef typename Id<TReadStore>::Type TId;
	typedef typename Size<TReadStore>::Type TSize;
	typedef typename Value<TFile>::Type TValue;
	typedef String<char> TCharString;

	// Write Reads
	for(TSize i = 0; i<length(readSt); ++i) {
		_streamWrite(target,"{FRG\n");
		_streamWrite(target,"act:");
		_streamPut(target, 'A');
		_streamPut(target, '\n');
		_streamWrite(target,"acc:");
		_streamPutInt(target, i + 1);
		_streamPut(target, '\n');
		_streamWrite(target,"typ:");
		_streamPut(target, 'R');
		_streamPut(target, '\n');
		_streamWrite(target,"src:\n");
		_streamWrite(target, value(readSt.data_names, i));
		_streamWrite(target, "\n.\n");
		_streamWrite(target,"etm:");
		_streamPut(target, '0');
		_streamPut(target, '\n');
		_streamWrite(target,"seq:\n");
		String<TAlphabet> read;
		loadRead(readSt, i, read);
		for(TSize k = 0;k<length(read); k+=70) {
			TSize endK = k + 70;
			if (endK > length(read)) endK = length(read);
			_streamWrite(target, infix(read, k, endK));
			_streamPut(target, '\n');
		}
		_streamWrite(target, ".\n");
		_streamWrite(target,"qlt:\n");
		String<char> qlt;
		loadQuality(readSt, i, qlt);
		for(TSize k = 0;k<length(qlt); k+=70) {
			TSize endK = k + 70;
			if (endK > length(qlt)) endK = length(qlt);
			_streamWrite(target, infix(qlt, k, endK));
			_streamPut(target, '\n');
		}
		_streamWrite(target, ".\n");
		_streamWrite(target,"clr:");
		_streamPutInt(target, value(readSt.data_clr, i).i1);
		_streamPut(target, ',');
		_streamPutInt(target, value(readSt.data_clr, i).i2);
		_streamPut(target, '\n');
		_streamWrite(target,"}\n");
	}
}


//////////////////////////////////////////////////////////////////////////////


template<typename TFile, typename TAlphabet, typename TSpec, typename TFragmentStore, typename TLibraryStore, typename TContigStore>
inline void 
write(TFile& target,
	  ReadStore<TAlphabet, TSpec>& readSt,
	  TFragmentStore&,
	  TLibraryStore&,
	  TContigStore& ctgSt,
	  CeleraCgb) 
{
	SEQAN_CHECKPOINT
	typedef ReadStore<TAlphabet, TSpec> TReadStore;
	typedef typename Id<TReadStore>::Type TId;
	typedef typename Size<TReadStore>::Type TSize;
	typedef typename Value<TFile>::Type TValue;
	typedef String<char> TCharString;
	String<GappedRead<> >& origGappedReads = value(ctgSt.data_reads, 0);
	String<char> contig;
	loadContig(ctgSt, 0, contig);

	// Find smallest offset
	TSize offsetLeft = length(contig);
	for(TSize k = 0; k<length(origGappedReads); ++k) if ((value(origGappedReads, k)).data_offset < offsetLeft) offsetLeft = (value(origGappedReads, k)).data_offset;

	// Sort the reads
	typedef std::set<std::pair<TSize, TSize> > TOffsetIndexPair;
	TOffsetIndexPair offsetIndexSet;
	for(TSize k = 0; k<length(origGappedReads); ++k) {
		TSize clr1 = (value(origGappedReads, k)).data_clr.i1;
		TSize clr2 = (value(origGappedReads, k)).data_clr.i2;
		clr1 += ((value(origGappedReads, k)).data_offset - offsetLeft);
		clr2 += ((value(origGappedReads, k)).data_offset - offsetLeft);
		if (clr1 > clr2) { TSize tmp = clr1; clr1 = clr2; clr2 = tmp; }
		offsetIndexSet.insert(std::make_pair(clr1, k));
	}
	String<GappedRead<> > gappedReads;
	for(typename TOffsetIndexPair::const_iterator itPos = offsetIndexSet.begin(); itPos != offsetIndexSet.end(); ++itPos) {
		appendValue(gappedReads, value(origGappedReads, (*itPos).second));
	}

	//// Write IAF record for all reads
	//for(TSize k = 0; k<length(gappedReads); ++k) {
	//	_streamWrite(target,"{IAF\n");
	//	_streamWrite(target,"acc:");
	//	_streamPutInt(target, (value(gappedReads, k)).data_source + 1);
	//	_streamPut(target, '\n');
	//	_streamWrite(target,"typ:");
	//	_streamPut(target, 'R');
	//	_streamPut(target, '\n');
	//	_streamWrite(target,"chi:0\ncha:0\nclr:-1,-1\nmst:U\n}\n");
	//}

	// Write Header
	_streamWrite(target,"{IUM\nacc:0\nsrc:\ngen> @@ [0,0]\n.\ncov:0.000\nsta:X\nfur:X\nabp:0\nbbp:0\n");
	_streamWrite(target,"len:");
	_streamPutInt(target, length(contig));
	_streamPut(target, '\n');
	_streamWrite(target,"cns:\n.\nqlt:\n.\nfor:0\n");
	_streamWrite(target,"nfr:");
	_streamPutInt(target, length(readSt));
	_streamPut(target, '\n');
	
	// Write gapped reads
	for(TSize k = 0; k<length(gappedReads); ++k) {
		_streamWrite(target,"{IMP\n");
		_streamWrite(target,"typ:");
		_streamPut(target, 'R');
		_streamPut(target, '\n');
		_streamWrite(target,"mid:");
		_streamPutInt(target, (value(gappedReads, k)).data_source + 1);
		_streamPut(target, '\n');
		TSize clr1 = (value(gappedReads, k)).data_clr.i1;
		TSize clr2 = (value(gappedReads, k)).data_clr.i2;
		clr1 += ((value(gappedReads, k)).data_offset - offsetLeft);
		clr2 += ((value(gappedReads, k)).data_offset - offsetLeft);
		_streamWrite(target,"con:");
		//TSize orig1 = clr1;
		//TSize orig2 = clr2;
		//if (orig1 > orig2) {TSize tmp = orig1; orig1 = orig2; orig2 = tmp; }
		//TSize best = 0;
		//TSize bestDist = 0;
		//for(TSize other = 0; other<length(gappedReads); ++other) {
		//	if (other == k) continue;
		//	TSize thisRead1 = (value(gappedReads, other)).data_clr.i1;
		//	TSize thisRead2 = (value(gappedReads, other)).data_clr.i2;
		//	thisRead1 += ((value(gappedReads, other)).data_offset - offsetLeft);
		//	thisRead2 += ((value(gappedReads, other)).data_offset - offsetLeft);
		//	if (thisRead1 > thisRead2) {TSize tmp = thisRead1; thisRead1 = thisRead2; thisRead2 = tmp; }
		//	if ((orig1 > thisRead1) && (orig2 < thisRead2)) {
		//		if ((best == 0) ||
		//			(bestDist > (thisRead2 - thisRead1))) {
		//			bestDist = (thisRead2 - thisRead1);
		//			best = other + 1;
		//		}
		//	}
		//}

		//_streamPutInt(target, best);
		_streamPut(target, '0');
		_streamPut(target, '\n');
		_streamWrite(target,"pos:");
		_streamPutInt(target, clr1);
		_streamPut(target, ',');
		_streamPutInt(target, clr2);
		_streamPut(target, '\n');
		_streamWrite(target,"dln:0\n");
		_streamWrite(target,"del:\n");
		_streamWrite(target,"}\n");
	}
	_streamWrite(target,"}\n");
}

*/




//////////////////////////////////////////////////////////////////////////////

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
