#ifndef SEQAN_HEADER_GRAPH_ALIGN_TCOFFEE_DISTANCE_H
#define SEQAN_HEADER_GRAPH_ALIGN_TCOFFEE_DISTANCE_H

namespace SEQAN_NAMESPACE_MAIN
{



//////////////////////////////////////////////////////////////////////////////
// T-Coffee - Determining sequence similarity
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////


template<typename TAlignmentGraph, typename TAlphabet>
inline double 
getSequenceSimilarity(TAlignmentGraph& g,
					  TAlphabet)
{
	SEQAN_CHECKPOINT
	typedef typename Size<TAlignmentGraph>::Type TSize;
	typedef typename Host<TAlignmentGraph>::Type TStringSet;
	typedef typename Value<TStringSet>::Type TString;
	typedef typename Infix<TString>::Type TInfix;

	TSize len1 = length(stringSet(g)[0]);
	TSize len2 = length(stringSet(g)[1]);

	double sim = 0.0;
	typedef typename Iterator<TAlignmentGraph, EdgeIterator>::Type TEdgeIterator;
	TEdgeIterator it(g);
	for(;!atEnd(it);++it) {
		typedef typename Iterator<TInfix>::Type TInfixIter;
		TInfix inf1 = label(g,sourceVertex(it));
		TInfix inf2 = label(g,targetVertex(it));
		TInfixIter sIt1 = begin(inf1);
		TInfixIter sIt2 = begin(inf2);
		while((!atEnd(sIt1)) || (!atEnd(sIt2))) {
			if ( (TAlphabet) *sIt1  == (TAlphabet) *sIt2) ++sim;
			goNext(sIt1); goNext(sIt2);
		}
	}
	
	// Normalize by distance
	if (len1 > len2) return sim / len2;
	else return sim / len1;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlignmentGraph, typename TValue, typename TSpec, typename TAlphabet>
inline void 
getSequenceSimilarity(TAlignmentGraph& g,
					  String<TValue, TSpec>& sim,
					  TAlphabet)
{
	SEQAN_CHECKPOINT
	typedef String<TValue, TSpec> TSimilarityMatrix;
	typedef typename Size<TAlignmentGraph>::Type TSize;
	typedef typename VertexDescriptor<TAlignmentGraph>::Type TVertexDescriptor;
	typedef typename Host<TAlignmentGraph>::Type TStringSet;
	typedef typename Value<TStringSet>::Type TString;
	typedef typename Infix<TString>::Type TInfix;

	TStringSet& str = stringSet(g);
	TSize numSeqs = length(str);
	clear(sim);
	fill(sim, numSeqs * numSeqs, 0);

	typedef typename Iterator<TAlignmentGraph, EdgeIterator>::Type TEdgeIterator;
	typedef typename Iterator<TInfix>::Type TInfixIter;
	TEdgeIterator it(g);
	for(;!atEnd(it);++it) {
		TVertexDescriptor sV = sourceVertex(it);
		TVertexDescriptor tV = targetVertex(it);
		TSize pos1 = idToPosition(str, sequenceId(g, sV));
		TSize pos2 = idToPosition(str, sequenceId(g, tV));
		TInfix inf1 = label(g,sV);
		TInfix inf2 = label(g,tV);
		TInfixIter sIt1 = begin(inf1);
		TInfixIter sIt2 = begin(inf2);
		while((!atEnd(sIt1)) || (!atEnd(sIt2))) {
			if ( (TAlphabet) *sIt1  == (TAlphabet) *sIt2) {
				value(sim, pos1 * numSeqs + pos2) += 1;
				value(sim, pos2 * numSeqs + pos1) += 1;
			}
			goNext(sIt1); goNext(sIt2);
		}
	}
	for(TSize i=0; i<numSeqs; ++i) {
		for(TSize j=i+1; j<numSeqs; ++j) {
			TSize len1 = length(str[i]);
			TSize len2 = length(str[j]);

			if (len1 > len2) {
				value(sim, i * numSeqs + j) /= len2;
				value(sim, j * numSeqs + i) /= len2;
				//value(sim, i * numSeqs + j) /= (len2 * len2);
				//value(sim, j * numSeqs + i) /= (len2 * len2);
			} else {
				value(sim, i * numSeqs + j) /= len1;
				value(sim, j * numSeqs + i) /= len1;
				//value(sim, i * numSeqs + j) /= (len1 * len1);
				//value(sim, j * numSeqs + i) /= (len1 * len1);
			}
		}
	}
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlignmentGraph, typename TValue, typename TSpec>
inline void 
getSequenceSimilarity(TAlignmentGraph& g,
					  String<TValue, TSpec>& sim)
{
	SEQAN_CHECKPOINT
	typedef typename Host<TAlignmentGraph>::Type TStringSet;
	typedef typename Value<TStringSet>::Type TString;
	getSequenceSimilarity(g, sim, typename Value<TString>::Type());
}

//////////////////////////////////////////////////////////////////////////////

template<typename TFragmentMatches, typename TStringSet, typename TSize1,  typename TSize2, typename TSize3, typename TAlphabet>
inline double 
getSequenceSimilarity(TFragmentMatches& matches,
					  TStringSet& str,
					  TSize1 const from,
					  TSize2 const to,
					  TSize3& alignLength,
					  TAlphabet)
{
	SEQAN_CHECKPOINT
	typedef typename Size<TFragmentMatches>::Type TSize;
	typedef typename Id<TFragmentMatches>::Type TId;
	typedef typename Value<TStringSet>::Type TString;
	typedef typename Infix<TString>::Type TInfix;

	TSize sim = 0;
	TSize match_length = 0;
	TSize len1 = length(str[0]);
	TSize len2 = length(str[1]);

	typedef typename Iterator<TInfix>::Type TInfixIter;
	for(TSize i = from;i<to;++i) {
		TInfix inf1 = label(matches[i], str, sequenceId(matches[i], 0));
		TInfix inf2 = label(matches[i], str, sequenceId(matches[i], 1));
		TInfixIter sIt1 = begin(inf1);
		TInfixIter sIt2 = begin(inf2);
		while((!atEnd(sIt1)) || (!atEnd(sIt2))) {
			if ( (TAlphabet) *sIt1  == (TAlphabet) *sIt2) {
				++sim;
			}
			goNext(sIt1); goNext(sIt2);
			++match_length;
		}
	}

	alignLength = match_length + (len1 - match_length) + (len2 - match_length);
	if (len1 > len2) return (double) sim / (double) len2;
	else return (double) sim / (double) len1;
}


//////////////////////////////////////////////////////////////////////////////
// T-Coffee - Similarity to distance matrix
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template<typename TMatrix>
inline void
similarityToDistanceMatrix(TMatrix& mat, 
						   KimuraDistance) 
{
	SEQAN_CHECKPOINT
	
	typedef typename Value<TMatrix>::Type TValue;
	typedef typename Size<TMatrix>::Type TSize;

	// Initialize the mat matrix
	TSize nseq = (TSize) sqrt((double)length(mat));

	// Calculate the mat
	for (TSize row=0;row<nseq;++row) {
		for(TSize col=row+1;col<nseq;++col) {
			// Kimura correction
			TValue val = 0.8 * getValue(mat, row*nseq+col) + 0.2;
			//if (val < 0.2) val = 0.2;
			val = (TValue) log((TValue) val - (1.0 - val * val) / 5.0);
			
			// Assign values
			assignValue(mat, row*nseq+col, val);
			assignValue(mat, col*nseq+row, val);
		}
		assignValue(mat, row*nseq+row, 0);
	}

	//// Debug code
	//for (TSize row=0;row<nseq;++row) {
	//	for(TSize col=0;col<nseq;++col) {
	//		std::cout << getValue(mat, row*nseq+col) << ",";
	//	}
	//	std::cout << std::endl;
	//}
	//std::cout << std::endl;
}


//////////////////////////////////////////////////////////////////////////////

template<typename TMatrix>
inline void
similarityToDistanceMatrix(TMatrix& mat, 
						   FractionalDistance) 
{
	SEQAN_CHECKPOINT
	
	typedef typename Value<TMatrix>::Type TValue;
	typedef typename Size<TMatrix>::Type TSize;

	// Initialize the mat matrix
	TSize nseq = (TSize) sqrt((double)length(mat));

	// Calculate the mat
	for (TSize row=0;row<nseq;++row) {
		for(TSize col=row+1;col<nseq;++col) {
			// Distance must be between 0 and 1
			TValue val = (1 - getValue(mat, row*nseq+col)) * 100;

			// Assign values
			assignValue(mat, row*nseq+col, val);
			assignValue(mat, col*nseq+row, val);
		}
		assignValue(mat, row*nseq+row, 0);
	}

	//// Debug code
	//for (TSize row=0;row<nseq;++row) {
	//	for(TSize col=0;col<nseq;++col) {
	//		std::cout << getValue(mat, row*nseq+col) << ",";
	//	}
	//	std::cout << std::endl;
	//}
	//std::cout << std::endl;
}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
