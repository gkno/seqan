#ifndef SEQAN_HEADER_GRAPH_ALIGN_TCOFFEE_H
#define SEQAN_HEADER_GRAPH_ALIGN_TCOFFEE_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Graph: T-Coffee - Tags
//////////////////////////////////////////////////////////////////////////////

/**
.Tag.TCoffeeDistance
..summary:Tag to use the T-Coffee distance metric.
..value.TCoffeeDistance:Use the TCoffee Distance.
*/
struct TCoffeeDistance_;
typedef Tag<TCoffeeDistance_> const TCoffeeDistance;

/**
.Tag.FractionalDistance
..summary:Tag to use the fractional distance metric.
..value.FractionalDistance:Use the Fractional Distance.
*/
struct FractionalDistance_;
typedef Tag<FractionalDistance_> const FractionalDistance;

/**
.Tag.GlobalPairwise_Library
..summary:Tag to specify the type of the library.
..value.GlobalPairwise_Library:Use of a pairwise global library.
*/
struct GlobalPairwise_Library_;
typedef Tag<GlobalPairwise_Library_> const GlobalPairwise_Library;

/**
.Tag.MUM_Library
..summary:Tag to specify the type of the library.
..value.MUM_Library:Use of a maximal unique match library.
*/
struct MUM_Library_;
typedef Tag<MUM_Library_> const MUM_Library;


//////////////////////////////////////////////////////////////////////////////
// Graph: T-Coffee - Scoring Schema
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////

template<typename TGraphType>
struct ScoreAlignmentGraph;

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TGraphType>
class Score<TValue, ScoreAlignmentGraph<TGraphType> >
{
public:
	Holder<TGraphType> data_graph;

public:
	Score(TGraphType _graph) {
		data_graph = _graph;
	}

	Score(Score const & other) {
		data_graph = other.data_graph;
	}

	Score & operator = (Score const & other) {
		data_graph = other.data_graph;
		return *this;
	}
};

template <typename TValue, typename TGraphType>
inline TValue
scoreGapExtend(Score<TValue, ScoreAlignmentGraph<TGraphType> > & me)
{
	return 0;
}


template <typename TValue, typename TGraphType>
inline TValue const
scoreGapExtend(Score<TValue, ScoreAlignmentGraph<TGraphType> > const & me)
{
	return 0;
}

template <typename TValue, typename TGraphType>
inline TValue
scoreGapOpen(Score<TValue, ScoreAlignmentGraph<TGraphType> > & me)
{
	return 0;
}
template <typename TValue, typename TGraphType>
inline TValue const
scoreGapOpen(Score<TValue, ScoreAlignmentGraph<TGraphType> > const & me)
{
	return 0;
}

template <typename TValue, typename TGraphType, typename T>
inline TValue
score(Score<TValue, ScoreAlignmentGraph<TGraphType> > & me,
	  T const & left,
	  T const & right)
{
	typedef typename Size<T>::Type TSize;
	typedef typename VertexDescriptor<TGraphType>::Type TVertexDescriptor;
	typedef typename EdgeDescriptor<TGraphType>::Type TEdgeDescriptor;
	TSize len1 = length(left);
	TSize len2 = length(right);
	TSize divider = len1 * len2;
	TVertexDescriptor nilVertex = getNil<TVertexDescriptor>();
	TValue sum = 0;
	for(TSize i = 0;i<len1;++i) {
		for(TSize j = 0;j<len2;++j) {
			if ((left[i] != nilVertex) && (right[j] != nilVertex)) {
				TEdgeDescriptor e = findEdge(value(me.data_graph), left[i], right[j]);
				if (e != 0) sum += getCargo(e);
			}
		}
	}
	return sum / divider;
}






//////////////////////////////////////////////////////////////////////////////
// Graph: T-Coffee - Distance Matrix
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
template<typename TMatrix>
void
kmerToDistanceMatrix(TMatrix& mat, TCoffeeDistance) {
	SEQAN_CHECKPOINT
	
	typedef typename Value<TMatrix>::Type TValue;
	typedef typename Size<TMatrix>::Type TSize;

	// Initialize the mat matrix
	TSize nseq = length(mat, 0);

	// Calculate the mat
	for (TSize row=0;row<nseq;++row) {
		for(TSize col=row+1;col<nseq;++col) {
			TValue mat0 = min(getValue(mat, row*nseq+row),getValue(mat, col*nseq+col));
			assignValue(mat, row*nseq+col, ((mat0 - getValue(mat, row*nseq+col)) / mat0 * 3 * 10.0 + 0.5) );
			assignValue(mat, col*nseq+row, getValue(mat, row*nseq+col));
		}
		assignValue(mat, row*nseq+row, 0);
	}

	// Normalize matrix
	double limit = 1000;
	double max = 100000;
	double norm = 100;
	for (TSize row=0; row<nseq; ++row) {
		for (TSize col=row+1; col<nseq; ++col) {
			assignValue(mat, row*nseq+col, ((limit - (100 - getValue(mat,row*nseq+col))) * norm ) / max);
			assignValue(mat, col*nseq+row, getValue(mat, row*nseq+col));
		}
	}


	//for (TSize row=0;row<nseq;++row) {
	//	for(TSize col=0;col<nseq;++col) {
	//		std::cout << getValue(mat, row*nseq+col) << ",";
	//	}
	//	std::cout << std::endl;
	//}
}

//////////////////////////////////////////////////////////////////////////////

template<typename TMatrix>
void
kmerToDistanceMatrix(TMatrix& mat, FractionalDistance) {
	SEQAN_CHECKPOINT
	
	typedef typename Value<TMatrix>::Type TValue;
	typedef typename Size<TMatrix>::Type TSize;

	// Initialize the mat matrix
	TSize nseq = length(mat, 0);

	// Calculate the mat
	for (TSize row=0;row<nseq;++row) {
		for(TSize col=row+1;col<nseq;++col) {
			// First the fractional common kmer count
			TValue val = getValue(mat, row*nseq+col);
			TValue minVal = getValue(mat, row*nseq+row);
			TValue tmpVal;
			if ((tmpVal = getValue(mat, col*nseq+col)) < minVal) minVal = tmpVal;
			val /= minVal;
			
			// The transformed measure
			//val = log10(0.1 + val);

			// Kimura correction
			//TValue d_sq = 0.0;
			//for(TSize i=0;i<nseq;++i) {
			//	d_sq += (getValue(mat, row*nseq+i) * getValue(mat, i*nseq+col));
			//}
			//val = log(val - (1 - d_sq) / 5);

			// Assign the final value
			assignValue(mat, row*nseq+col, 1 - val);
			assignValue(mat, col*nseq+row, 1 - val);
		}
		assignValue(mat, row*nseq+row, 0);
	}

	//for (TSize row=0;row<nseq;++row) {
	//	for(TSize col=0;col<nseq;++col) {
	//		std::cout << getValue(mat, row*nseq+col) << ",";
	//	}
	//	std::cout << std::endl;
	//}
}

//////////////////////////////////////////////////////////////////////////////

template<typename TString, typename TSpec, typename THitMatrix, typename TSize, typename TAlphabet>
void
getCommonKmerMatrix(StringSet<TString, TSpec> const& strSet, THitMatrix& mat, TSize ktup, TAlphabet) {
	SEQAN_CHECKPOINT
	typedef unsigned int TWord;
	typedef StringSet<TString, TSpec> TStringSet;
	typedef typename Size<TStringSet>::Type TStringSetSize;
	typedef typename Size<THitMatrix>::Type TMatrixSize;

	// Number of sequences
	TStringSetSize nseq = length(strSet);

	// Initialization
	// Matrix for common k-tupels between sequence i and j
	setDimension(mat, 2);setLength(mat, 0, nseq);setLength(mat, 1, nseq);
	fill(host(mat), nseq*nseq, 0);

	// StringSet where each string is a sequence of amino acid groups identifiers
	typedef StringSet<String<TAlphabet>, Owner<> > TStringSetAA;
	// q-gram length = ktup
	typedef Index<TStringSetAA, Index_QGram<SimpleShape> > TIndex;
	TIndex index;
	resize(indexText(index), nseq);
	resize(indexShape(index), (unsigned) ktup);

	// Recode the strings into amino acid groups
	for(TStringSetSize k=0;k<length(strSet);++k) indexText(index)[k] = strSet[k];

	// Build index	
	indexCreate(index, QGram_SA());

	// Counter for each sequence
	String<TWord> counter;
	resize(counter, nseq);
	
	TWord nqgrams = length(indexDir(index)) - 1;
	if (nqgrams > 0) {
		TWord j1;
		TWord j2 = indexDir(index)[0];
		for(TWord i = 0; i < nqgrams; ++i) {  // Iterate over all q-grams
			j1 = j2;
			j2 = indexDir(index)[i+1];
			if (j1 == j2) continue;   //Empty hit-list for this q-gram

			// Clear the counters
			arrayFill(begin(counter, Standard()), end(counter, Standard()), 0);
			for(TWord j = j1; j < j2; ++j) ++counter[getValueI1(indexSA(index)[j])];

			// Add-up the values in the matrix
			for(TWord k = 0; k < nseq; ++k)
				for(TWord k2 = k; k2 < nseq; ++k2) {
					TWord minVal = counter[k];
					if (counter[k2] < minVal) minVal = counter[k2];
					assignValue(mat, k*nseq+k2, getValue(mat, k*nseq+k2) + minVal);
				}
		}
	}

	// Copy upper triangle to lower triangle
	for(TWord k = 0; k < nseq; ++k) {
		for(TWord k2 = k + 1; k2 < nseq; ++k2) {
			assignValue(mat, k2*nseq+k, getValue(mat, k*nseq+k2));
		}
	}

	//for (TWord row=0;row<nseq;++row) {
	//	for(TWord col=0;col<nseq;++col) {
	//		std::cout << getValue(mat, row*nseq+col) << ",";
	//	}
	//	std::cout << std::endl;
	//}
}


//////////////////////////////////////////////////////////////////////////////

template<typename TString, typename TSpec, typename THitMatrix>
void
getCommonKmerMatrix(StringSet<TString, TSpec> const& strSet, THitMatrix& mat) {
	SEQAN_CHECKPOINT
	getCommonKmerMatrix(strSet, mat, 3);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TString, typename TSpec, typename THitMatrix, typename TSize>
void
getCommonKmerMatrix(StringSet<TString, TSpec> const& strSet, THitMatrix& mat, TSize ktup) {
	SEQAN_CHECKPOINT
	
	typedef typename Value<TString>::Type TAlphabet;

	// Get the common q-grams / k-tupels
	getCommonKmerMatrix(strSet, mat, ktup, TAlphabet());
}

//////////////////////////////////////////////////////////////////////////////

template<typename TMatrixSpec, typename TCargo, typename TSpec>
void
slowNjTree(Matrix<double, TMatrixSpec>& mat, Graph<Tree<TCargo, TSpec> >& g) {
	SEQAN_CHECKPOINT
	
	typedef typename Size<Matrix<double, TMatrixSpec> >::Type TSize;
	typedef Graph<Tree<TCargo, TSpec> > TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	TVertexDescriptor nilVertex = getNil<TVertexDescriptor>();
	
	TSize nseq = length(mat, 0);

	//for(TSize i=0;i<nseq;++i) {
	//	for(TSize j=0;j<nseq;++j) {
	//		std::cout << getValue(mat, i*nseq+j) << ",";
	//	}
	//	std::cout << std::endl;
	//}

	// First initialization
	clearVertices(g);
	String<double> av;    // Average branch length to a combined node
	fill(av,nseq,0.0);

	String<TVertexDescriptor> connector;   // Nodes that need to be connected
	resize(connector, nseq);

	for(TSize i=0;i<nseq;++i) {
		addVertex(g);  // Add all the nodes that correspond to sequences
		assignValue(connector, i, i);
		assignValue(mat, i*nseq+i, 0.0);
	}

	// Main cycle
	double fnseqs=(double) nseq;
	for(TSize nc=0; nc<(nseq-3); ++nc) {
		double sumOfBranches = 0.0;

		// Copy upper triangle matrix to lower triangle
		for(TSize col=1; col<nseq; ++col) {
			for(TSize row=0; row<col; ++row) {
				assignValue(mat, col*nseq+row, getValue(mat, row*nseq+col));
				// Determine the sum of all branches
				sumOfBranches = sumOfBranches + getValue(mat, row*nseq+col);
			}
		}

		// Compute the sum of branch lengths for all possible pairs
		double tmin = 0.0;	
		TSize mini = 0;  // Next pair of seq i and j to join
		TSize minj = 0;
		for(TSize col=1; col<nseq; ++col)  {
			if (getValue(connector,col) != nilVertex) {
				for(TSize row=0; row<col; ++row) {
					if (getValue(connector,row) != nilVertex) {
						double diToAllOthers = 0.0;
						double djToAllOthers = 0.0;
						
						for(TSize i=0; i<nseq; ++i) {
							diToAllOthers += getValue(mat, i*nseq+row);
							djToAllOthers += getValue(mat, i*nseq+col);
						}

						double dij = getValue(mat, row*nseq+col);
						double total = diToAllOthers + djToAllOthers + (fnseqs - 2.0)*dij +2.0*(sumOfBranches - diToAllOthers - djToAllOthers);
						total /= (2.0*(fnseqs - 2.0));

						if ((tmin == 0) || (total < tmin)) {
							tmin = total;
							mini = row;
							minj = col;
						}
					}
				}
			}
		}

		// Print nodes that are about to be joined
		//std::cout << mini << std::endl;
		//std::cout << minj << std::endl;
		//std::cout << tmin << std::endl;
		//std::cout << std::endl;
		
		// Compute branch lengths
		double dMinIToOthers = 0.0;
		double dMinJToOthers = 0.0;
		for(TSize i=0; i<nseq; ++i) {
			dMinIToOthers += getValue(mat, i*nseq + mini);
			dMinJToOthers += getValue(mat, i*nseq + minj);
		}
		double dmin = getValue(mat, mini*nseq + minj);
		dMinIToOthers = dMinIToOthers / (fnseqs - 2.0);
		dMinJToOthers = dMinJToOthers / (fnseqs - 2.0);
		double iBranch = (dmin + dMinIToOthers - dMinJToOthers) * 0.5;
		double jBranch = dmin - iBranch;
		iBranch -= av[mini];
		jBranch -= av[minj];
		
		// Set negative branch length to zero
		if( fabs(iBranch) < 0.0001) iBranch = 0.0;
		if( fabs(jBranch) < 0.0001) jBranch = 0.0;
	
		// Print branch lengths
		//std::cout << iBranch << std::endl;
		//std::cout << jBranch << std::endl;
		//std::cout << std::endl;
		
		// Build tree
		TVertexDescriptor internalVertex = addVertex(g);
		addEdge(g, internalVertex, getValue(connector, mini), (TCargo) iBranch);
		addEdge(g, internalVertex, getValue(connector, minj), (TCargo) jBranch);

		// Remember the average branch length for the new combined node
		// Must be subtracted from all branches that include this node
		if(dmin <= 0.0) dmin = 0.000001;
		av[mini] = dmin * 0.5;


		// Re-initialisation
		// mini becomes the new combined node, minj is killed
		fnseqs = fnseqs - 1.0;
		assignValue(connector, minj, nilVertex);
		assignValue(connector, mini, internalVertex);

		for(TSize j=0; j<nseq; ++j) {
			if( getValue(connector, j) != nilVertex ) {
				double minIminJToOther = ( getValue(mat, mini*nseq+j) + getValue(mat, minj*nseq+j)) * 0.5;
				// Use upper triangle
				if((TSize) mini < j) assignValue(mat, mini*nseq+j, minIminJToOther);
				if((TSize) mini > j) assignValue(mat, j*nseq+mini, minIminJToOther);
			}
		}
		for(TSize j=0; j<nseq; ++j) {
			assignValue(mat, minj*nseq+j, 0.0);
			assignValue(mat, j*nseq+minj, 0.0);
		}
	}

	// Only three nodes left

	// Find the remaining nodes
	String<TSize> l;
	fill(l,3,0);
	TSize count = 0;
	for(TSize i=0; i<nseq; ++i) {
		if(getValue(connector, i) != nilVertex) {
			l[count] = i;
			++count;
		}
	}

	// Remaining nodes
	//std::cout << l[0] << std::endl;
	//std::cout << l[1] << std::endl;
	//std::cout << l[2] << std::endl;
	//std::cout << std::endl;

	String<double> branch;
	resize(branch, 3);
	branch[0] = (getValue(mat, l[0]*nseq+l[1]) + getValue(mat, l[0]*nseq+l[2]) - getValue(mat, l[1]*nseq+l[2])) * 0.5;
	branch[1] = (getValue(mat, l[1]*nseq+l[2]) + getValue(mat, l[0]*nseq+l[1]) - getValue(mat, l[0]*nseq+l[2])) * 0.5;
	branch[2] =  (getValue(mat, l[0]*nseq+l[2]) + getValue(mat, l[0]*nseq+l[1]) - getValue(mat, l[1]*nseq+l[2])) * 0.5;
    
	branch[0] -= av[l[0]];
	branch[1] -= av[l[1]];
	branch[2] -= av[l[2]];

	// Print branch lengths
	//std::cout << branch[0] << std::endl;
	//std::cout << branch[1] << std::endl;
	//std::cout << branch[2] << std::endl;
	//std::cout << std::endl;
    
	// Reset tiny negative and positive branch lengths to zero
	if( fabs(branch[0]) < 0.0001) branch[0] = 0.0;
	if( fabs(branch[1]) < 0.0001) branch[1] = 0.0;
	if( fabs(branch[2]) < 0.0001) branch[2] = 0.0;
    
	// Build tree
	TVertexDescriptor internalVertex = addVertex(g);
	addEdge(g, internalVertex, getValue(connector, l[0]), (TCargo) branch[0]);
	addEdge(g, internalVertex, getValue(connector, l[1]), (TCargo) branch[1]);
	addEdge(g, internalVertex, getValue(connector, l[2]), (TCargo) branch[2]);
	g.data_root = internalVertex;
}



//////////////////////////////////////////////////////////////////////////////
// T-Coffee: Primary Library Generation
//////////////////////////////////////////////////////////////////////////////

template<typename TAlignmentGraph, typename TAlphabet>
double 
_getSequenceSimilarity(TAlignmentGraph& g,
					   TAlphabet)
{
	SEQAN_CHECKPOINT
	typedef typename Size<TAlignmentGraph>::Type TSize;
	typedef typename Id<TAlignmentGraph>::Type TId;
	typedef typename VertexDescriptor<TAlignmentGraph>::Type TVertexDescriptor;
	typedef typename Host<TAlignmentGraph>::Type TStringSet;
	typedef typename Value<TStringSet>::Type TString;
	typedef typename Infix<TString>::Type TInfix;

	TSize len1 = length(stringSet(g)[0]);
	TSize len2 = length(stringSet(g)[1]);

	double sim = 0.0;
	typedef typename Iterator<TAlignmentGraph, EdgeIterator>::Type TEdgeIterator;
	TEdgeIterator it(g);
	for(;!atEnd(it);++it) {
		TVertexDescriptor sV = sourceVertex(it);
		TVertexDescriptor tV = targetVertex(it);
		typedef typename Iterator<TInfix>::Type TInfixIter;
		TInfix inf1 = label(g,sV);
		TInfix inf2 = label(g,tV);
		TInfixIter sIt1 = begin(inf1);
		TInfixIter sIt2 = begin(inf2);
		while((!atEnd(sIt1)) || (!atEnd(sIt2))) {
			if ( (TAlphabet) *sIt1  == (TAlphabet) *sIt2) ++sim;
			goNext(sIt1); goNext(sIt2);
		}
	}
	
	if (len1 > len2) return sim / len2;
	else return sim / len1;
}

template<typename TStringSet, typename TCargo, typename TSpec, typename TSize>
void 
generatePrimaryLibrary(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
					   TSize minLen,
					   MUM_Library)
{
	SEQAN_CHECKPOINT

	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename Id<TGraph>::Type TId;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename Value<TStringSet>::Type TString;
	typedef Index<StringSet<TString, Owner<> >, Index_ESA<> > TIndex;
	
	clearVertices(g);

	TStringSet& str = stringSet(g);	
	TSize nseq = length(str);
	TIndex index;
	resize(indexText(index), nseq);
	for(TSize i = 0;i<nseq;++i) indexText(index)[i] = str[i];
	indexRequire(index, ESA_SA());
	indexRequire(index, ESA_LCP());
	indexRequire(index, ESA_BWT());

	typename Iterator<TIndex, MUMs>::Type it(index, minLen);	// set minimum MUM length
	String< typename SAValue<TIndex>::Type > occs;			// temp. string storing the hit positions

	while (!atEnd(it)) 
	{
		occs = getOccurences(it);							// gives hit positions (seqNo,seqOfs)
		orderOccurences(occs);								// order them by seqNo
			
		TSize matchLen = repLength(it);
		for(TSize i = 0; i < (TSize) length(occs); ++i) {
			//std::cout << positionToId(str, i) << ',' << getValueI2(occs[i]) << ',' << matchLen << std::endl;
			TVertexDescriptor v = addVertex(g, positionToId(str, i), getValueI2(occs[i]), matchLen);
			if (i > 0) {
				for(TSize k = i; k>0;--k) {
					addEdge(g,v,v-k,matchLen);
				}
			}
		}
		++it;
	}
}

template<typename TStringSet, typename TCargo, typename TSpec, typename TAlphabet>
void 
generatePrimaryLibrary(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
					   TAlphabet,
					   GlobalPairwise_Library)
{
	SEQAN_CHECKPOINT

	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename Size<TGraph>::Type TSize;
	typedef typename Id<TGraph>::Type TId;
	typedef typename Iterator<TGraph, EdgeIterator>::Type TEdgeIterator;
	typedef std::map<std::pair<TId,TId>, TCargo> TSeqSimilarity;

	clearVertices(g);

	// Pairwise alignments for all pairs of sequences
	Score<double> score_type = Score<double>(2,-1,-0.5,-2);
	TStringSet& str = stringSet(g);	
	TSize nseq = length(stringSet(g));
	TSeqSimilarity seqSimMap;

	// String of fragments to combine all pairwise alignments into a multiple alignment
	typedef Fragment<> TFragment;
	typedef String<TFragment, External<> > TFragmentString;
	TFragmentString matches;

	for(TSize i=0; i<nseq; ++i) {
		for(TSize j=i+1; j<nseq; ++j) {
			TId id1 = positionToId(str, i);
			TId id2 = positionToId(str, j);
			// Pairwise alignment graph
			TStringSet pairSet;
			assignValueById(pairSet, str, id1);
			assignValueById(pairSet, str, id2);
			typedef Graph<Alignment<TStringSet, void> > TPairGraph;
			typedef typename VertexDescriptor<TPairGraph>::Type TVD;
			typedef typename Iterator<TPairGraph, EdgeIterator>::Type TEI;
			TPairGraph pGraph(pairSet);
			globalAlignment(pGraph, score_type, Hirschberg() );

			// Determine a sequence weight
			TCargo seqSim = (TCargo) (_getSequenceSimilarity(pGraph, TAlphabet() ) * 100);
			if (id1 < id2) seqSimMap.insert(std::make_pair(std::make_pair(id1,id2),seqSim));
			else seqSimMap.insert(std::make_pair(std::make_pair(id2,id1),seqSim));
			//std::cout << pairSet[0] << std::endl;
			//std::cout << pairSet[1] << std::endl;
			//std::cout << pGraph << std::endl;
			//std::cout << seqSim << std::endl;
			
			TEI it(pGraph);
			for(;!atEnd(it);++it) {
				TVD sV = sourceVertex(it);
				TVD tV = targetVertex(it);
				push_back(matches, TFragment( (unsigned int) sequenceId(pGraph, sV), (unsigned int) fragmentBegin(pGraph,sV), (unsigned int) sequenceId(pGraph, tV),  (unsigned int)  fragmentBegin(pGraph,tV),  (unsigned int)  fragmentLength(pGraph,tV)));
			}
		}
	}

	//// Debug Code
	//// Print all the matches
	//std::cout << "The four sequences:" << std::endl;
	//for(TSize i = 0;i<length(str);++i) {
	//	std::cout << positionToId(str,i) << ':' << str[i] << std::endl;
	//}
	//std::cout << "The matches:" << std::endl;
	//for(TSize i = 0;i<length(matches);++i) {
	//	for(TSize i = 0;i<length(matches);++i) {
	//		TId tmp_id1 = sequenceId(matches[i],0);
	//		std::cout << tmp_id1 << ',' << fragmentBegin(matches[i],tmp_id1) << ',';
	//		for(TSize j = fragmentBegin(matches[i],tmp_id1); j < fragmentBegin(matches[i],tmp_id1) + fragmentLength(matches[i],tmp_id1); ++j) {
	//			std::cout << str[idToPosition(str, tmp_id1)][j];
	//		}
	//		TId tmp_id2 = sequenceId(matches[i],1);
	//		std::cout << ',' <<	tmp_id2 << ',' << fragmentBegin(matches[i],tmp_id2) << ',';
	//		for(TSize j = fragmentBegin(matches[i],tmp_id2); j < fragmentBegin(matches[i],tmp_id2) + fragmentLength(matches[i],tmp_id2); ++j) {
	//			std::cout << str[idToPosition(str, tmp_id2)][j];
	//		}
	//		std::cout << std::endl;
	//	}
	//}

	// Refine all matches and create multiple alignment
	matchRefinement(matches,stringSet(g),g);

	//// Debug Code
	//std::cout << "Refined matches" << std::endl;
	//TEdgeIterator it_tmp(g);
	//for(;!atEnd(it_tmp);++it_tmp) {
	//	TId id1 = sequenceId(g,sourceVertex(it_tmp));
	//	TId id2 = sequenceId(g,targetVertex(it_tmp));
	//	std::cout << id1 << ',' << fragmentBegin(g,sourceVertex(it_tmp)) << ',';
	//	std::cout << label(g,sourceVertex(it_tmp));
	//	std::cout << ',' <<	id2 << ',' << fragmentBegin(g,targetVertex(it_tmp)) << ',';
	//	std::cout << label(g,targetVertex(it_tmp));
	//	std::cout << std::endl;	
	//}

	// Adapt edge weights
	TEdgeIterator it(g);
	for(;!atEnd(it);++it) {
		TId id1 = sequenceId(g,sourceVertex(it));
		TId id2 = sequenceId(g,targetVertex(it));
		if (id1<id2) cargo(*it) = (seqSimMap.find(std::make_pair(id1,id2)))->second;
		else cargo(*it) = (seqSimMap.find(std::make_pair(id2,id1)))->second;
	}
}


template<typename TStringSet, typename TCargo, typename TSpec>
void 
tripletLibraryExtension(Graph<Alignment<TStringSet, TCargo, TSpec> >& g)
{
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename Size<TGraph>::Type TSize;
	typedef typename Id<TGraph>::Type TId;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename Iterator<TGraph, VertexIterator>::Type TVertexIterator;
	typedef typename Iterator<TGraph, EdgeIterator>::Type TEdgeIterator;
	typedef typename Iterator<TGraph, OutEdgeIterator>::Type TOutEdgeIterator;
	typedef std::map<TVertexDescriptor, TCargo> TIdToWeightMap;
	TIdToWeightMap idToWeightMap;

	// First add edges for the case that a and c is aligned, b and c is aligned, but a and b are not
	// Give these edges a weight of 0 because extension follows
	String<TVertexDescriptor> edges;
	TVertexIterator itVertex(g);
	for(;!atEnd(itVertex);++itVertex) {
		TOutEdgeIterator outIt1(g, *itVertex);
		if (atEnd(outIt1)) continue;
		TOutEdgeIterator outIt2(g, *itVertex);
		goNext(outIt2);
		while (!atEnd(outIt2)) {
			TVertexDescriptor tV1 = targetVertex(outIt1);
			TVertexDescriptor tV2 = targetVertex(outIt2);
			if ((sequenceId(g, tV1) != sequenceId(g,tV2)) && 
				(findEdge(g, tV1, tV2) == 0)) {
				appendValue(edges, tV1);
				appendValue(edges, tV2);
			}
			goNext(outIt1);
			goNext(outIt2);
		}
	}
	for(TSize i=0; i<length(edges);i+=2) addEdge(g, edges[i], edges[i+1], 0);
	
	// Now augment all existing edges
	String<TCargo> newCargoMap;
	resizeEdgeMap(g, newCargoMap);
	TEdgeIterator it(g);
	for(;!atEnd(it);++it) {
		assignProperty(newCargoMap, *it, cargo(*it));
		TVertexDescriptor sourceV = sourceVertex(it);
		TVertexDescriptor targetV = targetVertex(it);
		TOutEdgeIterator outIt1(g, sourceV);
		for(;!atEnd(outIt1);++outIt1) {
			if (targetV == targetVertex(outIt1)) continue;
			idToWeightMap.insert(std::make_pair(targetVertex(outIt1), getCargo(*outIt1)));
		}
		TOutEdgeIterator outIt2(g, targetV);
		for(;!atEnd(outIt2);++outIt2) {
			if (sourceV == targetVertex(outIt2)) continue;
			typename TIdToWeightMap::const_iterator pos = idToWeightMap.find(targetVertex(outIt2));
			if (pos != idToWeightMap.end()) {
				// Add the minimum of the two alignment values
				if (getCargo(*outIt2) > pos->second) property(newCargoMap, *it) += pos->second;
				else property(newCargoMap, *it) += getCargo(*outIt2);
			}
		}
		idToWeightMap.clear();
	}
	// Assign the new weights
	goBegin(it);
	for(;!atEnd(it);++it) cargo(*it) = getProperty(newCargoMap, *it);
}


template<typename TStringSet, typename TCargo, typename TSpec, typename TStringSet2, typename TSequence, typename TTag>
inline void 
_alignStringSetAccordingToGraph(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
								TStringSet2& strSet,
								TSequence& alignSeq,
								TTag)
{
	SEQAN_CHECKPOINT
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename Size<TGraph>::Type TSize;
	
	Score<typename Cargo<TGraph>::Type, ScoreAlignmentGraph<TGraph> > score_type = Score<typename Cargo<TGraph>::Type, ScoreAlignmentGraph<TGraph> >(g);
	TSequence tmp;
	globalAlignment(tmp, strSet, score_type, TTag());
	TSize len = length(tmp);
	clear(alignSeq);
	resize(alignSeq, len);
	for(TSize i = len; i>0;--i) alignSeq[len-i] = tmp[i-1];
}

template<typename TStringSet, typename TCargo, typename TSpec, typename TGuideTree, typename TVertexDescriptor, typename TSequence, typename TTag>
void 
_recursiveProgressiveAlignment(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
							   TGuideTree& tree,
							   TVertexDescriptor const root,
							   TSequence& alignSeq,
							   TTag tag)
{
	SEQAN_CHECKPOINT

	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename Size<TGraph>::Type TSize;
	typedef typename Id<TGraph>::Type TId;
	typedef typename Iterator<TGuideTree, AdjacencyIterator>::Type TAdjacencyIterator;
	typedef typename Iterator<TGuideTree, DfsPreorder>::Type TDfsPreorderIterator;

	TStringSet& str = stringSet(g);
	TVertexDescriptor nilVertex = getNil<TVertexDescriptor>();

	if(isLeaf(tree, root)) {
		TId seqId = positionToId(str, root);
		TSize i = 0;
		while(i<length(str[root])) {
			TVertexDescriptor nextVertex = findVertex(g, seqId, i);
			if (nextVertex == nilVertex) {
				TSize j = i + 1;
				while ((j < length(str[root])) && (findVertex(g, seqId, j) == nilVertex)) ++j;
				nextVertex = addVertex(g, seqId, i, j-i);
			}
			appendValue(alignSeq, String<TVertexDescriptor>(nextVertex));
			i += fragmentLength(g, nextVertex);
		}
	} else {
		// Align the two children
		TAdjacencyIterator adjIt(tree, root);
		typedef String<String<TVertexDescriptor> > TSegmentString;
		typedef StringSet<TSegmentString, Dependent<> > TSegmentStringSet;
		TSegmentString seq1;
		TSegmentString seq2;
		TSegmentStringSet strSet;
		bool first = true;
		for(;!atEnd(adjIt);goNext(adjIt)) {
			if (first) {
				_recursiveProgressiveAlignment(g,tree, *adjIt, seq1, TTag());
				assignValueById(strSet, seq1);
				first = false;
			} else {
				_recursiveProgressiveAlignment(g,tree, *adjIt, seq2, TTag());
				assignValueById(strSet, seq2);
				_alignStringSetAccordingToGraph(g,strSet,alignSeq, TTag());
			}
		}

		//// Debug Code
		//for(unsigned int i = 0; i<length(alignSeq);++i) {
		//	std::cout << '(';
		//	for(unsigned int j=0; j<length(alignSeq[i]);++j) {
		//		std::cout << getValue(alignSeq[i], j) << ',';
		//	}
		//	std::cout << ')' << ',';
		//}
		//std::cout << std::endl;
	}
}



template<typename TStringSet, typename TCargo, typename TSpec, typename TGuideTree, typename TOutGraph, typename TTag>
void 
progressiveAlignment(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
					 TGuideTree& tree,
					 TOutGraph& gOut,
					 TTag)
{
	SEQAN_CHECKPOINT

	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename Size<TGraph>::Type TSize;
	typedef typename Id<TGraph>::Type TId;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename Iterator<TGuideTree, AdjacencyIterator>::Type TAdjacencyIterator;
	typedef typename Iterator<TGuideTree, DfsPreorder>::Type TDfsPreorderIterator;

	// Align all three subtrees
	TAdjacencyIterator adjIt(tree, getRoot(tree));
	typedef String<String<TVertexDescriptor> > TSegmentString;
	typedef StringSet<TSegmentString, Dependent<> > TSegmentStringSet;
	TSegmentString seq1;
	TSegmentString seq2;
	TSegmentString alignSeq;
	TSegmentString seq3;
	TSegmentStringSet strSet;
	unsigned int count = 0;
	for(;!atEnd(adjIt);goNext(adjIt)) {
		if (count == 0) {
			_recursiveProgressiveAlignment(g,tree, *adjIt, seq1, TTag());
			assignValueById(strSet, seq1);
			++count;
		} else if (count == 1) {
			_recursiveProgressiveAlignment(g,tree, *adjIt, seq2, TTag());
			assignValueById(strSet, seq2);
			_alignStringSetAccordingToGraph(g,strSet,alignSeq, TTag());
			clear(strSet);
			assignValueById(strSet, alignSeq);
			++count;
		} else {
			_recursiveProgressiveAlignment(g,tree, *adjIt, seq3, TTag());
			assignValueById(strSet, seq3);
			_alignStringSetAccordingToGraph(g,strSet,alignSeq, TTag());
		}
	}

	// Create the alignment graph
	TVertexDescriptor nilVertex = getNil<TVertexDescriptor>();
	for(TSize i = 0; i<length(alignSeq);++i) {
		for(TSize j=0; j<length(alignSeq[i]);++j) {
			TVertexDescriptor v = getValue(alignSeq[i], j);
			if (v == nilVertex) continue;
			TVertexDescriptor l = addVertex(gOut, sequenceId(g, v), fragmentBegin(g,v), fragmentLength(g,v));
			TSize count = 1;
			for(TSize k = j; k>0; --k) {
				if (getValue(alignSeq[i], k - 1) != nilVertex) {
					addEdge(gOut, l - count, l);
					++count;
				}
			}
		}
	}
}
   
//////////////////////////////////////////////////////////////////////////////
// T-Coffee Library Reading / Writing
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////
// FUNCTIONS
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template<typename TFile, typename TChar, typename TStringSet, typename TCargo, typename TSpec, typename TString>
inline void
_parse_readSequenceData(TFile & file,
						TChar & c,
						Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
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

//////////////////////////////////////////////////////////////////////////////

template<typename TFile, typename TStringSet, typename TCargo, typename TSpec>
inline void
_readLibrary(TFile & file,
			 Graph<Alignment<TStringSet, TCargo, TSpec> >& g)
{
	SEQAN_CHECKPOINT
	typedef unsigned int TWord;
	typedef typename Value<TFile>::Type TValue;
	
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef typename Id<TGraph>::Type TId;

	TValue c;
	bool seq1ToN = false;
	if (_streamEOF(file)) return;
	else c = _streamGet(file);

	typedef std::pair<unsigned int, unsigned int> TSeqRes;
	typedef std::map<TSeqRes, unsigned int> TNodeMap;
	TNodeMap node_map;
	TWord seq1 = 0;
	TWord seq2 = 0;
	while (!_streamEOF(file)) {
		_parse_skipWhitespace(file,c);
		if (_streamEOF(file)) break;
		if (c == '#') {
			c = _streamGet(file);
			_parse_skipWhitespace(file,c);
			seq1 = _parse_readNumber(file, c);
			seq2 = _parse_readNumber(file, c);
			if (empty(g))
				if ((seq1 != 0) || (seq2 != 0)) seq1ToN = true;
			if (seq1ToN) {
				--seq1;
				--seq2;
			}
		} else if (c == '!') {
			_parse_skipLine(file, c);
		} else {
			unsigned int res1 = _parse_readNumber(file, c);
			_parse_skipWhitespace(file,c);
			unsigned int res2 = _parse_readNumber(file, c);
			_parse_skipWhitespace(file,c);
			unsigned int weight = _parse_readNumber(file, c);
			_parse_skipLine(file,c);
		
			// Insert new vertex if necessary
			--res1;
			--res2;
			bool newEdge = false;
			TSeqRes key = std::make_pair(seq1, res1);
			TNodeMap::iterator nodePos = node_map.find(key);
			TId id1;
			if (nodePos == node_map.end()) {
				id1 = addVertex(g, seq1, res1, 1); 
				node_map.insert(std::make_pair(key, id1));
				newEdge = true;
			} else {
				id1 = nodePos->second;
			}

			key = std::make_pair(seq2, res2);
			nodePos = node_map.find(key);
			TId id2;
			if (nodePos == node_map.end()) {
				id2 = addVertex(g, seq2, res2, 1); 
				node_map.insert(std::make_pair(key, id2));
				newEdge = true;
			} else {
				id2 = nodePos->second;
			}

			// Insert a new edge or adapt the weight
			if (newEdge) addEdge(g,id1,id2,weight);
			else {
				TEdgeDescriptor e = findEdge(g,id1,id2);
				if( e == 0 ) addEdge(g,id1,id2,weight);
				else cargo(e) += weight;
			}
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
	typedef unsigned int TWord;
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename Position<TFile>::Type TPosition;
	typedef typename Value<TFile>::Type TValue;
	typedef typename Id<TGraph>::Type TIdType;

	TValue c;
	if (_streamEOF(file)) return;
	else c = _streamGet(file);

	// Ignore first line
	_parse_skipLine(file, c);
	_parse_skipWhitespace(file, c);
	
	// Read number of sequences
	TWord nSeq = (TWord) _parse_readNumber(file, c);
	resize(stringSet(g), nSeq);

	// Read sequences
	for(TWord i=0; i<nSeq; ++i) {
		_parse_skipWhitespace(file, c);
		std::cout << _parse_readIdentifier(file, c) << ", ";
		_parse_skipWhitespace(file, c);
		std::cout << _parse_readNumber(file, c) << ", ";
		_parse_skipWhitespace(file, c);
		_parse_readSequenceData(file,c,g,stringSet(g)[i]);
		//std::cout << id << ", ";
		std::cout << getValueById(stringSet(g), i) << std::endl;
		//SEQAN_ASSERT(id < nSeq)		
	}
	// Reinitialize the graph, because we changed the sequences
	clearVertices(g);

	// Read library
	_readLibrary(file,g);
}

template<typename TFile, typename TStringSet, typename TCargo, typename TSpec>
void write(TFile & file, 
		   Graph<Alignment<TStringSet, TCargo, TSpec> > const& g,
		   TCoffeeLib) 
{
	SEQAN_CHECKPOINT
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;\
	typedef typename Value<TStringSet>::Type TString;
	typedef typename Size<TStringSet>::Type TSize;

	_streamWrite(file, "! TC_LIB_FORMAT_01\n");
	TSize len = length(getStringSet(g));
	_streamPutInt(file, len);
	_streamPut(file, '\n');
	for(TSize i=0;i<len;++i) {
		_streamWrite(file, "seq");
		_streamPutInt(file, i);
		_streamPut(file, ' ');
		TString str = value(getStringSet(g), i);
		_streamPutInt(file, length(str));
		_streamPut(file, ' ');
		_streamWrite(file, str);
		_streamPut(file, '\n');
	}


	typedef std::pair<unsigned int, unsigned int> TSeq;
	typedef Triple<unsigned int, unsigned int, unsigned int> TData;
	typedef std::multimap<TSeq, TData> TMap;
	TMap m;

	typedef typename Iterator<TGraph, EdgeIterator>::Type TIter;
	TIter it(g);
	for(;!atEnd(it);++it) {
		TVertexDescriptor sV = sourceVertex(it);
		TVertexDescriptor tV = targetVertex(it);
		if (sequenceId(g,sV) > sequenceId(g,tV)) {
			TVertexDescriptor tmp = sV;
			sV = tV;
			tV = tmp;
		}
		m.insert(std::make_pair(TSeq(sequenceId(g,sV), sequenceId(g,tV)), 
								TData(segmentBegin(g,sV) + 1, segmentBegin(g,tV) + 1, getCargo(*it))));
	}
	TSeq old;
	for(TMap::iterator pos = m.begin();pos!=m.end();++pos) {
		if (old != pos->first) {
			old = pos->first; 
			_streamPut(file, '#');
			_streamPutInt(file, pos->first.first);
			_streamPut(file, ' ');
			_streamPutInt(file, pos->first.second);
			_streamPut(file, '\n');		
		}
		_streamPutInt(file, pos->second.i1);
		_streamPut(file, ' ');
		_streamPutInt(file, pos->second.i2);
		_streamPut(file, ' ');
		_streamPutInt(file, pos->second.i3);
		_streamPut(file, '\n');	
	}
	_streamWrite(file, "! SEQ_0_TO_N-1");
	_streamPut(file, '\n');	
}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
