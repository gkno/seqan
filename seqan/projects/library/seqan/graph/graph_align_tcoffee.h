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
.Tag.LocalPairwise_Library
..summary:Tag to specify the type of the library.
..value.LocalPairwise_Library:Use of a pairwise local library.
*/
struct LocalPairwise_Library_;
typedef Tag<LocalPairwise_Library_> const LocalPairwise_Library;

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
	TGraphType const* data_graph;
	TValue data_inf;
	typename VertexDescriptor<TGraphType>::Type data_nilVertex;

public:
	Score(TGraphType const& _graph) {
		data_graph = &_graph;
		data_inf = infimumValue<TValue>();
		data_nilVertex = getNil<typename VertexDescriptor<TGraphType>::Type>();
	}

	Score(Score const & other) {
		data_graph = other.data_graph;
		data_inf = other.data_inf;
		data_nilVertex = other.data_nilVertex;
	}

	Score & operator = (Score const & other) {
		if (this == &other) return *this;
		data_graph = other.data_graph;
		data_inf = other.data_inf;
		data_nilVertex = other.data_nilVertex;
		return *this;
	}
};

template <typename TValue, typename TGraphType>
inline TValue
scoreGapExtend(Score<TValue, ScoreAlignmentGraph<TGraphType> > &)
{
	return 0;
}


template <typename TValue, typename TGraphType>
inline TValue const
scoreGapExtend(Score<TValue, ScoreAlignmentGraph<TGraphType> > const &)
{
	return 0;
}

template <typename TValue, typename TGraphType>
inline TValue
scoreGapOpen(Score<TValue, ScoreAlignmentGraph<TGraphType> > &)
{
	return 0;
}
template <typename TValue, typename TGraphType>
inline TValue const
scoreGapOpen(Score<TValue, ScoreAlignmentGraph<TGraphType> > const &)
{
	return 0;
}

template <typename TValue, typename TGraphType, typename T>
inline TValue
score(Score<TValue, ScoreAlignmentGraph<TGraphType> > & me,
	  T const & left,
	  T const & right)
{
	typedef typename EdgeDescriptor<TGraphType>::Type TEdgeDescriptor;
	typedef typename VertexDescriptor<TGraphType>::Type TVertexDescriptor;
	typedef typename Cargo<TGraphType>::Type TCargo;
	typedef typename Iterator<TGraphType, OutEdgeIterator>::Type TOutEdgeIterator;
	typedef typename Iterator<T const>::Type TStringIter;

	TValue sum = 0;
	TStringIter itLeftEnd = end(left);
	TStringIter itRightEnd = end(right);
	for(TStringIter itLeft = begin(left);itLeft != itLeftEnd;++itLeft) {
		if (*itLeft != me.data_nilVertex) {
			for(TStringIter itRight = begin(right);itRight != itRightEnd;++itRight) {
				if (*itRight != me.data_nilVertex) {
					TEdgeDescriptor e = findEdge(*me.data_graph, *itLeft, *itRight);
					if (e != 0) sum += getCargo(e);
				}
			}
		}
	}

	// Did we found at least one edge?
	if (sum > 0) return sum / ((TValue) length(left) * (TValue) length(right));
	else return me.data_inf;
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
// Graph: T-Coffee - Simple q-gram counter
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template<typename TString, typename TTupelString, typename TAlphabet>
inline void
_getTupelString(TString const& str, TTupelString& tupelString, unsigned int const ktup, TAlphabet) {
	typedef unsigned int TWord;
	typedef typename Size<TString>::Type TSize;

	// Alphabet size
	TSize alphabet_size = ValueSize<TAlphabet>::VALUE;

	// Assign a unique number to each k-tupel
	String<TWord> prod;  // Scaling according to position in k-tupel
	resize(prod,ktup);
	for (TWord i=0; i<ktup;++i) prod[ktup-i-1]=(TWord)pow((double)alphabet_size,(double)i);

	TSize len = length(str);
	clear(tupelString);
	resize(tupelString, len-(ktup - 1)); 
	TSize tupelIndex = 0;
	TSize endTupel = 0;
	tupelString[tupelIndex] = 0;
	for(;endTupel<ktup;++endTupel) {
		tupelString[tupelIndex] += ((unsigned int) ((TAlphabet) str[endTupel])) * prod[endTupel];
	}
	++tupelIndex;
	for(;endTupel<len;++endTupel) {
		tupelString[tupelIndex] = tupelString[tupelIndex - 1];
		tupelString[tupelIndex] -= ((unsigned int) ((TAlphabet) str[endTupel - ktup])) * prod[0];
		tupelString[tupelIndex] *= alphabet_size;
		tupelString[tupelIndex] += ((unsigned int) ((TAlphabet) str[endTupel]));
		++tupelIndex;
	}
}


template<typename TString, typename TSpec, typename THitMatrix, typename TSize, typename TAlphabet>
void
getCommonKmerMatrix(StringSet<TString, TSpec> const& strSet, THitMatrix& mat, TSize ktup, TAlphabet) {
	SEQAN_CHECKPOINT
	typedef unsigned int TWord;
	typedef String<TWord> TTupelString;
	typedef String<TTupelString> TTupelStringSet;

	// Number of sequences
	TSize nseq = length(strSet);
	TSize alphabet_size = ValueSize<TAlphabet>::VALUE;

	// Initialization
	// Matrix for common k-tupels between sequence i and j
	setDimension(mat, 2);setLength(mat, 0, nseq);setLength(mat, 1, nseq);
	fill(host(mat), nseq*nseq, 0);

	// Transform the set of strings into a set of strings of k-tupels
	TTupelStringSet tupSet;
	resize(tupSet, length(strSet));
	for(TSize k=0;k<(TSize) length(strSet);++k) _getTupelString(strSet[k], tupSet[k], ktup, TAlphabet());

	// Build for each sequence the q-gram Index and count common hits
	String<TWord> qIndex;
	String<TWord> compareIndex;
	for(TSize k=0;k<nseq;++k) {
		clear(qIndex);
		fill(qIndex, (unsigned int) pow((double)alphabet_size, (double)ktup), (TWord) 0, Exact());
		for(TSize i = 0;i < (TSize) length(tupSet[k]);++i) ++qIndex[ tupSet[k][i] ];
		TWord value;
	    for (TSize k2=k; k2<nseq; ++k2) {
			clear(compareIndex);
			fill(compareIndex, (unsigned int) pow((double)alphabet_size, (double)ktup), (TWord) 0, Exact());
			value = 0;
			for(TSize i = 0;i < (TSize) length(tupSet[k2]);++i) {
				//std::cout << tupSet[k2][i] << "," << compareIndex[ tupSet[k2][i] ] << "," << qIndex[ tupSet[k2][i] ]<< std::endl;
				if (compareIndex[ tupSet[k2][i] ] < qIndex[ tupSet[k2][i] ]) ++value;
				++compareIndex[ tupSet[k2][i] ];
			}
			assignValue(mat, k*nseq+k2, value);
		}
	}

	// Copy upper triangle to lower triangle
	for(TWord k = 0; k < (TWord) nseq; ++k) {
		for(TWord k2 = k + 1; k2 < (TWord) nseq; ++k2) {
			assignValue(mat, k2*nseq+k, getValue(mat, k*nseq+k2));
		}
	}

	//// Debug code
	//for (TSize row=0;row<nseq;++row) {
	//	for(TSize col=0;col<nseq;++col) {
	//		std::cout << getValue(mat, row*nseq+col) << ",";
	//	}
	//	std::cout << std::endl;
	//}
}



/*
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
	

	// Recode the strings into amino acid groups
	for(TStringSetSize k=0;k<length(strSet);++k) {
		indexText(index)[k] = strSet[k];
		if ( (TSize) length(strSet[k]) < ktup) ktup = length(strSet[k]) - 1;
	}
	resize(indexShape(index), (unsigned) ktup);
	
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

	for (TWord row=0;row<nseq;++row) {
		for(TWord col=0;col<nseq;++col) {
			std::cout << getValue(mat, row*nseq+col) << ",";
		}
		std::cout << std::endl;
	}
}
*/

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


template<typename TStringSet, typename TCargo, typename TSpec, typename TScore>
void 
generatePrimaryLibrary(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
					   TScore const& score_type,
					   LocalPairwise_Library)
{
	SEQAN_CHECKPOINT
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename Size<TGraph>::Type TSize;
	typedef typename Id<TGraph>::Type TId;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef typename Iterator<TGraph, EdgeIterator>::Type TEdgeIterator;
	typedef typename Value<TScore>::Type TScoreValue;

	// Clear graph
	clearVertices(g);

	// Pairwise alignments for all pairs of sequences
	TStringSet& str = stringSet(g);	
	TSize nseq = length(str);

	// String of fragments to combine all pairwise alignments into a multiple alignment
	typedef Fragment<> TFragment;
	typedef String<TFragment, Block<> > TFragmentString;
	typedef typename Iterator<TFragmentString>::Type TFragmentStringIter;
	TFragmentString matches;
	String<TScoreValue, Block<> > score_values;

	for(TSize i=0; i<nseq; ++i) {
		for(TSize j=i+1; j<nseq; ++j) {
			TId id1 = positionToId(str, i);
			TId id2 = positionToId(str, j);
			// Pairwise alignment graph
			TStringSet pairSet;
			assignValueById(pairSet, str, id1);
			assignValueById(pairSet, str, id2);
			typedef Graph<Alignment<TStringSet, unsigned int> > TPairGraph;
			typedef typename VertexDescriptor<TPairGraph>::Type TVD;
			typedef typename Iterator<TPairGraph, EdgeIterator>::Type TEI;
			TPairGraph pGraph(pairSet);
			localAlignment(pGraph, score_type, SmithWatermanClump() );
			
			//fstream strm02; // Alignment graph as dot
			//strm02.open(TEST_PATH "my_test.dot", ios_base::out | ios_base::trunc);
			//write(strm02,pGraph,DotDrawing());
			//strm02.close();

			TEI it(pGraph);
			for(;!atEnd(it);++it) {
				TVD sV = sourceVertex(it);
				TVD tV = targetVertex(it);
				push_back(matches, TFragment( (unsigned int) sequenceId(pGraph, sV), (unsigned int) fragmentBegin(pGraph,sV), (unsigned int) sequenceId(pGraph, tV),  (unsigned int)  fragmentBegin(pGraph,tV),  (unsigned int)  fragmentLength(pGraph,tV)));
				push_back(score_values, cargo(*it));
			}
		}
	}

	// Refine all matches and create multiple alignment
	matchRefinement(matches,str,const_cast< TScore&>(score_type),g);

	// Adapt edge weights
	TFragmentStringIter endIt = end(matches);
	for(TFragmentStringIter it = begin(matches); it != endIt; ++it) {
		TId id1 = sequenceId(*it,0);
		TId id2 = sequenceId(*it,1);
		TSize pos1 = fragmentBegin(*it, id1);
		TSize pos2 = fragmentBegin(*it, id2);
		TSize end1 = pos1 + fragmentLength(*it, id1);
		while(pos1 < end1) {
			SEQAN_TASSERT(pos2 < pos2 + fragmentLength(*it, id2))
			TVertexDescriptor p1 = findVertex(g, id1, pos1);
			TVertexDescriptor p2 = findVertex(g, id2, pos2);
			TEdgeDescriptor e = findEdge(g, p1, p2);
			if (e != 0) cargo(e) *= (TCargo) getValue(score_values, position(it));
			else addEdge(g, p1, p2, (TCargo) getValue(score_values, position(it)));
			SEQAN_TASSERT(fragmentLength(g, p1) == fragmentLength(g, p2))
			pos1 += fragmentLength(g, p1);
			pos2 += fragmentLength(g, p2);
		}
	}
}


template<typename TStringSet, typename TCargo, typename TSpec, typename TScore>
void 
generatePrimaryLibrary(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
					   TScore const& score_type,
					   GlobalPairwise_Library)
{
	SEQAN_CHECKPOINT

	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename Size<TGraph>::Type TSize;
	typedef typename Id<TGraph>::Type TId;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef typename Iterator<TGraph, EdgeIterator>::Type TEdgeIterator;
	typedef typename Value<TScore>::Type TScoreValue;

	// Clear graph
	clearVertices(g);

	// Pairwise alignments for all pairs of sequences
	TStringSet& str = stringSet(g);	
	TSize nseq = length(str);
	
	// String of fragments to combine all pairwise alignments into a multiple alignment
	typedef Fragment<> TFragment;
	typedef String<TFragment, Block<> > TFragmentString;
	typedef typename Iterator<TFragmentString>::Type TFragmentStringIter;
	TFragmentString matches;
	String<TScoreValue, Block<> > score_values;

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
			globalAlignment(pGraph, score_type, Gotoh() );
			
			// Determine a sequence weight
			TCargo seqSim = (TCargo) (_getSequenceSimilarity(pGraph, typename Value<TStringSet>::Type() ) * 100);

			//// Debug code
			//std::cout << pairSet[0] << std::endl;
			//std::cout << pairSet[1] << std::endl;
			//std::cout << pGraph << std::endl;
			//std::cout << seqSim << std::endl;
			
			// Extract the matches		
			TEI it(pGraph);
			for(;!atEnd(it);++it) {
				TVD sV = sourceVertex(it);
				TVD tV = targetVertex(it);
				push_back(matches, TFragment( (unsigned int) sequenceId(pGraph, sV), (unsigned int) fragmentBegin(pGraph,sV), (unsigned int) sequenceId(pGraph, tV),  (unsigned int)  fragmentBegin(pGraph,tV),  (unsigned int)  fragmentLength(pGraph,tV)));
				push_back(score_values, seqSim);
			}
		}
	}

	//// Debug Code
	//// Print all the matches
	//std::cout << "The sequences:" << std::endl;
	//for(TSize i = 0;i<length(str);++i) {
	//	std::cout << positionToId(str,i) << ':' << str[i] << std::endl;
	//}
	//std::cout << "The matches:" << std::endl;
	//for(TSize i = 0;i<length(matches);++i) {
	//	TId tmp_id1 = sequenceId(matches[i],0);
	//	std::cout << tmp_id1 << ',' << fragmentBegin(matches[i],tmp_id1) << ',';
	//	for(TSize j = fragmentBegin(matches[i],tmp_id1); j < fragmentBegin(matches[i],tmp_id1) + fragmentLength(matches[i],tmp_id1); ++j) {
	//		std::cout << str[idToPosition(str, tmp_id1)][j];
	//	}
	//	TId tmp_id2 = sequenceId(matches[i],1);
	//	std::cout << ',' <<	tmp_id2 << ',' << fragmentBegin(matches[i],tmp_id2) << ',';
	//	for(TSize j = fragmentBegin(matches[i],tmp_id2); j < fragmentBegin(matches[i],tmp_id2) + fragmentLength(matches[i],tmp_id2); ++j) {
	//		std::cout << str[idToPosition(str, tmp_id2)][j];
	//	}
	//	std::cout << std::endl;
	//}

	// Refine all matches and create multiple alignment
	matchRefinement(matches,stringSet(g),const_cast<TScore&>(score_type),g);

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
	TFragmentStringIter endIt = end(matches);
	for(TFragmentStringIter it = begin(matches); it != endIt; ++it) {
		TId id1 = sequenceId(*it,0);
		TId id2 = sequenceId(*it,1);
		TSize pos1 = fragmentBegin(*it, id1);
		TSize pos2 = fragmentBegin(*it, id2);
		TSize end1 = pos1 + fragmentLength(*it, id1);
		while(pos1 < end1) {
			SEQAN_TASSERT(pos2 < pos2 + fragmentLength(*it, id2))
			TVertexDescriptor p1 = findVertex(g, id1, pos1);
			TVertexDescriptor p2 = findVertex(g, id2, pos2);
			TEdgeDescriptor e = findEdge(g, p1, p2);
			if (e != 0) cargo(e) *= (TCargo) getValue(score_values, position(it));
			else addEdge(g, p1, p2, (TCargo) getValue(score_values, position(it)));
			SEQAN_TASSERT(fragmentLength(g, p1) == fragmentLength(g, p2))
			pos1 += fragmentLength(g, p1);
			pos2 += fragmentLength(g, p2);
		}
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
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef typename Iterator<TGraph, VertexIterator>::Type TVertexIterator;
	typedef typename Iterator<TGraph, EdgeIterator>::Type TEdgeIterator;
	typedef typename Iterator<TGraph, OutEdgeIterator>::Type TOutEdgeIterator;

	// Two tasks:
	// 1) Add edges for the case that a and c is aligned, b and c is aligned, but a and b are not, give these edges the appropriate weight
	// 2) Augment all existing edges
	String<TCargo> newCargoMap;
	resize(newCargoMap, getIdUpperBound(_getEdgeIdManager(g)), Exact());
	TEdgeIterator it(g);
	for(;!atEnd(it);++it) assignProperty(newCargoMap, *it, cargo(*it));
	typedef std::map<std::pair<TVertexDescriptor, TVertexDescriptor>, TCargo> TNewEdgeMap;
	TNewEdgeMap edges;
	TVertexIterator itVertex(g);
	for(;!atEnd(itVertex);++itVertex) {
		TOutEdgeIterator outIt1(g, *itVertex);
		while (!atEnd(outIt1)) {
			TOutEdgeIterator outIt2 = outIt1;
			goNext(outIt2);
			while (!atEnd(outIt2)) {
				TVertexDescriptor tV1 = targetVertex(outIt1);
				TVertexDescriptor tV2 = targetVertex(outIt2);
				if (sequenceId(g, tV1) != sequenceId(g,tV2)) {
					TEdgeDescriptor e = findEdge(g, tV1, tV2);
					if (e == 0) {
						// New edge
						TCargo val = cargo(*outIt1);
						if (val > cargo(*outIt2)) val = cargo(*outIt2);
						if (tV1 < tV2) {
							typename TNewEdgeMap::iterator pos = edges.find(std::make_pair(tV1, tV2));
							if (pos == edges.end()) {
								edges.insert(std::make_pair(std::make_pair(tV1, tV2), val));
							} else {
								pos->second += val;
							}
						} else {
							typename TNewEdgeMap::iterator pos = edges.find(std::make_pair(tV2, tV1));
							if (pos == edges.end()) {
								edges.insert(std::make_pair(std::make_pair(tV2, tV1), val));
							} else {
								pos->second += val;
							}
						}
					} else {
						if (getCargo(*outIt2) > getCargo(*outIt1)) property(newCargoMap, e) += getCargo(*outIt1);
						else property(newCargoMap, e) += getCargo(*outIt2);	
					}
				}
				goNext(outIt2);
			}
			goNext(outIt1);
		}
	}
	// Assign the new weights and clean-up the cargo map
	goBegin(it);
	for(;!atEnd(it);++it) cargo(*it) = getProperty(newCargoMap, *it);
	clear(newCargoMap);

	// Finally add the new edges created by the triplet approach
	// If there are very many edges we have to thin the graph
	if (numEdges(g) < 1000000) {
		for(typename TNewEdgeMap::const_iterator pos = edges.begin();pos!=edges.end();++pos) {
			addEdge(g, pos->first.first, pos->first.second, pos->second);
		}
	} else {
		for(typename TNewEdgeMap::const_iterator pos = edges.begin();pos!=edges.end();++pos) {
			TVertexDescriptor tV1 = pos->first.first;
			TVertexDescriptor tV2 = pos->first.second;
			TId localSeqId = sequenceId(g, tV2);
			TCargo cargoSum = 0;
			unsigned int count = 0;
			TOutEdgeIterator outIt(g, tV1);
			for(;!atEnd(outIt);++outIt) {
				if (sequenceId(g, targetVertex(outIt)) == localSeqId) {
					cargoSum+=cargo(*outIt);
					++count;
				}
			}
			// If new edge weight is above average add it, otherwise discard
			if ((count == (unsigned int) 0) ||
				((TCargo)pos->second > (TCargo) (cargoSum / count) ) ) {
				addEdge(g, pos->first.first, pos->first.second, pos->second);
			}
		}
	}
	
	// Clean-up the edge map
	edges.clear();
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
	
	//Score<typename Cargo<TGraph>::Type, ScoreAlignmentGraph<TGraph> > score_type = Score<typename Cargo<TGraph>::Type, ScoreAlignmentGraph<TGraph> >(g);
	Score<int, ScoreAlignmentGraph<TGraph> > score_type = Score<int, ScoreAlignmentGraph<TGraph> >(g);
	TSequence tmp;
	globalAlignment(tmp, strSet, score_type, TTag());
	TSize len = length(tmp);
	clear(alignSeq);
	resize(alignSeq, len);

	// Traceback is backwards, so reverse everything
	for(TSize i = len; i>0;--i) alignSeq[len-i] = tmp[i-1];
}

template<typename TStringSet, typename TCargo, typename TSpec, typename TGuideTree, typename TVertexDescriptor, typename TSequence, typename TTag>
void 
_recursiveProgressiveAlignment(Graph<Alignment<TStringSet, TCargo, TSpec> >& g,
							   TGuideTree& tree,
							   TVertexDescriptor const root,
							   TSequence& alignSeq,
							   TTag)
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
		clear(strSet);clear(seq1);clear(seq2);

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
			//std::cout << "First subtree aligned" << std::endl;
			assignValueById(strSet, seq1);
			++count;
		} else if (count == 1) {
			_recursiveProgressiveAlignment(g,tree, *adjIt, seq2, TTag());
			assignValueById(strSet, seq2);
			_alignStringSetAccordingToGraph(g,strSet,alignSeq, TTag());
			//std::cout << "Second subtree aligned" << std::endl;
			clear(strSet);clear(seq1);clear(seq2);
			assignValueById(strSet, alignSeq);
			++count;
		} else {
			_recursiveProgressiveAlignment(g,tree, *adjIt, seq3, TTag());
			//std::cout << "Third subtree aligned" << std::endl;
			assignValueById(strSet, seq3);
			_alignStringSetAccordingToGraph(g,strSet,alignSeq, TTag());
			clear(strSet);clear(seq3);
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
						Graph<Alignment<TStringSet, TCargo, TSpec> >&,
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


/*
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
	typedef typename Iterator<TGraph, EdgeIterator>::Type TEdgeIterator;
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

	typedef Fragment<> TFragment;
	typedef String<TFragment, External<> > TFragmentString;
	TFragmentString matches;
	while (!atEnd(it)) 
	{
		occs = getOccurences(it);							// gives hit positions (seqNo,seqOfs)
		orderOccurences(occs);								// order them by seqNo
			
		TSize matchLen = repLength(it);
		for(TSize i = 0; i < (TSize) length(occs); ++i) {
			//std::cout << positionToId(str, i) << ',' << getValueI2(occs[i]) << ',' << matchLen << std::endl;
			if (i > 0) {
				for(TSize k = i; k>0;--k) {
					push_back(matches, TFragment( (unsigned int) positionToId(str, i), (unsigned int) getValueI2(occs[i]), (unsigned int) positionToId(str, i - k),  (unsigned int)  getValueI2(occs[i-k]),  (unsigned int)  matchLen));
				}
			}
		}
		++it;
	}

	// Refine all matches and create multiple alignment
	matchRefinement(matches,stringSet(g),g);

	// Adapt edge weights
	for(TEdgeIterator it(g);!atEnd(it);++it) {
		cargo(*it) = fragmentLength(g, sourceVertex(it));
	}
}
*/

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
