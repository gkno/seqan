#ifndef SEQAN_HEADER_GRAPH_NUSSINOV_H
#define SEQAN_HEADER_GRAPH_NUSSINOV_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// ScoreNussinov: Nussinov RNA folding score
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////

struct ScoreNussinov;

//////////////////////////////////////////////////////////////////////////////

template <typename TValue>
class Score<TValue, ScoreNussinov>
{
public:
	TValue data_map;

public:
	Score(TValue _map):
		data_map(_map)
	{
	}

	Score(Score const & other):
		data_map(other.data_map)
	{
	}

	~Score()
	{
	}

	Score & operator = (Score const & other)
	{
		data_map = other.data_map;
		return *this;
	}

//____________________________________________________________________________
};
//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename T>
inline unsigned int
score(Score<TValue, ScoreNussinov> const& me,
	  T const & left,
	  T const & right)
{
	typename TValue::const_iterator pos;
	if ((pos = me.data_map.find(std::make_pair(left,right))) == me.data_map.end()) return 0;
	else return pos->second;
}




//////////////////////////////////////////////////////////////////////////////
// Alignment: Nussinov RNA folding
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////

template <typename TCargo, typename TSpec, typename TScoreValue, typename TSpec2, typename TString, typename TScore, typename TPos>
void
_align_nussinov_trace(Graph<Undirected<TCargo, TSpec> >& g,
					  String<TScoreValue, TSpec2> const& mat,
					  TString const& str,
					  TScore const& sc,
					  TPos const i,
					  TPos const j) 
{
	SEQAN_CHECKPOINT

	typedef typename Size<TString>::Type TSize;
	typedef typename Value<TString>::Type TCharacter;

	TSize len = length(str);
	if (i<j) {
		if (getValue(mat, i*len+j) == getValue(mat, (i+1)*len+j)) _align_nussinov_trace(g,mat,str,sc,i+1,j);
		else if (getValue(mat, i*len+j) == getValue(mat, i*len+(j-1))) _align_nussinov_trace(g,mat,str,sc,i,j-1);
		else if (getValue(mat, i*len+j) == getValue(mat, (i+1)*len+(j-1))+score(sc,str[i],str[j])) {
			addEdge(g,i,j);
			_align_nussinov_trace(g,mat,str,sc,i+1,j-1);
		}
		else {
			for(TSize k=i+1;k<j-1;++k) {
				if (getValue(mat, i*len+j)==getValue(mat, i*len+k) + getValue(mat,(k+1)*len+j)) {
					_align_nussinov_trace(g,mat,str,sc,i,k);
					_align_nussinov_trace(g,mat,str,sc,k+1,j);
				}
			}
		}
	}
}

//////////////////////////////////////////////////////////////////////////////

template <typename TScoreValue, typename TSpec, typename TString, typename TScore, typename TLoopSize>
TScoreValue
_align_nussinov(String<TScoreValue, TSpec>& mat,
				TString const& str,
				TScore const& sc,
				TLoopSize const& lSize) 
{
	SEQAN_CHECKPOINT
	typedef typename Size<TString>::Type TSize;
	typedef typename Value<TString>::Type TCharacter;

	// Initialization
	TSize len = length(str);
	fill(mat, len * len, 0);
	TScoreValue maxVal = 0;
	TScoreValue tmp = 0;
	
	// Recursion
	for(TSize n = 1; n < len; ++n) {
		for(TSize j=n; j < len; ++j) {
			TSize i = j - n;

			// Get the new maximum	
			maxVal = getValue(mat, (i+1)*len+j);
			if ((tmp = getValue(mat, i*len + (j-1))) > maxVal) maxVal = tmp;
			if ((j-i> (TSize) lSize) && (tmp = getValue(mat, (i+1)*len + (j-1)) + score(sc, str[i], str[j])) > maxVal) maxVal = tmp;
			for(TSize k = i+1;k+1<j;++k) {
				if ((tmp = getValue(mat, i*len + k) + getValue(mat,(k+1)*len+j)) > maxVal) maxVal = tmp;
			}
			assignValue(mat, i*len + j, maxVal);
		}
	}

	//// Debug code
	//for(TSize i = 0; i < len; ++i) {
	//	for(TSize j=0; j < len; ++j) {
	//		std::cout << getValue(mat, i*len+j) << ',';
	//	}
	//	std::cout << std::endl;
	//}

	return getValue(mat, 0*len+(len-1));
}


//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec, typename TString, typename TScore, typename TLoopSize>
unsigned int
rnaFolding(Graph<Undirected<TCargo, TSpec> >& g,
		   TString const& str,
		   TScore const& sc,
		   TLoopSize const& lSize,
		   Nussinov)
{
	SEQAN_CHECKPOINT
	typedef typename Size<TString>::Type TSize;
	String<unsigned int> mat;
	unsigned int maxVal = _align_nussinov(mat, str, sc, lSize);
	addVertex(g);
	for(TSize i = 1; i<length(str);++i) {
		addVertex(g);
		addEdge(g, i-1, i);
	}
	//addEdge(g, (unsigned int) length(str) - 1, (unsigned int) 0);
	_align_nussinov_trace(g, mat, str, sc, (unsigned int) 0, (unsigned int) length(str)-1);
	return maxVal;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec, typename TScore, typename TString>
unsigned int
rnaFolding(Graph<Undirected<TCargo, TSpec> >& g,
		   TString const& str,
		   TScore const& sc,
		   Nussinov)
{
	SEQAN_CHECKPOINT
	return rnaFolding(g,str,sc,0,Nussinov());
}


}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
