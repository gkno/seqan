#ifndef SEQAN_HEADER_GRAPH_REFINE_SCORING_H
#define SEQAN_HEADER_GRAPH_REFINE_SCORING_H


namespace SEQAN_NAMESPACE_MAIN
{

	



//fake score function 
template<typename TScoreValue,typename TStringSet,typename TAlign,typename TValue, typename TSize>
TScoreValue
getScore(TScoreValue &,
		 TStringSet &,
		 TAlign &,
		 TValue,
		 TValue,
		 TSize,
		 TSize)
{
SEQAN_CHECKPOINT
	return 1;
}				




}
#endif //#ifndef SEQAN_HEADER_...
