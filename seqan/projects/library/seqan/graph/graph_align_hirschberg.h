#ifndef SEQAN_HEADER_GRAPH_ALIGN_HIRSCHBERG_H
#define SEQAN_HEADER_GRAPH_ALIGN_HIRSCHBERG_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Alignment: Hirschberg Alignment - TODO!!!
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template <typename TTrace, typename TStringSet, typename TScoreValue>
TScoreValue
_align_hirschberg(TTrace& trace,
		  TStringSet const& str,
		  Score<TScoreValue, Simple> const& sc,
		  Hirschberg) 
{
	SEQAN_CHECKPOINT
/*
		typedef typename Value<TTrace>::Type TTraceValue;
	typedef typename Value<TStringSet const>::Type TString;
	typedef typename Infix<TString const>::Type TInfix;
	typedef typename Size<TString>::Type TSize;
	
	// TraceBack values
	enum {Diagonal = 0, Horizontal = 1, Vertical = 2};

	TSize len1 = length(str[0]);
	TSize len2 = length(str[1]);
	typedef std::pair<TSize, TSize> TPoint;
	typedef std::list<TPoint> TMidPointQueue;
	TMidPointQueue midpoints;
	TSize x1 = 0;
	TSize x2 = 0;
	TSize y1 = 0;
	TSize y2 = 0;
	bool firstRun = true;
	bool lastRun = false;
	do {
		// Iterator until the number of midpoints equals the sequence length
		if (!(midpoints.size() < len1-1)) lastRun = true;
		firstRun = true;
		TMidPointQueue::iterator it = midpoints.begin();
		do {
			// Step1: Get the alignment window
			if (firstRun) {
				x1 = 0;
				y1 = 0;
				firstRun = false;
			} else {
				x1 = x2;
				y1 = y2;
				++it;
			}

			if (it == midpoints.end()) {
				x2 = len1;
				y2 = len2;
			} else {
				TPoint p = *it;
				x2 = p.first;
				y2 = p.second;
			}

			// Step2: Do the alignment
			TSize middle = x1 + ((x2-x1)/2);
			TInfix infix1 = infix(str[0], x1, x2);
			TInfix infix2 = infix(str[1], y1, y2);
		
			std::cout << infix1 << std::endl;
			std::cout << infix2 << std::endl;
			std::cout << "--" << std::endl;
			
			// Step3: Insert the new midpoint
			if ((x2 - x1) > 1) {
				it = midpoints.insert(it, std::make_pair(middle,y1 + ((y2-y1)/2)));
				++it;
			}
		} while (it != midpoints.end());
	} while(!lastRun);
	*/
	return 0;
}


//////////////////////////////////////////////////////////////////////////////

template<typename TAlign, typename TStringSet, typename TScoreValue, typename TTag>
TScoreValue
_globalAlignment(TAlign& align,
		 TStringSet const& str,
		 Score<TScoreValue, Simple> const& sc,
		 TTag)
{
	SEQAN_CHECKPOINT
	String<TraceBack> trace;
	return _align_hirschberg(trace, str, sc, TTag());
}


}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
