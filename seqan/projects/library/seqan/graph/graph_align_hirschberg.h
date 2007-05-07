#ifndef SEQAN_HEADER_GRAPH_ALIGN_HIRSCHBERG_H
#define SEQAN_HEADER_GRAPH_ALIGN_HIRSCHBERG_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Alignment: Hirschberg Alignment - TODO!!!
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template <typename TTrace, typename TInfix, typename TScoreValue, typename TTag>
TScoreValue
_align_recursive_hirschberg(TTrace& trace,
							TInfix const& infix1,
							TInfix const& infix2,
							Score<TScoreValue, Simple> const& sc,
							TTag) 
{
	SEQAN_CHECKPOINT

	// TraceBack values
	enum {Diagonal = 0, Horizontal = 1, Vertical = 2};

	typedef typename Size<TInfix>::Type TSize;
	TSize len1 = length(infix1);
	TSize len2 = length(infix2);

	TSize middle = len1 / 2;
	TSize cut = len2 / 2;

	if (len1 == 1) return 1;
	else if (len2 == 1) return 1;
	
	TInfix newInfix1a = infix(host(infix1), beginPosition(infix1), beginPosition(infix1) + middle);
	TInfix newInfix2a = infix(host(infix2), beginPosition(infix2), beginPosition(infix2) + cut);
	TInfix newInfix1b = infix(host(infix1), beginPosition(infix1) + middle, endPosition(infix1));
	TInfix newInfix2b = infix(host(infix2), beginPosition(infix2) + cut, endPosition(infix2) );

	/*
	std::cout << middle << std::endl;
	std::cout << cut << std::endl;
	std::cout << infix1 << std::endl;
	std::cout << infix2 << std::endl;
	std::cout << newInfix1a << std::endl;
	std::cout << newInfix2a << std::endl;
	std::cout << newInfix1b << std::endl;
	std::cout << newInfix2b << std::endl;
	std::cout << "--" << std::endl;
	*/

	return _align_recursive_hirschberg(trace, newInfix1a, newInfix2a, sc, TTag()) +
	  _align_recursive_hirschberg(trace, newInfix1b, newInfix2b, sc, TTag());
}

//////////////////////////////////////////////////////////////////////////////

template <typename TTrace, typename TStringSet, typename TScoreValue, typename TTag>
TScoreValue
_align_hirschberg(TTrace& trace,
		  TStringSet const& str,
		  Score<TScoreValue, Simple> const& sc,
		  TTag) 
{
	SEQAN_CHECKPOINT
	typedef typename Value<TTrace>::Type TTraceValue;
	typedef typename Value<TStringSet const>::Type TString;
	typedef typename Infix<TString const>::Type TInfix;
	typedef typename Size<TString>::Type TSize;


	TInfix infix1 = infix(str[0], 0, length(str[0]));
	TInfix infix2 = infix(str[1], 0, length(str[1]));
	return _align_recursive_hirschberg(trace, infix1, infix2, sc, TTag());
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
