#ifndef SEQAN_HEADER_SCORE_BASE_H
#define SEQAN_HEADER_SCORE_BASE_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////

struct Simple;

//////////////////////////////////////////////////////////////////////////////


/**
.Class.Score:
..cat:Miscellaneous
..summary:A scoring scheme.
..signature:Score<TValue, TSpec>
..param.TValue:The value type.
...default:int
..param.TSpec:The specializing type.
...default:@Spec.Simple@
*/
template <typename TValue = int, typename TSpec = Simple>
class Score;

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSpec>
struct Value< Score<TValue, TSpec> >
{
	typedef TValue Type;
};

//////////////////////////////////////////////////////////////////////////////

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
