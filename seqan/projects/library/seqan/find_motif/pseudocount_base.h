#ifndef SEQAN_HEADER_PSEUDOCOUNT_BASE_H
#define SEQAN_HEADER_PSEUDOCOUNT_BASE_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////

/**
.Class.Pseudocount:
..summary:Holds the pseudocounts for each residue of a given sequence alphabet.
..cat:Motif Search
..signature:Pseudocount<TValue, TSpec>
..param.TValue:The type of sequence which is considered.
...metafunction:Metafunction.Value
...type:Spec.Dna
...type:Spec.AminoAcid
..param.TSpec:Specialization tag for determining the pseudocount method.
...type:Spec.CMode 
...type:Spec.PMode
*/

template<typename TValue, typename TSpec>
class Pseudocount;

//////////////////////////////////////////////////////////////////////////////
//Metafunctions
//////////////////////////////////////////////////////////////////////////////

///.Metafunction.Value.param.T.type:Class.Pseudocount

template<typename TValue, typename TSpec>
struct Value< Pseudocount<TValue, TSpec> >
{
	typedef TValue Type;
};
template<typename TValue, typename TSpec>
struct Value< Pseudocount<TValue, TSpec> const>
{
	typedef TValue const Type;
};

/////////////////////////////////////////////////////////////////////////////

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...