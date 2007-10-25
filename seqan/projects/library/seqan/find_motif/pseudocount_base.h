#ifndef SEQAN_HEADER_PSEUDOCOUNT_BASE_H
#define SEQAN_HEADER_PSEUDOCOUNT_BASE_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////

/**
.Class.Pseudocount:
..summary:Holds the pseudocounts for each residue of a given sequence alphabet.
..cat:Motif Finding
..signature:Pseudocount<TValue, TSpec>
..param.TValue:The type of sequence which is considered.
...type:Spec.Dna
...type:Spec.AminoAcid
..param.TSpec:Specialization tag for determining the pseudocount method.
...type:Spec.CMode 
...type:Spec.PMode
*/

template<typename TValue, typename TSpec>
class Pseudocount;

/////////////////////////////////////////////////////////////////////////////

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...