/*
 *  modifier_shortcuts.h
 *  SeqAn
 *
 *  Created by David Weese on 26.01.06.
 *
 */

#ifndef SEQAN_HEADER_MODIFIER_SHORTCUTS_H
#define SEQAN_HEADER_MODIFIER_SHORTCUTS_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////

typedef ModView< FunctorComplement<Dna> >						ModComplementDna;
typedef ModView< FunctorComplement<Dna5> >						ModComplementDna5;

//////////////////////////////////////////////////////////////////////////////

typedef ModifiedString<DnaString, ModComplementDna>				DnaStringComplement;
typedef ModifiedString<Dna5String, ModComplementDna5>			Dna5StringComplement;

//////////////////////////////////////////////////////////////////////////////

typedef ModifiedString<DnaString, ModReverse>					DnaStringReverse;
typedef ModifiedString<Dna5String, ModReverse>					Dna5StringReverse;

//////////////////////////////////////////////////////////////////////////////
/*
typedef ModifiedString<DnaStringReverse, ModComplementDna>		DnaStringReverseComplement;
typedef ModifiedString<Dna5StringReverse, ModComplementDna5>	Dna5StringReverseComplement;
*/
typedef ModifiedString<DnaStringComplement, ModReverse>			DnaStringReverseComplement;
typedef ModifiedString<Dna5StringComplement, ModReverse>		Dna5StringReverseComplement;

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

template < typename TSequence >
inline void complementInPlace(TSequence & sequence) {
	convertInPlace(sequence, FunctorComplement<typename Value<TSequence>::Type>());
} 

template < typename TSequence >
inline void reverseComplementInPlace(TSequence & sequence) {
	convertInPlace(sequence, FunctorComplement<typename Value<TSequence>::Type>());
	reverseInPlace(sequence);
} 


//////////////////////////////////////////////////////////////////////////////

}

#endif
