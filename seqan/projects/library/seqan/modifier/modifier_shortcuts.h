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

typedef ModView< FunctorComplement<Dna> >	ModComplementDna;
typedef ModView< FunctorComplement<Dna5> >	ModComplementDna5;

//////////////////////////////////////////////////////////////////////////////

typedef ModifiedString<DnaString, ModView< FunctorComplement<Dna> > >		DnaStringComplement;
typedef ModifiedString<Dna5String, ModView< FunctorComplement<Dna5> > >		Dna5StringComplement;

//////////////////////////////////////////////////////////////////////////////

typedef ModifiedString<DnaString, ModReverse>		DnaStringReverse;
typedef ModifiedString<Dna5String, ModReverse>		Dna5StringReverse;

//////////////////////////////////////////////////////////////////////////////
/*
typedef ModifiedString<DnaStringReverse, ModComplementDna>		DnaStringReverseComplement;
typedef ModifiedString<Dna5StringReverse, ModComplementDna5>	Dna5StringReverseComplement;
*/

typedef ModifiedString<
			ModifiedString<DnaString, ModView< FunctorComplement<Dna> > >, 
			ModReverse
		>	DnaStringReverseComplement;

typedef ModifiedString<
			ModifiedString<	Dna5String, ModView< FunctorComplement<Dna5> > >, 
			ModReverse
		>	Dna5StringReverseComplement;

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
