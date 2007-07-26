 /*==========================================================================
                SeqAn - The Library for Sequence Analysis
                          http://www.seqan.de 
 ============================================================================
  Copyright (C) 2007

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  Lesser General Public License for more details.

 ============================================================================
  $Id$
 ==========================================================================*/

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
