/*==========================================================================
                SeqAn - The Library for object Analysis
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
  Author: David Weese <david.weese@fu-berlin.de>
 ==========================================================================*/

#ifndef SEQAN_HEADER_MODIFIER_SHORTCUTS_H
#define SEQAN_HEADER_MODIFIER_SHORTCUTS_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////

typedef ModView< FunctorComplement<Dna> >	ModComplementDna;
typedef ModView< FunctorComplement<Dna5> >	ModComplementDna5;
typedef ModView< FunctorComplement<Rna> >	ModComplementRna;
typedef ModView< FunctorComplement<Rna5> >	ModComplementRna5;

//////////////////////////////////////////////////////////////////////////////

typedef ModifiedString<DnaString, ModView< FunctorComplement<Dna> > >		DnaStringComplement;
typedef ModifiedString<Dna5String, ModView< FunctorComplement<Dna5> > >		Dna5StringComplement;
typedef ModifiedString<RnaString, ModView< FunctorComplement<Rna> > >		RnaStringComplement;
typedef ModifiedString<Rna5String, ModView< FunctorComplement<Rna5> > >		Rna5StringComplement;

//////////////////////////////////////////////////////////////////////////////

typedef ModifiedString<DnaString, ModReverse>		DnaStringReverse;
typedef ModifiedString<Dna5String, ModReverse>		Dna5StringReverse;
typedef ModifiedString<RnaString, ModReverse>		RnaStringReverse;
typedef ModifiedString<Rna5String, ModReverse>		Rna5StringReverse;

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

typedef ModifiedString<
			ModifiedString<RnaString, ModView< FunctorComplement<Rna> > >, 
			ModReverse
		>	RnaStringReverseComplement;

typedef ModifiedString<
			ModifiedString<	Rna5String, ModView< FunctorComplement<Rna5> > >, 
			ModReverse
		>	Rna5StringReverseComplement;

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

/**
.Function.complementInPlace:
..cat:Modifier
..summary:Complement a sequence or a @Class.StringSet@ in-place.
..signature:complementInPlace(sequence)
..param.sequence:The sequence to complement.
...type:Class.String
...type:Class.Segment
..include:seqan/modifier.h
..see:Function.reverseComplementInPlace
..see:Function.toLowerInPlace
..see:Function.toUpperInPlace
 */
template < typename TSequence >
inline void complementInPlace(TSequence & sequence) 
{
    SEQAN_CHECKPOINT;
	convertInPlace(sequence, FunctorComplement<typename Value<TSequence>::Type>());
} 

template < typename TSequence >
inline void complementInPlace(TSequence const & sequence) 
{
    SEQAN_CHECKPOINT;
	convertInPlace(sequence, FunctorComplement<typename Value<TSequence>::Type>());
} 

/**
.Function.complementInPlace:
..signature:complementInPlace(stringSet)
..param.stringSet:The @Class.StringSet@ to complement.
...type:Class.StringSet
..include:seqan/modifier.h
 */
template < typename TSequence, typename TSpec >
inline void complementInPlace(StringSet<TSequence, TSpec> & stringSet)
{
    SEQAN_CHECKPOINT;
	unsigned seqCount = length(stringSet);
	for(unsigned seqNo = 0; seqNo < seqCount; ++seqNo)
		complementInPlace(stringSet[seqNo]);
}

template < typename TSequence, typename TSpec >
inline void complementInPlace(StringSet<TSequence, TSpec> const & stringSet)
{
    SEQAN_CHECKPOINT;
	unsigned seqCount = length(stringSet);
	for(unsigned seqNo = 0; seqNo < seqCount; ++seqNo)
		complementInPlace(stringSet[seqNo]);
}

/**
.Function.reverseComplementInPlace:
..cat:Modifier
..summary:Reverse and complement a sequence or a @Class.StringSet@ in-place.
..signature:reverseComplementInPlace(sequence)
..param.sequence:The sequence to complement.
...type:Class.String
...type:Class.Segment
..include:seqan/modifier.h
..see:Function.complementInPlace
..see:Function.toLowerInPlace
..see:Function.toUpperInPlace
 */
template < typename TSequence >
inline void reverseComplementInPlace(TSequence & sequence) 
{
    SEQAN_CHECKPOINT;
	convertInPlace(sequence, FunctorComplement<typename Value<TSequence>::Type>());
	reverseInPlace(sequence);
} 

// TODO(holtgrew): How is doing anything in-place on a const value possible?
template < typename TSequence >
inline void reverseComplementInPlace(TSequence const & sequence) 
{
    SEQAN_CHECKPOINT;
	convertInPlace(sequence, FunctorComplement<typename Value<TSequence>::Type>());
	reverseInPlace(sequence);
} 

/**
.Function.reverseComplementInPlace:
..signature:reverseComplementInPlace(stringSet)
..param.stringSet:The @Class.StringSet@ to complement.
...type:Class.StringSet
..include:seqan/modifier.h
 */
template < typename TSequence, typename TSpec >
inline void reverseComplementInPlace(StringSet<TSequence, TSpec> & stringSet)
{
    SEQAN_CHECKPOINT;
	unsigned seqCount = length(stringSet);
	for(unsigned seqNo = 0; seqNo < seqCount; ++seqNo)
		reverseComplementInPlace(stringSet[seqNo]);
}

// TODO(holtgrew): How is doing anything in-place on a const value possible?
template < typename TSequence, typename TSpec >
inline void reverseComplementInPlace(StringSet<TSequence, TSpec> const & stringSet)
{
    SEQAN_CHECKPOINT;
	unsigned seqCount = length(stringSet);
	for(unsigned seqNo = 0; seqNo < seqCount; ++seqNo)
		reverseComplementInPlace(stringSet[seqNo]);
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.toLowerInPlace:
..cat:Modifier
..summary:Convert characters in sequence or @Class.StringSet@ to lower case in-place.
..signature:toLowerInPlace(sequence)
..param.sequence:The sequence to convert into lowercase.
...type:Class.String
...type:Class.Segment
..include:seqan/modifier.h
..see:Function.toUpperInPlace
..see:Function.reverseComplementInPlace
..see:Function.complementInPlace
 */
template < typename TSequence >
inline void toLowerInPlace(TSequence & sequence) 
{
    SEQAN_CHECKPOINT;
	convertInPlace(sequence, FunctorLowcase<typename Value<TSequence>::Type>());
} 

// TODO(holtgrew): How is doing anything in-place on a const value possible?
template < typename TSequence >
inline void toLowerInPlace(TSequence const & sequence) 
{
    SEQAN_CHECKPOINT;
	convertInPlace(sequence, FunctorLowcase<typename Value<TSequence>::Type>());
} 

/**
.Function.toLowerInPlace:
..signature:toLowerInPlace(stringSet)
..param.stringSet:The @Class.StringSet@ to convert into lowercase.
...type:Class.StringSet
..include:seqan/modifier.h
 */	
template < typename TSequence, typename TSpec >
inline void toLowerInPlace(StringSet<TSequence, TSpec> & stringSet)
{
    SEQAN_CHECKPOINT;
	unsigned seqCount = length(stringSet);
	for(unsigned seqNo = 0; seqNo < seqCount; ++seqNo)
		toLowerInPlace(stringSet[seqNo]);
}

// TODO(holtgrew): How is doing anything in-place on a const value possible?
template < typename TSequence, typename TSpec >
inline void toLowerInPlace(StringSet<TSequence, TSpec> const & stringSet)
{
    SEQAN_CHECKPOINT;
	unsigned seqCount = length(stringSet);
	for(unsigned seqNo = 0; seqNo < seqCount; ++seqNo)
		toLowerInPlace(stringSet[seqNo]);
}

/**
.Function.toUpperInPlace:
..cat:Modifier
..summary:Convert characters in sequence or @Class.StringSet@ to lower case in-place.
..signature:toUpperInPlace(sequence)
..param.sequence:The sequence to convert into uppercase.
...type:Class.String
...type:Class.Segment
..include:seqan/modifier.h
..see:Function.toLowerInPlace
..see:Function.reverseComplementInPlace
..see:Function.complementInPlace
 */
template < typename TSequence >
inline void toUpperInPlace(TSequence & sequence) 
{
    SEQAN_CHECKPOINT;
	convertInPlace(sequence, FunctorUpcase<typename Value<TSequence>::Type>());
} 

// TODO(holtgrew): How is doing anything in-place on a const value possible?
template < typename TSequence >
inline void toUpperInPlace(TSequence const & sequence) 
{
    SEQAN_CHECKPOINT;
	convertInPlace(sequence, FunctorUpcase<typename Value<TSequence>::Type>());
} 

/**
.Function.toUpperInPlace:
..signature:toUpperInPlace(stringSet)
..param.stringSet:The @Class.StringSet@ to convert into uppercase.
...type:Class.StringSet
..include:seqan/modifier.h
 */	
template < typename TSequence, typename TSpec >
inline void toUpperInPlace(StringSet<TSequence, TSpec> & stringSet)
{
    SEQAN_CHECKPOINT;
	unsigned seqCount = length(stringSet);
	for(unsigned seqNo = 0; seqNo < seqCount; ++seqNo)
		toUpperInPlace(stringSet[seqNo]);
}

// TODO(holtgrew): How is doing anything in-place on a const value possible?
template < typename TSequence, typename TSpec >
inline void toUpperInPlace(StringSet<TSequence, TSpec> const & stringSet)
{
    SEQAN_CHECKPOINT;
	unsigned seqCount = length(stringSet);
	for(unsigned seqNo = 0; seqNo < seqCount; ++seqNo)
		toUpperInPlace(stringSet[seqNo]);
}

//////////////////////////////////////////////////////////////////////////////

/**
.Shortcut.ModComplementDna:
..cat:Modifier
..summary:Modifier specialization type for the complement of @Spec.Dna@ alphabet sequences.
..signature:DnaStringComplement
..shortcutfor:Spec.ModView
...signature:ModView< FunctorComplement<Dna> >
..see:Spec.ModView
..see:Class.FunctorComplement
*/

/**
.Shortcut.ModComplementDna5:
..cat:Modifier
..summary:Modifier specialization type for the complement of @Spec.Dna5@ alphabet sequences.
..signature:Dna5StringComplement
..shortcutfor:Spec.ModView
...signature:ModView< FunctorComplement<Dna5> >
..see:Spec.ModView
..see:Class.FunctorComplement
*/

/**
.Shortcut.ModComplementRna:
..cat:Modifier
..summary:Modifier specialization type for the complement of @Spec.Rna@ alphabet sequences.
..signature:RnaStringComplement
..shortcutfor:Spec.ModView
...signature:ModView< FunctorComplement<Rna> >
..see:Spec.ModView
..see:Class.FunctorComplement
*/

/**
.Shortcut.ModComplementRna5:
..cat:Modifier
..summary:Modifier specialization type for the complement of @Spec.Rna5@ alphabet sequences.
..signature:Rna5StringComplement
..shortcutfor:Spec.ModView
...signature:ModView< FunctorComplement<Rna5> >
..see:Spec.ModView
..see:Class.FunctorComplement
*/


/**
.Shortcut.DnaStringComplement:
..cat:Modifier
..summary:Modifier for the complement of a @Shortcut.DnaString@.
..signature:DnaStringComplement
..shortcutfor:Class.ModifiedString
...signature:ModifiedString<DnaString, ModView< FunctorComplement<Dna> > >
..see:Shortcut.DnaString
..see:Class.ModifiedString
..see:Spec.ModView
..see:Class.FunctorComplement
*/

/**
.Shortcut.Dna5StringComplement:
..cat:Modifier
..summary:Modifier for the complement of a @Shortcut.Dna5String@.
..signature:Dna5StringComplement
..shortcutfor:Class.ModifiedString
...signature:ModifiedString<Dna5String, ModView< FunctorComplement<Dna5> > >
..see:Shortcut.Dna5String
..see:Class.ModifiedString
..see:Spec.ModView
..see:Class.FunctorComplement
*/

/**
.Shortcut.RnaStringComplement:
..cat:Modifier
..summary:Modifier for the complement of a @Shortcut.RnaString@.
..signature:RnaStringComplement
..shortcutfor:Class.ModifiedString
...signature:ModifiedString<RnaString, ModView< FunctorComplement<Rna> > >
..see:Shortcut.RnaString
..see:Class.ModifiedString
..see:Spec.ModView
..see:Class.FunctorComplement
*/

/**
.Shortcut.Rna5StringComplement:
..cat:Modifier
..summary:Modifier for the complement of a @Shortcut.Rna5String@.
..signature:Rna5StringComplement
..shortcutfor:Class.ModifiedString
...signature:ModifiedString<Rna5String, ModView< FunctorComplement<Rna5> > >
..see:Shortcut.Rna5String
..see:Class.ModifiedString
..see:Spec.ModView
..see:Class.FunctorComplement
*/


/**
.Shortcut.DnaStringReverse:
..cat:Modifier
..summary:Modifier for the reverse of a @Shortcut.DnaString@.
..signature:DnaStringReverse
..shortcutfor:Class.ModifiedString
...signature:ModifiedString<DnaString, ModReverse>
..see:Shortcut.DnaString
..see:Class.ModifiedString
..see:Spec.ModReverse
*/

/**
.Shortcut.Dna5StringReverse:
..cat:Modifier
..summary:Modifier for the reverse of a @Shortcut.Dna5String@.
..signature:Dna5StringReverse
..shortcutfor:Class.ModifiedString
...signature:ModifiedString<Dna5String, ModReverse>
..see:Shortcut.Dna5String
..see:Class.ModifiedString
..see:Spec.ModReverse
*/

/**
.Shortcut.RnaStringReverse:
..cat:Modifier
..summary:Modifier for the reverse of a @Shortcut.RnaString@.
..signature:RnaStringReverse
..shortcutfor:Class.ModifiedString
...signature:ModifiedString<RnaString, ModReverse>
..see:Shortcut.RnaString
..see:Class.ModifiedString
..see:Spec.ModReverse
*/

/**
.Shortcut.Rna5StringReverse:
..cat:Modifier
..summary:Modifier for the reverse of a @Shortcut.Rna5String@.
..signature:Rna5StringReverse
..shortcutfor:Class.ModifiedString
...signature:ModifiedString<Rna5String, ModReverse>
..see:Shortcut.Rna5String
..see:Class.ModifiedString
..see:Spec.ModReverse
*/


/**
.Shortcut.DnaStringReverseComplement:
..cat:Modifier
..summary:Modifier for the reverse complement of a @Shortcut.DnaString@.
..signature:DnaStringReverseComplement
..shortcutfor:Class.ModifiedString
...signature:ModifiedString<ModifiedString<DnaString, ModView< FunctorComplement<Dna> > >, ModReverse>
..see:Shortcut.DnaString
..see:Class.ModifiedString
..see:Spec.ModReverse
..see:Spec.ModView
..see:Class.FunctorComplement
*/

/**
.Shortcut.Dna5StringReverseComplement:
..cat:Modifier
..summary:Modifier for the reverse complement of a @Shortcut.Dna5String@.
..signature:Dna5StringReverseComplement
..shortcutfor:Class.ModifiedString
...signature:ModifiedString<ModifiedString<Dna5String, ModView< FunctorComplement<Dna> > >, ModReverse>
..see:Shortcut.Dna5String
..see:Class.ModifiedString
..see:Spec.ModReverse
..see:Spec.ModView
..see:Class.FunctorComplement
*/

/**
.Shortcut.RnaStringReverseComplement:
..cat:Modifier
..summary:Modifier for the reverse complement of a @Shortcut.RnaString@.
..signature:RnaStringReverseComplement
..shortcutfor:Class.ModifiedString
...signature:ModifiedString<ModifiedString<RnaString, ModView< FunctorComplement<Rna> > >, ModReverse>
..see:Shortcut.RnaString
..see:Class.ModifiedString
..see:Spec.ModReverse
..see:Spec.ModView
..see:Class.FunctorComplement
*/

/**
.Shortcut.Rna5StringReverseComplement:
..cat:Modifier
..summary:Modifier for the reverse complement of a @Shortcut.Rna5String@.
..signature:Rna5StringReverseComplement
..shortcutfor:Class.ModifiedString
...signature:ModifiedString<ModifiedString<Rna5String, ModView< FunctorComplement<Rna> > >, ModReverse>
..see:Shortcut.Rna5String
..see:Class.ModifiedString
..see:Spec.ModReverse
..see:Spec.ModView
..see:Class.FunctorComplement
*/

}

#endif
