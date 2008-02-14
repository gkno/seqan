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

#ifndef SEQAN_HEADER_GRAPH_UTILITY_ALPHABETS_H
#define SEQAN_HEADER_GRAPH_UTILITY_ALPHABETS_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// RNA Alphabet
//////////////////////////////////////////////////////////////////////////////

template <typename T = void>
struct _Translate_Table_Rna5_2_Ascii
{
	static char const VALUE[5];
};
template <typename T>
char const _Translate_Table_Rna5_2_Ascii<T>::VALUE[5] = {'A', 'C', 'G', 'U', 'N'}; 

//////////////////////////////////////////////////////////////////////////////

template <typename T = void>
struct _Translate_Table_Byte_2_Rna5
{
	static char const VALUE[256];
};
template <typename T>
char const _Translate_Table_Byte_2_Rna5<T>::VALUE[256] = 
{
	0,   1,   2,   3,   4,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //0
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //1
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //2
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //3
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //4
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //5
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //6
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //7
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //8
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //9
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //10
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //11
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //12
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //13
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //14
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0  //15
};

//////////////////////////////////////////////////////////////////////////////


template <typename T = void>
struct _Translate_Table_Ascii_2_Rna5
{
	static char const VALUE[256];
};
template <typename T>
char const _Translate_Table_Ascii_2_Rna5<T>::VALUE[256] = 
{
	4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4, //0
	4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4, //1
	4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4, //2
	4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4, //3

	4,   0,   4,   1,   4,   4,   4,   2,   4,   4,   4,   4,   4,   4,   4,   4, //4
//	 ,   A,   B,   C,   D,   E,   D,   G,   H,   I,   J,   K,   L,   M,   N,   O,

	4,   4,   4,   4,   4,   3,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4, //5
//	P,   Q,   R,   S,   T,   U,   V,   W,   X,   Y,   Z,    ,    ,    ,    ,    

	4,   0,   4,   1,   4,   4,   4,   2,   4,   4,   4,   4,   4,   4,   4,   4, //6
//   ,   a,   b,   c,   d,   e,   f,   g,   h,   i,   j,   k,   l,   m,   n,   o,

	4,   4,   4,   4,   4,   3,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4, //7
//  p,   q,   r,   s,   t,   u,   v,   w,   x,   y,   z,    ,    ,    ,    ,   

	4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4, //8
	4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4, //9
	4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4, //10
	4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4, //11
	4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4, //12
	4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4, //13
	4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4, //14
	4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4  //15
};

//////////////////////////////////////////////////////////////////////////////

/**
.Spec.Rna5:
..cat:Alphabets
..summary:Rna5 alphabet.
..general:Class.SimpleType
..signature:Rna5
..remarks:
...text:The @Metafunction.ValueSize@ of $Rna5$ is 5. 
...text:Objects of type $Rna5$ cannot be converted into other types.
..see:Metafunction.ValueSize
*/
struct _Rna5 {};
typedef SimpleType<unsigned char, _Rna5> Rna5;

template <> struct ValueSize< Rna5 > { enum { VALUE = 5 }; };
template <> struct BitsPerValue< Rna5 > { enum { VALUE = 3 }; };


//////////////////////////////////////////////////////////////////////////////
//Rna5 assignment
//////////////////////////////////////////////////////////////////////////////

inline void 
assign(Ascii& target,
	   Rna5 const & source)
{
	SEQAN_CHECKPOINT
	target = _Translate_Table_Rna5_2_Ascii<>::VALUE[source.value];
}

//////////////////////////////////////////////////////////////////////////////

template <>
struct CompareType<Rna5, Byte> { typedef Rna5 Type; };
inline void assign(Rna5 & target, Byte c_source)
{
	SEQAN_CHECKPOINT
	target.value = _Translate_Table_Byte_2_Rna5<>::VALUE[c_source];
}

//////////////////////////////////////////////////////////////////////////////

template <>
struct CompareType<Rna5, Ascii> { typedef Rna5 Type; };
inline void assign(Rna5 & target, Ascii c_source)
{
	SEQAN_CHECKPOINT
	target.value = _Translate_Table_Ascii_2_Rna5<>::VALUE[(unsigned char)c_source];
}

//////////////////////////////////////////////////////////////////////////////


template <>
struct CompareType<Rna5, Unicode> { typedef Rna5 Type; };
inline void assign(Rna5 & target, Unicode c_source)
{
	SEQAN_CHECKPOINT
	target.value = _Translate_Table_Ascii_2_Rna5<>::VALUE[(unsigned char) c_source];
}








//////////////////////////////////////////////////////////////////////////////
// Compressed amino acid alphabets
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////
// Dayhoff
//////////////////////////////////////////////////////////////////////////////

template <typename T = void>
struct _Translate_Table_AA_2_AAGroupsDayhoff
{
	static char const VALUE[24];
};
template <typename T>
char const _Translate_Table_AA_2_AAGroupsDayhoff<T>::VALUE[24] = 
{
	0, // A 0 Ala Alanine                 
	3, // R 1 Arg Arginine                
	2, // N 2 Asn Asparagine              
	2, // D 3 Asp Aspartic Acid           
	5, // C 4 Cys Cystine                 
	2, // Q 5 Gln Glutamine               
	2, // E 6 Glu Glutamic Acid           
	0, // G 7 Gly Glycine                 
	3, // H 8 His Histidine               
	1, // I 9 Ile Isoleucine              
	1, // L 10 Leu Leucine                 
	3, // K 11 Lys Lysine                  
	1, // M 12 Met Methionine              
	4, // F 13 Phe Phenylalanine           
	0, // P 14 Pro Proline                 
	0, // S 15 Ser Serine                  
	0, // T 16 Thr Threonine               
	4, // W 17 Trp Tryptophan              
	4, // Y 18 Tyr Tyrosine                
	1, // V 19 Val Valine                  
	2, // B 20 Aspartic Acid, Asparagine   
	2, // Z 21 Glutamic Acid, Glutamine    
	6, // X 22 Unknown                     
	6  // * 23 Terminator                  
};	              

//////////////////////////////////////////////////////////////////////////////

/**
.Spec.AAGroupsDayhoff:
..cat:Alphabets
..summary:Compressed amino acid alphabet (Dayhoff).
..general:Class.SimpleType
..signature:AAGroupsDayhoff
..remarks:
...text:The @Metafunction.ValueSize@ of $AAGroupsDayhoff$ is 7. 
The groups are defined in the following way: "agjopst"=0, "ilmv"=1, "bdenqz"=2, "hkr"=3, "fwy"=4, "c"=5, all others = 6
...text:Objects of type $AAGroupsDayhoff$ cannot be converted into other types.
..see:Metafunction.ValueSize
*/
struct _AAGroupsDayhoff {};
typedef SimpleType<unsigned char, _AAGroupsDayhoff> AAGroupsDayhoff;

template <> struct ValueSize< AAGroupsDayhoff > { enum { VALUE = 7 }; };
template <> struct BitsPerValue< AAGroupsDayhoff > { enum { VALUE = 4 }; };

//////////////////////////////////////////////////////////////////////////////
//AAGroupsDayhoff assignment

template <>
struct CompareType<AAGroupsDayhoff, AminoAcid> { typedef AAGroupsDayhoff Type; };
inline void assign(AAGroupsDayhoff & target, AminoAcid const c_source)
{
SEQAN_CHECKPOINT
	target.value = _Translate_Table_AA_2_AAGroupsDayhoff<>::VALUE[(unsigned char) c_source];
}

//////////////////////////////////////////////////////////////////////////////

template <>
struct CompareType<AAGroupsDayhoff, Byte> { typedef AAGroupsDayhoff Type; };
inline void assign(AAGroupsDayhoff & target, Byte const c_source)
{
SEQAN_CHECKPOINT
	target.value = _Translate_Table_AA_2_AAGroupsDayhoff<>::VALUE[(unsigned char) _Translate_Table_Byte_2_AA<>::VALUE[c_source] ];
}

//////////////////////////////////////////////////////////////////////////////

template <>
struct CompareType<AAGroupsDayhoff, Ascii> { typedef AAGroupsDayhoff Type; };
inline void assign(AAGroupsDayhoff & target, Ascii const c_source)
{
SEQAN_CHECKPOINT
	target.value = _Translate_Table_AA_2_AAGroupsDayhoff<>::VALUE[(unsigned char) _Translate_Table_Ascii_2_AA<>::VALUE[(unsigned char) c_source] ];
}

//////////////////////////////////////////////////////////////////////////////

template <>
struct CompareType<AAGroupsDayhoff, Unicode> { typedef AAGroupsDayhoff Type; };
inline void assign(AAGroupsDayhoff & target, Unicode const c_source)
{
SEQAN_CHECKPOINT
	target.value = _Translate_Table_AA_2_AAGroupsDayhoff<>::VALUE[(unsigned char) _Translate_Table_Ascii_2_AA<>::VALUE[(unsigned char) c_source] ];
}

//////////////////////////////////////////////////////////////////////////////







//////////////////////////////////////////////////////////////////////////////
// SeB6
//////////////////////////////////////////////////////////////////////////////

template <typename T = void>
struct _Translate_Table_AA_2_AAGroupsSeB6
{
	static char const VALUE[24];
};
template <typename T>
char const _Translate_Table_AA_2_AAGroupsSeB6<T>::VALUE[24] = 
{
	0, // A 0 Ala Alanine                 
	2, // R 1 Arg Arginine                
	2, // N 2 Asn Asparagine              
	2, // D 3 Asp Aspartic Acid           
	1, // C 4 Cys Cystine                 
	2, // Q 5 Gln Glutamine               
	2, // E 6 Glu Glutamic Acid           
	4, // G 7 Gly Glycine                 
	2, // H 8 His Histidine               
	5, // I 9 Ile Isoleucine              
	5, // L 10 Leu Leucine                 
	2, // K 11 Lys Lysine                  
	5, // M 12 Met Methionine              
	3, // F 13 Phe Phenylalanine           
	1, // P 14 Pro Proline                 
	0, // S 15 Ser Serine                  
	0, // T 16 Thr Threonine               
	3, // W 17 Trp Tryptophan              
	3, // Y 18 Tyr Tyrosine                
	5, // V 19 Val Valine                  
	2, // B 20 Aspartic Acid, Asparagine   
	2, // Z 21 Glutamic Acid, Glutamine    
	6, // X 22 Unknown                     
	6  // * 23 Terminator                  
};	              

//////////////////////////////////////////////////////////////////////////////

/**
.Spec.AAGroupsSeB6:
..cat:Alphabets
..summary:Compressed amino acid alphabet (SeB6).
..general:Class.SimpleType
..signature:AAGroupsSeB6
..remarks:
...text:The @Metafunction.ValueSize@ of $AAGroupsSeB6$ is 7. 
...text:Objects of type $AAGroupsSeB6$ cannot be converted into other types.
..see:Metafunction.ValueSize
*/
struct _AAGroupsSeB6 {};
typedef SimpleType<unsigned char, _AAGroupsSeB6> AAGroupsSeB6;

template <> struct ValueSize< AAGroupsSeB6 > { enum { VALUE = 7 }; };
template <> struct BitsPerValue< AAGroupsSeB6 > { enum { VALUE = 4 }; };

//////////////////////////////////////////////////////////////////////////////
//AAGroupsSeB6 assignment

template <>
struct CompareType<AAGroupsSeB6, AminoAcid> { typedef AAGroupsSeB6 Type; };
inline void assign(AAGroupsSeB6 & target, AminoAcid const c_source)
{
SEQAN_CHECKPOINT
	target.value = _Translate_Table_AA_2_AAGroupsSeB6<>::VALUE[(unsigned char) c_source];
}

//////////////////////////////////////////////////////////////////////////////

template <>
struct CompareType<AAGroupsSeB6, Byte> { typedef AAGroupsSeB6 Type; };
inline void assign(AAGroupsSeB6 & target, Byte const c_source)
{
SEQAN_CHECKPOINT
	target.value = _Translate_Table_AA_2_AAGroupsSeB6<>::VALUE[(unsigned char) _Translate_Table_Byte_2_AA<>::VALUE[c_source] ];
}

//////////////////////////////////////////////////////////////////////////////

template <>
struct CompareType<AAGroupsSeB6, Ascii> { typedef AAGroupsSeB6 Type; };
inline void assign(AAGroupsSeB6 & target, Ascii const c_source)
{
SEQAN_CHECKPOINT
	target.value = _Translate_Table_AA_2_AAGroupsSeB6<>::VALUE[(unsigned char) _Translate_Table_Ascii_2_AA<>::VALUE[(unsigned char) c_source] ];
}

//////////////////////////////////////////////////////////////////////////////

template <>
struct CompareType<AAGroupsSeB6, Unicode> { typedef AAGroupsSeB6 Type; };
inline void assign(AAGroupsSeB6 & target, Unicode const c_source)
{
SEQAN_CHECKPOINT
	target.value = _Translate_Table_AA_2_AAGroupsSeB6<>::VALUE[(unsigned char) _Translate_Table_Ascii_2_AA<>::VALUE[(unsigned char) c_source] ];
}

//////////////////////////////////////////////////////////////////////////////






//////////////////////////////////////////////////////////////////////////////
// SeB8
//////////////////////////////////////////////////////////////////////////////

template <typename T = void>
struct _Translate_Table_AA_2_AAGroupsSeB8
{
	static char const VALUE[24];
};
template <typename T>
char const _Translate_Table_AA_2_AAGroupsSeB8<T>::VALUE[24] = 
{
	0, // A 0 Ala Alanine                 
	3, // R 1 Arg Arginine                
	2, // N 2 Asn Asparagine              
	2, // D 3 Asp Aspartic Acid           
	1, // C 4 Cys Cystine                 
	3, // Q 5 Gln Glutamine               
	3, // E 6 Glu Glutamic Acid           
	5, // G 7 Gly Glycine                 
	2, // H 8 His Histidine               
	6, // I 9 Ile Isoleucine              
	6, // L 10 Leu Leucine                 
	3, // K 11 Lys Lysine                  
	6, // M 12 Met Methionine              
	4, // F 13 Phe Phenylalanine           
	7, // P 14 Pro Proline                 
	0, // S 15 Ser Serine                  
	0, // T 16 Thr Threonine               
	4, // W 17 Trp Tryptophan              
	4, // Y 18 Tyr Tyrosine                
	6, // V 19 Val Valine                  
	2, // B 20 Aspartic Acid, Asparagine   
	3, // Z 21 Glutamic Acid, Glutamine    
	8, // X 22 Unknown                     
	8  // * 23 Terminator                  
};	              

//////////////////////////////////////////////////////////////////////////////

/**
.Spec.AAGroupsSeB8:
..cat:Alphabets
..summary:Compressed amino acid alphabet (SeB8).
..general:Class.SimpleType
..signature:AAGroupsSeB8
..remarks:
...text:The @Metafunction.ValueSize@ of $AAGroupsSeB8$ is 9. 
...text:Objects of type $AAGroupsSeB8$ cannot be converted into other types.
..see:Metafunction.ValueSize
*/
struct _AAGroupsSeB8 {};
typedef SimpleType<unsigned char, _AAGroupsSeB8> AAGroupsSeB8;

template <> struct ValueSize< AAGroupsSeB8 > { enum { VALUE = 9 }; };
template <> struct BitsPerValue< AAGroupsSeB8 > { enum { VALUE = 4 }; };

//////////////////////////////////////////////////////////////////////////////
//AAGroupsSeB8 assignment

template <>
struct CompareType<AAGroupsSeB8, AminoAcid> { typedef AAGroupsSeB8 Type; };
inline void assign(AAGroupsSeB8 & target, AminoAcid const c_source)
{
SEQAN_CHECKPOINT
	target.value = _Translate_Table_AA_2_AAGroupsSeB8<>::VALUE[(unsigned char) c_source];
}

//////////////////////////////////////////////////////////////////////////////

template <>
struct CompareType<AAGroupsSeB8, Byte> { typedef AAGroupsSeB8 Type; };
inline void assign(AAGroupsSeB8 & target, Byte const c_source)
{
SEQAN_CHECKPOINT
	target.value = _Translate_Table_AA_2_AAGroupsSeB8<>::VALUE[(unsigned char) _Translate_Table_Byte_2_AA<>::VALUE[c_source] ];
}

//////////////////////////////////////////////////////////////////////////////

template <>
struct CompareType<AAGroupsSeB8, Ascii> { typedef AAGroupsSeB8 Type; };
inline void assign(AAGroupsSeB8 & target, Ascii const c_source)
{
SEQAN_CHECKPOINT
	target.value = _Translate_Table_AA_2_AAGroupsSeB8<>::VALUE[(unsigned char) _Translate_Table_Ascii_2_AA<>::VALUE[(unsigned char) c_source] ];
}

//////////////////////////////////////////////////////////////////////////////

template <>
struct CompareType<AAGroupsSeB8, Unicode> { typedef AAGroupsSeB8 Type; };
inline void assign(AAGroupsSeB8 & target, Unicode const c_source)
{
SEQAN_CHECKPOINT
	target.value = _Translate_Table_AA_2_AAGroupsSeB8<>::VALUE[(unsigned char) _Translate_Table_Ascii_2_AA<>::VALUE[(unsigned char) c_source] ];
}

//////////////////////////////////////////////////////////////////////////////






//////////////////////////////////////////////////////////////////////////////
// Murphy
//////////////////////////////////////////////////////////////////////////////

template <typename T = void>
struct _Translate_Table_AA_2_AAGroupsMurphy
{
	static char const VALUE[24];
};
template <typename T>
char const _Translate_Table_AA_2_AAGroupsMurphy<T>::VALUE[24] = 
{
	0, // A 0 Ala Alanine                 
	7, // R 1 Arg Arginine                
	2, // N 2 Asn Asparagine              
	2, // D 3 Asp Aspartic Acid           
	1, // C 4 Cys Cystine                 
	2, // Q 5 Gln Glutamine               
	2, // E 6 Glu Glutamic Acid           
	4, // G 7 Gly Glycine                 
	5, // H 8 His Histidine               
	6, // I 9 Ile Isoleucine              
	6, // L 10 Leu Leucine                 
	7, // K 11 Lys Lysine                  
	6, // M 12 Met Methionine              
	3, // F 13 Phe Phenylalanine           
	8, // P 14 Pro Proline                 
	9, // S 15 Ser Serine                  
	9, // T 16 Thr Threonine               
	3, // W 17 Trp Tryptophan              
	3, // Y 18 Tyr Tyrosine                
	6, // V 19 Val Valine                  
	2, // B 20 Aspartic Acid, Asparagine   
	2, // Z 21 Glutamic Acid, Glutamine    
	10, // X 22 Unknown                     
	10  // * 23 Terminator                  
};	              

//////////////////////////////////////////////////////////////////////////////

/**
.Spec.AAGroupsMurphy:
..cat:Alphabets
..summary:Compressed amino acid alphabet (Murphy).
..general:Class.SimpleType
..signature:AAGroupsMurphy
..remarks:
...text:The @Metafunction.ValueSize@ of $AAGroupsMurphy$ is 11. 
...text:Objects of type $AAGroupsMurphy$ cannot be converted into other types.
..see:Metafunction.ValueSize
*/
struct _AAGroupsMurphy {};
typedef SimpleType<unsigned char, _AAGroupsMurphy> AAGroupsMurphy;

template <> struct ValueSize< AAGroupsMurphy > { enum { VALUE = 11 }; };
template <> struct BitsPerValue< AAGroupsMurphy > { enum { VALUE = 4 }; };

//////////////////////////////////////////////////////////////////////////////
//AAGroupsMurphy assignment

template <>
struct CompareType<AAGroupsMurphy, AminoAcid> { typedef AAGroupsMurphy Type; };
inline void assign(AAGroupsMurphy & target, AminoAcid const c_source)
{
SEQAN_CHECKPOINT
	target.value = _Translate_Table_AA_2_AAGroupsMurphy<>::VALUE[(unsigned char) c_source];
}

//////////////////////////////////////////////////////////////////////////////

template <>
struct CompareType<AAGroupsMurphy, Byte> { typedef AAGroupsMurphy Type; };
inline void assign(AAGroupsMurphy & target, Byte const c_source)
{
SEQAN_CHECKPOINT
	target.value = _Translate_Table_AA_2_AAGroupsMurphy<>::VALUE[(unsigned char) _Translate_Table_Byte_2_AA<>::VALUE[c_source] ];
}

//////////////////////////////////////////////////////////////////////////////

template <>
struct CompareType<AAGroupsMurphy, Ascii> { typedef AAGroupsMurphy Type; };
inline void assign(AAGroupsMurphy & target, Ascii const c_source)
{
SEQAN_CHECKPOINT
	target.value = _Translate_Table_AA_2_AAGroupsMurphy<>::VALUE[(unsigned char) _Translate_Table_Ascii_2_AA<>::VALUE[(unsigned char) c_source] ];
}

//////////////////////////////////////////////////////////////////////////////

template <>
struct CompareType<AAGroupsMurphy, Unicode> { typedef AAGroupsMurphy Type; };
inline void assign(AAGroupsMurphy & target, Unicode const c_source)
{
SEQAN_CHECKPOINT
	target.value = _Translate_Table_AA_2_AAGroupsMurphy<>::VALUE[(unsigned char) _Translate_Table_Ascii_2_AA<>::VALUE[(unsigned char) c_source] ];
}

//////////////////////////////////////////////////////////////////////////////





//////////////////////////////////////////////////////////////////////////////
// SolisG10
//////////////////////////////////////////////////////////////////////////////

template <typename T = void>
struct _Translate_Table_AA_2_AAGroupsSolisG10
{
	static char const VALUE[24];
};
template <typename T>
char const _Translate_Table_AA_2_AAGroupsSolisG10<T>::VALUE[24] = 
{
	0, // A 0 Ala Alanine                 
	0, // R 1 Arg Arginine                
	5, // N 2 Asn Asparagine              
	2, // D 3 Asp Aspartic Acid           
	1, // C 4 Cys Cystine                 
	0, // Q 5 Gln Glutamine               
	0, // E 6 Glu Glutamic Acid           
	3, // G 7 Gly Glycine                 
	4, // H 8 His Histidine               
	0, // I 9 Ile Isoleucine              
	0, // L 10 Leu Leucine                 
	0, // K 11 Lys Lysine                  
	0, // M 12 Met Methionine              
	0, // F 13 Phe Phenylalanine           
	6, // P 14 Pro Proline                 
	7, // S 15 Ser Serine                  
	8, // T 16 Thr Threonine               
	0, // W 17 Trp Tryptophan              
	9, // Y 18 Tyr Tyrosine                
	0, // V 19 Val Valine                  
	5, // B 20 Aspartic Acid, Asparagine   
	0, // Z 21 Glutamic Acid, Glutamine    
	10, // X 22 Unknown                     
	10  // * 23 Terminator                  
};	              

//////////////////////////////////////////////////////////////////////////////

/**
.Spec.AAGroupsSolisG10:
..cat:Alphabets
..summary:Compressed amino acid alphabet (SolisG10).
..general:Class.SimpleType
..signature:AAGroupsSolisG10
..remarks:
...text:The @Metafunction.ValueSize@ of $AAGroupsSolisG10$ is 11. 
...text:Objects of type $AAGroupsSolisG10$ cannot be converted into other types.
..see:Metafunction.ValueSize
*/
struct _AAGroupsSolisG10 {};
typedef SimpleType<unsigned char, _AAGroupsSolisG10> AAGroupsSolisG10;

template <> struct ValueSize< AAGroupsSolisG10 > { enum { VALUE = 11 }; };
template <> struct BitsPerValue< AAGroupsSolisG10 > { enum { VALUE = 4 }; };

//////////////////////////////////////////////////////////////////////////////
//AAGroupsSolisG10 assignment

template <>
struct CompareType<AAGroupsSolisG10, AminoAcid> { typedef AAGroupsSolisG10 Type; };
inline void assign(AAGroupsSolisG10 & target, AminoAcid const c_source)
{
SEQAN_CHECKPOINT
	target.value = _Translate_Table_AA_2_AAGroupsSolisG10<>::VALUE[(unsigned char) c_source];
}

//////////////////////////////////////////////////////////////////////////////

template <>
struct CompareType<AAGroupsSolisG10, Byte> { typedef AAGroupsSolisG10 Type; };
inline void assign(AAGroupsSolisG10 & target, Byte const c_source)
{
SEQAN_CHECKPOINT
	target.value = _Translate_Table_AA_2_AAGroupsSolisG10<>::VALUE[(unsigned char) _Translate_Table_Byte_2_AA<>::VALUE[c_source] ];
}

//////////////////////////////////////////////////////////////////////////////

template <>
struct CompareType<AAGroupsSolisG10, Ascii> { typedef AAGroupsSolisG10 Type; };
inline void assign(AAGroupsSolisG10 & target, Ascii const c_source)
{
SEQAN_CHECKPOINT
	target.value = _Translate_Table_AA_2_AAGroupsSolisG10<>::VALUE[(unsigned char) _Translate_Table_Ascii_2_AA<>::VALUE[(unsigned char) c_source] ];
}

//////////////////////////////////////////////////////////////////////////////

template <>
struct CompareType<AAGroupsSolisG10, Unicode> { typedef AAGroupsSolisG10 Type; };
inline void assign(AAGroupsSolisG10 & target, Unicode const c_source)
{
SEQAN_CHECKPOINT
	target.value = _Translate_Table_AA_2_AAGroupsSolisG10<>::VALUE[(unsigned char) _Translate_Table_Ascii_2_AA<>::VALUE[(unsigned char) c_source] ];
}

//////////////////////////////////////////////////////////////////////////////






//////////////////////////////////////////////////////////////////////////////
// SolisD10
//////////////////////////////////////////////////////////////////////////////

template <typename T = void>
struct _Translate_Table_AA_2_AAGroupsSolisD10
{
	static char const VALUE[24];
};
template <typename T>
char const _Translate_Table_AA_2_AAGroupsSolisD10<T>::VALUE[24] = 
{
	0, // A 0 Ala Alanine                 
	3, // R 1 Arg Arginine                
	2, // N 2 Asn Asparagine              
	2, // D 3 Asp Aspartic Acid           
	1, // C 4 Cys Cystine                 
	3, // Q 5 Gln Glutamine               
	3, // E 6 Glu Glutamic Acid           
	5, // G 7 Gly Glycine                 
	6, // H 8 His Histidine               
	7, // I 9 Ile Isoleucine              
	8, // L 10 Leu Leucine                 
	3, // K 11 Lys Lysine                  
	0, // M 12 Met Methionine              
	4, // F 13 Phe Phenylalanine           
	5, // P 14 Pro Proline                 
	2, // S 15 Ser Serine                  
	6, // T 16 Thr Threonine               
	9, // W 17 Trp Tryptophan              
	8, // Y 18 Tyr Tyrosine                
	7, // V 19 Val Valine                  
	2, // B 20 Aspartic Acid, Asparagine   
	3, // Z 21 Glutamic Acid, Glutamine    
	10, // X 22 Unknown                     
	10  // * 23 Terminator                  
};	              

//////////////////////////////////////////////////////////////////////////////

/**
.Spec.AAGroupsSolisD10:
..cat:Alphabets
..summary:Compressed amino acid alphabet (SolisD10).
..general:Class.SimpleType
..signature:AAGroupsSolisD10
..remarks:
...text:The @Metafunction.ValueSize@ of $AAGroupsSolisD10$ is 11. 
...text:Objects of type $AAGroupsSolisD10$ cannot be converted into other types.
..see:Metafunction.ValueSize
*/
struct _AAGroupsSolisD10 {};
typedef SimpleType<unsigned char, _AAGroupsSolisD10> AAGroupsSolisD10;

template <> struct ValueSize< AAGroupsSolisD10 > { enum { VALUE = 11 }; };
template <> struct BitsPerValue< AAGroupsSolisD10 > { enum { VALUE = 4 }; };

//////////////////////////////////////////////////////////////////////////////
//AAGroupsSolisD10 assignment

template <>
struct CompareType<AAGroupsSolisD10, AminoAcid> { typedef AAGroupsSolisD10 Type; };
inline void assign(AAGroupsSolisD10 & target, AminoAcid const c_source)
{
SEQAN_CHECKPOINT
	target.value = _Translate_Table_AA_2_AAGroupsSolisD10<>::VALUE[(unsigned char) c_source];
}

//////////////////////////////////////////////////////////////////////////////

template <>
struct CompareType<AAGroupsSolisD10, Byte> { typedef AAGroupsSolisD10 Type; };
inline void assign(AAGroupsSolisD10 & target, Byte const c_source)
{
SEQAN_CHECKPOINT
	target.value = _Translate_Table_AA_2_AAGroupsSolisD10<>::VALUE[(unsigned char) _Translate_Table_Byte_2_AA<>::VALUE[c_source] ];
}

//////////////////////////////////////////////////////////////////////////////

template <>
struct CompareType<AAGroupsSolisD10, Ascii> { typedef AAGroupsSolisD10 Type; };
inline void assign(AAGroupsSolisD10 & target, Ascii const c_source)
{
SEQAN_CHECKPOINT
	target.value = _Translate_Table_AA_2_AAGroupsSolisD10<>::VALUE[(unsigned char) _Translate_Table_Ascii_2_AA<>::VALUE[(unsigned char) c_source] ];
}

//////////////////////////////////////////////////////////////////////////////

template <>
struct CompareType<AAGroupsSolisD10, Unicode> { typedef AAGroupsSolisD10 Type; };
inline void assign(AAGroupsSolisD10 & target, Unicode const c_source)
{
SEQAN_CHECKPOINT
	target.value = _Translate_Table_AA_2_AAGroupsSolisD10<>::VALUE[(unsigned char) _Translate_Table_Ascii_2_AA<>::VALUE[(unsigned char) c_source] ];
}

//////////////////////////////////////////////////////////////////////////////





//////////////////////////////////////////////////////////////////////////////
// LiB10
//////////////////////////////////////////////////////////////////////////////

template <typename T = void>
struct _Translate_Table_AA_2_AAGroupsLiB10
{
	static char const VALUE[24];
};
template <typename T>
char const _Translate_Table_AA_2_AAGroupsLiB10<T>::VALUE[24] = 
{
	0, // A 0 Ala Alanine                 
	7, // R 1 Arg Arginine                
	5, // N 2 Asn Asparagine              
	2, // D 3 Asp Aspartic Acid           
	1, // C 4 Cys Cystine                 
	2, // Q 5 Gln Glutamine               
	2, // E 6 Glu Glutamic Acid           
	4, // G 7 Gly Glycine                 
	5, // H 8 His Histidine               
	6, // I 9 Ile Isoleucine              
	8, // L 10 Leu Leucine                 
	7, // K 11 Lys Lysine                  
	8, // M 12 Met Methionine              
	3, // F 13 Phe Phenylalanine           
	9, // P 14 Pro Proline                 
	0, // S 15 Ser Serine                  
	0, // T 16 Thr Threonine               
	3, // W 17 Trp Tryptophan              
	3, // Y 18 Tyr Tyrosine                
	6, // V 19 Val Valine                  
	5, // B 20 Aspartic Acid, Asparagine   
	2, // Z 21 Glutamic Acid, Glutamine    
	10, // X 22 Unknown                     
	10  // * 23 Terminator                  
};	              

//////////////////////////////////////////////////////////////////////////////

/**
.Spec.AAGroupsLiB10:
..cat:Alphabets
..summary:Compressed amino acid alphabet (LiB10).
..general:Class.SimpleType
..signature:AAGroupsLiB10
..remarks:
...text:The @Metafunction.ValueSize@ of $AAGroupsLiB10$ is 11. 
...text:Objects of type $AAGroupsLiB10$ cannot be converted into other types.
..see:Metafunction.ValueSize
*/
struct _AAGroupsLiB10 {};
typedef SimpleType<unsigned char, _AAGroupsLiB10> AAGroupsLiB10;

template <> struct ValueSize< AAGroupsLiB10 > { enum { VALUE = 11 }; };
template <> struct BitsPerValue< AAGroupsLiB10 > { enum { VALUE = 4 }; };

//////////////////////////////////////////////////////////////////////////////
//AAGroupsLiB10 assignment

template <>
struct CompareType<AAGroupsLiB10, AminoAcid> { typedef AAGroupsLiB10 Type; };
inline void assign(AAGroupsLiB10 & target, AminoAcid const c_source)
{
SEQAN_CHECKPOINT
	target.value = _Translate_Table_AA_2_AAGroupsLiB10<>::VALUE[(unsigned char) c_source];
}

//////////////////////////////////////////////////////////////////////////////

template <>
struct CompareType<AAGroupsLiB10, Byte> { typedef AAGroupsLiB10 Type; };
inline void assign(AAGroupsLiB10 & target, Byte const c_source)
{
SEQAN_CHECKPOINT
	target.value = _Translate_Table_AA_2_AAGroupsLiB10<>::VALUE[(unsigned char) _Translate_Table_Byte_2_AA<>::VALUE[c_source] ];
}

//////////////////////////////////////////////////////////////////////////////

template <>
struct CompareType<AAGroupsLiB10, Ascii> { typedef AAGroupsLiB10 Type; };
inline void assign(AAGroupsLiB10 & target, Ascii const c_source)
{
SEQAN_CHECKPOINT
	target.value = _Translate_Table_AA_2_AAGroupsLiB10<>::VALUE[(unsigned char) _Translate_Table_Ascii_2_AA<>::VALUE[(unsigned char) c_source] ];
}

//////////////////////////////////////////////////////////////////////////////

template <>
struct CompareType<AAGroupsLiB10, Unicode> { typedef AAGroupsLiB10 Type; };
inline void assign(AAGroupsLiB10 & target, Unicode const c_source)
{
SEQAN_CHECKPOINT
	target.value = _Translate_Table_AA_2_AAGroupsLiB10<>::VALUE[(unsigned char) _Translate_Table_Ascii_2_AA<>::VALUE[(unsigned char) c_source] ];
}

//////////////////////////////////////////////////////////////////////////////




//////////////////////////////////////////////////////////////////////////////
// LiA10
//////////////////////////////////////////////////////////////////////////////

template <typename T = void>
struct _Translate_Table_AA_2_AAGroupsLiA10
{
	static char const VALUE[24];
};
template <typename T>
char const _Translate_Table_AA_2_AAGroupsLiA10<T>::VALUE[24] = 
{
	0, // A 0 Ala Alanine                 
	6, // R 1 Arg Arginine                
	4, // N 2 Asn Asparagine              
	1, // D 3 Asp Aspartic Acid           
	0, // C 4 Cys Cystine                 
	6, // Q 5 Gln Glutamine               
	1, // E 6 Glu Glutamic Acid           
	3, // G 7 Gly Glycine                 
	4, // H 8 His Histidine               
	5, // I 9 Ile Isoleucine              
	7, // L 10 Leu Leucine                 
	6, // K 11 Lys Lysine                  
	7, // M 12 Met Methionine              
	2, // F 13 Phe Phenylalanine           
	8, // P 14 Pro Proline                 
	9, // S 15 Ser Serine                  
	9, // T 16 Thr Threonine               
	2, // W 17 Trp Tryptophan              
	2, // Y 18 Tyr Tyrosine                
	5, // V 19 Val Valine                  
	4, // B 20 Aspartic Acid, Asparagine   
	6, // Z 21 Glutamic Acid, Glutamine    
	10, // X 22 Unknown                     
	10  // * 23 Terminator                  
};	              

//////////////////////////////////////////////////////////////////////////////

/**
.Spec.AAGroupsLiA10:
..cat:Alphabets
..summary:Compressed amino acid alphabet (LiA10).
..general:Class.SimpleType
..signature:AAGroupsLiA10
..remarks:
...text:The @Metafunction.ValueSize@ of $AAGroupsLiA10$ is 11. 
...text:Objects of type $AAGroupsLiA10$ cannot be converted into other types.
..see:Metafunction.ValueSize
*/
struct _AAGroupsLiA10 {};
typedef SimpleType<unsigned char, _AAGroupsLiA10> AAGroupsLiA10;

template <> struct ValueSize< AAGroupsLiA10 > { enum { VALUE = 11 }; };
template <> struct BitsPerValue< AAGroupsLiA10 > { enum { VALUE = 4 }; };

//////////////////////////////////////////////////////////////////////////////
//AAGroupsLiA10 assignment

template <>
struct CompareType<AAGroupsLiA10, AminoAcid> { typedef AAGroupsLiA10 Type; };
inline void assign(AAGroupsLiA10 & target, AminoAcid const c_source)
{
SEQAN_CHECKPOINT
	target.value = _Translate_Table_AA_2_AAGroupsLiA10<>::VALUE[(unsigned char) c_source];
}

//////////////////////////////////////////////////////////////////////////////

template <>
struct CompareType<AAGroupsLiA10, Byte> { typedef AAGroupsLiA10 Type; };
inline void assign(AAGroupsLiA10 & target, Byte const c_source)
{
SEQAN_CHECKPOINT
	target.value = _Translate_Table_AA_2_AAGroupsLiA10<>::VALUE[(unsigned char) _Translate_Table_Byte_2_AA<>::VALUE[c_source] ];
}

//////////////////////////////////////////////////////////////////////////////

template <>
struct CompareType<AAGroupsLiA10, Ascii> { typedef AAGroupsLiA10 Type; };
inline void assign(AAGroupsLiA10 & target, Ascii const c_source)
{
SEQAN_CHECKPOINT
	target.value = _Translate_Table_AA_2_AAGroupsLiA10<>::VALUE[(unsigned char) _Translate_Table_Ascii_2_AA<>::VALUE[(unsigned char) c_source] ];
}

//////////////////////////////////////////////////////////////////////////////

template <>
struct CompareType<AAGroupsLiA10, Unicode> { typedef AAGroupsLiA10 Type; };
inline void assign(AAGroupsLiA10 & target, Unicode const c_source)
{
SEQAN_CHECKPOINT
	target.value = _Translate_Table_AA_2_AAGroupsLiA10<>::VALUE[(unsigned char) _Translate_Table_Ascii_2_AA<>::VALUE[(unsigned char) c_source] ];
}

//////////////////////////////////////////////////////////////////////////////





//////////////////////////////////////////////////////////////////////////////
// SeV10
//////////////////////////////////////////////////////////////////////////////

template <typename T = void>
struct _Translate_Table_AA_2_AAGroupsSeV10
{
	static char const VALUE[24];
};
template <typename T>
char const _Translate_Table_AA_2_AAGroupsSeV10<T>::VALUE[24] = 
{
	0, // A 0 Ala Alanine                 
	7, // R 1 Arg Arginine                
	2, // N 2 Asn Asparagine              
	2, // D 3 Asp Aspartic Acid           
	1, // C 4 Cys Cystine                 
	7, // Q 5 Gln Glutamine               
	2, // E 6 Glu Glutamic Acid           
	4, // G 7 Gly Glycine                 
	5, // H 8 His Histidine               
	6, // I 9 Ile Isoleucine              
	6, // L 10 Leu Leucine                 
	7, // K 11 Lys Lysine                  
	6, // M 12 Met Methionine              
	3, // F 13 Phe Phenylalanine           
	8, // P 14 Pro Proline                 
	0, // S 15 Ser Serine                  
	0, // T 16 Thr Threonine               
	9, // W 17 Trp Tryptophan              
	3, // Y 18 Tyr Tyrosine                
	6, // V 19 Val Valine                  
	2, // B 20 Aspartic Acid, Asparagine   
	7, // Z 21 Glutamic Acid, Glutamine    
	10, // X 22 Unknown                     
	10  // * 23 Terminator                  
};	              

//////////////////////////////////////////////////////////////////////////////

/**
.Spec.AAGroupsSeV10:
..cat:Alphabets
..summary:Compressed amino acid alphabet (SeV10).
..general:Class.SimpleType
..signature:AAGroupsSeV10
..remarks:
...text:The @Metafunction.ValueSize@ of $AAGroupsSeV10$ is 11. 
...text:Objects of type $AAGroupsSeV10$ cannot be converted into other types.
..see:Metafunction.ValueSize
*/
struct _AAGroupsSeV10 {};
typedef SimpleType<unsigned char, _AAGroupsSeV10> AAGroupsSeV10;

template <> struct ValueSize< AAGroupsSeV10 > { enum { VALUE = 11 }; };
template <> struct BitsPerValue< AAGroupsSeV10 > { enum { VALUE = 4 }; };

//////////////////////////////////////////////////////////////////////////////
//AAGroupsSeV10 assignment

template <>
struct CompareType<AAGroupsSeV10, AminoAcid> { typedef AAGroupsSeV10 Type; };
inline void assign(AAGroupsSeV10 & target, AminoAcid const c_source)
{
SEQAN_CHECKPOINT
	target.value = _Translate_Table_AA_2_AAGroupsSeV10<>::VALUE[(unsigned char) c_source];
}

//////////////////////////////////////////////////////////////////////////////

template <>
struct CompareType<AAGroupsSeV10, Byte> { typedef AAGroupsSeV10 Type; };
inline void assign(AAGroupsSeV10 & target, Byte const c_source)
{
SEQAN_CHECKPOINT
	target.value = _Translate_Table_AA_2_AAGroupsSeV10<>::VALUE[(unsigned char) _Translate_Table_Byte_2_AA<>::VALUE[c_source] ];
}

//////////////////////////////////////////////////////////////////////////////

template <>
struct CompareType<AAGroupsSeV10, Ascii> { typedef AAGroupsSeV10 Type; };
inline void assign(AAGroupsSeV10 & target, Ascii const c_source)
{
SEQAN_CHECKPOINT
	target.value = _Translate_Table_AA_2_AAGroupsSeV10<>::VALUE[(unsigned char) _Translate_Table_Ascii_2_AA<>::VALUE[(unsigned char) c_source] ];
}

//////////////////////////////////////////////////////////////////////////////

template <>
struct CompareType<AAGroupsSeV10, Unicode> { typedef AAGroupsSeV10 Type; };
inline void assign(AAGroupsSeV10 & target, Unicode const c_source)
{
SEQAN_CHECKPOINT
	target.value = _Translate_Table_AA_2_AAGroupsSeV10<>::VALUE[(unsigned char) _Translate_Table_Ascii_2_AA<>::VALUE[(unsigned char) c_source] ];
}

//////////////////////////////////////////////////////////////////////////////







//////////////////////////////////////////////////////////////////////////////
// SeB10
//////////////////////////////////////////////////////////////////////////////

template <typename T = void>
struct _Translate_Table_AA_2_AAGroupsSeB10
{
	static char const VALUE[24];
};
template <typename T>
char const _Translate_Table_AA_2_AAGroupsSeB10<T>::VALUE[24] = 
{
	0, // A 0 Ala Alanine                 
	8, // R 1 Arg Arginine                
	2, // N 2 Asn Asparagine              
	2, // D 3 Asp Aspartic Acid           
	1, // C 4 Cys Cystine                 
	3, // Q 5 Gln Glutamine               
	3, // E 6 Glu Glutamic Acid           
	5, // G 7 Gly Glycine                 
	6, // H 8 His Histidine               
	7, // I 9 Ile Isoleucine              
	7, // L 10 Leu Leucine                 
	8, // K 11 Lys Lysine                  
	7, // M 12 Met Methionine              
	4, // F 13 Phe Phenylalanine           
	9, // P 14 Pro Proline                 
	0, // S 15 Ser Serine                  
	0, // T 16 Thr Threonine               
	6, // W 17 Trp Tryptophan              
	4, // Y 18 Tyr Tyrosine                
	7, // V 19 Val Valine                  
	2, // B 20 Aspartic Acid, Asparagine   
	3, // Z 21 Glutamic Acid, Glutamine    
	10, // X 22 Unknown                     
	10  // * 23 Terminator                  
};	              

//////////////////////////////////////////////////////////////////////////////

/**
.Spec.AAGroupsSeB10:
..cat:Alphabets
..summary:Compressed amino acid alphabet (SeB10).
..general:Class.SimpleType
..signature:AAGroupsSeB10
..remarks:
...text:The @Metafunction.ValueSize@ of $AAGroupsSeB10$ is 11. 
...text:Objects of type $AAGroupsSeB10$ cannot be converted into other types.
..see:Metafunction.ValueSize
*/
struct _AAGroupsSeB10 {};
typedef SimpleType<unsigned char, _AAGroupsSeB10> AAGroupsSeB10;

template <> struct ValueSize< AAGroupsSeB10 > { enum { VALUE = 11 }; };
template <> struct BitsPerValue< AAGroupsSeB10 > { enum { VALUE = 4 }; };

//////////////////////////////////////////////////////////////////////////////
//AAGroupsSeB10 assignment

template <>
struct CompareType<AAGroupsSeB10, AminoAcid> { typedef AAGroupsSeB10 Type; };
inline void assign(AAGroupsSeB10 & target, AminoAcid const c_source)
{
SEQAN_CHECKPOINT
	target.value = _Translate_Table_AA_2_AAGroupsSeB10<>::VALUE[(unsigned char) c_source];
}

//////////////////////////////////////////////////////////////////////////////

template <>
struct CompareType<AAGroupsSeB10, Byte> { typedef AAGroupsSeB10 Type; };
inline void assign(AAGroupsSeB10 & target, Byte const c_source)
{
SEQAN_CHECKPOINT
	target.value = _Translate_Table_AA_2_AAGroupsSeB10<>::VALUE[(unsigned char) _Translate_Table_Byte_2_AA<>::VALUE[c_source] ];
}

//////////////////////////////////////////////////////////////////////////////

template <>
struct CompareType<AAGroupsSeB10, Ascii> { typedef AAGroupsSeB10 Type; };
inline void assign(AAGroupsSeB10 & target, Ascii const c_source)
{
SEQAN_CHECKPOINT
	target.value = _Translate_Table_AA_2_AAGroupsSeB10<>::VALUE[(unsigned char) _Translate_Table_Ascii_2_AA<>::VALUE[(unsigned char) c_source] ];
}

//////////////////////////////////////////////////////////////////////////////

template <>
struct CompareType<AAGroupsSeB10, Unicode> { typedef AAGroupsSeB10 Type; };
inline void assign(AAGroupsSeB10 & target, Unicode const c_source)
{
SEQAN_CHECKPOINT
	target.value = _Translate_Table_AA_2_AAGroupsSeB10<>::VALUE[(unsigned char) _Translate_Table_Ascii_2_AA<>::VALUE[(unsigned char) c_source] ];
}

//////////////////////////////////////////////////////////////////////////////





//////////////////////////////////////////////////////////////////////////////
// SeB14
//////////////////////////////////////////////////////////////////////////////

template <typename T = void>
struct _Translate_Table_AA_2_AAGroupsSeB14
{
	static char const VALUE[24];
};
template <typename T>
char const _Translate_Table_AA_2_AAGroupsSeB14<T>::VALUE[24] = 
{
	0, // A 0 Ala Alanine                 
	8, // R 1 Arg Arginine                
	10, // N 2 Asn Asparagine              
	2, // D 3 Asp Aspartic Acid           
	1, // C 4 Cys Cystine                 
	3, // Q 5 Gln Glutamine               
	3, // E 6 Glu Glutamic Acid           
	5, // G 7 Gly Glycine                 
	6, // H 8 His Histidine               
	7, // I 9 Ile Isoleucine              
	9, // L 10 Leu Leucine                 
	8, // K 11 Lys Lysine                  
	9, // M 12 Met Methionine              
	4, // F 13 Phe Phenylalanine           
	11, // P 14 Pro Proline                 
	12, // S 15 Ser Serine                  
	12, // T 16 Thr Threonine               
	13, // W 17 Trp Tryptophan              
	4, // Y 18 Tyr Tyrosine                
	7, // V 19 Val Valine                  
	2, // B 20 Aspartic Acid, Asparagine   
	3, // Z 21 Glutamic Acid, Glutamine    
	14, // X 22 Unknown                     
	14  // * 23 Terminator                  
};	              

//////////////////////////////////////////////////////////////////////////////

/**
.Spec.AAGroupsSeB14:
..cat:Alphabets
..summary:Compressed amino acid alphabet (SeB14).
..general:Class.SimpleType
..signature:AAGroupsSeB14
..remarks:
...text:The @Metafunction.ValueSize@ of $AAGroupsSeB14$ is 15. 
...text:Objects of type $AAGroupsSeB14$ cannot be converted into other types.
..see:Metafunction.ValueSize
*/
struct _AAGroupsSeB14 {};
typedef SimpleType<unsigned char, _AAGroupsSeB14> AAGroupsSeB14;

template <> struct ValueSize< AAGroupsSeB14 > { enum { VALUE = 15 }; };
template <> struct BitsPerValue< AAGroupsSeB14 > { enum { VALUE = 4 }; };

//////////////////////////////////////////////////////////////////////////////
//AAGroupsSeB14 assignment

template <>
struct CompareType<AAGroupsSeB14, AminoAcid> { typedef AAGroupsSeB14 Type; };
inline void assign(AAGroupsSeB14 & target, AminoAcid const c_source)
{
SEQAN_CHECKPOINT
	target.value = _Translate_Table_AA_2_AAGroupsSeB14<>::VALUE[(unsigned char) c_source];
}

//////////////////////////////////////////////////////////////////////////////

template <>
struct CompareType<AAGroupsSeB14, Byte> { typedef AAGroupsSeB14 Type; };
inline void assign(AAGroupsSeB14 & target, Byte const c_source)
{
SEQAN_CHECKPOINT
	target.value = _Translate_Table_AA_2_AAGroupsSeB14<>::VALUE[(unsigned char) _Translate_Table_Byte_2_AA<>::VALUE[c_source] ];
}

//////////////////////////////////////////////////////////////////////////////

template <>
struct CompareType<AAGroupsSeB14, Ascii> { typedef AAGroupsSeB14 Type; };
inline void assign(AAGroupsSeB14 & target, Ascii const c_source)
{
SEQAN_CHECKPOINT
	target.value = _Translate_Table_AA_2_AAGroupsSeB14<>::VALUE[(unsigned char) _Translate_Table_Ascii_2_AA<>::VALUE[(unsigned char) c_source] ];
}

//////////////////////////////////////////////////////////////////////////////

template <>
struct CompareType<AAGroupsSeB14, Unicode> { typedef AAGroupsSeB14 Type; };
inline void assign(AAGroupsSeB14 & target, Unicode const c_source)
{
SEQAN_CHECKPOINT
	target.value = _Translate_Table_AA_2_AAGroupsSeB14<>::VALUE[(unsigned char) _Translate_Table_Ascii_2_AA<>::VALUE[(unsigned char) c_source] ];
}

//////////////////////////////////////////////////////////////////////////////

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
