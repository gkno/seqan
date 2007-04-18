#ifndef SEQAN_HEADER_GRAPH_UTILITY_ALPHABETS_H
#define SEQAN_HEADER_GRAPH_UTILITY_ALPHABETS_H

namespace SEQAN_NAMESPACE_MAIN
{
	
//////////////////////////////////////////////////////////////////////////////
// Compressed amino acid alphabets
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////


template <typename T = void>
struct _Translate_Table_AA_2_AAGroups
{
	static char const VALUE[24];
};
template <typename T>
char const _Translate_Table_AA_2_AAGroups<T>::VALUE[24] = 
{
	0, // 0 Ala Alanine                 
	3, // 1 Arg Arginine                
	2, // 2 Asn Asparagine              
	2, // 3 Asp Aspartic Acid           
	5, // 4 Cys Cystine                 
	2, // 5 Gln Glutamine               
	2, // 6 Glu Glutamic Acid           
	0, // 7 Gly Glycine                 
	3, // 8 His Histidine               
	1, // 9 Ile Isoleucine              
	1, //10 Leu Leucine                 
	3, //11 Lys Lysine                  
	1, //12 Met Methionine              
	4, //13 Phe Phenylalanine           
	0, //14 Pro Proline                 
	0, //15 Ser Serine                  
	0, //16 Thr Threonine               
	4, //17 Trp Tryptophan              
	4, //18 Tyr Tyrosine                
	1, //19 Val Valine                  
	2, //20 Aspartic Acid, Asparagine   
	2, //21 Glutamic Acid, Glutamine    
	6, //22 Unknown                     
	6  //23 Terminator                  
};

//////////////////////////////////////////////////////////////////////////////

/**
.Spec.AAGroups:
..cat:Alphabets
..summary:Alphabet for Amino Acid Groups.
..general:Class.SimpleType
..signature:AAGroups
..remarks:
...text:The @Metafunction.ValueSize@ of $AAGroups$ is 7. 
The groups are defined in the following way: "agjopst"=0, "ilmv"=1, "bdenqz"=2, "hkr"=3, "fwy"=4, "c"=5, all others = 6
...text:Objects of type $AAGroups$ cannot be converted into other types.
...text:$AAGroups$ is typedef for $SimpleType<char,_AAGroups>$, while $_AAGroups$ is an auxilliary specialization tag class.
..see:Metafunction.ValueSize
*/
struct _AAGroups {};
typedef SimpleType<unsigned char,_AAGroups> AAGroups;

template <> struct ValueSize< AAGroups > { enum { VALUE = 7 }; };
template <> struct BitsPerValue< AAGroups > { enum { VALUE = 4 }; };

//////////////////////////////////////////////////////////////////////////////
//AAGroups assignment

template <>
struct CompareType<AAGroups, AminoAcid> { typedef AAGroups Type; };
inline void assign(AAGroups & target, AminoAcid const c_source)
{
SEQAN_CHECKPOINT
	target.value = _Translate_Table_AA_2_AAGroups<>::VALUE[(unsigned char) c_source];
}

//////////////////////////////////////////////////////////////////////////////

template <>
struct CompareType<AAGroups, Byte> { typedef AAGroups Type; };
inline void assign(AAGroups & target, Byte const c_source)
{
SEQAN_CHECKPOINT
	target.value = _Translate_Table_AA_2_AAGroups<>::VALUE[(unsigned char) _Translate_Table_Byte_2_AA<>::VALUE[c_source] ];
}

//////////////////////////////////////////////////////////////////////////////

template <>
struct CompareType<AAGroups, Ascii> { typedef AAGroups Type; };
inline void assign(AAGroups & target, Ascii const c_source)
{
SEQAN_CHECKPOINT
	target.value = _Translate_Table_AA_2_AAGroups<>::VALUE[(unsigned char) _Translate_Table_Ascii_2_AA<>::VALUE[(unsigned char) c_source] ];
}

//////////////////////////////////////////////////////////////////////////////

template <>
struct CompareType<AAGroups, Unicode> { typedef AAGroups Type; };
inline void assign(AAGroups & target, Unicode const c_source)
{
SEQAN_CHECKPOINT
	target.value = _Translate_Table_AA_2_AAGroups<>::VALUE[(unsigned char) _Translate_Table_Ascii_2_AA<>::VALUE[(unsigned char) c_source] ];
}

//////////////////////////////////////////////////////////////////////////////



}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
