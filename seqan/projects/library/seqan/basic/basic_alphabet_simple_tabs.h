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

#ifndef SEQAN_HEADER_BASIC_ALPHABET_SIMPLE_TABS_H
#define SEQAN_HEADER_BASIC_ALPHABET_SIMPLE_TABS_H


namespace SEQAN_NAMESPACE_MAIN
{
//////////////////////////////////////////////////////////////////////////////

//use same tables for Dna

template <typename T = void>
struct _Translate_Table_Dna5_2_Ascii
{
	static char const VALUE[5];
};
template <typename T>
char const _Translate_Table_Dna5_2_Ascii<T>::VALUE[5] = {'A', 'C', 'G', 'T', 'N'};

//____________________________________________________________________________

template <typename T = void>
struct _Translate_Table_Dna5_2_Iupac
{
	static char const VALUE[5];
};
template <typename T>
char const _Translate_Table_Dna5_2_Iupac<T>::VALUE[5] = {0x02, 0x04, 0x08, 0x01, 0x0f};

//____________________________________________________________________________

// use same table for Rna

template <typename T = void>
struct _Translate_Table_Rna5_2_Ascii
{
	static char const VALUE[5];
};
template <typename T>
char const _Translate_Table_Rna5_2_Ascii<T>::VALUE[5] = {'A', 'C', 'G', 'U', 'N'}; 

//____________________________________________________________________________

template <typename T = void>
struct _Translate_Table_Iupac_2_Ascii
{
	static char const VALUE[16];
};
template <typename T>
char const _Translate_Table_Iupac_2_Ascii<T>::VALUE[16] = 
{
	'U', //0000=0
	'T', //0001=1 //T=1: change between U and T is just inc/dec
	'A', //0010=2
	'W', //0011=3 TA
	'C', //0100=4 
	'Y', //0101=5 TC (pyrimidine)
	'M', //0110=6 AC
	'H', //0111=7 not-G
	'G', //1000=8
	'K', //1001=9 TG
	'R', //1010=A AG (purine)
	'D', //1011=B not-C
	'S', //1100=C CG
	'B', //1101=D non-A
	'V', //1110=E non-T
	'N'  //1111=F any
};

//____________________________________________________________________________

template <typename T = void>
struct _Translate_Table_Iupac_2_Dna
{
	static char const VALUE[16];
};
template <typename T>
char const _Translate_Table_Iupac_2_Dna<T>::VALUE[16] = 
{
	3, //'U'
	3, //'T'
	0, //'A'
	0, //'W' = TA
	1, //'C' 
	1, //'Y' = TC
	0, //'M' = AC
	0, //'H' = not-G
	2, //'G'
	2, //'K' = TG
	0, //'R' = AG
	0, //'D' = not-C
	1, //'S' = CG
	1, //'B' = non-A
	0, //'V' = non-T
	0  //'N' = any
};

//____________________________________________________________________________


template <typename T = void>
struct _Translate_Table_Iupac_2_Dna5
{
	static char const VALUE[16];
};
template <typename T>
char const _Translate_Table_Iupac_2_Dna5<T>::VALUE[16] = 
{
	3, //'U'
	3, //'T'
	0, //'A'
	4, //'W' = TA
	1, //'C' 
	4, //'Y' = TC
	4, //'M' = AC
	4, //'H' = not-G
	2, //'G'
	4, //'K' = TG
	4, //'R' = AG
	4, //'D' = not-C
	4, //'S' = CG
	4, //'B' = non-A
	4, //'V' = non-T
	4  //'N' = any
};

//____________________________________________________________________________

template <typename T = void>
struct _Translate_Table_Ascii_2_Dna
{
	static char const VALUE[256];
};
template <typename T>
char const _Translate_Table_Ascii_2_Dna<T>::VALUE[256] = 
{
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //0
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //1
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //2
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //3

	0,   0,   0,   1,   0,   0,   0,   2,   0,   0,   0,   0,   0,   0,   0,   0, //4
//	 ,   A,   B,   C,   D,   E,   D,   G,   H,   I,   J,   K,   L,   M,   N,   O,

	0,   0,   0,   0,   3,   3,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //5
//	P,   Q,   R,   S,   T,   U,   V,   W,   X,   Y,   Z,    ,    ,    ,    ,    

	0,   0,   0,   1,   0,   0,   0,   2,   0,   0,   0,   0,   0,   0,   0,   0, //6
//   ,   a,   b,   c,   d,   e,   f,   g,   h,   i,   j,   k,   l,   m,   n,   o,

	0,   0,   0,   0,   3,   3,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //7
//  p,   q,   r,   s,   t,   u,   v,   w,   x,   y,   z,    ,    ,    ,    ,   

	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //8
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //9
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //10
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //11
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //12
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //13
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //14
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0  //15
};

//____________________________________________________________________________


template <typename T = void>
struct _Translate_Table_Ascii_2_Dna5
{
	static char const VALUE[256];
};
template <typename T>
char const _Translate_Table_Ascii_2_Dna5<T>::VALUE[256] = 
{
	4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4, //0
	4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4, //1
	4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4, //2
	4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4, //3

	4,   0,   4,   1,   4,   4,   4,   2,   4,   4,   4,   4,   4,   4,   4,   4, //4
//	 ,   A,   B,   C,   D,   E,   D,   G,   H,   I,   J,   K,   L,   M,   N,   O,

	4,   4,   4,   4,   3,   3,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4, //5
//	P,   Q,   R,   S,   T,   U,   V,   W,   X,   Y,   Z,    ,    ,    ,    ,    

	4,   0,   4,   1,   4,   4,   4,   2,   4,   4,   4,   4,   4,   4,   4,   4, //6
//   ,   a,   b,   c,   d,   e,   f,   g,   h,   i,   j,   k,   l,   m,   n,   o,

	4,   4,   4,   4,   3,   3,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4, //7
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

//____________________________________________________________________________

template <typename T = void>
struct _Translate_Table_Ascii_2_Rna
{
	static char const VALUE[256];
};
template <typename T>
char const _Translate_Table_Ascii_2_Rna<T>::VALUE[256] = 
{
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //0
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //1
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //2
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //3

	0,   0,   0,   1,   0,   0,   0,   2,   0,   0,   0,   0,   0,   0,   0,   0, //4
//	 ,   A,   B,   C,   D,   E,   D,   G,   H,   I,   J,   K,   L,   M,   N,   O,

	0,   0,   0,   0,   0,   3,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //5
//	P,   Q,   R,   S,   T,   U,   V,   W,   X,   Y,   Z,    ,    ,    ,    ,    

	0,   0,   0,   1,   0,   0,   0,   2,   0,   0,   0,   0,   0,   0,   0,   0, //6
//   ,   a,   b,   c,   d,   e,   f,   g,   h,   i,   j,   k,   l,   m,   n,   o,

	0,   0,   0,   0,   0,   3,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //7
//  p,   q,   r,   s,   t,   u,   v,   w,   x,   y,   z,    ,    ,    ,    ,   

	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //8
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //9
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //10
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //11
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //12
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //13
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //14
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0  //15
};

//____________________________________________________________________________

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

//____________________________________________________________________________

template <typename T = void>
struct _Translate_Table_Ascii_2_Iupac
{
	static char const VALUE[256];
};
template <typename T>
char const _Translate_Table_Ascii_2_Iupac<T>::VALUE[256] = 
{
	15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15, //0
	15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15, //1
	15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15, //2
	15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15, //3

	15,   2,  13,   4,  11,  15,  15,   8,   7,  15,  15,   9,  15,   6,  15,  15, //4
	//   ,   A,   B,   C,   D,   E,   F,   G,   H,   I,   J,   K,   L,   M,   N,   O,

	15,  15,  10,  12,   1,   0,  14,   3,  15,   5,  15,  15,  15,  15,  15,  15, //5
	//  P,   Q,   R,   S,   T,   U,   V,   W,   X,   Y,   Z,    ,    ,    ,    ,   

	15,   2,  13,   4,  11,  15,  15,   8,   7,  15,  15,   9,  15,   6,  15,  15, //6
	//   ,   a,   b,   c,   d,   e,   f,   g,   h,   i,   j,   k,   l,   m,   n,   o,

	15,  15,  10,  12,   1,   0,  14,   3,  15,   5,  15,  15,  15,  15,  15,  15, //7
	//  p,   q,   r,   s,   t,   u,   v,   w,   x,   y,   z,    ,    ,    ,    ,   

	15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15, //8
	15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15, //9
	15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15, //10
	15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15, //11
	15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15, //12
	15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15, //13
	15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15, //14
	15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15  //15
};

//____________________________________________________________________________

template <typename T = void>
struct _Translate_Table_Byte_2_Dna
{
	static char const VALUE[256];
};
template <typename T>
char const _Translate_Table_Byte_2_Dna<T>::VALUE[256] = 
{
	0,   1,   2,   3,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //0
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

//____________________________________________________________________________


template <typename T = void>
struct _Translate_Table_Byte_2_Dna5
{
	static char const VALUE[256];
};
template <typename T>
char const _Translate_Table_Byte_2_Dna5<T>::VALUE[256] = {
	0,   1,   2,   3,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4, //0
	4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4, //1
	4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4, //2
	4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4, //3
	4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4, //4
	4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4, //5
	4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4, //6
	4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4, //7
	4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4, //8
	4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4, //9
	4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4, //10
	4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4, //11
	4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4, //12
	4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4, //13
	4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4, //14
	4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4  //15
};

//____________________________________________________________________________

template <typename T = void>
struct _Translate_Table_Byte_2_Rna
{
	static char const VALUE[256];
};
template <typename T>
char const _Translate_Table_Byte_2_Rna<T>::VALUE[256] = 
{
	0,   1,   2,   3,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //0
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

//____________________________________________________________________________

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

//____________________________________________________________________________

template <typename T = void>
struct _Translate_Table_Byte_2_Iupac
{
	static char const VALUE[256];
};
template <typename T>
char const _Translate_Table_Byte_2_Iupac<T>::VALUE[256] = 
{
	0,   1,   2,   3,   4,   5,   6,   7,   8,   9,  10,  11,  12,  13,  14,  15, //0
	15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15, //1
	15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15, //2
	15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15, //3
	15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15, //4
	15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15, //5
	15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15, //6
	15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15, //7
	15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15, //8
	15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15, //9
	15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15, //10
	15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15, //11
	15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15, //12
	15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15, //13
	15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15, //14
	15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,  15  //15
};

//____________________________________________________________________________


template <typename T = void>
struct _Translate_Table_AA_2_Ascii
{
	static char const VALUE[24];
};
template <typename T>
char const _Translate_Table_AA_2_Ascii<T>::VALUE[24] = 
{
	'A', // 0 Ala Alanine                 
	'R', // 1 Arg Arginine                
	'N', // 2 Asn Asparagine              
	'D', // 3 Asp Aspartic Acid           
	'C', // 4 Cys Cystine                 
	'Q', // 5 Gln Glutamine               
	'E', // 6 Glu Glutamic Acid           
	'G', // 7 Gly Glycine                 
	'H', // 8 His Histidine               
	'I', // 9 Ile Isoleucine              
	'L', //10 Leu Leucine                 
	'K', //11 Lys Lysine                  
	'M', //12 Met Methionine              
	'F', //13 Phe Phenylalanine           
	'P', //14 Pro Proline                 
	'S', //15 Ser Serine                  
	'T', //16 Thr Threonine               
	'W', //17 Trp Tryptophan              
	'Y', //18 Tyr Tyrosine                
	'V', //19 Val Valine                  
	'B', //20 Aspartic Acid, Asparagine   
	'Z', //21 Glutamic Acid, Glutamine    
	'X', //22 Unknown                     
	'*'  //23 Terminator                  
};

//____________________________________________________________________________


template <typename T = void>
struct _Translate_Table_Ascii_2_AA
{
	static char const VALUE[256];
};
template <typename T>
char const _Translate_Table_Ascii_2_AA<T>::VALUE[256] = 
{
	22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22, //0
	22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22, //1
	22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  23,  22,  22,  22,  22,  22, //2
//													   *	
	22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22, //3
	22,   0,  20,   4,   3,   6,  13,   7,   8,   9,  22,  11,  10,  12,   2,  22, //4
//    ,   A,   B,   C,   D,   E,   F,   G,   H,   I,   J,   K,   L,   M,   N,   O,

	14,   5,   1,  15,  16,  22,  19,  17,  22,  18,  21,  22,  22,  22,  22,  22, //5
//   P,   Q,   R,   S,   T,   U,   V,   W,   X,   Y,   Z,    ,    ,    ,    ,    ,

	22,   0,  20,   4,   3,   6,  13,   7,   8,   9,  22,  11,  10,  12,   2,  22, //6
//    ,   a,   b,   c,   d,   e,   f,   g,   h,   i,   j,   k,   l,   m,   n,   o,

	14,   5,   1,  15,  16,  22,  19,  17,  22,  18,  21,  22,  22,  22,  22,  22, //7
//   p,   q,   r,   s,   t,   u,   v,   w,   x,   y,   z,    ,    ,    ,    ,    ,

	22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22, //8
	22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22, //9
	22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22, //10
	22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22, //11
	22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22, //12
	22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22, //13
	22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22, //14
	22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22  //15
};

//____________________________________________________________________________


template <typename T = void>
struct _Translate_Table_Byte_2_AA
{
	static char const VALUE[256];
};
template <typename T>
char const _Translate_Table_Byte_2_AA<T>::VALUE[256] = 
{
	0,   1,   2,   3,   4,   5,   6,   7,   8,   9,  10,  11,  12,  13,  14,  15, //0
	16,  17,  18,  19,  20,  21,  22,  23,  22,  22,  22,  22,  22,  22,  22,  22, //1
	22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22, //2
	22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22, //3
	22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22, //4
	22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22, //5
	22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22, //6
	22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22, //7
	22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22, //8
	22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22, //9
	22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22, //10
	22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22, //11
	22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22, //12
	22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22, //13
	22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22, //14
	22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22,  22  //15
};


//____________________________________________________________________________




}//namespace SEQAN_NAMESPACE_MAIN
#endif //#ifndef SEQAN_HEADER_...
