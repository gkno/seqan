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

#ifndef SEQAN_HEADER_GRAPH_ALIGN_TCOFFEE_KMER_H
#define SEQAN_HEADER_GRAPH_ALIGN_TCOFFEE_KMER_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// T-Coffee - Simple k-mer counter
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template<typename TString, typename TTupelString, typename TKTup, typename TAlphabet>
inline void
_getTupelString(TString const& str, 
				TTupelString& tupelString,
				TKTup const ktup, 
				TAlphabet) 
{
	SEQAN_CHECKPOINT
	typedef typename Value<typename Value<TTupelString>::Type>::Type TWord;
	
	// Alphabet size
	TWord alphabet_size = ValueSize<TAlphabet>::VALUE;

	// Assign a unique number to each k-tupel
	String<TWord> prod;  // Scaling according to position in k-tupel
	resize(prod,ktup);
	for (TWord i=0; i< (TWord) ktup;++i) {
		prod[ktup-i-1] = 1;
		for(TWord j=0;j<i;++j) prod[ktup-i-1] *= alphabet_size;
	}

	TWord len = length(str);
	clear(tupelString);
	resize(tupelString, len-(ktup - 1)); 
	TWord tupelIndex = 0;
	TWord endTupel = 0;
	tupelString[tupelIndex] = 0;
	for(;endTupel< (TWord) ktup;++endTupel) {
		tupelString[tupelIndex] += (TWord) ((Byte) ((TAlphabet) str[endTupel])) * prod[endTupel];
	}
	++tupelIndex;
	for(;endTupel<len;++endTupel) {
		tupelString[tupelIndex] = tupelString[tupelIndex - 1];
		tupelString[tupelIndex] -= (TWord) ((Byte) ((TAlphabet) str[endTupel - ktup])) * prod[0];
		tupelString[tupelIndex] *= alphabet_size;
		tupelString[tupelIndex] += (TWord) ((Byte) ((TAlphabet) str[endTupel]));
		++tupelIndex;
	}
}

//////////////////////////////////////////////////////////////////////////////

template<typename TString, typename TTupelString, typename TKTup, typename TAlphabet>
inline void
_getNonOverlappingTupelString(TString const& str, 
							  TTupelString& tupelString,
							  TKTup const ktup, 
							  TAlphabet) 
{
	SEQAN_CHECKPOINT
	typedef typename Value<typename Value<TTupelString>::Type>::Type TWord;
	
	// Alphabet size
	TWord alphabet_size = ValueSize<TAlphabet>::VALUE;

	// Assign a unique number to each k-tupel
	String<TWord> prod;  // Scaling according to position in k-tupel
	resize(prod,ktup);
	for (TWord i=0; i< (TWord) ktup;++i) {
		prod[ktup-i-1] = 1;
		for(TWord j=0;j<i;++j) prod[ktup-i-1] *= alphabet_size;
	}
	TWord len = length(str);
	clear(tupelString);
	TWord lenOfTup = (TWord) std::floor((double) len / (double) ktup);
	fill(tupelString, lenOfTup, 0, Exact()); 
	TWord tupelIndex = 0;
	TWord endTupel = 0;
	tupelString[tupelIndex] = 0;
	while (tupelIndex < lenOfTup) {
		tupelString[tupelIndex] += (TWord) ( (Byte) ( (TAlphabet) str[endTupel] ) ) * prod[endTupel % 3];
		++endTupel;
		if (endTupel % 3 == 0) ++tupelIndex;
	}
}

//////////////////////////////////////////////////////////////////////////////

template<typename TString, typename TSpec, typename THitMatrix, typename TSize, typename TAlphabet>
inline void
getKmerSimilarityMatrix(StringSet<TString, TSpec> const& strSet, 
						THitMatrix& mat, 
						TSize ktup, 
						TAlphabet) 
{
	SEQAN_CHECKPOINT
	typedef __int64 TWord;
	typedef String<TWord> TTupelString;
	typedef String<TTupelString> TTupelStringSet;
	typedef typename Value<THitMatrix>::Type TValue;

	// Number of sequences
	TSize nseq = length(strSet);
	TSize alphabet_size = ValueSize<TAlphabet>::VALUE;

	// Initialization
	// Matrix for common k-tupels between sequence i and j
	clear(mat);
	fill(mat, nseq*nseq, 0);

	// Transform the set of strings into a set of strings of k-tupels
	TTupelStringSet tupSet;
	resize(tupSet, length(strSet));
	for(TSize k=0;k<(TSize) length(strSet);++k) _getTupelString(strSet[k], tupSet[k], ktup, TAlphabet());

	// Build for each sequence the q-gram Index and count common hits
	String<TWord> qIndex;
	String<TWord> compareIndex;
	for(TSize k=0;k<nseq;++k) {
		clear(qIndex);
		fill(qIndex, (TWord) (double) pow((double)alphabet_size, (double)ktup), (TWord) 0, Exact());
		for(TSize i = 0;i < (TSize) length(tupSet[k]);++i) ++qIndex[ tupSet[k][i] ];
		TWord value;
	    for (TSize k2=k; k2<nseq; ++k2) {
			clear(compareIndex);
			fill(compareIndex, (TWord) (double) pow((double)alphabet_size, (double)ktup), (TWord) 0, Exact());
			value = 0;
			for(TSize i = 0;i < (TSize) length(tupSet[k2]);++i) {
				//std::cout << tupSet[k2][i] << "," << compareIndex[ tupSet[k2][i] ] << "," << qIndex[ tupSet[k2][i] ]<< std::endl;
				if (compareIndex[ tupSet[k2][i] ] < qIndex[ tupSet[k2][i] ]) ++value;
				++compareIndex[ tupSet[k2][i] ];
			}
			assignValue(mat, k*nseq+k2, (TValue) value);
		}
	}

	// Copy upper triangle to lower triangle and scale
	for(TWord row = 0; row < (TWord) nseq; ++row) {
		for(TWord col = row+1; col < (TWord) nseq; ++col) {
			// Fractional common kmer count
			// = Number of common q-grams / Number of possible common q-grams
			TValue val = getValue(mat, row*nseq+col);
			TValue minVal = getValue(mat, row*nseq+row);
			if (getValue(mat, col*nseq+col) < minVal) minVal = getValue(mat, col*nseq+col);
			val /= minVal;

			// Assign the values
			assignValue(mat, row*nseq+col, val);
			assignValue(mat, col*nseq+row, val);
		}
		assignValue(mat, row*nseq+row, 1);
	}

	//// Debug code
	//for (TSize row=0;row<nseq;++row) {
	//	for(TSize col=0;col<nseq;++col) {
	//		std::cout << getValue(mat, row*nseq+col) << ",";
	//	}
	//	std::cout << std::endl;
	//}
	//std::cout << std::endl;
}




}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
