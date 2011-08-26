
// ==========================================================================
//                                  FMIndex
// ==========================================================================
// Copyright (c) 2006-2011, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Jochen Singer <your.email@example.net>
// ==========================================================================

#ifndef SANDBOX_MY_SANDBOX_APPS_FMINDEX_COMPRESSEDSA_H_
#define SANDBOX_MY_SANDBOX_APPS_FMINDEX_COMPRESSEDSA_H_

namespace seqan
{

	template< typename TOccTable, typename TPrefixSumTable >
	struct LFTable
	{
		TOccTable occTable;
		TPrefixSumTable prefixSumTable;

		LFTable() {};
		LFTable(TOccTable occTable, TPrefixSumTable prefixSumTable) :
			occTable(occTable),
			prefixSumTable(prefixSumTable)
		{};
	};

	struct FibreOccTable_;
	struct FibrePrefixSumTable_;

	typedef Tag< FibreOccTable_> const FibreOccTable;
	typedef Tag< FibrePrefixSumTable_> const FibrePrefixSumTable;

	template< typename TOccTable, typename TPrefixSumTable >
	struct Fibre< LFTable< TOccTable, TPrefixSumTable >, FibreOccTable >
	{
		typedef TOccTable Type;
	};
	
	template< typename TOccTable, typename TPrefixSumTable >
	struct Fibre< LFTable< TOccTable, TPrefixSumTable >, FibrePrefixSumTable >
	{
		typedef TPrefixSumTable Type;
	};

	template< typename TOccTable, typename TPrefixSumTable >
	typename Fibre< LFTable< TOccTable, TPrefixSumTable >, FibrePrefixSumTable >::Type &
	getFibre(LFTable< TOccTable, TPrefixSumTable > &lfTable, FibrePrefixSumTable)
	{
		return lfTable.prefixSumTable;
	}

	template< typename TOccTable, typename TPrefixSumTable >
	typename Fibre< LFTable< TOccTable, TPrefixSumTable >, FibreOccTable >::Type &
	getFibre(LFTable< TOccTable, TPrefixSumTable > &lfTable, FibreOccTable)
	{
		return lfTable.occTable;
	}
	
	template< typename TSparseString, typename TLfTable, typename TSpec >
	struct CompressedSA
	{
		TSparseString							compressedSA;
		TLfTable * 								lfTable;

		CompressedSA(){};

		typedef typename Size< typename Fibre< TSparseString, FibreSparseString >::Type >::Type TCompressedSaValue;
		template < typename TPos >
		TCompressedSaValue const operator[](TPos pos)
		{ 
			typedef typename Fibre< TSparseString, FibreIndicatorString>:: Type TIndicatorString;
			TIndicatorString const &indicatorString = getFibre(compressedSA, FibreIndicatorString());
			TPos counter = 0;
			while(!getBit(indicatorString, pos))
			{
				pos = lfMapping(*lfTable, pos);
				++counter;
			}
			return compressedSA[getRank(indicatorString, pos) - 1] + counter;
		}

		template < typename TPos >
		TCompressedSaValue operator[](TPos pos) const
		{
			typedef typename Fibre< TSparseString, FibreIndicatorString>:: Type TIndicatorString;
			TIndicatorString const &indicatorString = getFibre(compressedSA, FibreIndicatorString());
			TPos counter = 0;
			while(!getBit(indicatorString, pos))
			{
				pos = lfMapping(*lfTable, pos);
				++counter;
			}
			return getValue(compressedSA, getRank(indicatorString, pos) - 1) + counter;
		}
	};


	template< typename TSparseString, typename TLfTable, typename TSpec >
	struct Fibre< CompressedSA < TSparseString, TLfTable, TSpec >, FibreSA >
	{
		typedef TSparseString Type;
	};

	template< typename TSparseString, typename TLfTable, typename TSpec >
	typename Fibre< CompressedSA < TSparseString, TLfTable, TSpec >, FibreSA >::Type &
	getFibre(CompressedSA < TSparseString, TLfTable, TSpec > &compressedSA, FibreSA)
	{
		return compressedSA.compressedSA;
	}

/*	template< typename TSparseString, typename TLfTable, typename TSpec, typename TSize >
	void initCompressedSA(CompressedSA< TSparseString, TLfTable, TSpec > &compressedSA, 
			TLfTable &lfTable, 
			TSize size);*/

	template< typename TSparseString, typename TLfTable, typename TSpec >
	struct Iterator< CompressedSA < TSparseString, TLfTable, TSpec > const, Standard>
	{        
        typedef Iter<CompressedSA < TSparseString, TLfTable, TSpec > const, PositionIterator > Type;
	};
	
	template< typename TSparseString, typename TLfTable, typename TSpec >
	struct Iterator< CompressedSA < TSparseString, TLfTable, TSpec >, Standard >
	{        
        typedef Iter<CompressedSA < TSparseString, TLfTable, TSpec > , PositionIterator > Type;
	};
	
	template< typename TSparseString, typename TLfTable, typename TSpec >
	struct Iterator< CompressedSA < TSparseString, TLfTable, TSpec >, Rooted> :
		 Iterator< CompressedSA < TSparseString, TLfTable, TSpec >, Standard> {};
	
	template< typename TSparseString, typename TLfTable, typename TSpec >
	struct Iterator< CompressedSA < TSparseString, TLfTable, TSpec > const, Rooted> :
		 Iterator< CompressedSA < TSparseString, TLfTable, TSpec > const, Standard> {};

	template< typename TSparseString, typename TLfTable, typename TSpec >
	struct Value< CompressedSA< TSparseString, TLfTable, TSpec > > 
	{
        typedef typename Value< TSparseString >::Type Type;
	};

	template< typename TSparseString, typename TLfTable, typename TSpec >
	struct Value< CompressedSA < TSparseString, TLfTable, TSpec > const>
	{
        typedef typename Value< TSparseString >::Type const Type;
	};

}


#endif // SANDBOX_MY_SANDBOX_APPS_FMINDEX_COMPRESSEDSA_H_
