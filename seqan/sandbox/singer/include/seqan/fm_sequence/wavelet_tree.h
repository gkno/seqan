// ==========================================================================
//                                  wavelet tree
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

#ifndef SANDBOX_MY_SANDBOX_APPS_FMINDEX_WAVELETTREE_H_
#define SANDBOX_MY_SANDBOX_APPS_FMINDEX_WAVELETTREE_H_

namespace seqan{
	
	template < typename TText, typename TFreq >
	void getNumCharsImpl(TText &text,
		TFreq &freq)
	{
		typedef typename Size< TText >::Type TSize;
		for (TSize i = 0; i < length(text); ++i)
		{
			++freq[ordValue(text[i])];
		}
	}
	template < typename TText, typename TFreq >
	void getNumChars(TText &text,
		TFreq &freq)
	{
		typedef typename Value< TText >::Type TValue;
		unsigned int sigmaSize = ValueSize< TValue >::VALUE;

		typedef typename Value< TFreq >::Type TFreqValue;
		resize(freq, sigmaSize, 0);
		getNumCharsImpl(text, freq);
	}
	
	template < typename TSequence, typename TFreq >
	void getNumChars(StringSet< TSequence > &text,
		TFreq &freq)
	{
		typedef typename Value< TSequence >::Type TValue;
		unsigned int sigmaSize = ValueSize< TValue >::VALUE;

		typedef typename Value< TFreq >::Type TFreqValue;
		resize(freq, sigmaSize, 0);
		for (unsigned i = 0; i < length(text); ++i)
		{
			getNumCharsImpl(getValue(text,i), freq);
		}
	}

	//function to create the prefixSum array
	template< typename TPrefixSumTable, typename TFreqTable, typename TAlphabetSize >
	void createPrefixTable(TPrefixSumTable &prefixSumTable, TFreqTable &freqTable, TAlphabetSize sigmaSize)
	{
		typedef typename Value< TPrefixSumTable >::Type TCounterValue;	
		clear(prefixSumTable);
		resize(prefixSumTable, sigmaSize + 1, 0);
						
		TCounterValue temp = 0;
		TCounterValue sum = 1;
		for(TAlphabetSize i = 0; i < sigmaSize; ++i)
		{
			temp = freqTable[i];
			prefixSumTable[i] = sum;
			sum += temp;
		}
		prefixSumTable[sigmaSize] = sum;
	}

	template< typename TText, typename TSpec >
	struct WaveletTree;

	struct SingleString_;
	struct MultiString_;
	struct FibreBitStrings_;
	struct FibreSplitValues_;

	typedef Tag< SingleString_ > SingleString;
	typedef Tag< MultiString_ > MultiString;
	typedef Tag< FibreBitStrings_ > const FibreBitStrings;
	typedef Tag< FibreSplitValues_ > const FibreSplitValues;

	template< typename TText, typename TSpec >
	struct Fibre< WaveletTree< TText, TSpec >, FibreBitStrings >
	{
		typedef StringSet< RankSupportBitString< void > > Type;
	};
	
	template< typename TText, typename TSpec >
	struct Fibre< WaveletTree< TText, TSpec >, FibreSplitValues >
	{
		typedef WaveletTreeStructure< TText > Type;
	};

	template< typename TText, typename TSpec >
	struct WaveletTree
	{
		typedef typename Fibre< WaveletTree< TText, TSpec >, FibreBitStrings >::Type TBitStrings;
		TBitStrings bitStrings;
		WaveletTreeStructure< TText > splitValues;

		WaveletTree(){}

		WaveletTree(TText const &text)
		{
			createWaveletTree(*this, text);	
		}

		template< typename TFreqTable >
		WaveletTree(TText const &text, TFreqTable const &freqTable)
		{
			createWaveletTree(*this, text, freqTable);	
		}
		
		template< typename TFreqTable, typename TPrefixSumTable >
		WaveletTree(TText const &text, TFreqTable const &freqTable, TPrefixSumTable const &prefixSumTable)
		{
			createWaveletTree(*this, text, freqTable, prefixSumTable);	
		}
	};

	template< typename TText >
	struct WaveletTree< TText, SingleString >
	{
		typedef typename BitVector_< BitsPerValue< typename Value< TText >::Type >::VALUE >::Type TValue;
		typedef typename Fibre< WaveletTree, FibreBitStrings >::Type 	TBitStrings;
		typedef typename Size< TText >::Type							TSize;
		
		TBitStrings 						bitStrings;
		WaveletTreeStructure< TText >		splitValues;
		TSize 								dollarPosition;
		TValue		 						dollarSub;

		WaveletTree(){}

		WaveletTree(TText const &text)
		{
			dollarPosition = (TSize)-1;
			dollarSub = 0;
			createWaveletTree(*this, text);	
		}

		template< typename TFreqTable >
		WaveletTree(TText const &text, TFreqTable const &freqTable)
		{
			dollarPosition = (TSize)-1;
			dollarSub = 0;
			createWaveletTree(*this, text, freqTable);	
		}
		
		template< typename TFreqTable, typename TPrefixSumTable >
		WaveletTree(TText const &text, TFreqTable const &freqTable, TPrefixSumTable const &prefixSumTable)
		{
			dollarPosition = (TSize)-1;
			dollarSub = 0;
			createWaveletTree(*this, text, freqTable, prefixSumTable);	
		}
		
		bool operator==(const WaveletTree &b) const
		{

			typedef typename Size< TText >::Type							TSize;
			bool test = true;
			if(length(bitStrings) == length(b.bitStrings))
			{				
				for(TSize i = 0; i < length(bitStrings); ++i)
				{
					if(!(bitStrings[i] == b.bitStrings[i]))
					{
						test = false;
					}
				}
			}
			else
			{
				test = false;
			}

			return (test && 
					splitValues == b.splitValues && 
					dollarPosition == b.dollarPosition && 
					dollarSub == b.dollarSub);
		}

	};

	template< typename TText >
	struct WaveletTree< TText, MultiString >
	{
		typedef typename BitVector_< BitsPerValue< typename Value< TText >::Type >::VALUE >::Type 	TValue;
		typedef typename Fibre< WaveletTree, FibreBitStrings >::Type 								TBitStrings;
		typedef typename Size< TText >::Type														TSize;
		
		TBitStrings 						bitStrings;
		WaveletTreeStructure< TText >		splitValues;
		TBitStrings							dollarPosition;
		TValue		 						dollarSub;

		WaveletTree(){}

		WaveletTree(TText const &text)
		{
			dollarSub = 0;
			createWaveletTree(*this, text);	
		}

		template< typename TFreqTable >
		WaveletTree(TText const &text, TFreqTable const &freqTable)
		{
			dollarSub = 0;
			createWaveletTree(*this, text, freqTable);	
		}
		
		template< typename TFreqTable, typename TPrefixSumTable >
		WaveletTree(TText const &text, TFreqTable const &freqTable, TPrefixSumTable const &prefixSumTable)
		{
			dollarSub = 0;
			createWaveletTree(*this, text, freqTable, prefixSumTable);	
		}
		
		bool operator==(const WaveletTree &b) const
		{
			typedef typename Size< TText >::Type TSize;
			bool test = true;
			if(length(bitStrings) == length(b.bitStrings))
			{				
				for(TSize i = 0; i < length(bitStrings); ++i)
				{
					if(!(bitStrings[i] == b.bitStrings[i]))
					{
						test = false;
					}
				}
			}
			else
			{
				test = false;
			}

			return (test && 
					splitValues == b.splitValues && 
					dollarPosition == b.dollarPosition && 
					dollarSub == b.dollarSub);
		}

	};
	
	template< typename TText, typename TSpec >
	typename Fibre< WaveletTree< TText, TSpec >, FibreBitStrings >::Type &
	getFibre(WaveletTree< TText, TSpec > &tree, const FibreBitStrings)
	{
		return tree.bitStrings;
	}

	template< typename TText, typename TSpec >
	typename Fibre< WaveletTree< TText,TSpec >, FibreBitStrings >::Type const &
	getFibre(const WaveletTree< TText, TSpec > &tree, const FibreBitStrings)
	{
		return tree.bitStrings;
	}

	template< typename TText, typename TSpec >
	typename Fibre< WaveletTree< TText, TSpec >, FibreSplitValues >::Type &
	getFibre(WaveletTree< TText, TSpec > &tree, FibreSplitValues)
	{
		return tree.splitValues;
	}

	template< typename TText, typename TSpec >
	typename Fibre< WaveletTree< TText, TSpec >, FibreSplitValues >::Type const &
	getFibre(WaveletTree< TText, TSpec > const &tree, const FibreSplitValues)
	{
		return tree.splitValues;
	}
	
	template<  
		typename TText, 
		typename TWaveletTreeSpec,
		typename TTreeSplitValue, 
		typename TPos >
	unsigned getOccImpl(const WaveletTree< TText, TWaveletTreeSpec > &tree,
		const TTreeSplitValue ordChar,
		const TPos pos)
	{
		typedef typename BitVector_< BitsPerValue< typename Value< TText >::Type >::VALUE >::Type TValue;
		TPos sum = pos + 1;
		TValue treePos = 0;
		
		typename Iterator< const WaveletTreeStructure< TText > >::Type iter(tree.splitValues, treePos);
		do{
			TPos addValue = getRank(tree.bitStrings[treePos], sum - 1);
			if(ordChar <= getFibre(tree, FibreSplitValues()).treeNodes[treePos].i1)
			{
				sum -= addValue;
				goLeft(iter);
			}
			else
			{	
				sum = addValue;
				goRight(iter);
			}
			treePos = getPosition(iter);
			if(sum == 0)
			{
				return 0;
			}
		}while(treePos);
				
		return sum;
	}

	template< 
		typename TText, 
		typename TChar, 
		typename TPos >
	unsigned getOcc(const WaveletTree< TText, SingleString > &tree,
		const TChar character,
		const TPos pos)
	{
		typedef typename BitVector_< BitsPerValue< TChar >::VALUE >::Type TTreeSplitValue;
		TTreeSplitValue ordChar = ordValue(character);
		unsigned occ = getOccImpl(tree, ordChar, pos);
		if(ordChar == tree.dollarSub && pos >= tree.dollarPosition)
		{
			return occ - 1;
		}
		return occ;
	}
	
	template< 
		typename TText, 
		typename TChar, 
		typename TPos >
	unsigned getOcc(const WaveletTree< TText, MultiString > &tree,
		const TChar character,
		const TPos pos)
	{
		typedef typename BitVector_< BitsPerValue< TChar >::VALUE >::Type TTreeSplitValue;
		TTreeSplitValue ordChar = ordValue(character);
		unsigned occ = getOccImpl(tree, ordChar, pos);
		unsigned numDollar = getRank(tree.dollarPosition, pos);
		if(ordChar == tree.dollarSub)
		{
			return occ - numDollar;
		}
		return occ;
	}
	
	template<  
		typename TText, 
		typename TWaveletTreeSpec,
		typename TPos >
	typename BitVector_< BitsPerValue< typename Value< TText >::Type >::VALUE >::Type
	getCharacter(const WaveletTree< TText, TWaveletTreeSpec > &tree,
		const TPos pos)
	{
		typedef typename BitVector_< BitsPerValue< typename Value< TText >::Type >::VALUE >::Type TValue;
		TPos sum = pos + 1;
		TValue treePos = 0;
		typename Iterator< const WaveletTreeStructure< TText > >::Type iter(tree.splitValues, treePos);
		bool direction;
		TValue character;
		do
		{
			direction = getBit(tree.bitStrings[treePos], sum - 1);
			TPos addValue = getRank(tree.bitStrings[treePos], sum - 1);
			if(direction)
			{
				character = getCharacter(iter) + 1;
				sum = addValue;
				goRight(iter);
			}
			else
			{	
				character = getCharacter(iter);
				sum -= addValue;
				goLeft(iter);
			}
			treePos = getPosition(iter);
		}while(treePos);
		
		return character;
	}
	
	template < typename TContainer, typename TPos, typename TBitsPerNumber, typename TValue >
	void numberToBits(TContainer &container, TPos pos, TBitsPerNumber bpn, TValue value)
	{
		//determine the underlying type of the container (e.g.: char, int, unsigned, long, ...)
		typedef typename Value< TContainer >::Type TContainerValue;
		unsigned bitsPerContainerValue = BitsPerValue< TContainerValue >::VALUE;
		TContainerValue insertValue = value;

		//determine the right position in the container
		pos *= bpn;
		Pair< unsigned, unsigned > containerPos = Pair<unsigned, unsigned>(pos/bitsPerContainerValue, pos%bitsPerContainerValue); 
	
		int numShifts = bitsPerContainerValue - containerPos.i2 - bpn;

		//the block is fully contained within one container
		if(numShifts >= 0)
		{
			if(length(container) == containerPos.i1)
			{
				container += 0;
			}
			container[containerPos.i1] &= ~(((1 << bpn) -1) << numShifts);
			container[containerPos.i1] |=  (insertValue << numShifts);

		}
		else
		{
			container[containerPos.i1] &= (-1u << (bpn + numShifts));
			container[containerPos.i1] |=  (insertValue >> (-1*numShifts));

			if(length(container) == containerPos.i1+1)
			{
				container += 0;
			}
			container[containerPos.i1+1] &= (-1u >> (bpn + numShifts));
			container[containerPos.i1+1] |=  (insertValue << (bitsPerContainerValue + numShifts));
		}
	}


	template < typename TReturn, typename TContainer, typename TPos, typename TBitsPerNumber >
	TReturn bitsToNumber(TContainer container, TPos pos, TBitsPerNumber bpn)
	{
		//determine the underlying type of the container (e.g.: char, int, unsigned, long, ...)
		typedef typename Value< TContainer >::Type TContainerValue;
		unsigned bitsPerContainerValue = BitsPerValue< TContainerValue >::VALUE;
		TContainerValue bitMask;

		//determine the right position in the container
		pos *= bpn;
		Pair< unsigned, unsigned > containerPos = Pair<unsigned, unsigned>(pos/bitsPerContainerValue, pos%bitsPerContainerValue); 
	
		//determine wheter a border between two block is crossed
		int numShifts = bitsPerContainerValue - containerPos.i2 - bpn;

		//the block is fully contained within one container
		TReturn value;// = 0;
		if(numShifts >= 0)
		{
			bitMask =  (1 << bpn) -1;
			value = ((container[containerPos.i1] >> numShifts) & bitMask);
		}
		else
		{
			bitMask =  (1 << bpn) -1;
			value = ((container[containerPos.i1] << (-1 * numShifts)) & bitMask);
			//did it so complicated because sometimes 1 and sometimes 0 were shifted from the left (depends ont the return type signed/unsigned)
			value |= (container[containerPos.i1 + 1] >> (bitsPerContainerValue + numShifts) & ((1 << (-1*numShifts)) -1));
		}

		return value;
	}

	//this function determines 2^value	
	unsigned pow(unsigned value)
	{
		return (1 << value);
	}

	unsigned getTreeLevel(unsigned treePosition)
	{
		return floor(log(treePosition + 1)/log(2));
	}

	//this function is used to fill all bit strings in a wavelet tree AND counting
	template < typename TText, typename TWaveletTreeSpec >// typename TTreeNodeWidth >
	void fillWaveletTree(
			WaveletTree< TText, TWaveletTreeSpec > &tree, 		//the set of bit strings to be filled
			typename Iterator< WaveletTreeStructure< TText > >::Type iter,
			const TText &text,		 					//the original string
			const typename BitVector_< BitsPerValue< typename Value< TText >::Type >::VALUE >::Type lowerBound, 			//the lowest value to be considered
			const typename BitVector_< BitsPerValue< typename Value< TText >::Type >::VALUE >::Type upperBound)			//the highest value to be considered
	{

		typedef typename BitVector_< BitsPerValue< typename Value< TText >::Type >::VALUE >::Type TValue;
		typename Iterator< WaveletTreeStructure< TText > >::Type iterRight = iter;

		fillBitString( 
			tree,
			iter,
			text,
			lowerBound,
			upperBound);

		TValue leftUpperBound = getCharacter(iter);
		goLeft(iter);
		if(getPosition(iter))
		{
			fillWaveletTree(
				   	tree,
					iter,	
					text,
					lowerBound,
					leftUpperBound);
		}

		TValue rightLowerBound = ordValue(getCharacter(iterRight)) + 1;
		goRight(iterRight);
		if(getPosition(iterRight))
		{
			fillWaveletTree(
				   	tree,
					iterRight,	
					text,
					rightLowerBound,
					upperBound);
		}
	}

/*	//determine the start bucket
	template< typename TSpec >
	TSpec getFinalWidth(TSpec &width)
	{
		TSpec bitsPerValue = BitsPerValue< TSpec >::VALUE;
		if(width % bitsPerValue)
		{
			return width / bitsPerValue + 1;
		}
		return width / bitsPerValue;
	}

	template< typename TValue, typename TBPV >
	TValue getFinalWidth(TValue width, TBPV bitsPerValue)
	{
		if(width % bitsPerValue)
		{
			return width / bitsPerValue + 1;
		}
		return width / bitsPerValue;
	}*/

	
	template < 
			 typename TBWT, 
			 typename TWaveletTreeSpec,
			 typename TFreqTable,
			 typename TPrefixSumTable >// typename TTreeNodeWidth >
	void createWaveletTree(WaveletTree< 
							TBWT, 
							TWaveletTreeSpec
						> &tree, 
						TBWT const &bwt,
						TFreqTable const &freq,
						TPrefixSumTable const &prefixSumTable)
	{
		typedef typename BitVector_< BitsPerValue< typename Value< TBWT >::Type >::VALUE >::Type TValue;
		
		//generate the tree structure
		typedef typename Size< typename Value< TBWT >::Type >::Type TSize;
		TSize sigmaSize = ValueSize< typename Value< TBWT >::Type >::VALUE;
		TValue numberOfTreeNodes = sigmaSize - 1; 
		resize(tree.splitValues, numberOfTreeNodes, Pair<TValue, TValue >(0, 0));
		
		TValue smallestValue = 0;
		TValue highestValue = sigmaSize - 1;
		TValue numChildNodes = 0;

		typename Iterator< WaveletTreeStructure< TBWT > >::Type iter(tree.splitValues, 0);
		computeTreeEntries(freq, 
			iter, 
			smallestValue, 
			highestValue, 
			numChildNodes);

		numberOfTreeNodes = numChildNodes;
		resize(tree.splitValues, numberOfTreeNodes);
		//writeGraph(tree.splitValues);

		String< unsigned long long > lengthString;
		computeStringLengthFromTree(lengthString,
			prefixSumTable,
	   		tree.splitValues);

		
		resize(tree.bitStrings, numberOfTreeNodes);
		for(unsigned i = 0; i < numberOfTreeNodes; ++i)
		{
			resize(tree.bitStrings[i], lengthString[i]);
		}

		setPosition(iter, (TValue)0);

		fillWaveletTree(	
				tree,
				iter,
				bwt,
				(TValue)0,
				(TValue)-1);
	}

	template < 
			 typename TBWT, 
			 typename TWaveletTreeSpec,
			 typename TFreqTable >
	void createWaveletTree(WaveletTree< 
							TBWT, 
							TWaveletTreeSpec
						> &tree, 
						TBWT const &bwt,
						TFreqTable const &freqTable)
	{
		
		typedef typename Size< typename Value< TBWT >::Type >::Type TValue;
		TValue sigmaSize = ValueSize< typename Value< TBWT >::Type >::VALUE;
		typedef TFreqTable TPrefixSumTable;
		TPrefixSumTable prefixSumTable;
		createPrefixTable(prefixSumTable, freqTable, sigmaSize);
	
		createWaveletTree(tree, bwt, freqTable, prefixSumTable);
	}
	//init the CounterWavelet-tree 	
	template < 
			 typename TBWT, 
			 typename TWaveletTreeSpec >
	void createWaveletTree(WaveletTree< 
							TBWT, 
							TWaveletTreeSpec
						> &tree, 
						TBWT const &bwt)
	{

		typedef String< typename Size< TBWT >::Type > TFreqTable;
		TFreqTable freqTable;
		getNumChars(bwt, freqTable);

		createWaveletTree(tree, bwt, freqTable);
	}

	template< typename TText >
	inline bool open(
		WaveletTree< TText, SingleString > &tree, 
		const char *fileName,
		int openMode)
	{
		String<char> name;
		
		typedef typename BitVector_< BitsPerValue< typename Value< TText >::Type >::VALUE >::Type TValue;
		typedef typename Size< TText >::Type							TSize;

		String< Pair< TValue, TSize > > dollarValues;

		name = fileName;	append(name, ".dollar");
		if(!open(dollarValues, toCString(name), openMode))
		{
			return false;
		}
		tree.dollarSub = dollarValues[0].i1;
		tree.dollarPosition = dollarValues[0].i2;
		name = fileName;	append(name, ".tree");	open(getFibre(tree, FibreBitStrings()), toCString(name), openMode);
		name = fileName;	append(name, ".split");	open(getFibre(tree, FibreSplitValues()), toCString(name), openMode);
		return true;
	}

	template < typename TText >
	inline bool open(
		WaveletTree< TText, SingleString > &tree, 
		const char *fileName)
	{
		return open(tree, fileName, DefaultOpenMode< WaveletTree< TText, SingleString > >::VALUE);
	}

	template< typename TText >
	inline bool save(
		WaveletTree< TText, SingleString > const &tree, 
		const char *fileName,
		int openMode)
	{
		String<char> name;

		typedef typename BitVector_< BitsPerValue< typename Value< TText >::Type >::VALUE >::Type TValue;
		typedef typename Size< TText >::Type							TSize;

		String< Pair< TValue, TSize > > dollarValues;
		append(dollarValues, Pair<TValue, TSize >(tree.dollarSub, tree.dollarPosition));

		name = fileName;	append(name, ".dollar");save(dollarValues, toCString(name), openMode);
		name = fileName;	append(name, ".tree");	save(getFibre(tree, FibreBitStrings()), toCString(name), openMode);
		name = fileName;	append(name, ".split");	save(getFibre(tree, FibreSplitValues()), toCString(name), openMode);
		return true;
	}

	template < typename TText >
	inline bool save(
		WaveletTree< TText, SingleString > const &tree, 
		const char *fileName)
	{
		return save(tree, fileName, DefaultOpenMode< WaveletTree< TText, SingleString > >::VALUE);
	}
}
#endif  // #ifndef SANDBOX_MY_SANDBOX_APPS_FMINDEX_WAVELETTREE_H_

