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

#ifndef SANDBOX_MY_SANDBOX_APPS_FMINDEX_WAVELETTREESTRUCTURE_H_
#define SANDBOX_MY_SANDBOX_APPS_FMINDEX_WAVELETTREESTRUCTURE_H_

namespace seqan{

	template< typename TText >
	struct WaveletTreeStructure;
	
	struct FibreTreeNodes_;
	typedef Tag< FibreTreeNodes_ > const FibreTreeNodes;

	template< typename TText >
	struct Fibre< WaveletTreeStructure< TText >, FibreTreeNodes >
	{
		typedef typename BitVector_< BitsPerValue< typename Value< TText >::Type >::VALUE >::Type TValue;
		typedef String< Pair < TValue, TValue > > Type;
	};

	template< typename TText >
	struct WaveletTreeStructure
	{
		typename Fibre< WaveletTreeStructure, FibreTreeNodes >::Type treeNodes;
	
		WaveletTreeStructure(){};

		template< typename TFreqString >
		explicit WaveletTreeStructure(TFreqString &string);
		
		bool operator==(const WaveletTreeStructure &b) const
		{
			return (treeNodes == b.treeNodes);
		}
	};
	
	template< typename TText >
	template< typename TFreqString >	
	WaveletTreeStructure< TText >::WaveletTreeStructure(TFreqString &freqString)
	{
		typedef typename BitVector_< BitsPerValue< typename Value< TText >::Type >::VALUE >::Type TValue;
		TValue sigmaSize = length(freqString);

		resize(treeNodes, sigmaSize - 1);
		TValue numChildNodes = 0;
		TValue lowest = 0;
		TValue highest = sigmaSize - 1;
		TValue posInTree = 0;
		computeTreeEntries(freqString, 
				*this,
				lowest,
				highest,
				posInTree,
				numChildNodes);
	}

	template< typename TText >
	typename Fibre< WaveletTreeStructure< TText >, FibreTreeNodes >::Type &
	getFibre(WaveletTreeStructure< TText > &treeStructure, FibreTreeNodes)
	{
		return treeStructure.treeNodes;
	}

	template< typename TText >
	typename Fibre< WaveletTreeStructure< TText >, FibreTreeNodes >::Type const &
	getFibre(WaveletTreeStructure< TText > const &treeStructure, FibreTreeNodes)
	{
		return treeStructure.treeNodes;
	}

	/*template< typename TBitString, typename TText, typename TSpec >
	typename Fibre< WaveletTree< TBitString, TSplitValue, TPosInSubTree, TSpec >, FibreSplitValues >::Type const &
	getFibre(WaveletTree< TBitString, TSplitValue, TPosInSubTree, TSpec > const &tree, const FibreSplitValues)
	{
		return tree.splitValues;
	}*/

	template< typename TText >
	unsigned length(WaveletTreeStructure< TText > &tree)
	{
		return length(tree.treeNodes);
	}
	
	template<typename TText, typename TSize >
	void resize(WaveletTreeStructure< TText > &treeStructure, TSize size)
	{
		resize(treeStructure.treeNodes, size);
	}
	
	template<typename TText, typename TSize >
	void resize(WaveletTreeStructure< TText > &treeStructure, TSize size, 
			Pair< typename BitVector_< BitsPerValue< typename Value< TText >::Type >::VALUE >::Type,
				typename BitVector_< BitsPerValue< typename Value< TText >::Type >::VALUE >::Type > value)
	{
		resize(treeStructure.treeNodes, size, value);
	}
	
	template<typename TText >
	void clear(WaveletTreeStructure< TText > &treeStructure)
	{
		clear(treeStructure.treeNodes);
	}

	template<typename TText >
	struct Iter< WaveletTreeStructure < TText > const, Standard>
	{        
		typedef typename BitVector_< BitsPerValue< typename Value< TText >::Type >::VALUE >::Type TValue;
		TValue position;
	    const WaveletTreeStructure< TText > * waveletTreeStructure;
	
		Iter(const WaveletTreeStructure< TText > &treeStructure, const TValue pos):
			position(pos),
			waveletTreeStructure(&treeStructure)
		{}
	};

	template< typename TText >
	struct Iterator< WaveletTreeStructure < TText > const, Standard>
	{        
		typedef Iter< WaveletTreeStructure < TText > const, Standard> Type;
	};
	
	template<typename TText >
	struct Iter< WaveletTreeStructure < TText >, Standard>
	{        
		typedef typename BitVector_< BitsPerValue< typename Value< TText >::Type >::VALUE >::Type TValue;
		TValue position;
		WaveletTreeStructure< TText > * waveletTreeStructure;

		Iter(WaveletTreeStructure< TText > &treeStructure, TValue pos):
			position(pos),
			waveletTreeStructure(&treeStructure)
		{}
	};

	template<typename TText >
	struct Iterator< WaveletTreeStructure < TText >, Standard>
	{        
		typedef Iter< WaveletTreeStructure < TText >, Standard> Type;
	};
	
	template<typename TText >
	struct Iterator< WaveletTreeStructure < TText >, Rooted> :
		 Iterator< WaveletTreeStructure < TText >, Standard > {};
	
	template<typename TText >
	struct Iterator< WaveletTreeStructure < TText > const, Rooted> :
		 Iterator< WaveletTreeStructure < TText > const, Standard> {};


	template<typename TText >
	struct Value< WaveletTreeStructure< TText > >
	{
		typedef typename BitVector_< BitsPerValue< typename Value< TText >::Type >::VALUE >::Type TValue;
		typedef Pair< TValue, TValue > Type;
	};

	template< typename TText >
	struct Value< WaveletTreeStructure< TText > const>
	{
		typedef typename BitVector_< BitsPerValue< typename Value< TText >::Type >::VALUE >::Type TValue;
		typedef Pair< TValue, TValue > const Type;
	};

	template< typename TText >
	struct Reference< WaveletTreeStructure< TText > >
	{
		typedef typename Value< WaveletTreeStructure< TText > >::Type Type;
	};

	template<typename TText >
	struct Reference< const WaveletTreeStructure< TText > >
	{
		typedef typename Value< WaveletTreeStructure< TText > >::Type const Type;
	};


	template<typename TText >
	inline typename Iterator< WaveletTreeStructure < TText > const >::Type 
	begin(WaveletTreeStructure< TText > const &waveletTreeStructure)
	{
		return typename Iterator< WaveletTreeStructure < TText > >::Type(waveletTreeStructure, 0);
	}
	
	template<typename TText >
	inline typename Iterator< WaveletTreeStructure < TText > const >::Type 
	end(WaveletTreeStructure< TText > const &waveletTreeStructure)
	{
		return typename Iterator< WaveletTreeStructure < TText > >::Type(waveletTreeStructure, length(waveletTreeStructure.treeNodes));
	}

	template<typename TText >
	inline typename Iterator< WaveletTreeStructure < TText > >::Type 
	begin(WaveletTreeStructure< TText > &waveletTreeStructure)
	{
		return typename Iterator< WaveletTreeStructure < TText > >::Type(waveletTreeStructure.treeNodes, 0);
	}
	
	template<typename TText >
	inline typename Iterator< WaveletTreeStructure < TText > >::Type 
	end(WaveletTreeStructure< TText > &waveletTreeStructure)
	{
		return typename Iterator< WaveletTreeStructure < TText > >::Type(waveletTreeStructure.treeNodes, length(waveletTreeStructure.treeNodes));
	}

	template<typename TText, typename TPos >
	void setPosition(Iter< WaveletTreeStructure< TText >, Standard > &iter, TPos pos)
	{
		iter.position = pos;
	}

	template<typename TText >
	typename BitVector_< BitsPerValue< typename Value< TText >::Type >::VALUE >::Type
	getPosition(Iter< WaveletTreeStructure< TText >, Standard > &iter)
	{
		return iter.position;
	}

	template<typename TText >
	typename BitVector_< BitsPerValue< typename Value< TText >::Type >::VALUE >::Type
	getPosition(Iter< const WaveletTreeStructure< TText >, Standard > &iter)
	{
		return iter.position;
	}

	template<typename TText >
	void setCharacter(Iter< WaveletTreeStructure< TText >, Standard > &iter, 
		typename BitVector_< BitsPerValue< typename Value< TText >::Type >::VALUE >::Type character)
	{
		iter.waveletTreeStructure->treeNodes[iter.position].i1 = character;
	}

	template<typename TText >
	typename BitVector_< BitsPerValue< typename Value< TText >::Type >::VALUE >::Type
	getCharacter(Iter< WaveletTreeStructure< TText >, Standard > &iter)
	{
		return iter.waveletTreeStructure->treeNodes[iter.position].i1;
	}

	template<typename TText >
	typename BitVector_< BitsPerValue< typename Value< TText >::Type >::VALUE >::Type
	getCharacter(Iter< const WaveletTreeStructure< TText >, Standard > &iter)
	{
		return iter.waveletTreeStructure->treeNodes[iter.position].i1;
	}

	template<typename TText >
	void setLeftChildPos(Iter< WaveletTreeStructure< TText >, Standard > &iter)
	{
		if(iter.waveletTreeStructure->treeNodes[iter.position].i2 == 0)
		{
			iter.waveletTreeStructure->treeNodes[iter.position].i2 = 2;
			return;
		}
		if(iter.waveletTreeStructure->treeNodes[iter.position].i2 == 2)
		{
			return;
		}
		std::cerr << "ERROR: The right child has just been deleted!" << std::endl;
	}


	template<typename TText >
	typename BitVector_< BitsPerValue< typename Value< TText >::Type >::VALUE >::Type
	getLeftChildPos(Iter< WaveletTreeStructure< TText >, Standard > &iter)
	{
		if(iter.waveletTreeStructure->treeNodes[iter.position].i2 > 1)
		{
			return iter.position + 1;
		}
		return 0;
	}

	template<typename TText >
	typename BitVector_< BitsPerValue< typename Value< TText >::Type >::VALUE >::Type
	getLeftChildPos(Iter< const WaveletTreeStructure< TText >, Standard > &iter)
	{
		if(iter.waveletTreeStructure->treeNodes[iter.position].i2 > 1)
		{
			return iter.position + 1;
		}
		return 0;
	}

	template<typename TText >
	void setRightChildPosOnly(Iter< WaveletTreeStructure< TText >, Standard > &iter)
	{
		iter.waveletTreeStructure->treeNodes[iter.position].i2 = 1;
	}

	template<typename TText, typename TPos >
	void setRightChildPos(Iter< WaveletTreeStructure< TText >, Standard > &iter, TPos rightChildPosition)
	{
		if(iter.waveletTreeStructure->treeNodes[iter.position].i2 == 0)
		{
			iter.waveletTreeStructure->treeNodes[iter.position].i2 = 1;
			return;
		}
		if(iter.waveletTreeStructure->treeNodes[iter.position].i2 == 2)
		{
			iter.waveletTreeStructure->treeNodes[iter.position].i2 = rightChildPosition + 2;
			return ;
		}
		if(iter.waveletTreeStructure->treeNodes[iter.position].i2 == 1)
		{
			return;
		}
		iter.waveletTreeStructure->treeNodes[iter.position].i2 = rightChildPosition + 2;
	}

	template<typename TText >
	typename BitVector_< BitsPerValue< typename Value< TText >::Type >::VALUE >::Type
	getRightChildPos(Iter< WaveletTreeStructure< TText >, Standard > &iter)
	{
		if(iter.waveletTreeStructure->treeNodes[iter.position].i2 > 2)
		{
			return iter.waveletTreeStructure->treeNodes[iter.position].i2 - 2;
		}
		if(iter.waveletTreeStructure->treeNodes[iter.position].i2 == 1)
		{
			return iter.position + 1;
		}
		return 0;
	}

	template<typename TText >
	typename BitVector_< BitsPerValue< typename Value< TText >::Type >::VALUE >::Type
	getRightChildPos(Iter< const WaveletTreeStructure< TText >, Standard > &iter)
	{
		if(iter.waveletTreeStructure->treeNodes[iter.position].i2 > 2)
		{
			return iter.waveletTreeStructure->treeNodes[iter.position].i2 - 2;
		}
		if(iter.waveletTreeStructure->treeNodes[iter.position].i2 == 1)
		{
			return iter.position + 1;
		}
		return 0;
	}

	//////////////////////////////////////////////////////////////////////////////////////////////
	template<typename TText >
	void setNodeToLeaf(Iter< WaveletTreeStructure< TText >, Standard > &iter)
	{
		if(iter.waveletTreeStructure->treeNodes[iter.position].i2 != 0)
		{
			std::cerr << "You just deleted ";
			if(iter.waveletTreeStructure->treeNodes[iter.position].i2 == 1)
			{
				std::cerr << "the right sub tree!" << std::endl;
			}
			if(iter.waveletTreeStructure->treeNodes[iter.position].i2 == 2)
			{
				std::cerr << "the left sub tree!" << std::endl;
			}
			else
			{
				std::cerr << "both sub trees!" << std::endl;
			}
		}
		iter.waveletTreeStructure->treeNodes[iter.position].i2 = 0;
	}

	//////////////////////////////////////////////////////////////////////////////////////////////
	template<typename TText >
	void goLeft(Iter< WaveletTreeStructure< TText >, Standard > &iter)
	{
		iter.position = getLeftChildPos(iter);
	}

	template<typename TText >
	void goLeft(Iter< const WaveletTreeStructure< TText >, Standard > &iter)
	{
		iter.position = getLeftChildPos(iter);
	}

	//////////////////////////////////////////////////////////////////////////////////////////////
	template<typename TText >
	void goRight(Iter< WaveletTreeStructure< TText >, Standard > &iter)
	{
		iter.position = getRightChildPos(iter);
	}

	template<typename TText >
	void goRight(Iter< const WaveletTreeStructure< TText >, Standard > &iter)
	{
		iter.position = getRightChildPos(iter);
	}

	//////////////////////////////////////////////////////////////////////////////////////////////
	template< typename TText, typename TPosInSubTree >
	void printTreeLevel(WaveletTreeStructure< TText > &tree, TPosInSubTree posInTree)
	{
		TPosInSubTree leftChild = tree.treeNodes[posInTree].i2;
		TPosInSubTree rightChild = tree.treeNodes[posInTree].i3;
		std::cout << "(" << posInTree << ": " << tree.treeNodes[posInTree].i1 << ", " << leftChild << ", " << rightChild << ") | ";
		if(leftChild)
		{
			printTreeLevel(tree, leftChild);
		}
		if(rightChild)
		{
			printTreeLevel(tree, rightChild);
		}
	}

	template< typename TText >
	void printTree(WaveletTreeStructure< TText > &tree)
	{
		printTreeLevel(tree, 0);
	}

	template< typename TText, typename TString >
	void writeGraphImpl(Iter< WaveletTreeStructure< TText >, Standard > &iter, TString name)
	{
		typedef typename BitVector_< BitsPerValue< typename Value< TText >::Type >::VALUE >::Type TValue;
		typename Iterator< WaveletTreeStructure< TText > >::Type iter2 = iter;
		std::ofstream stream(toCString(name), std::ios::app);
		TValue pos = getLeftChildPos(iter);
		if(pos)		
		{
			stream << ordValue(iter.waveletTreeStructure->treeNodes[iter.position].i1) << " -> " << ordValue(iter.waveletTreeStructure->treeNodes[pos].i1) << ";" << std::endl;
			goLeft(iter);
			writeGraphImpl(iter, name);
		}
		else
		{
			stream << ordValue(iter.waveletTreeStructure->treeNodes[iter.position].i1) << " -> " << "leave1" <<ordValue(iter.position)<< ";" << std::endl;
		}
	
		pos = getRightChildPos(iter2);
		if(pos)
		{
			stream << ordValue(iter2.waveletTreeStructure->treeNodes[iter2.position].i1) << " -> " << ordValue(iter2.waveletTreeStructure->treeNodes[pos].i1) << ";" << std::endl;
			goRight(iter2);
			writeGraphImpl(iter2, name);
		}
		else
		{
			stream << ordValue(iter2.waveletTreeStructure->treeNodes[iter2.position].i1) << " -> " << "leave2" <<ordValue(iter2.position)<< ";" << std::endl;
		}
	    stream.close();
	}

	template< typename TText >
	void writeGraph(WaveletTreeStructure< TText > &tree)
	{

		typename Iterator< WaveletTreeStructure< TText > >::Type iter(tree, 0);

		String< char > name = "testfile.dot";
		std::ofstream stream(toCString(name), std::ios::out);
		stream << "digraph G {" << std::endl;
		stream.close();
		writeGraphImpl(iter, name);
		
		stream.open(toCString(name), std::ios::app);
		stream << "}" << std::endl;
	    stream.close();
	}
	
	template< typename TText, typename TString >
	void writeGraph(WaveletTreeStructure< TText > &tree, TString name)
	{

		typename Iterator< WaveletTreeStructure< TText > >::Type iter(tree, 0);
		std::ofstream stream(toCString(name), std::ios::out);
		stream << "digraph G {" << std::endl;
		stream.close();
		writeGraphImpl(iter, name);
		
		stream.open(toCString(name), std::ios::app);
		stream << "}" << std::endl;
	    stream.close();
	}
	
	template< typename TStringLengthString, typename TPrefixSumTable, typename TIter, typename TTreeSplitValue >
	void computeSingleStringLengthFromTree(TStringLengthString &lengthString,
			TPrefixSumTable &prefixSumTable,
			TIter &iter,
			TTreeSplitValue lowerBound,
			TTreeSplitValue upperBound) 
	{

		lengthString[iter.position] = prefixSumTable[upperBound + 1] - prefixSumTable[lowerBound];
		TIter iter2 = iter;

		TTreeSplitValue newSplit = getCharacter(iter);
		goLeft(iter);
		if(getPosition(iter))
		{
			computeSingleStringLengthFromTree(lengthString,
				prefixSumTable,
		   		iter,
				lowerBound,
				newSplit); 
		}

		newSplit = getCharacter(iter2) + 1;
		goRight(iter2);
		if(getPosition(iter2))
		{
			computeSingleStringLengthFromTree(lengthString,
				prefixSumTable,
		   		iter2,
				newSplit,
				upperBound); 
		}
	}
	
	template< typename TStringLengthString, typename TPrefixSumTable, typename TText >
	void computeStringLengthFromTree(TStringLengthString &lengthString,
		TPrefixSumTable &prefixSumTable,
	   	WaveletTreeStructure< TText > &tree)
	{
		resize(lengthString, length(tree.treeNodes), 0);
		typename Iterator< WaveletTreeStructure< TText > >::Type iter(tree, 0);
		typedef typename BitVector_< BitsPerValue< typename Value< TText >::Type >::VALUE >::Type TValue;

		computeSingleStringLengthFromTree(lengthString,
				prefixSumTable,
		   		iter,
				(TValue)0,
				(TValue)(length(prefixSumTable)-2));

	}

/*
	template< typename TFreqString, 
		typename TText, 
		typename TTreeNodeWidthString,
	   	typename TValue,
		typename TPosInTree,
		typename TNumOfChildNodes >
	void computeTreeEntries_(
		TFreqString &freq,
		String< Pair<  > > &treeNodes,
		TTreeNodeWidthString &treeNodeWidth,
		TValue smallestValue,
		TValue biggestValue,
		TPosInTree posInTree,
		TNumOfChildNodes &numChildNodes)
	{	
		TValue oldSmallestValue = smallestValue;
		TValue oldBiggestValue = biggestValue;
		typedef typename Value< TFreqString >::Type TSize;
  		TSize leftSize = freq[smallestValue];
		TSize rightSize = freq[biggestValue]; 

		if(smallestValue == biggestValue)
		{
			return;
		}
		
		if(smallestValue == biggestValue - 1)
		{
			treeNodes[posInTree].i1 = smallestValue;
			treeNodes[posInTree].i2 = 0;
			treeNodes[posInTree].i3 = 0;
			treeNodeWidth[posInTree] = freq[smallestValue] + freq[biggestValue];
			++numChildNodes;
			return;
		}

		TValue leftCounter = 0;
		TValue rightCounter = 0;
		while(biggestValue - 1 > smallestValue)
		{
			if(leftSize < rightSize)
			{
				leftSize += freq[smallestValue + 1];
				++smallestValue;
				++leftCounter;
			}
			else if(leftSize > rightSize)
			{
				rightSize += freq[biggestValue - 1];
				--biggestValue;
				++rightCounter;
			}
			else
			{
				if(leftCounter < rightCounter)
				{
					leftSize += freq[smallestValue + 1];
					++smallestValue;
					++leftCounter;
				}
				else
				{
					rightSize += freq[biggestValue - 1];
					--biggestValue;
					++rightCounter;
				}
			}
		}

		treeNodeWidth[posInTree] = leftSize + rightSize;
		treeNodes[posInTree].i1 = smallestValue;
		
		computeTreeEntries_(freq, treeNodes, treeNodeWidth, oldSmallestValue, smallestValue, posInTree + 1, numChildNodes);
	//	treeNodes[posInTree].i2 = posInTree + 1;
		treeNodes[posInTree].i2 = posInTree + numChildNodes + 1;
		++numChildNodes;
		
		TNumOfChildNodes numChildNodes2 = 0;
		computeTreeEntries_(freq, treeNodes, treeNodeWidth, biggestValue, oldBiggestValue, treeNodes[posInTree].i3, numChildNodes2);
		numChildNodes += numChildNodes2;
		return;
	}
*/

	template< typename TFreqString, 
		typename TText, 
	   	typename TValue,
		typename TNumOfChildNodes >
	void computeTreeEntries(
		TFreqString &freq,
		Iter< WaveletTreeStructure< TText >, Standard > &iter,
		TValue smallestValue,
		TValue biggestValue,
		TNumOfChildNodes &numChildNodes)
	{	
		
		typedef typename BitVector_< BitsPerValue< typename Value< TText >::Type >::VALUE >::Type TTreeValue;

		TValue oldSmallestValue = smallestValue;
		TValue oldBiggestValue = biggestValue;
		typedef typename Value< TFreqString >::Type TSize;
  		TSize leftSize = freq[smallestValue];
		TSize rightSize = freq[biggestValue]; 

		if(smallestValue == biggestValue - 1)
		{
			setCharacter(iter, (TTreeValue)smallestValue);
			setNodeToLeaf(iter);
			++numChildNodes;
			return;
		}

		TValue leftCounter = 0;
		TValue rightCounter = 0;
		while(biggestValue - 1 > smallestValue)
		{
			if(leftSize < rightSize)
			{
				++smallestValue;
				leftSize += freq[smallestValue];
				++leftCounter;
			}
			else if(leftSize > rightSize)
			{
				--biggestValue;
				rightSize += freq[biggestValue];
				++rightCounter;
			}
			else
			{
				if(leftCounter < rightCounter)
				{
					++smallestValue;
					leftSize += freq[smallestValue];
					++leftCounter;
				}
				else
				{
					--biggestValue;
					rightSize += freq[biggestValue];
					++rightCounter;
				}
			}
		}

		setCharacter(iter, (TTreeValue)smallestValue);
		//if the elements of the node do not appear in the text make a leaf here
		if(!leftSize && !rightSize)
		{
			setNodeToLeaf(iter);
			++numChildNodes;
			return;
		}

		//there must be either a left or right subtree or both
		typename Iterator< WaveletTreeStructure< TText > >::Type iter2 = iter;
		
		//is there only one element on the left side?
		if(oldSmallestValue == smallestValue){
			setRightChildPos(iter, iter.position + numChildNodes + 1);
			++numChildNodes;
		}
		else
		{
			setLeftChildPos(iter);
			goLeft(iter);
			computeTreeEntries(freq, iter, oldSmallestValue, smallestValue, numChildNodes);
		}
		
		//is there only one element on the right side?
		if(biggestValue == oldBiggestValue){
			setLeftChildPos(iter2);
			++numChildNodes;
		}
		else
		{		
			setRightChildPos(iter2, iter2.position + numChildNodes + 1);
			TNumOfChildNodes numChildNodes2 = 0;
			goRight(iter2);
			computeTreeEntries(freq, iter2, biggestValue, oldBiggestValue, numChildNodes2);
			numChildNodes += numChildNodes2 + 1;
			return;
		}
	}

	template< typename TText >
	inline bool open(
		WaveletTreeStructure< TText > &structure, 
		const char *fileName,
		int openMode)
	{
		String<char> name;
		name = fileName;	append(name, ".treestruct");	open(getFibre(structure, FibreTreeNodes()), toCString(name), openMode);
		return true;
	}

	template < typename TText >
	inline bool open(
		WaveletTreeStructure< TText > &structure, 
		const char *fileName)
	{
		return open(structure, fileName, DefaultOpenMode< WaveletTree< TText, SingleString > >::VALUE);
	}

	template< typename TText >
	inline bool save(
		WaveletTreeStructure< TText > const &structure, 
		const char *fileName,
		int openMode)
	{
		String<char> name;
		name = fileName;	append(name, ".treestruct");	save(getFibre(structure, FibreTreeNodes()), toCString(name), openMode);
		return true;
	}

	template < typename TText >
	inline bool save(
		WaveletTreeStructure< TText > const &structure, 
		const char *fileName)
	{
		return save(structure, fileName, DefaultOpenMode< WaveletTree< TText, SingleString > >::VALUE);
	}
}

#endif
