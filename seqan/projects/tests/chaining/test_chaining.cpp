#include <vector>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <cmath>

//#define SEQAN_DEBUG
#define SEQAN_TEST
#define SEQAN_NOSRAN //suppress srand

#include <seqan/sequence.h>
#include <seqan/chaining.h>

using namespace seqan;
using namespace std;

//////////////////////////////////////////////////////////////////////////////

template< typename TContainer >
void
_generateRandomFrags(TContainer & dest,
					 int num,
					 int min,
					 int max,
					 int minwidth,
					 int maxwidth,
					 int dim )
{
	typedef typename Value< TContainer >::Type FragType;
	typename Key< FragType >::Type * left_pos;
	typename Key< FragType >::Type * right_pos;
	reserve( dest, num );
	double d_dim = static_cast< double >( dim );
	for( int i = 0; i < num; ++i )
	{
		
		double width_sum = 0;
		left_pos = new typename Key< FragType >::Type[ dim ];
		right_pos = new typename Key< FragType >::Type[ dim ];
		for( int d = 0; d < dim; ++d )
		{
			left_pos[ d ] = ( rand() % ( max - min ) ) + min;
			int width = ( rand() % ( maxwidth - minwidth ) )+ minwidth;
			width_sum += width;
			right_pos[ d ] = left_pos[ d ] + width;			
			
		}
		FragType frag( left_pos, right_pos, dim, static_cast< typename Weight< FragType >::Type >(exp(log(width_sum)/d_dim)) * 100 );

		delete[] left_pos;
		delete[] right_pos;

		appendValue( dest, frag );
	}
}

//////////////////////////////////////////////////////////////////////////////
//helper function

template <typename TChain, typename TScoring>
void _showChain(TChain & ch,
				TScoring scoring)
{
	typedef typename Iterator<TChain>::Type TIterator;
	typedef typename Value<TChain>::Type TFragment;

	TIterator it = begin(ch);
	TIterator it_end = end(ch);
	unsigned int dim = dimension(*it); //chain must not be empty!
	while (it < it_end)
	{
		TFragment & frag = *it;
		printf("%4i: ", weight(frag));
		cout << "(";
		for (unsigned int i = 0; i < dim; ++i)
		{
			printf("%2i", leftPosition(frag, i));
			if (i < (dim-1)) cout << ", ";
		}
		cout << ")(";
		for (unsigned int i = 0; i < dim; ++i)
		{
			printf("%2i", rightPosition(frag, i));
			if (i < (dim-1)) cout << ", ";
		}
		cout << ")\n";

		++it;
		if (it == it_end) break;

		printf("%4i\n", scoreChainGap(scoring, frag, *it));
	}
}

//////////////////////////////////////////////////////////////////////////////

template <typename TScoring>
void testChainer(int count, 
				 int dim,
				 TScoring scoring)
{
	String< Fragment<int> > fragments;
	reserve(fragments, count);
	
	_generateRandomFrags(fragments, count, 1, 3 * count, 10, 20, dim);
//_showChain(fragments, scoring);

	String< Fragment<int> > ch;
	reserve(ch, count+2);


	//build chain
	int chain_score = chain(fragments, ch, scoring);
std::cout << chain_score << "\n";
//_showChain(ch, scoring);

	//verify validity of chain
	SEQAN_TASSERT(length(ch) > 0)
	int sum = weight(ch[0]);
	for (unsigned int i = 1; i < length(ch); ++i)
	{
		SEQAN_TASSERT(_chain_generic_chainable(ch[i-1], ch[i]))
		sum += scoreChainGap(scoring, ch[i-1], ch[i]) + weight(ch[i]);
	}
	//verify score of chain
	SEQAN_TASSERT(sum == chain_score)

	//build generic chain
	int chain_score2 = chain(fragments, ch, scoring, GenericChaining());
std::cout << chain_score2 << "\n";
//_showChain(ch, scoring);

	//verify validity of generic chain
	SEQAN_TASSERT(length(ch) > 0)
	sum = weight(ch[0]);
	for (unsigned int i = 1; i < length(ch); ++i)
	{
		SEQAN_TASSERT(_chain_generic_chainable(ch[i-1], ch[i]))
		sum += scoreChainGap(scoring, ch[i-1], ch[i]) + weight(ch[i]);
	}
	//verify score of generic chain
	SEQAN_TASSERT(sum == chain_score2)

	//compare results of two chaining algorithms
	SEQAN_TASSERT(chain_score2 == chain_score)
}

//////////////////////////////////////////////////////////////////////////////

void createTest()
{//Problembeispiel mit neuer Fragmentsemantik
	Score<int, ChainSoP> scoring;

	String< Fragment<int> > frags;

	Fragment<int> f(2);

	setWeight(f, 100);

	_setLeftPosition(f, 0, 28);
	_setLeftPosition(f, 1, 18);
	_setRightPosition(f, 0, 29);
	_setRightPosition(f, 1, 23);
	appendValue(frags, f); //+

	_setLeftPosition(f, 0, 6);
	_setLeftPosition(f, 1, 4);
	_setRightPosition(f, 0, 19);
	_setRightPosition(f, 1, 17);
	appendValue(frags, f); //!

	_setLeftPosition(f, 0, 12);
	_setLeftPosition(f, 1, 4);
	_setRightPosition(f, 0, 26);
	_setRightPosition(f, 1, 17);
	appendValue(frags, f); //!


	String< Fragment<int> > ch;
	reserve(ch, 5);
	int chain_score = chain(frags, ch, scoring);
	std::cout << chain_score << ", " << length(ch) << "\n";
	_showChain(ch, scoring);

	chain_score = chain(frags, ch, scoring, GenericChaining());
	std::cout << chain_score << ", " << length(ch) << "\n";
	_showChain(ch, scoring);
}

void createTest2()
{//Problembeispiel mit alter Fragmentsemantik
	Score<int, ChainSoP> scoring;
	String< Fragment<int> > frags;

	Fragment<int> f(2);

	setWeight(f, 100);

	_setLeftPosition(f, 0, 7);
	_setLeftPosition(f, 1, 5);
	_setRightPosition(f, 0, 19);
	_setRightPosition(f, 1, 17);
	appendValue(frags, f); //fragment 1

	_setLeftPosition(f, 0, 13);
	_setLeftPosition(f, 1, 5);
	_setRightPosition(f, 0, 26);
	_setRightPosition(f, 1, 17);
	appendValue(frags, f); //fragment 2

	_setLeftPosition(f, 0, 29);
	_setLeftPosition(f, 1, 19);
	_setRightPosition(f, 0, 29);
	_setRightPosition(f, 1, 23);
	appendValue(frags, f); //fragment 3


	String< Fragment<int> > ch;
	reserve(ch, 5);
	int chain_score = compute_chain(frags, ch, scoring, SemiDeferred());
	std::cout << chain_score << "\n";
	_showChain(ch, scoring);
}

//////////////////////////////////////////////////////////////////////////////

int main()
{
/*
	createTest2();
	return 0;

	srand(787);
	testChainer(10, 2, Score<int, ChainSoP>());

	for (int i = 0; i < 1000000; ++i)
	{
		cout << i << "\n";
		srand(i);
//	testChainer(1000, 2, Score<int, Zero>());
//	testChainer(1000, 2, Score<int, Manhattan>());
		testChainer(10, 2, Score<int, ChainSoP>());
	}
*/
	testChainer(1000, 2, Score<int, Zero>());
	testChainer(1000, 2, Score<int, Manhattan>());
	testChainer(1000, 2, Score<int, ChainSoP>());

	return 0;
}

