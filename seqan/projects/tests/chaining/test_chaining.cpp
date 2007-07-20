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

template <typename TScoring>
void testChainer(int count, 
				 int dim,
				 TScoring scoring)
{
	String< Fragment<int> > fragments;
	reserve(fragments, count);
	
	_generateRandomFrags(fragments, count, 1, 3 * count, 10, 20, dim);

for (unsigned int i = 0; i < length(fragments); ++i)
{
	for (unsigned int j = 0; j < dim; ++j)
		_setLeftPosition(fragments[i], j, leftPosition(fragments[i], j) + 1);
}

	String< Fragment<int> > ch;
	reserve(ch, count);


	//build chainer chain
	int chain_score = chain(fragments, ch, scoring, Chainer());
	std::cout << chain_score << ", " << length(ch) << "\n";

for (unsigned int i = 0; i < length(fragments); ++i)
{
	for (unsigned int j = 0; j < dim; ++j)
		_setLeftPosition(fragments[i], j, leftPosition(fragments[i], j) - 1);
}
for (unsigned int i = 1; i < length(ch); ++i)
{
	for (unsigned int j = 0; j < dim; ++j)
		_setLeftPosition(ch[i], j, leftPosition(ch[i], j) - 1);
}

	//verify validity of chainer chain
	SEQAN_TASSERT(length(ch) > 0)
	int sum = weight(ch[0]);
	for (unsigned int i = 1; i < length(ch); ++i)
	{
		SEQAN_TASSERT(_chain_generic_chainable(ch[i-1], ch[i]))
		sum += scoreChainGap(scoring, ch[i-1], ch[i]) + weight(ch[i]);
	}

	//verify score of chainer chain
	std::cout << sum << "\n";
	SEQAN_TASSERT(sum == chain_score)

	//build generic chain
	int chain_score2 = chain(fragments, ch, scoring, GenericChaining());
	std::cout << chain_score2 << "\n";

	//verify validity of generic chain
	SEQAN_TASSERT(length(ch) > 0)
	sum = weight(ch[0]);
	for (unsigned int i = 1; i < length(ch); ++i)
	{
		SEQAN_TASSERT(_chain_generic_chainable(ch[i-1], ch[i]))
		sum += scoreChainGap(scoring, ch[i-1], ch[i]) + weight(ch[i]);
	}

	//verify score of generic chain
	std::cout << sum << "\n";
	SEQAN_TASSERT(sum == chain_score2)

	//compare it with generic chaining
	SEQAN_TASSERT(chain_score2 == chain_score)

/*
	clear(ch);
	chain_score = chain(fragments, ch, Score<int, Manhattan>(), Chainer());
	std::cout << chain_score << "\n";

	clear(ch);
	chain_score = chain(fragments, ch, Score<int, ChainSoP>(), Chainer());
	std::cout << chain_score << "\n";
*/
}

//////////////////////////////////////////////////////////////////////////////


int main()
{
//	testChainer(1000, 2, Score<int, Zero>());
	testChainer(1000, 2, Score<int, Manhattan>());

	return 0;
}

