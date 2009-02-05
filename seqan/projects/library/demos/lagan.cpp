#include <fstream>
#include <seqan/seeds.h>
#include <seqan/file.h>

using namespace seqan;
using namespace std;

template <typename TSeed, typename TSegment, typename TSize, typename TScoring, typename TScore>
void laganChaining(list<TSeed> & chain,
				   TSegment const & a,
				   TSegment const & b,
				   TSize q,
				   TSize q_min,
				   TSize limit,
				   TSize gaps_max,
				   TScoring scoring_scheme,
				   TScore score_min)
{
	if (((TSize)length(a) <= gaps_max) && ((TSize)length(b) <= gaps_max)) return;

	//Step 1: find seeds
	typedef typename Value<TSeed>::Type TPosition;
	typedef SeedSet<TPosition, SimpleSeed, DefaultScore> TSeedSet;
	TSeedSet seedset(limit, score_min, scoring_scheme);

	typedef Index< TSegment, Index_QGram<SimpleShape > > TQGramIndex;
	TQGramIndex index_qgram(b);

	typedef Finder<TQGramIndex> TFinder;
	TFinder finder(index_qgram);

	while (length(seedset) == 0)
	{
		if (q <= q_min) return;

		resize(indexShape(index_qgram), q);
		for (int i = 0; i < length(a)-q+1; ++i)
		{
			while (find(finder, infix(a, i, i+q)))
			{
				typedef typename Position<TFinder>::Type TPosition;
				TPosition a_pos = beginPosition(a)+i;
				TPosition b_pos = beginPosition(b)+position(finder);
				if (!addSeed(seedset, a_pos, b_pos, q, 0, Merge()))
				if (!addSeed(seedset, a_pos, b_pos, q, host(a), host(b), 5, Chaos()))
					 addSeed(seedset, a_pos, b_pos, q, Single());
			}
			clear(finder);
		}
		--q;
	}

	//Step 2: global chaining
	globalChaining(seedset, chain);
	clear(seedset);

	//Step 3: recursively fill gaps
	if (q > q_min)
	{
		list<TSeed> subchain;
		typedef typename list<TSeed>::iterator TIterator;

		TIterator it = chain.begin();
		TIterator it2 = it; 
		++it2;

		laganChaining(subchain, 
			infix(host(a), beginPosition(a), leftDim0(*it)), 
			infix(host(b), beginPosition(b), leftDim1(*it)), 
			q, q_min, limit, gaps_max, scoring_scheme, score_min);
		chain.splice(it, subchain);

		while(it2 != chain.end())
		{
			laganChaining(subchain, 
				infix(host(a), rightDim0(*it), leftDim0(*it2)), 
				infix(host(b), rightDim1(*it), leftDim1(*it2)), 
				q, q_min, limit, gaps_max, scoring_scheme, score_min);
			chain.splice(it2, subchain);

			it = it2;
			++it2;
		}

		laganChaining(subchain, 
			infix(host(a), rightDim0(*it), endPosition(a)), 
			infix(host(b), rightDim1(*it), endPosition(b)),
			q, q_min, limit, gaps_max, scoring_scheme, score_min);
		chain.splice(it2, subchain);
	}
}

int main()
{
	//load sequences
	typedef String<Dna> TString;
	TString a = String<Dna, FileReader<Fasta> >("lagan1.fasta");
	TString b = String<Dna, FileReader<Fasta> >("lagan2.fasta");

	//call lagan
	if ((length(a) > 0) && (length(b) > 0))
	{
		typedef Seed<int, SimpleSeed> TSeed;
		list<TSeed> chain;

		//Step 1 to 3
		laganChaining(chain, 
			infix(a, 0, length(a)), 
			infix(b, 0, length(b)),
			13, 8, 6, 200, SimpleScore(3,-2,-1,-3), 30);

		//Step 4: banded alignment
		Align<TString, ArrayGaps> alignment;
		resize(rows(alignment), 2);
		setSource(row(alignment, 0), a);
		setSource(row(alignment, 1), b);
		int score = bandedChainAlignment(chain, 7, alignment, SimpleScore(3,-2,-1,-3));

		cout << "Score: " << score << endl;
		cout << alignment << endl;
	}
	else
	{
		cout << "Error - File problem" << endl;
	}
	return 0;
}
