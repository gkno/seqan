#include <seqan/seeds.h>
#include <seqan/file.h>
#include <fstream>

using namespace seqan;
using namespace std;

template <typename TSeedSet, typename TSegment, typename TSize>
TSize
findSeeds(TSeedSet & seedset,
		  TSegment const & seg_0,
		  TSegment const & seg_1,
		  TSize q,
		  TSize q_min)
{
	//add all common q-grams into seetset

	typedef Index< TSegment, Index_QGram<SimpleShape > > TQGramIndex;
	TQGramIndex index_qgram(seg_1);

	while ((length(seedset) == 0) && (q >= q_min))
	{
		resize(indexShape(index_qgram), q);

		typedef Finder<TQGramIndex> TFinder;
		TFinder finder(index_qgram);

		int q_gram_count = length(seg_0)-q+1;
		for (int i = 0; i < q_gram_count; ++i)
		{
			while (find(finder, infix(seg_0, i, i+q)))
			{
				typedef typename Position<TFinder>::Type TPosition;
				TPosition pos = position(finder);
				if (!addSeed(seedset, beginPosition(seg_0)+i, beginPosition(seg_1)+pos, q, 0, Merge()))
				if (!addSeed(seedset, beginPosition(seg_0)+i, beginPosition(seg_1)+pos, q, host(seg_0), host(seg_1), 5, Chaos()))
					 addSeed(seedset, beginPosition(seg_0)+i, beginPosition(seg_1)+pos, q, Single());
			}
			clear(finder);
		}
		--q;
	}
	return q;
}

template <typename TSeed, typename TSegment, typename TSize, typename TScoring, typename TScore>
void
lagan(list<TSeed> & chain,
	  TSegment const & seg_0,
	  TSegment const & seg_1,
	  TSize q,
	  TSize q_min,
	  TSize limit,
	  TSize gaps_max,
	  TScoring scoring_scheme,
	  TScore score_min)
{
	if (((TSize)length(seg_0) <= gaps_max) && ((TSize)length(seg_1) <= gaps_max)) return;

	//Step 1: find seeds
	typedef typename Value<TSeed>::Type TPosition;
	typedef SeedSet<TPosition, SimpleSeed, DefaultScore> TSeedSet;
	TSeedSet seedset(limit, score_min, scoring_scheme);
	q = findSeeds(seedset, seg_0, seg_1, q, q_min);
	if (!length(seedset)) return;

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

		lagan(subchain, 
			infix(host(seg_0), beginPosition(seg_0), leftDim0(*it)), 
			infix(host(seg_1), beginPosition(seg_1), leftDim1(*it)), 
			q, q_min, limit, gaps_max, scoring_scheme, score_min);
		chain.splice(it, subchain);

		while(it2 != chain.end())
		{
			lagan(subchain, 
				infix(host(seg_0), rightDim0(*it), leftDim0(*it2)), 
				infix(host(seg_1), rightDim1(*it), leftDim1(*it2)), 
				q, q_min, limit, gaps_max, scoring_scheme, score_min);
			chain.splice(it2, subchain);

			it = it2;
			++it2;
		}

		lagan(subchain, 
			infix(host(seg_0), rightDim0(*it), endPosition(seg_0)), 
			infix(host(seg_1), rightDim1(*it), endPosition(seg_1)),
			q, q_min, limit, gaps_max, scoring_scheme, score_min);
		chain.splice(it2, subchain);
	}
}

template <typename TString>
void 
loadFile(char const * filename,
		 TString & str)
{
	//load file into str
	fstream fstrm;
	fstrm.open(filename, ios_base::in | ios_base::binary);
	if (fstrm.fail())
	{
		cout << "error: " << filename << " not found" << endl;
		return;
	}
    read(fstrm, str, Fasta());
    fstrm.close(); 
}

int
main()
{
	//load sequences
	typedef String<Dna> TString;
	TString str_0;
	loadFile("lagan1.fasta", str_0);
	TString str_1;
	loadFile("lagan2.fasta", str_1);

	//call lagan
	if ((length(str_0) > 0) && (length(str_1) > 0))
	{
		typedef Seed<int, SimpleSeed> TSeed;
		list<TSeed> chain;

		//Step 1 to 3
		lagan(chain, 
			infix(str_0, 0, length(str_0)), 
			infix(str_1, 0, length(str_1)),
			13, 8, 6, 200, SimpleScore(3,-2,-1,-3), 30);

		//Step 4: banded alignment
		Align<TString, ArrayGaps> alignment;
		resize(rows(alignment), 2);
		setSource(row(alignment, 0), str_0);
		setSource(row(alignment, 1), str_1);
		int score = bandedChainAlignment(chain, 7, alignment, SimpleScore(3,-2,-1,-3));

		//cout << "Score: " << score << endl;
		cout << alignment << endl;
	}
	else
	{
		cout << "Error - File problem" << endl;
	}
}
