#include <iostream>
#include <seqan/seeds.h>

using namespace seqan;
using namespace std;

int main()
{
///Example 1: three algorithms for seed extension
	String<char> a = "SEEDabXcdXefXXX";
	String<char> b = "SEEDabYcdefYYYY";
	Seed<> seed1(0, 0, 4);          //left=0; length=4
	extendSeed(seed1, a, b, 1, MatchExtend());
	cout << rightPosition(seed1, 0) << endl;  //output: 6
	cout << rightPosition(seed1, 1) << endl;  //output: 6

	Seed<> seed2(0, 0, 4);          //left=0; length=4
	Score<> scoring(1, -1, -1);
	extendSeed(seed2, 2, scoring, a, b, 1, UngappedXDrop());
	cout << rightPosition(seed2, 0) << endl;  //output: 9
	cout << rightPosition(seed2, 1) << endl;  //output: 9

	Seed<> seed3(0, 0, 4);          //left=0; length=4
	extendSeed(seed3, 2, scoring, a, b, 1, GappedXDrop());
	cout << rightPosition(seed3, 0) << endl;  //output: 12
	cout << rightPosition(seed3, 1) << endl;  //output: 11

///Example 2: global chaining
	String< Seed<int, MultiSeed> > fragments;
	for (int i = 0; i < 1000; ++i)
	{
		Seed<int, MultiSeed> seed(3);
		for (int d = 0; d < 3; ++d)
		{
			int pos = rand();
			setLeftPosition(seed, d, pos);
			setRightPosition(seed, d, pos+100);
		}
		setWeight(seed, 200);
		appendValue(fragments, seed); 
	}
	String< Seed<int, MultiSeed> > global_chain;
	Score<int, Manhattan> scoring2;
	int chain_score = globalChaining(fragments, global_chain, scoring2);
	std::cout << chain_score << "\n";

	return 0;
}
