#include <iostream>
#include <cstdio>
#include <vector>

#define SEQAN_TEST

#include <seqan/basic.h>

void Main_TestGaps(); //test_align_gaps.cpp
void Main_TestAlign(); //test_align_align.cpp
void Main_TestLocalAlign();//test_align_local.cpp
void Main_TestMyers(); // test_align_myers.cpp

//////////////////////////////////////////////////////////////////////////////

int main() 
{
	SEQAN_TREPORT("TEST BEGIN")

	Main_TestLocalAlign();
	Main_TestGaps();
	Main_TestAlign();
	Main_TestMyers();




/*
	typedef Align< String<char>, ArrayGaps> TAlign;
	TAlign a;
	resize(a, 2);
	assignSource(row(a, 0), "abcdefg", 1, 6);
	assignSource(row(a, 1), "abfcxef", 0, 7);

	int score = globalAlignment(a, SimpleScore() );

	cout << a;

	Gaps<String<char >, ArrayGaps > gps;
	assignSource(gps, "hallo");

	//Iterator<Gaps<String<char >, ArrayGaps > >::Type it = iter(gps, 0);
	Iter<Gaps<String<char >, ArrayGaps>, GapsIterator<void> > it = iter(gps, 0);
	++it;
	cout << *it << "\n";
*/
	SEQAN_TREPORT("TEST END");

	return 0;
}

/*
int main()
{
	Gaps<String<char>, ArrayGaps > gaps1;
	Gaps<String<char>, ArrayGaps > gaps2;
	String<char> str1 = "abcdefghijklmn";
	setSource(gaps1, str1, 2, 12);
	setSource(gaps2, str1, 0, 10);
	insertGaps(gaps1, 1, 1);
	insertGaps(gaps1, 5, 3);
	insertGaps(gaps1, 11, 2);
	cout << gaps1 << "\n";
	cout << gaps2 << "\n";

	copyGaps(gaps2, 2, gaps1, 1, 6);

	cout << gaps1 << "\n";
	cout << gaps2 << "\n";


	Align<String<char>, ArrayGaps > ali;
	resize(rows(ali), 2);
	row(ali, 0) = gaps1;
	row(ali, 1) = gaps2;

	SEQAN_TASSERT(id(row(ali, 0)) == id(gaps1));
	detach(ali);
	SEQAN_TASSERT(id(row(ali, 0)) != id(gaps1));

}
*/
