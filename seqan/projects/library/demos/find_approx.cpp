#include <iostream>
#include <seqan/find.h>

using namespace seqan;
using namespace std;

///This program finds all occurrences of $CCT$ in $AACTTAACCTAA$ with $\leq 1$ error using the @Spec.MyersUkkonen@ approximate search algorithm.
int main() 
{
	String<char> haystk("AACTTAACCTAA");
	String<char> ndl("CCT");

	Finder<String<char> > fnd(haystk);
	Pattern<String<char>, MyersUkkonen> pat(ndl);
///The function @Function.setScoreLimit@ sets the limit score an occurrence must reach.
///Since the used scoring scheme is a distance measure (edit distance), all scores are negative.
///A score limit of $\geq -1$ therefore means a edit distance $\leq 1$.
	setScoreLimit(pat, -1);
	while (find(fnd, pat))
	{
		cout << position(fnd) << ": " << getScore(pat) << "\n";
///Note that @Function.position@ returns the end position of the found occurrence.
	}

	return 0;
}

