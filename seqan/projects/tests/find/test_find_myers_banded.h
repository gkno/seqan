/*
   Helper functions for testing banded edit distance aproximate
   string search algorithms.
 */

using namespace std;
using namespace seqan;

//////////////////////////////////////////////////////////////////////////////

String<int> mat;

template <typename TMatrix, typename TSeq1, typename TSeq2>
void dumpMat(TMatrix &mat, TSeq1 &seq1, TSeq2 &seq2)
{
    // diagonal DP
    for(unsigned j=0;j<=length(seq2);++j) {
        std::cerr << "mat  ";
        for(int i=length(seq1);i>=0;--i) {
            int val = mat[   j  * (length(seq1)+1) +   i];
            if (val != -9)
                std::cerr << std::setw(5) << val;
            else
                std::cerr << std::setw(5) << ' ';
        }
        if (j > 0)
            std::cerr << std::setw(5) << seq2[j-1];
        std::cerr << std::endl;
    }
}

template <typename TString>
bool testMyersUkkonen(TString seq1, TString seq2, bool dump = true) 
{
	Finder<TString> finder(seq2);
//	Pattern<TString, MyersUkkonenBanded> pattern(seq1);
	PatternState_<TString,  Myers<AlignTextBanded<NMatchesN_,NMatchesN_>, True, void> > state;

	bool equal = true;
	int delta = length(seq2) - length(seq1);
	
	clear(mat);
	fill(mat, (length(seq1)+1) * (length(seq2)+1), -9, Exact());

    for (unsigned i = 0; i <= length(seq1); ++i)
        mat[i] = i;
    for (unsigned i = 0; i <= length(seq2); ++i)
        mat[i * (length(seq1) + 1)] = (i < length(seq1))? 0 : i - length(seq1) + 1;

    while (!atEnd(finder))
        find(finder, seq1, state, -1000);


#ifndef SEQAN_TEST_MYERS_STRICTBANDED

	delta = 7;

	// banded DP alignment
	for(unsigned j = 1; j <= length(seq2); ++j)
		for(unsigned i = 1; i <= length(seq1); ++i) 
        {
			int diag = j-i;
			if (0 <= diag && (diag <= delta || i+delta >= length(seq1)))
			{
				int d = mat[(j-1) * (length(seq1)+1) + i-1];
				if (seq1[i-1] != seq2[j-1]) ++d;
				int h = mat[(j-1) * (length(seq1)+1) +   i] + 1;
				int v = mat[   j  * (length(seq1)+1) + i-1] + 1;
				int min = 99;
				if (0 < diag)
					if (min > h) min = h;
				if (diag < delta || i+delta > length(seq1))
					if (min > v) min = v;
				if (diag <= delta || i+delta > length(seq1))
					if (min > d) min = d;
				mat[j * (length(seq1)+1) + i] = min;
			}
		}

#else

	if (delta >= length(seq1)) return true;

	// real banded DP alignment
	for(int j = 1; j <= length(seq2); ++j)
		for(int i = 1; i <= length(seq1); ++i) {
			int diag = j-i;
			if (0 <= diag && diag <= delta)
			{
				int d = mat[(j-1) * (length(seq1)+1) + i-1];
				if (seq1[i-1] != seq2[j-1]) ++d;
				int h = mat[(j-1) * (length(seq1)+1) +   i] + 1;
				int v = mat[   j  * (length(seq1)+1) + i-1] + 1;
				int min = 99;
				if (0 < diag) if (min > h) min = h;
				if (diag < delta) if (min > v) min = v;
				if (min > d) min = d;
				mat[j * (length(seq1)+1) + i] = min;
			}
		}

#endif

	for(unsigned j = 1; j <= length(seq2); ++j)
		for(unsigned i = 1; i <= length(seq1); ++i) 
        {
            unsigned pos = j * (length(seq1)+1) + i;
            if (mat[pos] != state.DPMat[pos])
                equal = false;
        }
    
    if (equal) return true;

	if (dump)
    {
        dumpMat(state.DPMat, seq1, seq2);
        dumpMat(mat, seq1, seq2);

		for(unsigned j = 1; j <= length(seq2); ++j)
			for(unsigned i = 1; i <= length(seq1); ++i) 
			{
				unsigned pos = j * (length(seq1)+1) + i;
				SEQAN_ASSERT_EQ(mat[pos], state.DPMat[pos]);
			}
    }    

	return equal;
}

template <typename TFinderCSP, typename TPatternCSP, typename TText, typename TNeedle>
void testCSPImpl(TText &text, TNeedle &needle, int errors)
{
	PatternState_<TNeedle, Myers<AlignTextBanded<TFinderCSP,TPatternCSP>, True, void> > state;
	Finder<TText> finder(text);
    SEQAN_ASSERT_TRUE(find(finder, needle, state, -1000));
    SEQAN_ASSERT_EQ(position(finder), length(text) - 1) ;
    SEQAN_ASSERT_EQ(_getMatchScore(state), -errors);
}

SEQAN_DEFINE_TEST(test_myers_find_banded_csp)
{
    {
        Dna5String text   = "ACAGTNNTAAGNNNNA";
        Dna5String needle = "ACNGTACTAAGNNNNG";
        testCSPImpl<NMatchesAll_,  NMatchesAll_>  (text, needle, 1);
        testCSPImpl<NMatchesAll_,  NMatchesN_>    (text, needle, 2);
        testCSPImpl<NMatchesN_,    NMatchesAll_>  (text, needle, 3);
        testCSPImpl<NMatchesN_,    NMatchesN_>    (text, needle, 4);
        testCSPImpl<NMatchesNone_, NMatchesNone_> (text, needle, 8);
    }
    {
        CharString text   = "ACAGTNNTAAGNNNNA";
        CharString needle = "ACNGTACTAAGNNNNG";
        testCSPImpl<NMatchesAll_,  NMatchesAll_>  (text, needle, 1);
        testCSPImpl<NMatchesAll_,  NMatchesN_>    (text, needle, 2);
        testCSPImpl<NMatchesN_,    NMatchesAll_>  (text, needle, 3);
        testCSPImpl<NMatchesN_,    NMatchesN_>    (text, needle, 4);
        testCSPImpl<NMatchesNone_, NMatchesNone_> (text, needle, 8);
    }
}

//////////////////////////////////////////////////////////////////////////////

SEQAN_DEFINE_TEST(test_myers_find_banded)
{
	typedef Dna TValue;
	String<Dna> seq1 = "tgtaaaggagt";
	String<Dna> seq2 = "tgtgtaaaggagttgtggagttgtaaaaaggagt";

	for(unsigned li=2; li<=length(seq1); ++li)
		for(unsigned lj=li; lj<=length(seq2) && lj-li<8; ++lj)
			for(unsigned i=0; i+li<=length(seq1); ++i) 
            if (li > 8)
			{
				String<TValue> s1 = infix(seq1, i, i+li);
				for(unsigned j=0; j+lj<=length(seq2); ++j) 
				{
					String<TValue> s2 = infix(seq2, j, j+lj);
					if (!testMyersUkkonen(s1, s2, false)) {
						std::cerr << "DIFFERENCE("<<right<<") for " << s1 << "," << s2 << std::endl;
						testMyersUkkonen(s1, s2, true);
					}
				}
			}
}
