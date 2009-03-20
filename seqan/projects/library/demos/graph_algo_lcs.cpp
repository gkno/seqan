///A tutorial about the longest common subsequence algorithm.
#include <iostream>
#include <seqan/graph_algorithms.h>

using namespace seqan;

int main() {
///Create two sequences
	String<char> seq1("abacx");
	String<char> seq2("baabca");
///Out-parameter: A string of positions belonging to the longest common subsequence
	String<std::pair<unsigned int, unsigned int>, Block<> > pos;
///Longest common subsequence
	longestCommonSubsequence(seq1, seq2, pos);
///Console ouptut
	::std::cout << seq1 << ::std::endl;
	::std::cout << seq2 << ::std::endl;
	::std::cout << "Lcs:" << ::std::endl;
	for(int i = length(pos)-1; i>=0; --i) {
		::std::cout << seq1[pos[i].first] <<  ',';
	}
	::std::cout << ::std::endl;
	for(int i = length(pos)-1; i>=0; --i) {
		::std::cout << seq2[pos[i].second] <<  ',';
	}
	::std::cout << ::std::endl;
	return 0;
}
