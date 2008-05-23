/// This code example illustrates a consistency-based multiple sequence alignment using an amino acid alphabet
#include <seqan/graph_msa.h>
#include <iostream>

using namespace seqan;


int main() {
	typedef String<AminoAcid> TString;
	typedef StringSet<TString, Dependent<> > TStringSet;
	typedef Graph<Alignment<TStringSet, unsigned int, Default> > TGraph;
/// Alignments are carried out on a StringSet that holds the sequences
	TString str1 = "GARFIELDTHELASTFATCAT";
	TString str2 = "GARFIELDTHEFASTCAT";
	TString str3 = "GARFIELDTHEVERYFASTCAT";
	TString str4 = "THEFATCAT";
	TStringSet strSet;
	assignValueById(strSet, str1);
	assignValueById(strSet, str2);
	assignValueById(strSet, str3);
	assignValueById(strSet, str4);
/// Out-parameter: An alignment graph of multiple sequences
	Graph<Alignment<TStringSet, void, WithoutEdgeId> > gOut(strSet);
/// Consistency-based multiple sequence alignment
	globalAlignment(strSet, gOut, MSA_Protein() );
/// Console output
	std::cout << gOut << std::endl;
	return 0;
}
