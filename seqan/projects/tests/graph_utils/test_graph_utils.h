#ifndef SEQAN_HEADER_TEST_GRAPH_UTILS_H
#define SEQAN_HEADER_TEST_GRAPH_UTILS_H

using namespace std;
using namespace seqan;

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////

void  Test_Alphabets() {
	// Test Rna
	Rna5 r = 'U';
	SEQAN_TASSERT(r == 3)
	r = Byte(0);
	SEQAN_TASSERT(r == 'a')
	char c = 'c';
	r = c;
	c = r;
	SEQAN_TASSERT(r == Rna5('c'))
	SEQAN_TASSERT((c == 'C' || c == 'c'))
	Unicode u = 'g';
	r = u;
	SEQAN_TASSERT(r == Rna5('g'))
	

	// Test Dayhoff
	AAGroupsDayhoff gr;
	gr = AminoAcid('T');
	SEQAN_TASSERT(gr == 0)
	gr = Byte(3);
	SEQAN_TASSERT(gr == 2)
	gr = char('c');
	SEQAN_TASSERT(gr == 5)
	gr = Unicode('j');
	SEQAN_TASSERT(gr == 6)

	// Test SeB6
	AAGroupsSeB6 gr1;
	gr1 = AminoAcid('T');
	SEQAN_TASSERT(gr1 == 0)
	gr1 = Byte(3);
	SEQAN_TASSERT(gr1 == 2)
	gr1 = char('c');
	SEQAN_TASSERT(gr1 == 1)
	gr1 = Unicode('j');
	SEQAN_TASSERT(gr1 == 6)

	// Test SeB8
	AAGroupsSeB8 gr2;
	gr2 = AminoAcid('T');
	SEQAN_TASSERT(gr2 == 0)
	gr2 = Byte(3);
	SEQAN_TASSERT(gr2 == 2)
	gr2 = char('c');
	SEQAN_TASSERT(gr2 == 1)
	gr2 = Unicode('j');
	SEQAN_TASSERT(gr2 == 8)

	// Test Murphy
	AAGroupsMurphy gr3;
	gr3 = AminoAcid('T');
	SEQAN_TASSERT(gr3 == 9)
	gr3 = Byte(3);
	SEQAN_TASSERT(gr3 == 2)
	gr3 = char('c');
	SEQAN_TASSERT(gr3 == 1)
	gr3 = Unicode('j');
	SEQAN_TASSERT(gr3 == 10)

	// Test SolisG10
	AAGroupsSolisG10 gr4;
	gr4 = AminoAcid('T');
	SEQAN_TASSERT(gr4 == 8)
	gr4 = Byte(3);
	SEQAN_TASSERT(gr4 == 2)
	gr4 = char('c');
	SEQAN_TASSERT(gr4 == 1)
	gr4 = Unicode('j');
	SEQAN_TASSERT(gr4 == 10)

	// Test SolisD10
	AAGroupsSolisD10 gr5;
	gr5 = AminoAcid('T');
	SEQAN_TASSERT(gr5 == 6)
	gr5 = Byte(3);
	SEQAN_TASSERT(gr5 == 2)
	gr5 = char('c');
	SEQAN_TASSERT(gr5 == 1)
	gr5 = Unicode('j');
	SEQAN_TASSERT(gr5 == 10)

	// Test LiB10
	AAGroupsLiB10 gr6;
	gr6 = AminoAcid('T');
	SEQAN_TASSERT(gr6 == 0)
	gr6 = Byte(3);
	SEQAN_TASSERT(gr6 == 2)
	gr6 = char('c');
	SEQAN_TASSERT(gr6 == 1)
	gr6 = Unicode('j');
	SEQAN_TASSERT(gr6 == 10)

	// Test LiA10
	AAGroupsLiA10 gr7;
	gr7 = AminoAcid('T');
	SEQAN_TASSERT(gr7 == 9)
	gr7 = Byte(3);
	SEQAN_TASSERT(gr7 == 1)
	gr7 = char('c');
	SEQAN_TASSERT(gr7 == 0)
	gr7 = Unicode('j');
	SEQAN_TASSERT(gr7 == 10)

	// Test SeV10
	AAGroupsSeV10 gr8;
	gr8 = AminoAcid('T');
	SEQAN_TASSERT(gr8 == 0)
	gr8 = Byte(3);
	SEQAN_TASSERT(gr8 == 2)
	gr8 = char('c');
	SEQAN_TASSERT(gr8 == 1)
	gr8 = Unicode('j');
	SEQAN_TASSERT(gr8 == 10)

	// Test SeB10
	AAGroupsSeB10 gr9;
	gr9 = AminoAcid('T');
	SEQAN_TASSERT(gr9 == 0)
	gr9 = Byte(3);
	SEQAN_TASSERT(gr9 == 2)
	gr9 = char('c');
	SEQAN_TASSERT(gr9 == 1)
	gr9 = Unicode('j');
	SEQAN_TASSERT(gr9 == 10)

	// Test SeB14
	AAGroupsSeB14 gr10;
	gr10 = AminoAcid('T');
	SEQAN_TASSERT(gr10 == 12)
	gr10 = Byte(3);
	SEQAN_TASSERT(gr10 == 2)
	gr10 = char('c');
	SEQAN_TASSERT(gr10 == 1)
	gr10 = Unicode('j');
	SEQAN_TASSERT(gr10 == 14)
}


//////////////////////////////////////////////////////////////////////////////


void Test_GraphUtils() {
	Test_Alphabets();

	debug::verifyCheckpoints("projects/library/seqan/graph_utils/graph_drawing.h");
	debug::verifyCheckpoints("projects/library/seqan/graph_utils/graph_utility_alphabets.h");
	debug::verifyCheckpoints("projects/library/seqan/graph_utils/graph_utility_parsing.h");
	debug::verifyCheckpoints("projects/library/seqan/graph_utils/graph_utility_match_parsing.h");
}


}

#endif

