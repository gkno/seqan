#include <seqan/align.h>
#include <seqan/seeds.h>

using namespace seqan;
using namespace std;

int main() 
{
	typedef String<Dna>					TSequence;
	typedef Infix<TSequence>::Type		TInfix;
	typedef StringSet<TSequence>		TSequenceSet;
	typedef StringSet<TInfix>			TInfixSet;

	typedef int							TScoreValue;
	typedef Position<TSequence>::Type	TPos;
	typedef Seed<int, SimpleSeed>		TSeed;
	
	// read sequences
	TSequenceSet seqs;
	appendValue(seqs, String<Dna, FileReader<Fasta> >("../../../demos/sequence_1.fa"));
	appendValue(seqs, String<Dna, FileReader<Fasta> >("../../../demos/sequence_2.fa"));

	// set input infixes
	TInfixSet infixes;
	appendValue(infixes, TInfix(seqs[0], 2958, 3158));
	appendValue(infixes, TInfix(seqs[1], 3354, 3574));

	// set the scoring scheme
	TScoreValue penalty = -19;
	Score<TScoreValue> scoreMatrix(1, penalty, penalty);

	// align object for the complete sequences
	Align<TSequence> align;
	resize(rows(align), 2);
	assignSource(row(align, 0), seqs[0]);
	assignSource(row(align, 1), seqs[1]);

	// ---------- banded local alignment on infixes ----------
	Align<TInfix> localAlign;
    resize(rows(localAlign), 2);
    assignSource(row(localAlign, 0), infixes[0]);
    assignSource(row(localAlign, 1), infixes[1]);

	localAlignment(localAlign, scoreMatrix); // TODO: specify band as param3, parma4
	integrateAlign(align, localAlign);

	// ---------- seed extension with local alignment as seed ----------

	// begin and end position of local alignment in sequences
	TPos locAliBegin_0 = sourceBeginPosition(row(localAlign, 0)) + beginPosition(infixes[0]);
	TPos locAliBegin_1 = sourceBeginPosition(row(localAlign, 1)) + beginPosition(infixes[1]);
	TPos locAliEnd_0 = sourceEndPosition(row(localAlign, 0)) + beginPosition(infixes[0]);
	TPos locAliEnd_1 = sourceEndPosition(row(localAlign, 1)) + beginPosition(infixes[1]);

	TSeed seed(locAliBegin_0, locAliBegin_1, locAliEnd_0 - 1, locAliEnd_1 - 1);
	extendSeed(seed, 5*(-penalty), scoreMatrix, seqs[0], seqs[1], 2, GappedXDrop());

	// alignment on left extension
	if (leftPosition(seed, 0) < (int)locAliBegin_0 || leftPosition(seed, 1) < (int)locAliBegin_1) {
		TInfixSet leftExtension;
		appendValue(leftExtension, TInfix(seqs[0], leftPosition(seed, 0), locAliBegin_0));
		appendValue(leftExtension, TInfix(seqs[1], leftPosition(seed, 1), locAliBegin_1));

		Align<TInfix> leftAlign;
		resize(rows(leftAlign), 2);
		assignSource(row(leftAlign, 0), leftExtension[0]);
		assignSource(row(leftAlign, 1), leftExtension[1]);

		globalAlignment(leftAlign, scoreMatrix);
		integrateAlign(align, leftAlign);
	}

	// alignment on rigth extension
	if (rightPosition(seed, 0) + 1 > (int)locAliEnd_0 || rightPosition(seed, 1) + 1 > (int)locAliEnd_1) {
		TInfixSet rightExtension;
		appendValue(rightExtension, TInfix(seqs[0], locAliEnd_0, rightPosition(seed, 0) + 1));
		appendValue(rightExtension, TInfix(seqs[1], locAliEnd_1, rightPosition(seed, 1) + 1));

		Align<TInfix> rightAlign;
		resize(rows(rightAlign), 2);
		assignSource(row(rightAlign, 0), rightExtension[0]);
		assignSource(row(rightAlign, 1), rightExtension[1]);

		globalAlignment(rightAlign, scoreMatrix);
		integrateAlign(align, rightAlign);
	}

	// set extended begin and end positions of align
	setSourceBeginPosition(row(align, 0), leftPosition(seed, 0));
	setSourceBeginPosition(row(align, 1), leftPosition(seed, 1));
	setBeginPosition(row(align, 0), 0);
	setBeginPosition(row(align, 1), 0);
	setSourceEndPosition(row(align, 0), rightPosition(seed, 0) + 1);
	setSourceEndPosition(row(align, 1), rightPosition(seed, 1) + 1);

	std::cout << align;
	return 0;
}
