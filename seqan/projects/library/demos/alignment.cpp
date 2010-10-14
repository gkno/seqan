///A tutorial about global alignments.
#include <iostream>
#include <seqan/align.h>
#include <seqan/graph_align.h>

#include <seqan/seeds.h>

using namespace seqan;

//////////////////////////////////////////////////////////////////////////////
// insert gaps into alignment accoring to the traceback stored in trace
template <typename TSource, typename TSpec, typename TPosOffset, typename TTrace> 
void
_traceToAlign(Align<TSource, TSpec> & align_,
			  TPosOffset offset0,
			  TPosOffset offset1,
			  TTrace const & trace)
{
	typedef Align<TSource, TSpec> TAlign;
	typedef typename Size<TAlign>::Type TSize;

	typedef typename Row<TAlign>::Type TRow;
	typedef typename Iterator<TRow>::Type TRowIterator;

	TSize i = length(trace.sizes); //scan trace backwards
	TRowIterator it0 = begin(row(align_, 0)) + offset0;
	TRowIterator it1 = begin(row(align_, 1)) + offset1;
	while (i > 0)
	{
		--i;
		TSize siz = trace.sizes[i];
		switch ((int) trace.tvs[i])
		{
		case 1: //horizontal:
			insertGaps(it1, siz);
			break;

		case 2: //vertical:
			insertGaps(it0, siz);
			break;
		}
		goFurther(it0, siz);
		goFurther(it1, siz);
	}
}

//void birtesTestStichAlign() {
//	DnaString seq1 = "aaaaaaaagtcacagggactaatttttt";
//	DnaString seq2 = "cccccccgtcagtggaccctaagggggg";
//
//	std::cout << "Sequences:" << std::endl;
//	std::cout << seq1 << std::endl;
//	std::cout << seq2 << std::endl;
//	std::cout << std::endl;
//
//	Align<DnaString> align;
//	resize(rows(align), 2);
//	assignSource(row(align, 0), seq1);
//	assignSource(row(align, 1), seq2);
//
//	Score<int> scoring(1, -2, -2);
//
//	StringSet<Infix<DnaString>::Type> infixes;
//	appendValue(infixes, infix(seq1, 3, 23));
//	appendValue(infixes, infix(seq2, 3, 23));
//
//	localAlignment(align, infixes, LocalAlignmentFinder<>(), scoring, 4, -5, 5, BandedWatermanEggert());
//
//	std::cout << "local alignment:" << std::endl;
//	std::cout << align;
//	std::cout << std::endl;
//
//    int seedBegin1 = sourceBeginPosition(row(align, 0));
//    int seedBegin2 = sourceBeginPosition(row(align, 1));
//    int seedEnd1 = sourceEndPosition(row(align, 0));
//    int seedEnd2 = sourceEndPosition(row(align, 1));
//	Seed<int> seed(seedBegin1, seedBegin2, seedEnd1-1, seedEnd2-1);
//
//	std::cout << "Seed before extension:" << std::endl;
//	std::cout << infix(seq1, leftPosition(seed, 0), rightPosition(seed, 0)+1) << std::endl;
//	std::cout << infix(seq2, leftPosition(seed, 1), rightPosition(seed, 1)+1) << std::endl;
//	std::cout << std::endl;
//
//    extendSeed(seed, 4, scoring, seq1, seq2, 2, GappedXDrop());
//
//	std::cout << "Seed after extension:" << std::endl;
//	std::cout << infix(seq1, leftPosition(seed, 0), rightPosition(seed, 0)+1) << std::endl;
//	std::cout << infix(seq2, leftPosition(seed, 1), rightPosition(seed, 1)+1) << std::endl;
//	std::cout << std::endl;
//
//	setSourceBeginPosition(row(align, 0), leftPosition(seed, 0));
//	setSourceBeginPosition(row(align, 1), leftPosition(seed, 1));
//	setSourceEndPosition(row(align, 0), rightPosition(seed, 0)+1);
//	setSourceEndPosition(row(align, 1), rightPosition(seed, 1)+1);
//
//	std::cout << "After setting new source begin and end positions:" << std::endl;
//	std::cout << align;
//
//	StringSet<Infix<DnaString>::Type> infixStrLeft;
//	appendValue(infixStrLeft, infix(seq1, leftPosition(seed, 0), seedBegin1));
//	appendValue(infixStrLeft, infix(seq2, leftPosition(seed, 1), seedBegin2));
//
//	std::cout << "Left extension seq1: " << infix(seq1, leftPosition(seed, 0), seedBegin1) << std::endl;
//	std::cout << "Left extension seq2: " << infix(seq2, leftPosition(seed, 1), seedBegin2) << std::endl;
//	std::cout << std::endl;
//
//	_Align_Traceback<unsigned> trace;
//	globalAlignment(trace, infixStrLeft, scoring, -3, 3, BandedNeedlemanWunsch());
//	_traceToAlign(align,
//		toViewPosition(row(align, 0), leftPosition(seed, 0)) - beginPosition(row(align, 0)),
//		toViewPosition(row(align, 1), leftPosition(seed, 1)) - beginPosition(row(align, 1)),
//		trace);
//
//	std::cout << "After banded alignment on left extension:" << std::endl;
//	std::cout << align;
//
//	StringSet<Infix<DnaString>::Type> infixStrRight;
//	appendValue(infixStrRight, infix(seq1, seedEnd1, rightPosition(seed, 0)+1));
//	appendValue(infixStrRight, infix(seq2, seedEnd2, rightPosition(seed, 1)+1));
//
//	std::cout << "Right extension seq1: " << infix(seq1, seedEnd1, rightPosition(seed, 0)+1) << std::endl;
//	std::cout << "Right extension seq2: " << infix(seq2, seedEnd2, rightPosition(seed, 1)+1) << std::endl;
//	std::cout << std::endl;
//
//	clear(trace.sizes); clear(trace.tvs);
//	globalAlignment(trace, infixStrRight, scoring, -3, 3, BandedNeedlemanWunsch());
//	_traceToAlign(align, toViewPosition(row(align, 0), seedEnd1)-beginPosition(row(align, 0)), toViewPosition(row(align, 1), seedEnd2)-beginPosition(row(align, 1)), trace);
//
//	std::cout << "After banded alignment on right extension:" << std::endl;
//	std::cout << align;
//}

void birtesTestStichAlign2() {
	DnaString seq1 = "aaaaaaaagtcacagggactaattttttccatacagccttgctaaaaaaa";
	DnaString seq2 = "cccccccgtcagtggaccctaaggggggccagataagccgctccccc";

	std::cout << "Sequences:" << std::endl;
	std::cout << seq1 << std::endl;
	std::cout << seq2 << std::endl;
	std::cout << std::endl;

	Align<DnaString> align;
	resize(rows(align), 2);
	assignSource(row(align, 0), seq1);
	assignSource(row(align, 1), seq2);

	Score<int> scoring(1, -2, -2);

	Align<Infix<DnaString>::Type > infixAlign;
	resize(rows(infixAlign), 2);
	assignSource(row(infixAlign, 0), infix(seq1, 3, 45));
	assignSource(row(infixAlign, 1), infix(seq2, 3, 45));

	LocalAlignmentFinder<> finder = LocalAlignmentFinder<>();
	while (localAlignment(infixAlign, finder, scoring, 4, -5, 5, BandedWatermanEggert())) {

		std::cout << "local alignment:" << std::endl;
		std::cout << infixAlign;
		std::cout << std::endl;

		int seedBegin1 = sourceBeginPosition(row(infixAlign, 0)) + beginPosition(source(row(infixAlign, 0)));
		int seedBegin2 = sourceBeginPosition(row(infixAlign, 1)) + beginPosition(source(row(infixAlign, 1)));
		int seedEnd1 = sourceEndPosition(row(infixAlign, 0)) + beginPosition(source(row(infixAlign, 0)));
		int seedEnd2 = sourceEndPosition(row(infixAlign, 1)) + beginPosition(source(row(infixAlign, 1)));

		String<unsigned> pos;
		resize(pos, 2);
		value(pos, 0) = seedBegin1;
		value(pos, 1) = seedBegin2;
		integrateAlign(align, infixAlign, pos);
		std::cout << "align with local alignment:" << std::endl;
		std::cout << align;
		std::cout << std::endl;

		Seed<int> seed(seedBegin1, seedBegin2, seedEnd1-1, seedEnd2-1);

		std::cout << "Seed before extension:" << std::endl;
		std::cout << infix(seq1, leftPosition(seed, 0), rightPosition(seed, 0)+1) << std::endl;
		std::cout << infix(seq2, leftPosition(seed, 1), rightPosition(seed, 1)+1) << std::endl;
		std::cout << std::endl;

		extendSeed(seed, 4, scoring, seq1, seq2, 2, GappedXDrop());

		std::cout << "Seed after extension:" << std::endl;
		std::cout << infix(seq1, leftPosition(seed, 0), rightPosition(seed, 0)+1) << std::endl;
		std::cout << infix(seq2, leftPosition(seed, 1), rightPosition(seed, 1)+1) << std::endl;
		std::cout << std::endl;

		setSourceBeginPosition(row(align, 0), leftPosition(seed, 0));
		setSourceBeginPosition(row(align, 1), leftPosition(seed, 1));
		setBeginPosition(row(align, 0), 0);
		setBeginPosition(row(align, 1), 0);
		setSourceEndPosition(row(align, 0), rightPosition(seed, 0)+1);
		setSourceEndPosition(row(align, 1), rightPosition(seed, 1)+1);

		std::cout << "After setting new source begin and end positions:" << std::endl;
		std::cout << align;

		std::cout << "Left extension seq1: " << infix(seq1, leftPosition(seed, 0), seedBegin1) << std::endl;
		std::cout << "Left extension seq2: " << infix(seq2, leftPosition(seed, 1), seedBegin2) << std::endl;
		std::cout << std::endl;

		StringSet<Infix<DnaString>::Type> infixStrLeft;
		appendValue(infixStrLeft, infix(seq1, leftPosition(seed, 0), seedBegin1));
		appendValue(infixStrLeft, infix(seq2, leftPosition(seed, 1), seedBegin2));

		_Align_Traceback<unsigned> trace;
		globalAlignment(trace, infixStrLeft, scoring, -3, 3, BandedNeedlemanWunsch());

		Align<Infix<DnaString>::Type > infixAlign2;
		resize(rows(infixAlign2), 2);
		assignSource(row(infixAlign2, 0), infix(seq1, leftPosition(seed, 0), seedBegin1));
		assignSource(row(infixAlign2, 1), infix(seq2, leftPosition(seed, 1), seedBegin2));
		_pump_trace_2_Align(infixAlign2, trace);

		std::cout << "Banded alignment on left extension:" << std::endl;
		std::cout << infixAlign2;
		std::cout << std::endl;

		value(pos, 0) = toViewPosition(row(align, 0), leftPosition(seed, 0));
		value(pos, 1) = toViewPosition(row(align, 1), leftPosition(seed, 1));
		std::cout << leftPosition(seed, 0) << " " << sourceBeginPosition(row(align, 0)) << std::endl;
		std::cout << value(pos, 0) << " " << value(pos, 1) << std::endl;
		integrateAlign(align, infixAlign2, pos);
		std::cout << "align with banded alignment on left extension:" << std::endl;
		std::cout << align;

		std::cout << "Right extension seq1: " << infix(seq1, seedEnd1, rightPosition(seed, 0)+1) << std::endl;
		std::cout << "Right extension seq2: " << infix(seq2, seedEnd2, rightPosition(seed, 1)+1) << std::endl;
		std::cout << std::endl;

		StringSet<Infix<DnaString>::Type> infixStrRight;
		appendValue(infixStrRight, infix(seq1, seedEnd1, rightPosition(seed, 0)+1));
		appendValue(infixStrRight, infix(seq2, seedEnd2, rightPosition(seed, 1)+1));

		clear(trace.sizes); clear(trace.tvs);
		globalAlignment(trace, infixStrRight, scoring, -3, 3, BandedNeedlemanWunsch());

		Align<Infix<DnaString>::Type > infixAlign3;
		resize(rows(infixAlign3), 2);
		assignSource(row(infixAlign3, 0), infix(seq1, seedEnd1, rightPosition(seed, 0)+1));
		assignSource(row(infixAlign3, 1), infix(seq2, seedEnd2, rightPosition(seed, 1)+1));
		_pump_trace_2_Align(infixAlign3, trace);

		std::cout << "Banded alignment on right extension:" << std::endl;
		std::cout << infixAlign3;
		std::cout << std::endl;

		value(pos, 0) = toViewPosition(row(align, 0), leftPosition(seed, 0));
		value(pos, 1) = toViewPosition(row(align, 1), leftPosition(seed, 1));
		integrateAlign(align, infixAlign3, pos);
		std::cout << "align with banded alignment on right extension:" << std::endl;
		std::cout << align;

		clearGaps(align);
		setSourceBeginPosition(row(align, 0), 0);
		setSourceBeginPosition(row(align, 1), 0);
		setBeginPosition(row(align, 0), 0);
		setBeginPosition(row(align, 1), 0);
		setSourceEndPosition(row(align, 0), length(seq1));
		setSourceEndPosition(row(align, 1), length(seq2));
	}
}

void birtesTestViewPos() {
    Align<DnaString> align;
    resize(rows(align), 2);
    assignSource(row(align, 0), "acgtttaccg");
    assignSource(row(align, 1), "agtttatcg");

    globalAlignment(align, Score<int>(1, -1, -1));
    std::cout << align;
	std::cout << "begin1: " << beginPosition(row(align, 0)) << std::endl;
    std::cout << "begin2: " << beginPosition(row(align, 1)) << std::endl;
    std::cout << "sourceBegin1: " << sourceBeginPosition(row(align, 0)) << std::endl;
    std::cout << "sourceBegin2: " << sourceBeginPosition(row(align, 1)) << std::endl;

    setSourceBeginPosition(row(align, 0), toSourcePosition(row(align, 0), 2));
    setSourceBeginPosition(row(align, 1), toSourcePosition(row(align, 1), 2));
    std::cout << align;

	std::cout << "begin1: " << beginPosition(row(align, 0)) << std::endl;
    std::cout << "begin2: " << beginPosition(row(align, 1)) << std::endl;
    std::cout << "sourceBegin1: " << sourceBeginPosition(row(align, 0)) << std::endl;
    std::cout << "sourceBegin2: " << sourceBeginPosition(row(align, 1)) << std::endl;
    std::cout << "View Position of Position 2 in seq1: " << toViewPosition(row(align, 0), 2) << std::endl;
    std::cout << "View Position of Position 2 in seq2: " << toViewPosition(row(align, 1), 2) << std::endl;

	std::cout << "My view position of pos 2 in seq1: " << toViewPosition(row(align, 0), 2 + sourceBeginPosition(row(align, 0))) - toViewPosition(row(align, 0), sourceBeginPosition(row(align, 0))) << std::endl;
	std::cout << "My view position of pos 2 in seq2: " << toViewPosition(row(align, 1), 2 + sourceBeginPosition(row(align, 1))) - toViewPosition(row(align, 1), sourceBeginPosition(row(align, 1))) << std::endl;
}

int main()
{
/////Two DNA sequences that shall be aligned.
//	typedef String<Dna> TSequence;
//	TSequence seq1 = "atcgaatgcgga";
//	TSequence seq2 = "actcgttgca";
/////Scoring objects are used to define a scoring scheme.
/////In this case, affine gap costs with match = 0, mismatch = -1, gapextend = -1 and gapopen = -2.
//	Score<int> score(0, -1, -1, -2);
/////Example 1: We use @Class.Align@ to align the two sequences. 
/////Since we do not specify an @Tag.Global Alignment Algorithms|algorithm tag@ when we call @Function.globalAlignment@, 
/////a suitable algorithm (@Tag.Global Alignment Algorithms|Gotoh@) is automatically choosen.
//	Align<TSequence, ArrayGaps> align;
//	resize(rows(align), 2);
//	assignSource(row(align, 0), seq1);
//	assignSource(row(align, 1), seq2);
//
//	::std::cout << "Score = " << globalAlignment(align, score) << ::std::endl;
//	::std::cout << align << ::std::endl;
/////Example 2: We now choose explicitely the algorithm @Tag.Global Alignment Algorithms|MyersHirschberg@.
/////Since this algorithm always works on Levenshtein distance, $score$ is ignored here.
/////Therefore, this algorithm computes a different alignment and returns a different score.
//	::std::cout << "Score = " << globalAlignment(align, score, MyersHirschberg()) << ::std::endl;
//	::std::cout << align << ::std::endl;
/////Example 3: We now do the same as in case 1, but now we use an @Spec.Alignment Graph@ for storing the alignment.
/////Here we use @Tag.Global Alignment Algorithms|Gotoh's algorithm@.
//	typedef StringSet<TSequence, Dependent<> > TStringSet;
//	typedef Graph<Alignment<TStringSet, void> > TAlignmentGraph;
//
//	TStringSet string_set;
//	appendValue(string_set, seq1);
//	appendValue(string_set, seq2);
//	TAlignmentGraph alignment_graph(string_set);
//
//	::std::cout << "Score = " << globalAlignment(alignment_graph, score, Gotoh()) << ::std::endl;
//	::std::cout << alignment_graph << ::std::endl;

    //birtesTestViewPos();
	birtesTestStichAlign2();

	return 0;
}
