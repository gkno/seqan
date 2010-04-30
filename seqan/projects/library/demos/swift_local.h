#include <iostream>
#include <seqan/index.h>
#include <seqan/seeds.h>

using namespace seqan;

// returns true if align has a match at pos, otherwise false
template<typename TSource, typename TSize>
inline bool
isMatch(Align<TSource> const & align, TSize pos) {
SEQAN_CHECKPOINT
    if(isGap(row(align, 0), pos)) {
        return false;
    } else if(isGap(row(align, 1), pos)) {
        return false;
    } else if(row(align, 0)[pos] != row(align, 1)[pos]) {
        return false;
    } else {
        return true;
    }
}

template<typename TSource, typename TPos>
void
_fillGapsString(Align<TSource> const & align,
                String<Triple<TPos, TPos, TPos> > & gaps) {
SEQAN_CHECKPOINT
    typedef Triple<TPos, TPos, TPos> TGapInfo;
    TPos totalErrors = 0;
    TPos gapBegin = beginPosition(row(align, 0));
    TPos i = gapBegin;

    // append gap starting at beginPosition (also if its length is 0!)
    while(i < endPosition(row(align, 0)) && !isMatch(align, i)) {
        ++i;
        ++totalErrors;
    }
    appendValue(gaps, TGapInfo(gapBegin, i, totalErrors));

    // iterate over alignment and append gaps
    while (i < endPosition(row(align, 0))) {
        // skip matches
        while(i < endPosition(row(align, 0)) && isMatch(align, i)) {
            ++i;
        }
        gapBegin = i;
        // skip and count mismatches/indels
        while(i < endPosition(row(align, 0)) && !isMatch(align, i)) {
            ++i;
            ++totalErrors;
        }
        appendValue(gaps, TGapInfo(gapBegin, i, totalErrors));
    }
    /*for(unsigned l = 0; l < length(gaps); ++l) {
        std::cout << gaps[l].i1 << "  " << gaps[l].i2 << "  " << gaps[l].i3 << std::endl;
    }*/
}

// checks the error rate of the fragment between end of left and start of right
template<typename TPos, typename TFloat>
inline bool
_isEpsMatch(Triple<TPos, TPos, TPos> const & left,
           Triple<TPos, TPos, TPos> const & right,
           TFloat eps) {
SEQAN_CHECKPOINT
    // compute mismatches/indels and length
    TPos errors = right.i3 - left.i3 - (right.i2 - right.i1);
    TPos length = right.i1 - left.i2;

    // check error rate
    return errors/(TFloat)(length) <= eps;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Identifies the longest epsilon match in align and sets the view positions of
// align to start and end position of the longest epsilon match
template<typename TSource, typename TSize, typename TFloat>
void 
longestEpsMatch(Align<TSource> & align,
                    TSize matchMinLength,
                    TFloat epsilon) {
SEQAN_CHECKPOINT
    // Preprocessing: compute and store gaps and lengths
    // A gap is a triple of gap begin position, gap end position, and total number of errors of sequence from begin
    //   to end position of this gap.
    typedef typename Position<Align<TSource> >::Type TPosition;
    typedef String<Triple<TPosition, TPosition, TPosition> > TGapsString;
    TGapsString gaps;
    _fillGapsString(align, gaps);

    // Identify longest eps match by iterating over combinations of left and right positions
    typename Iterator<TGapsString >::Type rightIt = end(gaps) - 1;
    typename Iterator<TGapsString >::Type leftIt = begin(gaps);

    TPosition beginPos = 0;
    TPosition endPos = 0;
    TSize minLength = matchMinLength - 1;
    
    while ((*leftIt).i2 + minLength < (*rightIt).i1) {
        while ((*leftIt).i2 + minLength < (*rightIt).i1) {
            if(_isEpsMatch(*leftIt, *rightIt, epsilon)) {
                beginPos = (*leftIt).i2;
                endPos = (*rightIt).i1;
                minLength = endPos - beginPos;
                break;
            }
            --rightIt;
        }
        rightIt = end(gaps) - 1;
        ++leftIt;
    }

    // Set view positions to the eps-match
	setSourceBeginPosition(row(align, 0), toSourcePosition(row(align, 0), beginPos));
	setSourceBeginPosition(row(align, 1), toSourcePosition(row(align, 1), beginPos));
	setBeginPosition(row(align, 0), beginPos);
	setBeginPosition(row(align, 1), beginPos);
	setSourceEndPosition(row(align, 0), toSourcePosition(row(align, 0), endPos));
	setSourceEndPosition(row(align, 1), toSourcePosition(row(align, 1), endPos));
}

//template<typename TInfix>
//void
//insertMatch(String<Align<TInfix> > & matches, Align<TInfix> & match) {
//    typedef typename Size<String<Align<TInfix> > >::Type TSize;
//    TSize len = length(matches);
//    if(len == 0) {
//        appendValue(matches, match);
//    } else {
//        typename Iterator<String<Align<TInfix> > >::Type iter = begin(matches);
//        // TODO
//    }
//}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// banded chain alignment and X-drop extension for all local alignments with a min score
template<typename TInfixA, typename TInfixB, typename TEpsilon, typename TSize, typename TDelta, typename TDrop>
int 
verifySwiftHit(TInfixA const & a,
               TInfixB const & b,
               TEpsilon eps,
               TSize minLength,
               TDelta delta,
               TDrop xDrop,
               String<Align<TInfixB> > & matches) {
    typedef Seed<int, SimpleSeed> TSeed;
    typedef typename Position<TInfixB>::Type TPos;

    TPos maxLength = 1000000000;
    if (length(a)*length(b) > maxLength) {
        std::cerr << "Warning: SWIFT hit <" << beginPosition(a) << "," << endPosition(a);
        std::cerr << "> , <" << beginPosition(b) << "," << endPosition(b) << "> too long... verification skipped.\n" << std::flush;
        return 0;
    }

    // define a scoring scheme
    typedef int TScore;
    TScore match = 1;
    TScore mismatchIndel = (int)(-1/eps) + 1;
    Score<TScore> scoreMatrix(match, mismatchIndel, mismatchIndel);
    TScore scoreDropOff = xDrop * (-mismatchIndel);
    
    // create alignment object for the infixes
    Align<TInfixB> alignment;
    resize(rows(alignment), 2);
    assignSource(row(alignment, 0), a);
    assignSource(row(alignment, 1), b);

    // calculate minimal score for local alignments
    TEpsilon e = floor(eps*minLength);
    TSize minLength1 = (TSize)ceil((e+1) / eps);
    TEpsilon e1 = floor(eps*minLength1);
    TSize minScore = _min((TSize)ceil((minLength-e) / (e+1)), (TSize)ceil((minLength1-e1) / (e1+1)));

    // local alignment
    __int64 upperDiag = 0;
    __int64 lowerDiag = endPosition(a) - (__int64)endPosition(b) - beginPosition(a) + beginPosition(b);
    if (beginPosition(b) == 0) upperDiag = lowerDiag + delta;
    if (endPosition(b) == endPosition(host(b))) lowerDiag = upperDiag - delta;
    LocalAlignmentFinder<> finder = LocalAlignmentFinder<>();
    int count = 0;
    while (localAlignment(alignment, finder, scoreMatrix, minScore, lowerDiag, upperDiag, BandedWatermanEggert())) {
        //splitXDrop

        // gapped X-drop extension
        ++count;
        TSeed seed(sourceBeginPosition(row(alignment, 0)) + beginPosition(a),
                   sourceBeginPosition(row(alignment, 1)) + beginPosition(b),
                   sourceEndPosition(row(alignment, 0)) + beginPosition(a)-1,
                   sourceEndPosition(row(alignment, 1)) + beginPosition(b)-1);
        //TSize errors = (score - length(row(alignment, 0)) * match) / (mismatchIndel - match);
        //TSize scoreDropOff = 2*score;//(int)(errors - 2 * eps * length(row(alignment, 0))) * mismatchIndel + score;//
        extendSeed(seed, scoreDropOff, scoreMatrix, host(a), host(b), 2, GappedXDrop());

        // banded alignment
        Align<TInfixB> extAlign;
        resize(rows(extAlign), 2);
        assignSource(row(extAlign, 0), infix(host(a), leftDim0(seed), rightDim0(seed)+1));
        assignSource(row(extAlign, 1), infix(host(b), leftDim1(seed), rightDim1(seed)+1));

        TSize band = xDrop;
        bandedAlignment(extAlign, seed, band, scoreMatrix);
        //StringSet<TInfixB> str;
        //appendValue(str, infix(host(a), leftDim0(seed), rightDim0(seed)+1));
        //appendValue(str, infix(host(b), leftDim1(seed), rightDim1(seed)+1));
        //_Align_Traceback<TSize> trace;
        //globalAlignment(trace, str, scoreMatrix, rightDiagonal(seed), leftDiagonal(seed), BandedNeedlemanWunsch()); // funktioniert nicht auf Infixen!!!
        //_pump_trace_2_Align(extAlign, trace);

        if ((TSize)length(row(extAlign, 0)) < minLength) continue;

        // cut ends to obtain longest contained epsilon-match
        longestEpsMatch(extAlign, minLength, eps);

        if ((TSize)length(row(extAlign, 0)) < minLength) continue;

        // insert e-match in matches string
        //insertMatch(matches, extAlign);
        typename Iterator<String<Align<TInfixB> > >::Type iter = begin(matches);
        while (iter != end(matches)) {
            if (beginPosition(sourceSegment(row(extAlign,0))) == beginPosition(sourceSegment(row(*iter,0)))
                && endPosition(sourceSegment(row(extAlign,0))) == endPosition(sourceSegment(row(*iter,0)))
                && beginPosition(sourceSegment(row(extAlign,1))) == beginPosition(sourceSegment(row(*iter,1)))
                && endPosition(sourceSegment(row(extAlign,1))) == endPosition(sourceSegment(row(*iter,1)))) {
                    break;
            }
            ++iter;
        }
        if (iter == end(matches)) {
            //std::cout << alignment << std::endl;
            //std::cout << score << std::endl;
            //std::cout << extAlign << std::endl;
            appendValue(matches, extAlign);
        }
    }
    return count;
}

// calls swift, verifies swift hits, outputs eps-matches
template<typename TText, typename TIndex, typename TSize, typename TDrop, typename TInfix>
int localSwift(Finder<TText, Swift<SwiftLocal> > & finder,
                Pattern<TIndex, Swift<SwiftLocal> > & pattern,
                double epsilon,
                TSize minLength,
                TDrop xDrop,
                StringSet<String<Align<TInfix> > > & matches) {
    resize(matches, countSequences(needle(pattern)));
    TSize numSwiftHits = 0;

    int count = 0;
	while (find(finder, pattern, epsilon, minLength)) {
        ++numSwiftHits;
        //std::cout << positionRange(finder) << " ; " << positionRange(pattern) << std::endl;
        // verification
        count += verifySwiftHit(range(finder), range(pattern), epsilon, minLength,
                              pattern.bucketParams[0].delta + pattern.bucketParams[0].overlap,
                              xDrop, value(matches, pattern.curSeqNo));
	}
    std::cout << count << std::endl;
    return numSwiftHits;
}
