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

//// returns true if align has a mismatch at pos, otherwise false
//template<typename TSource, typename TSize>
//inline bool
//isMismatch(Align<TSource> const & align, TSize pos) {
//SEQAN__CHECKPOINT
//    if(isGap(row(align, 0), pos)) {
//        return false;
//    } else if(isGap(row(align, 1), pos)) {
//        return false;
//    } else if(row(align, 0)[pos] != row(align, 1)[pos]) {
//        return true;
//    } else {
//        return false;
//    }
//}

///////////////////////////////////////////////////////////////////////////

template<typename TAlign, typename TPos, typename TScoreValue>
inline void
_appendNegativeSegment(TAlign const & align, 
                       TPos & pos, TPos len, 
                       Score<TScoreValue> const & scoreMatrix, 
                       String<Triple<TPos, TPos, TScoreValue> > & queue) {
SEQAN_CHECKPOINT
    typedef Triple<TPos, TPos, TScoreValue> TMerger;
    TPos beginPos = pos;

    TScoreValue score = 0;
    while (pos < len) {
        if (isGap(row(align, 0), pos) || isGap(row(align, 1), pos)) {
            score += scoreGap(scoreMatrix);
        } else if (value(row(align, 0), pos) != value(row(align, 1), pos)) {
            score += scoreMismatch(scoreMatrix);
        } else {
            break;
        }
        ++pos;
    }
    if (pos == len) appendValue(queue, TMerger(beginPos, pos, infimumValue<TScoreValue>()));
    else appendValue(queue, TMerger(beginPos, pos, score));
}

template<typename TAlign, typename TPos, typename TScoreValue>
inline void
_appendPositiveSegment(TAlign const & align, 
                       TPos & pos, TPos len, 
                       Score<TScoreValue> const & scoreMatrix, 
                       String<Triple<TPos, TPos, TScoreValue> > & queue) {
SEQAN_CHECKPOINT
    if (pos == len) return;
    typedef Triple<TPos, TPos, TScoreValue> TMerger;
    TPos beginPos = pos;

    TScoreValue score = 0;
    while ((pos < len) && 
           (!isGap(row(align, 0), pos) &&
            !isGap(row(align, 1), pos) &&
            (value(row(align, 0), pos) == value(row(align, 1), pos)))) {
        score += scoreMatch(scoreMatrix);
        ++pos;
    }
    appendValue(queue, TMerger(beginPos, pos, score));
}

template<typename TMerger>
inline bool
_negativeMerge(String<TMerger> & queue) {
SEQAN_CHECKPOINT
    typedef typename TMerger::T1 TPos;
    TPos len = length(queue);
    if (len < 3) return false;

    TMerger cd = value(queue, len-1);
    TMerger bc = value(queue, len-2);
    TMerger ab = value(queue, len-3);

    if ((bc.i3 < 0) || (bc.i3 >= abs(_max(ab.i3, cd.i3)))) {
        return false;
    } else {
        String<TMerger> newMerger;
        appendValue(newMerger, TMerger(ab.i1, cd.i2, ab.i3 + bc.i3 + cd.i3)); 
        replace(queue, len-3, len, newMerger);

        return true;
    }
}

template<typename TMerger>
inline bool
_positiveMerge(String<TMerger> & queue) {
SEQAN_CHECKPOINT
    typedef typename TMerger::T1 TPos;
    TPos len = length(queue);
    if (len < 5) return false;

    TMerger ef = value(queue, len-1);
    TMerger de = value(queue, len-2);
    TMerger cd = value(queue, len-3);
    TMerger bc = value(queue, len-4);
    TMerger ab = value(queue, len-5);

    if ((cd.i3 >= 0) || (cd.i3 < _max(ab.i3, ef.i3))) {
        return false;
    } else {
        String<TMerger> newMerger;
        appendValue(newMerger, TMerger(bc.i1, de.i2, bc.i3 + cd.i3 + de.i3)); 
        replace(queue, len-4, len-1, newMerger);

        return true;
    }
}

///////////////////////////////////////////////////////////////////////////

// Implements the algorithm from Zhang et al. in Bioinformatics, 1999: "Post-processing long pairwise alignments".
// Splits an alignment into sub-alignments that contain no x-Drop.
template<typename TAlign, typename TScoreValue, typename TScoreValue1, typename TScoreValue2>
void
_splitAtXDrops(TAlign const & align,
              Score<TScoreValue> & scoreMatrix,
              TScoreValue1 scoreDropOff,
              TScoreValue2 minScore,
              String<TAlign> & alignmentString) {
SEQAN_CHECKPOINT
    typedef typename Position<Row<TAlign> >::Type TPos;
    typedef Triple<TPos, TPos, TScoreValue> TMerger;
    
    // initialization
    String<TMerger> queue;
    TPos pos = _min(toViewPosition(row(align, 0), sourceBeginPosition(row(align, 0))),
                    toViewPosition(row(align, 1), sourceBeginPosition(row(align, 1))));
    appendValue(queue, TMerger(pos, pos, infimumValue<TScoreValue1>()));

    TPos aliLength = _max(toViewPosition(row(align, 0), sourceEndPosition(row(align, 0))),
                          toViewPosition(row(align, 1), sourceEndPosition(row(align, 1))));
    TPos len;
    while ((pos < aliLength) || (length(queue) > 1)) {
        // construct useful tree
        if (!_negativeMerge(queue)) {
            if (!_positiveMerge(queue)) {
                _appendPositiveSegment(align, pos, aliLength, scoreMatrix, queue);
                _appendNegativeSegment(align, pos, aliLength, scoreMatrix, queue);
            }
        }

        // check for x-Drop
        len = length(queue);
        if ((len == 3) && (value(queue, 2).i3 < scoreDropOff * (-1))) {
            if (value(queue, 1).i3 >= minScore) {
                // append new sub-alignment
                TPos begin = value(queue, 1).i1;
                TPos end = value(queue, 1).i2;

                TAlign ali(align);
                setSourceBeginPosition(row(ali, 0), toSourcePosition(row(ali, 0), begin));
                setSourceBeginPosition(row(ali, 1), toSourcePosition(row(ali, 1), begin));
                setSourceEndPosition(row(ali, 0), toSourcePosition(row(ali, 0), end));
                setSourceEndPosition(row(ali, 1), toSourcePosition(row(ali, 1), end));
                appendValue(alignmentString, ali);
            }
            replace(queue, 0, 2, String<TMerger>());
        }
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

template<typename TRow>
inline bool
_isUpstream(TRow & row1, TRow & row2) {
SEQAN_CHECKPOINT
    typedef typename Position<TRow>::Type TPos;
    TPos end1 = endPosition(sourceSegment(row1));
    TPos begin2 = beginPosition(sourceSegment(row2));
    if (end1 <= begin2) return true;
    else return false;
}

//template<TRow>
//inline bool
//_isDownstream(TRow & row1, TRow & row2) {
//QAN_CHECKPOINT
//    _isUpstream(row2, row1);
//}

template<typename TSource>
inline void
_insertMatch(String<Align<TSource> > & matches, Align<TSource> & match) {
SEQAN_CHECKPOINT
    typedef typename Position<String<Align<TSource> > >::Type TPosString;
    typedef typename Row<Align<TSource> >::Type TRow;

    TPosString maxPos = length(matches);
    TPosString minPos = 0;

    if(maxPos == 0) {
        // matches string is empty
        appendValue(matches, match);
        return;
    }

    TRow matchRow0 = row(match, 0);
    TRow matchRow1 = row(match, 1);

    // determine insertion position
    TPosString pos = maxPos / 2;
    while (pos < maxPos) {
        TRow row0 = row(value(matches, pos), 0);
        if (_isUpstream(row0, matchRow0)) {
            // match is downstream of matches[pos] in first row
            if (minPos == pos) {
                ++pos;
                break;
            }
            minPos = pos;
        } else if (_isUpstream(matchRow0, row0)) {
            // match is upstream of matches[pos] in first row
            maxPos = pos;
        } else {
            // match overlaps in first row with another match
            TRow row1 = row(value(matches, pos), 1);
            if (_isUpstream(row1, matchRow1)) {
                // match is downstream of matches[pos] in second row
                if (minPos == pos) {
                    ++pos;
                    break;
                }
                minPos = pos;
            } else if (_isUpstream(matchRow1, row1)) {
                // match is upstream of matches[pos] in second row
                maxPos = pos;
            } else {
                // match overlaps in both rows with another match -> keep longer match
                // TODO: If non-overlapping part is greater than minLength, keep the match!!!
                if (length(matchRow0) > length(row0)) {
                    // new match is longer
                    replace(matches, pos, pos+1, String<Align<TSource> >());
                    --maxPos;
                } else {
                    return;
                }
            }
        }
        pos = minPos + ((maxPos-minPos) / 2);
    }
    // insert match in matches at position pos
    String<Align<TSource> > str;
    appendValue(str, match);
    replace(matches, pos, pos, str);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// banded chain alignment and X-drop extension for all local alignments with a min score
template<typename TInfixA, typename TInfixB, typename TEpsilon, typename TSize, typename TDelta, typename TDrop>
void // int // 
verifySwiftHit(TInfixA const & a,
               TInfixB const & b,
               TEpsilon eps,
               TSize minLength,
               TDelta delta,
               TDrop xDrop,
               String<Align<TInfixB> > & matches) {
SEQAN_CHECKPOINT
    typedef Seed<int, SimpleSeed> TSeed;
    typedef typename Position<TInfixB>::Type TPos;

    TPos maxLength = 1000000000;
    if (length(a) > maxLength) {
        std::cerr << "Warning: SWIFT hit <" << beginPosition(a) << "," << endPosition(a);
        std::cerr << "> , <" << beginPosition(b) << "," << endPosition(b) << "> too long... verification skipped.\n" << std::flush;
        return;
    }
    //int count = length(matches);

    TSize minLengthWithoutErrors = minLength - (int)floor(minLength*eps);

    // define a scoring scheme
    typedef int TScore;
    TScore match = 1;
    TScore mismatchIndel = (int)(-1/eps) + 1;
    Score<TScore> scoreMatrix(match, mismatchIndel, mismatchIndel);
    TScore scoreDropOff = xDrop * (-mismatchIndel);
    
    // create alignment object for the infixes
    Align<TInfixB> localAlign;
    resize(rows(localAlign), 2);
    assignSource(row(localAlign, 0), a);
    assignSource(row(localAlign, 1), b);

    // calculate minimal score for local alignments
    TEpsilon e = floor(eps*minLength);
    TSize minLength1 = (TSize)ceil((e+1) / eps);
    TEpsilon e1 = floor(eps*minLength1);
    TSize minScore = _min((TSize)ceil((minLength-e) / (e+1)), (TSize)ceil((minLength1-e1) / (e1+1)));

    // banded local alignment
    __int64 upperDiag = 0;
    __int64 lowerDiag = endPosition(a) - (__int64)endPosition(b) - beginPosition(a) + beginPosition(b);
    if (beginPosition(b) == 0) upperDiag = lowerDiag + delta;
    if (endPosition(b) == endPosition(host(b))) lowerDiag = upperDiag - delta;
    LocalAlignmentFinder<> finder = LocalAlignmentFinder<>();
    while (localAlignment(localAlign, finder, scoreMatrix, minScore, lowerDiag, upperDiag, BandedWatermanEggert())) {

        // split local alignments containing an X-drop
        String<Align<TInfixB> > alignmentString;
        _splitAtXDrops(localAlign, scoreMatrix, scoreDropOff, minScore, alignmentString);

        typename Iterator<String<Align<TInfixB> > >::Type aliIt = begin(alignmentString);
        while (aliIt != end(alignmentString)) {

            // determine extension direction
            char direction;
            if (length(alignmentString) == 1) direction = 2;
            else if (aliIt == begin(alignmentString)) direction = 0;
            else if (aliIt == end(alignmentString)-1) direction = 1;
            else direction = 3;

            Align<TInfixB> extAlign;
            if (direction == 3) {
                extAlign = Align<TInfixB>(*aliIt);
            }
            else {
                // gapped X-drop extension
                TSeed seed(sourceBeginPosition(row(*aliIt, 0)) + beginPosition(a),
                           sourceBeginPosition(row(*aliIt, 1)) + beginPosition(b),
                           sourceEndPosition(row(*aliIt, 0)) + beginPosition(a) - 1,
                           sourceEndPosition(row(*aliIt, 1)) + beginPosition(b) - 1);
                extendSeed(seed, scoreDropOff, scoreMatrix, host(a), host(b), direction, GappedXDrop());

                if (length(seed) < minLengthWithoutErrors) {
                    ++aliIt;
                    continue;
                }

                // banded alignment
                StringSet<TInfixB> str;
                appendValue(str, infix(host(a), leftPosition(seed, 0), rightPosition(seed, 0)+1));
                appendValue(str, infix(host(b), leftPosition(seed, 1), rightPosition(seed, 1)+1));
                __int64 startDiag = leftPosition(seed, 1) - leftPosition(seed, 0);
                _Align_Traceback<TSize> trace;
                globalAlignment(trace, str, scoreMatrix, startDiag - leftDiagonal(seed), startDiag - rightDiagonal(seed), BandedNeedlemanWunsch());

                resize(rows(extAlign), 2);
                assignSource(row(extAlign, 0), infix(host(a), leftPosition(seed, 0), rightPosition(seed, 0)+1));
                assignSource(row(extAlign, 1), infix(host(b), leftPosition(seed, 1), rightPosition(seed, 1)+1));
                _pump_trace_2_Align(extAlign, trace);
            }

            if ((TSize)length(row(extAlign, 0)) < minLength) {
                ++aliIt;
                continue;
            }

            // cut ends to obtain longest contained epsilon-match
            longestEpsMatch(extAlign, minLength, eps);

            if ((TSize)length(row(extAlign, 0)) < minLength) {
                ++aliIt;
                continue;
            }

            // insert e-match in matches string
            _insertMatch(matches, extAlign);
            ++aliIt;
        }
    }
    //return length(matches) - count;
}

// calls swift, verifies swift hits, outputs eps-matches
template<typename TText, typename TIndex, typename TSize, typename TDrop, typename TInfix>
int localSwift(Finder<TText, Swift<SwiftLocal> > & finder,
                Pattern<TIndex, Swift<SwiftLocal> > & pattern,
                double epsilon,
                TSize minLength,
                TDrop xDrop,
                StringSet<String<Align<TInfix> > > & matches) {
SEQAN_CHECKPOINT
    resize(matches, countSequences(needle(pattern)));
    TSize numSwiftHits = 0;

    TSize maxLength = 0;
    TSize totalLength = 0;

    //String<unsigned> counts;
    //for (unsigned i = 0; i < 6; ++i) {
    //    appendValue(counts, 0);
    //}
    //int c = 0;
    //int maxCount = 0;
	while (find(finder, pattern, epsilon, minLength)) {
        ++numSwiftHits;

        // verification
        /*c = */verifySwiftHit(infix(finder), infix(pattern), epsilon, minLength,
                              pattern.bucketParams[0].delta + pattern.bucketParams[0].overlap,
                              xDrop, value(matches, pattern.curSeqNo));
        totalLength += length(infix(finder));
        if ((TSize)length(infix(finder)) > maxLength) maxLength = length(infix(finder));

        //if (c < 6 && c >= 0) ++counts[c];
        //if (c < 0) std::cout << c << std::endl;
        //if (c > maxCount) maxCount = c;
	}
    //for (unsigned i = 0; i < 6; ++i) {
    //    std::cout << counts[i] << "  ";
    //}
    //std::cout << "\nmax: " << maxCount << std::endl;

    std::cout << "Longest hit: " << maxLength << std::endl;
    if (numSwiftHits > 0) std::cout << "Avg hit length: " << totalLength/numSwiftHits << std::endl;
    
    return numSwiftHits;
}
