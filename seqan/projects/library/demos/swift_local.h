#include <iostream>
#include <seqan/index.h>
#include <seqan/seeds.h>

using namespace seqan;

// returns true if align has a match at pos, otherwise false
template<typename TSource, typename TSize>
inline bool
isMatch(Align<TSource> const & align, TSize pos) {
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
    typedef Triple<TPos, TPos, TPos> TGapInfo;
    TPos totalErrors = 0;
    TPos gapBegin = beginPosition(row(align, 0));
    TPos i = gapBegin;

    // append starting gap also if its length is 0
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
    // compute mismatches and length
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
shrinkToMaxEpsMatch(Align<TSource> & align,
                    TSize matchMinLength,
                    TFloat epsilon) {
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

template<typename TSource, typename TPosition, typename TFloat>
inline bool
_isEpsMatch(Align<TSource> & align, TPosition left, TPosition right, TFloat eps) {
    // count mismatches
    typename Size<Align<TSource> >::Type errors = 0;
    for (TPosition i = left; i < right; ++i) {
        if (!isMatch(align, i)) {
            ++errors;
        }
    }
    
    // check error rate
    return errors/(TFloat)(right-left) <= eps;
}

template<typename TSource, typename TSize, typename TFloat>
void 
shrinkToMaxEpsMatchSimple(Align<TSource> & align,
                    TSize matchMinLength,
                    TFloat epsilon) {
    typedef typename Position<Align<TSource> >::Type TPosition;
    TPosition begin = 0;
    TPosition end = 0;
    TSize minLength = matchMinLength - 1;

    // iterate over all right positions
    TPosition right = endPosition(row(align,0));
    TPosition left;
    while (right > beginPosition(row(align, 0)) + minLength) {
        if (!isMatch(align, right-1)) {
            --right;
        } else {
            // iterate over all left positions until eps match found 
            // or fragment is shorter than minLength
            left = beginPosition(row(align, 0));
            while (right > left + minLength) {
                if (!isMatch(align, left)) {
                    ++left;
                } else {
                    // check for new longest epsilon match
                    if (_isEpsMatch(align, left, right, epsilon)) {
                        begin = left;
                        end = right;
                        minLength = right - left;
                        break;
                    }
                    // skip matches at left end
                    while (right > left + minLength && isMatch(align, left))
                        ++left;
                }
            }
            // skip matches at right end
            while (right > beginPosition(row(align, 0)) + minLength && isMatch(align, right-1))
                --right;
        }
    }
    
    // set view positions to the eps-match
	setSourceBeginPosition(row(align, 0), toSourcePosition(row(align, 0), begin));
	setSourceBeginPosition(row(align, 1), toSourcePosition(row(align, 1), begin));
	setBeginPosition(row(align, 0), begin);
	setBeginPosition(row(align, 1), begin);
	setSourceEndPosition(row(align, 0), toSourcePosition(row(align, 0), end));
	setSourceEndPosition(row(align, 1), toSourcePosition(row(align, 1), end));
    //std::cout << "begin:" << begin << " end:" << end << std::endl;
    //std::cout << align << std::endl;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// banded chain alignment and X-drop extension for all local alignments with a min score
template<typename TInfixA, typename TInfixB, typename TEpsilon, typename TSize, typename TMatches>
void 
verifySwiftHitByLocalAlign(TInfixA const & a,
               TInfixB const & b,
               TEpsilon eps,
               TSize minLength,
               TMatches & matches) {
    typedef Seed<int, SimpleSeed> TSeed;
        
    // define a scoring scheme
    typedef int TScore;
    TScore match = 1;
    TScore mismatchIndel = (int)(-1/eps) + 1;
    Score<TScore> scoreMatrix(match, mismatchIndel, mismatchIndel);
    
    // create alignment object for the infixes
    Align<TInfixB> alignment;
    resize(rows(alignment), 2);
    assignSource(row(alignment, 0), a);
    assignSource(row(alignment, 1), b);

    // calculate minimal score for local alignments
    TSize minLength1 = (TSize)ceil((floor(eps*minLength)+1) / eps);
    TSize minScore = (TSize)ceil((minLength - floor(eps*minLength)) / (floor(eps*minLength)+1));
    minScore = _min(minScore, (TSize)ceil((minLength1 - floor(eps*minLength1)) / (floor(eps*minLength1)+1)));

    // local alignment
    LocalAlignmentFinder<> finder(alignment);
    while (localAlignment(alignment, finder, scoreMatrix, minScore-1)) {
        // gapped X-drop extension
        TSeed seed(sourceBeginPosition(row(alignment, 0)) + beginPosition(a),
                   sourceBeginPosition(row(alignment, 1)) + beginPosition(b),
                   sourceEndPosition(row(alignment, 0)) + beginPosition(a),
                   sourceEndPosition(row(alignment, 1)) + beginPosition(b));
        TSize errors = (getScore(finder) - length(row(alignment, 0)) * match) / (mismatchIndel - match);
        TSize scoreDropOff = -(int)(2 * length(row(alignment, 0)) * eps * mismatchIndel) + errors * mismatchIndel;
        extendSeed(seed, scoreDropOff, scoreMatrix, host(a), host(b), 2, GappedXDrop());

        // banded alignment
        Align<TInfixB> extAlign;
        resize(rows(extAlign), 2);
        assignSource(row(extAlign, 0), infix(host(a), leftDim0(seed), rightDim0(seed)+1));
        assignSource(row(extAlign, 1), infix(host(b), leftDim1(seed), rightDim1(seed)+1));
        bandedAlignment(extAlign, seed, errors, scoreMatrix);

        if ((TSize)length(row(extAlign, 0)) < minLength) continue;

        // cut ends to obtain longest contained epsilon-match
        shrinkToMaxEpsMatch(extAlign, minLength, eps);

        if ((TSize)length(row(extAlign, 0)) < minLength) continue;

        // insert e-match in matches string
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
        if (iter == end(matches) &&
            (TSize)length(row(extAlign,0)) >= minLength) {
            appendValue(matches, extAlign);
        }
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// create all chains of threshold q grams that where the starting positions 
// of no two q-grams are further than limit apart
template<typename TSeed, typename TInfixA, typename TInfixB, typename TThreshold, typename TLimit>
void
qGramChains(String< String<TSeed> > & outputChains,
                 TInfixA const & a,
                 TInfixB const & b,
                 TThreshold q,
                 TThreshold threshold,
                 TLimit limit) {
    typedef String<TSeed> TChain;
    typedef String<TChain> TChainList;
    TChainList chainList;

    typedef Index<DnaString, Index_QGram<SimpleShape> > TQGramIndex;
	TQGramIndex index_qgram(b); 
    resize(indexShape(index_qgram), q);

	typedef Finder<TQGramIndex> TFinder;
    for (unsigned i = beginPosition(a); i < endPosition(a)-q+1; i++) {
        TFinder finder(index_qgram);
		
        while (find(finder, infix(host(a), i, i+q))) {
		    // for all common q-grams

		    // create a new seed from the q-gram
		    typedef typename Position<TFinder>::Type TFinderPos;
		    TFinderPos posA = i;
		    TFinderPos posB = beginPosition(b) + position(finder);
		    TSeed newSeed(posA, posB, q);

		    unsigned numChains = length(chainList);
            unsigned i = 0;     // current position in the chain list
            for (unsigned k = 0; k < numChains; k++) {
                // for all chains in the current chainList
                TSeed lastSeed = *(end(chainList[i]) - 1);
                if(posB <= leftDim1(lastSeed)) {
                    // do nothing if seed is before last seed in chain
                    i++;
                    continue;
                }

                TFinderPos rightPosA = rightDim0(lastSeed) + 1;
                TFinderPos rightPosB = rightDim1(lastSeed) + 1;
                TFinderPos leftPosA = leftDim0(lastSeed);
                TFinderPos leftPosB = leftDim1(lastSeed);
                // ::std::cout << "last seed: " << leftPosA << "," << leftPosB;
                // ::std::cout << " (" << rightPosA << "," << rightPosB << ")" << ::std::endl;
                
                if ((__int64)posA-(__int64)leftDim0(*(begin(chainList[i]))) > limit) {
                    // remove chain if new seed is further than limit apart from first seed
                    erase(chainList, i);
                } else {
                    // append q-gram to copy of the chain if possible
                    if((posA >= rightPosA && posB >= rightPosB) || 
                        (endDiagonal(lastSeed) == endDiagonal(newSeed) && posA > leftPosA && posB > leftPosB)) {
                        TChain chainCopy(chainList[i]);
                        appendValue(chainCopy, newSeed);

                        // append copy of chain to output if its length reached threshold else add it to chainList
                        if (length(chainCopy) == threshold) {
                            appendValue(outputChains, chainCopy);
                        } else {
                            appendValue(chainList, chainCopy);
                        }
                    }
                    i++;
                }
            }

            // add a new chain that starts with this q-gram
            TChain newChain;
            appendValue(newChain, newSeed);
            appendValue(chainList, newChain);
        }
    }
}


// banded chain alignment and gapped X-drop extension for all chains of threshold q-grams
template<typename TInfixA, typename TInfixB, typename TFloat, typename TSize, typename TShortSize, typename TMatches>
void 
verifySwiftHitByChains(TInfixA const & a,
               TInfixB const & b,
               TFloat eps,
               TSize minLength,
               TShortSize q,
               TShortSize threshold,
               TShortSize distanceCut,
               TSize bandwidth,
               TMatches & matches) {
	typedef Seed<int, SimpleSeed> TSeed;
	typedef typename Value<TSeed>::Type TPosition;
    typedef String<TSeed> TChain;
    typedef String<TChain> TChainList;

    TChainList outputChains;
    TShortSize limit = distanceCut - threshold - q + 1;

    // create all chains of threshold q grams with no two q-grams further than limit apart
    qGramChains(outputChains, a, b, q, threshold, limit);

    typename Iterator<TChainList>::Type oIt = begin(outputChains);
    while (oIt != end(outputChains)) {
        // merge overlapping q-grams on same diagonal to one seed
        typename Iterator<TChain>::Type chainIt1 = begin(*oIt);
        typename Iterator<TChain>::Type chainIt2 = begin(*oIt) + 1;
        unsigned i = 0;    // counts the position of chainIt1
        while (chainIt2 != end(*oIt)) {
            if (endDiagonal(*chainIt1) == endDiagonal(*chainIt2)
                && leftDim0(*chainIt2) > leftDim0(*chainIt1) && leftDim0(*chainIt2) <= rightDim0(*chainIt1)+1) {
                // merge the two seeds, set chain iterators to next positions
                setLeftDim0(*chainIt2, leftDim0(*chainIt1));
                setLeftDim1(*chainIt2, leftDim1(*chainIt1));
                erase(*oIt, i);
            } else {
                chainIt1++;
                chainIt2++;
                i++;
            }
        }
        
        // define a scoring scheme
        typedef int TScore;
        TScore match = 1;
        TScore mismatchIndel = -1;
        Score<TScore> scoreMatrix(match, mismatchIndel, mismatchIndel);

        TShortSize totalSeedLength = 0;
        for (chainIt1 = begin(*oIt); chainIt1 != end(*oIt); ++chainIt1) {
            totalSeedLength += length(*chainIt1);
        }

        // gapped X-drop extension
        TSeed seed(leftDim0(*begin(*oIt)),
                   leftDim1(*begin(*oIt)),
                   rightDim0(*(end(*oIt)-1)),
                   rightDim1(*(end(*oIt)-1)));
        TSize errors = (length(seed) - totalSeedLength) / 2;
        TScore scoreDropOff = (TScore)(2* length(seed) * eps) - errors;/*
        std::cout << host(a) << "  " << leftDim0(seed) << "  " << rightDim0(seed) << std::endl;
        std::cout << host(b) << "  " <<  leftDim1(seed) << "  " << rightDim1(seed) << std::endl;
        std::cout << scoreDropOff << "  " << errors << std::endl;*/
        extendSeed(seed, scoreDropOff, scoreMatrix, 
                   host(a), host(b), 2, GappedXDrop());/*
        std::cout << host(a) << "  " << leftDim0(seed) << "  " << rightDim0(seed) << std::endl;
        std::cout << host(b) << "  " <<  leftDim1(seed) << "  " << rightDim1(seed) << std::endl;*/

        // create alignment object for the infixes
        Align<TInfixB> alignment;
        resize(rows(alignment), 2);
        assignSource(row(alignment, 0), infix(host(a), leftDim0(seed), rightDim0(seed)+1));
        assignSource(row(alignment, 1), infix(host(b), leftDim1(seed), rightDim1(seed)+1));

        // banded chain alignment
        TScore score = bandedChainAlignment(*oIt, bandwidth, alignment, scoreMatrix);

        if (score >= length(sourceSegment(row(alignment, 1))) - 2*bandwidth) {
            // cut ends to obtain longest contained epsilon-match
            shrinkToMaxEpsMatch(alignment, minLength, eps);

            // insert new e-match in matches string
            typename Iterator<String<Align<TInfixB> > >::Type iter = begin(matches);
            while (iter != end(matches)) {
                if (beginPosition(sourceSegment(row(alignment,0))) == beginPosition(sourceSegment(row(*iter,0)))
                    && endPosition(sourceSegment(row(alignment,0))) == endPosition(sourceSegment(row(*iter,0)))
                    && beginPosition(sourceSegment(row(alignment,1))) == beginPosition(sourceSegment(row(*iter,1)))
                    && endPosition(sourceSegment(row(alignment,1))) == endPosition(sourceSegment(row(*iter,1)))) {
                        break;
                }
                ++iter;
            }
            if (iter == end(matches)) {
                appendValue(matches, alignment);
            }
        }

        oIt++;
    }
}

// banded chain alignment and gapped X-drop extension for all maximal chains of q-grams
template<typename TInfixA, typename TInfixB, typename TFloat, typename TSize, typename TShortSize, typename TMatches>
void 
verifySwiftHitByMaxChain(TInfixA const & a,
               TInfixB const & b,
               TFloat eps,
               TSize minLength,
               TShortSize q,
               TShortSize threshold,
               TShortSize distanceCut,
               TSize bandwidth,
               TMatches & matches) {
    typedef Seed<int, SimpleSeed> TSeed;
    SeedSet<int, SimpleSeed, DefaultScore> seeds(1, q);
    String<TSeed> chain;//, seeds;

    typedef Index<DnaString, Index_QGram<SimpleShape> > TQGramIndex;
	TQGramIndex index_qgram(b); 
    resize(indexShape(index_qgram), q);

    // create seed set of all common q-grams
	typedef Finder<TQGramIndex> TFinder;
    for (unsigned i = beginPosition(a); i < endPosition(a)-q+1; i++) {
        TFinder finder(index_qgram);
        while (find(finder, infix(host(a), i, i+q))) {
            //TSeed newSeed(i, beginPosition(b) + position(finder), q);
            //appendValue(seeds, newSeed);
            unsigned begin = beginPosition(b) + position(finder);
            std::cout << i << "  " << i+q << "    " << begin << "  " << begin+q << std::endl;
            //if(!addSeed(seeds, i, begin, i+q, begin+q, 0, Merge())) {
                addSeed(seeds, i, begin, q, q, Single());
            //}
        }
    }
    
    // compute chain of q-grams from seed set
    int score = globalChaining(seeds, chain, -1, length(a), length(b));

    if (score > -10/*TODO: score???*/) {
        typename Iterator<String<TSeed> >::Type it = begin(chain);
        while (it != end(chain)) {
            std::cout << rightDim0(*it) << "  " << leftDim0(*it) << "    " << rightDim1(*it) << "  " << leftDim1(*it) << std::endl;
            ++it;
        }
        // TODO:
        // seed extension
        // bandedAlignment
        // cutEnds
        // insert into matches
    }
}