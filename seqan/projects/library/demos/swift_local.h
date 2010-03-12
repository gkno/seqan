#include <iostream>
#include <seqan/index.h>
#include <seqan/seeds.h>

using namespace seqan;


template<typename TSource, typename TPosition, typename TSize, typename TFloat>
bool
isEpsMatch(Align<TSource> & align, TPosition left, TPosition right, TSize minLength, TFloat eps) {
    // count mismatches
    TSize mismatches = 0;
    for (TPosition i = left; i < right; ++i) {
        if (row(align, 0)[i] != row(align, 1)[i]) {
            ++mismatches;
        }
    }
    
    // check error rate
    if (mismatches/(TFloat)(right-left) <= eps)
        return true;
    return false;
}


template<typename TSource, typename TSize, typename TFloat>
void 
shrinkToMaxEpsMatch(Align<TSource> & align, TSize matchMinLength, TFloat epsilon) {
    typedef typename Position<Align<TSource> >::Type TPosition;
    TPosition begin = 0;
    TPosition end = 0;
    TSize minLength = matchMinLength - 1;

    // iterate over all right positions
    TPosition right = endPosition(row(align,0));
    TPosition left;
    while (right > beginPosition(row(align, 0)) + minLength) {
        if (row(align, 0)[right-1] != row(align, 1)[right-1]) {
            --right;
        } else {
            // iterate over all left positions until eps match found 
            // or fragment is shorter than minLength
            left = beginPosition(row(align, 0));
            while (left < right - minLength) {
                if (row(align, 0)[left] != row(align, 1)[left]) {
                    ++left;
                } else {
                    // check for new longest epsilon match
                    if (isEpsMatch(align, left, right, matchMinLength, epsilon)) {
                        if ((TSize)(right - left) > minLength) {
                            begin = left;
                            end = right;
                            minLength = right - left;
                        }
                        break;
                    }
                    // skip matches at left end
                    while (left < right - minLength && row(align, 0)[left] == row(align, 1)[left])
                        ++left;
                }
            }
            // skip matches at right end
            while (right > beginPosition(row(align, 0)) + minLength && row(align, 0)[right-1] == row(align, 1)[right-1])
                --right;
        }
    }
    
    // set view positions to the alignment
	setSourceBeginPosition(row(align, 0), toSourcePosition(row(align, 0), begin));
	setSourceBeginPosition(row(align, 1), toSourcePosition(row(align, 1), begin));
	setBeginPosition(row(align, 0), 0);//begin); // TODO: ask David what this function is needed for!
	setBeginPosition(row(align, 1), 0);//begin);
	setSourceEndPosition(row(align, 0), toSourcePosition(row(align, 0), end));
	setSourceEndPosition(row(align, 1), toSourcePosition(row(align, 1), end));
    //std::cout << "begin:" << begin << " end:" << end << std::endl;
    //std::cout << align << std::endl;
}


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
        TSize scoreDropOff = (int)(2 * length(row(alignment, 0)) * eps)*match - errors*mismatchIndel;
        extendSeed(seed, scoreDropOff, scoreMatrix, host(a), host(b), 2, GappedXDrop());

        // banded alignment
        Align<TInfixB> extAlign;
        resize(rows(extAlign), 2);
        assignSource(row(extAlign, 0), infix(host(a), leftDim0(seed), rightDim0(seed)));
        assignSource(row(extAlign, 1), infix(host(b), leftDim1(seed), rightDim1(seed)));
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
        if (iter == end(matches) && (TSize)length(row(extAlign,0)) >= minLength) {
            appendValue(matches, extAlign);
        }
    }
}


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
        Iterator<String<TSeed> >::Type it = begin(chain);
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