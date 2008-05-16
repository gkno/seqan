 /*==========================================================================
                SeqAn - The Library for Sequence Analysis
                          http://www.seqan.de 
 ============================================================================
  Copyright (C) 2007

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  Lesser General Public License for more details.

 ============================================================================
  $Id$
 ==========================================================================*/

#ifndef SEQAN_HEADER_FIND_SWIFT_H
#define SEQAN_HEADER_FIND_SWIFT_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// SWIFT to search a text for
// - semi-global alignments of one/multiple short sequences
// - local epsilon matches of one/multiple short sequences
//////////////////////////////////////////////////////////////////////////////

template < typename TObject, typename TSpec > 
class Index;

struct _SwiftLocal;
typedef Tag<_SwiftLocal> SwiftLocal;

struct _SwiftSemiGlobal;
typedef Tag<_SwiftSemiGlobal> SwiftSemiGlobal;


template <typename TSpec = SwiftSemiGlobal>
struct Swift {
	enum { SEMIGLOBAL = 1 };		// 0..match eps-match of min.length n0; 1..match the whole read
	enum { DIAGONAL = 1 };			// 0..use rectangular buckets (QUASAR); 1..use diagonal buckets (SWIFT)
	enum { LOG2DELTA_MIN = 3 };		// set minimal delta to 8
	enum { THRESHOLD_MIN = 1 };	    // set minimal threshold to 1
#ifdef WITH_1HULL
	enum { QGRAM_ERRORS = 1 };		// allow 1 error per q-gram
#else
	enum { QGRAM_ERRORS = 0 };		// allow 0 error per q-gram
#endif
	enum { HAMMING_ONLY = 1 };
	enum { PARAMS_BY_LENGTH = 1 };	// params are determined only by seq.length
};

template <>
struct Swift<SwiftLocal> {
	enum { SEMIGLOBAL = 0 };	// 0..match eps-match of min.length n0; 1..match the whole read
	enum { DIAGONAL = 1 };		// 0..use rectangular buckets (QUASAR); 1..use diagonal buckets (SWIFT)
	enum { LOG2DELTA_MIN = 3 };	// set minimal delta to 8
	enum { THRESHOLD_MIN = 1 };	// set minimal threshold to 1
	enum { QGRAM_ERRORS = 0 };	// allow 0 errors per q-gram
	enum { HAMMING_ONLY = 0 };
	enum { PARAMS_BY_LENGTH = 0 };
};

//////////////////////////////////////////////////////////////////////////////

	template <typename TSpec, typename TSize, typename TShortSize = unsigned short>
	struct _SwiftBucket {
		TSize			firstIncrement;
		TSize			lastIncrement;
		TShortSize		counter;
#ifdef SEQAN_DEBUG_SWIFT
		TSize			_lastIncDiag;
#endif
	};

	template <typename TSize, typename TShortSize>
	struct _SwiftBucket<SwiftSemiGlobal, TSize, TShortSize> {
		TSize			lastIncrement;
		TShortSize		counter;
#ifdef SEQAN_DEBUG_SWIFT
		TSize			_lastIncDiag;
#endif
	};

	template <typename TSpec, typename TSize, typename TShortSize = unsigned short>
	struct _SwiftBucketParams {
		TSize			firstBucket;	// first _SwiftBucket entry in pattern.buckets
		TSize			reuseMask;		// 2^ceil(log2(x)) reuse every x-th bucket)
		TShortSize		distanceCut;	// if lastIncrement is this far or farer away, threshold can't be reached
		TShortSize		delta;			// buckets begin at multiples of delta
		TShortSize		threshold;		// at least threshold q-gram hits induce an approx match
		TShortSize		overlap;		// number of diagonals/columns a bucket shares with its neighbor
		unsigned char	logDelta;		// log2(delta)
	};

//____________________________________________________________________________


	template <typename THstkPos>
	struct _SwiftHit {
		THstkPos	hstkPos;			// parallelogram begin in haystack 
		unsigned	ndlSeqNo;			// needle sequence number
		unsigned	bucketWidth;		// (non-diagonal) bucket width (bktHeight + delta + overlap (for diagonals))
	};

//____________________________________________________________________________


	template < typename THaystack, typename TSpec >
	class Finder< THaystack, Swift<TSpec> >
	{
	public:
		typedef typename Iterator<THaystack, Rooted>::Type		TIterator;
		typedef typename Position<THaystack>::Type				THstkPos;
		typedef _SwiftHit<THstkPos>								TSwiftHit;
		typedef String<TSwiftHit>								THitString;
		typedef typename Iterator<THitString, Standard>::Type	THitIterator;

		TIterator		data_iterator;
		TIterator		haystackEnd;
		bool			_needReinit;	// if true, the Pattern needs to be reinitialized
		bool			_atEnd;
		THitString		hits;
		THitIterator	curHit, endHit;
		THstkPos		curPos;

		Finder():
			_needReinit(true) {}

		Finder(THaystack &haystack):
			data_iterator(begin(haystack, Rooted())),
			_needReinit(true) {}

		Finder(TIterator &iter):
			data_iterator(iter),
			_needReinit(true) {}

		Finder(TIterator const &iter):
			data_iterator(iter),
			_needReinit(true) {}

		Finder(Finder const &orig):
			data_iterator(orig.data_iterator),
			haystackEnd(orig.haystackEnd),
			hits(orig.hits),
			_needReinit(orig._needReinit) 
		{
			curHit = begin(hits, Standard()) + (orig.curHit - begin(orig.hits, Standard()));
			endHit = end(hits, Standard());
		};

		inline typename Reference<TIterator>::Type 
		operator* () { return value(hostIterator(*this)); }

		inline typename Reference<TIterator const>::Type 
		operator* () const { return value(hostIterator(*this)); }

		operator TIterator () const	{ return data_iterator;	}
	};


//____________________________________________________________________________

	// forward
    template < typename TInput, typename TSpec >
    struct Pipe;

	template < typename TTuples, typename TPipeSpec, typename TSpec >
	class Finder< Pipe<TTuples, TPipeSpec>, Swift<TSpec> >
	{
	public:
		typedef Pipe<TTuples, TPipeSpec>						TInput;
		typedef typename Size<TInput>::Type						THstkPos;
		typedef _SwiftHit<THstkPos>								TSwiftHit;
		typedef String<TSwiftHit>								THitString;
		typedef typename Iterator<THitString, Standard>::Type	THitIterator;

		TInput			&in;
		bool			_needReinit;	// if true, the Pattern needs to be reinitialized
		THitString		hits;
		THitIterator	curHit, endHit;
		THstkPos		curPos, dotPos, dotPos2;

		Finder(TInput &_in):
			in(_in),
			_needReinit(true) {}

		Finder(Finder const &orig):
			in(orig.in),
			hits(orig.hits),
			_needReinit(orig._needReinit) 
		{
			curHit = begin(hits, Standard()) + (orig.curHit - begin(orig.hits, Standard()));
			endHit = end(hits, Standard());
		};
	};


//____________________________________________________________________________

	
	template <typename THaystack, typename TSpec>
	inline bool
	atEnd(Finder<THaystack, Swift<TSpec> > & me)
	{
		return hostIterator(hostIterator(me)) == hostIterator(me.haystackEnd);
	}

	template <typename THaystack, typename TSpec>
	inline void
	goEnd(Finder<THaystack, Swift<TSpec> > & me)
	{
		hostIterator(me) = me.haystackEnd;
	}


//____________________________________________________________________________


	template < typename TNeedle, typename TIndexSpec, typename TSpec >
	class Pattern< Index<TNeedle, TIndexSpec>, Swift<TSpec> >
	{
	public:
		typedef Index<TNeedle, TIndexSpec>								TIndex;
		typedef typename Size<TIndex>::Type								TSize;
		typedef unsigned short											TShortSize;
		typedef typename Fibre<TIndex, Tag<_Fibre_SA> const >::Type		TSA;
		typedef typename Fibre<TIndex, Tag<_Fibre_Shape> const >::Type	TShape;
		typedef typename Iterator<TSA const, Standard>::Type			TIterator;
		
		typedef _SwiftBucket<TSpec, TSize, TShortSize>					TBucket;
		typedef String<TBucket>											TBucketString;
		typedef _SwiftBucketParams<TSpec, TSize, TShortSize>			TBucketParams;
		typedef String<TBucketParams>									TBucketParamsString;
		
		TShape					shape;
		TBucketString			buckets;
		TBucketParamsString		params;
		unsigned				curSeqNo;

		Holder<TIndex>	data_host;

		Pattern() 
		{
			clear(*this);
		}
		Pattern(TIndex &_index): data_host(_index) 
		{
			clear(*this);
		}
		Pattern(TIndex const &_index): data_host(_index)
		{
			clear(*this);
		}
	};
	
//____________________________________________________________________________


template <typename TParams>
inline void _printSwiftParams(TParams &params)
{
	::std::cout << "  firstBucket: " << params.firstBucket << ::std::endl;
	::std::cout << "  reuseMask:   " << params.reuseMask << ::std::endl;
	::std::cout << "  distanceCut: " << params.distanceCut << ::std::endl;
	::std::cout << "  delta:       " << params.delta << ::std::endl;
	::std::cout << "  threshold:   " << params.threshold << ::std::endl;
	::std::cout << "  overlap:     " << params.overlap << ::std::endl;
	::std::cout << "  logDelta:    " << (int)params.logDelta << ::std::endl << ::std::endl;
}

template < typename TNeedle, typename TIndexSpec, typename TSpec >
inline void _printSwiftBuckets(Pattern< Index<TNeedle, TIndexSpec>, Swift<TSpec> > &p)
{
	typedef Index<TNeedle, TIndexSpec> TIndex;
	typedef typename Pattern<TIndex, Swift<TSpec> >::TBucketParams TParams;

	unsigned j = 0;
	TParams *params = &_getSwiftParams(p, 0);

	for(unsigned i=0; i<length(p.buckets) && i<10; ++i) 
	{
		if ((i & params->reuseMask) == 0)
		{
			::std::cout << ::std::endl << "ReadBucket #" << j << "    " << '"';
			::std::cout << indexText(host(p))[j] << '"' << ::std::endl;
			::std::cout << "  length:      " << sequenceLength(j, host(p)) << ::std::endl;
			params = &_getSwiftParams(p, j++);
			_printSwiftParams(*params);
		}

		::std::cout << "    lastInc: " << (int)p.buckets[i].lastIncrement;
		::std::cout << "  \tCounter: " << p.buckets[i].counter << ::std::endl;
	}
}

template <typename TIndex, typename TSpec, typename TSize>
inline typename Pattern<TIndex, Swift<TSpec> >::TBucketParams &
_getSwiftParams(Pattern<TIndex, Swift<TSpec> > & pattern, TSize seqNo) 
{
	if (Swift<TSpec>::PARAMS_BY_LENGTH)
		return pattern.params[sequenceLength(seqNo, host(pattern))];
	else
		return pattern.params[seqNo];
}

template <typename TIndex, typename TSpec, typename TParams, typename TSize>
inline unsigned
_getSwiftBucketNo(Pattern<TIndex, Swift<TSpec> > const &, TParams & params, TSize seqNo) 
{
	if (Swift<TSpec>::PARAMS_BY_LENGTH)
		return (params.reuseMask + 1) * seqNo;
	else
		return params.firstBucket;
}

template <typename TIndex, typename TFloat, typename _TSize, typename TSpec>
inline void _patternInit(Pattern<TIndex, Swift<TSpec> > & pattern, TFloat errorRate, _TSize minLengthForAll) 
{
	typedef Pattern<TIndex, Swift<TSpec> >			TPattern;
	typedef typename Size<TIndex>::Type				TSize;
	typedef typename Fibre<TIndex, QGram_SA>::Type	TSA;
	typedef typename Iterator<TSA, Standard>::Type	TSAIter;
	typedef typename TPattern::TBucket				TBucket;
	typedef typename TPattern::TBucketParams		TBucketParams;
	
	indexRequire(host(pattern), QGram_SADir());
	TIndex const &index = host(pattern);
	TSize span = length(pattern.shape);
	TSize count = 0;
	TSize seqCount = countSequences(index);
	TSize bucketsPerCol2Max = 0;
	TSize maxLength = 0;
	
	if (Swift<TSpec>::PARAMS_BY_LENGTH) {
		for(unsigned seqNo = 0; seqNo < seqCount; ++seqNo) {
			TSize length = sequenceLength(seqNo, host(pattern));
			if (maxLength < length)
				maxLength = length;
		}
		resize(pattern.params, maxLength + 1);
	} else
		resize(pattern.params, seqCount);
	
	if (minLengthForAll != 0) 
	{
		// global matches
		TSize minLength = minLengthForAll;
		for(unsigned seqNo = 0; seqNo < seqCount; ++seqNo) 
		{
			// swift q-gram lemma
			TBucketParams &params = _getSwiftParams(pattern, seqNo);
			// n..next length that could decrease threshold
			TSize n = (TSize) ceil((floor(errorRate * minLength) + 1) / errorRate);
			// minimal threshold is minimum errors of minLength and n
			int threshold = (TSize) min(
				(n + 1) - span * (floor(errorRate * n) + 1),
				(minLength + 1) - span * (floor(errorRate * minLength) + 1));
			if (threshold < Swift<TSpec>::THRESHOLD_MIN) threshold = Swift<TSpec>::THRESHOLD_MIN;
			params.threshold = threshold;

			TSize errors = (TSize) floor((2 * params.threshold + span - 3) / (1 / errorRate - span));
			params.overlap = errors;
			params.distanceCut = (params.threshold - 1) + span * (errors + span);
			params.logDelta = (TSize) ceil(log((double)errors) / log(2.0));
			params.delta = 1 << params.logDelta;
			// TODO: classical swift for rectangular buckets
		}
	} else
		for(unsigned seqNo = 0; seqNo < seqCount; ++seqNo) 
		{
			// get pattern length and max. allowed errors
			TBucketParams &params = _getSwiftParams(pattern, seqNo);
			TSize length = sequenceLength(seqNo, host(pattern));
			TSize errors = (TSize) floor(errorRate * length);
			TSize errorsWC = errors / (1 + Swift<TSpec>::QGRAM_ERRORS);

			// q-gram lemma: How many conserved q-grams we see at least?
			// (define a minimal threshold of 1)
			int threshold = length + 1 - span * (errorsWC + 1);
			if (threshold < Swift<TSpec>::THRESHOLD_MIN) threshold = Swift<TSpec>::THRESHOLD_MIN;
			params.threshold = threshold;
			
			if (Swift<TSpec>::HAMMING_ONLY == 0)
				errors = 0;			

			// a bucket has distanceCut different positions of q-grams
			// if a q-gram is this far or farer away it can't belong to the
			// same bucket
			params.distanceCut = length - (span - 1) + errors;

			TSize bucketsPerCol2;
			if (Swift<TSpec>::DIAGONAL == 1)
			{
				// Use overlapping parallelograms				
				params.overlap = errors;
				
				// delta must be a power of 2 greater then errors (define a minimal delta of 8)
				params.logDelta = (TSize) ceil(log((double)(errors + 1)) / log(2.0));
				if (params.logDelta < Swift<TSpec>::LOG2DELTA_MIN) params.logDelta = Swift<TSpec>::LOG2DELTA_MIN;
				params.delta = 1 << params.logDelta;

				// the formula for bucketsPerCol is (worst-case):
				// (height-(q-1) - 1 - (delta+1-e))/delta + 3
				//    ^-- full paral. in the middle --^     ^-- 2 at the bottom, 1 at the top
				TSize bucketsPerCol = (length - span + 2 * params.delta + errors - 1) / params.delta;
				bucketsPerCol2 = 1 << (TSize) ceil(log((double)bucketsPerCol) / log(2.0));
			}
			else
			{
				// Use overlapping rectangles
				params.overlap = length - span + errors;

				// delta must be a power of 2 greater then seq.length + errors (define a minimal delta of 32)
				params.logDelta = (TSize) ceil(log((double)(length - span + 1 + errors)) / log(2.0));
				if (params.logDelta < Swift<TSpec>::LOG2DELTA_MIN) params.logDelta = Swift<TSpec>::LOG2DELTA_MIN;
				params.delta = 1 << params.logDelta;

				bucketsPerCol2 = 2;
			}

			SEQAN_ASSERT(distanceCut <= bucketsPerCol * (TSize) delta);

			params.firstBucket = count;
			params.reuseMask = bucketsPerCol2 - 1;
			
			if (Swift<TSpec>::PARAMS_BY_LENGTH) {
				++count;
				if (bucketsPerCol2Max < bucketsPerCol2)
					bucketsPerCol2Max = bucketsPerCol2;
			} else
				count += bucketsPerCol2;
			
/*			if (seqNo<3)
				_printSwiftParams(params);
*/		}

	if (Swift<TSpec>::PARAMS_BY_LENGTH) {
		count *= bucketsPerCol2Max;
		for(unsigned i = 0; i < length(pattern.params); ++i)
			pattern.params[i].reuseMask = bucketsPerCol2Max - 1;
	}
			
	resize(pattern.buckets, count);
	for(unsigned seqNo = 0, i = 0, next = 0; seqNo < seqCount; ++seqNo)
	{
		TBucketParams &params = _getSwiftParams(pattern, seqNo);
		next = i + params.reuseMask + 1;
		for(; i < next; ++i) {
			pattern.buckets[i].lastIncrement = 0 - params.distanceCut;
			pattern.buckets[i].counter = 0;
		}
	}
}

template <
	typename THaystack,
	typename TIndex, 
	typename TSpec,
	typename THValue
>
inline bool _swiftMultiProcessQGram(
	Finder<THaystack, Swift<TSpec> > &finder, 
	Pattern<TIndex, Swift<TSpec> > &pattern,
	THValue hash)
{
	typedef Finder<THaystack, Swift<TSpec> >					TFinder;
	typedef Pattern<TIndex, Swift<TSpec> >						TPattern;

	typedef typename Size<TIndex>::Type							TSize;
	typedef typename Fibre<TIndex, QGram_SA>::Type				TSA;
	typedef typename Iterator<TSA, Standard>::Type				TSAIter;
	typedef typename TPattern::TBucketString					TBucketString;
	typedef typename Iterator<TBucketString, Standard>::Type	TBucketIter;
	typedef typename TPattern::TBucketParams					TBucketParams;
	typedef typename TFinder::TSwiftHit							THit;
	
	TIndex const &index = host(pattern);
	
	TSAIter saBegin = begin(indexSA(index), Standard());
	TSAIter occ =  saBegin + indexDir(index)[hash];
	TSAIter occEnd = saBegin + indexDir(index)[hash + 1];
	TBucketIter bktBegin = begin(pattern.buckets, Standard());
	Pair<unsigned> ndlPos;

	for(; occ != occEnd; ++occ) 
	{
		posLocalize(ndlPos, *occ, stringSetLimits(index));
		TBucketParams &params = _getSwiftParams(pattern, getSeqNo(ndlPos));

		__int64 diag = finder.curPos;
		if (Swift<TSpec>::DIAGONAL == 1) diag -= getSeqOffset(ndlPos);
		
		unsigned bktNo = (diag >> params.logDelta) & params.reuseMask;
		unsigned bktOfs = diag & (params.delta - 1);

		TBucketIter bkt = bktBegin + (_getSwiftBucketNo(pattern, params, getSeqNo(ndlPos)) + bktNo);
		
		do {
			if ((*bkt).lastIncrement + params.distanceCut <= finder.curPos)
			{
				// too far away - could be:
				// 1. too far to reach the treshold
				// 2. too far to be in the same bucket (we must ensure that bucketIdx doesn't collide) 
				if ((*bkt).counter >= params.threshold)
				{
					// upper bucket no. of lastIncr. q-gram
					unsigned upperBktNo = (*bkt).lastIncrement >> params.logDelta;

					TSize height = 0;
					if (Swift<TSpec>::DIAGONAL == 1)
						height = sequenceLength(getSeqNo(ndlPos), host(pattern)) - 1;

					// we must decrement bucket no. until (no. mod reuse == bktNo)
					__int64 bktBeginHstk = 
						(__int64) (upperBktNo - ((upperBktNo - bktNo) & params.reuseMask)) << params.logDelta;

					if (bktBeginHstk >= 0) {
						THit hit = {
							bktBeginHstk,							// bucket begin in haystack
							getSeqNo(ndlPos),						// needle seq. number
							height + params.delta + params.overlap	// bucket width (non-diagonal)
						};
						appendValue(finder.hits, hit);
					} else {
						// match begins left of haystack begin
						THit hit = {
							0,										// bucket begin in haystack
							getSeqNo(ndlPos),						// needle seq. number
							height + params.delta + params.overlap	// bucket width (non-diagonal)
							+ bktBeginHstk
						};
						appendValue(finder.hits, hit);
					}
				}
				(*bkt).lastIncrement = finder.curPos;
				(*bkt).counter = 1;
#ifdef SEQAN_DEBUG_SWIFT
				(*bkt)._lastIncDiag = diag;
#endif
			}
			else
			{
				if ((*bkt).lastIncrement == finder.curPos) 
					break;	// increment only once per sequence			
				(*bkt).lastIncrement = finder.curPos;
				if (++(*bkt).counter <= 0)
					(*bkt).counter = params.threshold;
#ifdef SEQAN_DEBUG_SWIFT
				(*bkt)._lastIncDiag = diag;
#endif
			}

			if (bktOfs >= params.overlap) break;
			
			// repeat with the previous overlapping bucket
			bktOfs = params.overlap;
			if (!bktNo) {
				bktNo = params.reuseMask;
				bkt += bktNo;
			} else {
				--bktNo;
				--bkt;
			}
		} while (true);
	}

	finder.curHit = begin(finder.hits, Standard());
	finder.endHit = end(finder.hits, Standard());

	return !empty(finder.hits);
}

template <
	typename THaystack,
	typename TIndex, 
	typename TSpec
>
inline bool _swiftMultiFlushBuckets(
	Finder<THaystack, Swift<TSpec> > &finder, 
	Pattern<TIndex, Swift<TSpec> > &pattern)
{
	typedef Finder<THaystack, Swift<TSpec> >					TFinder;
	typedef Pattern<TIndex, Swift<TSpec> >						TPattern;

	typedef typename TPattern::TBucketString					TBucketString;
	typedef typename Iterator<TBucketString, Standard>::Type	TBucketIterator;
	typedef typename TPattern::TBucketParams					TBucketParams;

	typedef typename Size<TIndex>::Type							TSize;
	typedef typename TFinder::TSwiftHit							THit;
	
	TBucketIterator	bkt = begin(pattern.buckets, Standard());
	TBucketIterator	bktEnd;
	TSize seqCount = countSequences(host(pattern));
	
	for(TSize ndlSeqNo = 0; ndlSeqNo < seqCount; ++ndlSeqNo) 
	{
		TBucketParams &params = _getSwiftParams(pattern, ndlSeqNo);
		bktEnd = bkt + (params.reuseMask + 1);
		for(unsigned bktNo = 0; bkt != bktEnd; ++bkt, ++bktNo)
		{
			if ((*bkt).counter >= params.threshold)
			{
				// upper bucket no. of lastIncr. q-gram
				unsigned upperBktNo = (*bkt).lastIncrement >> params.logDelta;

				TSize height = 0;
				if (Swift<TSpec>::DIAGONAL == 1)
					height = sequenceLength(ndlSeqNo, host(pattern)) - 1;

				// we must decrement bucket no. until (no. mod reuse == bktNo)
				__int64 bktBeginHstk = 
					(__int64) (upperBktNo - ((upperBktNo - bktNo) & params.reuseMask)) << params.logDelta;

				if (bktBeginHstk >= 0) {
					THit hit = {
						bktBeginHstk,							// bucket begin in haystack
						ndlSeqNo,								// needle seq. number
						height + params.delta + params.overlap	// bucket width (non-diagonal)
					};
					appendValue(finder.hits, hit);
				} else {
					// match begins left of haystack begin
					THit hit = {
						0,										// bucket begin in haystack
						ndlSeqNo,								// needle seq. number
						height + params.delta + params.overlap	// bucket width (non-diagonal)
						+ bktBeginHstk
					};
					appendValue(finder.hits, hit);
				}

			}
			(*bkt).lastIncrement = 0 - params.distanceCut;
			(*bkt).counter = 0;
		}
	}

	finder.curHit = begin(finder.hits, Standard());
	finder.endHit = end(finder.hits, Standard());

	return !empty(finder.hits);
}

template <typename TNeedle, typename TIndexSpec, typename TSpec>
inline bool 
empty(Pattern<Index<TNeedle, TIndexSpec>, Swift<TSpec> > & me) 
{
	return empty(me.params);
}

template <typename TNeedle, typename TIndexSpec, typename TSpec>
inline void 
clear(Pattern<Index<TNeedle, TIndexSpec>, Swift<TSpec> > & me) 
{
	clear(me.params);
	clear(me.buckets);
}

//____________________________________________________________________________

template <typename THaystack, typename TSpec>
inline typename Position<Finder<THaystack, Swift<TSpec> > >::Type
position(Finder<THaystack, Swift<TSpec> > & finder)
{
	typename Finder<THaystack, Swift<TSpec> >::TSwiftHit &hit = *finder.curHit;
	return hit.hstkPos + hit.bucketWidth;
}

template <typename THaystack, typename TSpec>
inline typename Position<Finder<THaystack, Swift<TSpec> > >::Type
position(Finder<THaystack, Swift<TSpec> > const & finder)
{
	typename Finder<THaystack, Swift<TSpec> >::TSwiftHit &hit = *finder.curHit;
	return hit.hstkPos + hit.bucketWidth;
}

//____________________________________________________________________________

template <typename THaystack, typename TSpec>
inline typename Position<Finder<THaystack, Swift<TSpec> > >::Type
beginPosition(Finder<THaystack, Swift<TSpec> > & finder)
{
	typename Finder<THaystack, Swift<TSpec> >::TSwiftHit &hit = *finder.curHit;
	return hit.hstkPos;
}

template <typename THaystack, typename TSpec>
inline typename Position<Finder<THaystack, Swift<TSpec> > >::Type
beginPosition(Finder<THaystack, Swift<TSpec> > const & finder)
{
	typename Finder<THaystack, Swift<TSpec> >::TSwiftHit &hit = *finder.curHit;
	return hit.hstkPos;
}

//____________________________________________________________________________

template <typename THaystack, typename TSpec>
inline typename Position<Finder<THaystack, Swift<TSpec> > >::Type
endPosition(Finder<THaystack, Swift<TSpec> > & finder)
{
	typename Finder<THaystack, Swift<TSpec> >::TSwiftHit &hit = *finder.curHit;
	return hit.hstkPos + hit.bucketWidth;
}

template <typename THaystack, typename TSpec>
inline typename Position<Finder<THaystack, Swift<TSpec> > >::Type
endPosition(Finder<THaystack, Swift<TSpec> > const & finder)
{
	typename Finder<THaystack, Swift<TSpec> >::TSwiftHit &hit = *finder.curHit;
	return hit.hstkPos + hit.bucketWidth;
}

//____________________________________________________________________________

template <typename THaystack, typename TSpec>
inline typename Infix<THaystack>::Type
range(Finder<THaystack, Swift<TSpec> > &finder)
{
	typedef typename Size<THaystack>::Type TSize;

	typename Finder<THaystack, Swift<TSpec> >::TSwiftHit &hit = *finder.curHit;
	TSize hitEnd = hit.hstkPos + hit.bucketWidth;
	TSize textEnd = length(haystack(finder));

	if (hitEnd > textEnd)
		return infix(haystack(finder), hit.hstkPos, textEnd);
	else
		return infix(haystack(finder), hit.hstkPos, hitEnd);
}

template <typename THaystack, typename TSpec, typename TText>
inline typename Infix<TText>::Type
range(Finder<THaystack, Swift<TSpec> > &finder, TText &text)
{
	typedef  typename Size<TText>::Type TSize;

	typename Finder<THaystack, Swift<TSpec> >::TSwiftHit &hit = *finder.curHit;
	TSize hitEnd = hit.hstkPos + hit.bucketWidth;
	TSize textEnd = length(text);

	if (hitEnd > textEnd)
		return infix(text, hit.hstkPos, textEnd);
	else
		return infix(text, hit.hstkPos, hitEnd);
}

template <typename TNeedle, typename TIndexSpec, typename TSpec>
inline typename Value<TNeedle>::Type &
range(Pattern<Index<TNeedle, TIndexSpec>, Swift<TSpec> > &pattern)
{
	return indexText(needle(pattern))[pattern.curSeqNo];
}


template <typename THaystack, typename TNeedle, typename TIndexSpec, typename TSpec>
inline bool 
find(
	 Finder<THaystack, Swift<TSpec> > &finder,
	 Pattern<Index<TNeedle, TIndexSpec>, Swift<TSpec> > &pattern, 
	 double errorRate)
{
	if (empty(finder)) 
	{
		_patternInit(pattern, errorRate, 0);
		_finderSetNonEmpty(finder);
		finder.haystackEnd =
			hostIterator(finder) + (length(host(finder)) - length(pattern.shape) + 1);

		if (atEnd(finder)) return false;
		finder.curPos = 0;
		if (_swiftMultiProcessQGram(finder, pattern, hash(pattern.shape, hostIterator(hostIterator(finder)))))
		{
			pattern.curSeqNo = (*finder.curHit).ndlSeqNo;
			return true;
		}
	} else
		if (++finder.curHit != finder.endHit) 
		{
			pattern.curSeqNo = (*finder.curHit).ndlSeqNo;
			return true;
		}

	clear(finder.hits);
	do {
		++finder;
		if (atEnd(finder)) 
			if (_swiftMultiFlushBuckets(finder, pattern))
				break;
			else
				return false;
		++finder.curPos;
	} while (!_swiftMultiProcessQGram(finder, pattern, hashNext(pattern.shape, hostIterator(hostIterator(finder)))));

	pattern.curSeqNo = (*finder.curHit).ndlSeqNo;
	return true;
}

template <typename THashes, typename TPipeSpec, typename TNeedle, typename TIndexSpec, typename TSpec>
inline bool 
find(
	 Finder<Pipe<THashes, TPipeSpec>, Swift<TSpec> > &finder,
	 Pattern<Index<TNeedle, TIndexSpec>, Swift<TSpec> > &pattern, 
	 double errorRate,
	 bool printDots)
{
	if (empty(finder)) 
	{
		_patternInit(pattern, errorRate, 0);
		_finderSetNonEmpty(finder);
		finder.dotPos = 100000;
		finder.dotPos2 = 10 * finder.dotPos;

		beginRead(finder.in);
		if (eof(finder.in)) 
		{
			endRead(finder.in);
			return false;
		}
		finder.curPos = (*finder.in).i1;
		if (_swiftMultiProcessQGram(finder, pattern, hash(pattern.shape, (*finder.in).i2)))
		{
			pattern.curSeqNo = (*finder.curHit).ndlSeqNo;
			return true;
		}
	} else
		if (++finder.curHit != finder.endHit) 
		{
			pattern.curSeqNo = (*finder.curHit).ndlSeqNo;
			return true;
		}

	clear(finder.hits);
	if (eof(finder.in)) return false;

	do {
		++finder.in;
		if (eof(finder.in)) {
			endRead(finder.in);
			_printSwiftBuckets(pattern);
			if (_swiftMultiFlushBuckets(finder, pattern))
				break;
			else
				return false;
		}
		finder.curPos = (*finder.in).i1;
		if (printDots)
			if (finder.curPos == finder.dotPos) 
			{
				finder.dotPos += 100000;
				if (finder.curPos == finder.dotPos2)
				{
					finder.dotPos2 += 1000000;
					::std::cerr << (finder.curPos / 1000000) << "MBp" << ::std::flush;
				} else
					::std::cerr << "." << ::std::flush;
			}
	} while (!_swiftMultiProcessQGram(finder, pattern, hash(pattern.shape, (*finder.in).i2)));

	pattern.curSeqNo = (*finder.curHit).ndlSeqNo;
	return true;
}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_FIND_SHIFTAND_H
