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
  Author: David Weese <weese@fu-berlin.de>
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

// TODO(bkehr): Is this documentatin right? Should the Specializations be Tags?
// weese: of minimal length?
/**
.Spec.Swift:
..summary:Provides a fast filter alogrithm that guarantees to find all regions overlapping with potential \epsilon-matches.
An \epsilon-match is a matching region of minimal length and an error rate of at most \epsilon.
..general:Class.Pattern
..general:Class.Finder
..cat:Searching
..signature:Finder<TIndex, Swift<TSpec> >
..signature:Pattern<TIndex, Swift<TSpec> >
..param.TIndex: A q-gram index of needle(s).
...type:Spec.Index_QGram
..param.TSpec: Specifies the type of Swift filter.
*/
///.Class.Pattern.param.TSpec.type:Spec.Swift

/**
.Spec.SwiftLocal:
..summary:The specialization for the general swift filter that finds epsilon matches between haystack and needle.
..general:Spec.Swift
..cat:Searching
..signature:Finder<TIndex, Swift<SwiftLocal> >
..signature:Pattern<TIndex, Swift<SwiftLocal> >
..param.TIndex: A q-gram index of needle(s).
...type:Spec.Index_QGram
*/
///.Spec.Swift.param.TSpec.type:Spec.SwiftLocal
/**
.Spec.SwiftSemiGlobal:
..summary:The specialization for the semi-global swift filter that finds regions of the haystack where a needle matches with an error rate less than \epsilon.
..general:Spec.Swift
..cat:Searching
..signature:Finder<TIndex, Swift<SwiftSemiGlobal> >
..signature:Pattern<TIndex, Swift<SwiftSemiGlobal> >
..param.TIndex: A q-gram index of needle(s).
...type:Spec.Index_QGram
*/
///.Spec.Swift.param.TSpec.type:Spec.SwiftSemiGlobal

template < typename TObject, typename TSpec > class Index;
template < typename TObject > struct SAValue;

struct _SwiftLocal;
typedef Tag<_SwiftLocal> SwiftLocal;

template <typename TSpec = void>
struct _SwiftSemiGlobal;
typedef Tag<_SwiftSemiGlobal<void> > SwiftSemiGlobal;

struct _Hamming;
typedef Tag<_SwiftSemiGlobal<_Hamming> > SwiftSemiGlobalHamming;


template <typename TSpec = SwiftSemiGlobal>
struct Swift;

template <>
struct Swift<SwiftSemiGlobal> {
	enum { SEMIGLOBAL = 1 };		// 0..match eps-match of min.length n0; 1..match the whole read
	enum { DIAGONAL = 1 };			// 0..use rectangular buckets (QUASAR); 1..use diagonal buckets (SWIFT)
	enum { QGRAM_ERRORS = 0 };		// q-gram must match exactly
	enum { HAMMING_ONLY = 0 };		// 0..no indels; 1..allow indels
	enum { PARAMS_BY_LENGTH = 1 };	// params are determined only by seq.length
};

template <>
struct Swift<SwiftSemiGlobalHamming> {
	enum { SEMIGLOBAL = 1 };		// 0..match eps-match of min.length n0; 1..match the whole read
	enum { DIAGONAL = 1 };			// 0..use rectangular buckets (QUASAR); 1..use diagonal buckets (SWIFT)
	enum { QGRAM_ERRORS = 0 };		// q-gram must match exactly
	enum { HAMMING_ONLY = 1 };		// 0..no indels; 1..allow indels
	enum { PARAMS_BY_LENGTH = 1 };	// params are determined only by seq.length
};

template <>
struct Swift<SwiftLocal> {
	enum { SEMIGLOBAL = 0 };		// 0..match eps-match of min.length n0; 1..match the whole read
	enum { DIAGONAL = 1 };			// 0..use rectangular buckets (QUASAR); 1..use diagonal buckets (SWIFT)
	enum { QGRAM_ERRORS = 0 };		// allow 0 errors per q-gram
	enum { HAMMING_ONLY = 0 };		// 0..no indels; 1..allow indels
	enum { PARAMS_BY_LENGTH = 0 };	// params are determined by seq.no.
};

struct SwiftParameters {
	int minThreshold;
	int minLog2Delta;
	int tabooLength;
	bool printDots;

	SwiftParameters():
		minThreshold(1),		// set minimal threshold to 1
		minLog2Delta(4),		// set minimal delta to 16
		tabooLength(1),			// minimal genomic distance between q-gram hits
		printDots(false) {}		// print a . for every 100kbps mapped genome
};

//////////////////////////////////////////////////////////////////////////////

	template <typename TSpec, typename _TSize, typename _TShortSize = unsigned short>
	struct _SwiftBucket 
	{
		typedef _TSize			TSize;
		typedef _TShortSize		TShortSize;

		TSize					firstIncrement;
		TSize					lastIncrement;
		TShortSize				counter;		// q-gram hits
		TShortSize				threshold;		// at least threshold q-gram hits induce an approx match
#ifdef SEQAN_DEBUG_SWIFT
		TSize					_lastIncDiag;
#endif
	};

	template <typename _TSpec, typename _TSize, typename _TShortSize>
	struct _SwiftBucket<_SwiftSemiGlobal<_TSpec>, _TSize, _TShortSize>
	{
		typedef _TSize			TSize;
		typedef _TShortSize		TShortSize;

		TSize					lastIncrement;
		TShortSize				counter;		// q-gram hits
		TShortSize				threshold;		// at least threshold q-gram hits induce an approx match
#ifdef SEQAN_DEBUG_SWIFT
		int						_lastIncDiag;
#endif
	};

	template <typename TSpec, typename _TSize, typename _TShortSize = unsigned short>
	struct _SwiftBucketParams 
	{
		typedef _TSize			TSize;
		typedef _TShortSize		TShortSize;

		TSize			firstBucket;	// first _SwiftBucket entry in pattern.buckets
		TSize			reuseMask;		// 2^ceil(log2(x)) reuse every x-th bucket)
		TShortSize		threshold;		// at least threshold q-gram hits induce an approx match
		TShortSize		distanceCut;	// if lastIncrement is this far or farer away, threshold can't be reached
		TShortSize		delta;			// buckets begin at multiples of delta
		TShortSize		overlap;		// number of diagonals/columns a bucket shares with its neighbor
		TShortSize		tabooLength;	// minimal genomic distance between q-gram hits
		unsigned char	logDelta;		// log2(delta)
	};

	template <typename _TSpec, typename _TSize, typename _TShortSize>
	struct _SwiftBucketParams< Swift<Tag<_SwiftSemiGlobal<_TSpec> > >, _TSize, _TShortSize>
	{
		typedef _TSize			TSize;
		typedef _TShortSize		TShortSize;

		TSize			firstBucket;	// first _SwiftBucket entry in pattern.buckets
		TSize			reuseMask;		// 2^ceil(log2(x)) reuse every x-th bucket)
		TShortSize		threshold;		// at least threshold q-gram hits induce an approx match
		TShortSize		delta;			// buckets begin at multiples of delta
		TShortSize		overlap;		// number of diagonals/columns a bucket shares with its neighbor
		TShortSize		tabooLength;	// minimal genomic distance between q-gram hits
		unsigned char	logDelta;		// log2(delta)
	};

//____________________________________________________________________________


	template <typename TSpec, typename THstkPos>
	struct _SwiftHit 
	{
		THstkPos	hstkPos;			// parallelogram begin in haystack 
		unsigned	ndlSeqNo;			// needle sequence number
		THstkPos	ndlPos;				// begin position of hit in needle
		unsigned	bucketWidth;		// (non-diagonal) bucket width (hitLengthNeedle + delta + overlap (for diagonals))
		unsigned	hitLengthNeedle;	// length of the hit in needle
	};

	template <typename TSpec, typename THstkPos>
	struct _SwiftHit<Tag<_SwiftSemiGlobal<TSpec> >, THstkPos>
	{
		THstkPos	hstkPos;			// parallelogram begin in haystack 
		unsigned	ndlSeqNo;			// needle sequence number
		unsigned	bucketWidth;		// (non-diagonal) bucket width (bktHeight + delta + overlap (for diagonals))
	};

//____________________________________________________________________________


	template <typename THaystack, typename TSpec>
	class Finder<THaystack, Swift<TSpec> >
	{
	public:
		typedef typename Iterator<THaystack, Rooted>::Type			TIterator;
		typedef typename Position<THaystack>::Type					THstkPos;
		typedef _SwiftHit<TSpec, __int64>							TSwiftHit;
		typedef String<TSwiftHit>									THitString;
		typedef typename Iterator<THitString, Standard>::Type		THitIterator;
		typedef typename SAValue<THaystack>::Type					TSAValue;
		typedef Repeat<TSAValue, unsigned>							TRepeat;
		typedef String<TRepeat>										TRepeatString;
		typedef typename Iterator<TRepeatString, Standard>::Type	TRepeatIterator;

		TIterator		data_iterator;
		TIterator		haystackEnd;
		bool			_needReinit;	// if true, the Pattern needs to be reinitialized
		THitString		hits;
		THitIterator	curHit, endHit;
		THstkPos		startPos, curPos, endPos;
		THstkPos		dotPos, dotPos2;
		TRepeatString	data_repeats;
		TRepeatIterator	curRepeat, endRepeat;

		Finder():
			_needReinit(true) { }

		Finder(THaystack &haystack):
			data_iterator(begin(haystack, Rooted())),
			_needReinit(true) { }

		template <typename TRepeatSize, typename TPeriodSize>
		Finder(THaystack &haystack, TRepeatSize minRepeatLen, TPeriodSize maxPeriod):
			data_iterator(begin(haystack, Rooted())),
			_needReinit(true) 
		{
			findRepeats(data_repeats, haystack, minRepeatLen, maxPeriod);
		}

		Finder(TIterator &iter):
			data_iterator(iter),
			_needReinit(true) { }

		Finder(TIterator const &iter):
			data_iterator(iter),
			_needReinit(true) { }

		Finder(Finder const &orig):
			data_iterator(orig.data_iterator),
			haystackEnd(orig.haystackEnd),
			_needReinit(orig._needReinit),
			hits(orig.hits),
            startPos(orig.startPos),
            curPos(orig.curPos),
            endPos(orig.endPos),
            dotPos(orig.dotPos),
            dotPos2(orig.dotPos2),
			data_repeats(orig.data_repeats)
		{
			curHit = begin(hits, Standard()) + (orig.curHit - begin(orig.hits, Standard()));
			endHit = end(hits, Standard());
			curRepeat = begin(data_repeats, Standard()) + (orig.curRepeat - begin(orig.data_repeats, Standard()));
			endRepeat = end(data_repeats, Standard());
		};

		inline typename Reference<TIterator>::Type 
		operator* () { return value(hostIterator(*this)); }

		inline typename Reference<TIterator const>::Type 
		operator* () const { return value(hostIterator(*this)); }

		operator TIterator () const	{ return data_iterator;	}
        
        Finder & operator = (Finder const &orig) 
        {
            data_iterator = orig.data_iterator;
            haystackEnd = orig.haystackEnd;
            _needReinit = orig._needReinit;
            hits = orig.hits;
            startPos = orig.startPos;
            curPos = orig.curPos;
            endPos = orig.endPos;
            dotPos = orig.dotPos;
            dotPos2 = orig.dotPos2;
            data_repeats = orig.data_repeats;
            curHit = begin(hits, Standard()) + (orig.curHit - begin(orig.hits, Standard()));
            endHit = end(hits, Standard());
            curRepeat = begin(data_repeats, Standard()) + (orig.curRepeat - begin(orig.data_repeats, Standard()));
            endRepeat = end(data_repeats, Standard());
            return *this;
        }
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
		typedef _SwiftHit<TSpec, __int64>						TSwiftHit;
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


	template <typename TIndex, typename TSpec>
	class Pattern<TIndex, Swift<TSpec> >
	{
	public:
		typedef typename Size<TIndex>::Type								TSize;
		typedef unsigned												TShortSize;
		typedef typename Fibre<TIndex, Tag<_Fibre_SA> const >::Type		TSA;
		typedef typename Fibre<TIndex, Tag<_Fibre_Shape> const >::Type	TShape;
		typedef typename Iterator<TSA const, Standard>::Type			TIterator;
		
		typedef _SwiftBucket<TSpec, TSize, TShortSize>					TBucket;
		typedef String<TBucket>											TBucketString;
		typedef _SwiftBucketParams<TSpec, TSize, TShortSize>			TBucketParams;
		typedef String<TBucketParams>									TBucketParamsString;
		
		TShape					shape;
		TBucketString			buckets;
		TBucketParamsString		bucketParams;
		SwiftParameters			params;
		unsigned				curSeqNo;
        __int64                 curBeginPos, curEndPos;
		__int64					finderPosOffset;
		__int64					finderPosNextOffset;
		__int64					finderLength;

		double					_currentErrorRate;
		int						_currentMinLengthForAll;

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


template <typename _TSpec, typename TSize, typename TShortSize>
inline void _printSwiftParams(_SwiftBucketParams<_TSpec, TSize, TShortSize > &bucketParams)
{
	::std::cout << "  firstBucket: " << bucketParams.firstBucket << ::std::endl;
	::std::cout << "  reuseMask:   " << bucketParams.reuseMask << ::std::endl;
	::std::cout << "  distanceCut: " << bucketParams.distanceCut << ::std::endl;
	::std::cout << "  delta:       " << bucketParams.delta << ::std::endl;
	::std::cout << "  threshold:   " << bucketParams.threshold << ::std::endl;
	::std::cout << "  overlap:     " << bucketParams.overlap << ::std::endl;
	::std::cout << "  logDelta:    " << (int)bucketParams.logDelta << ::std::endl << ::std::endl;
}

template <typename _TSpec, typename TSize, typename TShortSize>
inline void _printSwiftParams(_SwiftBucketParams<Tag<_SwiftSemiGlobal<_TSpec> >, TSize, TShortSize > &bucketParams)
{
	::std::cout << "  firstBucket: " << bucketParams.firstBucket << ::std::endl;
	::std::cout << "  reuseMask:   " << bucketParams.reuseMask << ::std::endl;
	::std::cout << "  delta:       " << bucketParams.delta << ::std::endl;
	::std::cout << "  threshold:   " << bucketParams.threshold << ::std::endl;
	::std::cout << "  overlap:     " << bucketParams.overlap << ::std::endl;
	::std::cout << "  logDelta:    " << (int)bucketParams.logDelta << ::std::endl << ::std::endl;
}

template <typename TIndex, typename TSpec>
inline void _printSwiftBuckets(Pattern< TIndex, Swift<TSpec> > &p)
{
	typedef typename Pattern<TIndex, Swift<TSpec> >::TBucketParams TParams;

	unsigned j = 0;
	TParams *bucketParams = &_swiftBucketParams(p, 0);

	for(unsigned i=0; i<length(p.buckets) && i<10; ++i) 
	{
		if ((i & bucketParams->reuseMask) == 0)
		{
			::std::cout << ::std::endl << "ReadBucket #" << j << "    " << '"';
			::std::cout << indexText(host(p))[j] << '"' << ::std::endl;
			::std::cout << "  length:      " << sequenceLength(j, host(p)) << ::std::endl;
			bucketParams = &_swiftBucketParams(p, j++);
			_printSwiftParams(*bucketParams);
		}

		::std::cout << "    lastInc: " << (int)p.buckets[i].lastIncrement;
		::std::cout << "  \tCounter: " << p.buckets[i].counter << ::std::endl;
	}
}

template <typename TIndex, typename TSpec, typename TSize>
inline typename Pattern<TIndex, Swift<TSpec> >::TBucketParams &
_swiftBucketParams(Pattern<TIndex, Swift<TSpec> > & pattern, TSize seqNo) 
{
	if (Swift<TSpec>::PARAMS_BY_LENGTH)
		return pattern.bucketParams[sequenceLength(seqNo, host(pattern))];
	else
		return pattern.bucketParams[seqNo];
}

template <typename TIndex, typename TSpec, typename TParams, typename TSize>
inline unsigned
_swiftBucketNo(Pattern<TIndex, Swift<TSpec> > const &, TParams &bucketParams, TSize seqNo) 
{
	if (Swift<TSpec>::PARAMS_BY_LENGTH)
		return (bucketParams.reuseMask + 1) * seqNo;	// assumes the same reuseMask for all reads
	else
		return bucketParams.firstBucket;
}

template <typename TIndex, typename TSpec, typename TSeqNo>
inline int
_qgramLemma(Pattern<TIndex, Swift<TSpec> > const & pattern, TSeqNo seqNo, int errors)
{
	// q-gram lemma: How many conserved q-grams we see at least?
	// each error destroys at most <weight> many (gapped) q-grams
	return 
		sequenceLength(seqNo, host(pattern)) - length(indexShape(host(pattern))) + 1 
		- errors * weight(indexShape(host(pattern)));
}

template <typename TIndex, typename TSpec, typename TSeqNo, typename TThreshold>
inline void
setMinThreshold(Pattern<TIndex, Swift<TSpec> > & pattern, TSeqNo seqNo, TThreshold thresh) 
{
	typedef Pattern<TIndex, Swift<TSpec> >						TPattern;
	typedef typename TPattern::TBucketParams					TBucketParams;
	typedef typename TPattern::TBucketString					TBucketString;
	typedef typename Iterator<TBucketString, Standard>::Type	TBucketIterator;

	TBucketParams &bucketParams = _swiftBucketParams(pattern, seqNo);
	TBucketIterator it = begin(pattern.buckets, Standard()) + _swiftBucketNo(pattern, bucketParams, seqNo);
	TBucketIterator itEnd = it + (bucketParams.reuseMask + 1);

	for (; it != itEnd; ++it)
		if ((*it).threshold < thresh)
			(*it).threshold = thresh;
}


template <typename TIndex, typename TFloat, typename _TSize, typename TSpec>
inline void _patternInit(Pattern<TIndex, Swift<TSpec> > &pattern, TFloat errorRate, _TSize minLengthForAll) 
{
	typedef Pattern<TIndex, Swift<TSpec> >						TPattern;
	typedef typename Size<TIndex>::Type							TSize;
	typedef typename Fibre<TIndex, QGram_SA>::Type				TSA;
	typedef typename Iterator<TSA, Standard>::Type				TSAIter;
	typedef typename TPattern::TBucket							TBucket;
	typedef typename TBucket::TSize								TBucketSize;
	typedef typename TPattern::TBucketParams					TBucketParams;
	typedef typename TPattern::TBucketString					TBucketString;
	typedef typename Iterator<TBucketString, Standard>::Type	TBucketIterator;

	double _newErrorRate = errorRate;
	TSize seqCount = countSequences(host(pattern));
	
	if (pattern._currentErrorRate != _newErrorRate || pattern._currentMinLengthForAll != minLengthForAll)
	{
		// settings have been changed -> initialize bucket parameters
	
		pattern._currentErrorRate = _newErrorRate;
		pattern._currentMinLengthForAll = minLengthForAll;
		pattern.finderPosOffset = 0;
		pattern.finderPosNextOffset = pattern.finderLength;

		indexRequire(host(pattern), QGram_SADir());
		pattern.shape = indexShape(host(pattern));

		TSize span = length(pattern.shape);
		TSize count = 0;
		TSize bucketsPerCol2Max = 0;
		TSize maxLength = 0;
		
		if (Swift<TSpec>::PARAMS_BY_LENGTH) {
			for(unsigned seqNo = 0; seqNo < seqCount; ++seqNo) {
				TSize length = sequenceLength(seqNo, host(pattern));
				if (maxLength < length)
					maxLength = length;
			}
			resize(pattern.bucketParams, maxLength + 1);
		} else
			resize(pattern.bucketParams, seqCount);
		
    	if (Swift<TSpec>::SEMIGLOBAL == 0) 
	    {
		    // global matches
    		TSize minLength = minLengthForAll;
	    	for(unsigned seqNo = 0; seqNo < seqCount; ++seqNo) 
		    {
    			// swift q-gram lemma
	    		TBucketParams &bucketParams = _swiftBucketParams(pattern, seqNo);
		    	// n..next length that could decrease threshold
			    TSize n = (TSize) ceil((floor(errorRate * minLength) + 1) / errorRate);
    			// minimal threshold is minimum errors of minLength and n
	    		int threshold = (TSize) _min(
		    		(n + 1) - span * (floor(errorRate * n) + 1),
			    	(minLength + 1) - span * (floor(errorRate * minLength) + 1));

				if (threshold > pattern.params.minThreshold)
					bucketParams.threshold = threshold;
				else
					bucketParams.threshold = pattern.params.minThreshold;

			    TSize errors = (TSize) floor((2 * bucketParams.threshold + span - 3) / (1 / errorRate - span));
			
			
    			// a bucket has distanceCut different positions of q-grams
	    		// if a q-gram is this far or further away it can't belong to the
		    	// same bucket
			    bucketParams.distanceCut = (bucketParams.threshold - 1) + span * errors + span;

    			TSize bucketsPerCol2;
	    		if(Swift<TSpec>::DIAGONAL == 1) 
		    	{
			    	// Use overlapping parallelograms
				    bucketParams.overlap = errors;
    
	    			// delta must be a power of 2 and greater than errors
		    		bucketParams.logDelta = (TSize) ceil(log((double)errors + 1) / log(2.0));
			    	if (bucketParams.logDelta < pattern.params.minLog2Delta) 
				    	bucketParams.logDelta = pattern.params.minLog2Delta;
    				bucketParams.delta = 1 << bucketParams.logDelta;
	    			bucketParams.tabooLength = pattern.params.tabooLength;

		    		// maximal number of buckets in one column
			    	TSize bucketsPerCol = (sequenceLength(seqNo, host(pattern)) - span + 2 * bucketParams.delta + errors - 1) / bucketParams.delta;
				    bucketsPerCol2 = 1 << (TSize) ceil(log((double)bucketsPerCol) / log(2.0)); // next greater or equal power of 2
    			}
	    		else
		    	{
			    	// TODO: classical swift for rectangular buckets
				    // Use overlapping rectangles
    				//bucketParams.overlap = ;

	    			// delta must be a power of 2 greater than seq.length + errors (define a minimal delta of 32)
		    		//bucketParams.logDelta = ;
			    	//if (bucketParams.logDelta < pattern.params.minLog2Delta) 
				    //	bucketParams.logDelta = pattern.params.minLog2Delta;
    				//bucketParams.delta = 1 << bucketParams.logDelta;
	    			bucketsPerCol2 = 1;
		    	}
			
    			bucketParams.firstBucket = count; // firstBucket is only used if Swift<TSpec>::PARAMS_BY_LENGTH == 0
	    		bucketParams.reuseMask = bucketsPerCol2 - 1;
		    	bucketParams.tabooLength = pattern.params.tabooLength;
			
    			if (Swift<TSpec>::PARAMS_BY_LENGTH) {
	    			++count;
		    		if (bucketsPerCol2Max < bucketsPerCol2)
			    		bucketsPerCol2Max = bucketsPerCol2;
    			} else
	    			count += bucketsPerCol2;
		
		    	if (seqNo<3)
			    	_printSwiftParams(bucketParams);
		    }
		} else
			for(unsigned seqNo = 0; seqNo < seqCount; ++seqNo) 
			{
				// get pattern length and max. allowed errors
				TBucketParams &bucketParams = _swiftBucketParams(pattern, seqNo);
				TSize length = sequenceLength(seqNo, host(pattern));
				TSize errors = (TSize) floor(errorRate * length);
				TSize errorsWC = errors / (1 + Swift<TSpec>::QGRAM_ERRORS);

				// q-gram lemma: How many conserved q-grams we see at least?
				// (define a minimal threshold of 1)
				int threshold = length - span + 1 - errorsWC * weight(pattern.shape);
				if (threshold > pattern.params.minThreshold)
					bucketParams.threshold = threshold;
				else
					bucketParams.threshold = pattern.params.minThreshold;
				
				if (Swift<TSpec>::HAMMING_ONLY != 0)
					errors = 0;			

				// a bucket has distanceCut different positions of q-grams
				// if a q-gram is this far or farer away it can't belong to the
				// same bucket
	//			bucketParams.distanceCut = length - (span - 1) + errors;

				TSize bucketsPerCol2;
				if (Swift<TSpec>::DIAGONAL == 1)
				{
					// Use overlapping parallelograms				
					bucketParams.overlap = errors;
					
					// delta must be a power of 2 greater then errors (define a minimal delta of 8)
					bucketParams.logDelta = (TSize) ceil(log((double)(errors + 1)) / log(2.0));
					if (bucketParams.logDelta < pattern.params.minLog2Delta) 
						bucketParams.logDelta = pattern.params.minLog2Delta;
					bucketParams.delta = 1 << bucketParams.logDelta;

					// the formula for bucketsPerCol is (worst-case):
					// (height-(q-1) - 1 - (delta+1-e))/delta + 3
					//    ^-- full paral. in the middle --^     ^-- 2 at the bottom, 1 at the top
					TSize bucketsPerCol = (length - span + 2 * bucketParams.delta + errors - 1) / bucketParams.delta;
					bucketsPerCol2 = 1 << (TSize) ceil(log((double)bucketsPerCol) / log(2.0));
				}
				else
				{
					// Use overlapping rectangles
					bucketParams.overlap = length - span + errors;

					// delta must be a power of 2 greater then seq.length + errors (define a minimal delta of 32)
					bucketParams.logDelta = (TSize) ceil(log((double)(length - span + 1 + errors)) / log(2.0));
					if (bucketParams.logDelta < pattern.params.minLog2Delta) 
						bucketParams.logDelta = pattern.params.minLog2Delta;
					bucketParams.delta = 1 << bucketParams.logDelta;

					bucketsPerCol2 = 2;
				}

	//			SEQAN_ASSERT(distanceCut <= bucketsPerCol * (TSize) delta);

				bucketParams.firstBucket = count;
				bucketParams.reuseMask = bucketsPerCol2 - 1;
				bucketParams.tabooLength = pattern.params.tabooLength;
				
				if (Swift<TSpec>::PARAMS_BY_LENGTH) {
					++count;
					if (bucketsPerCol2Max < bucketsPerCol2)
						bucketsPerCol2Max = bucketsPerCol2;
				} else
					count += bucketsPerCol2;
				
	/*			if (seqNo<3)
					_printSwiftParams(bucketParams);
	*/		}

		if (Swift<TSpec>::PARAMS_BY_LENGTH) {
			count *= bucketsPerCol2Max;
			for(unsigned i = 0; i < length(pattern.bucketParams); ++i)
				pattern.bucketParams[i].reuseMask = bucketsPerCol2Max - 1;
		}
		resize(pattern.buckets, count);

		TBucketIterator	bktEnd, bkt = begin(pattern.buckets, Standard());
		for(unsigned seqNo = 0; seqNo < seqCount; ++seqNo)
		{
			TBucketParams &bucketParams = _swiftBucketParams(pattern, seqNo);
			TBucketSize lastIncrement = (TBucketSize)0 - (TBucketSize)bucketParams.tabooLength;
			for(bktEnd = bkt + bucketParams.reuseMask + 1; bkt != bktEnd; ++bkt) 
			{
				(*bkt).lastIncrement = lastIncrement;
				(*bkt).counter = 0;
				(*bkt).threshold = bucketParams.threshold;
			}
		}
	}
	else
	{
		// settings are unchanged -> reset buckets

		// finderPosOffset is used to circumvent expensive resetting of all buckets
		pattern.finderPosOffset = pattern.finderPosNextOffset;
		pattern.finderPosNextOffset += pattern.finderLength;
		
		// reset buckets only if an overflow of the finder position would occur
		if (pattern.finderPosNextOffset <= pattern.finderPosOffset)
		{
			pattern.finderPosOffset = 0;
			pattern.finderPosNextOffset = pattern.finderLength;
			
			TBucketIterator	bktEnd, bkt = begin(pattern.buckets, Standard());
			for (TSize ndlSeqNo = 0; ndlSeqNo < seqCount; ++ndlSeqNo) 
			{
				TBucketParams &bucketParams = _swiftBucketParams(pattern, ndlSeqNo);
				TBucketSize lastIncrement = (TBucketSize)0 - (TBucketSize)bucketParams.tabooLength;
				for (bktEnd = bkt + (bucketParams.reuseMask + 1); bkt != bktEnd; ++bkt)
				{
					(*bkt).lastIncrement = lastIncrement;
					(*bkt).counter = 0;
				}
			}
		}
	}

/*
	std::cerr << "Swift bucket params: " << length(pattern.bucketParams) << std::endl;
	std::cerr << "Swift buckets:       " << length(pattern.buckets) << std::endl;
	std::cerr << "Buckets per read:    " << bucketsPerCol2Max << std::endl;
*/
}


/////////////////////////////////////////////////////////////
// Creates a new hit and appends it to the finders hit list
template <
	typename THaystack,
	typename TIndex,
	typename TSpec,
	typename TBucket,
	typename TBucketParams,
	typename TSize
>
inline void _createHit(
	Finder<THaystack, Swift<TSpec> > & finder,
	Pattern<TIndex, Swift<TSpec> > & pattern,
	TBucket & bkt,
	TBucketParams & bucketParams,
	__int64 diag,
	TSize ndlSeqNo)
{
	typedef typename Finder<THaystack, Swift<TSpec> >::TSwiftHit	THit;
	__int64 lastInc = (__int64)(*bkt).lastIncrement - pattern.finderPosOffset;
	__int64 firstInc = (__int64)(*bkt).firstIncrement - pattern.finderPosOffset;

    if(diag > lastInc) 
	{
		// bucket is reused since last increment 
		TSize reusePos = (bucketParams.reuseMask + 1) << bucketParams.logDelta;
		diag -= (__int64)ceil((diag-lastInc)/(double)reusePos) * reusePos;
	}

    // determine width, height, and begin position in needle
	TSize width = lastInc - firstInc + length(pattern.shape);
	TSize height = width + bucketParams.delta + bucketParams.overlap;
	__int64 ndlBegin = lastInc + length(pattern.shape) - diag - height;

    // create the hit
	THit hit = {                //                              *
		firstInc,               // bucket begin in haystack     * *
		ndlSeqNo,               // needle seq. number           *   *
		ndlBegin,               // bucket begin in needle       *     *
		width,                  // bucket width (non-diagonal)    *   *
		height                  // bucket height                    * * 
	};                          //                                    *

    // append the hit to the finders hit list
    appendValue(finder.hits, hit);
}

//////////////////////////////////////////////////////////////////////
// Updates the counters of the buckets in which the q-gram with hash value hash occurs.
// Assures that those updated bucket counters are set to one
//    - that exceeded the reuse mask since last increment
//    - for which the last increment lies more than distanceCut away.
// If a bucket counter reaches threshold a hit is appended to the hit-list of the finder.
// Returns true if the hit-list of the finder is not empty after this calls.
template <
	typename TFinder,
	typename TIndex,
	typename TSpec,
	typename THashValue
>
inline bool _swiftMultiProcessQGram(
	TFinder & finder,
	Pattern<TIndex, Swift<TSpec> > & pattern,
	THashValue hash)
{
	typedef Pattern<TIndex, Swift<TSpec> >						TPattern;

	typedef typename Size<TIndex>::Type							TSize;
	typedef typename Fibre<TIndex, QGram_SA>::Type				TSA;
	typedef typename Iterator<TSA, Standard>::Type				TSAIter;
	typedef typename TPattern::TBucketString					TBucketString;
	typedef typename Iterator<TBucketString, Standard>::Type	TBucketIter;
	typedef typename Value<TBucketString>::Type					TBucket;
	typedef typename TBucket::TShortSize						TShortSize;
	typedef typename TPattern::TBucketParams					TBucketParams;
	typedef typename TFinder::TSwiftHit							THit;
	
	TIndex const &index = host(pattern);
	
	// create an iterator over the positions of the q-gram occurences in pattern
	TSAIter saBegin = begin(indexSA(index), Standard());
	TSAIter occ = saBegin + indexDir(index)[getBucket(index.bucketMap, hash)];
    TSAIter occEnd = saBegin + indexDir(index)[getBucket(index.bucketMap, hash) + 1];
    TBucketIter bktBegin = begin(pattern.buckets, Standard());
	Pair<unsigned> ndlPos;
	
/*	std::cerr<<"\t["<<(occEnd-occ)<<"]"<< std::flush;
	
	if ((occEnd-occ)>100)
	{
		std::cerr<<" ";
		for(int i=0;i<length(indexShape(host(pattern)));++i)
			std::cerr<<*(hostIterator(hostIterator(finder))+i);
	}
*/	
	// iterate over all q-gram occurences and do the processing
	__int64 curPos = finder.curPos + pattern.finderPosOffset;
	for(; occ != occEnd; ++occ)
	{
		posLocalize(ndlPos, *occ, stringSetLimits(index)); // get pair of SeqNo and Pos in needle
		TBucketParams &bucketParams = _swiftBucketParams(pattern, getSeqNo(ndlPos));

		// begin position of the diagonal of q-gram occurence in haystack (possibly negative)
		__int64 diag = finder.curPos;
		if (Swift<TSpec>::DIAGONAL == 1) diag -= getSeqOffset(ndlPos);

		unsigned bktNo = (diag >> bucketParams.logDelta) & bucketParams.reuseMask; // bucket no of diagonal
		unsigned bktOfs = diag & (bucketParams.delta - 1); // offset of diagonal to bucket begin
		__int64  bktBeginHstk = diag & ~(__int64)(bucketParams.delta - 1); // haystack position of bucket begin diagonal
		
		// global (over all pattern sequences) number of current bucket
		TBucketIter bkt = bktBegin + (_swiftBucketNo(pattern, bucketParams, getSeqNo(ndlPos)) + bktNo);
		
		TShortSize hitCount;

		do {
			if ((__int64)(*bkt).lastIncrement < bktBeginHstk + pattern.finderPosOffset
				|| (__int64)((*bkt).lastIncrement + bucketParams.distanceCut) < (__int64)(curPos + length(pattern.shape)))
			{
				// last increment was before the beginning of the current bucket => bucket is reused
				// OR last increment was in the same bucket but lies more than distanceCut away
                
                if ((*bkt).counter >= (*bkt).threshold)
                {
	                // create a new hit and append it to the finders hit list
				    _createHit(finder, pattern, bkt, bucketParams, bktBeginHstk, getSeqNo(ndlPos));
                }

				// reuse bucket
				hitCount = 1;
				(*bkt).firstIncrement = curPos;
			}
			else if((*bkt).lastIncrement + bucketParams.tabooLength > curPos)
			{
				// bkt counter was already incremented for another q-gram at
				//   a haystack position that is closer than tabooLength
				// we jump directly to
				//   where we check whether the q-gram falls into another overlapping bucket or not
				goto checkOverlap;
			}
			else
			{
				if((*bkt).counter == 0) (*bkt).firstIncrement = curPos;
				hitCount = (*bkt).counter + 1;
			}

			(*bkt).lastIncrement = curPos;
			(*bkt).counter = hitCount;
#ifdef SEQAN_DEBUG_SWIFT
			(*bkt)._lastIncDiag = diag;
#endif

checkOverlap:
			// check if q-gram falls into another overlapping bucket
			if(bktOfs >= bucketParams.overlap) break;

			// set to previous overlapping bucket for next iteration
			bktBeginHstk -= bucketParams.delta;
			bktOfs += bucketParams.delta;
			if(bktNo) {
				--bktNo;
				--bkt;
			} else {
				bktNo = bucketParams.reuseMask;
				bkt += bktNo;
			}
		}
		while(true);
	}

	finder.curHit = begin(finder.hits, Standard());
	finder.endHit = end(finder.hits, Standard());

	return !empty(finder.hits);
}

///////////////////////////////////////////////////////////////////
// Updates the counters of the buckets in which the q-gram with hash value hash occurs.
// Assures that those updated bucket counters that exceeded the reuse mask since last increment are set to one.
// If a bucket counter reaches threshold a hit is appended to the hit-list of the finder.
// Returns true if the hit-list of the finder is not empty after this call.
template <
	typename TFinder,
	typename TIndex, 
	typename _TSpec,
	typename THValue
>
inline bool _swiftMultiProcessQGram(
	TFinder &finder, 
	Pattern<TIndex, Swift<Tag<_SwiftSemiGlobal<_TSpec> > > > &pattern,
	THValue hash)
{
	typedef Pattern<TIndex, Swift<Tag<_SwiftSemiGlobal<_TSpec> > > >	TPattern;

	typedef typename Size<TIndex>::Type							TSize;
	typedef typename Fibre<TIndex, QGram_SA>::Type				TSA;
	typedef typename Iterator<TSA, Standard>::Type				TSAIter;
	typedef typename TPattern::TBucketString					TBucketString;
	typedef typename Iterator<TBucketString, Standard>::Type	TBucketIter;
	typedef typename Value<TBucketString>::Type					TBucket;
	typedef typename TBucket::TShortSize						TShortSize;
	typedef typename TPattern::TBucketParams					TBucketParams;
	typedef typename TFinder::TSwiftHit							THit;
	
	TIndex const &index = host(pattern);	
	
	// create an iterator over the positions of the q-gram occurences in pattern
	TSAIter saBegin = begin(indexSA(index), Standard());
	TSAIter occ = saBegin + indexDir(index)[getBucket(index.bucketMap, hash)];
    TSAIter occEnd = saBegin + indexDir(index)[getBucket(index.bucketMap, hash) + 1];
    TBucketIter bktBegin = begin(pattern.buckets, Standard());
	Pair<unsigned> ndlPos;
	
/*	std::cerr<<"\t["<<(occEnd-occ)<<"]"<< std::flush;
	
	if ((occEnd-occ)>100)
	{
		std::cerr<<" ";
		for(int i=0;i<length(indexShape(host(pattern)));++i)
			std::cerr<<*(hostIterator(hostIterator(finder))+i);
	}
*/	
	// iterate over all q-gram occurences and do the processing
	__int64 curPos = finder.curPos + pattern.finderPosOffset;
	for(; occ != occEnd; ++occ) 
	{
		posLocalize(ndlPos, *occ, stringSetLimits(index));
		TBucketParams &bucketParams = _swiftBucketParams(pattern, getSeqNo(ndlPos));

		__int64 diag = finder.curPos;
		if (Swift<Tag<_SwiftSemiGlobal<_TSpec> > >::DIAGONAL == 1) diag -= getSeqOffset(ndlPos);
		
		unsigned bktNo = (diag >> bucketParams.logDelta) & bucketParams.reuseMask;
		unsigned bktOfs = diag & (bucketParams.delta - 1);
		__int64  bktBeginHstk = diag & ~(__int64)(bucketParams.delta - 1);

		TBucketIter bkt = bktBegin + (_swiftBucketNo(pattern, bucketParams, getSeqNo(ndlPos)) + bktNo);		
		TShortSize hitCount;

		do 
		{
			if ((__int64)(*bkt).lastIncrement < bktBeginHstk + pattern.finderPosOffset)
			{
				// last increment was before the beginning of the current bucket
				// (we must ensure that bucketIdx doesn't collide)
				hitCount = 1;
			}
			else
			{
				if ((__int64)((*bkt).lastIncrement + bucketParams.tabooLength) > curPos)
					goto checkOverlap;	// increment only once per sequence			
				hitCount = (*bkt).counter + 1;
			}

			(*bkt).lastIncrement = curPos;
			(*bkt).counter = hitCount;
#ifdef SEQAN_DEBUG_SWIFT
			(*bkt)._lastIncDiag = diag;
#endif

			if (hitCount == (*bkt).threshold)
			{

				TSize height = 0;
				if (Swift<Tag<_SwiftSemiGlobal<_TSpec> > >::DIAGONAL == 1)
					height = sequenceLength(getSeqNo(ndlPos), host(pattern)) - 1;

#ifdef SEQAN_DEBUG_SWIFT
				// upper bucket no. of lastIncr. q-gram
				__int64 upperBktNo = ((*bkt).lastIncrement - pattern.finderPosOffset) >> bucketParams.logDelta;

				// we must decrement bucket no. until (no. mod reuse == bktNo)
				__int64 _bktBeginHstk = 
					 (upperBktNo - ((upperBktNo - bktNo) & bucketParams.reuseMask)) << bucketParams.logDelta;

				if ((*bkt)._lastIncDiag - _bktBeginHstk >= bucketParams.delta + bucketParams.overlap || (*bkt)._lastIncDiag < _bktBeginHstk) {
					::std::cerr << "qgram stored in wrong bucket (diag:" << (*bkt)._lastIncDiag << ", begin:" << _bktBeginHstk;
					::std::cerr << ", delta:" << bucketParams.delta << ", overlap:" << bucketParams.overlap << ")" << ::std::endl;
				}
#endif
//				if (bktBeginHstk >= 0) 
//				{
					THit hit = {
						bktBeginHstk,										// bucket begin in haystack
						getSeqNo(ndlPos),									// needle seq. number
						height + bucketParams.delta + bucketParams.overlap	// bucket width (non-diagonal)
					};
					appendValue(finder.hits, hit);
//				} else {
//					// match begins left of haystack begin
//					THit hit = {
//						0,													// bucket begin in haystack
//						getSeqNo(ndlPos),									// needle seq. number
//						height + bucketParams.delta + bucketParams.overlap	// bucket width (non-diagonal)
//						+ (diag & ~(__int64)(bucketParams.delta - 1))
//					};
//					appendValue(finder.hits, hit);
//				}
			}

		checkOverlap:
			if (bktOfs >= bucketParams.overlap) break;

			// repeat with the previous overlapping bucket
			bktBeginHstk -= bucketParams.delta;
			bktOfs += bucketParams.delta;
			if (bktNo) {
				--bktNo;
				--bkt;
			} else {
				bktNo = bucketParams.reuseMask;
				bkt += bktNo;
			}
		} while (true);
	}

	finder.curHit = begin(finder.hits, Standard());
	finder.endHit = end(finder.hits, Standard());

	return !empty(finder.hits);
}

//////////////////////////////////////////////////////
// resets counter and lastIncrement of all buckets
template <
	typename TFinder,
	typename TIndex, 
	typename TSpec
>
inline bool _swiftMultiFlushBuckets(
	TFinder & finder,
	Pattern<TIndex, Swift<TSpec> > & pattern
	)
{
	typedef Pattern<TIndex, Swift<TSpec> >						TPattern;

	typedef typename TPattern::TBucket							TBucket;
	typedef typename TBucket::TSize								TBucketSize;
	typedef typename TPattern::TBucketString					TBucketString;
	typedef typename Iterator<TBucketString, Standard>::Type	TBucketIterator;
	typedef typename TPattern::TBucketParams					TBucketParams;

	typedef typename Size<TIndex>::Type							TSize;

	TBucketIterator	bkt = begin(pattern.buckets, Standard());
	TBucketIterator	bktEnd;
	TSize seqCount = countSequences(host(pattern));
	__int64 hstkLength = length(haystack(finder));

	for(TSize ndlSeqNo = 0; ndlSeqNo < seqCount; ++ndlSeqNo) 
	{
		TBucketParams &bucketParams = _swiftBucketParams(pattern, ndlSeqNo);
		bktEnd = bkt + (bucketParams.reuseMask + 1);
		for(unsigned bktNo = 0; bkt != bktEnd; ++bkt, ++bktNo)
		{
            if ((*bkt).counter >= (*bkt).threshold)
            {
		    	// hstkPos / delta: gives the number of the bucket that is at the top of this column (modulo reuseMask missing)
			    TSize topBucket = (TSize)(hstkLength >> bucketParams.logDelta);
    			// number of buckets in last column above the bucket with the number bktNo
	    		TSize bucketNoInCol = (topBucket + bucketParams.reuseMask + 1 - bktNo) & bucketParams.reuseMask;
		    	// begin position of lower diagonal of this bucket in haystack (possibly negative)
			    __int64 diag = (hstkLength & ~(__int64)(bucketParams.delta - 1)) - (bucketNoInCol << bucketParams.logDelta);
			
    	        // create a new hit and append it to the finders hit list
	    		_createHit(finder, pattern, bkt, bucketParams, diag, ndlSeqNo);

		    	(*bkt).lastIncrement = (TBucketSize)0 - (TBucketSize)bucketParams.tabooLength;
			    (*bkt).counter = 0;
            }
		}
	}
	finder.curHit = begin(finder.hits, Standard());
	finder.endHit = end(finder.hits, Standard());

	return !empty(finder.hits);
}

//////////////////////////////////////////////////////
// no resetting is needed for the semiglobal version
template <
	typename TFinder,
	typename TIndex, 
	typename _TSpec
>
inline bool _swiftMultiFlushBuckets(
	TFinder &, 
	Pattern<TIndex, Swift<Tag<_SwiftSemiGlobal<_TSpec> > > > &)
{
    // there is nothing to be done here as we dump matches immediately after reaching the threshold
    return false;
}

template <typename TIndex, typename TSpec>
inline bool 
empty(Pattern<TIndex, Swift<TSpec> > & me) 
{
	return empty(me.bucketParams);
}

template <typename TIndex, typename TSpec>
inline void 
clear(Pattern<TIndex, Swift<TSpec> > & me) 
{
	me.finderPosOffset = 0;
	me.finderPosNextOffset = 0;
	me.finderLength = 0;
	me._currentErrorRate = -1;
	me._currentMinLengthForAll = -1;
	clear(me.bucketParams);
	clear(me.buckets);
}

//____________________________________________________________________________

template <typename THaystack, typename TSpec>
inline typename Position<Finder<THaystack, Swift<TSpec> > >::Type
position(Finder<THaystack, Swift<TSpec> > const & finder)
{
	typename Finder<THaystack, Swift<TSpec> >::TSwiftHit &hit = *finder.curHit;
	return hit.hstkPos + hit.bucketWidth;
}

template <typename THaystack, typename TSpec>
inline typename Position<Finder<THaystack, Swift<TSpec> > >::Type
position(Finder<THaystack, Swift<TSpec> > & finder)
{
	return position(const_cast<Finder<THaystack, Swift<TSpec> > const &>(finder));
}

//____________________________________________________________________________

template <typename TIndex, typename TSpec>
inline typename SAValue<TIndex>::Type
position(Pattern<TIndex, Swift<TSpec> > const & pattern)
{
	__int64 hitEnd = pattern.curEndPos;
	__int64 textLength = sequenceLength(pattern.curSeqNo, needle(pattern));
	if(hitEnd > textLength) hitEnd = textLength;

	typename SAValue<TIndex >::Type pos;
	posLocalToX(pos, Pair<unsigned, __int64>(pattern.curSeqNo, hitEnd), stringSetLimits(host(pattern)));
	return pos;
}

template <typename TIndex, typename TSpec>
inline typename SAValue<TIndex>::Type
position(Pattern<TIndex, Swift<Tag<_SwiftSemiGlobal<TSpec> > > > const & pattern)
{
	typedef typename Size<TIndex>::Type TSize;
    typename SAValue<TIndex >::Type pos;
	posLocalToX(pos, Pair<unsigned, TSize>(pattern.curSeqNo, length(needle(pattern))), stringSetLimits(host(pattern)));
	return pos;
}

template <typename TIndex, typename TSpec>
inline typename SAValue<TIndex>::Type
position(Pattern<TIndex, Swift<TSpec> > & pattern)
{
	return position(const_cast<Pattern<TIndex, Swift<TSpec> > const &>(pattern));
}

//____________________________________________________________________________

template <typename THaystack, typename TSpec>
inline typename Position<Finder<THaystack, Swift<TSpec> > >::Type
beginPosition(Finder<THaystack, Swift<TSpec> > const & finder)
{
	return (*finder.curHit).hstkPos;
}

template <typename THaystack, typename TSpec>
inline typename Position<Finder<THaystack, Swift<TSpec> > >::Type
beginPosition(Finder<THaystack, Swift<TSpec> > & finder)
{
	return beginPosition(const_cast<Finder<THaystack, Swift<TSpec> > const &>(finder));
}

//____________________________________________________________________________

template <typename TIndex, typename TSpec>
inline typename SAValue<TIndex>::Type
beginPosition(Pattern<TIndex, Swift<TSpec> > const & pattern)
{
	__int64 hitBegin = pattern.curBeginPos;
	if (hitBegin < 0) hitBegin = 0;
	
	typename SAValue<TIndex >::Type pos;
	posLocalToX(pos, Pair<unsigned, __int64>(pattern.curSeqNo, hitBegin), stringSetLimits(host(pattern)));
	return pos;
}

template <typename TIndex, typename TSpec>
inline typename SAValue<TIndex >::Type
beginPosition(Pattern<TIndex, Swift<Tag<_SwiftSemiGlobal<TSpec> > > > const & pattern)
{
    typename SAValue<TIndex >::Type pos;
	posLocalToX(pos, Pair<unsigned>(pattern.curSeqNo, 0), stringSetLimits(host(pattern)));
	return pos;
}

template <typename TIndex, typename TSpec>
inline typename SAValue<TIndex>::Type
beginPosition(Pattern<TIndex, Swift<TSpec> > & pattern)
{
	return beginPosition(const_cast<Pattern<TIndex, Swift<TSpec> > const &>(pattern));
}

//____________________________________________________________________________

template <typename THaystack, typename TSpec>
inline typename Position<Finder<THaystack, Swift<TSpec> > >::Type
endPosition(Finder<THaystack, Swift<TSpec> > const & finder)
{
	typename Finder<THaystack, Swift<TSpec> >::TSwiftHit &hit = *finder.curHit;
	return hit.hstkPos + hit.bucketWidth;
}

template <typename THaystack, typename TSpec>
inline typename Position<Finder<THaystack, Swift<TSpec> > >::Type
endPosition(Finder<THaystack, Swift<TSpec> > & finder)
{
	return endPosition(const_cast<Finder<THaystack, Swift<TSpec> > const &>(finder));
}

//____________________________________________________________________________

template <typename TIndex, typename TSpec>
inline typename SAValue<TIndex>::Type
endPosition(Pattern<TIndex, Swift<TSpec> > const & pattern)
{
	__int64 hitEnd = pattern.curEndPos;
	__int64 textLength = sequenceLength(pattern.curSeqNo, needle(pattern));
	if(hitEnd > textLength) hitEnd = textLength;

	typename SAValue<TIndex >::Type pos;
	posLocalToX(pos, Pair<unsigned, __int64>(pattern.curSeqNo, hitEnd), stringSetLimits(host(pattern)));
	return pos;
}

template <typename TIndex, typename TSpec>
inline typename SAValue<TIndex >::Type
endPosition(Pattern<TIndex, Swift<Tag<_SwiftSemiGlobal<TSpec> > > > const & pattern)
{
	typedef typename Size<TIndex>::Type TSize;
	typename SAValue<TIndex >::Type pos;
	posLocalToX(pos, Pair<unsigned, TSize>(pattern.curSeqNo, length(needle(pattern))), stringSetLimits(host(pattern)));
	return pos;
}

template <typename TIndex, typename TSpec>
inline typename SAValue<TIndex>::Type
endPosition(Pattern<TIndex, Swift<TSpec> > & pattern)
{
	return endPosition(const_cast<Pattern<TIndex, Swift<TSpec> > const &>(pattern));
}

//____________________________________________________________________________
/**
.Function.positionRangeNoClip
..cat:Searching
..summary:Returns a pair of the begin and end position in or beyond the haystack or needle for the last hit found.
..signature:positionRangeNoClip(finder)
..signature:positionRangeNoClip(pattern)
..param.finder:A @Class.Finder@ object.
..param.pattern:A @Class.Pattern@ object.
..returns:A pair of the begin and end position in the haystack or needle for the last hit found. These positions could
be negative or beyond the end of $finder$ or $pattern$ when using filter algorithms.
...remarks:The return type is $Pair<typename SAValue<THost>::Type>$ if $THost$ is the type of haystack or needle.
..see:Function.positionRange
*/
///.Function.positionRangeNoClip.param.finder.type:Spec.Swift
///.Function.positionRangeNoClip.param.pattern.type:Spec.Swift

template <typename THaystack, typename TSpec>
inline Pair<typename Position<Finder<THaystack, Swift<TSpec> > >::Type>
positionRangeNoClip(Finder<THaystack, Swift<TSpec> > const & finder)
{
	typedef typename Position<Finder<THaystack, Swift<TSpec> > >::Type TPosition;
	typedef Pair<TPosition> TPair;
	typename Finder<THaystack, Swift<TSpec> >::TSwiftHit &hit = *finder.curHit;
	return TPair((TPosition)hit.hstkPos, (TPosition)(hit.hstkPos + hit.bucketWidth));
}

template <typename THaystack, typename TSpec>
inline Pair<typename Position<Finder<THaystack, Swift<TSpec> > >::Type>
positionRangeNoClip(Finder<THaystack, Swift<TSpec> > & finder)
{
	return positionRangeNoClip(const_cast<Finder<THaystack, Swift<TSpec> > const &>(finder));
}

//____________________________________________________________________________
/**
.Function.positionRange
..cat:Searching
..summary:Returns a pair of the begin and end position in the haystack or needle for the last hit found.
..signature:positionRange(finder)
..signature:positionRange(pattern)
..param.finder:A @Class.Finder@ object.
..param.pattern:A @Class.Pattern@ object.
..returns:A pair of the begin and end position in the haystack or needle for the last hit found.
...remarks:The return type is $Pair<typename SAValue<THost>::Type>$ if $THost$ is the type of haystack or needle.
..see:Function.beginPosition
..see:Function.endPosition
*/
///.Function.positionRange.param.finder.type:Spec.Swift
///.Function.positionRange.param.pattern.type:Spec.Swift

template <typename THaystack, typename TSpec>
inline Pair<typename Position<Finder<THaystack, Swift<TSpec> > >::Type>
positionRange(Finder<THaystack, Swift<TSpec> > const & finder)
{
	typedef typename Position<Finder<THaystack, Swift<TSpec> > >::Type TPosition;
	typedef Pair<TPosition> TPair;
	typename Finder<THaystack, Swift<TSpec> >::TSwiftHit &hit = *finder.curHit;

	__int64 hitBegin = hit.hstkPos;
	__int64 hitEnd = hit.hstkPos + hit.bucketWidth;
	__int64 textEnd = length(haystack(finder));

	if (hitBegin < 0) hitBegin = 0;
	if (hitEnd > textEnd) hitEnd = textEnd;
	return TPair((TPosition)hitBegin, (TPosition)hitEnd);
}

template <typename THaystack, typename TSpec>
inline Pair<typename Position<Finder<THaystack, Swift<TSpec> > >::Type>
positionRange(Finder<THaystack, Swift<TSpec> > & finder)
{
	return positionRange(const_cast<Finder<THaystack, Swift<TSpec> > const &>(finder));
}

//____________________________________________________________________________

template <typename TIndex, typename TSpec>
inline Pair<typename SAValue<TIndex>::Type>
positionRange(Pattern<TIndex, Swift<TSpec> > & pattern)
{
	return Pair<typename SAValue<TIndex>::Type> (beginPosition(pattern), endPosition(pattern));
}

//____________________________________________________________________________

template <typename TSwiftHit, typename TText>
inline typename Infix<TText>::Type
swiftInfixNoClip(TSwiftHit const &hit, TText &text)
{
	return infix(text, hit.hstkPos, hit.hstkPos + hit.bucketWidth);
}

template <typename TSwiftHit, typename TText>
inline typename Infix<TText>::Type
swiftInfix(TSwiftHit const &hit, TText &text)
{
	__int64 hitBegin = hit.hstkPos;
	__int64 hitEnd = hit.hstkPos + hit.bucketWidth;
	__int64 textEnd = length(text);

	if (hitBegin < 0) hitBegin = 0;
	if (hitEnd > textEnd) hitEnd = textEnd;
	return infix(text, hitBegin, hitEnd);
}

//____________________________________________________________________________

///.Function.infix.remarks:For finders or patterns of filtering algorithms (e.g. @Spec.Swift@) the returned infix is a potential match.
///.Function.infix.param.finder.type:Spec.Swift
///.Function.infix.param.pattern.type:Spec.Swift

template <typename THaystack, typename TSpec>
inline typename Infix<THaystack>::Type
infix(Finder<THaystack, Swift<TSpec> > &finder)
{
	return swiftInfix(*finder.curHit, haystack(finder));
}

template <typename THaystack, typename TSpec, typename TText>
inline typename Infix<TText>::Type
infix(Finder<THaystack, Swift<TSpec> > &finder, TText &text)
{
	return swiftInfix(*finder.curHit, text);
}

//____________________________________________________________________________

template <typename THaystack, typename TSpec>
inline typename Infix<THaystack>::Type
infixNoClip(Finder<THaystack, Swift<TSpec> > &finder)
{
	return swiftInfixNoClip(*finder.curHit, haystack(finder));
}

template <typename THaystack, typename TSpec, typename TText>
inline typename Infix<TText>::Type
infixNoClip(Finder<THaystack, Swift<TSpec> > &finder, TText &text)
{
	return swiftInfixNoClip(*finder.curHit, text);
}

//____________________________________________________________________________

template <typename TIndex, typename TSpec, typename TText>
inline typename Infix<TText>::Type
infix(Pattern<TIndex, Swift<TSpec> > const & pattern, TText &text)
{
    __int64 hitBegin = pattern.curBeginPos;
	__int64 hitEnd = pattern.curEndPos;
	__int64 textLength = sequenceLength(pattern.curSeqNo, needle(pattern));

	if (hitEnd > textLength) hitEnd = textLength;
    if (hitBegin < 0) hitBegin = 0;

	return infix(text, hitBegin, hitEnd);
}

template <typename TIndex, typename TSpec>
inline typename Infix< typename GetSequenceByNo< TIndex const >::Type >::Type
infix(Pattern<TIndex, Swift<TSpec> > const & pattern)
{
	return infix(pattern, getSequenceByNo(pattern.curSeqNo, needle(pattern)));
}

template <typename TIndex, typename TSpec>
inline typename Infix< typename GetSequenceByNo< TIndex const >::Type >::Type
infix(Pattern<TIndex, Swift<Tag<_SwiftSemiGlobal<TSpec> > > > const & pattern)
{
	return infix(getSequenceByNo(pattern.curSeqNo, needle(pattern)), 0, sequenceLength(pattern.curSeqNo, needle(pattern)));
}

template <typename TIndex, typename TSpec>
inline typename Infix< typename GetSequenceByNo< TIndex const >::Type >::Type
infix(Pattern<TIndex, Swift<TSpec> > & pattern)
{
	return infix(const_cast<Pattern<TIndex, Swift<TSpec> > const &>(pattern));
}

//____________________________________________________________________________

template <typename THaystack, typename TSpec>
inline void 
_printDots(Finder<THaystack, Swift<TSpec> > &finder)
{
	while (finder.curPos >= finder.dotPos) 
	{
		finder.dotPos += 100000;
		if (finder.dotPos >= finder.dotPos2)
		{
			::std::cerr << (finder.dotPos2 / 1000000) << "M" << ::std::flush;
			finder.dotPos2 += 1000000;
		} else
			::std::cerr << "." << ::std::flush;
	}
}

template <typename TFinder, typename TIndex, typename TSpec>
inline bool 
_nextNonRepeatRange(
	TFinder &finder,
	Pattern<TIndex, Swift<TSpec> > &pattern)
{
	typedef typename TFinder::TRepeat		TRepeat;
	typedef typename Value<TRepeat>::Type	TPos;

	if (finder.curRepeat == finder.endRepeat) return false;

	do 
	{
		finder.startPos = (*finder.curRepeat).endPosition;
		if (++finder.curRepeat == finder.endRepeat) 
		{
			finder.endPos = length(host(finder));
			if (finder.startPos + length(pattern.shape) > finder.endPos)
				return false;
			else
				break;
		} else
			finder.endPos = (*finder.curRepeat).beginPosition;
		// repeat until the shape fits in non-repeat range
	} while (finder.startPos + length(pattern.shape) > finder.endPos);

	finder.curPos = finder.startPos;
	hostIterator(finder) = begin(host(finder)) + finder.startPos;
	finder.haystackEnd = begin(host(finder)) + (finder.endPos - length(pattern.shape) + 1);

//	if (pattern.params.printDots)
//		::std::cerr << ::std::endl << "  scan range (" << finder.startPos << ", " << finder.endPos << ") " << std::flush;

	return true;
}

template <typename TFinder, typename TIndex, typename TSpec>
inline bool 
_firstNonRepeatRange(
	TFinder &finder,
	Pattern<TIndex, Swift<TSpec> > &pattern)
{
	typedef typename TFinder::TRepeat		TRepeat;
	typedef typename Value<TRepeat>::Type	TPos;

	finder.curRepeat = begin(finder.data_repeats, Standard());
	finder.endRepeat = end(finder.data_repeats, Standard());

	if (finder.curRepeat == finder.endRepeat)
		finder.endPos = length(host(finder));
	else
		finder.endPos = (*finder.curRepeat).beginPosition;

	if (length(pattern.shape) > finder.endPos)
		return _nextNonRepeatRange(finder, pattern);

	finder.curPos = finder.startPos = 0;
	hostIterator(finder) = begin(host(finder));
	finder.haystackEnd = begin(host(finder)) + (finder.endPos - length(pattern.shape) + 1);

//	if (pattern.params.printDots)
//		::std::cerr << ::std::endl << "  scan range (" << finder.startPos << ", " << finder.endPos << ") " << std::flush;

	return true;
}

template <typename TFinder, typename TIndex, typename TSpec>
inline void
_copySwiftHit(
	TFinder &finder,
	Pattern<TIndex, Swift<TSpec> > &pattern)
{
	pattern.curSeqNo = (*finder.curHit).ndlSeqNo;
	pattern.curBeginPos = (*finder.curHit).ndlPos;
	pattern.curEndPos = (*finder.curHit).ndlPos + (*finder.curHit).hitLengthNeedle;
}

template <typename TFinder, typename TIndex, typename TSpec>
inline void 
_copySwiftHit(
	TFinder &finder,
	Pattern<TIndex, Swift<Tag<_SwiftSemiGlobal<TSpec> > > > &pattern)
{
	pattern.curSeqNo = (*finder.curHit).ndlSeqNo;
	pattern.curBeginPos = 0;
	pattern.curEndPos = length(indexText(needle(pattern))[pattern.curSeqNo]);
}

template <typename TFinder, typename TIndex, typename TSpec>
inline bool 
find(
	TFinder &finder,
	Pattern<TIndex, Swift<Tag<_SwiftSemiGlobal<TSpec> > > > &pattern, 
	double errorRate)
{
	return find(finder, pattern, errorRate, 0);
}

template <typename THaystack, typename TIndex, typename TSpec, typename TSize>
inline bool 
find(
	Finder<THaystack, Swift<TSpec> > &finder,
	Pattern<TIndex, Swift<TSpec> > &pattern, 
	double errorRate,
	TSize minLength)
{
	typedef typename Fibre<TIndex, QGram_Shape>::Type	TShape;
	typedef	typename Value<TShape>::Type				THashValue;

	if (empty(finder)) 
	{
		pattern.finderLength = pattern.params.tabooLength + length(container(finder));
		_patternInit(pattern, errorRate, minLength);
		_finderSetNonEmpty(finder);
		finder.dotPos = 100000;
		finder.dotPos2 = 10 * finder.dotPos;

		if (!_firstNonRepeatRange(finder, pattern)) return false;
		if (_swiftMultiProcessQGram(finder, pattern, hash(pattern.shape, hostIterator(hostIterator(finder)))))
		{
			_copySwiftHit(finder, pattern);
			return true;
		}
	} 
	else
	{
		if (++finder.curHit < finder.endHit) 
		{
			_copySwiftHit(finder, pattern);
			return true;
		}
	}

	// all previous matches reported -> search new ones
	clear(finder.hits);

	// are we at the end of the text?
	if (atEnd(finder) && finder.curRepeat == finder.endRepeat) 
	{
		finder.curHit = finder.endHit;
		return false;
	}

	do 
	{
		if (pattern.params.printDots) _printDots(finder);
		if (atEnd(++finder)) 
		{
			if (!_nextNonRepeatRange(finder, pattern)) 
			{
				if(_swiftMultiFlushBuckets(finder, pattern))
				{
					_copySwiftHit(finder, pattern);
					return true;
				}
				else
					return false;
			}
			hash(pattern.shape, hostIterator(hostIterator(finder)));
		}
		else
		{
			++finder.curPos;
			hashNext(pattern.shape, hostIterator(hostIterator(finder)));
		}
		
		if (_swiftMultiProcessQGram(finder, pattern, value(pattern.shape)))
		{
			_copySwiftHit(finder, pattern);
			return true;
		}

	} while (true);
}

template <typename THashes, typename TPipeSpec, typename TIndex, typename TSpec>
inline bool 
find(
	Finder<Pipe<THashes, TPipeSpec>, Swift<TSpec> > &finder,
	Pattern<TIndex, Swift<TSpec> > &pattern, 
	double errorRate)
{
	if (empty(finder)) 
	{
		pattern.finderLength = 0;
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
			_copySwiftHit(finder, pattern);
			return true;
		}
	} else
		if (++finder.curHit != finder.endHit) 
		{
			_copySwiftHit(finder, pattern);
			return true;
		}

	clear(finder.hits);
	if (eof(finder.in)) return false;

	do 
	{
		++finder.in;
		if (eof(finder.in)) 
		{
			endRead(finder.in);
#ifdef SEQAN_DEBUG_SWIFT
			_printSwiftBuckets(pattern);
#endif
			if(_swiftMultiFlushBuckets(finder, pattern))
			{
				_copySwiftHit(finder, pattern);
				return true;
			}
			else 
				return false;
		}
		finder.curPos = (*finder.in).i1;
		if (pattern.params.printDots) _printDots(finder);

	} while (!_swiftMultiProcessQGram(finder, pattern, hash(pattern.shape, (*finder.in).i2)));

	_copySwiftHit(finder, pattern);
	return true;
}

/**
.Function.windowFindBegin:
..cat:Searching
..summary:Initializes the pattern. Sets the finder on the begin position.
 Gets the first non-repeat range and sets it in the finder.
 Used together with @Function.windowFindBegin@ and @Function.windowFindEnd@.
..signature:windowFindBegin(finder, pattern, errorRate)
..param.finder:A SWIFT finder.
..param.pattern: A SWIFT pattern.
..param.errorRate:Error rate that is allowed between reads and reference.
 Schould be the same in as in @Function.windowFindNext@.
...type:Class.Double
*/
template <typename THaystack, typename TIndex, typename TSpec>
inline bool 
windowFindBegin(
	Finder<THaystack, Swift<TSpec> > &finder,
	Pattern<TIndex, Swift<TSpec> > &pattern, 
	double errorRate)
{
	SEQAN_CHECKPOINT
	
	pattern.finderLength = pattern.params.tabooLength + length(container(finder));
	_patternInit(pattern, errorRate, 0);
	_finderSetNonEmpty(finder);
	finder.dotPos = 100000;
	finder.dotPos2 = 10 * finder.dotPos;

	if (!_firstNonRepeatRange(finder, pattern)) return false;
    
    return true;
}


/**
.Function.windowFindNext:
..cat:Searching
..summary:Searches over the next window with the finder. The found hits can be retrieved with @Function.getSwiftHits@
 Used together with @Function.windowFindBegin@ and @Function.windowFindEnd@.
..signature:windowFindNext(finder, pattern, finderWindowLength)
..param.finder:A SWIFT finder.
..param.pattern: A SWIFT pattern.
..param.finderWindowLength:Number of bases that are scanned beginning from the position the finder is at.
 Including bases that are marked as repeats and that are skipped.
...type:nolink:unsigned int
..returns:true, if there are bases that can be scanned. false, otherwise
..see:Function.windowFindBegin
..see:Function.windowFindEnd
..see:Function.getSwiftHits
*/
template <typename THaystack, typename TIndex, typename TSpec, typename TSize>
inline bool 
windowFindNext(
	Finder<THaystack, Swift<TSpec> > &finder,
	Pattern<TIndex, Swift<TSpec> > &pattern, 
	TSize finderWindowLength
#ifdef RAZERS_TIMER
    ,
    Pair<int, int>& posLength	// weese: ifdefs in signatures are eval, and should be removed
#endif
               )
{
	SEQAN_CHECKPOINT
	
	typedef typename Fibre<TIndex, QGram_Shape>::Type	TShape;
	typedef	typename Value<TShape>::Type				THashValue;
	
	typedef Finder<THaystack, Swift<TSpec> >			TFinder;
	typedef typename TFinder::THstkPos					THstkPos;
    
	// all previous matches reported -> search new ones
	clear(finder.hits);

#ifdef RAZERS_TIMER
    posLength.i1 = finder.curPos;
#endif
    
	THstkPos windowEnd = finder.curPos + finderWindowLength;
	// iterate over all non-repeat regions within the window
	for (; finder.curPos < windowEnd; )
	{
		THstkPos localEnd = finder.endPos;
		if (localEnd > windowEnd) localEnd = windowEnd;

#ifdef RAZERS_TIMER
        posLength.i2 += localEnd - finder.curPos;
#endif
        
		// filter a non-repeat region within the window
		TShape &shape = pattern.shape;
		_swiftMultiProcessQGram(finder, pattern, hash(shape, hostIterator(hostIterator(finder))));
        
		for (++finder.curPos, ++finder; finder.curPos < localEnd; ++finder.curPos, ++finder){
			_swiftMultiProcessQGram(finder, pattern, hashNext(shape, hostIterator(hostIterator(finder))));			
        }
            
		if (pattern.params.printDots) _printDots(finder);

		if (finder.curPos == finder.endPos)
			if (!_nextNonRepeatRange(finder, pattern))
                return false;
	}
	return true;
}

/**
.Function.windowFindEnd:
..cat:Searching
..summary:Flushes the pattern. Used together with @Function.windowFindBegin@ and @Function.windowFindNext@.
..signature:windowFindNext(finder, pattern)
..param.finder:A SWIFT finder.
..param.pattern: A SWIFT pattern.
..see:Function.windowFindBegin
*/
template <typename THaystack, typename TIndex, typename TSpec>
inline void 
windowFindEnd(
	Finder<THaystack, Swift<TSpec> > & finder,
	Pattern<TIndex, Swift<TSpec> > &pattern)
{
	SEQAN_CHECKPOINT
	
	_swiftMultiFlushBuckets(finder, pattern);
}


/**
.Function.getSwiftHits:
..cat:Searching
..summary:Gets the string of hits from the finder
..signature:getSwiftHits(finder)
..param.finder:A SWIFT finder.
..returns:@Class.String@ of Hits (use Finder<...>::THitString as Type).
*/
template <typename THaystack, typename TSpec>
inline typename Finder<THaystack, Swift<TSpec> >::THitString &
getSwiftHits(Finder<THaystack, Swift<TSpec> > &finder)
{
	SEQAN_CHECKPOINT
	
	return finder.hits;
}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_FIND_SHIFTAND_H
