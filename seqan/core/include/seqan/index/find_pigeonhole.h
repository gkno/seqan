// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2010, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: David Weese <david.weese@fu-berlin.de>
// ==========================================================================
// This filter uses the pigeonhole principle to search k-error matches 
// between strings in a stringset and a text.
// The pigeonhole principle states:
// If a string is split into k+1 pieces then in every k-error match of the
// whole string one piece matches without error.
// ==========================================================================


#ifndef CORE_INCLUDE_SEQAN_FIND_PIGEONHOLE_H_
#define CORE_INCLUDE_SEQAN_FIND_PIGEONHOLE_H_

namespace seqan 
{

// ==========================================================================
// Forwards
// ==========================================================================

    template < typename TObject, typename TSpec > class Index;
    template < typename TObject > struct SAValue;

// ==========================================================================
// Tags, Classes, Enums
// ==========================================================================

    template <typename TSpec = void>
    struct Pigeonhole {
        enum { ONE_PER_DIAGONAL = 1 };	// 1..record last seed diagonal and ignore seeds on the same diag. (heuristic)
        enum { HAMMING_ONLY = 0 };		// 0..no indels; 1..allow indels
    };

    template <>
    struct Pigeonhole<Hamming_> {
        enum { ONE_PER_DIAGONAL = 1 };	// 1..record last seed diagonal and ignore seeds on the same diag. (heuristic)
        enum { HAMMING_ONLY = 1 };		// 0..no indels; 1..allow indels
    };
	
//////////////////////////////////////////////////////////////////////////////


//____________________________________________________________________________


	template <typename TSpec, typename THstkPos>
	struct PigeonholeHit_ 
	{
		THstkPos	hstkPos;			// seed diagonal begin in haystack 
		unsigned	ndlSeqNo;			// needle sequence number
	};
    
    template <typename THaystack, typename TPattern, typename TPigeonholeSpec>
    struct FindResult<Finder<THaystack, Pigeonhole<TPigeonholeSpec> >, TPattern>
    {
        typedef SwiftHitSemiGlobal_<__int64> Type;
    };
        
    struct PigeonholeParameters {
        unsigned overlap;           // overlap length of adjacent q-grams, default is to use non-overlapping q-grams (=0)
        bool printDots;
        bool debug;
        
        PigeonholeParameters():
            overlap(0),
            printDots(false),		// print a . for every 100kbps mapped genome
            debug(false) {}
    };
    

//____________________________________________________________________________


/**
.Spec.Pigeonhole:
..summary:Provides a fast filter alogrithm that uses the pigeonhole lemma, i.e. if a pattern matches with k errors in the text, every partition into k+1 parts contains one part that matches without error.
..general:Class.Pattern
..general:Class.Finder
..cat:Searching
..signature:Finder<THaystack, Pigeonhole<TSpec> >
..signature:Pattern<TIndex, Pigeonhole<TSpec> >
..param.THaystack:The type of the sequence that should be searched.
..param.TIndex: A q-gram index of needle(s) that should be searched for.
...type:Spec.IndexQGram
..param.TSpec: Specifies the type of Swift filter.
..include:seqan/index.h
..remarks:
 The @Class.Pattern@ must be a q-gram index over multiple patterns. The tolerated error rate must be given when @Function.find@ or @Function.windowFindBegin@ is called.
 In these functions the length of the index @Class.Shape@ is set automatically thus it must be modifiable at runtime, e.g. @Spec.OneGappedShape@.

*/
///.Class.Pattern.param.TSpec.type:Spec.Pigeonhole
///.Class.Finder.param.TSpec.type:Spec.Pigeonhole

	template <typename THaystack, typename TSpec>
	class Finder<THaystack, Pigeonhole<TSpec> >
	{
	public:
		typedef typename Iterator<THaystack, Rooted>::Type			TIterator;
		typedef typename Position<THaystack>::Type					THstkPos;
//		typedef PigeonholeHit_<TSpec, __int64>						TPigeonholeHit;
		typedef typename FindResult<Finder>::Type                   TPigeonholeHit;
		typedef typename WindowFindResult<Finder>::Type				THitString;
		typedef typename Iterator<THitString, Standard>::Type		THitIterator;
		typedef typename SAValue<THaystack>::Type					TSAValue;
		typedef Repeat<TSAValue, unsigned>							TRepeat;
        // TODO(holtgrew): Make this a holder so we can share the actual string when copying.
		typedef String<TRepeat>										TRepeatString;
		typedef typename Iterator<TRepeatString, Standard>::Type	TRepeatIterator;

		TIterator		data_iterator;
		TIterator		haystackEnd;
		bool			_needReinit;	// if true, the Pattern needs to be reinitialized
		THitString		hits;
		THitIterator	curHit, endHit;
		THstkPos		startPos, curPos, endPos;
        THstkPos        windowStart;
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
            windowStart(orig.windowStart),
            dotPos(orig.dotPos),
            dotPos2(orig.dotPos2),
			data_repeats(orig.data_repeats)
		{
			curHit = begin(hits, Standard()) + (orig.curHit - begin(orig.hits, Standard()));
			endHit = end(hits, Standard());
			curRepeat = begin(data_repeats, Standard()) + (orig.curRepeat - begin(orig.data_repeats, Standard()));
			endRepeat = end(data_repeats, Standard());
		}

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
            windowStart = orig.windowStart;
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

	
	template <typename THaystack, typename TSpec>
	inline bool
	atEnd(Finder<THaystack, Pigeonhole<TSpec> > & me)
	{
		return hostIterator(hostIterator(me)) == hostIterator(me.haystackEnd);
	}

	template <typename THaystack, typename TSpec>
	inline void
	goEnd(Finder<THaystack, Pigeonhole<TSpec> > & me)
	{
		hostIterator(me) = me.haystackEnd;
	}


//____________________________________________________________________________


	template <typename TIndex, typename TSpec>
	class Pattern<TIndex, Pigeonhole<TSpec> >
	{
	public:
		typedef typename Size<TIndex>::Type								TSize;
		typedef unsigned												TShortSize;
		typedef typename Fibre<TIndex, Tag<FibreSA_> const >::Type		TSA;
		typedef typename Fibre<TIndex, Tag<Fibre_Shape_> const >::Type	TShape;
		typedef typename Iterator<TSA const, Standard>::Type			TIterator;
				
		TShape					shape;
        __int64                 curBeginPos, curEndPos;
		__int64					finderLength;
		__int64					seqDisabled;
		TSize					finderPosOffset;						// these must be of
		TSize					finderPosNextOffset;					// type TSize of TBucket
        TSize                   maxSeqLen;
		unsigned				curSeqNo;
		String<__int64>         lastSeedDiag;
		PigeonholeParameters    params;
		double					_currentErrorRate;

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


template <typename TValue, typename TSpec>
inline unsigned
_pigeonholeMaxShapeWeight(Shape<TValue, TSpec> const &)
{
    typedef typename Value<Shape<TValue, SimpleShape> >::Type THashValue;
    return (unsigned)(BitsPerValue<THashValue>::VALUE - 1)/ (unsigned)BitsPerValue<TValue>::VALUE;
}

template <typename TShape, typename TSize>
inline bool
_pigeonholeUpdateShapeLength(TShape, TSize)
{
    return false;
}

template <typename TValue, typename TSize>
inline bool
_pigeonholeUpdateShapeLength(Shape<TValue, SimpleShape> &shape, TSize newLength)
{
    TSize weight = _pigeonholeMaxShapeWeight(shape);

    if (weight > newLength) 
		weight = newLength;
		
    if (weight == length(shape))
        return false;

    resize(shape, weight);
    return true;
}

template <typename TValue, typename TSize>
inline bool
_pigeonholeUpdateShapeLength(Shape<TValue, OneGappedShape> &shape, TSize newLength)
{
    if (length(shape) == newLength)
        return false;
    
    TSize weight = _pigeonholeMaxShapeWeight(shape);

    if (weight > newLength)
		weight = newLength;

	CharString str;
	resize(str, weight / 2, '1');
	resize(str, newLength - (weight + 1) / 2, '0');
	resize(str, newLength, '1');
	stringToShape(shape, str);
    return true;
}

template <typename TValue, typename TSize>
inline bool
_pigeonholeUpdateShapeLength(Shape<TValue, GenericShape> &shape, TSize newLength)
{
    if (length(shape) == newLength)
        return false;
    
    TSize weight = _pigeonholeMaxShapeWeight(shape);
    
    if (weight >= newLength)
        weight = newLength;
    
    CharString str;
    resize(str, weight / 2, '1');
    resize(str, newLength - (weight + 1) / 2, '0');
    resize(str, newLength, '1');
    stringToShape(shape, str);

    return true;
}

template <typename TIndex, typename TSpec, typename TSeqNo>
inline void
maskPatternSequence(Pattern<TIndex, Pigeonhole<TSpec> > & pattern, TSeqNo seqNo, bool enable) 
{
    SEQAN_ASSERT_NEQ_MSG((int)Pigeonhole<TSpec>::ONE_PER_DIAGONAL, 0, "Disabling sequences in the pigeonhole filter requires the ONE_PER_DIAGONAL heuristic.");
	
	if (enable)
		pattern.lastSeedDiag[seqNo] = -pattern.maxSeqLen;
	else
		pattern.lastSeedDiag[seqNo] = pattern.seqDisabled;
}

template <typename TIndex, typename TFloat, typename TSpec>
inline void _patternInit(Pattern<TIndex, Pigeonhole<TSpec> > &pattern, TFloat errorRate) 
{
	typedef typename Size<TIndex>::Type TSize;

	double _newErrorRate = errorRate;
	TSize seqCount = countSequences(host(pattern));

	if (pattern._currentErrorRate != _newErrorRate)
	{
		// settings have been changed -> initialize bucket parameters
	
		pattern._currentErrorRate = _newErrorRate;
		pattern.finderPosOffset = 0;
		pattern.finderPosNextOffset = pattern.finderLength;

        TSize minQ = MaxValue<TSize>::VALUE;
        TSize maxQ = 3;
        TSize maxSeqLen = 0;
        for(unsigned seqNo = 0; seqNo < seqCount; ++seqNo) 
        {
            // get pattern length and max. allowed errors
            TSize length = sequenceLength(seqNo, host(pattern));
            if (maxSeqLen < length) maxSeqLen = length;
			
			// sequence must have sufficient length
			if (length <= pattern.params.overlap) continue;
			
			// cut overlap many characters from the end
			length -= pattern.params.overlap;
            TSize errors = (TSize) floor(errorRate * length);
            TSize q = length / (errors + 1);
			
			
			// ignore too short q-grams
			if (q < 3) continue;
            if (minQ > q) minQ = q;
            if (maxQ < q) maxQ = q;
        }
        pattern.maxSeqLen = maxSeqLen;
        if (minQ < 3) minQ = maxQ;
		
        TIndex &index = host(pattern);
		
		if (minQ == MaxValue<TSize>::VALUE)
		{
			// disable index
			minQ = pattern.maxSeqLen + 1;
		}
		
        if (_pigeonholeUpdateShapeLength(pattern.shape, minQ + pattern.params.overlap) || getStepSize(index) != minQ)
        {
            clear(index);
            setStepSize(index, minQ);
            CharString str;
            shapeToString(str, pattern.shape);
#pragma omp critical
			{
				std::cout << std::endl << "Pigeonhole settings:" << std::endl;
				std::cout << "  shape:    " << length(str) << '\t' << str << std::endl;
				std::cout << "  stepsize: " << getStepSize(index) << std::endl;
			}
         }
        indexShape(host(pattern)) = pattern.shape;
        indexRequire(host(pattern), QGramSADir());

        clear(pattern.lastSeedDiag);
        if (Pigeonhole<TSpec>::ONE_PER_DIAGONAL)
            resize(pattern.lastSeedDiag, seqCount, -maxSeqLen);
		pattern.seqDisabled = -maxSeqLen - 1;
    }
}

template <
	typename TFinder,
	typename TIndex, 
	typename TSpec,
	typename THValue
>
inline bool _pigeonholeProcessQGram(
	TFinder &finder, 
	Pattern<TIndex, Pigeonhole<TSpec> > &pattern,
	THValue hash)
{
	typedef Pattern<TIndex, Pigeonhole<TSpec> >         TPattern;
	typedef typename TFinder::THstkPos					THstkPos;
    
    typedef typename Fibre<TIndex, QGramSA>::Type       TSA;
    typedef typename Iterator<TSA, Standard>::Type      TSAIter;
    typedef typename TFinder::TPigeonholeHit            THit;
    
	TIndex const &index = host(pattern);	
    
	// all previous matches reported -> search new ones
	clear(finder.hits);
    
    TSAIter saBegin = begin(indexSA(index), Standard());
    Pair<unsigned> ndlPos;
    THit hit;

    register unsigned bktNo = getBucket(index.bucketMap, hash);
    register TSAIter occ = saBegin + indexDir(index)[bktNo];
    register TSAIter occEnd = saBegin + indexDir(index)[bktNo + 1];
    
    for(; occ != occEnd; ++occ) 
    {
        posLocalize(ndlPos, *occ, stringSetLimits(index));
        hit.hstkPos = finder.curPos - getSeqOffset(ndlPos);		// bucket begin in haystack
        hit.ndlSeqNo = getSeqNo(ndlPos);						// needle seq. number
        if (Pigeonhole<TSpec>::ONE_PER_DIAGONAL)
        {
            __int64 diag = hit.hstkPos + (__int64)pattern.finderPosOffset;
            if (pattern.lastSeedDiag[hit.ndlSeqNo] == diag)
                continue;
            pattern.lastSeedDiag[hit.ndlSeqNo] = diag;
        }
		unsigned ndlLength = sequenceLength(hit.ndlSeqNo, host(pattern));

		if (Pigeonhole<TSpec>::HAMMING_ONLY != 0)
		{
			hit.bucketWidth = ndlLength;
		}
		else
		{
			unsigned indels = (unsigned)floor(pattern._currentErrorRate * ndlLength);
			hit.bucketWidth = ndlLength + (indels << 1);
			hit.hstkPos -= indels;
		}
        appendValue(finder.hits, hit);
    }

	finder.curHit = begin(finder.hits, Standard());
	finder.endHit = end(finder.hits, Standard());

	return !empty(finder.hits);
}


template <typename TIndex, typename TSpec>
inline bool 
empty(Pattern<TIndex, Pigeonhole<TSpec> > & me) 
{
	return me._currentErrorRate == -1;
}

template <typename TIndex, typename TSpec>
inline void 
clear(Pattern<TIndex, Pigeonhole<TSpec> > & me) 
{
	me.finderPosOffset = 0;
	me.finderPosNextOffset = 0;
	me.finderLength = 0;
	me._currentErrorRate = -1;
	clear(me.lastSeedDiag);
}

//____________________________________________________________________________

template <typename THaystack, typename TSpec>
inline typename Position<Finder<THaystack, Pigeonhole<TSpec> > >::Type
position(Finder<THaystack, Pigeonhole<TSpec> > const & finder)
{
	typename Finder<THaystack, Pigeonhole<TSpec> >::TPigeonholeHit &hit = *finder.curHit;
	return hit.hstkPos + hit.bucketWidth;
}

template <typename THaystack, typename TSpec>
inline typename Position<Finder<THaystack, Pigeonhole<TSpec> > >::Type
position(Finder<THaystack, Pigeonhole<TSpec> > & finder)
{
	return position(const_cast<Finder<THaystack, Pigeonhole<TSpec> > const &>(finder));
}

//____________________________________________________________________________

template <typename TIndex, typename TSpec>
inline typename SAValue<TIndex>::Type
position(Pattern<TIndex, Pigeonhole<TSpec> > const & pattern)
{
	typedef typename Size<TIndex>::Type TSize;
    typename SAValue<TIndex >::Type pos;
	posLocalToX(pos, Pair<unsigned, TSize>(pattern.curSeqNo, length(needle(pattern))), stringSetLimits(host(pattern)));
	return pos;
}

template <typename TIndex, typename TSpec>
inline typename SAValue<TIndex>::Type
position(Pattern<TIndex, Pigeonhole<TSpec> > & pattern)
{
	return position(const_cast<Pattern<TIndex, Pigeonhole<TSpec> > const &>(pattern));
}

//____________________________________________________________________________

template <typename THaystack, typename TSpec>
inline __int64
beginPosition(Finder<THaystack, Pigeonhole<TSpec> > const & finder)
{
	return (*finder.curHit).hstkPos;
}

template <typename THaystack, typename TSpec>
inline __int64
beginPosition(Finder<THaystack, Pigeonhole<TSpec> > & finder)
{
	return beginPosition(const_cast<Finder<THaystack, Pigeonhole<TSpec> > const &>(finder));
}

//____________________________________________________________________________

template <typename TIndex, typename TSpec>
inline typename SAValue<TIndex >::Type
beginPosition(Pattern<TIndex, Pigeonhole<TSpec> > const & pattern)
{
    typename SAValue<TIndex >::Type pos;
	posLocalToX(pos, Pair<unsigned>(pattern.curSeqNo, 0), stringSetLimits(host(pattern)));
	return pos;
}

template <typename TIndex, typename TSpec>
inline typename SAValue<TIndex>::Type
beginPosition(Pattern<TIndex, Pigeonhole<TSpec> > & pattern)
{
	return beginPosition(const_cast<Pattern<TIndex, Pigeonhole<TSpec> > const &>(pattern));
}

//____________________________________________________________________________

template <typename THaystack, typename TSpec>
inline typename Position<Finder<THaystack, Pigeonhole<TSpec> > >::Type
endPosition(Finder<THaystack, Pigeonhole<TSpec> > const & finder)
{
	typename Finder<THaystack, Pigeonhole<TSpec> >::TPigeonholeHit &hit = *finder.curHit;
	return hit.hstkPos + hit.bucketWidth;
}

template <typename THaystack, typename TSpec>
inline typename Position<Finder<THaystack, Pigeonhole<TSpec> > >::Type
endPosition(Finder<THaystack, Pigeonhole<TSpec> > & finder)
{
	return endPosition(const_cast<Finder<THaystack, Pigeonhole<TSpec> > const &>(finder));
}

//____________________________________________________________________________

template <typename TIndex, typename TSpec>
inline typename SAValue<TIndex >::Type
endPosition(Pattern<TIndex, Pigeonhole<TSpec> > const & pattern)
{
	typedef typename Size<TIndex>::Type TSize;
	typename SAValue<TIndex >::Type pos;
	posLocalToX(pos, Pair<unsigned, TSize>(pattern.curSeqNo, length(needle(pattern))), stringSetLimits(host(pattern)));
	return pos;
}

template <typename TIndex, typename TSpec>
inline typename SAValue<TIndex>::Type
endPosition(Pattern<TIndex, Pigeonhole<TSpec> > & pattern)
{
	return endPosition(const_cast<Pattern<TIndex, Pigeonhole<TSpec> > const &>(pattern));
}

///.Function.positionRangeNoClip.param.finder.type:Spec.Pigeonhole
///.Function.positionRangeNoClip.param.pattern.type:Spec.Pigeonhole

template <typename THaystack, typename TSpec>
inline Pair<typename Position<Finder<THaystack, Pigeonhole<TSpec> > >::Type>
positionRangeNoClip(Finder<THaystack, Pigeonhole<TSpec> > const & finder)
{
	typedef typename Position<Finder<THaystack, Pigeonhole<TSpec> > >::Type TPosition;
	typedef Pair<TPosition> TPair;
	typename Finder<THaystack, Pigeonhole<TSpec> >::TPigeonholeHit &hit = *finder.curHit;
	return TPair((TPosition)hit.hstkPos, (TPosition)(hit.hstkPos + hit.bucketWidth));
}

template <typename THaystack, typename TSpec>
inline Pair<typename Position<Finder<THaystack, Pigeonhole<TSpec> > >::Type>
positionRangeNoClip(Finder<THaystack, Pigeonhole<TSpec> > & finder)
{
	return positionRangeNoClip(const_cast<Finder<THaystack, Pigeonhole<TSpec> > const &>(finder));
}

///.Function.positionRange.param.finder.type:Spec.Pigeonhole
///.Function.positionRange.param.pattern.type:Spec.Pigeonhole

template <typename THaystack, typename TSpec>
inline Pair<typename Position<Finder<THaystack, Pigeonhole<TSpec> > >::Type>
positionRange(Finder<THaystack, Pigeonhole<TSpec> > const & finder)
{
	typedef typename Position<Finder<THaystack, Pigeonhole<TSpec> > >::Type TPosition;
	typedef Pair<TPosition> TPair;
	typename Finder<THaystack, Pigeonhole<TSpec> >::TPigeonholeHit &hit = *finder.curHit;

	__int64 hitBegin = hit.hstkPos;
	__int64 hitEnd = hit.hstkPos + hit.bucketWidth;
	__int64 textEnd = length(haystack(finder));

	if (hitBegin < 0) hitBegin = 0;
	if (hitEnd > textEnd) hitEnd = textEnd;
	return TPair((TPosition)hitBegin, (TPosition)hitEnd);
}

template <typename THaystack, typename TSpec>
inline Pair<typename Position<Finder<THaystack, Pigeonhole<TSpec> > >::Type>
positionRange(Finder<THaystack, Pigeonhole<TSpec> > & finder)
{
	return positionRange(const_cast<Finder<THaystack, Pigeonhole<TSpec> > const &>(finder));
}

//____________________________________________________________________________

template <typename TIndex, typename TSpec>
inline Pair<typename SAValue<TIndex>::Type>
positionRange(Pattern<TIndex, Pigeonhole<TSpec> > & pattern)
{
	return Pair<typename SAValue<TIndex>::Type> (beginPosition(pattern), endPosition(pattern));
}

//____________________________________________________________________________

template <typename TPigeonholeHit, typename TText>
inline typename Infix<TText>::Type
pigeonholeInfixNoClip(TPigeonholeHit const &hit, TText &text)
{
	return infix(text, hit.hstkPos, hit.hstkPos + hit.bucketWidth);
}

template <typename TPigeonholeHit, typename TText>
inline typename Infix<TText>::Type
pigeonholeInfix(TPigeonholeHit const &hit, TText &text)
{
	__int64 hitBegin = hit.hstkPos;
	__int64 hitEnd = hit.hstkPos + hit.bucketWidth;
	__int64 textEnd = length(text);

	if (hitBegin < 0) hitBegin = 0;
	if (hitEnd > textEnd) hitEnd = textEnd;
	SEQAN_ASSERT_LEQ(hitBegin, hitEnd);
	return infix(text, hitBegin, hitEnd);
}

//____________________________________________________________________________

///.Function.infix.param.finder.type:Spec.Pigeonhole

template <typename THaystack, typename TSpec>
inline typename Infix<THaystack>::Type
infix(Finder<THaystack, Pigeonhole<TSpec> > &finder)
{
	return swiftInfix(*finder.curHit, haystack(finder));
}

template <typename THaystack, typename TSpec, typename TText>
inline typename Infix<TText>::Type
infix(Finder<THaystack, Pigeonhole<TSpec> > &finder, TText &text)
{
	return swiftInfix(*finder.curHit, text);
}

//____________________________________________________________________________

template <typename THaystack, typename TSpec>
inline typename Infix<THaystack>::Type
infixNoClip(Finder<THaystack, Pigeonhole<TSpec> > &finder)
{
	return swiftInfixNoClip(*finder.curHit, haystack(finder));
}

template <typename THaystack, typename TSpec, typename TText>
inline typename Infix<TText>::Type
infixNoClip(Finder<THaystack, Pigeonhole<TSpec> > &finder, TText &text)
{
	return swiftInfixNoClip(*finder.curHit, text);
}

//____________________________________________________________________________

template <typename TIndex, typename TSpec, typename TText>
inline typename Infix<TText>::Type
infix(Pattern<TIndex, Pigeonhole<TSpec> > const & pattern, TText &text)
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
infix(Pattern<TIndex, Pigeonhole<TSpec> > const & pattern)
{
	return infix(getSequenceByNo(pattern.curSeqNo, needle(pattern)), 0, sequenceLength(pattern.curSeqNo, needle(pattern)));
}

template <typename TIndex, typename TSpec>
inline typename Infix< typename GetSequenceByNo< TIndex const >::Type >::Type
infix(Pattern<TIndex, Pigeonhole<TSpec> > & pattern)
{
	return infix(const_cast<Pattern<TIndex, Pigeonhole<TSpec> > const &>(pattern));
}

//____________________________________________________________________________

template <typename THaystack, typename TSpec>
inline void 
_printDots(Finder<THaystack, Pigeonhole<TSpec> > &finder)
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
	Pattern<TIndex, Pigeonhole<TSpec> > &pattern)
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
	Pattern<TIndex, Pigeonhole<TSpec> > &pattern)
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
_copyPigeonholeHit(
	TFinder &finder,
	Pattern<TIndex, Pigeonhole<TSpec> > &pattern)
{
	pattern.curSeqNo = (*finder.curHit).ndlSeqNo;
	pattern.curBeginPos = 0;
	pattern.curEndPos = length(indexText(needle(pattern))[pattern.curSeqNo]);
}

template <typename THaystack, typename TIndex, typename TSpec>
inline bool 
find(
	Finder<THaystack, Pigeonhole<TSpec> > &finder,
	Pattern<TIndex, Pigeonhole<TSpec> > &pattern, 
	double errorRate)
{
	if (empty(finder)) 
	{
		_patternInit(pattern, errorRate);
		pattern.finderLength = pattern.maxSeqLen + length(container(finder));
		_finderSetNonEmpty(finder);
		finder.dotPos = 100000;
		finder.dotPos2 = 10 * finder.dotPos;

		if (!_firstNonRepeatRange(finder, pattern)) return false;
		if (_pigeonholeProcessQGram(finder, pattern, hash(pattern.shape, hostIterator(hostIterator(finder)))))
		{
			_copyPigeonholeHit(finder, pattern);
			return true;
		}
	} 
	else
	{
		if (++finder.curHit < finder.endHit) 
		{
			_copyPigeonholeHit(finder, pattern);
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
                return false;
			hash(pattern.shape, hostIterator(hostIterator(finder)));
		}
		else
		{
			++finder.curPos;
			hashNext(pattern.shape, hostIterator(hostIterator(finder)));
		}
		
		if (_pigeonholeProcessQGram(finder, pattern, value(pattern.shape)))
		{
			_copyPigeonholeHit(finder, pattern);
			return true;
		}

	} while (true);
}

/**
.Function.windowFindBegin:
..param.finder.type:Spec.Pigeonhole
..param.pattern.type:Spec.Pigeonhole
*/

template <typename THaystack, typename TIndex, typename TSpec>
inline bool 
windowFindBegin(
	Finder<THaystack, Pigeonhole<TSpec> > &finder,
	Pattern<TIndex, Pigeonhole<TSpec> > &pattern, 
	double errorRate)
{
	SEQAN_CHECKPOINT
	
	_patternInit(pattern, errorRate);
	pattern.finderLength = pattern.maxSeqLen + length(container(finder));
	_finderSetNonEmpty(finder);
	finder.dotPos = 100000;
	finder.dotPos2 = 10 * finder.dotPos;
    finder.windowStart = 0;

	if (!_firstNonRepeatRange(finder, pattern)) return false;
    
    return true;
}


/**
.Function.windowFindNext:
..param.finder.type:Spec.Pigeonhole
..param.pattern.type:Spec.Pigeonhole
*/

template <typename THaystack, typename TIndex, typename TSpec, typename TSize>
inline bool 
windowFindNext(
	Finder<THaystack, Pigeonhole<TSpec> > &finder,
	Pattern<TIndex, Pigeonhole<TSpec> > &pattern, 
	TSize finderWindowLength)
{
	SEQAN_CHECKPOINT
	
	typedef typename Fibre<TIndex, QGramShape>::Type	TShape;
	typedef Finder<THaystack, Pigeonhole<TSpec> >		TFinder;
	typedef typename TFinder::THstkPos					THstkPos;
    
    typedef typename Fibre<TIndex, QGramSA>::Type       TSA;
    typedef typename Iterator<TSA, Standard>::Type      TSAIter;
    typedef typename TFinder::TPigeonholeHit            THit;

    TIndex const &index = host(pattern);	

	// all previous matches reported -> search new ones
	clear(finder.hits);
    
	THstkPos windowEnd = finder.windowStart + finderWindowLength;
    TSAIter saBegin = begin(indexSA(index), Standard());
    Pair<unsigned> ndlPos;
    THit hit;
	double errorRate = pattern._currentErrorRate;
	__int64 seqDisabled = pattern.seqDisabled;
        
    // iterate over all non-repeat regions within the window
    for (; finder.curPos < windowEnd; )
    {
        THstkPos nonRepeatEnd = finder.endPos - length(pattern.shape) + 1;
        THstkPos localEnd = _min(windowEnd, nonRepeatEnd);
        
        // filter a non-repeat region within the window
        if (finder.curPos < localEnd)
        {
            TShape &shape = pattern.shape;
            
            hashInit(shape, hostIterator(hostIterator(finder)));
            for (; finder.curPos < localEnd; ++finder.curPos, ++finder)
            {
                hashNext(shape, hostIterator(hostIterator(finder)));
                
                register unsigned bktNo = getBucket(index.bucketMap, value(shape));
                register TSAIter occ = saBegin + indexDir(index)[bktNo];
                register TSAIter occEnd = saBegin + indexDir(index)[bktNo + 1];

				for(; occ != occEnd; ++occ) 
				{
					posLocalize(ndlPos, *occ, stringSetLimits(index));
					hit.hstkPos = finder.curPos - getSeqOffset(ndlPos);	// bucket begin in haystack
					hit.ndlSeqNo = getSeqNo(ndlPos);						// needle seq. number

					if (Pigeonhole<TSpec>::ONE_PER_DIAGONAL)
					{
						register __int64 lastDiag = pattern.lastSeedDiag[hit.ndlSeqNo];
						if (lastDiag == seqDisabled) continue;
						
						register __int64 diag = hit.hstkPos + (__int64)pattern.finderPosOffset;
						if (lastDiag == diag)
						{
//							std::cout<<"double hit"<<std::endl;
							continue;
						}
						pattern.lastSeedDiag[hit.ndlSeqNo] = diag;
					}
//                     std::cout<<"hit at:\thpos="<<finder.curPos<<"\tnpos="<<getSeqOffset(ndlPos)<<"\treadNo="<<hit.ndlSeqNo<<"\thash="<<value(shape)<<std::endl;
					unsigned ndlLength = sequenceLength(hit.ndlSeqNo, host(pattern));
					unsigned errors = (unsigned)floor(errorRate * ndlLength);
					hit.bucketWidth = ndlLength + (errors << 1) + 1;
					hit.hstkPos -= errors;
					appendValue(finder.hits, hit);
				}
            }
        }
        
        if (pattern.params.printDots) _printDots(finder);
        
        if (finder.curPos >= nonRepeatEnd)
            if (!_nextNonRepeatRange(finder, pattern))
            {
                finder.windowStart = windowEnd;
                return false;
            }
    }
    finder.windowStart = windowEnd;
	return true;
}

/**
.Function.windowFindEnd:
..param.finder.type:Spec.Pigeonhole
..param.pattern.type:Spec.Pigeonhole
*/

template <typename THaystack, typename TIndex, typename TSpec>
inline void 
windowFindEnd(
	Finder<THaystack, Pigeonhole<TSpec> > &,
	Pattern<TIndex, Pigeonhole<TSpec> > &)
{
}

///.Function.getWindowFindHits.param.finder.type:Spec.Pigeonhole

template <typename THaystack, typename TSpec>
inline typename Finder<THaystack, Pigeonhole<TSpec> >::THitString &
getWindowFindHits(Finder<THaystack, Pigeonhole<TSpec> > &finder)
{
	SEQAN_CHECKPOINT
	
	return finder.hits;
}

///.Function.getMaxDeviationOfOrder.param.pattern.type:Spec.Pigeonhole

template <typename TIndex, typename TSpec>
inline typename Size<TIndex>::Type
getMaxDeviationOfOrder(Pattern<TIndex, Pigeonhole<TSpec> > &pattern)
{
	SEQAN_CHECKPOINT;
   
	return (pattern.maxSeqLen <= length(indexShape(host(pattern))))? 0: pattern.maxSeqLen - length(indexShape(host(pattern)));
}

}// namespace seqan

#endif //#ifndef CORE_INCLUDE_SEQAN_FIND_PIGEONHOLE_H_
