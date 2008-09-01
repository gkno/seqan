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

#ifndef SEQAN_HEADER_RECOGNITION_RATE_DP_H
#define SEQAN_HEADER_RECOGNITION_RATE_DP_H

namespace SEQAN_NAMESPACE_MAIN
{

//____________________________________________________________________________
// Extensions to SeqAn

	struct _ErrorAlphabet {};
	typedef SimpleType<unsigned char,_ErrorAlphabet> ErrorAlphabet;

	template <> struct ValueSize< ErrorAlphabet > { enum { VALUE = 4 }; };
	template <> struct BitsPerValue< ErrorAlphabet > { enum { VALUE = 2 }; };

	template <typename T = void>
	struct _Translate_Table_Error_2_Ascii
	{
		static char const VALUE[4];
	};
	template <typename T>
	char const _Translate_Table_Error_2_Ascii<T>::VALUE[4] = {'.', 'M', 'I', 'D'};

	inline void assign(Ascii & c_target, 
					   ErrorAlphabet const & source)
	{
	SEQAN_CHECKPOINT
		c_target = _Translate_Table_Error_2_Ascii<>::VALUE[source.value];
	}


	struct ErrorPackedString;

	template <typename TValue>
	struct Host<String<TValue, Packed<ErrorPackedString> > >
	{
		typedef String<__int64, Array<1> > Type;
	};

	template <typename TValue>
	struct Host<String<TValue, Packed<ErrorPackedString> > const >
	{
		typedef String<__int64, Array<1> > const Type;
	};


//____________________________________________________________________________

	enum ErrorType {
		SEQAN_MATCH    = 0,
		SEQAN_MISMATCH = 1,
		SEQAN_INSERT   = 2,
		SEQAN_DELETE   = 3
	};

	// descriptor of the modification pattern
	// in the recursion it modifies the last q-gram of a read sequence
	template <typename TFloat>
	struct State
	{
		TFloat prob;				// probability of this state
		int transition[4];			// returns previous state
		unsigned char errors;		// errors in this state
		unsigned char len;			// length of this pattern (shapeSpan-errors <= this value <= shapeSpan+errors)
		bool skipFirst:1;			// skip this pattern if it is the first
		bool skipLast:1;			// skip this pattern if it is the last
		bool intermediate:1;		// this is an intermediate result (beginning with INSERT)
		bool qgramHit:1;			// is this a q-gram hit? (result of the former delta function)
	};


#ifdef PLATFORM_WINDOWS

	template<typename TValue>
	inline bool isnan(TValue value)
	{
		return value != value;
	}

	template<typename TValue>
	inline bool isinf(TValue value)
	{
		return value == log(0.0);
	}

	template<typename TValue>
	inline TValue log1p(TValue value)
	{
		return log(1 + value);
	}

#endif

#ifdef USE_LOGVALUES

	template <typename TValue>
	inline TValue
	_transform(TValue a)
	{
		return log(a);
	}

	template <typename TValue>
	inline TValue
	_transformBack(TValue a)
	{
		return exp(a);
	}

	//////////////////////////////////////////////////////////////////////////////
	// Returns the sum of two probability values in log space
	template <typename TValue>
	inline void
	_probAdd(TValue &a, TValue b)
	{
		if (isinf(a)) {
			a = b;
			return;
		}
		if (isinf(b)) return;
		if (isnan(a + log(1 + exp(b - a)))) return;
		a += log(1 + exp(b - a));
	}

	template <typename TValue>
	inline TValue
	_probMul(TValue a, TValue b)
	{
		return a + b;
	}

	template <typename TValue>
	inline TValue
	_probDiv(TValue a, TValue b)
	{
		return a - b;
	}

#else

	template <typename TValue>
	inline TValue
	_transform(TValue a)
	{
		return a;
	}

	template <typename TValue>
	inline TValue
	_transformBack(TValue a)
	{
		return a;
	}

	template <typename TValue>
	inline void
	_probAdd(TValue &a, TValue b)
	{
		a += b;
	}

	template <typename TValue>
	inline TValue
	_probMul(TValue a, TValue b)
	{
		return a * b;
	}

	template <typename TValue>
	inline TValue
	_probDiv(TValue a, TValue b)
	{
		return a / b;
	}

#endif

struct ErrorPatternLess
{
	template <typename TPattern>
	bool operator() (TPattern const &a, TPattern const &b) const
	{
		typedef typename Iterator<TPattern const>::Type TIter;
		TIter itA = end(a, Standard());
		TIter itB = end(b, Standard());
		TIter itEnd;
		if (length(a) <= length(b))
		{
			itEnd = begin(a, Standard());
			for (; itA != itEnd;) 
			{
				--itA;
				--itB;
				if (*itA < *itB) return true;
				if (*itA > *itB) return false;
			}
			return false;
		} else 
		{
			itEnd = begin(b, Standard());
			for (; itB != itEnd;) 
			{
				--itA;
				--itB;
				if (*itA < *itB) return true;
				if (*itA > *itB) return false;
			}
			return true;
		}
	}
};

template <typename TPatternStore, typename TPattern>
inline int 
_getErrorPatternIndex(TPatternStore const &patternStore, TPattern const &pattern)
{
	typedef typename Iterator<TPatternStore const>::Type TIter;
	TIter lb = lower_bound(begin(patternStore, Standard()), end(patternStore, Standard()), pattern, ErrorPatternLess());
	if (*lb == pattern) {
//		::std::cout << pattern;
		return lb - begin(patternStore, Standard());
	} else {
//		::std::cerr << "  !Pattern Not Found! " << pattern << "\tnext is " << *lb << ::std::endl;
		return -1;
	}
}

// Cut 1 read character and trailing INSERTs of the pattern
template <typename TPattern>
inline int 
_cutErrorPattern(TPattern &_pattern)
{
	typedef typename Iterator<TPattern const, Standard>::Type TIter;
	TPattern const & pattern = const_cast<TPattern const&>(_pattern);
	TIter it = end(pattern, Standard());
	int cuttedErrors = -2;

	// cut trailing INSERTs
	do {
		--it;
		++cuttedErrors;
	} while ((int)getValue(it) == SEQAN_INSERT);

	// cut non INSERT
	if ((int)getValue(it) != SEQAN_MATCH)
		++cuttedErrors;

	//  and all adjacent INSERTs
	do {
		--it;
		++cuttedErrors;
	} while ((int)getValue(it) == SEQAN_INSERT);

	resize(_pattern, 1 + (it - begin(pattern, Standard())));
	return cuttedErrors;
}

template < typename TLogErrorDistr >
typename Value<TLogErrorDistr>::Type 
_getProb(TLogErrorDistr const &logError, int errorType, int readPos)
{
	int maxN = length(logError) / 4;
	SEQAN_ASSERT(readPos >= 0 && readPos < maxN);
	return logError[maxN * (int)errorType + readPos];
}

//////////////////////////////////////////////////////////////////////////////
// Returns log probability of q-gram-configuration q ending at position pos in sequence
template < typename TPattern, typename TLogErrorDistr >
typename Value<TLogErrorDistr>::Type 
_getPatternProb(TLogErrorDistr const &logError, TPattern const &pattern, int readPos)
{
	typename Value<TLogErrorDistr>::Type prob = _transform(1.0);
	for (int i = 0, j = 0; j < (int)length(pattern); ++j)
	{
		prob = _probMul(prob, _getProb(logError, getValue(pattern, j), readPos + i));
		if ((int)getValue(pattern, j) != SEQAN_INSERT)
			++i;
	}
	return prob;
}



//////////////////////////////////////////////////////////////////////////////
// Initialize states-string for edit/hamming-distance filters
template <
	typename TStateString,
	typename TShape,
	typename TLogErrorDistr >
void initPatterns(
	TStateString &states,				// resulting states-string
	TShape const &bitShape,				// bit-string of the shape
	int maxErrors,						// allowed errors per pattern
	TLogErrorDistr const &logError,		// error distribution
	bool optionHammingOnly = false,		// enumerate only hamming patterns (disabled by default)
	bool optionMinOutput = true)		// omit output
{
#ifndef DEBUG_RECOG_DP
//	typedef String<ErrorAlphabet, Packed<ErrorPackedString> >	TPattern;
	typedef String<ErrorAlphabet>								TPattern;
#endif

	typedef typename Iterator<TPattern, Standard>::Type			TIter;
	typedef typename Value<TStateString>::Type					TState;
	
	ErrorType lastErrorType = (optionHammingOnly)? SEQAN_MISMATCH: SEQAN_DELETE;

	SEQAN_ASSERT(SEQAN_MATCH == 0);
	SEQAN_ASSERT((length(logError) % 4) == 0);

#ifndef DEBUG_RECOG_DP
	String<TPattern> patternStore;
#endif

	// a modifier is a pair of position and error type
	String<Pair<int, ErrorType> > mods;
	fill(mods, maxErrors, Pair<int, ErrorType> (0, SEQAN_MATCH));

	TPattern pattern;
	int span = length(bitShape);
	int maxN = length(logError) / 4;

	//////////////////////////////////////////////////////////////////////////////
	// Enumerate all edit-modification patterns with up to k errors
	if (maxErrors == 0) 
	{
		fill(pattern, span, (ErrorAlphabet)SEQAN_MATCH);
		appendValue(patternStore, pattern);
	}
	else
	do 
	{
		clear(pattern);
		fill(pattern, span, (ErrorAlphabet)SEQAN_MATCH);

		// place errors in the pattern
		bool skip = false;
		for (int i = 0; (i < maxErrors) && !skip; ++i)
		{
//			::std::cout << mods[i].i1 << " " << (ErrorAlphabet)mods[i].i2 << "\t";
			switch (mods[i].i2)
			{
			case SEQAN_MISMATCH:
			case SEQAN_DELETE:
				if (pattern[mods[i].i1] != (ErrorAlphabet)SEQAN_MATCH)
				{
					skip = true;
					break;
				}
				pattern[mods[i].i1] = (ErrorAlphabet)mods[i].i2;
				break;

			case SEQAN_INSERT:
				insertValue(pattern, mods[i].i1, (ErrorAlphabet)SEQAN_INSERT);
				break;
				
			case SEQAN_MATCH:
				break;
			}
		}

		// remove redundant patterns
		if (!skip) 
		{
			TIter it = begin(pattern, Standard());
			TIter itEnd = end(pattern, Standard());
			int left = getValue(it);
			int right;
			for (++it; (it != itEnd) && !skip; ++it, left = right) 
			{
				right = getValue(it);
				
				if (left == SEQAN_MISMATCH && right == SEQAN_DELETE) 
					skip = true;	// MISMATCH before DELETE is DELETE before MISMATCH (already enumerated)

				if (left == SEQAN_MISMATCH && right == SEQAN_INSERT) 
					skip = true;	// MISMATCH before INSERT is INSERT before MISMATCH (already enumerated)

				if (left == SEQAN_INSERT && right == SEQAN_DELETE) 
					skip = true;	// INSERT before DELETE is one MISMATCH (already enumerated)
				
				if (left == SEQAN_DELETE && right == SEQAN_INSERT) 
					skip = true;	// DELETE before INSERT is one MISMATCH (already enumerated)
			}
			if (left == SEQAN_INSERT)
				skip = true;		// no trailing INSERT allowed
		}

		if (!skip)
		{
			appendValue(patternStore, pattern);
//			::std::cout << pattern << ::std::endl;
		}

		// reposition modifiers
		int i = 0;
		for (; i < maxErrors; ++i)
		{
			if (mods[i].i2 == SEQAN_MATCH) continue;
			int endPos = (mods[i].i2 == SEQAN_INSERT)? span + 1: span;
			if (++mods[i].i1 < endPos) 
			{
				for(--i; i >= 0; --i)
					mods[i].i1 = mods[i + 1].i1;
				break;
			}
		}

		if (i < maxErrors) continue;

		for (i = 0; i < maxErrors; ++i)
			mods[i].i1 = 0;
		
		// next state combination
		for (i = 0; i < maxErrors; ++i)
		{
			if (mods[i].i2 == lastErrorType) continue;
			mods[i].i2 = (ErrorType)(mods[i].i2 + 1);
			for(--i; i >= 0; --i)
				mods[i].i2 = SEQAN_MISMATCH;
			break;
		}
		
		if (i == maxErrors) break;

	} while (true);
	
	if (!optionMinOutput) 
		::std::cout << "Stored " << length(patternStore) << " modification patterns" << ::std::flush;

	::std::sort(begin(patternStore, Standard()), end(patternStore, Standard()), ErrorPatternLess());
	for (int p = 1; p < (int)length(patternStore); ++p)
	{
		if (patternStore[p-1] == patternStore[p])
			::std::cerr << "  !Found duplicate! " << patternStore[p] << ::std::endl;
	}

	if (!optionMinOutput) 
		::std::cout << " and sorted them." << ::std::endl;

	//////////////////////////////////////////////////////////////////////////////
	// Calculate transitions
	resize(states, length(patternStore));
	for (int p = 0; p < (int)length(patternStore); ++p)
	{
		pattern = patternStore[p];
		TState &state = states[p];

//		::std::cout << pattern << "\t";

		// count errors of current pattern
		int errors = 0;
		for (int i = 0; i < (int)length(pattern); ++i)
			if ((int)getValue(pattern, i) != SEQAN_MATCH)
				++errors;
				
		state.len = length(pattern);
		state.errors = errors;
		state.intermediate = (int)getValue(pattern, 0) == SEQAN_INSERT;
		state.prob = _getPatternProb(logError, pattern, maxN - span);
//		::std::cout << pattern << "\t";

		state.skipFirst = false;
		state.skipLast = false;

		int err = 0, del = 0;
		for (int j = 0; j < (int)length(pattern); ++j)
		{
			switch ((int)getValue(pattern, j)) {
				case SEQAN_MATCH:
					++del;
					break;

				case SEQAN_DELETE:
					++del;
	
				case SEQAN_INSERT:
					++err;
					break;

				default:;
			}
			if (del > 0 && del <= err)
				state.skipFirst = true;
		}
		err = del = 0;
		for (int j = (int)length(pattern) - 1; j >= 0; --j)
		{
			switch ((int)getValue(pattern, j)) {
				case SEQAN_MATCH:
					++del;
					break;

				case SEQAN_DELETE:
					++del;
	
				case SEQAN_INSERT:
					++err;
					break;
			
				default:;
			}
			if (del > 0 && del <= err)
				state.skipLast = true;
		}

		// apply pattern to read q-gram
		// and check if shape is recognized in the genome
		state.qgramHit = false;
		int delta = 0;
		for (int j = 0, readPos = 0, genomePos = 0; j < (int)length(pattern); ++j) 
		{
			switch ((int)getValue(pattern, j))
			{
				case SEQAN_MATCH:
					if (readPos == 0) {
						// assert(bitShape[0] == '1')
						delta = genomePos;
						state.qgramHit = true;
					} else
						if (bitShape[readPos] == '1')
							state.qgramHit &= (readPos + delta == genomePos);
//					::std::cout << readPos;
					++readPos; ++genomePos;
					break;
				case SEQAN_MISMATCH:
					// was it a relevant read position?
					if (bitShape[readPos] == '1')
						state.qgramHit = false;
//					::std::cout << 'x';
					++readPos; ++genomePos;
					break;
				case SEQAN_DELETE:
					// was it a relevant read position?
					if (bitShape[readPos] == '1')
						state.qgramHit = false;
					++readPos;
					break;
				case SEQAN_INSERT:
					++genomePos;
//					::std::cout << 'x';
			}				
		}
//		::std::cout << ::std::endl;

		// prepend INSERT
		++errors;
		insertValue(pattern, 0, SEQAN_INSERT);
		if (errors <= maxErrors)
			state.transition[SEQAN_INSERT] = _getErrorPatternIndex(patternStore, pattern);
		else
			state.transition[SEQAN_INSERT] = -1;

		// prepend MISMATCH and cut INSERTS
		errors -= _cutErrorPattern(pattern);
		pattern[0] = SEQAN_MISMATCH;
		if (errors <= maxErrors)
			state.transition[SEQAN_MISMATCH] = _getErrorPatternIndex(patternStore, pattern);
		else
			state.transition[SEQAN_MISMATCH] = -1;

		// prepend DELETE
		pattern[0] = SEQAN_DELETE;
		if (errors <= maxErrors)
			state.transition[SEQAN_DELETE] = _getErrorPatternIndex(patternStore, pattern);
		else
			state.transition[SEQAN_DELETE] = -1;

		// prepend MATCH
		--errors;
		pattern[0] = SEQAN_MATCH;
		if (errors <= maxErrors)
			state.transition[SEQAN_MATCH] = _getErrorPatternIndex(patternStore, pattern);
		else
			state.transition[SEQAN_MATCH] = -1;
/*		
		::std::cout << "\t" << state.errors;
		::std::cout << "\t" << state.qgramHit;
		::std::cout << "\t" << state.leftError;
		::std::cout << "\t" << state.rightError;
		::std::cout << "\t" << state.transition[0];
		::std::cout << "\t" << state.transition[1];
		::std::cout << "\t" << state.transition[2];
		::std::cout << "\t" << state.transition[3];
		::std::cout << ::std::endl;
*/	}
	if (!optionMinOutput) 
		::std::cout << "Preprocessing finished." << ::std::endl;
}


//////////////////////////////////////////////////////////////////////////////
// Compute filtering loss of any q-gram filter (given a states-string)
template <
	typename TLossString, 
	typename TLogErrorDistr, 
	typename TStateString >
void computeFilteringLoss(
	TLossString &found,
	TStateString const &states,
	int span,
	int maxT,
	int maxErrors,
	TLogErrorDistr const &logError,
	bool optionAbsolute = false,
	bool optionMinOutput = true)
{
	typedef typename Value<TLossString>::Type		TFloat;
	typedef typename Value<TLogErrorDistr>::Type	TProbValue;
	typedef typename Value<TStateString>::Type		TState;

	typedef String<TFloat>							TMatrixCol;
	typedef String<int>								TIntCol;

	SEQAN_ASSERT((length(logError) % 4) == 0);

	int maxN = length(logError) / 4;
	int statesCount = length(states);
//	int span = length(bitShape);

	// columns n-1 and n for recursion 
	TMatrixCol col0;
	TMatrixCol col1;
	fill(col0, maxErrors * statesCount * maxT, _transform(0.0));
	resize(col1, maxErrors * statesCount * maxT);

#ifdef COUNT_LOSSES
	TFloat positive = _transform(0.0);
	TFloat negative = _transform(1.0);
#else
	TFloat positive = _transform(1.0);
	TFloat negative = _transform(0.0);
#endif

	// RECURSION BEGIN
	for (int s = 0; s < statesCount; ++s) 
	{
		TState const &state = states[s];

		if (state.skipFirst) continue;

		// we miss no match if threshold t is 0
		col0[s*maxT] = positive;

		// for n==0
		if (state.qgramHit)
		{
			// we miss no match if read q-gram is recognized
			// --> probability of finding this MMP is 1, if t=1
			col0[s*maxT+1] = positive;
			// --> probability of finding this MMP is 0, if t>1
			for (int t = 2; t < maxT; ++t)
				col0[s*maxT+t] = negative;
		} else
		{
			// we miss 1 match if t>0 and read q-gram is not recognized
			// --> probability of finding this MMP is 0, if t>=1
			for (int t = 1; t < maxT; ++t)
				col0[s*maxT+t] = negative;
		}
	}

	// iterate over sequence length n
	TMatrixCol *col = &col1;
	TMatrixCol *colPrev = &col0;

#ifdef DEBUG_RECOG_DP
	::std::cout << span << ":0";
	dump(col0, 0,statesCount);
	::std::cout << " :1";
	dump(col0, 1,statesCount);
#endif
	

	// RECURSION
	//
	// found(n,q,t,e) = (1-errorProb[n-span]) * found(n-1,0|(q>>1),t-delta,e) delta=1/0 <-> q hat 0/>0 fehler
	//               + errorProb[n-span] * found(n-1,1|(q>>1),t-delta,e-1)
	
	// rekursion (fuer q-gram matches <=1 fehler)
	// found(n,q,t,e) = (1-errorProb[n-span]) * found(n-1,0|(q>>1),t-delta,e) delta=1/0 <-> q hat <=1/>1 fehler
	//               + errorProb[n-span] * found(n-1,1|(q>>1),t-delta,e-1)
	
	for (int n = span; n < maxN; ++n)
	{
		for (int e = 0; e < maxErrors * statesCount; e += statesCount)
		{		
			for (int s = 0; s < statesCount; ++s)
			{
				TState const &state = states[s];				
				for (int t = 0; t < maxT; ++t)
				{
					int _t = t;
					if (_t > 0 && state.qgramHit) --_t;

					// MATCH
					TFloat recovered = _probMul(
						_getProb(logError, SEQAN_MATCH, n-span),
						(*colPrev)[(e+state.transition[SEQAN_MATCH])*maxT+_t]);

					// MISMATCH, INSERT, DELETE
					for (int m = SEQAN_MISMATCH; m < 4; ++m)
						if (e > 0)
						{
							int prevState = state.transition[m];
							if (prevState >= 0)
								if (m == SEQAN_INSERT)
									_probAdd(recovered, _probMul(_getProb(logError,m,n-span), (*col)[((e-statesCount)+prevState)*maxT+t]));
								else
									_probAdd(recovered, _probMul(_getProb(logError,m,n-span), (*colPrev)[((e-statesCount)+prevState)*maxT+_t]));
						}
					(*col)[(e+s)*maxT+t] = recovered;
				}
			}
			if (!optionMinOutput)
				::std::cout << '.' << ::std::flush;
		}

		TMatrixCol *tmp = col;
		col = colPrev;
		colPrev = tmp;

#ifdef DEBUG_RECOG_DP
		::std::cout << n+1 << ":0";
		dump(*colPrev, 0,statesCount);
		::std::cout << " :1";
		dump(*colPrev, 1,statesCount);
		::std::cout << " :2";
		dump(*colPrev, 2,statesCount);
#endif
	}
	
	if (!optionMinOutput)
		::std::cout << ::std::endl;

	// RECURSION END
	for (int eSum = 0; eSum < maxErrors; ++eSum)
		for (int t = 0; t < maxT; ++t) 
		{
			TFloat recovered = _transform(0.0);
			for (int s = 0; s < statesCount; ++s)
			{
				TState const &state = states[s];

				// skip intermediate results
				if (state.intermediate || state.skipLast) continue;
				if (state.errors <= eSum)
				{
					int e = eSum - state.errors;
					// multiply probability for the trailing pattern
					_probAdd(recovered, _probMul(state.prob, (*colPrev)[(e*statesCount+s)*maxT+t]));
				}
			}

#ifndef COUNT_LOSSES
			// we can only normalize probs if t==0 contains all k-pattern probs
			if (t > 0 && !optionAbsolute)
				recovered = _probDiv(recovered, found[eSum*maxT]);
#endif

			found[eSum*maxT+t] = recovered;
		}
}

}

#endif
