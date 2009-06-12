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
==========================================================================*/

#ifndef SEQAN_HEADER_SEQAN_CONSENSUS_PROFILE_H
#define SEQAN_HEADER_SEQAN_CONSENSUS_PROFILE_H


namespace SEQAN_NAMESPACE_MAIN
{


template<typename TValue, typename TCount = unsigned int, typename TSpec = Default>
class ProfileType;


template<typename TValue, typename TCount, typename TSpec>
class ProfileType {
	public:
		typedef typename Size<ProfileType>::Type TSize;

		TCount count[ValueSize<TValue>::VALUE];

		ProfileType() {
			for(TSize i = 0; i<ValueSize<TValue>::VALUE; ++i) count[i] = 0;
		}

		~ProfileType() {}

		ProfileType(ProfileType const& _other) {
			for(TSize i = 0; i<ValueSize<TValue>::VALUE; ++i) count[i] = _other.count[i];
		}

		template <typename TOther> 
		ProfileType(TOther const& other_data) {
			for(TSize i = 0; i<ValueSize<TValue>::VALUE; ++i) count[i] = 0;
			count[ordValue(TValue(other_data))] = 1;
		}

		ProfileType const& 
		operator = (ProfileType const& other_data) 
		{
			if (this == &other_data) return *this;
			for(TSize i = 0; i<ValueSize<TValue>::VALUE; ++i) count[i] = other_data.count[i];
			return *this;
		}

		template <typename TOther> 
		ProfileType const& 
		operator = (TOther const& other_data) {
			for(TSize i = 0; i<ValueSize<TValue>::VALUE; ++i) count[i] = 0;
			count[ordValue(TValue(other_data))] = 1;
			return *this;
		}

		operator char() { 
			return (char) TValue(*this); 
		}

		bool operator==(ProfileType const& other_data) const {
			for(TSize i = 0; i<ValueSize<TValue>::VALUE; ++i) {
				if (count[i] != other_data.count[i]) return false;
			}
			return true;
		}

		bool operator!=(ProfileType const& other_data) const {
			for(TSize i = 0; i<ValueSize<TValue>::VALUE; ++i) {
				if (count[i] != other_data.count[i]) return true;
			}
			return false;
		}
};


//////////////////////////////////////////////////////////////////////////////

template<typename TValue, typename TCount, typename TSpec>
struct ValueSize<ProfileType<TValue, TCount, TSpec> >
{
	enum { VALUE = ValueSize<TValue>::VALUE };
};

////////////////////////////////////////////////////////////////////////////////

template <typename TSourceValue, typename TSourceCount, typename TSourceSpec>
inline bool
empty(ProfileType<TSourceValue, TSourceCount, TSourceSpec> const & source)
{
	typedef typename Size<ProfileType<TSourceValue, TSourceCount, TSourceSpec> const >::Type TSize;
	for(TSize i = 0; i<ValueSize<TSourceValue>::VALUE; ++i) {
		if (source.count[i]) return false;
	}
	return true;
}

// Empty iff no letters or only gaps
template <typename TSourceValue, typename TSourceCount, typename TSourceSpec, typename TGapOrdValue>
inline bool
empty(ProfileType<TSourceValue, TSourceCount, TSourceSpec> const & source,
	  TGapOrdValue gapOrdVal)
{
	typedef typename Size<ProfileType<TSourceValue, TSourceCount, TSourceSpec> const >::Type TSize;
	for(TSize i = 0; i<ValueSize<TSourceValue>::VALUE; ++i) {
		if ((i != gapOrdVal) && (source.count[i])) return false;
	}
	return true;
}

////////////////////////////////////////////////////////////////////////////////

template <typename TSourceValue, typename TSourceCount, typename TSourceSpec>
inline typename Size<ProfileType<TSourceValue, TSourceCount, TSourceSpec> const >::Type
_getMaxIndex(ProfileType<TSourceValue, TSourceCount, TSourceSpec> const & source)
{
	typedef typename Size<ProfileType<TSourceValue, TSourceCount, TSourceSpec> const >::Type TSize;
	TSize maxIndex = 0;
	TSourceCount maxCount = source.count[0];
	for(TSize i = 1; i<ValueSize<TSourceValue>::VALUE; ++i) {
		if (source.count[i] > maxCount) {
			maxIndex = i;
			maxCount = source.count[i];
		}
	}
	return maxIndex;
}


////////////////////////////////////////////////////////////////////////////////

template <typename TTargetValue, typename TTargetSpec, typename TSourceValue, typename TSourceCount, typename TSourceSpec>
inline void 
assign(SimpleType<TTargetValue, TTargetSpec> & target, 
	   ProfileType<TSourceValue, TSourceCount, TSourceSpec> const & source)
{
	target.value = _getMaxIndex(source);
}

////////////////////////////////////////////////////////////////////////////////

template <typename THost, char CHAR, typename TSpec, typename T, typename TSourceValue, typename TSourceCount, typename TSourceSpec>
inline typename Convert<ModifiedAlphabet<THost, ModExpand<CHAR, TSpec> >, ProfileType<TSourceValue, TSourceCount, TSourceSpec> >::Type
convertImpl(Convert<ModifiedAlphabet<THost, ModExpand<CHAR, TSpec> >, T> const,
			ProfileType<TSourceValue, TSourceCount, TSourceSpec> const & source)
{
	ModifiedAlphabet<THost, ModExpand<CHAR, TSpec> > my;
	my.data = _getMaxIndex(source);
	return my;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TTarget, typename T, typename TSourceValue, typename TSourceCount, typename TSourceSpec>
inline typename Convert<TTarget, ProfileType<TSourceValue, TSourceCount, TSourceSpec> >::Type
convertImpl(Convert<TTarget, T> const,
			ProfileType<TSourceValue, TSourceCount, TSourceSpec> const & source)
{
	return convertImpl(Convert<TTarget, T>(), TSourceValue(source));
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStream, typename TValue, typename TCount, typename TSpec>
TStream& operator<<(TStream& os, ProfileType<TValue, TCount, TSpec> const& rhs) {
	typedef typename Size<ProfileType<TValue, TCount, TSpec> const>::Type TSize;
	for(TSize i = 0; i<ValueSize<TValue>::VALUE; ++i) {
		os << i << ':' << ' ' << rhs.count[i] << std::endl;
	}
	return os;
}



//////////////////////////////////////////////////////////////////////////////

struct ConsensusScore_;
typedef Tag<ConsensusScore_> const ConsensusScore;

//////////////////////////////////////////////////////////////////////////////

template <typename TValue>
class Score<TValue, ConsensusScore>
{
public:
	TValue data_gap_extend;
	TValue data_gap_open;

public:
	Score():
		data_gap_extend(-1),
		data_gap_open(-1)
	{}

};


template <typename TValue, typename TSourceValue, typename TSourceCount, typename TSourceSpec>
inline TValue
score(Score<TValue, ConsensusScore> const &,
	  ProfileType<TSourceValue, TSourceCount, TSourceSpec> left,   // Consensus profile
	  ProfileType<TSourceValue, TSourceCount, TSourceSpec> right)  // Single string
{
	typedef ProfileType<TSourceValue, TSourceCount, TSourceSpec> TProfileType;
	typedef typename Size<TProfileType>::Type TSize;
	TSize indexRight = 0;
	for(TSize i = 0; i<ValueSize<TProfileType>::VALUE; ++i) {
		if (right.count[i] > 0) {
			indexRight = i;
			break;
		}
	}
	TSourceCount maxCount = left.count[indexRight];
	for(TSize i = 0; i<ValueSize<TProfileType>::VALUE; ++i) {
		if (left.count[i] > maxCount) return -1;
	}
	return 0;
}

}

#endif

