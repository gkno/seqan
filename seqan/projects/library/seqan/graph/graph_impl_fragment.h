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

#ifndef SEQAN_HEADER_GRAPH_IMPL_FRAGMENT_H
#define SEQAN_HEADER_GRAPH_IMPL_FRAGMENT_H

namespace SEQAN_NAMESPACE_MAIN
{


	
template<typename TId = unsigned int, typename TPos = unsigned int, typename TSize = unsigned int, typename TSpec = Default()>
class Fragment;



template<typename TId, typename TPos, typename TSize, typename TSpec>
class Fragment {
 public:
  TId seqId1;
  TPos begin1;
  TId seqId2;
  TPos begin2;
  TSize len;
  
  Fragment()
  {
  }

  Fragment(TId sqId1, TPos beg1, TId sqId2, TPos beg2, TSize l) :
    seqId1(sqId1), begin1(beg1), seqId2(sqId2), begin2(beg2), len(l) {
  }

};

template<typename TId, typename TPos, typename TSize, typename TSpec>
struct Position<Fragment<TId, TPos, TSize, TSpec> > {
	typedef TPos Type;
};

template<typename TId, typename TPos, typename TSize, typename TSpec>
struct Position<Fragment<TId, TPos, TSize, TSpec> const> {
	typedef TPos const Type;
};


//////////////////////////////////////////////////////////////////////////////

template<typename TId, typename TPos, typename TSize, typename TSpec, typename TStringSet, typename TVal>
inline typename Infix<typename Value<TStringSet>::Type>::Type
label(Fragment<TId, TPos, TSize, TSpec> const& f,
      TStringSet& str,
      TVal const seqId)
{
	SEQAN_CHECKPOINT
	typedef Fragment<TId, TPos, TSize, TSpec> TFragment;

	if ((TId) seqId == f.seqId1) {
		return infix(getValueById(str, (TId) seqId), f.begin1, f.begin1 + f.len);
	} else {
		return infix(getValueById(str, (TId) seqId), f.begin2, f.begin2 + f.len);
	}
}

//////////////////////////////////////////////////////////////////////////////

template<typename TId, typename TPos, typename TSize, typename TSpec, typename TVal>
inline TId
sequenceId(Fragment<TId, TPos, TSize, TSpec> const& f,
		   TVal const seqId)
{
	SEQAN_CHECKPOINT
	if ((TId) seqId == 0) {
		return f.seqId1;
	} else {
		return f.seqId2;
	}
}

//////////////////////////////////////////////////////////////////////////////

template<typename TId, typename TPos, typename TSize, typename TSpec, typename TVal>
inline typename Position<Fragment<TId, TPos, TSize, TSpec> >::Type&
fragmentBegin(Fragment<TId, TPos, TSize, TSpec> const& f,
			  TVal const seqId)
{
	SEQAN_CHECKPOINT
	if ((TId) seqId == f.seqId1) {
	  return const_cast<typename Position<Fragment<TId, TPos, TSize, TSpec> >::Type&>(f.begin1);
	} else {
	  return const_cast<typename Position<Fragment<TId, TPos, TSize, TSpec> >::Type&>(f.begin2);
	}
}

//////////////////////////////////////////////////////////////////////////////

template<typename TId, typename TPos, typename TSize, typename TSpec, typename TVal>
inline TSize&
fragmentLength(Fragment<TId, TPos, TSize, TSpec> const& f,
			   TVal const)
{
	SEQAN_CHECKPOINT
	return const_cast<TSize&>(f.len);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TId, typename TPos, typename TSize, typename TSpec>
inline TSize&
fragmentLength(Fragment<TId, TPos, TSize, TSpec> const& f)
{
	SEQAN_CHECKPOINT
	return const_cast<TSize&>(f.len);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TId, typename TPos, typename TSize, typename TSpec, typename TId1, typename TPosition1, typename TId2, typename TPosition2>
inline void
getProjectedPosition(Fragment<TId, TPos, TSize, TSpec> const& f,
					 TId1 const seqId,
					 TPosition1 const pos,
					 TId2& seqId2,
					 TPosition2& pos2)
{
	SEQAN_CHECKPOINT
	if ((TId) seqId == f.seqId1) {
		SEQAN_TASSERT((TPosition1)f.begin1<=pos)
		SEQAN_TASSERT(pos - f.begin1 < f.len)	
		pos2 = f.begin2 + (pos - f.begin1);
		seqId2 = f.seqId2;
		return;
	} else {
		SEQAN_TASSERT((TPosition1)f.begin2<=pos)
		SEQAN_TASSERT(pos - f.begin2 < f.len)
		pos2 = f.begin1 + (pos - f.begin2);
		seqId2 = f.seqId1;
		return;
	}
}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
