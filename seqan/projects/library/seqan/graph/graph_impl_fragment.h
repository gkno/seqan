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
  
  Fragment(TId sqId1, TPos beg1, TId sqId2, TPos beg2, TSize l) :
    seqId1(sqId1), begin1(beg1), seqId2(sqId2), begin2(beg2), len(l) {
  }

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
	return seqId;
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
inline typename Size<Fragment<TId, TPos, TSize, TSpec> >::Type&
fragmentLength(Fragment<TId, TPos, TSize, TSpec> const& f,
	      TVal const seqId)
{
	SEQAN_CHECKPOINT
	return const_cast<typename Size<Fragment<TId, TPos, TSize, TSpec> >::Type&>(f.len);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TId, typename TPos, typename TSize, typename TSpec, typename TVal, typename TPosition>
inline typename Position<Fragment<TId, TPos, TSize, TSpec> >::Type
getProjectedPosition(Fragment<TId, TPos, TSize, TSpec> const& f,
					 TVal const seqId,
					 TPosition const pos)
{
	SEQAN_CHECKPOINT
	if ((TId) seqId == f.seqId1) {
		SEQAN_TASSERT((TPosition)f.begin1<=pos)
		SEQAN_TASSERT(pos - f.begin1 < f.len)	
		return f.begin2 + (pos - f.begin1);
	} else {
		SEQAN_TASSERT((TPosition)f.begin2<=pos)
		SEQAN_TASSERT(pos - f.begin2 < f.len)
		return f.begin1 + (pos - f.begin2);
	}
}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
