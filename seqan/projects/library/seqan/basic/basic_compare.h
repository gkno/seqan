#ifndef SEQAN_HEADER_BASIC_COMPARE_H
#define SEQAN_HEADER_BASIC_COMPARE_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////

template <typename TLeft, typename TRight>
struct CompareType;

template <typename TLeft, typename TRight>
struct CompareType<TLeft const, TRight>
{
	typedef typename CompareType<TLeft, TRight>::Type const Type;
};
template <typename TLeft, typename TRight>
struct CompareType<TLeft, TRight const>
{
	typedef typename CompareType<TLeft, TRight>::Type const Type;
};
template <typename TLeft, typename TRight>
struct CompareType<TLeft const, TRight const>
{
	typedef typename CompareType<TLeft, TRight>::Type const Type;
};

//////////////////////////////////////////////////////////////////////////////

template<typename _T> inline
bool lexLess(const _T& _Left, const _T& _Right)
{	// return lexicographical _Left < _Right
    return 
        reinterpret_cast<typename _MakeUnsigned<_T>::Type const&>(_Left) <
        reinterpret_cast<typename _MakeUnsigned<_T>::Type const&>(_Right);
}

//////////////////////////////////////////////////////////////////////////////
}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
