#ifndef SEQAN_HEADER_BASIC_ITERATOR_ADAPT_STD_H
#define SEQAN_HEADER_BASIC_ITERATOR_ADAPT_STD_H

//////////////////////////////////////////////////////////////////////////////

namespace std
{
	template<typename TContainer, typename TSpec>
	struct iterator_traits<seqan::Iter<TContainer, TSpec> >
	{
		typedef seqan::Iter<TContainer, TSpec> TIter;

		typedef random_access_iterator_tag iterator_category;
		typedef typename seqan::Value<TIter>::Type value_type;
		typedef typename seqan::Difference<TIter>::Type difference_type;
		typedef typename seqan::Value<TIter>::Type * pointer;
		typedef typename seqan::Reference<TIter>::Type reference;
	};
}



namespace SEQAN_NAMESPACE_MAIN
{
//////////////////////////////////////////////////////////////////////////////
// ??? TODO: adaptors for std-iterators to seqan-iterators


//////////////////////////////////////////////////////////////////////////////

} //namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
