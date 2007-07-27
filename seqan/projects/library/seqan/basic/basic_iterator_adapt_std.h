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
