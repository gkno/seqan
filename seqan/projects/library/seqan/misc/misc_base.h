/*
 *  misc_base.h
 *  SeqAn
 *
 *  Created by David Weese on 17.07.05.
 *
 */

#ifndef SEQAN_HEADER_MISC_BASE_H
#define SEQAN_HEADER_MISC_BASE_H

namespace SEQAN_NAMESPACE_MAIN
{

	//////////////////////////////////////////////////////////////////////////////
	// In SeqAn sets and maps store elements as pairs of (key,object) 
	// the elements of sets without objects are the keys.
	//////////////////////////////////////////////////////////////////////////////

	//////////////////////////////////////////////////////////////////////////////
	// Key/Object meta-functions
	//////////////////////////////////////////////////////////////////////////////

	template <typename TElement>
	struct Key {
		typedef TElement Type;
	};

	template <typename TKey, typename TObject, typename TSpec>
	struct Key< Pair<TKey, TObject, TSpec> > {
		typedef TKey Type;
	};

	template <typename TElement>
	struct Object {
		typedef Nothing Type;
	};

	template <typename TKey, typename TObject, typename TSpec>
	struct Object< Pair<TKey, TObject, TSpec> > {
		typedef TObject Type;
	};


	//////////////////////////////////////////////////////////////////////////////
	// keyOf function
	//////////////////////////////////////////////////////////////////////////////

	template <typename TElement>
	inline TElement & 
	keyOf(TElement &element) {
		return element;
	}
	template <typename TElement>
	inline TElement const & 
	keyOf(TElement const &element) {
		return element;
	}

	template <typename TKey, typename TObject, typename TSpec>
	inline TKey & 
	keyOf(Pair<TKey, TObject, TSpec> &element) {
		return element.i1;
	}
	template <typename TKey, typename TObject, typename TSpec>
	inline TKey const &
	keyOf(Pair<TKey, TObject, TSpec> const &element) {
		return element.i1;
	}

	//////////////////////////////////////////////////////////////////////////////
	// objectOf
	//////////////////////////////////////////////////////////////////////////////

	template <typename TKey, typename TObject, typename TSpec>
	inline TObject & 
	objectOf(Pair<TKey, TObject, TSpec> &element) {
		return element.i2;
	}
	template <typename TKey, typename TObject, typename TSpec>
	inline TObject const &
	objectOf(Pair<TKey, TObject, TSpec> const &element) {
		return element.i2;
	}

}

#endif

