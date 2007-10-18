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

#ifndef SEQAN_HEADER_INDEX_ESA_BASE_H
#define SEQAN_HEADER_INDEX_ESA_BASE_H

namespace SEQAN_NAMESPACE_MAIN
{

	// dfs order
	struct _Preorder;
	struct _Postorder;

	template <typename TDFSOrder = _Postorder, typename THideEmptyEdges = True>
	struct VSTreeIteratorTraits {
		typedef TDFSOrder DFSOrder;
		typedef THideEmptyEdges HideEmptyEdges;
	};

/**
.Tag.Preorder:
..summary:Preorder depth-first-search.
..cat:Index
..signature:Preorder
..remarks:When given as a second parameter in @Function.goNext@ the Suffix Tree is traversed in a preorder fashion (visit the node before its children).
..see:Tag.Postorder
*/

/**
.Tag.Postorder:
..summary:Postorder depth-first-search.
..cat:Index
..signature:Postorder
..remarks:When given as a second parameter in @Function.goNext@ the Suffix Tree is traversed in a postorder fashion (visit the node after its children).
..see:Tag.Preorder
*/

	// predefined iterator traits
	struct Preorder:			VSTreeIteratorTraits<_Preorder,  True> {};
	struct Postorder:			VSTreeIteratorTraits<_Postorder, True> {};
	struct PreorderEmptyEdges:	VSTreeIteratorTraits<_Preorder,  False> {};	// also iterate over
	struct PostorderEmptyEdges:	VSTreeIteratorTraits<_Postorder, False> {};	// empty edges (with $-label)

	// traits for TopDown iterators (w/o ParentLinks) for which postorder/preorder is ignored
	struct HideEmptyEdges:		VSTreeIteratorTraits<_Postorder, True> {};
	struct EmptyEdges:			VSTreeIteratorTraits<_Postorder, False> {};	// empty edges (with $-label)
	
	// MultiMEMs are more specialized MaxRepeats
	template <typename TSpec = void>
	struct _MaxRepeats;	// base class
	struct _MultiMEMs;	// subclass tag



	// virtual suffix tree iterators
	template <typename TSpec = void>
	struct VSTree;

		// top down traversal iterators
		template <typename TSpec = Preorder>
		struct TopDown;						// starts in the suffix tree root and can go down and go right

			// allows an top-down iterator to go up
			template < typename TSpec = Preorder >
			struct ParentLinks {};			// .. can also go up

		// bottom up traversal iterators
		template <typename TSpec = Postorder>
		struct BottomUp;					// starts in the first node of a depth-first-search and can go next

			struct	SuperMaxRepeats;					// maximal repeat and not part of a longer repeat
			struct	SuperMaxRepeatsFast;
			struct	MUMs;								// Maximal Unique Match (unique in every sequence)

			typedef _MaxRepeats<void>		MaxRepeats;	// maximal repeat
			struct	MaxRepeatOccurences;
			typedef _MaxRepeats<_MultiMEMs> MultiMEMs;	// Multiple Maximal Exact Match
			struct	MultiMEMOccurences;					// i.e. maximal match over different sequences


/**
.Metafunction.GetVSTreeIteratorTraits:
..cat:Index
..summary:Default behaviour of @Function.goNext@ when no second parameter is given.
..signature:GetVSTreeIteratorTraits<TIterator>::Type
..param.TIterator:A @Spec.VSTree Iterator@.
..returns:$Tag.Postorder$ by default and $Tag.Preorder$ if $TIterator$ is $VSTree<TopDown<ParentLinks<> > >$ or $VSTree<TopDown<ParentLinks<Preorder> > >$.
*/

	template <typename TIterator>
	struct GetVSTreeIteratorTraits:
		DeepestSpec<TIterator> {};

//////////////////////////////////////////////////////////////////////////////

	template < typename TText, typename TSpec >
	struct VertexDescriptor< Index<TText, Index_ESA<TSpec> > > {
		typedef typename Size< Index<TText, Index_ESA<TSpec> > >::Type TSize;
		typedef Pair<
					Pair<TSize>,	// node range uniquely identifies every single suffix tree node
					TSize			// right boundary of parent node's range (allows to go right)
				> Type;
	};


//////////////////////////////////////////////////////////////////////////////
// needful forward declarations

	struct ArrayGaps;

	template <typename TSource, typename TSpec>
	class Align;


//////////////////////////////////////////////////////////////////////////////
// ESA fibres

/**
.Tag.ESA_Text:
..summary:The original text the index should be based on.
..general:Metafunction.Fibre
..cat:Index
..signature:Fibre<TIndex, ESA_Text>::Type
..param.TIndex:The ESA index type.
...type:Spec.Index_ESA
*/
/**
.Tag.ESA_RawText:
..summary:The raw text the index is really based on.
..general:Metafunction.Fibre
..cat:Index
..signature:Fibre<TIndex, ESA_RawText>::Type
..param.TIndex:The ESA index type.
...type:Spec.Index_ESA
..remarks:@Tag.ESA_Text@ and @Tag.ESA_RawText@ fibres are equal by default. 
They differ if the index text is a set of strings. Then, raw text is the concatenation of all strings in this set.
*/
/**
.Tag.ESA_SA:
..summary:The suffix array.
..general:Metafunction.Fibre
..cat:Index
..signature:Fibre<TIndex, ESA_SA>::Type
..param.TIndex:The ESA index type.
...type:Spec.Index_ESA
..remarks:The suffix array contains the indices of all suffices of @Tag.ESA_RawText@ in lexicographical order.
..remarks:A @Class.String@ over the alphabet of a size type.
*/
/**
.Tag.ESA_LCP:
..summary:The lcp table.
..general:Metafunction.Fibre
..cat:Index
..signature:Fibre<TIndex, ESA_LCP>::Type
..param.TIndex:The ESA index type.
...type:Spec.Index_ESA
..remarks:The lcp table contains the lcp-value of two adjacent suffices in the suffix array @Tag.ESA_SA@.
..returns:A @Class.String@ over the alphabet of a size type.
*/
/**
.Tag.ESA_ChildTab:
..summary:The child table.
..general:Metafunction.Fibre
..cat:Index
..signature:Fibre<TIndex, ESA_ChildTab>::Type
..param.TIndex:The ESA index type.
...type:Spec.Index_ESA
..remarks:The child table contains structural information of the suffix tree (see Abhouelda et al.).
..returns:A @Class.String@ over the alphabet of a size type.
*/
/**
.Tag.ESA_BWT:
..summary:The Burrows-Wheeler table.
..general:Metafunction.Fibre
..cat:Index
..signature:Fibre<TIndex, ESA_BWT>::Type
..param.TIndex:The ESA index type.
...type:Spec.Index_ESA
..remarks:The Burrows-Wheeler table contains the Burrows-Wheeler transformation of @Tag.ESA_RawText@.
..remarks:The entries are the characters left of the corresponding suffix in the suffix array @Tag.ESA_SA@.
..returns:The type @Tag.ESA_RawText@ returns.
*/

///.Metafunction.Fibre.param.TSpec.type:Tag.ESA_Text
///.Metafunction.Fibre.param.TSpec.type:Tag.ESA_RawText
///.Metafunction.Fibre.param.TSpec.type:Tag.ESA_SA
///.Metafunction.Fibre.param.TSpec.type:Tag.ESA_LCP
///.Metafunction.Fibre.param.TSpec.type:Tag.ESA_LCPE
///.Metafunction.Fibre.param.TSpec.type:Tag.ESA_ChildTab
///.Metafunction.Fibre.param.TSpec.type:Tag.ESA_BWT


	typedef Tag<_Fibre_Text> const		ESA_Text;
	typedef Tag<_Fibre_RawText> const	ESA_RawText;
	typedef Tag<_Fibre_SA> const		ESA_SA;
	typedef Tag<_Fibre_RawSA> const		ESA_RawSA;
	typedef Tag<_Fibre_SAE> const		ESA_SAE;
	typedef Tag<_Fibre_LCP> const		ESA_LCP;
	typedef Tag<_Fibre_LCPE> const		ESA_LCPE;
	typedef Tag<_Fibre_ChildTab> const	ESA_ChildTab;
	typedef Tag<_Fibre_BWT> const		ESA_BWT;


//////////////////////////////////////////////////////////////////////////////
// ESA index

/**
.Spec.Index_ESA:
..summary:An index based on an enhanced suffix array.
..cat:Index
..general:Class.Index
..signature:Index<TText, Index_ESA<> >
..param.TText:The text type.
...type:Class.String
..remarks:The fibres (see @Class.Index@ and @Metafunction.Fibre@) of this index are a suffix array (see @Tag.ESA_SA@), a lcp table (see @Tag.ESA_LCP@), etc.
..remarks:This index can be accessed as a Suffix Tree using the @Spec.VSTree Iterator@ classes.
*/

/*
	already defined in index_base.h

	template <typename TSpec = void>
	struct Index_ESA;
*/

	template < typename TText, typename TSpec >
	class Index<TText, Index_ESA<TSpec> > {
	public:
		Holder<typename Fibre<Index, ESA_Text>::Type>	text;
		typename Fibre<Index, ESA_SA>::Type				sa;
		typename Fibre<Index, ESA_LCP>::Type			lcp;
		typename Fibre<Index, ESA_LCPE>::Type			lcpe;
		typename Fibre<Index, ESA_ChildTab>::Type		childtab;
		typename Fibre<Index, ESA_BWT>::Type			bwt;

		Index() {}

		template <typename _TText>
		Index(_TText &_text):
			text(_text) {}

		template <typename _TText>
		Index(_TText const &_text):
			text(_text) {}
	};

//////////////////////////////////////////////////////////////////////////////
///.Function.clear.param.object.type:Class.Index

	template <typename TText, typename TSpec>
	inline void clear(Index<TText, Index_ESA<TSpec> > &index) {
		clear(getFibre(index, ESA_SA()));
		clear(getFibre(index, ESA_LCP()));
		clear(getFibre(index, ESA_LCPE()));
		clear(getFibre(index, ESA_ChildTab()));
		clear(getFibre(index, ESA_BWT()));
	}


//////////////////////////////////////////////////////////////////////////////
// open

	template < typename TObject, typename TSpec >
	inline bool open(
		Index< TObject, Index_ESA<TSpec> > &index, 
		const char *fileName,
		int openMode)
	{
		String<char> name;
		name = fileName;	append(name, ".txt");	
		if ((!open(getFibre(index, ESA_Text()), toCString(name), openMode)) && 
			(!open(getFibre(index, ESA_Text()), fileName), openMode)) return false;

		name = fileName;	append(name, ".sa");	open(getFibre(index, ESA_SA()), toCString(name), openMode);
		name = fileName;	append(name, ".lcp");	open(getFibre(index, ESA_LCP()), toCString(name), openMode);
		name = fileName;	append(name, ".child");	open(getFibre(index, ESA_ChildTab()), toCString(name), openMode);
		name = fileName;	append(name, ".bwt");	open(getFibre(index, ESA_BWT()), toCString(name), openMode);
		return true;
	}
	template < typename TObject, typename TSpec >
	inline bool open(
		Index< TObject, Index_ESA<TSpec> > &index, 
		const char *fileName) 
	{
		return open(index, fileName, OPEN_RDONLY);
	}


//////////////////////////////////////////////////////////////////////////////
// save

	template < typename TObject, typename TSpec >
	inline bool save(
		Index< TObject, Index_ESA<TSpec> > &index, 
		const char *fileName,
		int openMode)
	{
		String<char> name;
		name = fileName;	append(name, ".txt");	
		if ((!save(getFibre(index, ESA_Text()), toCString(name), openMode)) && 
			(!save(getFibre(index, ESA_Text()), fileName), openMode)) return false;

		name = fileName;	append(name, ".sa");	save(getFibre(index, ESA_SA()), toCString(name), openMode);
		name = fileName;	append(name, ".lcp");	save(getFibre(index, ESA_LCP()), toCString(name), openMode);
		name = fileName;	append(name, ".child");	save(getFibre(index, ESA_ChildTab()), toCString(name), openMode);
		name = fileName;	append(name, ".bwt");	save(getFibre(index, ESA_BWT()), toCString(name), openMode);
		return true;
	}
	template < typename TObject, typename TSpec >
	inline bool save(
		Index< TObject, Index_ESA<TSpec> > &index, 
		const char *fileName)
	{
		return save(index, fileName, OPEN_WRONLY | OPEN_CREATE);
	}

}

#endif
