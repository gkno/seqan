/*
 *  index_esa_base.h
 *  SeqAn
 *
 *  Created by David Weese on 17.07.05.
 *
 */

#ifndef SEQAN_HEADER_INDEX_ESA_BASE_H
#define SEQAN_HEADER_INDEX_ESA_BASE_H

namespace SEQAN_NAMESPACE_MAIN
{
	
	// virtual suffix tree iterators
	template <typename TSpec = void>
	struct VSTree;


		// top down traversal iterators
		template <typename TSpec = void>
		struct TopDown;

				struct _Preorder;
				typedef Tag<_Preorder> Preorder;

			// allows an top-down iterator to go up
			template < typename TSpec = Preorder >
			struct ParentLinks {};

				// dfs order
				struct _Postorder;
				//struct _Preorder;

				typedef Tag<_Postorder> Postorder;
				//typedef Tag<_Preorder> Preorder;


		// bottom up traversal iterators
		template <typename TSpec = void>
		struct BottomUp;

			// bottom up repeat search iterators
			struct SuperMaxRepeats;
			struct SuperMaxRepeatsFast;
			struct MaxRepeats;
			struct MUMs;
			struct MaxRepeatOccurences;


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

/**
.Metafunction.DefaultDFSOrder:
..cat:Index
..summary:Default behaviour of @Function.goNext@ when no second parameter is given.
..signature:DefaultDFSOrder<TIterator>::Type
..param.TIterator:A @Spec.VSTree Iterator@.
..returns:$Tag.Postorder$ by default and $Tag.Preorder$ if $TIterator$ is $VSTree<TopDown<ParentLinks<> > >$ or $VSTree<TopDown<ParentLinks<Preorder> > >$.
*/
    template < typename TIterator >
    struct DefaultDFSOrder {
        typedef Postorder Type;
    };

    template < typename TIndex >
    struct DefaultDFSOrder< Iter< TIndex, VSTree< TopDown< ParentLinks<Preorder> > > > > {
        typedef Preorder Type;
    };

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


	typedef Tag<_Fibre_Text>		ESA_Text;
	typedef Tag<_Fibre_RawText>		ESA_RawText;
	typedef Tag<_Fibre_SA>			ESA_SA;
	typedef Tag<_Fibre_RawSA>		ESA_RawSA;
	typedef Tag<_Fibre_SAE>			ESA_SAE;
	typedef Tag<_Fibre_LCP>			ESA_LCP;
	typedef Tag<_Fibre_LCPE>		ESA_LCPE;
	typedef Tag<_Fibre_ChildTab>	ESA_ChildTab;
	typedef Tag<_Fibre_BWT>			ESA_BWT;


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
// ESA finders

	struct _Finder_MLR;		// simple Suffix Array finder with mlr-heuristic
	struct _Finder_LCPE;	// Suffix Array finder using an enhanced LCP-Table
	struct _Finder_LCPH;	// hybrid of Suffix Array and Enhanced LCP (to reduce random accesses)

/**
.Tag.ESA_FIND_MLR:
..summary:Exact string matching using a suffix array binary search with the mlr-heuristic.
..general:Class.Finder
..cat:Index
..signature:Finder<TIndex, ESA_FIND_MLR>
..param.TIndex:The index type.
...type:Spec.Index_ESA
*/

	typedef Tag<_Finder_MLR>	ESA_FIND_MLR;

/**
.Tag.ESA_FIND_LCPE:
..summary:Exact string matching using a suffix array binary search and a lcp-interval tree.
..general:Class.Finder
..cat:Index
..signature:Finder<TIndex, ESA_FIND_LCPE>
..param.TIndex:The index type.
...type:Spec.Index_ESA
*/

	typedef Tag<_Finder_LCPE>	ESA_FIND_LCPE;
	typedef Tag<_Finder_LCPH>	ESA_FIND_LCPH;

	template < typename TText, typename TSpec >
	struct DefaultFinder<Index<TText, Index_ESA<TSpec> > > {
        typedef Tag<_Finder_MLR> Type;	// standard suffix array finder is mlr-heuristic
    };


//////////////////////////////////////////////////////////////////////////////

	template < typename TText, typename TSpecESA, typename TSpecFinder >
	class Finder< Index<TText, Index_ESA<TSpecESA> >, TSpecFinder >
	{
	public:
		typedef Index<TText, Index_ESA<TSpecESA> > TIndex;
		typedef typename Iterator< typename Fibre<TIndex, ESA_SA>::Type const >::Type TIterator;

		Holder<TIndex>	index;
		Pair<TIterator>	range;
		TIterator		data_iterator;

		Finder() 
		{
			clear(*this);
		}
		Finder(TIndex &_index): index(_index) 
		{
			clear(*this);
		}
		Finder(TIndex const &_index): index(_index)
		{
			clear(*this);
		}


//____________________________________________________________________________

		friend inline typename _Parameter<TIndex>::Type 
		host(Finder & me)
		{
SEQAN_CHECKPOINT
			return value(me.index);
		}

		friend inline typename _Parameter<TIndex>::Type 
		host(Finder const & me)
		{
SEQAN_CHECKPOINT
			return value(me.index);
		}

		friend inline typename _Parameter<TIndex>::Type 
		container(Finder & me)
		{
SEQAN_CHECKPOINT
			return value(me.index);
		}

		friend inline typename _Parameter<TIndex>::Type 
		container(Finder const & me)
		{
SEQAN_CHECKPOINT
			return value(me.index);
		}

//____________________________________________________________________________

		friend inline void
		setHost(Finder & me, typename _Parameter<TIndex>::Type container_)
		{
SEQAN_CHECKPOINT
			me.index = container;
		}

		friend inline void
		setContainer(Finder & me, typename _Parameter<TIndex>::Type container_)
		{
SEQAN_CHECKPOINT
			me.index = container;
		}

//____________________________________________________________________________

		friend inline TIterator &
		hostIterator(Finder & me)
		{
SEQAN_CHECKPOINT
			return me.data_iterator;
		}

		friend inline TIterator const &
		hostIterator(Finder const & me)
		{
SEQAN_CHECKPOINT
			return me.data_iterator;
		}


//____________________________________________________________________________

		friend inline bool
		empty(Finder & me)
		{
SEQAN_CHECKPOINT
			return me.range.i1 == me.range.i2;
		}

		friend inline void
		clear(Finder & me)
		{
SEQAN_CHECKPOINT
			me.range.i1 = me.range.i2 = TIterator();
		}

//____________________________________________________________________________

		friend inline bool
		atBegin(Finder & me)
		{
SEQAN_CHECKPOINT
			return (empty(me) || hostIterator(me) == me.range.i1);
		}

		friend inline bool
		atEnd(Finder & me)
		{
SEQAN_CHECKPOINT
			return (empty(me) || hostIterator(me) == me.range.i2);
		}

//____________________________________________________________________________

		friend inline void
		goBegin(Finder & me)
		{
SEQAN_CHECKPOINT
			hostIterator(me) = me.range.i1;
		}

		friend inline void
		goEnd(Finder & me)
		{
SEQAN_CHECKPOINT
			hostIterator(me) = me.range.i2;
		}

//____________________________________________________________________________
/*
		template <typename TPosition>
		friend inline void 
		setPosition(Finder & me, TPosition pos_)
		{
SEQAN_CHECKPOINT
			hostIterator(me) = me.range.i1 + pos_;
		}
*/
//____________________________________________________________________________

		friend inline typename Position<Finder>::Type
		position(Finder & me)
		{
SEQAN_CHECKPOINT
			if (empty(me)) return 0;
			return *me.data_iterator;
		}

		friend inline typename Position<Finder>::Type
		position(Finder const & me)
		{
SEQAN_CHECKPOINT
			if (empty(me)) return 0;
			return hostIterator(me) - begin(container(me), Rooted());
		}

	};


//////////////////////////////////////////////////////////////////////////////
// find implementation

	template < typename TFinder, typename TPattern >
	inline void _findFirstESAIndex(
		TFinder &finder,
		TPattern const &pattern,
		ESA_FIND_MLR const)
	{
		typename Haystack<TFinder>::Type &index = haystack(finder);
		indexRequire(index, ESA_SA());
		finder.range = equalRangeSA(indexRawText(index), indexSA(index), pattern);
	}

	template < typename TFinder, typename TPattern >
	inline void _findFirstESAIndex(
		TFinder &finder,
		TPattern const &pattern,
		ESA_FIND_LCPE const)
	{
		typename Haystack<TFinder>::Type &index = haystack(finder);
		indexRequire(index, ESA_SA());
		indexRequire(index, ESA_LCPE());
		finder.range = equalRangeLCPE(indexRawText(index), indexSA(index), indexLCPE(index), pattern);
	}

//////////////////////////////////////////////////////////////////////////////
// find

	template < typename TText, typename TSpecESA, typename TSpecFinder, typename TPattern >
	inline bool find(
		Finder<Index<TText, Index_ESA<TSpecESA> >, TSpecFinder> &finder,
		TPattern const &pattern)
	{
		if (empty(finder)) {
			_findFirstESAIndex(finder, needle(pattern), TSpecFinder());
			hostIterator(finder) = finder.range.i1;
		} else
			++hostIterator(finder);
		return !atEnd(finder);
	}

	template < typename TText, typename TSpecESA, typename TSpecFinder >
	inline bool find(Finder<Index<TText, Index_ESA<TSpecESA> >, TSpecFinder> &finder)
	{
		if (empty(finder)) return false;
		++hostIterator(finder);
		return !atEnd(finder);
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
