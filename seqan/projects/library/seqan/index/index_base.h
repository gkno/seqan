/*
 *  index_base.h
 *  genindex
 *
 *  Created by David Weese on 17.07.05.
 *
 */

#ifndef SEQAN_HEADER_INDEX_BASE_H
#define SEQAN_HEADER_INDEX_BASE_H

//#define SEQAN_TEST_INDEX

namespace SEQAN_NAMESPACE_MAIN
{

	// index construction specs
	struct Skew3;
	struct Skew7;
	struct ManberMyers;
	struct QGram1;
	struct QGram2;
	struct QGram3;

	// lcp table construction algorithms
	struct Kasai;
	struct KasaiInPlace;

	// enhanced suffix array construction algorithms
	struct ChildTab;
	struct BWT;

	// Index/ExternalIndex specs
	// different types of suffix array based indexes ordered by their complexity/size
	struct Index_QGram;			// suffix array sorted by the first q chars

	template <typename TSpec = void>
	struct Index_ESA;

/**
.Metafunction.DefaultIndexSpec:
..cat:Index
..summary:Default @Class.Index@ specialization type.
..signature:DefaultIndexSpec<TText>::Type
..param.TText:The given text type.
..returns:Can be @Spec.Index_ESA@ or $Index_QGram$, etc.
..remarks:Currently @Spec.Index_ESA@ is default if $TText$ is a @Class.String@.
*/
    template < typename TObject >
    struct DefaultIndexSpec {
        typedef Index_ESA<> Type;
    };

/**
.Metafunction.DefaultIndexStringSpec:
..cat:Index
..summary:Default @Class.String@ specialization type of the @Metafunction.Fibre@ of an @Class.Index@. 
..signature:DefaultIndexStringSpec<TIndex>::Type
..param.TIndex:An @Class.Index@ Type.
..returns:If the underlying text is a @Class.String@ or a set of Strings (see @Class.StringSet@) the String's spec. type is returned.
..remarks:Most of the @Class.Index@ fibres are strings. The @Class.String@ specialization type is chosen by this meta-function.
*/
    template < typename TIndex >
    struct DefaultIndexStringSpec {
        typedef Alloc<> Type;
    };

    template < typename TValue, typename TSpec >
    struct DefaultIndexStringSpec< String<TValue, TSpec> > {
        typedef TSpec Type;
    };

	template < typename TString, typename TSpec >
    struct DefaultIndexStringSpec< StringSet<TString, TSpec> > {
		typedef typename DefaultIndexStringSpec<TString>::Type Type;
    };


/**
.Class.Index:
..summary:Contains preprocessing data of a fixed text. Allows fast dictionary look-up and advanced computations.
..cat:Pattern Matching
..signature:Index<TText[, TSpec]>
..param.TText:The text type.
...type:Class.String
..param.TSpec:The index type.
...default:The result of @Metafunction.DefaultIndexSpec@
..remarks:An index contains various arrays or objects, also called fibres (see @Metafunction.Fibre@).
..remarks:These fibres are created on demand depending on the requirements of an algorithm.
*/

///.Function.setHaystack.param.haystack.type:Class.Index

	// index as a haystack
	template < 
        typename TObject, 
        typename TSpec = typename DefaultIndexSpec<TObject>::Type > 
	class Index;


	template <typename _TText>
	struct BWTStruct {
		typedef _TText						TText;
		typedef typename Size<_TText>::Type	TSize;

		TText tab;
		TSize undefined;	// index of undefined entry (sa[undefined] == 0)
		
		BWTStruct() {}
		template <typename T1>
		BWTStruct(T1 &t1): tab(t1) {}
		template <typename T1, typename T2>
		BWTStruct(T1 const &t1): tab(t1) {}
		template <typename T1, typename T2>
		BWTStruct(T1 &t1, T2 &t2): tab(t1, t2) {}
		template <typename T1, typename T2>
		BWTStruct(T1 const &t1, T2 const &t2): tab(t1, t2) {}
	};

	template <typename TRawText>
	struct Value< BWTStruct<TRawText> > {
		typedef typename Value<TRawText>::Type Type;
	};



/**
.Metafunction.Fibre:
..summary:Type of a specific bundle member (fibre).
..signature:Fibre<TIndex, TSpec>::Type
..cat:Index
..param.TIndex:The fibre container type.
..param.TSpec:Type to specify the fibre.
..returns:Fibre type.
..remarks:An @Class.Index@ can be seen as a bundle consisting of various fibres. In most cases this type is $String<Size<TIndex>::Type>$.
..remarks:A @Metafunction.Fibre@ need not to be a real container. It can also be view (see @Tag.ESA_RawText@).
*/
	// meta function to get the type of a bundle fibre
	template < typename TIndex, typename TSpec >
	struct Fibre {
		typedef String< typename Size<TIndex>::Type > Type;
	};

	struct FibreRecord {
		unsigned	id;
		void*		ptr;
		bool		owner;
	};

	// less function to search in sorted list for fibre id
	struct FibreLess: public ::std::binary_function<FibreRecord, unsigned, bool>
	{	// functor for operator>
		inline bool operator()(const FibreRecord& _Left, const unsigned _Right) const
		{	// apply operator> to operands
			return (_Left.id < _Right);
		}
	};

/**
.Metafunction.DefaultIndexCreator:
..cat:Index
..summary:Default algorithm to create a demanded and not yet existing @Metafunction.Fibre@.
..signature:DefaultIndexCreator<TIndex, TFibre>::Type
..param.TIndex:An @Class.Index@ Type.
..param.TFibre:A tag specifying the fibre (e.g. @Tag.ESA_SA@).
..returns:A tag specifying the default algorithm to create the fibre with.
*/
    // standard algorithm for indices creation
    template < typename TIndex, typename TFibre >
    struct DefaultIndexCreator;


/**
	.Class.Bundle:
	..summary:General purpose container of various members.
	..signature:Bundle<TValue, TSize>
	..param.TValue:The value type, that is the type of the items/characters stored in the string.
	...remarks:Use @Metafunction.Value@ to get the value type for a given class.
	..param.TSpec:The specializing type.
	...default:$Alloc<>$, see @Spec.Alloc String@.
*/
/*
	template < typename TSpec = void >
	struct Bundle {
		typedef ::std::vector<FibreRecord>	TFibreRecords;
		TFibreRecords						fibres;
	};

	template < typename TBundleSpec, typename TFibreSpec >
	inline FibreRecord& getRecord(Bundle<TBundleSpec> &bundle, TFibreSpec const) {
		unsigned id = (unsigned)_ClassIdentifier<TFibreSpec>::getID();

		typename Bundle<TBundleSpec>::TFibreRecords::iterator first = lower_bound(bundle.fibres.begin(), bundle.fibres.end(), id, FibreLess());
		if (!first->id != id) {
			FibreRecord rec;
			rec.id = id;
			rec.ptr = NULL;
			rec.owner = true;
			bundle.fibres.insert(first, rec);
		} else
			return *first;
	}

	template < typename TBundleSpec, typename TFibreSpec >
	inline typename Fibre<Bundle<TBundleSpec>, TFibreSpec>::Type & getFibre(Bundle<TBundleSpec> &bundle, TFibreSpec const) {
		typedef typename Fibre<Bundle<TBundleSpec>, TFibreSpec>::Type Type;
		unsigned id = (unsigned)_ClassIdentifier<TFibreSpec>::getID();

		FibreRecord &rec = getRecord(bundle, TFibreSpec());
		if (!rec.ptr)
			rec.ptr = new Type();
		return *reinterpret_cast<Type*>(rec.ptr);
	}

	template < typename TBundleSpec, typename TFibreSpec >
	inline typename Fibre<Bundle<TBundleSpec>, TFibreSpec>::Type const & getFibre(Bundle<TBundleSpec> const &bundle, TFibreSpec const) {
		typedef typename Fibre<Bundle<TBundleSpec>, TFibreSpec>::Type Type;
		unsigned id = (unsigned)_ClassIdentifier<TFibreSpec>::getID();

		FibreRecord &rec = getRecord(bundle, TFibreSpec());
		return *reinterpret_cast<Type*>(rec.ptr);
	}
*/

	// various fibre specs for enhanced suffix arrays

	struct _ESA_Text;		// Original text. Can be a String or a StringSet
	struct _ESA_RawText;	// Concatenation of the strings above
	struct _ESA_SA;			// suffix array of raw text
	struct _ESA_SAE;		// suffix array reordered in a b-tree
	struct _ESA_LCP;		// lcp table of raw text
	struct _ESA_LCPE;		// lcp interval tree
	struct _ESA_ChildTab;	// childtab (Kurtz et al.) of raw text
	struct _ESA_BWT;		// burrows wheeler table of raw text

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


	typedef Tag<_ESA_Text>		ESA_Text;
	typedef Tag<_ESA_RawText>	ESA_RawText;
	typedef Tag<_ESA_SA>		ESA_SA;
	typedef Tag<_ESA_SAE>		ESA_SAE;
	typedef Tag<_ESA_LCP>		ESA_LCP;
	typedef Tag<_ESA_LCPE>		ESA_LCPE;
	typedef Tag<_ESA_ChildTab>	ESA_ChildTab;
	typedef Tag<_ESA_BWT>		ESA_BWT;

	template <typename TObject>
	struct SAValue {
		typedef typename Size<TObject>::Type Type;
	};
	
	template < typename TText, typename TSpec >
	struct SAValue< Index<TText, TSpec> > {
		typedef typename SAValue<TText>::Type Type;
	};

	// to speed up sequence number computation
	// we use a pair of seqNo and localPosition
/*	template < typename TString, typename TSpec, typename TIndexSpec >
	struct SAValue< Index<StringSet<TString, TSpec>, TIndexSpec> > {
		typedef Pair<
			typename Size< StringSet<TString, TSpec> >::Type,
			typename SAValue<TString>::Type
		> Type;
	};
*/


	// specialization for ESA indices

	template < typename TObject, typename TSpec >
    struct DefaultIndexStringSpec< Index<TObject, TSpec> > {
		typedef typename DefaultIndexStringSpec<TObject>::Type Type;
    };


	template < typename TObject, typename TSpec, typename TFibreSpec >
	struct Fibre< Index<TObject, Index_ESA<TSpec> >, TFibreSpec > {
		typedef String< 
			typename Size< Index<TObject, Index_ESA<TSpec> > >::Type,
			typename DefaultIndexStringSpec< Index<TObject, Index_ESA<> > >::Type > Type;
	};

	template < typename TText, typename TSpec >
	struct Fibre< Index<TText, Index_ESA<TSpec> >, ESA_Text > {
		typedef TText Type;
	};

	template < typename TText, typename TSpec >
	struct Fibre< Index<TText, Index_ESA<TSpec> >, ESA_RawText > {
		typedef typename Concatenator<TText>::Type Type;
	};

	template < typename TText, typename TSpec >
	struct Fibre< Index<TText, Index_ESA<TSpec> >, ESA_SA > {
		typedef String<
			typename SAValue< Index<TText, Index_ESA<TSpec> > >::Type,
			typename DefaultIndexStringSpec< Index<TText, Index_ESA<> > >::Type > Type;
	};

	template < typename TText, typename TSpec >
	struct Fibre< Index<TText, Index_ESA<TSpec> >, ESA_BWT > {
		typedef BWTStruct< typename Fibre< Index<TText, Index_ESA<TSpec> >, ESA_RawText >::Type > Type;
	};


// compression speed hack
/*	template < typename TValue, typename TConf >
	struct Fibre< Bundle<Index_ESA< String<TValue, External<TConf> > > >, ESA_LCP > {
		typedef String<unsigned char, External<TConf> > Type;
	};
*/	

	
	// ESA creators

	template < typename TText, typename TSpec >
	struct DefaultIndexCreator<Index<TText, TSpec>, ESA_SA> {
        typedef Skew7 Type;							// standard suffix array creator is skew7
    };

	template < typename TText, typename TSpec >
	struct DefaultIndexCreator<Index<TText, TSpec>, ESA_LCP> {
        typedef KasaiInPlace Type;
    };

	template < typename TText, typename TSpec >
	struct DefaultIndexCreator<Index<TText, TSpec>, ESA_BWT> {
        typedef BWT Type;
    };

	template < typename TText, typename TSpec >
	struct DefaultIndexCreator<Index<TText, TSpec>, ESA_ChildTab> {
        typedef ChildTab Type;
    };


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

    template < typename TText, typename TSpec >
    struct Value< Index<TText, Index_ESA<TSpec> > > {
		typedef typename Value< typename Fibre< Index<TText, Index_ESA<TSpec> >, ESA_RawText >::Type >::Type Type;
    };

	template < typename TText, typename TSpec >
    struct Size< Index<TText, Index_ESA<TSpec> > > {
		typedef typename Size< typename Fibre< Index<TText, Index_ESA<TSpec> >, ESA_RawText >::Type >::Type Type;
    };


	// ESA finders

	struct _ESA_MLR;		// simple Suffix Array finder with mlr-heuristic
	struct _ESA_LCPE;		// Suffix Array finder using an enhanced LCP-Table
	struct _ESA_LCPH;	    // hybrid of Suffix Array and Enhanced LCP (to reduce random accesses)

/**
.Tag.ESA_MLR:
..summary:Exact string matching using a suffix array binary search with the mlr-heuristic.
..general:Class.Finder
..cat:Index
..signature:Finder<TIndex, ESA_MLR>
..param.TIndex:The index type.
...type:Spec.Index_ESA
*/

	typedef Tag<_ESA_MLR>	ESA_MLR;

/**
.Tag.ESA_LCPE:
..summary:Exact string matching using a suffix array binary search and a lcp-interval tree.
..general:Class.Finder
..cat:Index
..signature:Finder<TIndex, ESA_LCPE>
..param.TIndex:The index type.
...type:Spec.Index_ESA
*/

	typedef Tag<_ESA_LCPE>	ESA_LCPE;
	typedef Tag<_ESA_LCPH>	ESA_LCPH;

	template < typename TText, typename TSpec >
	struct DefaultFinder<Index<TText, Index_ESA<TSpec> > > {
        typedef ESA_MLR Type;	// standard suffix array finder is mlr-heuristic
    };


	// direct access to the enhanced suffix array tables

/**
.Function.getFibre:
..summary:Returns a specific @Metafunction.Fibre@ of an @Class.Index@ object.
..cat:Index
..signature:getFibre(index, fibre_tag)
..param.index:The @Class.Index@ object holding the fibre.
...type:Class.Index
..param.fibre_tag:A tag that identifies the @Metafunction.Fibre@ (e.g. @Tag.ESA_SA@).
..returns:A reference to the @Metafunction.Fibre@ object.
*/

	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, Index_ESA<TSpec> >, ESA_Text>::Type & getFibre(Index<TText, Index_ESA<TSpec> > &index, ESA_Text const) {
		return value(index.text);
	}
	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, Index_ESA<TSpec> >, ESA_Text>::Type const & getFibre(Index<TText, Index_ESA<TSpec> > const &index, ESA_Text const) {
		return value(index.text);
	}

	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, Index_ESA<TSpec> >, ESA_RawText>::Type & getFibre(Index<TText, Index_ESA<TSpec> > &index, ESA_RawText const) {
		return concat(value(index.text));
	}
	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, Index_ESA<TSpec> >, ESA_RawText>::Type const & getFibre(Index<TText, Index_ESA<TSpec> > const &index, ESA_RawText const) {
		return concat(value(index.text));
	}

	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, Index_ESA<TSpec> >, ESA_SA>::Type & getFibre(Index<TText, Index_ESA<TSpec> > &index, ESA_SA const) {
		return index.sa;
	}
	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, Index_ESA<TSpec> >, ESA_SA>::Type const & getFibre(Index<TText, Index_ESA<TSpec> > const &index, ESA_SA const) {
		return index.sa;
	}

	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, Index_ESA<TSpec> >, ESA_LCP>::Type & getFibre(Index<TText, Index_ESA<TSpec> > &index, ESA_LCP const) {
		return index.lcp;
	}
	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, Index_ESA<TSpec> >, ESA_LCP>::Type const & getFibre(Index<TText, Index_ESA<TSpec> > const &index, ESA_LCP const) {
		return index.lcp;
	}

	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, Index_ESA<TSpec> >, ESA_LCPE>::Type & getFibre(Index<TText, Index_ESA<TSpec> > &index, ESA_LCPE const) {
		return index.lcpe;
	}
	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, Index_ESA<TSpec> >, ESA_LCPE>::Type const & getFibre(Index<TText, Index_ESA<TSpec> > const &index, ESA_LCPE const) {
		return index.lcpe;
	}

	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, Index_ESA<TSpec> >, ESA_ChildTab>::Type & getFibre(Index<TText, Index_ESA<TSpec> > &index, ESA_ChildTab const) {
		return index.childtab;
	}
	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, Index_ESA<TSpec> >, ESA_ChildTab>::Type const & getFibre(Index<TText, Index_ESA<TSpec> > const &index, ESA_ChildTab const) {
		return index.childtab;
	}

	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, Index_ESA<TSpec> >, ESA_BWT>::Type & getFibre(Index<TText, Index_ESA<TSpec> > &index, ESA_BWT const) {
		return index.bwt;
	}
	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, Index_ESA<TSpec> >, ESA_BWT>::Type const & getFibre(Index<TText, Index_ESA<TSpec> > const &index, ESA_BWT const) {
		return index.bwt;
	}



///.Function.length.param.object.type:Class.Index

	template <typename TText, typename TSpec>
	inline typename Size<Index<TText, Index_ESA<TSpec> > >::Type length(Index<TText, Index_ESA<TSpec> > const &index) {
		return length(getFibre(index, ESA_RawText()));
	}

/**
.Function.rawtextAt:
..summary:Shortcut for $value(indexRawText(..), ..)$.
..cat:Index
..signature:rawtextAt(position, index)
..param.position:A position in the array on which the value should be accessed.
..param.index:The @Class.Index@ object holding the fibre.
...type:Spec.Index_ESA
..returns:A reference or proxy to the value.
*/

	template <typename TSize, typename TIndex>
	inline typename Reference<typename Fibre<TIndex, ESA_RawText>::Type>::Type rawtextAt(TSize i, TIndex &index) {
		return value(getFibre(index, ESA_RawText()), i);
	}
	template <typename TSize, typename TIndex>
	inline typename Reference<typename Fibre<TIndex, ESA_RawText>::Type const>::Type rawtextAt(TSize i, TIndex const &index) {
		return value(getFibre(index, ESA_RawText()), i);
	}

/**
.Function.saAt:
..summary:Shortcut for $value(indexSA(..), ..)$.
..cat:Index
..signature:saAt(position, index)
..param.position:A position in the array on which the value should be accessed.
..param.index:The @Class.Index@ object holding the fibre.
...type:Spec.Index_ESA
..returns:A reference or proxy to the value.
*/

	template <typename TSize, typename TIndex>
	inline typename Reference<typename Fibre<TIndex, ESA_SA>::Type>::Type saAt(TSize i, TIndex &index) {
		return value(getFibre(index, ESA_SA()), i);
	}
	template <typename TSize, typename TIndex>
	inline typename Reference<typename Fibre<TIndex, ESA_SA>::Type const>::Type saAt(TSize i, TIndex const &index) {
		return value(getFibre(index, ESA_SA()), i);
	}

/**
.Function.lcpAt:
..summary:Shortcut for $value(indexLCP(..), ..)$.
..cat:Index
..signature:lcpAt(position, index)
..param.position:A position in the array on which the value should be accessed.
..param.index:The @Class.Index@ object holding the fibre.
...type:Spec.Index_ESA
..returns:A reference or proxy to the value.
*/

	template <typename TSize, typename TIndex>
	inline typename Reference<typename Fibre<TIndex, ESA_LCP>::Type>::Type lcpAt(TSize i, TIndex &index) {
		return value(getFibre(index, ESA_LCP()), i);
	}
	template <typename TSize, typename TIndex>
	inline typename Reference<typename Fibre<TIndex, ESA_LCP>::Type const>::Type lcpAt(TSize i, TIndex const &index) {
		return value(getFibre(index, ESA_LCP()), i);
	}

/**
.Function.lcpeAt:
..summary:Shortcut for $value(indexLCPE(..), ..)$.
..cat:Index
..signature:lcpeAt(position, index)
..param.position:A position in the array on which the value should be accessed.
..param.index:The @Class.Index@ object holding the fibre.
...type:Spec.Index_ESA
..returns:A reference or proxy to the value.
*/

	template <typename TSize, typename TIndex>
	inline typename Reference<typename Fibre<TIndex, ESA_LCPE>::Type>::Type lcpeAt(TSize i, TIndex &index) {
		return value(getFibre(index, ESA_LCPE()), i);
	}
	template <typename TSize, typename TIndex>
	inline typename Reference<typename Fibre<TIndex, ESA_LCPE>::Type const>::Type lcpeAt(TSize i, TIndex const &index) {
		return value(getFibre(index, ESA_LCPE()), i);
	}

/**
.Function.childAt:
..summary:Shortcut for $value(indexChildTab(..), ..)$.
..cat:Index
..signature:childAt(position, index)
..param.position:A position in the array on which the value should be accessed.
..param.index:The @Class.Index@ object holding the fibre.
...type:Spec.Index_ESA
..returns:A reference or proxy to the value.
*/

	template <typename TSize, typename TIndex>
	inline typename Reference<typename Fibre<TIndex, ESA_ChildTab>::Type>::Type childAt(TSize i, TIndex &index) {
		return value(getFibre(index, ESA_ChildTab()), i);
	}
	template <typename TSize, typename TIndex>
	inline typename Reference<typename Fibre<TIndex, ESA_ChildTab>::Type const>::Type childAt(TSize i, TIndex const &index) {
		return value(getFibre(index, ESA_ChildTab()), i);
	}

/**
.Function.bwtAt:
..summary:Shortcut for $value(indexBWT(..), ..)$.
..cat:Index
..signature:bwtAt(position, index)
..param.position:A position in the array on which the value should be accessed.
..param.index:The @Class.Index@ object holding the fibre.
...type:Spec.Index_ESA
..returns:A reference or proxy to the value.
*/

	template <typename TSize, typename TIndex>
	inline typename Reference<typename Fibre<TIndex, ESA_BWT>::Type>::Type bwtAt(TSize i, TIndex &index) {
		return value(getFibre(index, ESA_BWT()).tab, i);
	}
	template <typename TSize, typename TIndex>
	inline typename Reference<typename Fibre<TIndex, ESA_BWT>::Type const>::Type bwtAt(TSize i, TIndex const &index) {
		return value(getFibre(index, ESA_BWT()).tab, i);
	}

	template <typename TSize, typename TIndex>
	inline bool bwtValidAt(TSize i, TIndex const &index) {
		return i != getFibre(index, ESA_BWT()).undefined;
	}

///.Function.clear.param.object.type:Class.Index

	template <typename TValue>
	struct _SizeInval {
		enum	{ VALUE = ~(TValue)0 };
	};

	template <typename TValue>
	inline void _setSizeInval(TValue &v) {
		v = _SizeInval<TValue>::VALUE;
	}

	template <typename TValue>
	inline bool _isSizeInval(TValue const &v) {
		return v == _SizeInval<TValue>::VALUE;
	}

	template <typename TText, typename TSpec>
	inline void clear(Index<TText, Index_ESA<TSpec> > &index) {
		clear(getFibre(index, ESA_SA()));
		clear(getFibre(index, ESA_LCP()));
		clear(getFibre(index, ESA_LCPE()));
		clear(getFibre(index, ESA_ChildTab()));
		clear(getFibre(index, ESA_BWT()).tab);
		_setSizeInval(getFibre(index, ESA_BWT()).undefined);
	}


//////////////////////////////////////////////////////////////////////////////

/**
.Function.indexText:
..summary:Shortcut for $getFibre(.., ESA_Text)$.
..cat:Index
..signature:indexText(index)
..param.index:The @Class.Index@ object holding the fibre.
...type:Spec.Index_ESA
..returns:A reference to the @Tag.ESA_Text@ fibre (original text).
*/

	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, Index_ESA<TSpec> >, ESA_Text>::Type & indexText(Index<TText, Index_ESA<TSpec> > &index) { return getFibre(index, ESA_Text()); }
	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, Index_ESA<TSpec> >, ESA_Text>::Type const & indexText(Index<TText, Index_ESA<TSpec> > const &index) { return getFibre(index, ESA_Text()); }

/**
.Function.indexRawText:
..summary:Shortcut for $getFibre(.., ESA_RawText)$.
..cat:Index
..signature:indexRawText(index)
..param.index:The @Class.Index@ object holding the fibre.
...type:Spec.Index_ESA
..returns:A reference to the @Tag.ESA_RawText@ fibre (concatenated input text).
*/

	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, Index_ESA<TSpec> >, ESA_RawText>::Type & indexRawText(Index<TText, Index_ESA<TSpec> > &index) { return getFibre(index, ESA_RawText()); }
	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, Index_ESA<TSpec> >, ESA_RawText>::Type const & indexRawText(Index<TText, Index_ESA<TSpec> > const &index) { return getFibre(index, ESA_RawText()); }

/**
.Function.indexSA:
..summary:Shortcut for $getFibre(.., ESA_SA)$.
..cat:Index
..signature:indexSA(index)
..param.index:The @Class.Index@ object holding the fibre.
...type:Spec.Index_ESA
..returns:A reference to the @Tag.ESA_SA@ fibre (suffix array).
*/

	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, Index_ESA<TSpec> >, ESA_SA>::Type & indexSA(Index<TText, Index_ESA<TSpec> > &index) { return getFibre(index, ESA_SA()); }
	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, Index_ESA<TSpec> >, ESA_SA>::Type const & indexSA(Index<TText, Index_ESA<TSpec> > const &index) { return getFibre(index, ESA_SA()); }

/**
.Function.indexLCP:
..summary:Shortcut for $getFibre(.., ESA_LCP)$.
..cat:Index
..signature:indexLCP(index)
..param.index:The @Class.Index@ object holding the fibre.
...type:Spec.Index_ESA
..returns:A reference to the @Tag.ESA_LCP@ fibre (lcp table).
*/

	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, Index_ESA<TSpec> >, ESA_LCP>::Type & indexLCP(Index<TText, Index_ESA<TSpec> > &index) { return getFibre(index, ESA_LCP()); }
	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, Index_ESA<TSpec> >, ESA_LCP>::Type const & indexLCP(Index<TText, Index_ESA<TSpec> > const &index) { return getFibre(index, ESA_LCP()); }

/**
.Function.indexLCPE:
..summary:Shortcut for $getFibre(.., ESA_LCPE)$.
..cat:Index
..signature:indexLCPE(index)
..param.index:The @Class.Index@ object holding the fibre.
...type:Spec.Index_ESA
..returns:A reference to the @Tag.ESA_LCPE@ fibre (enhanced lcp table).
*/

	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, Index_ESA<TSpec> >, ESA_LCPE>::Type & indexLCPE(Index<TText, Index_ESA<TSpec> > &index) { return getFibre(index, ESA_LCPE()); }
	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, Index_ESA<TSpec> >, ESA_LCPE>::Type const & indexLCPE(Index<TText, Index_ESA<TSpec> > const &index) { return getFibre(index, ESA_LCPE()); }

/**
.Function.indexBWT:
..summary:Shortcut for $getFibre(.., ESA_BWT)$.
..cat:Index
..signature:indexBWT(index)
..param.index:The @Class.Index@ object holding the fibre.
...type:Spec.Index_ESA
..returns:A reference to the @Tag.ESA_BWT@ fibre (Burrows-Wheeler table).
*/

	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, Index_ESA<TSpec> >, ESA_BWT>::Type & indexBWT(Index<TText, Index_ESA<TSpec> > &index) { return getFibre(index, ESA_BWT()); }
	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, Index_ESA<TSpec> >, ESA_BWT>::Type const & indexBWT(Index<TText, Index_ESA<TSpec> > const &index) { return getFibre(index, ESA_BWT()); }

/**
.Function.indexChildTab:
..summary:Shortcut for $getFibre(.., ESA_ChildTab)$.
..cat:Index
..signature:indexChildTab(index)
..param.index:The @Class.Index@ object holding the fibre.
...type:Spec.Index_ESA
..returns:A reference to the @Tag.ESA_ChildTab@ fibre (child table).
*/

	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, Index_ESA<TSpec> >, ESA_ChildTab>::Type & indexChildTab(Index<TText, Index_ESA<TSpec> > &index) { return getFibre(index, ESA_ChildTab()); }
	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, Index_ESA<TSpec> >, ESA_ChildTab>::Type const & indexChildTab(Index<TText, Index_ESA<TSpec> > const &index) { return getFibre(index, ESA_ChildTab()); }

}

#endif
