/*
 *  index_shims.h
 *  genindex
 *
 *  Created by David Weese on 17.07.05.
 *
 */

#ifndef SEQAN_HEADER_INDEX_SHIMS_H
#define SEQAN_HEADER_INDEX_SHIMS_H

namespace SEQAN_NAMESPACE_MAIN
{

	//////////////////////////////////////////////////////////////////////////////
	// explicit external creation (used as a wrapper-end by algorithms)
	//////////////////////////////////////////////////////////////////////////////

/**
.Function.createSuffixArray:
..summary:Creates a suffix array from a given text.
..cat:Index
..signature:createSuffixArray(suffixArray, text[, algo_tag])
..param.suffixArray:A reference to the resulting suffix array.
..param.text:A given text.
..param.algo_tag:A tag that identifies the algorithm which is used for creation.
..remarks:The size of $suffixArray$ must be at least $length(text)$ before calling this function.
*/

	// build suffix array with an external pipeling algorithm (skew3, skew7, ...)
	template < 
		typename TSA, 
		typename TObject, 
		typename ConstrSpec >
	void createSuffixArrayExt(
		TSA &suffixArray,
		TObject &text,
		ConstrSpec const &)
	{
        // signed characters behave different than unsigned when compared
        // to get the same index with signed or unsigned chars we simply cast them to unsigned
        // before feeding them into the pipeline
        typedef typename _MakeUnsigned< typename Value<TObject>::Type >::Type TUValue;

        // specialization
		typedef Pipe< TObject, Source<> >				src_t;
        typedef Pipe< src_t, Caster<TUValue> >          unsigner_t;
		typedef Pipe< unsigner_t, ConstrSpec >	        creator_t;

		// instantiation
		src_t		src(text);
        unsigner_t  unsigner(src);
		creator_t	creator;

		// processing
		creator << unsigner;
		suffixArray << creator;
		#ifdef SEQAN_TEST_INDEX
			isSuffixArray(suffixArray, text);
		#endif
	}


	// build suffix array (external) for mutliple sequences
	template < 
		typename TSA, 
		typename TString, 
		typename TSpec,
		typename ConstrSpec >
	void createSuffixArrayExt(
		TSA &suffixArray,
		StringSet<TString, TSpec> &stringSet,
		ConstrSpec const &)
	{
        // signed characters behave different than unsigned when compared
        // to get the same index with signed or unsigned chars we simply cast them to unsigned
        // before feeding them into the pipeline
		typedef typename Concatenator<StringSet<TString, TSpec> >::Type TConcat;
        typedef typename _MakeUnsigned< typename Value<TConcat>::Type >::Type TUValue;
		typedef Multi<
			ConstrSpec, 
			typename Value<TSA>::Type, 
			typename StringSetLimits<StringSet<TString, TSpec> >::Type > MultiConstrSpec;

        // specialization
		typedef Pipe< TConcat, Source<> >				src_t;
        typedef Pipe< src_t, Caster<TUValue> >          unsigner_t;
		typedef Pipe< unsigner_t, MultiConstrSpec >	    creator_t;

		// instantiation
		src_t		src(concat(stringSet));
        unsigner_t  unsigner(src);
		creator_t	creator(stringSetLimits(stringSet));

		// processing
		creator << unsigner;
		suffixArray << creator;
		#ifdef SEQAN_TEST_INDEX
			isSuffixArray(suffixArray, stringSet);
		#endif
	}


	// build lcp table with an external pipelining algorithm (ext kasai, ...)
	template < 
        typename TLCPTable,
		typename TObject, 
        typename TSA,
		typename ConstrSpec >
	void createLCPTableExt(
		TLCPTable &LCP,
		TObject &text,
		TSA &suffixArray,
		ConstrSpec const &spec = Kasai())
	{
		// specialization
		typedef Pipe< TObject, Source<> >							srcText_t;
		typedef Pipe< TSA, Source<> >   							srcSA_t;
	    typedef Pipe< Bundle2< srcText_t, srcSA_t >, ConstrSpec >	creator_t;

		// instantiation
		srcText_t	srcText(text);
		srcSA_t		srcSA(suffixArray);
		creator_t	creator;

		// processing
	    creator << bundle2(srcText, srcSA);
		LCP << creator;
		#ifdef SEQAN_TEST_INDEX
			isLCPTable(LCP, suffixArray, text);
		#endif
	}


	// build lcp table (external) for mutliple sequences
	template < 
        typename TLCPTable,
		typename TString,
		typename TSpec,
        typename TSA,
		typename ConstrSpec >
	void createLCPTableExt(
		TLCPTable &LCP,
		StringSet<TString, TSpec> &stringSet,
		TSA &suffixArray,
		ConstrSpec const &spec = Kasai())
	{
		typedef typename Concatenator<StringSet<TString, TSpec> >::Type TConcat;
		typedef Multi<
			ConstrSpec, 
			typename Value<TSA>::Type, 
			typename StringSetLimits<StringSet<TString, TSpec> >::Type > MultiConstrSpec;

		// specialization
		typedef Pipe< TConcat, Source<> >								srcText_t;
		typedef Pipe< TSA, Source<> >   								srcSA_t;
	    typedef Pipe< Bundle2< srcText_t, srcSA_t >, MultiConstrSpec >	creator_t;

		// instantiation
		srcText_t	srcText(concat(stringSet));
		srcSA_t		srcSA(suffixArray);
		creator_t	creator(stringSetLimits(stringSet));

		// processing
	    creator << bundle2(srcText, srcSA);
		LCP << creator;
		#ifdef SEQAN_TEST_INDEX
			isLCPTable(LCP, suffixArray, text);
		#endif
	}


    // build enhanced LCP table with an external pipelining algorithm (ext kasai, ...)
	// and a dynamic programming tree construction alg
	template < 
		typename TValue, 
        typename TSpec,
		typename TObject, 
        typename TSA,
		typename ConstrSpec >
	void createLCPETableExt(
		String< TValue, TSpec > &LCPE,
		TObject &text,
		TSA &suffixArray,
		ConstrSpec const &)
	{
		// specialization
		typedef Pipe< TObject, Source<> >						    srcText_t;
		typedef Pipe< TSA, Source<> >					    	    srcSA_t;
	    typedef Pipe< Bundle2< srcText_t, srcSA_t >, ConstrSpec >	creator_t;

		// instantiation
		srcText_t	srcText(text);
		srcSA_t		srcSA(suffixArray);
		creator_t	creator;

		// processing
	    creator << bundle2(srcText, srcSA);
		#ifdef SEQAN_TEST_INDEX
			isLCPTable(creator, suffixArray, text);
		#endif
		createLCPBinTree(LCPE, creator);
	}

	// build enhanced LCP+SuffixArray table with an external version of Kasai
	// and a dynamic programming tree construction alg
	template < 
		typename TValue, 
        typename TSpec,
		typename TObject, 
        typename TSA,
		typename ConstrSpec >
	void createLCPHTableExt(
		String< TValue, TSpec > &LCPH,
		TObject &text,
		TSA &suffixArray,
		ConstrSpec const &spec) 
	{
		// specialization
		typedef Pipe< TObject, Source<> >					        srcText_t;
		typedef Pipe< TSA, Source<> >								srcSA_t;
	    typedef Pipe< Bundle2< srcText_t, srcSA_t >, ConstrSpec >	creator_t;

		// instantiation
		srcText_t	srcText(text);
		srcSA_t		srcSA(suffixArray);
		creator_t	creator;

		// processing
	    creator << bundle2(srcText, srcSA);
		#ifdef SEQAN_TEST_INDEX
			isLCPTable(creator, srcSA, text);
		#endif
		createHybridBinTree(LCPH, creator, suffixArray);
	}


	//////////////////////////////////////////////////////////////////////////////
	// various enhanced table creators
	//////////////////////////////////////////////////////////////////////////////

    // build enhanced LCP table with an lcp algorithm
	// and a dynamic programming tree construction alg
    template <
        typename TValue,
        typename TSpec,
        typename TText,
        typename TSuffixArray,
		typename ConstrSpec >
    void createLCPETable(
		String< TValue, TSpec > &LCPE,
		TText &s,
		TSuffixArray &SA,
		ConstrSpec const &spec)
	{
        //TSuffixArray LCP;
        //resize(LCP, length(s), Exact());
		// we use LCPE[n-lcpSize..n-1] as a temporary buffer instead of allocating one
		typename Suffix<String< TValue, TSpec > >::Type LCP = suffix(LCPE, length(LCPE) - length(s));

		createLCPTable(LCP, s, SA, spec);
		#ifdef SEQAN_TEST_INDEX
			isLCPTable(LCP, SA, s);
		#endif
        createLCPBinTree(LCPE, LCP);
    }

    template <
        typename TValue,
        typename TConfig,
        typename TText,
        typename TSuffixArray,
		typename ConstrSpec >
    void createLCPETable(
		String< TValue, External<TConfig> > &LCPE,
		TText &s,
		TSuffixArray &SA,
		ConstrSpec const &spec)
	{
        createLCPETableExt(LCPE, s, SA, spec);
    }

    template <
        typename TValue,
        typename TSpec,
        typename TText,
        typename TSuffixArray>
    inline void createLCPETable(
		String< TValue, TSpec > &LCPE,
		TText &s,
		TSuffixArray &SA)
	{
		CreateLCPETable(LCPE, s, SA, Kasai());
    }


	// build enhanced LCP+SuffixArray table with an lcp algorithm
	// and a dynamic programming tree construction alg

    template <
        typename TValue,
        typename TSpec,
        typename TText,
        typename TSuffixArray,
		typename ConstrSpec >
    void createLCPHTable(
		String< TValue, TSpec > &LCPH,
		TText &s,
		TSuffixArray &SA,
		ConstrSpec const &spec = Kasai())
	{
        TSuffixArray LCP;
        resize(LCP, length(s), Exact());
        createLCPTable(LCP, s, SA, spec);
		#ifdef SEQAN_TEST_INDEX
			isLCPTable(LCP, SA, s);
		#endif
        createHybridBinTree(LCPH, LCP, SA);
    }

    template <
        typename TValue,
        typename TConfig,
        typename TText,
        typename TSuffixArray,
		typename ConstrSpec >
    void createLCPHTable(
		String< TValue, External<TConfig> > &LCPH,
		TText &s,
		TSuffixArray &SA,
		ConstrSpec const &spec)
	{
        createLCPHTableExt(LCPH, s, SA, spec);
    }

	template < 
		typename TValue, 
        typename TSpec,
        typename TText,
        typename TSuffixArray >
	void createLCPHTable(
		String< TValue, TSpec > &LCPH,
		TText &s,
		TSuffixArray &SA)
	{
		createLCPHTable(LCPH, s, SA, Kasai());
	}


	// ESA finders

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

//____________________________________________________________________________

	template < typename TFinder, typename TPattern >
	inline void _findFirstESAIndex(
		TFinder &finder,
		TPattern const &pattern,
		ESA_MLR const)
	{
		typename Haystack<TFinder>::Type &index = haystack(finder);
		indexRequire(index, ESA_SA());
		finder.range = equalRangeSA(indexRawText(index), indexSA(index), pattern);
	}

	template < typename TFinder, typename TPattern >
	inline void _findFirstESAIndex(
		TFinder &finder,
		TPattern const &pattern,
		ESA_LCPE const)
	{
		typename Haystack<TFinder>::Type &index = haystack(finder);
		indexRequire(index, ESA_SA());
		indexRequire(index, ESA_LCPE());
		finder.range = equalRangeLCPE(indexRawText(index), indexSA(index), indexLCPE(index), pattern);
	}

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
	// index creation

    template <typename TObject>
    inline void setKeepFirst(TObject &obj, bool b) {}

	template <typename TValue, typename TConfig>
    inline void setKeepFirst(String<TValue, External<TConfig> > &obj, bool b) { obj.keepFirst = b; }



/**
.Function.indexCreate:
..summary:Creates a specific @Metafunction.Fibre@.
..cat:Index
..signature:indexCreate(index, fibre_tag[, algo_tag])
..param.index:The @Class.Index@ object holding the fibre.
...type:Class.Index
..param.fibre_tag:A tag that identifies the @Metafunction.Fibre@ (e.g. @Tag.ESA_SA@).
..param.algo_tag:A tag that identifies the algorithm which is used to create the fibre.
...default:The result of @Metafunction.DefaultIndexCreator@.
..returns:A $bool$ which is $true$ on a successful creation.
..remarks:$indexCreate$ calls the fibre corresponding $createXXX(..)$ function (e.g. @Function.createSuffixArray@).
*/

	template <typename TText, typename TSpec, typename TFibre>
	inline bool indexCreate(Index<TText, Index_ESA<TSpec> > &index, TFibre const fibre) {
		return indexCreate(index, fibre, typename DefaultIndexCreator<Index<TText, Index_ESA<TSpec> >, TFibre>::Type());
	}

	// table creators

	template <typename TText, typename TSpec, typename TSpecAlg>
	inline bool indexCreate(Index<TText, Index_ESA<TSpec> > &index, ESA_SA const, TSpecAlg const alg) {
		resize(indexSA(index), length(indexRawText(index)), Exact());
		createSuffixArray(indexSA(index), indexText(index), alg);
		return true;
	}

	template <typename TText, typename TSpec, typename TSpecAlg>
	inline bool indexCreate(Index<TText, Index_ESA<TSpec> > &index, ESA_LCP const, TSpecAlg const alg) {
		resize(indexLCP(index), length(indexRawText(index)), Exact());
		createLCPTable(indexLCP(index), indexText(index), indexSA(index), alg);
		return true;
	}

	template <typename TText, typename TSpec, typename TSpecAlg>
	inline bool indexCreate(Index<TText, Index_ESA<TSpec> > &index, ESA_LCPE const, TSpecAlg const alg) {
	//TODO: separate LCP from LCPE (for now LCPE = LCP + extra)
//		resize(indexLCP(index), length(indexRawText(index)), Exact());
//		createLCPETable(indexLCPE(index), indexRawText(index), indexSA(index), alg);
		return false;
	}

	template <typename TText, typename TSpec>
	inline bool indexCreate(Index<TText, Index_ESA<TSpec> > &index, ESA_BWT const, BWT const) {
		resize(indexBWT(index), length(indexRawText(index)), Exact());
		createBWTable(indexBWT(index), indexRawText(index), indexRawSA(index));
		return true;
	}

	template <typename TText, typename TSpec>
	inline bool indexCreate(Index<TText, Index_ESA<TSpec> > &index, ESA_ChildTab const, ChildTab const) {
		resize(indexChildTab(index), length(indexRawText(index)), Exact());
		createChildTable(indexChildTab(index), indexLCP(index));
		return true;
	}

	// automatic creation

/**
.Function.indexSupplied:
..summary:Returns whether a specific @Metafunction.Fibre@ is present.
..cat:Index
..signature:indexSupplied(index, fibre_tag)
..param.index:The @Class.Index@ object holding the fibre.
...type:Class.Index
..param.fibre_tag:A tag that identifies the @Metafunction.Fibre@ (e.g. @Tag.ESA_SA@).
..returns:A $bool$ which is $true$, iff the fibre is present.
*/

	template <typename TText, typename TSpec, typename TFibre>
	inline bool indexSupplied(Index<TText, Index_ESA<TSpec> > &index, TFibre const fibre) {
		return !empty(getFibre(index, fibre));
	}

	template <typename TText, typename TSpec>
	inline bool indexSupplied(Index<TText, Index_ESA<TSpec> > &index, ESA_BWT) {
		return !empty(getFibre(index, ESA_BWT()));
	}

/**
.Function.indexRequire:
..summary:On-demand creation of a specific @Metafunction.Fibre@.
..cat:Index
..signature:indexRequire(index, fibre_tag)
..param.index:The @Class.Index@ object holding the fibre.
...type:Class.Index
..param.fibre_tag:A tag that identifies the @Metafunction.Fibre@ (e.g. @Tag.ESA_SA@).
..returns:A $bool$ which is $true$ on a successful creation.
..remarks:If the fibre already exists (@Function.indexSupplied@ is true) then $indexRequire$ does nothing.
If the fibre doesn't exist then @Function.indexCreate@ is called to create it.
*/

	template <typename TText, typename TSpec, typename TFibre>
	inline bool indexRequire(Index<TText, Index_ESA<TSpec> > &index, TFibre const fibre) {
		if (indexSupplied(index, fibre)) return true;				// if the table doesn't exist,
		if (!indexSolveDependencies(index, fibre)) return false;	// fulfill requirements
		return indexCreate(index, fibre);							// and create table
	}

	// dependencies

	template <typename TText, typename TSpec, typename TFibre>
	inline bool indexSolveDependencies(Index<TText, Index_ESA<TSpec> > &index, TFibre const fibre) {
		return true;
	}

	template <typename TText, typename TSpec, typename TFibre>
	inline bool indexSolveDependencies(Index<TText, Index_ESA<TSpec> > &index, ESA_LCP) {
		return indexRequire(index, ESA_SA());
	}

	template <typename TText, typename TSpec, typename TFibre>
	inline bool indexSolveDependencies(Index<TText, Index_ESA<TSpec> > &index, ESA_LCPE) {
		return indexRequire(index, ESA_LCP());
	}

	template <typename TText, typename TSpec, typename TFibre>
	inline bool indexSolveDependencies(Index<TText, Index_ESA<TSpec> > &index, ESA_ChildTab) {
		return indexRequire(index, ESA_LCP());
	}

	template <typename TText, typename TSpec, typename TFibre>
	inline bool indexSolveDependencies(Index<TText, Index_ESA<TSpec> > &index, ESA_BWT) {
		return indexRequire(index, ESA_SA());
	}

//____________________________________________________________________________

	template < typename TValue, typename TSpec >
	inline bool open(String<TValue, TSpec> &string, const char *fileName, int dummy = OPEN_RDONLY) {
		String<TValue, External<> > extString;
		if (!open(extString, fileName, OPEN_RDONLY)) return false;
		string = extString;
		return true;
	}

	template < typename TValue, typename TSpec >
	inline bool open(StringSet<String<TValue, TSpec> > &multi, const char *fileName, int dummy = OPEN_RDONLY) {
		String<TValue, External<> > extString;
		bool more = true;
		char id[11]; // 2^32 has 10 decimal digits + 1 (0x00)
		unsigned i = 0;
		for(; more; ++i) {
			sprintf(id, "%u", i);
			String<char> name;
			name = fileName;	append(name, id);
			if (more = open(extString, toCString(name), OPEN_RDONLY)) {
				getValue(multi, i) = extString;
				close(extString);
			}
		}
		return i > 1;
	}

	template < typename TValue, typename TSpec >
	inline bool save(String<TValue, TSpec> &string, const char *fileName, int dummy = OPEN_WRONLY | OPEN_CREATE) {
		if (length(string) == 0) return true;
		String<TValue, External<> > extString;
		if (!open(extString, fileName, OPEN_WRONLY | OPEN_CREATE)) return false;
		extString = string;
		return true;
	}

	template < typename TValue, typename TSpec, typename TSetSpec >
	inline bool save(StringSet<String<TValue, TSpec>, TSetSpec > &multi, const char *fileName, int dummy = OPEN_WRONLY | OPEN_CREATE) {
		if (length(multi) == 0) return true;
		String<TValue, External<> > extString;
		bool result = true;
		char id[12]; // 2^32 has 10 decimal digits + 2 ('.' and 0x00)
		for(unsigned i = 0; i < length(multi); ++i) {
			sprintf(id, ".%u", i);
			String<char> name;
			name = fileName;	
			append(name, &(id[0]));
			if (open(extString, toCString(name), OPEN_WRONLY | OPEN_CREATE)) {
				extString = getValue(multi, i);
				close(extString);
			} else
				result = false;
		}
		return result;
	}

	template < typename TObject, typename TSpec >
	inline bool open(
		Index< TObject, Index_ESA<TSpec> > &index, 
		const char *fileName,
		int openMode = OPEN_RDONLY)
	{
		String<char> name;
		name = fileName;	append(name, ".txt");	
		if (!open(getFibre(index, ESA_Text()), toCString(name), openMode) && 
			!open(getFibre(index, ESA_Text()), fileName), openMode) return false;

		name = fileName;	append(name, ".sa");	open(getFibre(index, ESA_SA()), toCString(name), openMode);
		name = fileName;	append(name, ".lcp");	open(getFibre(index, ESA_LCP()), toCString(name), openMode);
		name = fileName;	append(name, ".child");	open(getFibre(index, ESA_ChildTab()), toCString(name), openMode);
		name = fileName;	append(name, ".bwt");	open(getFibre(index, ESA_BWT()), toCString(name), openMode);
		return true;
	}

	template < typename TObject, typename TSpec >
	inline bool save(
		Index< TObject, Index_ESA<TSpec> > &index, 
		const char *fileName,
		int openMode = OPEN_RDONLY)
	{
		String<char> name;
		name = fileName;	append(name, ".txt");	
		if (!save(getFibre(index, ESA_Text()), toCString(name)) && 
			!save(getFibre(index, ESA_Text()), fileName)) return false;

		name = fileName;	append(name, ".sa");	save(getFibre(index, ESA_SA()), toCString(name));
		name = fileName;	append(name, ".lcp");	save(getFibre(index, ESA_LCP()), toCString(name));
		name = fileName;	append(name, ".child");	save(getFibre(index, ESA_ChildTab()), toCString(name));
		name = fileName;	append(name, ".bwt");	save(getFibre(index, ESA_BWT()), toCString(name));
		return true;
	}

}

#endif
