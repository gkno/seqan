/*
 *  index_shims.h
 *  SeqAn
 *
 *  Created by David Weese on 17.07.05.
 *
 */

#ifndef SEQAN_HEADER_INDEX_SHIMS_H
#define SEQAN_HEADER_INDEX_SHIMS_H

namespace SEQAN_NAMESPACE_MAIN
{

	//////////////////////////////////////////////////////////////////////////////
	// Suffix Array creation wrappers
	//////////////////////////////////////////////////////////////////////////////

	// build suffix array with an external pipeling algorithm (skew3, skew7, ...)
	template < 
		typename TSA, 
		typename TObject, 
		typename ConstrSpec >
	void createSuffixArrayExt(
		TSA &suffixArray,
		TObject const &text,
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
		StringSet<TString, TSpec> const &stringSet,
		ConstrSpec const &)
	{
        // signed characters behave different than unsigned when compared
        // to get the same index with signed or unsigned chars we simply cast them to unsigned
        // before feeding them into the pipeline
		typedef typename Concatenator<StringSet<TString, TSpec> >::Type			TConcat;
        typedef typename _MakeUnsigned< typename Value<TConcat>::Type >::Type	TUValue;
		typedef Multi<
			ConstrSpec, 
			typename Value<TSA>::Type, 
			typename StringSetLimits<StringSet<TString, TSpec> >::Type >		MultiConstrSpec;

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
			//isSuffixArray(suffixArray, stringSet);
		#endif
	}


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

    template < typename TSA,
               typename TText,
			   typename TAlgSpec >
    inline void createSuffixArray(
		TSA &SA,
		TText &s,
		TAlgSpec const &alg)
	{
		// -> call internal memory algorithm with an extended interface (+ alphabet size, max_depth)
		createSuffixArray(SA, s, alg, ValueSize< typename Value<TText>::Type >::VALUE, 0);
	}

	template < 
		typename TSA, 
		typename TValue, 
		typename TConfig,
		typename TAlgSpec >
	inline void createSuffixArray(
		TSA &SA,
		String< TValue, External<TConfig> > &s,
		TAlgSpec const &alg)
	{
		// -> explicitly create SA using external memory
        createSuffixArrayExt(SA, s, alg);
	}

	template < 
		typename TSA, 
		typename TValue, 
		typename TSpec,
		typename TSSetSpec,
		typename TAlgSpec >
	inline void createSuffixArray(
		TSA &SA,
		StringSet< String<TValue, TSpec>, TSSetSpec > &s,
		TAlgSpec const &)
	{
        createSuffixArrayExt(SA, s, Skew7());
	}


//____________________________________________________________________________


	//////////////////////////////////////////////////////////////////////////////
	// LCP Table creation wrappers
	//////////////////////////////////////////////////////////////////////////////

	template < 
        typename TLCPTable,
		typename TObject, 
        typename TSA,
		typename ConstrSpec >
	void createLCPTableExt(
		TLCPTable &LCP,
		TObject const &text,
		TSA const &suffixArray,
		ConstrSpec const &)
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
		StringSet<TString, TSpec> const &stringSet,
		TSA const &suffixArray,
		ConstrSpec const &)
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


/**
.Function.createLCPTable:
..summary:Creates a lcp table from a given text and suffix array.
..cat:Index
..signature:createLCPTable(lcp, text, suffixArray[, algo_tag])
..param.lcp:A reference to the resulting lcp table.
..param.text:A given text.
..param.suffixArray:The suffix array of $text$.
..param.algo_tag:A tag that identifies the algorithm which is used for creation.
..remarks:The size of $lcp$ must be at least $length(text)$ before calling this function.
*/

	template < 
        typename TLCPTable,
        typename TValue,
		typename TConfig,
		typename TSA >
	void createLCPTable(
		TLCPTable &LCP,
		String< TValue, External<TConfig> > const &s,
		TSA const &SA,
		Kasai const)
	{
        createLCPTableExt(LCP, s, SA, Kasai());
	}

	template < 
        typename TLCPTable,
        typename TValue,
		typename TConfig,
		typename TSSetSpec,
		typename TSA >
	void createLCPTable(
		TLCPTable &LCP,
		StringSet< String<TValue, External<TConfig> >, TSSetSpec > const &s,
		TSA const &SA,
		Kasai const)
	{
        createLCPTableExt(LCP, s, SA, Kasai());
	}


//____________________________________________________________________________


	//////////////////////////////////////////////////////////////////////////////
	// Enhanced LCP Table creation wrappers
	//////////////////////////////////////////////////////////////////////////////

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

    // build enhanced LCP table with an lcp algorithm
	// and a dynamic programming tree construction alg
    template <
        typename TValue,
        typename TSpec,
        typename TText,
        typename TSA,
		typename ConstrSpec >
    void createLCPETable(
		String< TValue, TSpec > &LCPE,
		TText &s,
		TSA &SA,
		ConstrSpec const &spec)
	{
        //TSA LCP;
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
        typename TSA,
		typename ConstrSpec >
    void createLCPETable(
		String< TValue, External<TConfig> > &LCPE,
		TText &s,
		TSA &SA,
		ConstrSpec const &spec)
	{
        createLCPETableExt(LCPE, s, SA, spec);
    }

    template <
        typename TValue,
        typename TSpec,
        typename TText,
        typename TSA>
    inline void createLCPETable(
		String< TValue, TSpec > &LCPE,
		TText &s,
		TSA &SA)
	{
		CreateLCPETable(LCPE, s, SA, Kasai());
    }

/*
	// build enhanced LCP+SuffixArray table with an lcp algorithm
	// and a dynamic programming tree construction alg

    template <
        typename TValue,
        typename TSpec,
        typename TText,
        typename TSA,
		typename ConstrSpec >
    void createLCPHTable(
		String< TValue, TSpec > &LCPH,
		TText &s,
		TSA &SA,
		ConstrSpec const &spec = Kasai())
	{
        TSA LCP;
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
        typename TSA,
		typename ConstrSpec >
    void createLCPHTable(
		String< TValue, External<TConfig> > &LCPH,
		TText &s,
		TSA &SA,
		ConstrSpec const &spec)
	{
        createLCPHTableExt(LCPH, s, SA, spec);
    }

	template < 
		typename TValue, 
        typename TSpec,
        typename TText,
        typename TSA >
	void createLCPHTable(
		String< TValue, TSpec > &LCPH,
		TText &s,
		TSA &SA)
	{
		createLCPHTable(LCPH, s, SA, Kasai());
	}
*/

//////////////////////////////////////////////////////////////////////////////

	template <typename TOccValue>
	struct _SAValueLess:
		public ::std::less<TOccValue> {};

	template <typename T1, typename T2, typename TCompression>
	struct _SAValueLess< Pair<T1,T2,TCompression> >:
		public ::std::binary_function< Pair<T1,T2,TCompression>, Pair<T1,T2,TCompression>, bool> 
	{
		inline bool operator()(const Pair<T1,T2,TCompression> &a, const Pair<T1,T2,TCompression> &b) const {
			return	getValueI1(a) < getValueI1(b) ||
					getValueI1(a) == getValueI1(b) && getValueI2(a) < getValueI2(b);
		}
	};

	template <typename TValue, typename TSpec>
	inline void orderOccurences(String<TValue, TSpec> &occString)
	{
		::std::sort(begin(occString, Standard()), end(occString, Standard()), _SAValueLess<TValue>());
	}


//////////////////////////////////////////////////////////////////////////////
// table creators

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

	template <typename TText, typename TSpec, typename TSpecAlg>
	inline bool indexCreate(Index<TText, TSpec> &index, Tag<_Fibre_SA> const, TSpecAlg const alg) {
		resize(indexSA(index), length(indexRawText(index)), Exact());
		createSuffixArray(indexSA(index), indexText(index), alg);
		return true;
	}

	template <typename TText, typename TSpec, typename TSpecAlg>
	inline bool indexCreate(Index<TText, TSpec> &index, Tag<_Fibre_LCP> const, TSpecAlg const alg) {
		resize(indexLCP(index), length(indexRawText(index)), Exact());
		createLCPTable(indexLCP(index), indexText(index), indexSA(index), alg);
		return true;
	}

	template <typename TText, typename TSpec, typename TSpecAlg>
	inline bool indexCreate(Index<TText, TSpec> &index, Tag<_Fibre_LCPE> const, TSpecAlg const alg) {
	//TODO: separate LCP from LCPE (for now LCPE = LCP + extra)
//		resize(indexLCP(index), length(indexRawText(index)), Exact());
//		createLCPETable(indexLCPE(index), indexRawText(index), indexSA(index), alg);
		return false;
	}

	template <typename TText, typename TSpec>
	inline bool indexCreate(Index<TText, TSpec> &index, Tag<_Fibre_BWT> const, BWT const) {
		resize(indexBWT(index), length(indexRawText(index)), Exact());
		createBWTable(indexBWT(index), indexRawText(index), indexRawSA(index));
		return true;
	}

	template <typename TText, typename TSpec>
	inline bool indexCreate(Index<TText, TSpec> &index, Tag<_Fibre_ChildTab> const, ChildTab const) {
		resize(indexChildTab(index), length(indexRawText(index)), Exact());
		createChildTable(indexChildTab(index), indexLCP(index));
		return true;
	}

	template <typename TText, typename TSpec, typename TFibre>
	inline bool indexCreate(Index<TText, TSpec> &index, TFibre const fibre) {
		return indexCreate(index, fibre, typename DefaultIndexCreator<Index<TText, TSpec>, TFibre>::Type());
	}

//////////////////////////////////////////////////////////////////////////////
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
	inline bool indexSupplied(Index<TText, TSpec> &index, TFibre const fibre) {
		return !empty(getFibre(index, fibre));
	}


//////////////////////////////////////////////////////////////////////////////
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
	inline bool indexRequire(Index<TText, TSpec> &index, TFibre const fibre) {
		if (indexSupplied(index, fibre)) return true;				// if the table doesn't exist,
		if (!indexSolveDependencies(index, fibre)) return false;	// fulfill requirements
		return indexCreate(index, fibre);							// and create table
	}


//////////////////////////////////////////////////////////////////////////////
// solve dependencies

	template <typename TText, typename TSpec, typename TFibre>
	inline bool indexSolveDependencies(Index<TText, TSpec> &index, TFibre const fibre) {
		return true;
	}

	template <typename TText, typename TSpec, typename TFibre>
	inline bool indexSolveDependencies(Index<TText, TSpec> &index, Tag<_Fibre_LCP>) {
		return indexRequire(index, Tag<_Fibre_SA>());
	}

	template <typename TText, typename TSpec, typename TFibre>
	inline bool indexSolveDependencies(Index<TText, TSpec> &index, Tag<_Fibre_LCPE>) {
		return indexRequire(index, Tag<_Fibre_LCP>());
	}

	template <typename TText, typename TSpec, typename TFibre>
	inline bool indexSolveDependencies(Index<TText, TSpec> &index, Tag<_Fibre_ChildTab>) {
		return indexRequire(index, Tag<_Fibre_LCP>());
	}

	template <typename TText, typename TSpec, typename TFibre>
	inline bool indexSolveDependencies(Index<TText, TSpec> &index, Tag<_Fibre_BWT>) {
		return indexRequire(index, Tag<_Fibre_SA>());
	}


//////////////////////////////////////////////////////////////////////////////
// open

	template < typename TValue, typename TSpec >
	inline bool open(String<TValue, TSpec> &string, const char *fileName, int openMode) {
		String<TValue, External< ExternalConfigLarge<> > > extString;
		if (!open(extString, fileName, openMode)) return false;
		assign(string, extString);
		return true;
	}
	template < typename TValue, typename TSpec >
	inline bool open(String<TValue, TSpec> &string, const char *fileName) {
		return open(string, fileName, OPEN_RDONLY);
	}

	template < typename THost, typename TSpec >
	inline bool open(Segment<THost, TSpec> &string, const char *fileName, int openMode) {
		String<typename Value<THost>::Type, External< ExternalConfigLarge<> > > extString;
		if (!open(extString, fileName, openMode)) return false;
		assign(string, extString);
		return true;
	}
	template < typename THost, typename TSpec >
	inline bool open(Segment<THost, TSpec> &string, const char *fileName) {
		return open(string, fileName, OPEN_RDONLY);
	}

	template < typename TValue, typename TSpec, typename TSSSpec >
	inline bool open(StringSet<String<TValue, TSpec>, TSSSpec> &multi, const char *fileName, int openMode) {
		bool more = true;
		char id[11]; // 2^32 has 10 decimal digits + 1 (0x00)
		unsigned i = 0;
		for(; more; ++i) {
			sprintf(id, "%u", i);
			String<char> name;
			name = fileName;	append(name, id);
			more = open(value(multi, i), toCString(name), openMode);
		}
		return i > 1;
	}
	template < typename TValue, typename TSpec, typename TSSSpec>
	inline bool open(StringSet<String<TValue, TSpec>, TSSSpec> &multi, const char *fileName) {
		return open(multi, fileName, OPEN_RDONLY);
	}


//////////////////////////////////////////////////////////////////////////////
// save

	template < typename TValue, typename TSpec >
	inline bool save(String<TValue, TSpec> const &string, const char *fileName, int openMode) {
		if (length(string) == 0) return true;
		String<TValue, External< ExternalConfigLarge<> > > extString;
		if (!open(extString, fileName, openMode)) return false;
		assign(extString, string);
		return true;
	}
	template < typename TValue, typename TSpec >
	inline bool save(String<TValue, TSpec> &string, const char *fileName) {
		return save(string, fileName, OPEN_WRONLY | OPEN_CREATE);
	}

	template < typename THost, typename TSpec >
	inline bool save(Segment<THost, TSpec> const &string, const char *fileName, int openMode) {
		if (length(string) == 0) return true;
		String<typename Value<THost>::Type, External< ExternalConfigLarge<> > > extString;
		if (!open(extString, fileName, openMode)) return false;
		assign(extString, string);
		return true;
	}
	template < typename THost, typename TSpec >
	inline bool save(Segment<THost, TSpec> const &string, const char *fileName) {
		return save(string, fileName, OPEN_WRONLY | OPEN_CREATE);
	}

	template < typename TValue, typename TSpec, typename TSSSpec>
	inline bool save(StringSet<String<TValue, TSpec>, TSSSpec> const &multi, const char *fileName, int openMode) {
		if (length(multi) == 0) return true;
		bool result = true;
		char id[12]; // 2^32 has 10 decimal digits + 2 ('.' and 0x00)
		for(unsigned i = 0; i < length(multi); ++i) {
			sprintf(id, ".%u", i);
			String<char> name;
			name = fileName;	
			append(name, &(id[0]));
			if (!save(getValue(multi, i), toCString(name), openMode))
				result = false;
		}
		return result;
	}
	template < typename TValue, typename TSpec, typename TSSSpec>
	inline bool save(StringSet<String<TValue, TSpec>, TSSSpec> const &multi, const char *fileName) {
		return save(multi, fileName, OPEN_WRONLY | OPEN_CREATE);
	}

}

#endif
