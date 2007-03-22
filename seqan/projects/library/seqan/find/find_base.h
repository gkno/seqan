#ifndef SEQAN_HEADER_FIND_BASE_H
#define SEQAN_HEADER_FIND_BASE_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////

/**
.Metafunction.DefaultFinder:
..cat:Pattern Matching
..summary:Default @Class.Finder@ specialization type.
..signature:DefaultFinder<THaystack>::Type
..param.THaystack:The given haystack type.
..returns:Is $void$ by default and @Tag.ESA_FIND_MLR@ if $THaystack$ is an @Class.Index@.
*/
	template < typename TObject >
	struct DefaultFinder {
		typedef void Type;
	};

/**
.Metafunction.DefaultPattern:
..cat:Pattern Matching
..summary:Default @Class.Pattern@ specialization type.
..signature:DefaultPattern<TNeedle>::Type
..param.TNeedle:The given needle type.
..returns:Is $void$ by default.
*/
	template < typename TObject >
	struct DefaultPattern {
		typedef void Type;
	};

/**
.Metafunction.Haystack:
..summary:Returns the haystack type of a @Class.Finder@ type.
..cat:Pattern Matching
..signature:Haystack<TFinder>::Type
..param.TFinder:A @Class.Finder@ type.
...type:Class.Finder
..returns:The haystack type of $TFinder$, i.e. $THaystack$ for $Finder<THaystack, TSpec>$.
*/

	template <typename TFinder>
	struct Haystack {
		typedef typename Container<TFinder>::Type Type;
	};

/**
.Metafunction.Needle:
..summary:Returns the needle type of a @Class.Pattern@ type.
..cat:Pattern Matching
..signature:Needle<TPattern>::Type
..param.TPattern:A @Class.Pattern@ type.
...type:Class.Pattern
..returns:The needle type of $TPattern$, i.e. $TNeedle$ for $Pattern<TNeedle, TSpec>$.
*/

	template <typename TPattern>
	struct Needle {
		typedef typename Host<TPattern>::Type Type;
	};

//////////////////////////////////////////////////////////////////////////////

/**
.Class.Pattern:
..summary:Holds the needle and preprocessing data (depends on algorithm).
..cat:Pattern Matching
..signature:Pattern<TNeedle[, TSpec]>
..param.TNeedle:The needle type.
...type:Class.String
..param.TSpec:The online-algorithm to search with.
...remarks:Leave empty for index-based pattern matching (see @Class.Index@).
...default:The result of @Metafunction.DefaultPattern@
..remarks:If $TNeedle$ is a set of strings, then $position(pattern)$ returns the index of the currently matching needle.
*/

	template < typename TNeedle, typename TSpec = typename DefaultPattern<TNeedle>::Type >
	class Pattern
	{
	private:
		Holder<TNeedle> data_host;

	public:

		Pattern() {}

		template <typename _TNeedle>
		Pattern(_TNeedle & ndl):
			data_host(ndl) {}

		template <typename _TNeedle>
		Pattern(_TNeedle const & ndl):
			data_host(ndl) {}

		friend inline typename Host<Pattern>::Type & 
		host(Pattern &me) { 
			return value(me.data_host);
		}

		friend inline typename Host<Pattern>::Type const & 
		host(Pattern const & me) {
			return value(me.data_host);
		}

		friend inline void
		setHost(Pattern &me, TNeedle const & ndl) {
			me.data_host = ndl;
		}

		friend inline void
		setHost(Pattern &me, TNeedle & ndl) {
			me.data_host = ndl;
		}
	};

	template < typename TObject >
	inline typename Needle<TObject>::Type &
	needle(TObject &obj) {
		return obj;
	}

	template < typename TObject >
	inline typename Needle<TObject const>::Type &
	needle(TObject const &obj) {
		return obj;
	}

/**
.Function.needle:
..summary:Returns the needle of a @Class.Pattern@ object (not implemented for some online-algorithms).
..cat:Pattern Matching
..signature:needle(pattern)
..param.pattern:The @Class.Pattern@ object to search with.
...type:Class.Pattern
..returns:The needle object to search for.
..remarks:The result type is @Metafunction.Needle@$<TPattern>::Type$ for pattern of type $TPattern$.
*/

///.Function.position.param.iterator.type:Class.Pattern

	template < typename TNeedle, typename TSpec >
	inline typename Needle< Pattern<TNeedle, TSpec> >::Type &
	needle(Pattern<TNeedle, TSpec> &obj) {
		return host(obj);
	}

	template < typename TNeedle, typename TSpec >
	inline typename Needle< Pattern<TNeedle, TSpec> const>::Type &
	needle(Pattern<TNeedle, TSpec> const &obj) {
		return host(obj);
	}

/**
.Function.setNeedle:
..summary:Sets the needle of a @Class.Pattern@ object and optionally induces preprocessing.
..cat:Pattern Matching
..signature:setNeedle(pattern, needle)
..param.pattern:The @Class.Pattern@ object to search with.
...type:Class.Pattern
..param.needle:The needle object to search for.
...type:Class.String
*/

	template < typename TNeedle, typename TSpec >
	inline void
	setNeedle(Pattern<TNeedle, TSpec> &obj, TNeedle const &ndl) {
		setHost(obj, ndl);
	}


//////////////////////////////////////////////////////////////////////////////

/**
.Function.find:
..summary:Search for a @Class.Pattern@ in a @Class.Finder@ object.
..cat:Pattern Matching
..signature:find(finder, pattern)
..signature:find(finder, pattern, k)
..param.finder:The @Class.Finder@ object to search through.
...remarks:For online-algorithm $patterns$, finder can also be an arbitrary @Concept.Rooted Iterator@.
...type:Class.Finder
...type:Concept.Rooted Iterator
..param.pattern:The @Class.Pattern@ object to search for.
...remarks:For index $finders$, pattern can also be a Sequence.
...type:Class.Pattern
..param.k:Desired minimal score (for approximate matching).
...remarks:$k$ has to be a number <= 0.
...remarks:Differences are deletions, insertions and substitutions.
..returns:$boolean$ that indicates whether an occurence of $pattern$ was found or not.
..remarks:Repeated calls of this function iterate through all occurences of $pattern$.
*/

/**
.Class.Finder:
..summary:Holds the haystack and a current search context.
..cat:Pattern Matching
..signature:Finder<THaystack[, TSpec]>
..param.THaystack:The haystack type.
...type:Class.String
...type:Class.Index
..param.TSpec:The index-algorithm to search with (Optional).
...default:The result of @Metafunction.DefaultFinder@
...remarks:Leave empty for online pattern matching (see @Class.Pattern@).
...remarks:If $THaystack$ is an @Class.Index@, then $TSpec$ specifies the index search algorithm.
..remarks:$position(finder)$ returns the position of the current hit in the haystack.
If $THaystack$ is a set of strings or an index of a set of strings, then $position(finder)$ returns a @Class.Pair@ $(hayNo, pos)$,
in which $hayNo$ is the haystack index and $pos$ the local position of the hit.
..remarks:Use $clear(finder)$ to reset a finder object and search from the beginning.
*/

///.Function.clear.param.object.type:Class.Finder
///.Function.position.param.iterator.type:Class.Finder

	template < typename THaystack, typename TSpec = typename DefaultFinder<THaystack>::Type >
	class Finder
	{
		typedef typename Iterator<THaystack, Rooted>::Type TIterator;

	public:
		TIterator data_iterator;
		bool _empty;

		Finder():
			_empty(true) {}

		Finder(THaystack &haystack):
			data_iterator(begin(haystack, Rooted() )),
			_empty(true) {}

		Finder(TIterator &iter):
			data_iterator(iter),
			_empty(false) {}

		Finder(TIterator const &iter):
		data_iterator(iter),
			_empty(false) {}

		Finder(Finder const &orig):
			data_iterator(orig.data_iterator),
			_empty(orig._empty) {};

//____________________________________________________________________________

		inline typename Reference<TIterator>::Type 
		operator* () 
		{
SEQAN_CHECKPOINT
			return value(hostIterator(*this));
		}

		inline typename Reference<TIterator const>::Type 
		operator* () const
		{
SEQAN_CHECKPOINT
			return value(hostIterator(*this));
		}

//____________________________________________________________________________

		operator TIterator () const
		{
SEQAN_CHECKPOINT
			return data_iterator;
		}

//____________________________________________________________________________

	};



//____________________________________________________________________________

	template <typename THaystack, typename TSpec>
	inline typename _Parameter<THaystack>::Type 
	host(Finder<THaystack, TSpec> & me)
	{
SEQAN_CHECKPOINT
		return container(hostIterator(me));
	}

	template <typename THaystack, typename TSpec>
	inline typename _Parameter<THaystack>::Type 
	host(Finder<THaystack, TSpec> const & me)
	{
SEQAN_CHECKPOINT
		return container(hostIterator(me));
	}

	template <typename THaystack, typename TSpec>
	inline typename _Parameter<THaystack>::Type 
	container(Finder<THaystack, TSpec> & me)
	{
SEQAN_CHECKPOINT
		return container(hostIterator(me));
	}

	template <typename THaystack, typename TSpec>
	inline typename _Parameter<THaystack>::Type 
	container(Finder<THaystack, TSpec> const & me)
	{
SEQAN_CHECKPOINT
		return container(hostIterator(me));
	}

//____________________________________________________________________________

	template <typename THaystack, typename TSpec>
	inline void
	setHost(Finder<THaystack, TSpec> & me, typename _Parameter<THaystack>::Type container_)
	{
SEQAN_CHECKPOINT
		setContainer(hostIterator(me), container_);
	}

	template <typename THaystack, typename TSpec>
	inline void
	setContainer(Finder<THaystack, TSpec> & me, typename _Parameter<THaystack>::Type container_)
	{
SEQAN_CHECKPOINT
		setContainer(hostIterator(me), container_);
	}

//____________________________________________________________________________

	template <typename THaystack, typename TSpec>
	inline typename Iterator<THaystack, Rooted>::Type &
	hostIterator(Finder<THaystack, TSpec> & me)
	{
SEQAN_CHECKPOINT
		return me.data_iterator;
	}

	template <typename THaystack, typename TSpec>
	inline typename Iterator<THaystack, Rooted>::Type const &
	hostIterator(Finder<THaystack, TSpec> const & me)
	{
SEQAN_CHECKPOINT
		return me.data_iterator;
	}

//____________________________________________________________________________

	template <typename THaystack, typename TSpec>
	inline bool
	empty(Finder<THaystack, TSpec> & me)
	{
SEQAN_CHECKPOINT
		return me._empty;
	}

	template <typename THaystack, typename TSpec>
	inline void
	clear(Finder<THaystack, TSpec> & me)
	{
SEQAN_CHECKPOINT
		me._empty = true;
	}
/*
	template <typename THaystack, typename TSpec>
	inline void
	init(Finder<THaystack, TSpec> & me)
	{
SEQAN_CHECKPOINT
		if (me._empty) {
			me._empty = false;
			goBegin(hostIterator(me));
		}
	}
*/

//____________________________________________________________________________

	template <typename THaystack, typename TSpec>
	inline bool
	atBegin(Finder<THaystack, TSpec> & me)
	{
SEQAN_CHECKPOINT
		return (!empty(me) && atBegin(hostIterator(me)));
	}

	template <typename THaystack, typename TSpec>
	inline bool
	atEnd(Finder<THaystack, TSpec> & me)
	{
SEQAN_CHECKPOINT
		return (!empty(me) && atEnd(hostIterator(me)));
	}

//____________________________________________________________________________

	template <typename THaystack, typename TSpec>
	inline void
	goBegin(Finder<THaystack, TSpec> & me)
	{
SEQAN_CHECKPOINT
		me._empty = false;
		goBegin(hostIterator(me));
	}

	template <typename THaystack, typename TSpec>
	inline void
	goEnd(Finder<THaystack, TSpec> & me)
	{
SEQAN_CHECKPOINT
		me._empty = false;
		goEnd(hostIterator(me));
	}

//____________________________________________________________________________

	template <typename THaystack, typename TSpec>
	inline typename Position<Finder<THaystack, TSpec> >::Type
	position(Finder<THaystack, TSpec> & me)
	{
SEQAN_CHECKPOINT
		if (empty(me)) return 0;
		return hostIterator(me) - begin(container(me), Rooted());
	}

	template <typename THaystack, typename TSpec>
	inline typename Position<Finder<THaystack, TSpec> >::Type
	position(Finder<THaystack, TSpec> const & me)
	{
SEQAN_CHECKPOINT
		if (empty(me)) return 0;
		return hostIterator(me) - begin(container(me), Rooted());
	}

	template <typename THaystack, typename TSpec, typename TPosition>
	inline void 
	setPosition(Finder<THaystack, TSpec> & me, TPosition pos_)
	{
SEQAN_CHECKPOINT
		hostIterator(me) = begin(container(me), Rooted()) + pos_;
	}

//____________________________________________________________________________

	template <typename THaystack, typename TSpec>
	inline Finder<THaystack, TSpec> &
	operator--(Finder<THaystack, TSpec> & me)
	{
SEQAN_CHECKPOINT
		--hostIterator(me);
		return me;
	}

	template <typename THaystack, typename TSpec>
	inline Finder<THaystack, TSpec> &
	operator++(Finder<THaystack, TSpec> & me)
	{
SEQAN_CHECKPOINT
/*			if (beforeBegin()) {
			goBegin(hostIterator(me));
		} else*/
			++hostIterator(me);
		return me;
	}

//////////////////////////////////////////////////////////////////////////////
// operator +
//////////////////////////////////////////////////////////////////////////////

	template <typename THaystack, typename TSpec, typename TIntegral>
	inline Finder<THaystack, TSpec> const
	operator + (Finder<THaystack, TSpec> const & left, TIntegral right)
	{
SEQAN_CHECKPOINT
		return Finder<THaystack, TSpec>(hostIterator(left) + right);
	}

//////////////////////////////////////////////////////////////////////////////
// operator +=
//////////////////////////////////////////////////////////////////////////////

	template <typename THaystack, typename TSpec, typename TIntegral>
	inline Finder<THaystack, TSpec> &
	operator += (Finder<THaystack, TSpec> & left,
					TIntegral right)
	{
SEQAN_CHECKPOINT
		hostIterator(left) += right;
		return left;
	}

//////////////////////////////////////////////////////////////////////////////
// operator -
//////////////////////////////////////////////////////////////////////////////

	template <typename THaystack, typename TSpec, typename TIntegral>
	inline Finder<THaystack, TSpec> const
	operator - (Finder<THaystack, TSpec> const & left, TIntegral right)
	{
SEQAN_CHECKPOINT
		return Finder<THaystack, TSpec>(hostIterator(left) - right);
	}

	template <typename THaystack, typename TSpec, typename TIntegral>
	inline typename Difference<Finder<THaystack, TSpec> const>::Type
	operator - (Finder<THaystack, TSpec> const & left, Finder<THaystack, TSpec> const & right)
	{
SEQAN_CHECKPOINT
		return hostIterator(left) - hostIterator(right);
	}

//////////////////////////////////////////////////////////////////////////////
// operator -=
//////////////////////////////////////////////////////////////////////////////

	template <typename THaystack, typename TSpec, typename TIntegral>
	inline Finder<THaystack, TSpec> &
	operator -= (Finder<THaystack, TSpec> & left,
					TIntegral right)
	{
SEQAN_CHECKPOINT
		hostIterator(left) -= right;
		return left;
	}

//____________________________________________________________________________


/**
.Function.setHaystack:
..summary:Sets the haystack of a @Class.Finder@ object.
..cat:Pattern Matching
..signature:setHaystack(finder, haystack)
..param.finder:The @Class.Finder@ object to search with.
...type:Class.Finder
..param.haystack:The haystack object the finder searches through.
...type:Class.String
*/

	template < typename THaystack, typename TSpec >
	inline void
	setHaystack(Finder<THaystack, TSpec> &obj, THaystack const &hstk) {
		setHost(obj, hstk);
	}

/**
.Function.haystack:
..summary:Returns the haystack of a @Class.Finder@ object.
..cat:Pattern Matching
..signature:haystack(finder)
..param.finder:The @Class.Finder@ object to search through.
...type:Class.Finder
..returns:The haystack object.
..remarks:The result type is @Metafunction.Haystack@$<TFinder>::Type$ for finder of type $TFinder$.
*/

	template < typename TObject >
	inline typename Haystack<TObject>::Type &
	haystack(TObject &obj) {
		return container(obj);
	}

	template < typename TObject >
	inline typename Haystack<TObject const>::Type &
	haystack(TObject const &obj) {
		return container(obj);
	}


//////////////////////////////////////////////////////////////////////////////


	template <typename THaystack, typename TSpec>
	struct Container< Finder<THaystack, TSpec> > {
		typedef THaystack Type;
	};

	template <typename THaystack, typename TSpec>
	struct Container< Finder<THaystack, TSpec> const> {
		typedef THaystack const Type;
	};

	template <typename THaystack, typename TSpec>
	struct Host< Finder<THaystack, TSpec> > {
		typedef THaystack Type;
	};

	template <typename THaystack, typename TSpec>
	struct Host< Finder<THaystack, TSpec> const> {
		typedef THaystack const Type;
	};


	template <typename THaystack, typename TSpec>
	struct Value< Finder<THaystack, TSpec> > {
		typedef typename Value<THaystack>::Type Type;
	};

	template <typename THaystack, typename TSpec>
	struct Position< Finder<THaystack, TSpec> > {
		typedef typename Position<THaystack>::Type Type;
	};

	template <typename THaystack, typename TSpec>
	struct Difference< Finder<THaystack, TSpec> > {
		typedef typename Difference<THaystack>::Type Type;
	};

	template <typename THaystack, typename TSpec>
	struct Size< Finder<THaystack, TSpec> > {
		typedef typename Size<THaystack>::Type Type;
	};

//____________________________________________________________________________


	template <typename TNeedle, typename TSpec>
	struct Container< Pattern<TNeedle, TSpec> > {
		typedef TNeedle Type;
	};

	template <typename TNeedle, typename TSpec>
	struct Container< Pattern<TNeedle, TSpec> const > {
		typedef TNeedle const Type;
	};

	template <typename TNeedle, typename TSpec>
	struct Host< Pattern<TNeedle, TSpec> > {
		typedef TNeedle Type;
	};

	template <typename TNeedle, typename TSpec>
	struct Host< Pattern<TNeedle, TSpec> const > {
		typedef TNeedle const Type;
	};


	template <typename TPattern, typename TSpec>
	struct Value< Pattern<TPattern, TSpec> > {
		typedef typename Value<TPattern>::Type Type;
	};

	template <typename TPattern, typename TSpec>
	struct Position< Pattern<TPattern, TSpec> > {
		typedef typename Position<TPattern>::Type Type;
	};

	template <typename TPattern, typename TSpec>
	struct Difference< Pattern<TPattern, TSpec> > {
		typedef typename Difference<TPattern>::Type Type;
	};

	template <typename TPattern, typename TSpec>
	struct Size< Pattern<TPattern, TSpec> > {
		typedef typename Size<TPattern>::Type Type;
	};


//////////////////////////////////////////////////////////////////////////////

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
