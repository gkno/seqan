/*
 *  Created by Anne-Katrin Emde on 3.1.2007.
 */

#ifndef SEQAN_HEADER_INDEX_QGRAM_H
#define SEQAN_HEADER_INDEX_QGRAM_H

namespace SEQAN_NAMESPACE_MAIN
{


//////////////////////////////////////////////////////////////////////////////
// q-gram index fibres

	struct _Fibre_Dir;			// directory/hash table, contains start indices of buckets

	typedef Tag<_Fibre_Text>		QGram_Text;
	typedef Tag<_Fibre_RawText>		QGram_RawText;
	typedef Tag<_Fibre_SA>			QGram_SA;
	typedef Tag<_Fibre_RawSA>		QGram_RawSA;
	typedef Tag<_Fibre_Dir>			QGram_Dir;

	struct QGram_Alg {};

//////////////////////////////////////////////////////////////////////////////
// q-gram index

/**
.Spec.Index_QGram:
..summary:An index based on an array of sorted q-grams.
..cat:Index
..general:Class.Index
..signature:Index<TText, Index_QGram<> >
..param.TText:The text type.
...type:Class.String
..remarks:The fibres (see @Class.Index@ and @Metafunction.Fibre@) of this index are a suffix array sorted by the first q characters (see @Tag.QGram_SA@) and a q-gram directory (see @Tag.QGram_Dir@).
*/

	template < typename TShapeSpec >
	struct Index_QGram;


	template < typename TText, typename TShapeSpec >
	class Index<TText, Index_QGram<TShapeSpec> > {
	public:
		Holder<typename Fibre<Index, QGram_Text>::Type>	text;
		typename Fibre<Index, QGram_SA>::Type			sa;		// suffix array sorted by the first q chars
		typename Fibre<Index, QGram_Dir>::Type			dir;	// bucket directory

		typedef Shape<typename Value<Index>::Type, TShapeSpec>	TShape;

		TShape shape;

		Index() {}

		template <typename _TText>
		Index(_TText &_text):
			text(_text) {}

		template <typename _TText>
		Index(_TText const &_text):
			text(_text) {}

		template <typename _TText, typename _TShape>
		Index(_TText &_text, _TShape const &_shape):
			text(_text),
			shape(_shape) {}

		template <typename _TText, typename _TShape>
		Index(_TText const &_text, _TShape const &_shape):
			text(_text),
			shape(_shape) {}
	};

    template < typename TText, typename TSpec >
    struct Value< Index<TText, Index_QGram<TSpec> > > {
		typedef typename Value< typename Fibre< Index<TText, Index_QGram<TSpec> >, QGram_RawText >::Type >::Type Type;
    };

	template < typename TText, typename TSpec >
    struct Size< Index<TText, Index_QGram<TSpec> > > {
		typedef typename Size< typename Fibre< Index<TText, Index_QGram<TSpec> >, QGram_RawText >::Type >::Type Type;
    };


//////////////////////////////////////////////////////////////////////////////
// default fibre creators

	template < typename TText, typename TShapeSpec >
	struct DefaultIndexCreator<Index<TText, Index_QGram<TShapeSpec> >, Tag<_Fibre_SA> > {
        typedef QGram_Alg Type;
    };

	template < typename TText, typename TShapeSpec >
	struct DefaultIndexCreator<Index<TText, Index_QGram<TShapeSpec> >, QGram_Dir > {
        typedef QGram_Alg Type;
    };

//////////////////////////////////////////////////////////////////////////////

	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec>, QGram_Dir>::Type & 
	getFibre(Index<TText, TSpec> &index, QGram_Dir const) {
		return index.dir;
	}
	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec> const, QGram_Dir>::Type & 
	getFibre(Index<TText, TSpec> const &index, QGram_Dir const) {
		return index.dir;
	}

//////////////////////////////////////////////////////////////////////////////
/**
.Function.indexDir:
..summary:Shortcut for $getFibre(.., QGram_Dir)$.
..cat:Index
..signature:indexDir(index)
..param.index:The @Class.Index@ object holding the fibre.
...type:Spec.Index_QGram
..returns:A reference to the @Tag.QGram_Dir@ fibre (q-gram directory).
*/

	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec>, QGram_Dir>::Type & 
	indexDir(Index<TText, TSpec> &index) { 
		return getFibre(index, QGram_Dir()); 
	}
	template <typename TText, typename TSpec>
	inline typename Fibre<Index<TText, TSpec> const, QGram_Dir>::Type & 
	indexDir(Index<TText, TSpec> const &index) { 
		return getFibre(index, QGram_Dir()); 
	}

	template <typename TText, typename TShapeSpec>
	inline Shape<typename Value<Index<TText, Index_QGram<TShapeSpec> > >::Type, TShapeSpec> &
	indexShape(Index<TText, Index_QGram<TShapeSpec> > &index) {
		return index.shape; 
	}
	template <typename TText, typename TShapeSpec>
	inline Shape<typename Value<Index<TText, Index_QGram<TShapeSpec> > >::Type, TShapeSpec> const &
	indexShape(Index<TText, Index_QGram<TShapeSpec> > const &index) {
		return index.shape; 
	}

//////////////////////////////////////////////////////////////////////////////
/**
.Function.createQGramIndex:
..summary:Builds a q-gram index on a sequence. 
..cat:Index
..signature:createQGramIndex(sa, dir, text, shape)
..param.text:The sequence.
..param.shape:The shape to be used.
...type:Class.Shape
..param.sa:The resulting list in which all q-grams are sorted alphabetically.
..param.dir:The resulting array that indicates at which position in index the corresponding q-grams can be found.
..returns:Index contains the sorted list of qgrams. For each possible q-gram pos contains the first position in index that corresponds to this q-gram. 
*/

	template < typename TSA, 
			typename TDir,
			typename TText,
			typename TShape >
	void createQGramIndex(
		TSA &sa,
		TDir &dir,
		TText const &text,
		TShape &shape)
	{
	SEQAN_CHECKPOINT
		typedef typename Iterator<TText const, Standard>::Type	TIterator;
		typedef typename Iterator<TDir, Standard>::Type			TDirIterator;
		typedef typename Value<TShape>::Type					TValue;
		typedef typename Position<TText>::Type					TPosition;
		typedef typename Size<TText>::Type						TSize;
		
		// clear counters
		arrayFill(begin(dir, Standard()), end(dir, Standard()), 0);

		//count q-grams
		TIterator qgram_1_it = begin(text,Standard());
		TIterator qgram_2_it = qgram_1_it + 1;
		TValue x, y = 0;			// hashwert des q-Grams
		x = hash(shape, qgram_1_it);
		++dir[x];
		TPosition i = 1;
		TSize num_qgrams = length(text) - shapeSpan(shape) + 1;
		while( i < num_qgrams)
		{
			y = hashNext(shape, qgram_1_it, qgram_2_it, x);
			++qgram_1_it;
			++qgram_2_it;
			++dir[y];
			x = y;
			++i;
		}

		//cumulative sum
		{
			TDirIterator it = begin(dir, Standard()), itEnd = end(dir, Standard());
			TSize diff = 0, diff_prev = 0, sum = 0;
			while (it != itEnd) {
				sum += diff_prev;
				diff_prev = diff;
				diff = *it;
				*it = sum;
				++it;
			}
			sum += diff_prev;
			dir[0] = sum;

			SEQAN_ASSERT(sum + diff == num_qgrams);
		}
		
		//build sa
		qgram_1_it = begin(text,Standard());
		qgram_2_it = qgram_1_it + 1;
		TValue lastBucketIndex = length(dir) - 1;

		// first hash
		x = hash(shape, qgram_1_it);
		if (x != lastBucketIndex)
			sa[dir[x + 1]++] = 0;
		else
			sa[dir[0]++] = 0;
		i = 1;
		while(i < num_qgrams)
		{
			// next hash
			y = hashNext(shape, qgram_1_it, qgram_2_it, x);
			++qgram_1_it;
			++qgram_2_it;
			if (y != lastBucketIndex)
				sa[dir[y + 1]++] = i++;
			else
				sa[dir[0]++] = i++;

			// für die nächste Runde
			x = y;
		}
		dir[0] = 0;
	}

	template < 
		typename TSA, 
		typename TDir,
		typename TString,
		typename TSpec,
		typename TShape,
		typename TLimitsString >
	void createQGramIndex(
		TSA &sa,
		TDir &dir,
		StringSet<TString, TSpec> const &stringSet,
		TShape &shape,
		TLimitsString &limits)
	{
	SEQAN_CHECKPOINT
		typedef typename Iterator<TString const, Standard>::Type	TIterator;
		typedef typename Iterator<TDir, Standard>::Type				TDirIterator;
		typedef typename Value<TShape>::Type						TValue;
		typedef typename Size<TString>::Type						TSize;
		
		// 1. clear counters
		arrayFill(begin(dir, Standard()), end(dir, Standard()), 0);

		// 2. count q-grams
		for(unsigned seqNo = 0; seqNo < length(stringSet); ++seqNo) 
		{
			TString const &sequence = value(stringSet, seqNo);

			//count q-grams
			TIterator it_l = begin(sequence, Standard());
			TIterator it_r = it_l + 1;
			TValue x, y = 0;			// hashwert des q-Grams
			x = hash(shape, it_l);
			++dir[x];
			TSize	i = 1;
			TSize	num_qgrams = length(sequence) - shapeSpan(shape) + 1;
			while (i < num_qgrams)
			{
				y = hashNext(shape, it_l, it_r, x);
				++it_l;
				++it_r;
				++dir[y];
				x = y;
				++i;
			}
		}

		// 3. cumulative sum
		{
			TDirIterator it = begin(dir, Standard());
			TDirIterator itEnd = end(dir, Standard());
			TSize diff = 0, diff_prev = 0, sum = 0;
			while (it != itEnd) {
				sum += diff_prev;
				diff_prev = diff;
				diff = *it;
				*it = sum;
				++it;
			}
			sum += diff_prev;
			dir[0] = sum;

			SEQAN_ASSERT(sum + diff == length(sa));
		}
		
		// 4. fill suffix array
		for(unsigned seqNo = 0; seqNo < length(stringSet); ++seqNo) 
		{
			TString const &sequence = value(stringSet, seqNo);

			TIterator it_l = begin(sequence, Standard());
			TIterator it_r = it_l + 1;
			TValue lastBucketIndex = length(dir) - 1;

			typename Value<TSA>::Type localPos;
			setValueI1(localPos, seqNo);
			setValueI2(localPos, 0);

			// first hash
			TValue x = hash(shape, it_l), y = 0;
			if (x != lastBucketIndex)
				sa[dir[x + 1]++] = localPos;
			else
				sa[dir[0]++] = localPos;

			TSize	i = 1;
			TSize	num_qgrams = length(sequence) - shapeSpan(shape) + 1;

			while (i < num_qgrams)
			{
				// next hash
				y = hashNext(shape, it_l, it_r, x);
				++it_l;
				++it_r;
				setValueI2(getValueI2(localPos) + 1);
				if (y != lastBucketIndex)
					sa[dir[y + 1]++] = localPos;
				else
					sa[dir[0]++] = localPos;

				// für die nächste Runde
				x = y;
			}
		}
		dir[0] = 0;
	}


//////////////////////////////////////////////////////////////////////////////
// create q-gram index of *one* sequence in external memory 

	// *** COMPARATORS & MAPS ***
        
    template <typename InType, typename Result = int>
    struct _qgram_comp : public ::std::binary_function<InType,InType,Result> {
        inline Result operator()(const InType &a, const InType &b) const
        {
			typedef typename Value<InType, 2>::Type TQGram;
			typename Value<TQGram>::Type const *sa = a.i2.i;
            typename Value<TQGram>::Type const *sb = b.i2.i;
			typename Value<TQGram>::Type const *saEnd = sa + length(a.i2);

            for(; sa != saEnd; ++sa, ++sb) {
                if (*sa == *sb) continue;
                return (*sa < *sb)? -1 : 1;
            }
			return posCompare(a.i1, b.i1);
        }
    };

    // optimized for bitvectors
    template <typename T1, typename TValue, unsigned _size, typename Result>
    struct _qgram_comp< Pair<T1, Tuple<TValue, _size, Compressed>, Compressed >, Result > :
        public ::std::binary_function<
            Pair<T1, Tuple<TValue, _size, Compressed>, Compressed >,
            Pair<T1, Tuple<TValue, _size, Compressed>, Compressed >,
            Result> {       
        inline Result operator()(
            const Pair<T1, Tuple<TValue, _size, Compressed>, Compressed > &a,
            const Pair<T1, Tuple<TValue, _size, Compressed>, Compressed > &b) const
        {
            if (a.i2 < b.i2) return -1;
            if (a.i2 > b.i2) return 1;
			return posCompare(a.i1, b.i1);
        }
    };


    template <typename TValue, typename TResult = unsigned>
    struct _qgram_hash : public ::std::unary_function<TValue, TResult> {
        inline TResult operator()(TValue const &a) const
        {
			typedef typename Value<TValue, 2>::Type	TQGram;
			TResult hash = 0;
			unsigned len = length(a.i2);
            for(unsigned i = 0; i < len; ++i) {
				hash *= ValueSize< typename Value<TQGram>::Type >::VALUE;
				hash += (TResult)a.i2[i];
            }
            return hash;
        }
    };

	// TODO: replace fixed tuple size of 6 with q and add q to Shape template arguments
	template < 
		typename TSA, 
		typename TDir,
		typename TText,
		typename TShape >
	void createQGramIndexExt(
		TSA &suffixArray,
		TDir &dir,
		TText &text,
		TShape &shape)
	{
        // signed characters behave different than unsigned when compared
        // to get the same index with signed or unsigned chars we simply cast them to unsigned
        // before feeding them into the pipeline
        typedef typename _MakeUnsigned< typename Value<TText>::Type >::Type TUValue;

        // *** SPECIALIZATION ***

		typedef Pipe< TText, Source<> >				TSource;
        typedef Pipe< TSource, Caster<TUValue> >    TUnsigner;
		typedef Pipe< TUnsigner, Tupler<7> >	    TTupler;
						                typedef _qgram_comp<_TypeOf(TTupler)> qcomp_t;
        typedef Pool< 
					_TypeOf(TTupler), 
					SorterSpec< SorterConfigSize<qcomp_t, _TSizeOf(TTupler) > > 
				> TSortTuples;
										typedef _qgram_hash<_TypeOf(TTupler), typename Size<TDir>::Type> qhash_t;

        // *** INSTANTIATION ***

		TSource			src(text);
        TUnsigner		unsigner(src);
		TTupler			tupler(unsigner);
		TSortTuples		sorter;

		// sort q-grams
		sorter << tupler;

		// fill sa and dir
		if (!beginRead(sorter)) return;

		typename Iterator<TSA>::Type itSA = begin(suffixArray);
		typename Iterator<TDir>::Type itDir = begin(dir);

		qcomp_t	qcomp;
		qhash_t qhash;

		typename Value<TSortTuples>::Type	old_qgram;
		typename Size<TDir>::Type			hash, old_hash = 0;
        typename Size<TSortTuples>::Type	leftToRead = length(sorter);
		bool first = true;

		while (leftToRead) {
			// copy occurence position
			*itSA = (*sorter).i1;
			if (first || qcomp(old_qgram, *sorter) != 0) {
				old_qgram = *sorter;
				hash = qhash(old_qgram);

				SEQAN_ASSERT(old_hash < hash);

				// copy bucket begin
				typename Size<TSortTuples>::Type i = length(sorter) - leftToRead;
				for(; old_hash < hash; ++old_hash, ++itDir)
					*itDir = i;
				first = false;
			}
			++src;
			--leftToRead;
		}

		// fill bucket table
		typename Size<TSortTuples>::Type i = length(sorter);
		hash = length(dir);
		for(; old_hash < hash; ++old_hash, ++itDir)
			*itDir = i;

		endRead(sorter);
	}


//////////////////////////////////////////////////////////////////////////////
// create q-gram index of *multiple* sequences in external memory 

	template < 
		typename TSA, 
		typename TDir,
		typename TString,
		typename TSpec,
		typename TShape,
		typename TLimitsString >
	void createQGramIndexExt(
		TSA &suffixArray,
		TDir &dir,
		StringSet<TString, TSpec> const &stringSet,
		TShape &shape,
		TLimitsString &limits)
	{
        // signed characters behave different than unsigned when compared
        // to get the same index with signed or unsigned chars we simply cast them to unsigned
        // before feeding them into the pipeline
		typedef typename Concatenator<StringSet<TString, TSpec> >::Type			TConcat;
        typedef typename _MakeUnsigned< typename Value<TConcat>::Type >::Type	TUValue;
		typedef Multi<
			Tupler<7, true, Compressed>, 
			typename Value<TSA>::Type,
			typename StringSetLimits< StringSet<TString, TSpec> >::Type >		TTuplerSpec;

        // *** SPECIALIZATION ***

		typedef Pipe< TConcat, Source<> >			TSource;
        typedef Pipe< TSource, Caster<TUValue> >    TUnsigner;
		typedef Pipe< TUnsigner, TTuplerSpec >	    TTupler;
						                typedef _qgram_comp<_TypeOf(TTupler)> qcomp_t;
        typedef Pool< 
					_TypeOf(TTupler), 
					SorterSpec< SorterConfigSize<qcomp_t, _TSizeOf(TTupler) > > 
				> TSortTuples;
										typedef _qgram_hash<_TypeOf(TTupler), typename Size<TDir>::Type> qhash_t;

        // *** INSTANTIATION ***

		TSource			src(concat(stringSet));
        TUnsigner		unsigner(src);
		TTupler			tupler(unsigner, limits);
		TSortTuples		sorter;

		// sort q-grams
		sorter << tupler;

		// fill sa and dir
		if (!beginRead(sorter)) return;

		typename Iterator<TSA>::Type itSA = begin(suffixArray);
		typename Iterator<TDir>::Type itDir = begin(dir);

		qcomp_t	qcomp;
		qhash_t qhash;

		typename Value<TSortTuples>::Type	old_qgram;
		typename Size<TDir>::Type			hash, old_hash = 0;
        typename Size<TSortTuples>::Type	leftToRead = length(sorter);
		bool first = true;

		while (leftToRead) {
			// copy occurence position
			*itSA = (*sorter).i1;
			if (first || qcomp(old_qgram, *sorter) != 0) {
				old_qgram = *sorter;
				hash = qhash(old_qgram);

				SEQAN_ASSERT(old_hash < hash);

				// copy bucket begin
				typename Size<TSortTuples>::Type i = length(sorter) - leftToRead;
				for(; old_hash < hash; ++old_hash, ++itDir)
					*itDir = i;
				first = false;
			}
			++src;
			--leftToRead;
		}

		// fill bucket table
		typename Size<TSortTuples>::Type i = length(sorter);
		hash = length(dir);
		for(; old_hash < hash; ++old_hash, ++itDir)
			*itDir = i;

		endRead(sorter);
	}



//////////////////////////////////////////////////////////////////////////////
// interface for automatic index creation 

	template <typename TText, typename TShapeSpec>
	inline bool indexCreate(Index<TText, Index_QGram<TShapeSpec> > &index, Tag<_Fibre_SA> const, QGram_Alg const) 
	{
		typedef Index<TText, Index_QGram<TShapeSpec> >				TIndex;
		typedef Shape<typename Value<TIndex>::Type, SimpleShape>	TShape;

		TShape &shape = indexShape(index);
		stringToShape(shape, "xxxxxxx");

		// count all overlapping q-grams
		typename Size<TIndex>::Type qgram_count = 0;
		for(unsigned i = 0; i < countSequences(index); ++i)
			if (sequenceLength(i, index) >= shapeSpan(shape))
				qgram_count += sequenceLength(i, index) - (shapeSpan(shape) - 1);

		resize(indexSA(index), qgram_count, Exact());
		resize(indexDir(index), intPow(ValueSize<TShape>::VALUE, shapeSpan(shape) - shapeCountBlanks(shape)), Exact());
		createQGramIndex(indexSA(index), indexDir(index), indexText(index), indexShape(index));
		return true;
	}

	template <typename TText, typename TShapeSpec>
	inline bool indexCreate(Index<TText, Index_QGram<TShapeSpec> > &index, QGram_Dir const, QGram_Alg const alg)
	{
		return indexCreate(index, QGram_SA(), alg);
	}

//////////////////////////////////////////////////////////////////////////////
// open

	template < typename TObject, typename TSpec >
	inline bool open(
		Index< TObject, Index_QGram<TSpec> > &index, 
		const char *fileName,
		int openMode)
	{
		String<char> name;
		name = fileName;	append(name, ".txt");	
		if ((!open(getFibre(index, QGram_Text()), toCString(name), openMode)) && 
			(!open(getFibre(index, QGram_Text()), fileName), openMode)) return false;

		name = fileName;	append(name, ".sa");	open(getFibre(index, QGram_SA()), toCString(name), openMode);
		name = fileName;	append(name, ".dir");	open(getFibre(index, QGram_Dir()), toCString(name), openMode);
		return true;
	}
	template < typename TObject, typename TSpec >
	inline bool open(
		Index< TObject, Index_QGram<TSpec> > &index, 
		const char *fileName) 
	{
		return open(index, fileName, OPEN_RDONLY);
	}


//////////////////////////////////////////////////////////////////////////////
// save

	template < typename TObject, typename TSpec >
	inline bool save(
		Index< TObject, Index_QGram<TSpec> > &index, 
		const char *fileName,
		int openMode)
	{
		String<char> name;
		name = fileName;	append(name, ".txt");	
		if ((!save(getFibre(index, QGram_Text()), toCString(name), openMode)) && 
			(!save(getFibre(index, QGram_Text()), fileName), openMode)) return false;

		name = fileName;	append(name, ".sa");	save(getFibre(index, QGram_SA()), toCString(name), openMode);
		name = fileName;	append(name, ".dir");	save(getFibre(index, QGram_Dir()), toCString(name), openMode);
		return true;
	}
	template < typename TObject, typename TSpec >
	inline bool save(
		Index< TObject, Index_QGram<TSpec> > &index, 
		const char *fileName)
	{
		return save(index, fileName, OPEN_WRONLY | OPEN_CREATE);
	}

}

#endif //#ifndef SEQAN_HEADER_...
