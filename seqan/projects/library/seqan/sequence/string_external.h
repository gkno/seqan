/*
 *  string_external.h
 *  genindex
 *
 *  Created by David Weese on 17.07.05.
 *
 */

#ifndef SEQAN_HEADER_STRING_EXTERNAL_H
#define SEQAN_HEADER_STRING_EXTERNAL_H

namespace SEQAN_NAMESPACE_MAIN
{

    // standard external string
    // size is uint32
    template < typename _TFile = File<>,            // default file type
               unsigned _PageSize = 4 * 1024 * 1024,// 1MTypes per default
			   unsigned _Frames = 2,			    // simultanous frames
			   bool		_TempByDefault = true >		// open temp. file in default c'tor
    struct ExternalConfig {
        typedef _TFile TFile;
        typedef unsigned SizeType;
        enum { PageSize = _PageSize };
        enum { Frames = _Frames };
		enum { TempByDefault = _TempByDefault };
    };

    // the same as ExternalConfig
    // but size type is size type of _TFile (i.e. uint64)
    //
    // ATTENTION:
    // pipes use the size type 
    // uint64 blows up your suffix arrays, lcp-tables, ...
    template < typename _TFile = File<>,            // default file type
               unsigned _PageSize = 1 * 1024 * 1024,// 1MTypes per default
			   unsigned _Frames = 2,			    // simultanous frames
			   bool		_TempByDefault = true >		// open temp. file in default c'tor
    struct ExternalConfigLarge {
        typedef _TFile TFile;
        typedef typename Size<_TFile>::Type SizeType;
        enum { PageSize = _PageSize };
        enum { Frames = _Frames };
		enum { TempByDefault = _TempByDefault };
    };

    template < typename _TFile = File<>,            // default file type
               unsigned _PageSize = 1 * 1024 * 1024,// 1MTypes per default
			   unsigned _Frames = 2 >			    // simultanous frames
    struct ExternalConfigManualOpen {
        typedef _TFile TFile;
        typedef typename Size<_TFile>::Type SizeType;
        enum { PageSize = _PageSize };
        enum { Frames = _Frames };
		enum { TempByDefault = false };				// don't open temp. file in default c'tor
    };

    // custom size type
    template < typename TSize,
		       typename _TFile = File<>,            // default file type
               unsigned _PageSize = 1 * 1024 * 1024,// 1MTypes per default
			   unsigned _Frames = 2,			    // simultanous frames {
			   bool		_TempByDefault = true >		// open temp. file in default c'tor
    struct ExternalConfigSize {
		typedef TSize SizeType;
        typedef _TFile TFile;
        enum { PageSize = _PageSize };
        enum { Frames = _Frames };
		enum { TempByDefault = _TempByDefault };
    };

    template < typename TConfig = ExternalConfig<> >
    struct External;


	// calculates the smallest power of 2 bigger than i
	// and returns the exponent
	template < typename ST >
	ST ceilPower2(ST i) {
		ST e = 0;
		for(ST j = 1; j < i; j = j << 1)
			++e;
		return e;
	}

/*
	//////////////////////////////////////////////////////////////////////////////
	// simple vector based page table
	template < typename T >
	struct PageTable : public std::vector<T>
	{
		typedef std::vector<T> Base;

		PageMapper(SizeType _maxPages, T _default = 0):
			Base(_maxPages, _default) {}
	};
*/

/*	//////////////////////////////////////////////////////////////////////////////
	// dynamically growing page table
	template < typename T,
			   typename ST,
			   unsigned _maxPagesBits = 32>
	struct PageTableDynamic
	{
		typedef T	Type;
		typedef T*	TypePtr;
		typedef ST	SizeType;

		TypePtr		pageTbl[_maxPagesBits];	// page tables indexed by the exponent (2^32 pages max)

		PageMapper() {
			memset(pageTbl, 0, sizeof(pageTbl));
		}

		~PageMapper() {
			for(unsigned i = 0; i < _maxPageBits; ++i)
				delete[] pageTbl[i];
		}

		Type operator[] (SizeType _pageNo) {
			unsigned e = ceilPower2(_pageNo);
			if (e < 10) e = 10;					// begin partitioning with entry 4096
			TypePtr tbl = pageTbl[e];
			if (!tbl) {
				tbl = pageTbl[e] = new Type[1 << e];
				memset(tbl, 0, sizeof(Type) * (1 << e));
			}
			if (e) _pageNo -= 1 << (e - 1);
			return tbl[_pageNo];
		}
	};
*/

    //////////////////////////////////////////////////////////////////////////////
	// random vector iterator
	template < typename _Vector >
	struct VectorIterator
	{
		typedef _Vector		                        Vector;
		typedef VectorIterator				        iterator;
        typedef typename Vector::Type			    Type;
		typedef typename Vector::SizeType		    SizeType;
		typedef typename Vector::const_reference	const_reference;
		typedef typename Vector::volatile_ptr		volatile_ptr;
        enum { _PageSize = Vector::PageSize };

		SizeType	offset;
		Vector		*vector;

        typedef VectorIterator                      std_iterator;

		//////////////////////////////////////////////////////////////////////////////
		// public iterator interface

		typedef std::random_access_iterator_tag		iterator_category;
		typedef typename Vector::value_type			value_type;
		typedef typename Vector::difference_type	difference_type;
		typedef typename Vector::pointer			pointer;
		typedef typename Vector::reference			reference;
		
	    VectorIterator():
			offset(0),
			vector(NULL) {}

	    explicit VectorIterator(Vector *_vector, SizeType _offset):
			vector(_vector),
			offset(_offset) {}

		//////////////////////////////////////////////////////////////////////////////
		// iterator conversion interface

        VectorIterator(const std_iterator &I):
			offset(I.offset),
			vector(I.vector) {}

		//////////////////////////////////////////////////////////////////////////////
		// iterator arithmetic

        inline difference_type operator- (const iterator &I) const {
			return offset - I.offset;
		};
		
		inline iterator operator- (difference_type delta) const {
			return iterator(vector, offset - delta);
		};
		
		inline iterator& operator-= (difference_type delta) const {
			offset -= delta;
			return *this;
		};
		
		inline iterator operator+ (difference_type delta) const {
			return iterator(vector, offset + delta);
		};
		
		inline iterator& operator+= (difference_type delta) const {
			offset += delta;
			return *this;
		};
		
		inline reference operator* () {
			return (*vector)[offset];
		}
    
		inline const_reference operator* () const {
			return (*vector)[offset];
		}
    
		inline iterator& operator++ () {
			++offset; return *this;
		}

		inline iterator operator++ (int) {
			iterator before = *this;
			++offset; return before;
		}

		inline iterator& operator-- () {
			--offset; return *this;
		}

		inline iterator operator-- (int) {
			iterator before = *this;
			--offset; return before;
		}

		inline bool operator== (const iterator &I) const {
			SEQAN_ASSERT(vector == I.vector);
			return offset == I.offset;
		}

		inline bool operator!= (const iterator &I) const {
			SEQAN_ASSERT(vector == I.vector);
			return offset != I.offset;
		}

		inline bool operator< (const iterator &I) const {
			SEQAN_ASSERT(vector == I.vector);
			return offset < I.offset;
		}
	};


	//////////////////////////////////////////////////////////////////////////////
	// const random vector iterator
	template < typename _Vector >
	struct VectorConstIterator
	{
		typedef _Vector		                        Vector;
		typedef VectorConstIterator		            iterator;
        typedef typename Vector::Type			    Type;
		typedef typename Vector::SizeType		    SizeType;
		typedef typename Vector::const_reference	const_reference;
		typedef typename Vector::volatile_ptr	    volatile_ptr;
        enum { _PageSize = Vector::PageSize };

		SizeType	offset;
		Vector		*vector;
		
        friend struct VectorIterator<Vector>;
        typedef VectorIterator<Vector>              std_iterator;
        typedef VectorConstIterator                 std_const_iterator;

        //////////////////////////////////////////////////////////////////////////////
		// public iterator interface

		typedef std::random_access_iterator_tag		iterator_category;
		typedef typename Vector::value_type			value_type;
		typedef typename Vector::difference_type	difference_type;
		typedef typename Vector::const_pointer	    pointer;
		typedef typename Vector::const_reference	reference;
		
	    VectorConstIterator():
			offset(0),
			vector(NULL) {}

	    explicit VectorConstIterator(Vector *_vector, SizeType _offset):
			vector(_vector),
			offset(_offset) {}

		//////////////////////////////////////////////////////////////////////////////
		// iterator conversion interface

		VectorConstIterator(const std_iterator &I):
			offset(I.offset),
			vector(I.vector) {}

		VectorConstIterator(const std_const_iterator &I):
			offset(I.offset),
			vector(I.vector) {}

		//////////////////////////////////////////////////////////////////////////////
		// iterator arithmetic

		inline difference_type operator- (const iterator &I) const {
			return offset - I.offset;
		};
		
		inline iterator operator- (difference_type delta) const {
			return iterator(vector, offset - delta);
		};
		
		inline iterator& operator-= (difference_type delta) const {
			offset -= delta;
			return *this;
		};
		
		inline iterator operator+ (difference_type delta) const {
			return iterator(vector, offset + delta);
		};
		
		inline iterator& operator+= (difference_type delta) const {
			offset += delta;
			return *this;
		};
		
		inline const_reference operator* () const {
			return (*vector)[offset];
		}
    
		inline iterator& operator++ () {
			++offset; return *this;
		}

		inline iterator operator++ (int) {
			iterator before = *this;
			++offset; return before;
		}

		inline iterator& operator-- () {
			--offset; return *this;
		}

		inline iterator operator-- (int) {
			iterator before = *this;
			--offset; return before;
		}

		inline bool operator== (const iterator &I) const {
			SEQAN_ASSERT(vector == I.vector);
			return offset == I.offset;
		}

		inline bool operator!= (const iterator &I) const {
			SEQAN_ASSERT(vector == I.vector);
			return offset != I.offset;
		}

		inline bool operator< (const iterator &I) const {
			SEQAN_ASSERT(vector == I.vector);
			return offset < I.offset;
		}		
	};


	//////////////////////////////////////////////////////////////////////////////
	// forward vector iterator
	template < typename _Vector >
    struct VectorFwdIterator : public ::std::iterator <
                                    ::std::bidirectional_iterator_tag,
                                    typename _Vector::value_type,
                                    typename _Vector::size_type,
                                    typename _Vector::pointer,
                                    typename _Vector::reference >
	{
		typedef _Vector		                        Vector;   
		typedef VectorFwdIterator			        iterator;
        typedef typename Vector::Type			    Type;
		typedef typename Vector::SizeType		    SizeType;
		typedef typename Vector::const_reference	const_reference;
		typedef typename Vector::volatile_ptr	    volatile_ptr;
        enum { _PageSize = Vector::PageSize };

//        friend class Vector;
        friend struct VectorIterator<Vector>;
        typedef VectorIterator<Vector>              std_iterator;
        typedef VectorConstIterator<Vector>         std_const_iterator;

		//////////////////////////////////////////////////////////////////////////////
		// public iterator interface

		// in fact this is also a random iterator, but you better use
		// the VectorIterator class for *real* random access
//		typedef std::bidirectional_iterator_tag		iterator_category;
		typedef std::random_access_iterator_tag		iterator_category;
		typedef typename Vector::value_type			value_type;
		typedef typename Vector::difference_type	difference_type;
		typedef typename Vector::pointer			pointer;
		typedef typename Vector::reference			reference;
		
		Vector			*vector;

        bool            dirty;
		int     		pageNo;
		unsigned		pageOfs;
        int             prefetch;   // -n .. prefetch n pages downwards, n .. prefetch n pages upwards, 0 .. disabled
		volatile_ptr	begin;
		
	    VectorFwdIterator():
			vector(NULL),
			pageNo(0),
			pageOfs(0),
            prefetch(0),
			begin(NULL) {}

		VectorFwdIterator(const iterator &I):
			vector(I.vector),
			pageNo(I.pageNo),
			pageOfs(I.pageOfs),
            prefetch(I.prefetch),
			begin(NULL) {}

	    explicit VectorFwdIterator(Vector *_vector, SizeType _offset):
			vector(_vector),
			pageNo(_offset / _PageSize),
			pageOfs(_offset % _PageSize),
            prefetch(0),
			begin(NULL) {}

		explicit VectorFwdIterator(Vector *_vector, SizeType _pageNo, SizeType _pageOfs):
			vector(_vector),
			pageNo(_pageNo),
			pageOfs(_pageOfs),
            prefetch(0),
			begin(NULL) {}

		~VectorFwdIterator() {
			invalidate();
		}

		//////////////////////////////////////////////////////////////////////////////
		// iterator conversion interface

		VectorFwdIterator(const std_iterator &I):
			vector(I.vector),
			pageNo(I.offset / _PageSize),
            pageOfs(I.offset % _PageSize),
            prefetch(0),
			begin(NULL) {}

		inline iterator& operator=(std_iterator const & _Right) {
			invalidate();
			pageNo = _Right.offset / _PageSize;
			pageOfs = _Right.offset % _PageSize;
            vector = _Right.vector;
			return *this;
		}

        inline operator std_iterator() const {
            return std_iterator(vector, (SizeType)pageNo * (SizeType)_PageSize + pageOfs);
        }

        inline operator std_const_iterator() const {
            return std_const_iterator(vector, (SizeType)pageNo * (SizeType)_PageSize + pageOfs);
        }

		inline iterator& operator=(iterator const & _Right) {
			invalidate();
			vector = _Right.vector;
			pageNo = _Right.pageNo;
			pageOfs = _Right.pageOfs;
            prefetch = _Right.prefetch;
			return *this;
		}

		//////////////////////////////////////////////////////////////////////////////
		// iterator arithmetic

		inline difference_type operator- (const iterator &I) const {
			return (difference_type)(pageNo - I.pageNo) * (difference_type)_PageSize + (pageOfs - I.pageOfs);
		};
		
		inline iterator operator- (difference_type delta) const {
			difference_type dPNo  = delta / _PageSize;
			difference_type dPOfs = delta % _PageSize;
			if (pageOfs >= dPOfs)
				return iterator(vector, pageNo - dPNo, pageOfs - dPOfs);
			else
				return iterator(vector, pageNo - dPNo - 1, _PageSize + pageOfs - dPOfs);
		};
		
		inline iterator& operator-= (difference_type delta) {
			difference_type dPNo  = delta / _PageSize;
			difference_type dPOfs = delta % _PageSize;
			if (pageOfs < dPOfs) {
				++dPNo;
				pageOfs = _PageSize + pageOfs - dPOfs;
			} else
				pageOfs -= dPOfs;
			if (dPNo) invalidate(0);
			pageNo -= dPNo;
			return *this;
		};
		
		inline iterator operator+ (difference_type delta) const {
			difference_type dPNo  = delta / _PageSize;
			difference_type nPOfs = pageOfs + delta % _PageSize;
			if (nPOfs < _PageSize)
				return iterator(vector, pageNo + dPNo, nPOfs);
			else
				return iterator(vector, pageNo + dPNo + 1, nPOfs - _PageSize);
		};
		
		inline iterator& operator+= (difference_type delta) {
			difference_type dPNo  = delta / _PageSize;
			difference_type nPOfs = pageOfs + delta % _PageSize;
			if (nPOfs >= _PageSize) {
				++dPNo;
				nPOfs -= _PageSize;
			}
			if (dPNo) invalidate(0);
			pageNo += dPNo;
			pageOfs = nPOfs;
			return *this;
		};
		
		inline void validate() const {
			typename Vector::PageFrameRef pf = vector->getSharedPage(pageNo, prefetch);
            const_cast<iterator*>(this)->dirty = pf.dirty;
			const_cast<iterator*>(this)->begin = pf.begin;
		}

        inline void invalidate(int _prefetch = 0) const {
            if (begin) {
                const_cast<iterator*>(this)->begin = NULL;
				vector->releasePage(pageNo, (prefetch != 0) || (_prefetch != 0));
                const_cast<iterator*>(this)->prefetch = _prefetch;
            }
		}

		inline reference operator* () {
			if (!begin) validate();
            // synchronize PageFrame dirty flag on dirty false->true change
            if (!dirty) {
                dirty = true;
    			vector->getPage(pageNo).dirty = true;
            }
			return begin[pageOfs];
		}
    
		inline const_reference operator* () const {
			if (!begin) validate();
			return begin[pageOfs];
		}
    
		inline iterator& operator++ () {
			if (++pageOfs == _PageSize) {
				invalidate(1);
				pageOfs = 0;
				++pageNo;
			}
			return *this;
		}

		inline iterator operator++ (int) {
			iterator before = *this;
			if (++pageOfs == _PageSize) {
				invalidate(1);
				pageOfs = 0;
				++pageNo;
			}
			return before;
		}

		inline iterator& operator-- () {
			if (pageOfs)
				--pageOfs;
			else {
				invalidate(-1);
				pageOfs = _PageSize - 1;
				--pageNo;
			}
			return *this;
		}

		inline iterator operator-- (int) {
			iterator before = *this;
			if (pageOfs)
				--pageOfs;
			else {
				invalidate(-1);
				pageOfs = _PageSize - 1;
				--pageNo;
			}
			return before;
		}

		inline bool operator== (const iterator &I) const {
			SEQAN_ASSERT(vector == I.vector);
			return pageNo == I.pageNo && pageOfs == I.pageOfs;
		}

		inline bool operator!= (const iterator &I) const {
			SEQAN_ASSERT(vector == I.vector);
			return pageNo != I.pageNo || pageOfs != I.pageOfs;
		}

		inline bool operator< (const iterator &I) const {
			SEQAN_ASSERT(vector == I.vector);
			return pageNo < I.pageNo || (pageNo == I.pageNo && pageOfs < I.pageOfs);
		}
    };


	//////////////////////////////////////////////////////////////////////////////
	// const forward vector iterator
	template < typename _Vector >
    struct VectorFwdConstIterator : public ::std::iterator <
                                    ::std::bidirectional_iterator_tag,
                                    typename _Vector::value_type,
                                    typename _Vector::size_type,
                                    typename _Vector::pointer,
                                    typename _Vector::reference >
	{
		typedef _Vector		                        Vector;   
		typedef VectorFwdConstIterator		        iterator;
        typedef typename Vector::Type			    Type;
		typedef typename Vector::SizeType		    SizeType;
		typedef typename Vector::const_reference	const_reference;
		typedef typename Vector::volatile_ptr	    volatile_ptr;
        enum { _PageSize = Vector::PageSize };

//        friend class _Vector;
        friend struct VectorIterator<Vector>;
        friend struct VectorConstIterator<Vector>;
        typedef VectorIterator<Vector>              std_iterator;
        typedef VectorConstIterator<Vector>         std_const_iterator;
        typedef VectorFwdIterator<Vector>			fwd_iterator;

		//////////////////////////////////////////////////////////////////////////////
		// public iterator interface

		// in fact this is also a random iterator, but you better use
		// the VectorIterator class for *real* random access
//		typedef std::bidirectional_iterator_tag		iterator_category;
		typedef std::random_access_iterator_tag		iterator_category;
		typedef typename Vector::value_type			value_type;
		typedef typename Vector::difference_type	difference_type;
		typedef typename Vector::const_pointer		pointer;
		typedef typename Vector::const_reference	reference;
		
		Vector			*vector;

		int     		pageNo;
		unsigned		pageOfs;
        int             prefetch;   // -n .. prefetch n pages downwards, n .. prefetch n pages upwards, 0 .. disabled
		volatile_ptr	begin;
		

        VectorFwdConstIterator():
			vector(NULL),
			pageNo(0),
			pageOfs(0),
            prefetch(0),
			begin(NULL) {}

		VectorFwdConstIterator(const iterator &I):
			vector(I.vector),
			pageNo(I.pageNo),
			pageOfs(I.pageOfs),
            prefetch(I.prefetch),
			begin(NULL) {}

		VectorFwdConstIterator(const fwd_iterator &I):
			vector(I.vector),
			pageNo(I.pageNo),
			pageOfs(I.pageOfs),
            prefetch(I.prefetch),
			begin(NULL) {}

		~VectorFwdConstIterator() {
			invalidate();
		}

	    VectorFwdConstIterator(Vector *_vector, SizeType _offset):
			vector(_vector),
			pageNo(_offset / _PageSize),
			pageOfs(_offset % _PageSize),
            prefetch(0),
			begin(NULL) {}

		VectorFwdConstIterator(Vector *_vector, SizeType _pageNo, SizeType _pageOfs):
			vector(_vector),
			pageNo(_pageNo),
			pageOfs(_pageOfs),
            prefetch(0),
			begin(NULL) {}

		//////////////////////////////////////////////////////////////////////////////
		// iterator conversion interface

		VectorFwdConstIterator(std_iterator &I):
			vector(I.vector),
			pageNo(I.offset / _PageSize),
            pageOfs(I.offset % _PageSize),
            prefetch(0),
			begin(NULL) {}

        VectorFwdConstIterator(std_const_iterator const &I):
			vector(I.vector),
			pageNo(I.offset / _PageSize),
            pageOfs(I.offset % _PageSize),
            prefetch(0),
			begin(NULL) {}

		inline iterator& operator=(std_iterator const & _Right) {
			invalidate();
			pageNo = _Right.offset / _PageSize;
			pageOfs = _Right.offset % _PageSize;
            vector = _Right.vector;
			return *this;
		}

		inline iterator& operator=(std_const_iterator const & _Right) {
			invalidate();
			pageNo = _Right.offset / _PageSize;
			pageOfs = _Right.offset % _PageSize;
            vector = _Right.vector;
			return *this;
		}

        inline operator std_const_iterator() const {
            return std_const_iterator(vector, (SizeType)pageNo * (SizeType)_PageSize + pageOfs);
        }

		inline iterator& operator=(iterator const & _Right) {
			invalidate();
			vector = _Right.vector;
			pageNo = _Right.pageNo;
			pageOfs = _Right.pageOfs;
            prefetch = _Right.prefetch;
			return *this;
		}

		inline iterator& operator=(fwd_iterator const & _Right) {
			invalidate();
			vector = _Right.vector;
			pageNo = _Right.pageNo;
			pageOfs = _Right.pageOfs;
            prefetch = _Right.prefetch;
			return *this;
		}

		//////////////////////////////////////////////////////////////////////////////
		// iterator arithmetic

		inline difference_type operator- (const iterator &I) const {
			return (difference_type)(pageNo - I.pageNo) * (difference_type)_PageSize + (pageOfs - I.pageOfs);
		};
		
		inline iterator operator- (difference_type delta) const {
			difference_type dPNo  = delta / _PageSize;
			difference_type dPOfs = delta % _PageSize;
			if (pageOfs >= dPOfs)
				return iterator(vector, pageNo - dPNo, pageOfs - dPOfs);
			else
				return iterator(vector, pageNo - dPNo - 1, _PageSize + pageOfs - dPOfs);
		};
		
		inline iterator& operator-= (difference_type delta) {
			difference_type dPNo  = delta / _PageSize;
			difference_type dPOfs = delta % _PageSize;
			if (pageOfs < dPOfs) {
				++dPNo;
				pageOfs = _PageSize + pageOfs - dPOfs;
			} else
				pageOfs -= dPOfs;
			if (dPNo) invalidate(0);
			pageNo -= dPNo;
			return *this;
		};
		inline iterator operator+ (difference_type delta) const {
			difference_type dPNo  = delta / _PageSize;
			difference_type nPOfs = pageOfs + delta % _PageSize;
			if (nPOfs < _PageSize)
				return iterator(vector, pageNo + dPNo, nPOfs);
			else
				return iterator(vector, pageNo + dPNo + 1, nPOfs - _PageSize);
		};
		
		inline iterator& operator+= (difference_type delta) {
			difference_type dPNo  = delta / _PageSize;
			difference_type nPOfs = pageOfs + delta % _PageSize;
			if (nPOfs >= _PageSize) {
				++dPNo;
				nPOfs -= _PageSize;
			}
			if (dPNo) invalidate(0);
			pageNo += dPNo;
			pageOfs = nPOfs;
			return *this;
		};
		
		inline void validate() const {
			typename Vector::PageFrameRef pf = vector->getSharedPage(pageNo, prefetch);
			const_cast<iterator*>(this)->begin = pf.begin;
		}

        inline void invalidate(int _prefetch = 0) const {
            if (begin) {
                const_cast<iterator*>(this)->begin = NULL;
				vector->releasePage(pageNo, (prefetch != 0) || (_prefetch != 0));
                const_cast<iterator*>(this)->prefetch = _prefetch;
            }
		}

		inline const_reference operator* () const {
			if (!begin) validate();
			return begin[pageOfs];
		}
    
		inline iterator& operator++ () {
			if (++pageOfs == _PageSize) {
				invalidate(1);
				pageOfs = 0;
				++pageNo;
			}
			return *this;
		}

		inline iterator operator++ (int) {
			iterator before = *this;
			if (++pageOfs == _PageSize) {
				invalidate(1);
				pageOfs = 0;
				++pageNo;
			}
			return before;
		}

		inline iterator& operator-- () {
			if (pageOfs)
				--pageOfs;
			else {
				invalidate(-1);
				pageOfs = _PageSize - 1;
				--pageNo;
			}
			return *this;
		}

		inline iterator operator-- (int) {
			iterator before = *this;
			if (pageOfs)
				--pageOfs;
			else {
				invalidate(-1);
				pageOfs = _PageSize - 1;
				--pageNo;
			}
			return before;
		}

		inline bool operator== (const iterator &I) const {
			SEQAN_ASSERT(vector == I.vector);
			return pageNo == I.pageNo && pageOfs == I.pageOfs;
		}

		inline bool operator!= (const iterator &I) const {
			SEQAN_ASSERT(vector == I.vector);
			return pageNo != I.pageNo || pageOfs != I.pageOfs;
		}

		inline bool operator< (const iterator &I) const {
			SEQAN_ASSERT(vector == I.vector);
			return pageNo < I.pageNo || (pageNo == I.pageNo && pageOfs < I.pageOfs);
		}

	};



    //////////////////////////////////////////////////////////////////////////////
    // (Prototypes)
    //////////////////////////////////////////////////////////////////////////////

    template < typename TValue,
               typename TConfig >
	class String<TValue, External<TConfig> >
	{

	public:
        enum { _PageSize = TConfig::PageSize,
               _Frames   = TConfig::Frames,
               PageSize  = TConfig::PageSize,
			   TempByDefault = TConfig::TempByDefault };

		typedef TValue	                    Type;
        typedef typename TConfig::TFile     TFile;
        typedef typename TConfig::SizeType  SizeType;

		typedef std::vector<int>						        PageTable;
		typedef PageFrame<TValue, TFile, Fixed<_PageSize> >  	PageFrame;
		typedef PageContainer<PageFrame, _Frames>		        Cache;
		typedef PageFrame&								        PageFrameRef;

    public: // debug
		PageTable			pager;
		Cache				cache;
		TFile				file;
        bool                _temporary, _ownFile;
		SizeType			_size;
        int                 lastDiskPage;       // the last page on disk and in mem 
        unsigned            lastDiskPageSize;   // can be smaller than PageSize
        bool keepFirst;                         // true .. try to keep the lowest page frames

    public:

		//////////////////////////////////////////////////////////////////////////////
		// public iterator types

        friend struct VectorIterator<String>;
        friend struct VectorConstIterator<String>;
        friend struct VectorFwdIterator<String>;
        friend struct VectorFwdConstIterator<String>;

		typedef VectorIterator<String>					VectorIterator;
		typedef VectorConstIterator<String>		        VectorConstIterator;
		typedef VectorFwdIterator<String>				VectorFwdIterator;
		typedef VectorFwdConstIterator<String>	        VectorFwdConstIterator;

		typedef VectorFwdIterator		iterator;
		typedef VectorFwdConstIterator	const_iterator;
		typedef Type					value_type;
		typedef Type&					reference;
		typedef Type const &			const_reference;
		typedef Type*					pointer;
		typedef Type const *			const_pointer;
		typedef SizeType				size_type;
		typedef SizeType				difference_type;
		typedef VolatilePtr<Type>       volatile_ptr;

		String(SizeType __size = 0):
            file(NULL),
			_size(0)
        {
            _temporary = true;
            _ownFile = false;
            lastDiskPage = 0;       // actually, these values need not to be initialized
            lastDiskPageSize = 0;   // here, because of "write before read"

			if (TempByDefault)
				if (!openTemp())
					::std::cout << "External String couldn't open temporary file" << ::std::endl;

			resize(__size);
            keepFirst = false;
        }

		String(TFile &_file)
        {
			open(_file);
            keepFirst = false;
        }

		String(const char *fileName, int openMode = OPEN_RDWR + OPEN_CREATE | OPEN_APPEND):
			file(NULL)
        {
			open(fileName, openMode);
            keepFirst = false;
        }

		~String() {
			close();
		}

		//////////////////////////////////////////////////////////////////////////////
		// vector interface

        inline void clear() {
            keepFirst = false;
            pager.clear();
            resize(0);
        }

		inline size_type length() const {
			return _size;
		}

		inline size_type capacity() const {
			return (size_type)pager.capacity() * (size_type)_PageSize;
		}

		inline void reserve(size_type _newCapacity) {
			pager.reserve(enclosingBlocks(_newCapacity, (unsigned)_PageSize));
		}

		inline void resize(size_type _newSize) {
			pager.resize(enclosingBlocks(_newSize, (unsigned)_PageSize), PageFrame::UNINITIALIZED);
            if (_newSize < _size && file) {
                waitForAll();                                       // wait for all pending transfers
                ::seqan::resize(file, (size_type)_newSize * (size_type)sizeof(TValue));   // before shrinking the file size
                lastDiskPage = _newSize / _PageSize;
                lastDiskPageSize = _newSize % _PageSize;
            }
			_size = _newSize;
		}

		inline reference operator[] (size_type offset) {
			PageFrameRef pf = getPage(offset / _PageSize);
			pf.dirty = true;
			return pf[offset % _PageSize];
		}

		inline const_reference operator[] (size_type offset) const {
			return const_cast<String*>(this)->getPage(offset / _PageSize)[offset % _PageSize];
		}

		inline iterator begin() {
			return iterator(this, 0);
		}

		inline const_iterator begin() const {
			return const_iterator(const_cast<String*>(this), 0);
		}

		inline iterator end() {
			return iterator(this, _size);
		}

		inline const_iterator end() const {
            return const_iterator(const_cast<String*>(this), _size);
		}

		inline void push(const_reference obj) {
			resize(_size + 1);
			back() = obj;
		}

		inline void push_back(const_reference obj) {
			resize(_size + 1);
			back() = obj;
		}

		inline void pop_back()	{
			resize(_size - 1);
		}

		inline reference front() {
			return (*this)[0];
		}

		inline const_reference front() const {
			return (*this)[0];
		}

		inline reference back() {
			return (*this)[_size - 1];
		}

		inline const_reference back() const {
			return (*this)[_size - 1];
		}

	    template <typename TSource>
	    inline String & operator= (TSource const & source)
	    {
		    assign(*this, source);
		    return *this;
	    }

	    inline String & operator= (String const & source)
	    {
		    assign(*this, source);
		    return *this;
	    }

        inline operator bool() {
            return file;
        }

    protected:

		//////////////////////////////////////////////////////////////////////////////
		// swapping interface

        void _dumpCache() {
            for(int i = 0; i < cache.size(); ++i) {
                PageFrameRef pf = cache[i];
                std::cout << "[" << pf.pageNo << "]";
                if (pf.dirty)
                    std::cout << "*";
                else
                    std::cout << " ";

                if (pf.status == PageFrame::READY)
                    std::cout << "   ";
                else
                    std::cout << ".  ";
            }
            std::cout << std::endl;
        }


        // return a priority for a page frame (the higher is more persistent)
        inline typename PageFrame::Priority getPriority(int pageNo) const {
            if (keepFirst && pageNo < cache.size() - 10) // save 1 for random access
                return PageFrame::PERMANENT_LEVEL;
            else
                return PageFrame::NORMAL_LEVEL;
        }

		// write page to disk if dirty and remove from page table now or after finishing IO
		inline void flush(PageFrameRef pf) {
            if (pf.status == PageFrame::READY && pf.dirty) {    // write if dirty and not i/o transferring
				nukeCopies(pf.begin);				            // proceeding writes should wait and set dirty bit

                if (pf.priority > PageFrame::NORMAL_LEVEL && pf.priority <= PageFrame::ITERATOR_LEVEL)
					cache.upgrade(pf, PageFrame::PREFETCH_LEVEL);

                if (pf.pageNo != _size / _PageSize)
    				writePage(pf, pf.pageNo, file);
                else {
                    lastDiskPage = _size / _PageSize;
                    lastDiskPageSize = _size % _PageSize;;
				    writeLastPage(pf, pf.pageNo, file, lastDiskPageSize);
                }
                pf.dataStatus = PageFrame::ON_DISK;
			}
		}

		// write page synchronously to disk if dirty and remove from page table
		inline void swapOutAndWait(PageFrameRef pf) {
			nukeCopies(pf.begin);      				// proceeding writes should wait and set dirty bit

            if (pf.status != PageFrame::READY) {            
				pager[pf.pageNo] = PageFrame::ON_DISK;		// page is not dirty and on disk
				waitFor(pf);                                // after finishing i/o transfer
                pf.pageNo = -1;                             // cut back link
                return;
            }

			if (pf.dirty) {                                 // write if dirty
                if (pf.pageNo != _size / _PageSize) {
    				writePage(pf, pf.pageNo, file);
                    if (pf.pageNo >= lastDiskPage)
                        lastDiskPage = -1;       			// make lastDiskPage(Size) invalid because file size is aligned
                } else {
				    writeLastPage(pf, pf.pageNo, file, _size % _PageSize);
                    lastDiskPage = _size / _PageSize;
                    lastDiskPageSize = _size % _PageSize;;
                }
				pager[pf.pageNo] = PageFrame::ON_DISK;		// page is marked to be on disk
				waitFor(pf);
			} else
				pager[pf.pageNo] = pf.dataStatus;			// restore original data status

            pf.pageNo = -1;                                 // cut back link
		}

		// wait until IO of every page is finished
		void waitForAll() {
			for(typename Cache::iterator I = cache.begin(); I != cache.end(); ++I)
                waitFor(*I);
		}

		struct testIODone : public std::unary_function<PageFrameRef,bool> {
			String &me;
			testIODone(String &_me): me(_me) {}

			inline bool operator() (PageFrameRef pf) {
                if (waitFor(pf, 0)) {
                    if (pf.pageNo >= me.lastDiskPage)
                        me.lastDiskPage = -1;    // make lastDiskPage(Size) invalid because file size is aligned
                    return true;
                } else
                    return false;
			}
		};

        inline PageFrameRef getPage(
            int pageNo, 
            typename PageFrame::Priority maxLevel = PageFrame::NORMAL_LEVEL, 
            typename PageFrame::Priority newLevel = PageFrame::NORMAL_LEVEL,
            int prefetchPages = 0)
        {
			int frameNo = pager[pageNo];
			if (frameNo >= 0) {					// cache hit

				PageFrameRef pf = cache[frameNo];
				cache.upgrade(
                    pf, 
                    Max(pf.priority, newLevel));    		// update lru order

				if (waitFor(pf))    						// wait for i/o transfer to complete
                    if (pf.pageNo >= lastDiskPage) {
                        lastDiskPage = -1;       			// make lastDiskPage(Size) invalid because file size is aligned
                    }

                if (prefetchPages > 0) prefetch(pageNo + 1, pageNo + 1 + prefetchPages, frameNo);
                else if (prefetchPages < 0) prefetch(pageNo + prefetchPages, pageNo, frameNo);

				return pf;

			} else {							// cache miss

				typename PageFrame::DataStatus dataStatus = static_cast<typename PageFrame::DataStatus>(frameNo);
				frameNo = cache.mru(testIODone(*this), maxLevel);   // try to get an undirty and READY pageframe
				if (frameNo < 0)							// if there is none,
					frameNo = cache.mruDirty();				// get the most recently used dirty frame
				PageFrameRef pf = cache[frameNo];

				// *** frame is choosen ***

				if (pf.begin)
					swapOutAndWait(pf);						// write synchronously to disk, if page is dirty
				else
					allocPage(pf, file);                    // allocate memory if page is virgin

				// *** frame is free now ***

				pf.dataStatus = dataStatus;
				if (dataStatus == PageFrame::ON_DISK)
                    if (pageNo != lastDiskPage)
					    readPage(pageNo, pf, file);
                    else
                        readLastPage(pageNo, pf, file, lastDiskPageSize);
				
				pager[pageNo] = frameNo;					// assign new page to page table
				pf.pageNo = pageNo;							// set back link
				cache.upgrade(
                    pf,
                    Max(getPriority(pageNo), newLevel));    // update lru order

                if (prefetchPages > 0) prefetch(pageNo + 1, pageNo + 1 + prefetchPages, frameNo);
                else if (prefetchPages < 0) prefetch(pageNo + prefetchPages, pageNo, frameNo);
                
                waitFor(pf);    							// wait for i/o transfer to complete
				return pf;
			}
		}

    public:

        // prefetch is non-blocking and should speed up swapping
		inline void prefetch(int pageBegin, int pageEnd, int except = -1) {
            if (!file) return;
            if (pageBegin < 0)              pageBegin = 0;
            if (pageEnd >= pager.size())    pageEnd = pager.size() - 1;
            for(int pageNo = pageBegin; pageNo < pageEnd; ++pageNo) {
			    int frameNo = pager[pageNo];
				typename PageFrame::DataStatus dataStatus = static_cast<typename PageFrame::DataStatus>(frameNo);
                if (dataStatus == PageFrame::ON_DISK &&             // prefetch only if page is on disk
                    pageNo != lastDiskPage)                         // reading the last page is blocking
                {   
				    frameNo = cache.mru(
                        testIODone(*this),
                        PageFrame::NORMAL_LEVEL);                   // choose undirty and ready page

                    if (frameNo < 0 || frameNo == except) return;   // no lowlevel-page left for prefetching
				    PageFrameRef pf = cache[frameNo];
                    #ifdef SEQAN_VERBOSE
						::std::cout << "prefetch: page " << pageNo << ::std::endl;
                    #endif

                    // *** frame is choosen ***

				    if (pf.begin)
					    swapOutAndWait(pf);						    // write synchronously to disk, if page is dirty
				    else
					    allocPage(pf, file);                        // allocate memory if page is virgin

    				// *** frame is free now ***

    				pf.dataStatus = dataStatus;
                    readPage(pageNo, pf, file);
				    pager[pageNo] = frameNo;					    // assign new page to page table
    				pf.pageNo = pageNo;							    // set back link
                    cache.upgrade(pf, PageFrame::PREFETCH_LEVEL);   // update lru order
                }
            }
		}
		
	    template < typename _TDefault >
		inline static int _prefetchIffAsync(int prefetchPages, _TDefault const &) {
			return 0;
		}
		
	    template < typename _TConfig >
		inline static int _prefetchIffAsync(int prefetchPages, File<Async<_TConfig> > const &) {
			return prefetchPages;
		}

		inline PageFrameRef getSharedPage(int pageNo, int prefetchPages = 0) {
			return getPage(
                pageNo, 
                PageFrame::PREFETCH_LEVEL, 
                PageFrame::ITERATOR_LEVEL,
                _prefetchIffAsync(prefetchPages, file));
		}

		inline void releasePage(int pageNo, bool writeThrough = false) {
            if (pageNo==4)
                pageNo=4;
			int frameNo = pager[pageNo];
			if (frameNo >= 0) {								        // release only cached pages
				PageFrameRef pf = cache[frameNo];
				if (pf.begin.isLonely() && pf.priority <= PageFrame::ITERATOR_LEVEL) {
					cache.upgrade(pf, Max(getPriority(pageNo), PageFrame::NORMAL_LEVEL));
                    if (writeThrough) {
                        #ifdef SEQAN_VERBOSE
                            if (pf.dirty)
								::std::cout << "writeThrough: page " << pageNo << ::std::endl;
                        #endif
					    flush(pf);							        // write if dirty
                    }
				}
			}
		}
        
		// cancel all transactions
		inline void cancel() {
			if (file)
			    for(typename Cache::iterator I = cache.begin(); I != cache.end(); ++I)
					if (I->begin) ::seqan::cancel(*I, file);
		}

		// write all dirty pages to disk
		inline void flush() {
			if (file) {
			    for(typename Cache::iterator I = cache.begin(); I != cache.end(); ++I)
				    if (I->begin) flush(*I);
				waitForAll();
			}
		}

		// flush and free all allocated pages
		inline void free() {
            flush();
            for(typename Cache::iterator I = cache.begin(); I != cache.end(); ++I) {
				if (I->pageNo >= 0) {
                    pager[I->pageNo] = I->dataStatus;
					I->pageNo = PageFrame::UNINITIALIZED;
				}
//                ::std::cout << *I << ::std::endl;
                if (I->begin) freePage(*I, file);
            }
		}

		inline bool open(const char *fileName, int openMode = OPEN_RDWR + OPEN_CREATE | OPEN_APPEND) {
            _temporary = false;
			if (_ownFile = SEQAN_NAMESPACE_MAIN::open(file, fileName, openMode))
                _size = size(file) / sizeof(TValue);
            else
                _size = 0;
			pager.resize(enclosingBlocks(_size, (unsigned)_PageSize), (_size)? PageFrame::ON_DISK: PageFrame::UNINITIALIZED);
            lastDiskPage = _size / _PageSize;
            lastDiskPageSize = _size % _PageSize;
			return _ownFile;
		}

		inline bool open(TFile &_file) {
			file = _file;
            _temporary = false;
            _ownFile = false;
            if (file)
                _size = size(file) / sizeof(TValue);
            else
                _size = 0;
            pager.resize(enclosingBlocks(_size, (unsigned)_PageSize), (_size)? PageFrame::ON_DISK: PageFrame::UNINITIALIZED);
            lastDiskPage = _size / _PageSize;
            lastDiskPageSize = _size % _PageSize;
			return _file;
        }

		inline bool openTemp() {
            _temporary = true;
            lastDiskPage = 0;
            lastDiskPageSize = 0;
			pager.clear();
			return _ownFile = SEQAN_NAMESPACE_MAIN::openTemp(file);
		}

		// close associated file
		inline bool close() {
			if (_temporary)	cancel();
			else			flush();
			free();
			pager.clear();
			if (_ownFile) {
				_ownFile = false;
				return SEQAN_NAMESPACE_MAIN::close(file);
			} else
				return true;
		}

        inline void rename(int frameNo) {
			PageFrameRef pf = cache[frameNo];
            cache.rename(frameNo);                                  // update lru entry
            if (pf.pageNo >= 0)
                pager[pf.pageNo] = frameNo;					        // update back link
        }

        // change the number of in-mem pageframes
        // more pages mean less swapping, 
        // less pages mean more free mem
        inline void resizeCache(unsigned newFrames) {
            unsigned oldFrames = cache.size();
            if (_size)
                newFrames = Min(newFrames, enclosingBlocks(_size, (unsigned)_PageSize));
            if (newFrames < oldFrames) {
                flush();
                for(unsigned i = newFrames; i < oldFrames; ++i) {
    			    int frameNo = cache.mruDirty();             // get the most recently used frame (can be dirty)
                    if (frameNo < 0) break;
				    PageFrameRef pf = cache[frameNo];

				    // *** frame is choosen ***

                    if (pf.begin) {
					    swapOutAndWait(pf);						// write synchronously to disk, if page is dirty
        				freePage(pf, file);                     // free memory
                    }

                    cache.erase(frameNo);                       // erase page frame from cache

                    for(int j = frameNo; j < cache.size(); ++j)
                        rename(j);                              // update remaining pages
                }
            } else if (oldFrames < newFrames) {
                cache.resize(newFrames);
            }
        }

    };


    //////////////////////////////////////////////////////////////////////////////
	// handler that manages a simple memory buffer
/*    template < typename TValue,
               typename TConfig >
	struct BufferHandler< Pipe< String<TValue, External<TConfig> >, Source<ContainerSpec> > >
    {
        typedef TValue                                                      Type;
        typedef typename Size< String<TValue, External<TConfig> > >::Type   SizeType;
        typedef SimpleBuffer<TValue, SizeType>                              SimpleBuffer;

        typedef Pipe< String<TValue, External<TConfig> >, Source<ContainerSpec> > Pipe;

		Pipe			&pipe;
        int             pageNo;

		BufferHandler(Pipe &_pipe):
			pipe(_pipe) {}

        inline SimpleBuffer& begin() {
            return pipe.in.getPage(pageNo = 0);
        }

        inline SimpleBuffer& next() {
            return pipe.in.getPage(++pageNo);
        }

        inline void process() {}
        inline void end() {}
        inline void cancel() {}
    };*/

	template < typename TValue,
               typename TConfig >
	struct BufferHandler< Pipe< String<TValue, External<TConfig> >, Source<ContainerSpec> > >
    {
        typedef TValue                                                      Type;
        typedef typename Size< String<TValue, External<TConfig> > >::Type   SizeType;
        typedef SimpleBuffer<TValue, SizeType>                              Buffer;

        typedef Pipe< String<TValue, External<TConfig> >, Source<ContainerSpec> >   Pipe;
        typedef typename String<TValue, External<TConfig> >::VectorFwdConstIterator ISource;
		typedef typename Iterator<Buffer>::Type										ITarget;

		Pipe		&pipe;
		unsigned	bufferSize;
        SizeType    rest;
        Buffer		buffer;
		ISource		source;

		BufferHandler(Pipe &_pipe, unsigned requestedSize):
			pipe(_pipe),
			bufferSize(requestedSize),
            rest(0) {}

        inline Buffer& first() {
            rest = length(pipe.in);
			allocPage(buffer, Min(bufferSize, rest), *this);
			source = begin(pipe.in);
			for(ITarget target = buffer.begin; target != buffer.end; ++target) {
				*target = *source;
				++source;
			}
            if (!(rest -= size(buffer))) source = ISource();
			return buffer;
        }

        inline Buffer& next() {
			resize(buffer, Min(bufferSize, rest));
			ITarget _end = buffer.begin + size(buffer);
			for(ITarget target = buffer.begin; target != _end; ++target) {
				*target = *source;
				++source;
			}
            if (!(rest -= size(buffer))) source = ISource();
			return buffer;
        }

        inline void process() {}
        inline void end() { cancel(); }
        inline void cancel() { source = ISource(); freePage(buffer, *this); }
    };

    template < typename TValue, typename TConfig >
    struct Value< BufferHandler< Pipe< String<TValue, External<TConfig> >, Source<ContainerSpec> > > > {
        typedef SimpleBuffer< TValue > Type;
    };


    //////////////////////////////////////////////////////////////////////////////
    // global interface
    //////////////////////////////////////////////////////////////////////////////

    template < typename TValue, typename TConfig >
    struct Value< String<TValue, External<TConfig> > >
    {
	    typedef typename String<TValue, External<TConfig> >::Type Type;
    };

    template < typename TValue, typename TConfig >
    struct Size< String<TValue, External<TConfig> > >
    {
        typedef typename String<TValue, External<TConfig> >::SizeType Type;
    };

    template < typename TValue, typename TConfig >
    struct Position< String<TValue, External<TConfig> > >
    {
        typedef typename String<TValue, External<TConfig> >::SizeType Type;
    };

    template < typename TValue, typename TConfig >
    struct Difference< String<TValue, External<TConfig> > >
    {
		typedef typename _MakeSigned<typename String<TValue, External<TConfig> >::SizeType>::Type Type;
    };

    template < typename TValue, typename TConfig >
    struct Iterator< String<TValue, External<TConfig> > const, Standard >
    {
        typedef typename String<TValue, External<TConfig> >::const_iterator Type;
    };

    template < typename TValue, typename TConfig >
    struct Iterator< String<TValue, External<TConfig> >, Standard >
    {
        typedef typename String<TValue, External<TConfig> >::iterator Type;
    };

    template < typename TValue, typename TConfig >
	struct DefaultIteratorSpec< String<TValue, External<TConfig> > > {
		typedef Standard Type;
	};
	
    template < typename TValue, typename TConfig >
	struct DefaultIteratorSpec< String<TValue, External<TConfig> > const > {
		typedef Standard Type;
	};
	

	// iterator metafunctions
	template < typename TVector >
    struct Value< VectorIterator<TVector> >			{ typedef typename Value<TVector>::Type Type; };
	template < typename TVector >
    struct Value< VectorConstIterator<TVector> >	{ typedef typename Value<TVector>::Type Type; };
	template < typename TVector >
    struct Value< VectorFwdIterator<TVector> >		{ typedef typename Value<TVector>::Type Type; };
	template < typename TVector >
    struct Value< VectorFwdConstIterator<TVector> > { typedef typename Value<TVector>::Type Type; };

	template < typename TVector >
	struct Reference< VectorConstIterator<TVector> >:
		public Reference< typename Value<TVector>::Type const > {};

	template < typename TVector >
	struct Reference< VectorFwdConstIterator<TVector> >:
		public Reference< typename Value<TVector>::Type const > {};

	template < typename TVector >
    struct Size< VectorIterator<TVector> >			{ typedef typename Size<TVector>::Type Type; };
	template < typename TVector >
    struct Size< VectorConstIterator<TVector> >		{ typedef typename Size<TVector>::Type Type; };
	template < typename TVector >
    struct Size< VectorFwdIterator<TVector> >		{ typedef typename Size<TVector>::Type Type; };
	template < typename TVector >
    struct Size< VectorFwdConstIterator<TVector> >	{ typedef typename Size<TVector>::Type Type; };

	template < typename TVector >
    struct Position< VectorIterator<TVector> >		{ typedef typename Position<TVector>::Type Type; };
	template < typename TVector >
    struct Position< VectorConstIterator<TVector> >	{ typedef typename Position<TVector>::Type Type; };
	template < typename TVector >
    struct Position< VectorFwdIterator<TVector> >	{ typedef typename Position<TVector>::Type Type; };
	template < typename TVector >
    struct Position< VectorFwdConstIterator<TVector> > { typedef typename Position<TVector>::Type Type; };

	template < typename TVector >
    struct Difference< VectorIterator<TVector> >		{ typedef typename Difference<TVector>::Type Type; };
	template < typename TVector >
    struct Difference< VectorConstIterator<TVector> >	{ typedef typename Difference<TVector>::Type Type; };
	template < typename TVector >
    struct Difference< VectorFwdIterator<TVector> >		{ typedef typename Difference<TVector>::Type Type; };
	template < typename TVector >
    struct Difference< VectorFwdConstIterator<TVector> > { typedef typename Difference<TVector>::Type Type; };

    template < typename TValue, typename TConfig >
    inline void 
    clear(String<TValue, External<TConfig> > &me) {
        me.clear();
    }

    template < typename TValue, typename TConfig >
    inline void 
    flush(String<TValue, External<TConfig> > &me) {
        me.flush();
    }

	template < typename TValue, typename TConfig >
    inline bool 
    open(String<TValue, External<TConfig> > &me, const char *fileName, int openMode = OPEN_RDWR + OPEN_CREATE | OPEN_APPEND) {
		return me.open(fileName, openMode);
    }

	template < typename TValue, typename TConfig >
    inline bool 
    open(String<TValue, External<TConfig> > &me, typename TConfig::TFile file) {
		return me.open(file);
    }

	template < typename TValue, typename TConfig >
    inline bool 
    openTemp(String<TValue, External<TConfig> > &me) {
		return me.open();
    }

	template < typename TValue, typename TConfig >
    inline bool 
    save(String<TValue, External<TConfig> > &me, const char *fileName, int openMode = OPEN_RDWR + OPEN_CREATE | OPEN_APPEND) {
		// External Strings are persistent, thus there is no need to save them
		//ExtStringsDontNeedToBeSaved error;
	}

	template < typename TValue, typename TConfig >
    inline bool 
    save(String<TValue, External<TConfig> > &me, typename TConfig::TFile file) {
		// External Strings are persistent, thus there is no need to save them
		//ExtStringsDontNeedToBeSaved error;
	}

    template < typename TValue, typename TConfig >
    inline bool 
    close(String<TValue, External<TConfig> > &me) {
		return me.close();
    }



    template < typename TValue, typename TConfig >
    inline typename Size< String<TValue, External<TConfig> > >::Type
    size(String<TValue, External<TConfig> > &me)
    {
        return me.length();
    }

    template < typename TValue, typename TConfig >
    inline typename Size< String<TValue, External<TConfig> > >::Type
    length(String<TValue, External<TConfig> > &me)
    {
        return me.length();
    }

    template < typename TValue, typename TConfig, typename TSize >
    inline typename Size< String<TValue, External<TConfig> > >::Type
    resize(
	    String<TValue, External<TConfig> > &me,
	    TSize new_length)
    {
	    me.resize(new_length);
        return me.length();
    }

    template < typename TValue, typename TConfig, typename TSize, typename TExpand >
    inline typename Size< String<TValue, External<TConfig> > >::Type
    resize(
	    String<TValue, External<TConfig> > &me,
		TSize new_length,
		Tag<TExpand> const tag)
	{
		return resize(me, new_length);
	}

/*
    template < typename TValue, typename TConfig, typename TSize>
    inline typename Size< String<TValue, External<TConfig> > >::Type
    resize(
	    String<TValue, External<TConfig> > &me,
	    TSize new_length,
		TagLimit)
    {
	    me.resize(new_length);
        return me.length();
    }
*/
    template < typename TValue, typename TConfig, typename TSpec >
    inline typename Iterator<String<TValue, External<TConfig> >, Tag<TSpec> const>::Type
    begin(String<TValue, External<TConfig> > &me, Tag<TSpec> const) {
		return me.begin();
    }

    template < typename TValue, typename TConfig, typename TSpec >
    inline typename Iterator<String<TValue, External<TConfig> > const, Tag<TSpec> const>::Type
    begin(String<TValue, External<TConfig> > const &me, Tag<TSpec> const) {
		return me.begin();
    }

    template < typename TValue, typename TConfig, typename TSpec >
    inline typename Iterator<String<TValue, External<TConfig> >, Tag<TSpec> const>::Type
    end(String<TValue, External<TConfig> > &me, Tag<TSpec> const) {
		return me.end();
    }

    template < typename TValue, typename TConfig, typename TSpec >
    inline typename Iterator<String<TValue, External<TConfig> > const, Tag<TSpec> const>::Type
    end(String<TValue, External<TConfig> > const &me, Tag<TSpec> const) {
		return me.end();
    }

    template < typename TValue, typename TConfig, typename TPos >
    inline typename Reference<String<TValue, External<TConfig> > >::Type 
    at(String<TValue, External<TConfig> > &me, TPos pos)
    {
	    return me[pos];
    }

    template < typename TValue, typename TConfig >
    inline void
    push_back(String<TValue, External<TConfig> > &me, TValue const &_Val)
    {
	    return me.push_back(_Val);
    }
/*
    template < typename TSpec >
	std::ostream& operator<<(std::ostream &out, String<char, TSpec > &p) {

        typename Iterator< String<char, TSpec > >::Type _cur = begin(p), _end = end(p);
        while (_cur != _end) {
		    out << *_cur;
            ++_cur;
        }
		return out;
	}

    template < typename TValue, typename TSpec >
	std::ostream& operator<<(std::ostream &out, String<TValue, TSpec > &p) {

        typename Iterator< String<TValue, TSpec > >::Type _cur = begin(p), _end = end(p);
        while (_cur != _end) {
		    out << *_cur << " ";
            ++_cur;
        }
		return out;
	}
*/
    // sequence -> external string
    template < typename TValue,
               typename TConfig,
               typename TSource >
    inline void assign(String<TValue, External<TConfig> > &target, TSource const &source) {

        typedef typename Iterator<TSource const>::Type                          ISource;
        typedef typename String<TValue, External<TConfig> >::VectorFwdIterator  ITarget;

        resize(target, length(source));

        ISource it_source       = begin(source);
        ITarget it_target       = begin(target);
		ITarget it_target_end   = end(target);
		while (it_target != it_target_end)
		{
			*it_target = *it_source;
			++it_target;
			++it_source;
		}
    }

    template < typename TValue, typename TConfig >
    inline void const * 
    id(String<TValue, External<TConfig> > const &me)
    {
        return &(*me.pager.begin());
    }

}

#endif
