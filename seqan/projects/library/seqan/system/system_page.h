/*
 *  system_page.h
 *  genindex
 *
 *  Created by David Weese on 17.07.05.
 *
 */

#ifndef SEQAN_HEADER_SYSTEM_PAGE_H
#define SEQAN_HEADER_SYSTEM_PAGE_H

#include <cassert>
#include <limits>
#include <cstdio>
#include <list>
#include <vector>

//////////////////////////////////////////////////////////////////////////////

namespace SEQAN_NAMESPACE_MAIN
{


//////////////////////////////////////////////////////////////////////////////

	template <typename T1, typename T2> inline
	T1 enclosingBlocks(T1 _size, T2 _blockSize) {
		return (_size + _blockSize - 1) / _blockSize;
	}

	template <typename T1, typename T2> inline
	T1 alignSize(T1 _size, T2 _aligning) {
        if (_size < _aligning)
            return _aligning;
        else
		    return (_size / _aligning) * (T1)_aligning;
	}

	//////////////////////////////////////////////////////////////////////////////
	// simple class for fast random accessable containers buffers
    template < typename TIterator, typename TValue, typename TSize >
	struct IteratorBuffer {
        typedef TValue      Type;
        typedef TSize       SizeType;
		typedef TIterator	Iterator;

        TIterator           begin;
        TIterator           end;

        IteratorBuffer():
            begin(TIterator()),
            end(TIterator()) {}

		IteratorBuffer(TIterator _begin, TIterator _end):
            begin(_begin),
            end(_end) {}

        IteratorBuffer(TIterator _begin, SizeType _size):
            begin(_begin),
            end(_begin + _size) {}

//        inline Type& operator[](SizeType i) { return begin[i]; }
        inline Type const & operator[](SizeType i) const { return begin[i]; }
    };

    template < typename TIterator, typename TValue, typename TSize >
    struct Iterator< IteratorBuffer<TIterator, TValue, TSize> >
    {
        typedef TIterator Type;
    };

    template < typename TIterator, typename TValue, typename TSize >
    struct Value< IteratorBuffer<TIterator, TValue, TSize> >
    {
        typedef TValue Type;
    };

    template < typename TIterator, typename TValue, typename TSize >
    struct Size< IteratorBuffer<TIterator, TValue, TSize> >
    {
        typedef TSize Type;
    };

    template < typename TIterator, typename TValue, typename TSize >
    inline typename Size< IteratorBuffer<TIterator, TValue, TSize> >::Type
    size(IteratorBuffer<TIterator, TValue, TSize> const &me) {
        return me.end - me.begin;
    }

    template < typename TIterator, typename TValue, typename TSize >
    inline typename Size< IteratorBuffer<TIterator, TValue, TSize> >::Type
    length(IteratorBuffer<TIterator, TValue, TSize> const &me) {
        return me.end - me.begin;
    }


	//////////////////////////////////////////////////////////////////////////////
	// base class for memory buffers
    template < typename TValue, typename TSize = unsigned >
	struct SimpleBuffer {
        typedef TValue      Type;
        typedef TSize       SizeType;

		typedef	TValue&		TypeRef;
		typedef TValue*     TypePtr;
		typedef TValue*     Iterator;

		Iterator            begin;      // the beginning of the buffer
        Iterator            end;        // end of valid data
        TSize               pageSize;   // size of allocated memory

        SimpleBuffer():
            begin(NULL),
            end(NULL) {}

        SimpleBuffer(TypePtr _begin, TypePtr _end):
            begin(_begin),
            end(_end) {}

        SimpleBuffer(TypePtr _begin, SizeType _size):
            begin(_begin),
            end(_begin + _size) {}

        SimpleBuffer(SizeType _pageSize):
            begin(NULL),
            end(NULL),
            pageSize(_pageSize) {}

        inline Type& operator[](SizeType i) { return begin[i]; }
        inline Type const & operator[](SizeType i) const { return begin[i]; }
	};

    template < typename TValue >
    struct Iterator< SimpleBuffer<TValue> >
    {
        typedef TValue* Type;
    };

    template < typename TValue >
    struct Value< SimpleBuffer<TValue> >
    {
        typedef TValue Type;
    };

    template < typename TValue, typename TSize >
    struct Size< SimpleBuffer<TValue, TSize> >
    {
        typedef TSize Type;
    };

    template < typename TValue >
    inline typename Size<SimpleBuffer<TValue> >::Type
    pageSize(SimpleBuffer<TValue> &me) {
        return me.pageSize;
    }

    template < typename TValue, typename TSize >
    inline void setPageSize(SimpleBuffer<TValue> &me, TSize size) {
        me.pageSize = size;
    }

    template < typename TValue >
    inline typename Size<SimpleBuffer<TValue> >::Type
    size(SimpleBuffer<TValue> const &me) {
        return me.end - me.begin;
    }

    template < typename TValue >
    inline typename Size<SimpleBuffer<TValue> >::Type
    length(SimpleBuffer<TValue> const &me) {
        return me.end - me.begin;
    }

    template < typename TValue, typename TSize >
    inline void resize(SimpleBuffer<TValue> &me, TSize size) {
        me.end = me.begin + size;
    }

    template < typename TValue, typename TSize, typename T >
	inline void allocPage(SimpleBuffer<TValue> &pf, TSize size, T const & me) {
        setPageSize(pf, size);
        allocate(me, pf.begin, pageSize(pf));
        resize(pf, size);
	}

	template < typename TValue, typename T > inline
	void freePage(SimpleBuffer<TValue> &pf, T const & me) {
		deallocate(me, pf.begin, pageSize(pf));
		pf.begin = NULL;
        resize(pf, 0);
        setPageSize(pf, 0);
	}


    //////////////////////////////////////////////////////////////////////////////
	// a bucket is a structure to represent a small window of a page
    // used by algorithms which need a global view of all pages (merge sort, mapper)
    template < typename TValue >
    struct PageBucket {
        unsigned    pageOfs;                // begin of bucket window with relation to page begin
        TValue  	*begin, *cur, *end;     // begin/end of buckets memory buffer and a pointer
    };

    template < typename TValue >
    struct PageBucketExtended : public PageBucket< TValue > {
		int     	pageNo;		            // related page (needed by merger sort)
    };

	template < typename TValue >
    ::std::ostream& operator<<(::std::ostream &out, const PageBucketExtended<TValue> &pb) {
        for(TValue *cur = pb.begin; cur != pb.end; cur++)
            out << *cur << " ";
        return out;
    }


    template < typename TValue, typename TFile, typename TSpec >
    struct PageFrame {};    

	//////////////////////////////////////////////////////////////////////////////
	// page frame of dynamic size

    template < typename TSpec = void >
    struct Dynamic;

    // forward declaration
    template < typename TPageFrame >
    struct PageChain;

    template < typename TValue, typename TFile >
	struct PageFrame< TValue, TFile, Dynamic<> >: public SimpleBuffer< TValue >
	{
        typedef TValue                          Type;
        typedef unsigned                        SizeType;
		typedef TFile							File;

		typedef	TValue&							TypeRef;
		typedef TValue*           	            TypePtr;
        typedef PageChain<PageFrame>		    PageChain;
		typedef SimpleBuffer<TValue>	        Base;
        typedef typename aRequest<TFile>::Type  aRequest;

		bool			dirty;		// data needs to be written to disk before freeing
		unsigned   		pageNo;		// maps frames to pages (reverse vector mapper)
        aRequest        request;    // request structure of the async io process

        enum Status		{ READY, READING, WRITING };
		Status status;

        PageFrame       *next;      // next buffer in a chained list
//        PageChain       *chain;     // related chain 

        PageFrame(/*BufChain *_chain = NULL*/):
            dirty(false),
            pageNo(-1),
			status(READY),
            next(NULL),
//            chain(_chain),
			Base() { }
    };
    

	//////////////////////////////////////////////////////////////////////////////
	// page frame of static size

    template < unsigned _PageSize >
    struct Fixed;

    typedef ::std::list<int>		PageLRUList;    // least recently usage list
	typedef PageLRUList::iterator	PageLRUEntry;

    template < typename TValue,
               typename TFile,
               unsigned _PageSize >
	struct PageFrame<TValue, TFile, Fixed<_PageSize> >
	{
        typedef TValue                          Type;
        typedef unsigned                        SizeType;
		typedef TFile							File;
		enum { PageSize = _PageSize };

		typedef	TValue&							TypeRef;
		typedef VolatilePtr<TValue>	            TypePtr;
        typedef typename aRequest<TFile>::Type	aRequest;

		bool			dirty;		// data needs to be written to disk before freeing
		int     		pageNo;		// maps frames to pages (reverse vector mapper)
		TypePtr			begin;	    // start address of page memory
        aRequest        request;    // request structure of the async io process

		enum Status		{ READY, READING, WRITING };
		enum DataStatus	{ ON_DISK = -1, UNINITIALIZED = -2 };
		enum Priority	{ NORMAL_LEVEL = 0, PREFETCH_LEVEL = 1, ITERATOR_LEVEL = 2, PERMANENT_LEVEL = 3 };

		Status status;
		DataStatus dataStatus;

		PageLRUEntry	lruEntry;   // priority based lru
        Priority        priority;

		PageFrame():
            dirty(false),
            pageNo(-1),
			begin(NULL),
			status(READY),
            priority(NORMAL_LEVEL) {}

        inline Type& operator[](SizeType i) { return begin[i]; }
        inline Type const & operator[](SizeType i) const { return begin[i]; }
	};

    template < typename TValue, typename TFile, typename TSize >
    inline void resize(PageFrame<TValue, TFile, Dynamic<> > &me, TSize size) {
        me.end = me.begin + size;
    }

    template < typename TValue, typename TFile, unsigned _PageSize >
    inline typename Size<PageFrame<TValue, TFile, Fixed<_PageSize> > >::Type
    size(PageFrame<TValue, TFile, Fixed<_PageSize> > &me) {
        return _PageSize;
    }

    template < typename TValue, typename TFile, unsigned _PageSize >
    inline typename Size<PageFrame<TValue, TFile, Fixed<_PageSize> > >::Type
    length(PageFrame<TValue, TFile, Fixed<_PageSize> > &me) {
        return _PageSize;
    }

    template < typename TValue, typename TFile, unsigned _PageSize >
    inline typename Size<PageFrame<TValue, TFile, Fixed<_PageSize> > >::Type
    pageSize(PageFrame<TValue, TFile, Fixed<_PageSize> > &me) {
        return _PageSize;
    }

    template < typename TValue, typename TFile, unsigned _PageSize, typename TSize >
    inline void resize(PageFrame<TValue, TFile, Fixed<_PageSize> > &me, TSize size) { }

#ifndef qsort_r
    
    static void *qsort_context;

    template < typename TValue, typename TCompare >
    static int compQSort(const void *a, const void *b) {
        TCompare &C = *reinterpret_cast<TCompare*>(qsort_context);
        return C(*reinterpret_cast<const TValue*>(a), *reinterpret_cast<const TValue*>(b));
    }
    
    template < typename TValue, typename TCompare >
    static void quickSort(SimpleBuffer<TValue> &buf, TCompare &C) {
        // poor windows qsort doesn't support contexts and thus no reentrance
        qsort_context = &C;
        qsort(buf.begin, size(buf), sizeof(TValue), compQSort<TValue, TCompare>);
    }

#else

    template < typename TValue, typename TCompare >
    static int compQSort_r(void *p, const void *a, const void *b) {
        TCompare &C = *reinterpret_cast<TCompare*>(C);
        return C(*reinterpret_cast<const TValue*>(a), *reinterpret_cast<const TValue*>(b));
    }
    
    template < typename TValue, typename TCompare >
    static void quickSort(SimpleBuffer<TValue> &buf, TCompare &C) {
        qsort_r(buf.begin, size(buf), sizeof(TValue), &C, compQSort_r<TValue, TCompare>);
    }

#endif


	//////////////////////////////////////////////////////////////////////////////
	// various page frame methods

	template < typename TValue, typename TFile, typename TSpec >
    ::std::ostream& operator<<(::std::ostream &out, const PageFrame<TValue, TFile, TSpec > &pf) {
        out << "PageFrame @ " << pf.pageNo;
        if (pf.dirty)
            out << " DIRTY";
        else
            out << " CLEAN";

        switch (pf.status) {
			case PageFrame<TValue, TFile, TSpec >::READY:
                out << " READY";
                break;
			case PageFrame<TValue, TFile, TSpec >::READING:
                out << " READING";
                break;
			case PageFrame<TValue, TFile, TSpec >::WRITING:
                out << " WRITING";
        }

        if (pf.dataStatus == pf.ON_DISK)
            out << " ON_DISK";
        else
            out << " UNITIALIZED";

        out << " Prio:" << (int)pf.priority;
        out << " Buffer:" << (unsigned)(TValue*)pf.begin;

        return out;
	}

	template < typename TValue, typename TFile, typename TSpec, typename T > inline
	void allocPage(PageFrame<TValue, TFile, TSpec> &pf, T const & me) {
		TValue* tmp = NULL;
		allocate(me, tmp, pageSize(pf));
		pf.begin = tmp;
		#ifdef SEQAN_VVERBOSE
			printf("allocPage: %x\n", (unsigned)tmp);
		#endif
	}

	template < typename TValue, typename TFile, typename TSpec, typename T > inline
	void freePage(PageFrame<TValue, TFile, TSpec> &pf, T const & me) {
		#ifdef SEQAN_VVERBOSE
            if ((TValue*)pf.begin)
			    printf("freePage:  %x\n", (unsigned)(TValue*)pf.begin);
		#endif
        nukeCopies(pf.begin);
		deallocate(me, (TValue*)pf.begin, pageSize(pf));
		pf.begin = NULL;
        resize(pf, 0);
	}

	template < typename TValue, typename TFile, typename TSpec > inline
	bool readPage(int pageNo, PageFrame<TValue, TFile, TSpec> &pf, TFile &file) {
		typedef typename Position<TFile>::Type pos_t;
		#ifdef SEQAN_VVERBOSE
			printf("readPage:  %x from page %d\n", (unsigned)(TValue*)pf.begin, pageNo);
		#endif
		pf.dirty = false;
		pf.status = pf.READING;
//        resize(pf, pageSize(pf));
		return areadAt(file, (TValue*)pf.begin, size(pf), (pos_t)pageNo * (pos_t)pageSize(pf), pf.request);
	}

	template < typename TValue, typename TFile, typename TSpec > inline
	bool writePage(PageFrame<TValue, TFile, TSpec> &pf, int pageNo, TFile &file) {
		typedef typename Position<TFile>::Type pos_t;
		#ifdef SEQAN_VVERBOSE
			printf("writePage: %x from page %d\n", (unsigned)(TValue*)pf.begin, pageNo);
		#endif
		pf.status = pf.WRITING;
//        resize(pf, pageSize(pf));
		return awriteAt(file, (TValue*)pf.begin, size(pf), (pos_t)pageNo * (pos_t)pageSize(pf), pf.request);
	}

	template < typename TValue, typename TFile, typename TSpec, typename TSize> inline
    bool readLastPage(int pageNo, PageFrame<TValue, TFile, TSpec> &pf, TFile &file, TSize size) {
		typedef typename Position<TFile>::Type pos_t;
		#ifdef SEQAN_VVERBOSE
			printf("readPage:  %x from page %d size %d\n", (unsigned)(TValue*)pf.begin, pageNo, size);
		#endif
		pf.dirty = false;
		pf.status = pf.READY;
//        resize(pf, size);
		return readAt(file, (TValue*)pf.begin, size, (pos_t)pageNo * (pos_t)pageSize(pf));
	}

	template < typename TValue, typename TFile, typename TSpec, typename TSize > inline
	bool writeLastPage(PageFrame<TValue, TFile, TSpec> &pf, int pageNo, TFile &file, TSize size) {
		typedef typename Position<TFile>::Type pos_t;
		#ifdef SEQAN_VVERBOSE
			printf("writePage: %x from page %d size %d\n", (unsigned)(TValue*)pf.begin, pageNo, size);
		#endif
		pf.dirty = false;
		pf.status = pf.READY;
//        resize(pf, size);
		return writeAt(file, (TValue*)pf.begin, size, (pos_t)pageNo * (pos_t)pageSize(pf));
	}

	template < typename TValue, typename TFile, typename TSpec > inline
	bool waitFor(PageFrame<TValue, TFile, TSpec> &pf) {
		if ((pf.status != pf.READY) && waitFor(pf.request)) {
			pf.status = pf.READY;
			pf.dirty = false;
            return true;
		}
        return false;
	}

	template < typename TValue, typename TFile, typename TSpec, typename TTime > inline
	bool waitFor(PageFrame<TValue, TFile, TSpec> &pf, TTime timeOut) {
		if ((pf.status != pf.READY) && waitFor(pf.request, timeOut)) {
			pf.status = pf.READY;
			pf.dirty = false;
			return true;
		}
        return false;
	}

	template < typename TValue, typename TFile, typename TSpec > inline
	bool cancel(PageFrame<TValue, TFile, TSpec> &pf, TFile &file) {
        waitFor(pf, 0);
		if (pf.status != pf.READY) {
            if (!cancel(file, pf.request)) return false;
            pf.status = pf.READY;
        }
        return true;
	}

	//////////////////////////////////////////////////////////////////////////////
	// page based read/write methods used by Pool classes

    template < typename TValue, typename TFile, typename TSpec > inline
	bool readPage(PageFrame<TValue, TFile, Dynamic<TSpec> > &pf, TFile &file) {
        if (size(pf) == pageSize(pf))
            return readPage(pf.pageNo, pf, file);
        else
            return readLastPage(pf.pageNo, pf, file, size(pf));
	}

	template < typename TValue, typename TFile, typename TSpec > inline
	bool writePage(PageFrame<TValue, TFile, Dynamic<TSpec> > &pf, TFile &file) {
        if (size(pf) == pageSize(pf))
            return writePage(pf, pf.pageNo, file);
        else
            return writeLastPage(pf, pf.pageNo, file, size(pf));
	}

	template < typename TValue, typename TFile > inline
	unsigned readBucket(PageBucket<TValue> &b, int pageNo, unsigned pageSize, unsigned dataSize, TFile &file) {
		typedef typename Position<TFile>::Type pos_t;
        unsigned readSize = Min(dataSize - b.pageOfs, (unsigned)(b.end - b.begin));
		#ifdef SEQAN_VVERBOSE
			printf("readBucket:  %x from page %d at %d size %d\n", (unsigned)b.begin, pageNo, pageNo * pageSize + b.pageOfs, readSize);
		#endif
        if (readSize && readAt(file, b.begin, readSize, (pos_t)pageNo * (pos_t)pageSize + b.pageOfs)) {
            b.pageOfs += readSize;
            b.cur = b.begin;
            b.end = b.begin + readSize;
            return readSize;
        } else
            return 0;
	}

	template < typename TValue, typename TFile > inline
	bool writeBucket(PageBucket<TValue> &b, int pageNo, unsigned pageSize, TFile &file) {
		typedef typename Position<TFile>::Type pos_t;
		#ifdef SEQAN_VVERBOSE
			printf("writeBucket:  %x from page %d at %d size %d\n", (unsigned)b.begin, pageNo, pageNo * pageSize + b.pageOfs, b.cur - b.begin);
		#endif
        if ((b.cur == b.begin) || writeAt(file, b.begin, b.cur - b.begin, (pos_t)pageNo * (pos_t)pageSize + b.pageOfs)) {
            b.pageOfs += b.cur - b.begin;
            b.cur = b.begin;
            return true;
        } else
            return false;
	}

	template < typename TValue, typename TFile, typename TSpec > inline
	bool writeBucket(PageFrame<TValue, TFile, Dynamic<TSpec> > &pf, unsigned &pageOfs, TFile &file) {
		typedef typename Position<TFile>::Type pos_t;
		#ifdef SEQAN_VVERBOSE
			printf("writeBucket:  %x from page %d at %d size %d\n", (unsigned)pf.begin, pf.pageNo, pf.pageNo * pageSize(pf) + pageOfs, size(pf));
		#endif
        if (pf.end == pf.begin) return true;
        if (awriteAt(file, pf.begin, size(pf), (pos_t)pf.pageNo * (pos_t)pageSize(pf) + pageOfs, pf.request)) {
            pf.status = pf.WRITING;
            pageOfs += size(pf);
            return true;
        } else
            return false;
	}



    template < typename TPageFrame >
    struct PageChain {
        typedef typename TPageFrame::Type		Type;
        typedef typename TPageFrame::SizeType	SizeType;

		typedef	Type&							TypeRef;
		typedef Type*							TypePtr;
		typedef TPageFrame						PageFrame;

		//////////////////////////////////////////////////////////////////////////////
		// public iterator interface

		typedef TPageFrame * 		iterator;
		typedef TPageFrame const *	const_iterator;
		typedef TPageFrame      	value_type;
		typedef TPageFrame &     	reference;
		typedef TPageFrame const &  const_reference;
		typedef TPageFrame * 		pointer;
		typedef TPageFrame const *  const_pointer;
		typedef unsigned       		size_type;
		typedef int                 difference_type;

        pointer             first, last;
/*        pointer             firstFree, lastFree;
        Semaphore           freeBuffers;
        Mutex               dataLock;*/

        unsigned            frames, maxFrames;
        
        PageChain(unsigned _maxFrames = UINT_MAX):
            first(NULL),
            last(NULL),
            frames(0),
            maxFrames(_maxFrames)
//            freeBuffers(_buffers)
        {
            for(unsigned i = 0; i < _maxFrames; ++i)
                pushBack();
        }
        
        ~PageChain()
        {
            while (first)
                popFront();
        }
        
        inline reference operator[](int k) {
            pointer p = first;
            while (k) {
                p = p->next;
                --k;
            }
            return *p;
        }
        
        inline const_reference operator[](int k) const {
            pointer p = first;
            while (k) {
                p = p->next;
                --k;
            }
            return *p;
        }

        inline pointer getReadyPage() {
            if (!first || (frames < maxFrames && waitFor(*first, 0)))
                return pushBack();
            else {
                waitFor(*first);
                return firstToEnd();
            }
        }

        template < typename TFile >
        inline void cancelAll(TFile &file) {
            pointer p = first;
            while (p) {
                cancel(*p, file);
                p = p->next;
            }
        }

        inline void waitForAll() {
            pointer p = first;
            while (p) {
                waitFor(*p);
                p = p->next;
            }
        }

/*        BufHeader* getFreeBuffer() {
            freeBuffers.lock();
            dataLock.lock();
            BufHeader *h = firstFree;
            if (!(firstFree = firstFree->next)) lastFree = NULL;
            dataLock.unlock();
            return h;
        }
        
        void addFreeBuffer(BufHeader *h) {
            dataLock.lock();
            if (lastFree) {
                lastFree->next = h;
                lastFree = h;
            } else {
                firstFree = h;
                lastFree = h;
            }
            dataLock.unlock();
            freeBuffers.unlock();
        }*/
        
    private:

        inline pointer firstToEnd() {
            last->next = first;
            last = first;
            first = first->next;
            last->next = NULL;
            return last;
        }

        inline pointer pushBack() {
            pointer p = new PageFrame();
            if (p) {
                if (last)
                    last->next = p;
                else
                    first = p;
                last = p;
                ++frames;
            }
            return p;
        }

        inline pointer popFront() {
            pointer p = first;
            if (p) {
                first = first->next;
                if (!first) last = NULL;
                --frames;
                delete p;
            }
            return p;
        }
    };


/*    
    template < typename ST >
    struct BufArrayConfig
    {
        typedef ST SizeType;

        unsigned int    disks;
        unsigned int    buffersPerDisk;
        SizeType        bufferSize;         // prefetch/cache buffer size
        
        BufArrayConfig():
            disks(1),
            buffersPerDisk(3),
            bufferSize(2 * 1024 * 1024) { }
    };
    
    template < typename T, typename ST, typename _File >
    struct BufArray
    {
        typedef T Type;
        typedef ST SizeType;
        typedef BufChain<T,ST,_File> BufChain;
        typedef BufArrayConfig<ST> Config;
        
        Config		conf;
        BufChain**	disk;
        
        BufArray(const Config &_conf):
            conf(_conf)
        {
            disk = new BufChain*[conf.disks];
            for(unsigned int i = 0; i < conf.disks; i++)
                disk[i] = new BufChain(conf.buffersPerDisk, conf.bufferSize, i);
            BufChain *prev = disk[conf.disks - 1];
            for(unsigned int i = 0; i < conf.disks; i++)
                prev = prev->next = disk[i];
        }

        ~BufArray() {
            for(unsigned int i = 0; i < conf.disks; i++)
                delete disk[i];
            delete[] disk;
        }

        BufChain& operator[](const int k) {
            return *disk[k];
        }

        const BufChain& operator[](const int k) const {
            return *disk[k];
        }
    };*/



	//////////////////////////////////////////////////////////////////////////////
	// page container with lru mechanism
	// the lru strategy uses different priorities
	// the page with the least priority is used first
	// 0..random access pages
	// 1..forward iterator pages
	// 2..quasi permanent pages

    template < typename TPageFrame,
               unsigned _Frames,
               unsigned _PriorityLevels = TPageFrame::PERMANENT_LEVEL + 1 >
	struct PageContainer : public ::std::vector<TPageFrame>
	{
        typedef typename TPageFrame::Type			Type;
        typedef typename TPageFrame::SizeType		SizeType;
		typedef typename ::std::vector<TPageFrame>	Base;

		typedef	Type&							TypeRef;
		typedef Type*							TypePtr;
		typedef TPageFrame						PageFrame;

		enum { PriorityLevels = _PriorityLevels };

		//////////////////////////////////////////////////////////////////////////////
		// public iterator interface

		typedef typename Base::iterator			iterator;
		typedef typename Base::const_iterator	const_iterator;
		typedef TPageFrame      				value_type;
		typedef TPageFrame &     				reference;
		typedef TPageFrame const &				const_reference;
		typedef TPageFrame * 					pointer;
		typedef TPageFrame const *				const_pointer;
		typedef unsigned       					size_type;
		typedef int								difference_type;
		
        PageLRUList lruList[_PriorityLevels];

		PageContainer():
			Base(_Frames)
		{
            for(size_type i = 0; i < _Frames; ++i)
			    (*this)[i].lruEntry = lruList[0].insert(lruList[0].end(), i);
        }

		inline void reserve(unsigned _Count) {
			Base::reserve(_Count);
		}

		void resize(unsigned _Count) {
			unsigned _Size = Base::size();
			if (_Size < _Count)
                for(unsigned i = _Size; i < _Count; ++i)
                    push_back();
			else 
				if (_Size > _Count)
					for(unsigned i = _Count; i < _Size; ++i)
						pop_back();
		}

		inline void push_back() {
			Base::push_back(PageFrame());
			Base::back().lruEntry = lruList[0].insert(lruList[0].end(), Base::size() - 1);
		}

		inline void erase(int frameNo) {
			lruList[(*this)[frameNo].priority].erase((*this)[frameNo].lruEntry);
            Base::erase(Base::begin() + frameNo);
		}

        inline void rename(int frameNo) {
            *((*this)[frameNo].lruEntry) = frameNo;
        }

		inline void pop_back() {
			lruList[Base::back().priority].erase(Base::back().lruEntry);
			Base::pop_back();
		}


		//////////////////////////////////////////////////////////////////////////////
		// lru strategy interface

		inline void upgrade(const PageFrame &pf) {
			lruList[pf.priority].splice(lruList[pf.priority].begin(), lruList[pf.priority], pf.lruEntry);
			pf.lruEntry = lruList[pf.priority].begin();
		}

		inline void downgrade(const PageFrame &pf) {
			lruList[pf.priority].splice(lruList[pf.priority].end(), lruList[pf.priority], pf.lruEntry);
			pf.lruEntry = lruList[pf.priority].end() - 1;
		}

		inline void upgrade(PageFrame &pf, int newPriority) {
			lruList[newPriority].splice(lruList[newPriority].begin(), lruList[pf.priority], pf.lruEntry);
			pf.lruEntry = lruList[newPriority].begin();
			pf.priority = static_cast<typename PageFrame::Priority> (newPriority);
		}

		inline void downgrade(PageFrame &pf, int newPriority) {
			lruList[newPriority].splice(lruList[newPriority].end(), lruList[pf.priority], pf.lruEntry);
			pf.lruEntry = lruList[newPriority].end() - 1;
			pf.priority = static_cast<typename PageFrame::Priority> (newPriority);
		}

		inline void _dump() {
			for(unsigned i = 0; i < _PriorityLevels; ++i) {
                std::cout << "|";
                PageLRUList::const_iterator I = lruList[i].end();
                PageLRUList::const_iterator first = lruList[i].begin();
				while (I != first) {
					--I;
                    PageFrame &pf = (*this)[*I];
                    std::cout << pf.pageNo;
                    if (pf.dirty) std::cout << "*";
                    else          std::cout << " ";
                    if (pf.status == PageFrame::READY) std::cout << "  ";
                    else                               std::cout << ". ";
				};
            }
            std::cout << std::endl;
		}

        // Function is a functor which is called with a PageFrame object,
        // that is dirty or not READY (in an IO transfer)
		template <class Function>
		inline int mru(Function _Func, unsigned maxLevel = _PriorityLevels - 1) {
			for(unsigned i = 0; i <= maxLevel; ++i) {
                PageLRUList::const_iterator I = lruList[i].end();
                PageLRUList::const_iterator first = lruList[i].begin();
				while (I != first) {
					--I;
					PageFrame& pf =(*this)[*I];
					if (pf.status == PageFrame::READY && !pf.dirty)
						return *I;
					else
						if (_Func(pf)) return *I;
				};
            }
			#ifdef SEQAN_VVERBOSE
				printf("ALL PAGES DIRTY OR IN USE (try to use const iterators) :-(\n");
			#endif
			return -1;
		}

		inline int mruDirty() {
            for(unsigned i = 0; i < _PriorityLevels; ++i)
                if (!lruList[i].empty())
                    return lruList[i].back();
			return -1;
		}
	};




    template < typename TValue, typename TSize, typename T, class Function >
    inline bool equiDistantDistribution(
        SimpleBuffer<TValue> &_clusterBuffer, unsigned _bufferSize, T const &me,
        TSize _size, unsigned _pageSize,
        Function const &_Func)
    {
        unsigned _pages         = enclosingBlocks(_size, (unsigned)_pageSize);
        if (!_pages) {
            printf("equiDistantDistribution: _pages is null!\n");
            return false;
        }

        if (_bufferSize < _pages) {
            printf("equiDistantDistribution: clusterBufferSize is too small -> raised to %d\n", _pages);
            _bufferSize = _pages;
        }

        unsigned lastPageSize   = _size % _pageSize;
        unsigned pages          = _pages;

        if (_bufferSize > _size)
            _bufferSize = _size;

        allocPage(_clusterBuffer, _bufferSize, me);
        PageBucketExtended<TValue> pb;
        pb.begin = _clusterBuffer.begin;

        unsigned clusterSize = _bufferSize / pages;
        if (lastPageSize > 0 && clusterSize >= lastPageSize) {
            // last page bucket would get more memory than page would need
            // --> exclude from equi-size distribution
            if (--pages) {
                _bufferSize -= lastPageSize;
                clusterSize = _bufferSize / pages;
            }
        }

        if (pages) {
            unsigned remainder = _bufferSize % pages;
            for(unsigned i = 0, numerator = 0; i < pages; ++i) {
                pb.end = pb.begin + clusterSize;
                if ((numerator += remainder) >= pages) {    // simple bresenham for distribution
                    numerator -= pages;
                    ++pb.end;
                }
                pb.cur = pb.begin;
                pb.pageOfs = 0;
			    _Func(pb);
                pb.begin = pb.end;
            }
        }

        if (pages < _pages) {
            pb.end = pb.begin + lastPageSize;
            pb.cur = pb.begin;
            pb.pageOfs = 0;
			_Func(pb);
        }

        return true;
    }

    template < typename TValue, typename TSize, typename T, class Function >
    inline unsigned equiDistantAlignedDistribution(
        SimpleBuffer<TValue> &_clusterBuffer, unsigned aligning, unsigned _bufferSize, T const &me,
        TSize _size, unsigned _pageSize,
        Function const &_Func)
    {
        unsigned _pages         = enclosingBlocks(_size, (unsigned)_pageSize);
        if (!_pages) {
            printf("equiDistantDistribution: _pages is null!\n");
            return 0;
        }

        if (_bufferSize < _pages) {
            printf("equiDistantAlignedDistribution: clusterBufferSize is too small -> raised to %d\n", _pages);
            _bufferSize = _pages;
        }

        unsigned lastPageSize   = _size % _pageSize;
        unsigned pages          = _pages;

        if (_bufferSize > _size)
            _bufferSize = _size;

        unsigned clusterSize = _bufferSize / pages;
        unsigned aclusterSize = (clusterSize / aligning) * aligning;
        if (clusterSize - aclusterSize > aligning / 2)
            aclusterSize += aligning;

		if (aclusterSize != 0) {

			if (lastPageSize > 0 && aclusterSize > lastPageSize) {
				// last page bucket would get more memory than page would need
				// --> exclude from equi-size distribution
				--pages;
				allocPage(_clusterBuffer, aclusterSize * pages + lastPageSize, me);
			} else
				allocPage(_clusterBuffer, aclusterSize * pages, me);

			PageBucketExtended<TValue> pb;
			pb.begin = _clusterBuffer.begin;

			if (pages) {
				for(unsigned i = 0; i < pages; ++i) {
					pb.end = pb.begin + aclusterSize;
					pb.cur = pb.begin;
					pb.pageOfs = 0;
					_Func(pb);
					pb.begin = pb.end;
				}
			}

			if (pages < _pages) {
				pb.end = pb.begin + lastPageSize;
				pb.cur = pb.begin;
				pb.pageOfs = 0;
				_Func(pb);
			}
		}

        return aclusterSize;
    }

}

#endif
