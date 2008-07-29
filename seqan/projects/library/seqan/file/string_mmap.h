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
  $Id: string_external.h 2416 2008-06-17 15:30:38Z doering@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

#ifndef SEQAN_HEADER_STRING_MMAP_H
#define SEQAN_HEADER_STRING_MMAP_H


//////////////////////////////////////////////////////////////////////////////

namespace SEQAN_NAMESPACE_MAIN
{

    template < typename TConfig = ExternalConfig<> >
    struct MMap;
	
	
	//////////////////////////////////////////////////////////////////////////////
    // Memory Mapped String
    //////////////////////////////////////////////////////////////////////////////

    template < typename TValue,
               typename TConfig >
	class String<TValue, MMap<TConfig> >
	{
	public:

        typedef typename TConfig::TFile							TFile;
        typedef typename TConfig::TSize							TSize;

		TValue				*data_begin;
		TValue				*data_end;
		TSize				data_capacity;

		TFile				file;
		int					_openMode;
        bool                _temporary, _ownFile;

#ifdef PLATFORM_WINDOWS
        HANDLE				handle;
#endif

		String(TSize size = 0):
			data_begin(0),
			data_end(0),
			data_capacity(0),
            file(NULL)
        {
            _temporary = true;
            _ownFile = false;

			resize(*this, size);
        }

		String(TFile &_file):
			data_begin(0),
			data_end(0),
			data_capacity(0),
            file(NULL)
        {
			open(*this, _file);
        }

		String(const char *fileName, int openMode = DefaultOpenMode<TFile>::VALUE):
			data_begin(0),
			data_end(0),
			data_capacity(0),
            file(NULL)
        {
			open(*this, fileName, openMode);
        }

		template <typename TSource>
		String & operator =(TSource const & source)
		{
	SEQAN_CHECKPOINT
			assign(*this, source);
			return *this;
		}
		String & operator =(String const & source)
		{
	SEQAN_CHECKPOINT
			assign(*this, source);
			return *this;
		}

		~String() 
		{
			close(*this);
		}

//____________________________________________________________________________

		template <typename TPos>
		inline typename Reference<String>::Type
		operator [] (TPos pos)
		{
	SEQAN_CHECKPOINT
			return value(*this, pos);
		}

		template <typename TPos>
		inline typename Reference<String const>::Type 
		operator [] (TPos pos) const
		{
	SEQAN_CHECKPOINT
			return value(*this, pos);
		}

//____________________________________________________________________________

		friend inline typename Iterator<String, Standard>::Type
		begin(String & me,
			Standard)
		{
	SEQAN_CHECKPOINT
			return me.data_begin;
		}
		friend inline typename Iterator<String const, Standard>::Type
		begin(String const & me,
			Standard)
		{
	SEQAN_CHECKPOINT
			return me.data_begin;
		}

//____________________________________________________________________________

		friend inline typename Iterator<String, Standard>::Type
		end(String & me,
			Standard)
		{
	SEQAN_CHECKPOINT
			return me.data_end;
		}
		friend inline typename Iterator<String const, Standard>::Type
		end(String const & me,
			Standard)
		{
	SEQAN_CHECKPOINT
			return me.data_end;
		}

//____________________________________________________________________________

		friend inline typename Size<String>::Type
		capacity(String & me) 
		{
	SEQAN_CHECKPOINT
			return me.data_capacity;
		}

		friend inline typename Size<String>::Type
		capacity(String const & me) 
		{
	SEQAN_CHECKPOINT
			return me.data_capacity;
		}

        inline operator bool() 
		{
            return file;
        }
//____________________________________________________________________________

		friend inline void 
		_setLength(
			String & me, 
			size_t new_length)
		{
	SEQAN_CHECKPOINT
			me.data_end = me.data_begin + new_length;
		}

//____________________________________________________________________________

		friend inline void 
		_setCapacity(
			String & me, 
			size_t new_capacity)
		{
	SEQAN_CHECKPOINT
			me.data_capacity = new_capacity;
		}

//____________________________________________________________________________

};



    //////////////////////////////////////////////////////////////////////////////
    // meta-function interface

    template < typename TValue, typename TConfig >
    struct Size< String<TValue, MMap<TConfig> > >
    {
        typedef size_t Type;
    };

    template < typename TValue, typename TConfig >
    struct Difference< String<TValue, MMap<TConfig> > >
    {
		typedef typename _MakeSigned<size_t>::Type Type;
    };
//____________________________________________________________________________

    template < typename TValue, typename TConfig >
	struct DefaultOverflowExplicit<String<TValue, MMap<TConfig> > >
	{
		typedef Generous Type;
	};

    template < typename TValue, typename TConfig >
	struct DefaultOverflowImplicit<String<TValue, MMap<TConfig> > >
	{
		typedef Generous Type;
	};
//____________________________________________________________________________

    template < typename TValue, typename TConfig >
	struct IsContiguous< String<TValue, MMap<TConfig> > >
	{
		typedef True Type;
		enum { VALUE = true };
	};

    template < typename TValue, typename TConfig >
	struct AllowsFastRandomAccess< String<TValue, MMap<TConfig> > >
	{
		typedef False Type;
		enum { VALUE = false };
	};


	//////////////////////////////////////////////////////////////////////////////
    // global interface

//____________________________________________________________________________

    template < typename TValue, typename TConfig >
	inline void 
	waitForAll(String<TValue, MMap<TConfig> > &me)
	{
	}
	
    template < typename TValue, typename TConfig >
    inline void 
	flush(String<TValue, MMap<TConfig> > &me) 
	{
#ifndef PLATFORM_WINDOWS
		::msync(me.data_begin, length(me) * sizeof(TValue), MS_SYNC);
#endif
    }

	// cancel all transactions
    template < typename TValue, typename TConfig >
	inline void 
	cancel(String<TValue, MMap<TConfig> > &me)
	{
#ifndef PLATFORM_WINDOWS
		::msync(me.data_begin, capacity(me) * sizeof(TValue), MS_INVALIDATE);
#endif
	}

	// flush and free all allocated pages
    template < typename TValue, typename TConfig >
	inline void 
	free(String<TValue, MMap<TConfig> > &me)
	{
#ifndef PLATFORM_WINDOWS
		::madvise(me.data_begin, capacity(me) * sizeof(TValue), MADV_DONTNEED);
#endif
	}
//____________________________________________________________________________

#ifdef PLATFORM_WINDOWS

    static SECURITY_ATTRIBUTES MMapStringDefaultAttributes = {
        sizeof(SECURITY_ATTRIBUTES),
        NULL,
        true
    };

	template < typename TValue, typename TConfig, typename TSize >
    inline bool 
    _map(String<TValue, MMap<TConfig> > &me, TSize new_capacity) 
	{
		if (new_capacity > 0) 
		{
			resize(me.file, new_capacity * sizeof(TValue));
			DWORD prot = 0;
			DWORD access = 0;
			if ((me._openMode & OPEN_MASK) == OPEN_RDONLY) 
			{
				prot = PAGE_READONLY;
				access = FILE_MAP_READ;
			} else {
				prot = PAGE_READWRITE;
				access = FILE_MAP_ALL_ACCESS;
			}
            LARGE_INTEGER size;
			size.QuadPart = new_capacity;
			size.QuadPart *= sizeof(TValue);

			me.handle = CreateFileMapping(me.file.handle, &MMapStringDefaultAttributes, prot, size.HighPart, size.LowPart, NULL);
			if (me.handle == NULL)
			{
			#ifdef SEQAN_DEBUG
				::std::cerr << "CreateFileMapping failed. (ErrNo=" << GetLastError() << ")" << ::std::endl;
			#endif
				return false;
			}

			void *addr = MapViewOfFile(me.handle, access, 0, 0, 0);	
			if (addr == NULL)
			{
			#ifdef SEQAN_DEBUG
				::std::cerr << "MapViewOfFile failed. (ErrNo=" << GetLastError() << ")" << ::std::endl;
			#endif
				return false;
			}
				
			me.data_begin = (TValue *) addr;
			_setLength(me, new_capacity);
			_setCapacity(me, new_capacity);
		}
		return true;
	}

	template < typename TValue, typename TConfig, typename TCapSize >
    inline bool 
    _remap(String<TValue, MMap<TConfig> > &me, TCapSize new_capacity) 
	{
		typedef typename Size< String<TValue, MMap<TConfig> > >::Type TSize;

		if (me.data_begin) 
		{
			TSize seq_length = length(me);
			
			DWORD prot = 0;
			DWORD access = 0;
			if ((me._openMode & OPEN_MASK) == OPEN_RDONLY) 
			{
				prot = PAGE_READONLY;
				access = FILE_MAP_READ;
			} else {
				prot = PAGE_READWRITE;
				access = FILE_MAP_ALL_ACCESS;
			}
            LARGE_INTEGER size;
			size.QuadPart = new_capacity;
			size.QuadPart *= sizeof(TValue);

			bool result = true;
			result &= (UnmapViewOfFile(me.data_begin) != 0);
			result &= (CloseHandle(me.handle) != 0);

			HANDLE handle = CreateFileMapping(me.file.handle, &MMapStringDefaultAttributes, prot, size.HighPart, size.LowPart, NULL);
			if (handle == NULL)
			{
			#ifdef SEQAN_DEBUG
				::std::cerr << "CreateFileMapping failed. (ErrNo=" << GetLastError() << ")" << ::std::endl;
			#endif
				return false;
			}

			void *addr = MapViewOfFile(handle, access, 0, 0, 0);	
			if (addr == NULL)
			{
			#ifdef SEQAN_DEBUG
				::std::cerr << "MapViewOfFile failed. (ErrNo=" << GetLastError() << ")" << ::std::endl;
			#endif
				return false;
			}

			if (capacity(me) > new_capacity)
				resize(me.file, new_capacity * sizeof(TValue));

			me.handle = handle;
			me.data_begin = (TValue*) addr;
			_setLength(me, seq_length);
			_setCapacity(me, new_capacity);
			return true;
		} else
			return _map(me, new_capacity);
	}

	template < typename TValue, typename TConfig >
    inline bool 
    _unmap(String<TValue, MMap<TConfig> > &me) 
	{
		bool result = true;
		if (me.data_begin) 
		{
			if (UnmapViewOfFile(me.data_begin) == 0)
			{
			#ifdef SEQAN_DEBUG
				::std::cerr << "UnmapViewOfFile failed. (ErrNo=" << GetLastError() << ")" << ::std::endl;
			#endif
				result = false;
			}
			
			if (CloseHandle(me.handle) == 0)
			{
			#ifdef SEQAN_DEBUG
				::std::cerr << "CloseHandle failed. (ErrNo=" << GetLastError() << ")" << ::std::endl;
			#endif
				result = false;
			}
			
			resize(me.file, length(me) * sizeof(TValue));
			me.data_begin = NULL;
			me.data_end = NULL;
			me.data_capacity = 0;
		}
		return result;
	}

#else

	template < typename TValue, typename TConfig, typename TSize >
    inline bool 
    _map(String<TValue, MMap<TConfig> > &me, TSize new_capacity) 
	{
		if (new_capacity > 0) 
		{
			resize(me.file, new_capacity * sizeof(TValue));
			int prot = 0;
			if (me._openMode & OPEN_RDONLY) prot |= PROT_READ;
			if (me._openMode & OPEN_WRONLY) prot |= PROT_WRITE;
			void *addr = ::mmap(NULL, new_capacity * sizeof(TValue), prot, MAP_SHARED, me.file.handle, 0);
			
			if (addr == MAP_FAILED)
			{
			#ifdef SEQAN_DEBUG
				::std::cerr << "mmap failed. errno=" << errno << " (" << ::strerror(errno) << ")" << ::std::endl;
			#endif
				return false;
			}
				
			me.data_begin = (TValue *) addr;
			_setLength(me, new_capacity);
			_setCapacity(me, new_capacity);
		}
		return true;
	}

	template < typename TValue, typename TConfig, typename TCapSize >
    inline bool 
    _remap(String<TValue, MMap<TConfig> > &me, TCapSize new_capacity) 
	{
		typedef typename Size< String<TValue, MMap<TConfig> > >::Type TSize;

		if (me.data_begin) 
		{
			TSize seq_length = length(me);
			
			if (capacity(me) < new_capacity)
				resize(me.file, new_capacity * sizeof(TValue));

#ifdef MREMAP_MAYMOVE
			void *addr = ::mremap(me.data_begin, capacity(me) * sizeof(TValue), new_capacity * sizeof(TValue), MREMAP_MAYMOVE);
#else
			// for BSD systems without mremap(..) like Mac OS X ...
			int prot = 0;
			if (me._openMode & OPEN_RDONLY) prot |= PROT_READ;
			if (me._openMode & OPEN_WRONLY) prot |= PROT_WRITE;
//			void *addr = ::mmap(me.data_begin, new_capacity * sizeof(TValue), prot, MAP_SHARED | MAP_FIXED, me.file.handle, 0);
			::munmap(me.data_begin, capacity(me) * sizeof(TValue));
			void *addr = ::mmap(NULL, new_capacity * sizeof(TValue), prot, MAP_SHARED, me.file.handle, 0);
#endif

			if (addr == MAP_FAILED) 
			{
			#ifdef SEQAN_DEBUG
				::std::cerr << "mremap failed. errno=" << errno << " (" << ::strerror(errno) << ")" << ::std::endl;
			#endif
				return false;
			}

			if (capacity(me) > new_capacity)
				resize(me.file, new_capacity * sizeof(TValue));

			me.data_begin = (TValue*) addr;
			_setLength(me, seq_length);
			_setCapacity(me, new_capacity);
			return true;
		} else
			return _map(me, new_capacity);
	}

	template < typename TValue, typename TConfig >
    inline bool 
    _unmap(String<TValue, MMap<TConfig> > &me) 
	{
		if (me.data_begin) 
		{
			int error = ::munmap(me.data_begin, capacity(me) * sizeof(TValue));
			if (error != 0)
			{
			#ifdef SEQAN_DEBUG
				::std::cerr << "munmap failed. errno=" << errno << " (" << ::strerror(errno) << ")" << ::std::endl;
			#endif
				return false;
			}
			
			resize(me.file, length(me) * sizeof(TValue));
			me.data_begin = NULL;
			me.data_end = NULL;
			me.data_capacity = 0;
		}
		return true;
	}

#endif

	template < typename TValue, typename TConfig, typename TSize, typename TExpand >
    inline bool 
    _remap(String<TValue, MMap<TConfig> > &me, TSize new_capacity, Tag<TExpand> const expand) 
	{
		return _remap(me, new_capacity);
	}

	template < typename TValue, typename TConfig, typename TSize >
    inline bool 
    _remap(String<TValue, MMap<TConfig> > &me, TSize new_capacity, Generous) 
	{
		return _remap(me, computeGenerousCapacity(me, new_capacity));
	}

	template < typename TValue, typename TConfig >
    inline void 
    clear(String<TValue, MMap<TConfig> > &me) 
	{
		cancel(me);
		_unmap(me);
		resize(me.file, 0);
	}
//____________________________________________________________________________

	template < typename TValue, typename TConfig >
    inline bool 
    open(String<TValue, MMap<TConfig> > &me, const char *fileName, int openMode) 
	{
		close(me);
		me._temporary = false;
				
		if ((me._ownFile = open(me.file, fileName, openMode))) 
		{
			me._openMode = openMode;
			return _map(me, size(me.file) / sizeof(TValue));
		}

		return false;
    }

	template < typename TValue, typename TConfig >
    inline bool 
    open(String<TValue, MMap<TConfig> > &me, const char *fileName) 
	{
		typedef typename String<TValue, MMap<TConfig> >::Type TFile;
		return open(me, fileName, DefaultOpenMode<TFile>::VALUE);
    }

	template < typename TValue, typename TConfig >
    inline bool 
    open(String<TValue, MMap<TConfig> > &me, typename TConfig::TFile file) 
	{
		close(me);
		me.file = file;
        me._temporary = false;
        me._ownFile = false;

		if (me.file) 
		{
			me._openMode = OPEN_RDWR;
			return _map(me, size(me.file) / sizeof(TValue));
		}

		return false;
    }

	template < typename TValue, typename TConfig >
    inline bool 
    openTemp(String<TValue, MMap<TConfig> > &me) 
	{
		close(me);
        me._temporary = true;
		me._openMode = OPEN_RDWR;
		
		return me._ownFile = openTemp(me.file);
    }
//____________________________________________________________________________

	template < typename TValue, typename TConfig >
	inline void _ensureFileIsOpen(String<TValue, MMap<TConfig> > &me) 
	{
		if (!me.file)
		{
			me._temporary = true;
			if (!(me._ownFile = openTemp(me.file)))
				::std::cerr << "Memory Mapped String couldn't open temporary file" << ::std::endl;
		}
	}
//____________________________________________________________________________

    template < typename TValue, typename TConfig, typename TNewSize, typename TExpand >
    inline typename Size< String<TValue, MMap<TConfig> > >::Type
    reserve(
	    String<TValue, MMap<TConfig> > &me,
		TNewSize new_capacity,
		Tag<TExpand> const expand)
	{
		typedef typename Size< String<TValue, MMap<TConfig> > >::Type TSize;

		TSize old_capacity = capacity(me);
		if (old_capacity >= (TSize)new_capacity) return new_capacity;
			
		if (new_capacity > 0)
		{
			_ensureFileIsOpen(me);
			_remap(me, new_capacity, expand);
		} else
			_unmap(me);
		
		return capacity(me);
	}
//____________________________________________________________________________

	template < typename TValue, typename TConfig >
    inline bool 
    save(String<TValue, MMap<TConfig> > const &me, const char *fileName, int openMode) {
		// Memory Mapped Strings are persistent, thus there is no need to save them
		//MMapStringsDontNeedToBeSaved error;
		return true;
	}

	template < typename TValue, typename TConfig >
    inline bool 
    save(String<TValue, MMap<TConfig> > const &me, const char *fileName) {
		// Memory Mapped Strings are persistent, thus there is no need to save them
		//MMapStringsDontNeedToBeSaved error;
		return true;
	}

	template < typename TValue, typename TConfig >
    inline bool 
    save(String<TValue, MMap<TConfig> > const &me, typename TConfig::TFile file) {
		// Memory Mapped Strings are persistent, thus there is no need to save them
		//MMapStringsDontNeedToBeSaved error;
		return true;
	}
//____________________________________________________________________________

	template < typename TValue, typename TConfig >
    inline bool 
    close(String<TValue, MMap<TConfig> > &me) 
	{
		if (me.file) 
		{
			// close associated file
			if (me._temporary) 
			{
				me._temporary = false;
				cancel(me);
			}

			_unmap(me);

			if (me._ownFile) 
			{
				me._ownFile = false;
				return close(me.file);
			}
		}
		return true;
    }
//____________________________________________________________________________
// ugly by-passing the string_base interface as it was designed for sequential
// strings and their allocation/deallocation properties (no remap/reallocate)

	template < typename TValue, typename TConfig, typename TSource, typename TExpand >
	inline void
	assign(String<TValue, MMap<TConfig> > &target, 
				TSource const &source,
				Tag<TExpand> const expand)
	{
		typedef String<TValue, External<TConfig> >	TTarget;

		typename Size<TSource>::Type source_length	= length(source);
		typename Size<TTarget>::Type part_length	= reserve(target, source_length, expand);
		
		if (part_length > source_length)
			part_length = source_length;
		
		arrayConstructCopy(begin(source, Standard()), begin(source, Standard()) + part_length, begin(target, Standard()));
		_setLength(target, part_length);
	}

	template < typename TValue, typename TConfig, typename TSource, typename TExpand >
	inline void
	append(String<TValue, MMap<TConfig> > &target, 
				TSource const &source,
				Tag<TExpand> const expand)
	{
		typedef String<TValue, External<TConfig> >	TTarget;

		typename Size<TSource>::Type source_length	= length(source);
		typename Size<TTarget>::Type target_length	= length(target);
		typename Size<TTarget>::Type part_length	= reserve(target, source_length + target_length, expand) - target_length;
		
		if (part_length > source_length)
			part_length = source_length;
		
		arrayConstructCopy(begin(source, Standard()), begin(source, Standard()) + part_length, begin(target, Standard()) + target_length);
		_setLength(target, target_length + part_length);
	}


//////////////////////////////////////////////////////////////////////////////



} //namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
