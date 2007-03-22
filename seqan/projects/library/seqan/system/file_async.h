/*
 *  file_async.h
 *  SeqAn
 *
 *  Created by David Weese on 17.07.05.
 *
 */

//SEQAN_NO_GENERATED_FORWARDS: no forwards are generated for this file

#ifndef SEQAN_HEADER_FILE_ASYNC_H
#define SEQAN_HEADER_FILE_ASYNC_H

namespace SEQAN_NAMESPACE_MAIN
{

    //////////////////////////////////////////////////////////////////////////////
    // interface for typed files - deprecated
//     template < typename T = unsigned char >
//     struct AsyncConfig;
// 
//     template < typename TConfig = AsyncConfig<> >
//     struct Async;


#ifdef PLATFORM_WINDOWS


    static DWORD _transferedBytes;  // for reporting

    template < typename TConfig >
	class File<Async<TConfig> >
    {

    public:

        typedef LONGLONG    FilePtr;
        typedef ULONGLONG   SizeType;
        typedef DWORD       _SizeType;
        typedef HANDLE      Handle;
//        typedef IOQueue     Queue;

//        CountedPtr<Queue>   queueHolder;
//        Queue*              queue;
        Handle              hFile, hFileAsync;
//        FilePtr             position;
//        DWORD               aligning;
        bool                noBuffering;

        File(void *dummy = NULL): // to be compatible with the FILE*(NULL) constructor
            hFile(INVALID_HANDLE_VALUE) {}

/*        File(Queue *_queue, char const *fileName = NULL, int openMode = OPEN_RDWR + OPEN_CREATE | OPEN_APPEND):
            queue(_queue),
            hFile(INVALID_HANDLE_VALUE),
            position(0)
        {
            if (fileName)
                open(fileName, openMode);
            else
                openTemp();
        }

        File(char const *fileName = NULL, int openMode = OPEN_RDWR + OPEN_CREATE | OPEN_APPEND):
//            queueHolder(new Queue),
            hFile(INVALID_HANDLE_VALUE),
            position(0)
        {
//            queue = queueHolder;
            if (fileName)
                open(fileName, openMode);
            else
                openTemp();
        }
*/

        ~File() {
//            close();
        }

        bool open(char const *fileName, int openMode = OPEN_RDWR + OPEN_CREATE | OPEN_APPEND) {
			SEQAN_PROADD(PROOPENFILES, 1);
            noBuffering = getExtraFlags(openMode | OPEN_ASYNC) & (FILE_FLAG_NO_BUFFERING | FILE_FLAG_OVERLAPPED);
            hFileAsync = CreateFile(fileName,
                                getFileAccess(openMode | OPEN_ASYNC),
                                FILE_SHARE_DELETE | FILE_SHARE_READ | FILE_SHARE_WRITE,
                                NULL,
                                getCreationFlags(openMode | OPEN_ASYNC),
                                getExtraFlags(openMode | OPEN_ASYNC),
                                NULL);

            if (hFileAsync == INVALID_HANDLE_VALUE) {
				::std::cerr << "Open failed on file " << fileName << ". (ErrNo=" << GetLastError() << ")" << ::std::endl;
                return false;
            }
            #ifdef SEQAN_VERBOSE
                ::std::cout << "file opened asynchronously " << fileName << " handle " << ::std::hex << hFileAsync << ::std::dec << ::std::endl;
            #endif

            if (noBuffering) {
                hFile = CreateFile(fileName,                // in this case io must be sector aligned
                                getFileAccess(openMode),    // so we open a second file, for unaligned access
                                FILE_SHARE_DELETE | FILE_SHARE_READ | FILE_SHARE_WRITE,
                                NULL,
                                OPEN_EXISTING,
                                getExtraFlags(openMode & ~OPEN_ASYNC),
                                NULL);
                if (hFile == INVALID_HANDLE_VALUE) {
                	::std::cerr << "Open failed on secondary file " << fileName << ". (ErrNo=" << GetLastError() << ")" << ::std::endl;
                    return false;
                }
	            #ifdef SEQAN_VERBOSE
                	::std::cout << "async file opened  " << fileName << " handle " << ::std::hex << hFile << ::std::dec << ::std::endl;
                #endif
            } else
                hFile = hFileAsync;

            return true;
        }

        bool openTemp(int openMode = OPEN_RDWR + OPEN_CREATE) {
            char szTempName[MAX_PATH];
#ifdef SEQAN_DEFAULT_TMPDIR
            char szTempPath[MAX_PATH] = SEQAN_DEFAULT_TMPDIR;
#else
            char szTempPath[MAX_PATH];
            if (!GetTempPath(MAX_PATH, szTempPath)) {
				::std::cerr << "Couldn't get a temporary path name. (ErrNo=" << GetLastError() << ")" << ::std::endl;
                return false;
            }
#endif
            if (!GetTempFileName(szTempPath, "GNDX", 0, szTempName)) {
				::std::cerr << "Couldn't get a temporary file name. (ErrNo=" << GetLastError() << ")" << ::std::endl;
                return false;
            }
            return open(szTempName, openMode | OPEN_TEMPORARY);
        }

        inline bool close() {
            BOOL result = true;
            #ifdef SEQAN_VERBOSE
                ::std::cout << "files closed handles " << ::std::hex << hFileAsync << " and " << hFile << ::std::dec << ::std::endl;
            #endif
            if (hFile != hFileAsync)
                result &= CloseHandle(hFileAsync);
            result &= CloseHandle(hFile);
            hFileAsync = INVALID_HANDLE_VALUE;
            hFile = INVALID_HANDLE_VALUE;
			SEQAN_PROSUB(PROOPENFILES, 1);
            return result;
        }

        inline bool read(void *memPtr, _SizeType count) const {
            SEQAN_PROADD(PROIO, (count + PROPAGESIZE - 1) / PROPAGESIZE);
            SEQAN_PROTIMESTART(tw);
		    bool result = ReadFile(hFile, memPtr, count, &_transferedBytes, NULL);
            SEQAN_PROADD(PROCWAIT, SEQAN_PROTIMEDIFF(tw));
            return result;
        }

        inline bool write(void const *memPtr, _SizeType count) const {
            SEQAN_PROADD(PROIO, (count + PROPAGESIZE - 1) / PROPAGESIZE);
            SEQAN_PROTIMESTART(tw);
		    bool result = WriteFile(hFile, memPtr, count, &_transferedBytes, NULL);
            SEQAN_PROADD(PROCWAIT, SEQAN_PROTIMEDIFF(tw));
            return result;
        }

		inline FilePtr seek(FilePtr _pos, DWORD origin = FILE_BEGIN) {
//          LARGE_INTEGER li = _pos;
//			return SetFilePointer(hFileAsync, li.LowPart, &li.HighPart, MoveMethod);
            LARGE_INTEGER new_pos, pos;
            pos.QuadPart = _pos;
            SetFilePointerEx(hFile, pos, &new_pos, origin);
//            position = new_pos.QuadPart;
            return new_pos.QuadPart;
		}

		inline FilePtr tell() {
			return seek(0, FILE_CURRENT);
        }

		inline FilePtr size() const {
            LARGE_INTEGER result;
            DWORD dwError, high;
            result.LowPart = GetFileSize(hFile, &high);
            result.HighPart = high;
            if (result.LowPart == INVALID_FILE_SIZE && (dwError = GetLastError()) != NO_ERROR) {
				::std::cerr << "Couldn't get file size. (ErrNo=" << dwError << ")" << ::std::endl;
                return 0;
            }
            return result.QuadPart;
        }

        inline bool setEOF() const {
            return SetEndOfFile(hFile);
        }

		inline static DWORD error() {
			return GetLastError();
		}

        operator bool () const {
            return (hFile != INVALID_HANDLE_VALUE) && (hFileAsync != INVALID_HANDLE_VALUE);
        }

    protected:

        DWORD getFileAccess(int openMode) {
            DWORD result;
            switch (openMode & OPEN_MASK) {
                case OPEN_RDONLY:
                    result = GENERIC_READ;
                    break;
                case OPEN_WRONLY:
                    result = GENERIC_WRITE;
                    break;
                case OPEN_RDWR:
                    result = GENERIC_READ | GENERIC_WRITE;
            }
            return result;
        }

        DWORD getCreationFlags(int openMode) {
            if (openMode & OPEN_CREATE)
                if (openMode & OPEN_APPEND)
                    return OPEN_ALWAYS;
                else
                    return CREATE_ALWAYS;
            else
                return OPEN_EXISTING;
        }

        DWORD getExtraFlags(int openMode) {
            DWORD extra = FILE_ATTRIBUTE_NORMAL | FILE_FLAG_RANDOM_ACCESS;// | FILE_FLAG_WRITE_THROUGH;
            if (openMode & OPEN_ASYNC) {
                extra |= FILE_FLAG_OVERLAPPED;
                #ifdef SEQAN_DIRECTIO
                    extra |= FILE_FLAG_NO_BUFFERING;
                #endif
            }
            if (openMode & OPEN_TEMPORARY)  extra |= FILE_FLAG_DELETE_ON_CLOSE;
            return extra;
        }

    };


    //////////////////////////////////////////////////////////////////////////////
    // (SeqAn adaption)
    //////////////////////////////////////////////////////////////////////////////

    struct aiocb_win32 {
        OVERLAPPED  overlapped;
        Event       xmitDone;
    };

    template < typename TConfig >
    struct aRequest<File<Async<TConfig> > >
    {
        typedef aiocb_win32 Type;
    };
/*
    template < typename TConfig >
    struct aEvent<File<Async<TConfig> > >
    {
        typedef Event Type;
    };


    template < typename TConfig >
    struct aQueue<File<Async<TConfig> > >
    {
        typedef IOQueue Type;
    };

    template < typename TConfig >
    struct aHint<File<Async<TConfig> > >
    {
        typedef typename aQueue<File<Async<TConfig> > >::Type::aHint Type;
    };

    template < typename TConfig >
    struct aCallback<File<Async<TConfig> > >
    {
        typedef typename aQueue<File<Async<TConfig> > >::Type::aCallback Type;
    };*/


    template < typename TConfig >
    inline typename Size<File<Async<TConfig> > >::Type size(File<Async<TConfig> > &me) {
        return me.size();
    }

    template < typename TConfig >
    inline bool setEOF(File<Async<TConfig> > &me) {
        return me.setEOF();
    }

    template < typename TConfig >
    inline void reopen(File<Async<TConfig> > & me, int openMode) {
        me.reopen(openMode);
    }

    template < typename TConfig >
    inline unsigned sectorSize(File<Async<TConfig> > const &me) {
        DWORD SpC, nofC, tnoC, aligning;
        if (GetDiskFreeSpace(NULL, &SpC, &aligning, &nofC, &tnoC) == 0)  {
            ::std::cerr << "Error " << GetLastError() << " while querying cluster size" << ::std::endl;
            return 4096;
        }
        return aligning;
    }


    template < typename TConfig, typename TValue, typename TSize, typename TPos >
    inline bool areadAt(File<Async<TConfig> > & me, TValue *memPtr, TSize const count, TPos const fileOfs,
        aiocb_win32 &request)
    {
        SEQAN_PROTIMESTART(tw);
        LARGE_INTEGER ofs;
        ofs.QuadPart = fileOfs;
        ofs.QuadPart *= sizeof(TValue);
        request.overlapped.Offset = ofs.LowPart;
        request.overlapped.OffsetHigh = ofs.HighPart;
        if (!request.xmitDone) open(request.xmitDone);
        request.overlapped.hEvent = request.xmitDone.hEvent;
        if (ReadFile(
            me.hFileAsync, 
            memPtr, 
            count * sizeof(TValue),
            &ofs.LowPart,
            &request.overlapped) || (me.error() == ERROR_IO_PENDING))
        {
            SEQAN_PROADD(PROIO, (sizeof(TValue) * count + PROPAGESIZE - 1) / PROPAGESIZE);
            SEQAN_PROADD(PROIWAIT, SEQAN_PROTIMEDIFF(tw));
            return true;
        }
        if (me.error() == ERROR_NO_SYSTEM_RESOURCES) {  // read synchronoulsy instead
            #ifdef SEQAN_DEBUG_OR_TEST_
            	::std::cerr << "Warning: Falling back to sync. read. :( " << ::std::endl;
            #endif
			signal(request.xmitDone);
            return readAt(me, memPtr, count, fileOfs);
        }
        return false;
    }
    
    template < typename TConfig, typename TValue, typename TSize, typename TPos >
    inline bool awriteAt(File<Async<TConfig> > & me, TValue const *memPtr, TSize const count, TPos const fileOfs,
        aiocb_win32 &request)
    {
        SEQAN_PROTIMESTART(tw);
        LARGE_INTEGER ofs;
        ofs.QuadPart = fileOfs;
        ofs.QuadPart *= sizeof(TValue);
        request.overlapped.Offset = ofs.LowPart;
        request.overlapped.OffsetHigh = ofs.HighPart;
        if (!request.xmitDone) open(request.xmitDone);
        request.overlapped.hEvent = request.xmitDone.hEvent;
        if (WriteFile(
            me.hFileAsync, 
            memPtr, 
            count * sizeof(TValue),
            &ofs.LowPart,
            &request.overlapped) || (me.error() == ERROR_IO_PENDING))
        {
            SEQAN_PROADD(PROIO, (sizeof(TValue) * count + PROPAGESIZE - 1) / PROPAGESIZE);
            SEQAN_PROADD(PROIWAIT, SEQAN_PROTIMEDIFF(tw));
            return true;
        }
        if (me.error() == ERROR_NO_SYSTEM_RESOURCES) {  // write synchronoulsy instead
            #ifdef SEQAN_DEBUG_OR_TEST_
            	::std::cerr << "Warning: Falling back to sync. write. :( " << ::std::endl;
            #endif
			signal(request.xmitDone);
            return writeAt(me, memPtr, count, fileOfs);
        }
        return false;
    }

    //////////////////////////////////////////////////////////////////////
    // queue specific functions

    inline bool waitFor(aiocb_win32 &request) {
        SEQAN_PROTIMESTART(tw);
		if (!waitFor(request.xmitDone, 60000))
            ::std::cout << "waitFor timeout\n";
        SEQAN_PROADD(PROCWAIT, SEQAN_PROTIMEDIFF(tw));
        return true;
	}

    template < typename TTime >
    inline bool waitFor(aiocb_win32 &request, TTime timeout_millis) {
        SEQAN_PROTIMESTART(tw);
		bool result = waitFor(request.xmitDone, timeout_millis);
        SEQAN_PROADD(PROCWAIT, SEQAN_PROTIMEDIFF(tw));
        return result;
	}

	template < typename TSize >
	inline TSize waitForAny(aiocb_win32 const * const contexts[], TSize count, DWORD timeout_millis = Event::Infinite) {
        Event::Handle handles[count];
        for(TSize i = 0; i < count; ++i)
            handles[i] = contexts[i]->xmitDone.hEvent;

        SEQAN_PROTIMESTART(tw);
        DWORD result = WaitForMultipleObjects(count, &handles, false, timeout_millis);
        SEQAN_PROADD(PROCWAIT, SEQAN_PROTIMEDIFF(tw));
        if (result >= WAIT_OBJECT_0 && result < WAIT_OBJECT_0 + count)
    		return result - WAIT_OBJECT_0;
        return count;
	}

    template < typename TConfig >
    inline bool cancel(File<Async<TConfig> > & me, aiocb_win32 const &request) {
        return CancelIo(me.handleAsync);
    }

    template < typename TConfig >
    inline bool flush(File<Async<TConfig> > & me) {
		if (me.hFile != me.hFileAsync)	// in case of equality no direct access was done -> no flush needed
        	return FlushFileBuffers(me.hFile);
        else
            return true;
    }

    template < typename TConfig, typename aRequest >
    inline void release(File<Async<TConfig> > & me, aRequest & request) { }


/*        
    //////////////////////////////////////////////////////////////////////
    // callback based read/write

    template < typename TConfig, typename TValue, typename TSize,
               typename aCallback, typename aHint >
    inline typename aRequest<File<Async<TConfig> > >::Type
    aread(File<Async<TConfig> > & me, TValue *memPtr, TSize const count,
        aCallback* cb, aHint* hint)
    {
        DWORD bsize = (DWORD)(count * sizeof(TValue));
        typename aRequest<File<Async<TConfig> > >::Type request = 
            me.queue->areadAt(
                me.hFileAsync,
                me.position,
                memPtr,
                bsize,
                cb,
                hint);
        me.position += bsize;
        return request;
    }
    
    template < typename TConfig, typename TValue, typename TSize,
               typename aCallback, typename aHint >
    inline typename aRequest<File<Async<TConfig> > >::Type
    awrite(File<Async<TConfig> > & me, TValue const *memPtr, TSize const count,
        aCallback* cb, aHint* hint)
    {
        DWORD bsize = (DWORD)(count * sizeof(TValue));
        typename aRequest<File<Async<TConfig> > >::Type request = 
            me.queue->awriteAt(
                memPtr,
                me.hFileAsync,
                me.position,
                bsize,
                cb,
                hint);
        me.position += bsize;
        return request;
    }

    template < typename TConfig, typename TValue, typename TSize, typename TPos,
               typename aCallback, typename aHint >
    inline typename aRequest<File<Async<TConfig> > >::Type
    areadAt(File<Async<TConfig> > & me, TValue *memPtr, TSize const count, TPos const fileOfs,
        aCallback* cb, aHint* hint)
    {
        DWORD bsize = (DWORD)(count * sizeof(TValue));
        return me.queue->areadAt(
            me.hFileAsync,
            fileOfs * sizeof(TValue),
            memPtr,
            bsize,
            cb,
            hint);
    }
    
    template < typename TConfig, typename TValue, typename TSize, typename TPos,
               typename aCallback, typename aHint >
    inline typename aRequest<File<Async<TConfig> > >::Type
    awriteAt(File<Async<TConfig> > & me, TValue *memPtr, TSize const count, TPos const fileOfs,
        aCallback* cb, aHint* hint)
    {
        DWORD bsize = (DWORD)(count * sizeof(TValue));
        return me.queue->awriteAt(
            memPtr,
            me.hFileAsync,
            fileOfs * sizeof(TValue),
            bsize,
            cb,
            hint);
    }


    //////////////////////////////////////////////////////////////////////
    // event based read/write

    template < typename TConfig, typename TValue, typename TSize,
               typename aEvent >
    inline typename aRequest<File<Async<TConfig> > >::Type
    aread(File<Async<TConfig> > & me, TValue *memPtr, TSize const count,
        aEvent &event)
    {
        DWORD bsize = (DWORD)(count * sizeof(TValue));
        typename aRequest<File<Async<TConfig> > >::Type request = 
            me.queue->areadAt(
                me.hFileAsync,
                me.position,
                memPtr,
                bsize,
                event);
        me.position += bsize;
        return request;
    }
    
    template < typename TConfig, typename TValue, typename TSize,
               typename aEvent >
    inline typename aRequest<File<Async<TConfig> > >::Type
    awrite(File<Async<TConfig> > & me, TValue *memPtr, TSize const count,
        aEvent &event)
    {
        DWORD bsize = (DWORD)(count * sizeof(TValue));
        typename aRequest<File<Async<TConfig> > >::Type request =  
            me.queue->awriteAt(
                memPtr,
                me.hFileAsync,
                me.position,
                bsize,
                event);
        me.position += bsize;
        return request;
    }

    template < typename TConfig, typename TValue, typename TSize, typename TPos,
               typename aEvent >
    inline typename aRequest<File<Async<TConfig> > >::Type
    areadAt(File<Async<TConfig> > & me, TValue *memPtr, TSize const count, TPos const fileOfs,
        aEvent &event)
    {
        DWORD bsize = (DWORD)(count * sizeof(TValue));
        return me.queue->areadAt(
            me.hFileAsync,
            fileOfs * sizeof(TValue),
            memPtr,
            bsize,
            event);
    }
    
    template < typename TConfig, typename TValue, typename TSize, typename TPos,
               typename aEvent >
    inline typename aRequest<File<Async<TConfig> > >::Type
    awriteAt(File<Async<TConfig> > & me, TValue *memPtr, TSize const count, TPos const fileOfs,
        aEvent &event)
    {
        DWORD bsize = (DWORD)(count * sizeof(TValue));
        return me.queue->awriteAt(
            memPtr,
            me.hFileAsync,
            fileOfs * sizeof(TValue),
            bsize,
            event);
    }


    //////////////////////////////////////////////////////////////////////
    // queue specific functions

    template < typename TConfig >
    inline void flush(File<Async<TConfig> > & me) {
        me.queue->flush();
    }

    template < typename TConfig, typename aRequest >
    inline void release(File<Async<TConfig> > & me, aRequest & request) {
        me.queue->release(request);
    }
*/

	//////////////////////////////////////////////////////////////////////////////
	// page aligned allocate for direct file io

    struct TagAllocateAligned_;	//< allocate page aligned memory for direct i/o access
    typedef Tag<TagAllocateAligned_> const TagAllocateAligned;

	template <typename T, typename TValue, typename TSize>
	inline void
	allocate(T const & me, 
			 TValue * & data,
			 TSize count,
			 TagAllocateAligned const)
	{
		data = (TValue *) VirtualAlloc(NULL, count * sizeof(TValue), MEM_COMMIT, PAGE_READWRITE);
        if (!data)
			::std::cerr << "AlignAllocator: Could not allocate memory of size " << ::std::hex << count * sizeof(TValue) << ::std::dec << ". (ErrNo=" << GetLastError() << ")" << ::std::endl;
        else
            SEQAN_PROADD(PROMEMORY, sizeof(TValue) * count);
	}

	//////////////////////////////////////////////////////////////////////////////
	// page aligned deallocate for direct file io

	template <typename T, typename TValue, typename TSize>
	inline void 
	deallocate( T const & me,
				TValue * data, 
				TSize count,
				TagAllocateAligned const)
	{
        if (data) {
            SEQAN_PROSUB(PROMEMORY, sizeof(TValue) * count);
			VirtualFree(data, 0, MEM_RELEASE);
        }
	}

#else

    
    template < typename TConfig >
    class File<Async<TConfig> > : public File<Sync<TConfig> >
    {

    public:

        typedef File<Sync<TConfig> >  Base;

        typedef off_t			FilePtr;
		typedef off_t           SizeType;   // type of file size
        typedef size_t          _SizeType;  // type of transfer size (for read or write)
		typedef int				Handle;

        Handle handleAsync;
		using Base::handle;

		File(void *dummy = NULL): 	// to be compatible with the FILE*(NULL) constructor
			handleAsync(-1) {}

        bool open(char const *fileName, int openMode = OPEN_RDWR + OPEN_CREATE | OPEN_APPEND) {
            handle = ::open(fileName, Base::getOFlag(openMode & ~OPEN_ASYNC), S_IREAD | S_IWRITE);
			if (handle == -1) {
				::std::cerr << "Open failed on file " << fileName << ". (" << strerror(errno) << ")" << ::std::endl;
				return false;
			}

			if (Base::getOFlag(openMode | OPEN_ASYNC) & O_DIRECT) {
				handleAsync = ::open(fileName, Base::getOFlag(openMode | OPEN_ASYNC & ~OPEN_CREATE & ~OPEN_APPEND), S_IREAD | S_IWRITE);
				if (handleAsync == -1 || errno == EINVAL) {	// fall back to cached access
					#ifdef SEQAN_DEBUG_OR_TEST_
						::std::cerr << "Warning: Direct access openening failed. (" << strerror(errno) << ")" << ::std::endl;
					#endif
					handleAsync = handle;
				}
				#ifdef SEQAN_DEBUG_OR_TEST_
				    else
						::std::cerr << "Direct access successfully initiated" << ::std::endl;
				#endif
			} else
				handleAsync = handle;
			
			SEQAN_PROADD(PROOPENFILES, 1);
            return true;
        }

        bool openTemp(int openMode = OPEN_RDWR + OPEN_CREATE) {
#ifdef SEQAN_DEFAULT_TMPDIR
			char tmpFileName[] = SEQAN_DEFAULT_TMPDIR "/GNDXXXXXXX";
#else
			char tmpFileName[] = "/var/tmp/GNDXXXXXXX";
#endif
			if ((handle = handleAsync = ::mkstemp(tmpFileName)) == -1) {
				::std::cerr << "Couldn't create temporary file " << tmpFileName << ". (" << strerror(errno) << ")" << ::std::endl;
				return false;
			}
			if (!(close() && open(tmpFileName, openMode))) return false;
            #ifdef SEQAN_DEBUG
				if (::unlink(tmpFileName) == -1)
					::std::cerr << "Couldn't unlink temporary file " << tmpFileName << ". (" << strerror(errno) << ")" << ::std::endl;
            #else
				::unlink(tmpFileName);
			#endif
			return true;
        }

        inline bool close() {
			bool result = true;
			if (handle != handleAsync)
	            result &= (::close(handleAsync) == 0);
            result &= (::close(handle) == 0);
            handleAsync = -1;
            handle = -1;
			SEQAN_PROSUB(PROOPENFILES, 1);
            return result;
        }
    };


    //////////////////////////////////////////////////////////////////////////////
    // (SeqAn adaption)
    //////////////////////////////////////////////////////////////////////////////
/*
    template < typename TConfig >
    struct aQueue<File<Async<TConfig> > >
    {
        typedef void* Type;
    };
*/
    template < typename TConfig >
    struct aRequest<File<Async<TConfig> > >
    {
        typedef aiocb Type;
    };
/*
    template < typename TConfig >
    struct aEvent<File<Async<TConfig> > >
    {
        typedef aiocb Type;
    };
*/

    //////////////////////////////////////////////////////////////////////
    // event based read/write

//    enum { _AsyncIOSignal = SIGIO };

	inline void printRequest(aiocb &request) {
		::std::cout << ::std::hex;
		::std::cout << "fildes:  " << (unsigned int) request.aio_fildes << ::std::endl;
		::std::cout << "buffer:  " << (unsigned int) request.aio_buf << ::std::endl;
		::std::cout << "offset:  " << (unsigned int) request.aio_offset<< ::std::endl;
		::std::cout << "nbytes:  " << (unsigned int) request.aio_nbytes << ::std::endl;
		::std::cout << "event:   " << (unsigned int) request.aio_sigevent.sigev_notify << ::std::endl;
		::std::cout << ::std::dec;
	}

    template < typename TConfig, typename TValue, typename TSize, typename TPos >
    bool areadAt(File<Async<TConfig> > & me, TValue *memPtr, TSize const count, TPos const fileOfs,
        aiocb &request)
    {
        SEQAN_PROTIMESTART(tw);
        memset(&request, 0, sizeof(aiocb));
        request.aio_fildes = me.handleAsync;
        request.aio_buf = memPtr;
        request.aio_offset = fileOfs;
        request.aio_offset *= sizeof(TValue);
        request.aio_nbytes = count * sizeof(TValue);
        request.aio_sigevent.sigev_notify = SIGEV_NONE;
/*        request.aio_sigevent.sigev_notify = SIGEV_SIGNAL;
        request.aio_sigevent.sigev_signo = _AsyncIOSignal;
        request.aio_sigevent.sigev_value.sival_ptr = &request;*/
        SEQAN_PROADD(PROIO, (request.aio_nbytes + PROPAGESIZE - 1) / PROPAGESIZE);
		int result = aio_read(&request);
        SEQAN_PROADD(PROIWAIT, SEQAN_PROTIMEDIFF(tw));
        #ifdef SEQAN_DEBUG
			if (result) ::std::cerr << "areadAt returned " << result << ::std::endl;
		#endif
		return result == 0;
    }
    
    template < typename TConfig, typename TValue, typename TSize, typename TPos >
    bool awriteAt(File<Async<TConfig> > & me, const TValue *memPtr, TSize const count, TPos const fileOfs,
        aiocb &request)
    {
        SEQAN_PROTIMESTART(tw);
        memset(&request, 0, sizeof(aiocb));
        request.aio_fildes = me.handleAsync;
        request.aio_buf = const_cast<TValue*>(memPtr);
        request.aio_offset = fileOfs;
        request.aio_offset *= sizeof(TValue);
        request.aio_nbytes = count * sizeof(TValue);
        request.aio_sigevent.sigev_notify = SIGEV_NONE;
/*        request.aio_sigevent.sigev_notify = SIGEV_SIGNAL;
        request.aio_sigevent.sigev_signo = _AsyncIOSignal;
        request.aio_sigevent.sigev_value.sival_ptr = &request;*/
        SEQAN_PROADD(PROIO, (request.aio_nbytes + PROPAGESIZE - 1) / PROPAGESIZE);
		int result = aio_write(&request);
        SEQAN_PROADD(PROIWAIT, SEQAN_PROTIMEDIFF(tw));
        #ifdef SEQAN_DEBUG
			if (result) ::std::cerr << "awriteAt returned " << result << ::std::endl;
		#endif
        return result == 0;
    }

    template < typename TConfig >
    inline bool flush(File<Async<TConfig> > & me) {
		return me.handle == me.handleAsync || fdatasync(me.handle) == 0;
    }

    //////////////////////////////////////////////////////////////////////
    // queue specific functions

	inline bool waitFor(aiocb &request) {
        aiocb * cblist = &request;
        SEQAN_PROTIMESTART(tw);
		int result = aio_suspend(&cblist, 1, NULL);
        SEQAN_PROADD(PROCWAIT, SEQAN_PROTIMEDIFF(tw));
        #ifdef SEQAN_DEBUG
			if (result) {
	 			int eno = aio_error(&request);
				if (eno != EINPROGRESS)
					::std::cerr << "waitFor: aio_error returned " << strerror(eno) << " and errno is " << strerror(errno) << ::std::endl;
			}
		#endif
		return result == 0;
	}

	inline bool waitFor(aiocb &request, long timeout_millis) {
        aiocb * cblist = &request;
        timespec ts;
        ts.tv_sec = timeout_millis / 1000;
        ts.tv_nsec = (timeout_millis % 1000) * 1000;
        SEQAN_PROTIMESTART(tw);
		int result = aio_suspend(&cblist, 1, &ts);
        SEQAN_PROADD(PROCWAIT, SEQAN_PROTIMEDIFF(tw));
        #ifdef SEQAN_DEBUG
			if (result) {
	 			int eno = aio_error(&request);
				if (eno != EINPROGRESS)
					::std::cerr << "waitFor: aio_error returned " << strerror(eno) << " and errno is " << strerror(errno) << ::std::endl;
			}
		#endif
        return result == 0;
	}

	template < typename TSize >
	inline TSize waitForAny(aiocb const * const contexts[], TSize count) {
        SEQAN_PROTIMESTART(tw);
		bool result = aio_suspend(contexts, count, NULL);
        SEQAN_PROADD(PROCWAIT, SEQAN_PROTIMEDIFF(tw));
        return result == 0;
	}

	template < typename TSize >
	inline TSize waitForAny(aiocb const * const contexts[], TSize count, long timeout_millis) {
        timespec ts;
        ts.tv_sec = timeout_millis / 1000;
        ts.tv_nsec = (timeout_millis % 1000) * 1000;
        SEQAN_PROTIMESTART(tw);
		bool result = aio_suspend(contexts, count, &ts);
        SEQAN_PROADD(PROCWAIT, SEQAN_PROTIMEDIFF(tw));
        return result == 0;
	}

    template < typename TConfig >
    inline bool cancel(File<Async<TConfig> > & me, aiocb &request) {
        return aio_cancel(me.handleAsync, &request) == 0;
    }

    inline int error(aiocb const &request) {
        return aio_error(&request);
    }

    inline int return_value(aiocb &request) {
        return aio_return(&request);
    }

    template < typename TConfig, typename >
    inline void release(File<Async<TConfig> > & me, aiocb const &request) {
    }
/*
    typedef void (*sighandler_t)(int);
    static unsigned _AsyncIOHandlerRefCount = 0;
    static struct sigaction _AsyncIOOldSig;

    inline void _AsyncIOHandler(int sigNo, siginfo_t *info, void *hint) {
        SEQAN_ASSERT(sigNo == _AsyncIOSignal);
        // TODO: signal respective event
        // currently we don't need async IO handlers because
        // we only wait for single events
    }

    static sighandler_t _addAsyncIOHandler() {
        struct sigaction newSig, oldSig;
        newSig.sa_sigaction = _AsyncIOHandler;
        sigemptyset(&newSig.sa_mask);
        newSig.sa_flags = SA_RESTART + SA_SIGINFO;
        if (sigaction(_AsyncIOSignal, &newSig, &oldSig) < 0)
            return SIG_ERR;
        return oldSig.sa_handler;
    }
*/
    
	//////////////////////////////////////////////////////////////////////////////
	// page aligned allocate for direct file io

    struct TagAllocateAligned_;	//< allocate page aligned memory for direct i/o access
    typedef Tag<TagAllocateAligned_> const TagAllocateAligned;

	template <typename T, typename TValue, typename TSize>
	inline void
	allocate(T const & me, 
			 TValue * & data,
			 TSize count,
			 TagAllocateAligned const)
	{
        int error = posix_memalign(reinterpret_cast<void**>(&data), sysconf(_SC_PAGESIZE), count * sizeof(TValue));
        if (error) {
			::std::cerr << "AlignAllocator: Could not allocate memory of size " << ::std::hex << 
				count * sizeof(TValue) << " and an alignment of " << 
				(unsigned) sysconf(_SC_PAGESIZE) <<	::std::dec << ". (ErrNo=" << error << ")" << ::std::endl;
            data = NULL;
        } else
            SEQAN_PROADD(PROMEMORY, sizeof(TValue) * count);
	}

	//////////////////////////////////////////////////////////////////////////////
	// page aligned deallocate for direct file io

	template <typename T, typename TValue, typename TSize>
	inline void 
	deallocate( T const & me,
				TValue * data, 
				TSize count,
				TagAllocateAligned const)
	{
        if (data) {
        	SEQAN_PROSUB(PROMEMORY, sizeof(TValue) * count);
			free(data);
		}
	}

#endif

    //////////////////////////////////////////////////////////////////////////////
    // global functions

    template < typename TConfig >
    struct Value< File< Async<TConfig> > >
    {
            typedef typename TConfig::Type Type;
    };

    template < typename TConfig >
    struct Size< File< Async<TConfig> > >
    {
        typedef typename File< Async<TConfig> >::SizeType Type;
    };

    template < typename TConfig >
    struct Position< File< Async<TConfig> > >
    {
        typedef typename File< Async<TConfig> >::FilePtr Type;
    };

    template < typename TConfig >
    struct Difference< File< Async<TConfig> > >
    {
        typedef typename File< Async<TConfig> >::FilePtr Type;
    };



    template <typename TConfig, typename TValue, typename TSize>
	inline void
	allocate( File<Async<TConfig> > const & me,
			  TValue * & data, 
			  TSize count)
	{
		allocate(me, data, count, TagAllocateAligned());
	}

    template <typename TConfig, typename TValue, typename TSize>
	inline void
	deallocate( File<Async<TConfig> > const & me,
				TValue * data, 
				TSize count)
	{
		deallocate(me, data, count, TagAllocateAligned());
	}


}

#endif
