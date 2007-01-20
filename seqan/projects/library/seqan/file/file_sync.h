/*
 *  file_simple.h
 *  genindex
 *
 *  Created by David Weese on 17.07.05.
 *
 */

#ifndef SEQAN_HEADER_FILE_SIMPLE_H
#define SEQAN_HEADER_FILE_SIMPLE_H

#include <fcntl.h>          // O_CREAT ..
#include <sys/stat.h>       // 
#include <cstdio>           // tmpnam(..)

#ifdef PLATFORM_WINDOWS
# include <io.h>            // read(..) ..
#else
# include <cstdlib>
# include <cerrno>
# include <unistd.h>
#endif


//////////////////////////////////////////////////////////////////////////////

namespace SEQAN_NAMESPACE_MAIN
{

    //////////////////////////////////////////////////////////////////////////////
    // interface for typed files - deprecated
//     template < typename T = unsigned char >
//     struct SimpleConfig;
// 
//     template < typename TConfig = SimpleConfig<> >
//     struct Sync;


#ifdef PLATFORM_WINDOWS

    //////////////////////////////////////////////////////////////////////////////
    // Windows rtl file access
    template < typename TConfig >
	class File<Sync<TConfig> >
    {

    public:

		typedef __int64			FilePtr;
		typedef __int64         SizeType;   // type of file size
        typedef unsigned int    _SizeType;  // type of transfer size (for read or write)
		typedef int				Handle;

        File(void *dummy = NULL): // to be compatible with the FILE*(NULL) constructor
            handle(-1) {}

/*        File(char const *fileName = NULL, int openMode = OPEN_RDWR + OPEN_CREATE | OPEN_APPEND):
            handle(-1)
        {
            if (fileName)
                open(fileName, openMode);
            else
                openTemp();
        }

        File(Nothing *queue, char const *fileName = NULL, int openMode = OPEN_RDWR + OPEN_CREATE | OPEN_APPEND):
            handle(-1)
        {
            if (fileName)
                open(fileName, openMode);
            else
                openTemp();
        }*/

		~File() {
//			close();            
		}

        bool open(char const *fileName, int openMode = OPEN_RDWR + OPEN_CREATE | OPEN_APPEND) {
            handle = _open(fileName, getOFlag(openMode), _S_IREAD | _S_IWRITE);
			if (handle == -1) {
				printf("Open failed on file %s: %s\n", fileName, strerror(errno));
				return false;
			}
			SEQAN_PROADD(PROOPENFILES, 1);
            return true;
        }

        bool openTemp(int openMode = OPEN_RDWR + OPEN_CREATE) {
			char *fileName = _tempnam(SEQAN_DEFAULT_TMPDIR, "GNDX");
			if (!fileName) {
				printf("Cannot create a unique temp. filename\n");
				return false;
			}
            return open(fileName, openMode | OPEN_TEMPORARY);
        }

        inline bool close() {
            if (_close(handle) != 0)
                return false;
            handle = -1;
			SEQAN_PROSUB(PROOPENFILES, 1);
            return true;
        }

		inline int read(void *buffer, _SizeType count) const {
            SEQAN_PROADD(PROIO, (count + PROPAGESIZE - 1) / PROPAGESIZE);
            SEQAN_PROTIMESTART(tw);
		    int result = _read(handle, buffer, count);
            SEQAN_PROADD(PROCWAIT, SEQAN_PROTIMEDIFF(tw));
            return result;
		}

		inline int write(void const *buffer, _SizeType count) const {
            SEQAN_PROADD(PROIO, (count + PROPAGESIZE - 1) / PROPAGESIZE);
            SEQAN_PROTIMESTART(tw);
		    int result = _write(handle, buffer, count);
            SEQAN_PROADD(PROCWAIT, SEQAN_PROTIMEDIFF(tw));
            return result;
		}

		inline FilePtr seek(FilePtr pos, int origin = SEEK_SET) const {
			return _lseeki64(handle, pos, origin);
		}

		inline FilePtr tell() const {
			return _telli64(handle);
		}

		static int error() {
			return errno;
		}

        operator bool () const {
            return handle != -1;
        }

    protected:

        Handle handle;

        inline int getOFlag(int openMode) const {
			int result;

			switch (openMode & OPEN_MASK) {
                case OPEN_RDONLY:
                    result = _O_RDONLY;
					break;
                case OPEN_WRONLY:
                    result = _O_WRONLY;
					break;
                case OPEN_RDWR:
                    result = _O_RDWR;
					break;
			}

			if (openMode & OPEN_CREATE)     result |= _O_CREAT;
			if (openMode & OPEN_APPEND)     result |= _O_APPEND;
            if (openMode & OPEN_TEMPORARY)  result |= _O_TEMPORARY;
			return result | _O_BINARY;
        }

    };

	inline bool fileExists(const char *fileName) {
		struct _stat buf;
		return _stat(fileName, &buf) == 0;
	}

	inline bool fileUnlink(const char *fileName) {
		return _unlink(fileName) == 0;
	}

#else

    //////////////////////////////////////////////////////////////////////////////
    // Unix file access
    template < typename TConfig >
	class File<Sync<TConfig> >
    {

    public:

		typedef off_t			FilePtr;
		typedef off_t           SizeType;   // type of file size
        typedef size_t          _SizeType;  // type of transfer size (for read or write)
		typedef int				Handle;

        Handle handle;

        File(void *dummy = NULL): // to be compatible with the FILE*(NULL) constructor
            handle(-1) {}

/*        File(char const *fileName = NULL, int openMode = OPEN_RDWR + OPEN_CREATE | OPEN_APPEND):
            handle(-1)
        {
            if (fileName)
                open(fileName, openMode);
            else
                openTemp();
        }

        File(Nothing *queue, char const *fileName = NULL, int openMode = OPEN_RDWR + OPEN_CREATE | OPEN_APPEND):
            handle(-1)
        {
            if (fileName)
                open(fileName, openMode);
            else
                openTemp();
        }*/

		~File() {
//            close();            
		}

        bool open(char const *fileName, int openMode = OPEN_RDWR + OPEN_CREATE | OPEN_APPEND) {
            handle = ::open(fileName, getOFlag(openMode), S_IREAD | S_IWRITE);
			if (handle == -1 && errno == EINVAL) {	// fall back to cached access
	            #ifdef SEQAN_DEBUG_OR_TEST_
					printf("Warning: Direct access openening failed: %s.\n", fileName);
				#endif			
          	    handle = ::open(fileName, getOFlag(openMode & ~OPEN_ASYNC), S_IREAD | S_IWRITE);
			}
			
			if (handle == -1) {
				printf("Open failed on file %s: %s\n", fileName, strerror(errno));
				return false;
			}
			SEQAN_PROADD(PROOPENFILES, 1);
            return true;
        }

        bool openTemp(int openMode = OPEN_RDWR + OPEN_CREATE) {
			char tmpFileName[] = "/GNDXXXXXXX";
			if ((handle = ::mkstemp(tmpFileName)) == -1) {
				printf("Cannot create temp. file %s\n", strerror(errno));
				return false;
			}
			if (!(close() && open(tmpFileName, openMode))) return false;
			int result = ::unlink(tmpFileName);
            #ifdef SEQAN_DEBUG
				if (result == -1) printf("Cannot unlink temp. file: %s\n", strerror(errno));
            #endif
			return true;
        }

        inline bool close() {
            if (::close(handle) == -1) return false;
            handle = -1;
			SEQAN_PROSUB(PROOPENFILES, 1);
            return true;
        }

		inline ssize_t read(void *buffer, _SizeType count) const {
            SEQAN_PROADD(PROIO, (count + PROPAGESIZE - 1) / PROPAGESIZE);
            SEQAN_PROTIMESTART(tw);
		    ssize_t result = ::read(handle, buffer, count);
            SEQAN_PROADD(PROCWAIT, SEQAN_PROTIMEDIFF(tw));
            return result;
		}

		inline ssize_t write(void const *buffer, _SizeType count) const {
            SEQAN_PROADD(PROIO, (count + PROPAGESIZE - 1) / PROPAGESIZE);
            SEQAN_PROTIMESTART(tw);
		    ssize_t result = ::write(handle, buffer, count);
            SEQAN_PROADD(PROCWAIT, SEQAN_PROTIMEDIFF(tw));
            return result;
		}
/*
		inline ssize_t readAt(void *buffer, _SizeType count, FilePtr offset) const {
            return ::pread(handle, buffer, count, offset);
		}

		inline ssize_t writeAt(void const *buffer, _SizeType count, FilePtr offset) const {
    	    return ::pwrite(handle, buffer, count, offset);
		}
*/
		inline FilePtr seek(FilePtr pos, int origin = SEEK_SET) const {
            FilePtr result = ::lseek(handle, pos, origin);
//			#ifdef SEQAN_DEBUG
				if (result < 0) printf("seek returned %d %s\n", result, strerror(errno));
//			#endif
			return result;
		}

		inline FilePtr tell() const {
            return seek(0, SEEK_CUR);
        }

		static int error() {
            return errno;
		}

        operator bool () const {
            return handle != -1;
        }

    protected:

        inline int getOFlag(int openMode) const {
			int result = O_LARGEFILE;

			switch (openMode & OPEN_MASK) {
                case OPEN_RDONLY:
                    result |= O_RDONLY;
					break;
                case OPEN_WRONLY:
                    result |= O_WRONLY;
					break;
                case OPEN_RDWR:
                    result |= O_RDWR;
					break;
			}

			if (openMode & OPEN_CREATE)     result |= O_CREAT;
			if (openMode & OPEN_APPEND)     result |= O_APPEND;
//			if (openMode & OPEN_TEMPORARY)  result |= O_TEMPORARY;
        #ifdef SEQAN_DIRECTIO
    		if (openMode & OPEN_ASYNC)		result |= O_DIRECT;
        #endif
			return result;
        }

    };

	inline bool fileExists(const char *fileName) {
		struct stat buf;
		return stat(fileName, &buf) != -1;
	}

	inline bool fileUnlink(const char *fileName) {
		return unlink(fileName) == 0;
	}

#endif

    //////////////////////////////////////////////////////////////////////////////
    // global functions

    template < typename TConfig >
    struct Value< File<Sync<TConfig> > >
    {
	    typedef typename TConfig::Type Type;
    };

    template < typename TConfig >
    struct Size< File<Sync<TConfig> > >
    {
        typedef typename File<Sync<TConfig> >::SizeType Type;
    };

    template < typename TConfig >
    struct Position< File<Sync<TConfig> > >
    {
        typedef typename File<Sync<TConfig> >::FilePtr Type;
    };

    template < typename TConfig >
    struct Difference< File<Sync<TConfig> > >
    {
        typedef typename File<Sync<TConfig> >::FilePtr Type;
    };

    template < typename TConfig, typename TValue, typename TSize >
    inline bool read(File<Sync<TConfig> > & me, TValue *memPtr, TSize const count) {
		return me.read(memPtr, count * sizeof(TValue)) == count * sizeof(TValue);
    }
    
    template < typename TConfig, typename TValue, typename TSize >
    inline bool write(File<Sync<TConfig> > & me, TValue const *memPtr, TSize const count) {
		return me.write(memPtr, count * sizeof(TValue)) == count * sizeof(TValue);
    }

}

#endif
