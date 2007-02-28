/*
 *  file_page_raid0.h
 *  genindex
 *
 *  Created by David Weese on 17.07.05.
 *
 */

#ifndef SEQAN_HEADER_FILE_PAGE_RAID0_H
#define SEQAN_HEADER_FILE_PAGE_RAID0_H


//////////////////////////////////////////////////////////////////////////////

namespace SEQAN_NAMESPACE_MAIN
{


    //////////////////////////////////////////////////////////////////////////////
    // page based read/write for striped files

	template < typename TValue, unsigned _FileCount, typename TFile, typename TSpec > inline
	bool readPage(
		int pageNo, 
		PageFrame<TValue, File< Striped<_FileCount, TFile> >, TSpec> &pf, 
		File< Striped<_FileCount, TFile> > &file)
	{
		typedef typename Position<TFile>::Type pos_t;
		#ifdef SEQAN_VERBOSE
			printf("readPage:  %x from page %d\n", (unsigned)(TValue*)pf.begin, pageNo);
		#endif
		pf.dirty = false;
		pf.status = pf.READING;
		return areadAt(
			file[pageNo % _FileCount], 
			(TValue*)pf.begin, 
			size(pf), 
			(pos_t)(pageNo / _FileCount) * (pos_t)pageSize(pf), 
			pf.request);
	}

	template < typename TValue, unsigned _FileCount, typename TFile, typename TSpec > inline
	bool writePage(
		PageFrame<TValue, File< Striped<_FileCount, TFile> >, TSpec> &pf, 
		int pageNo, 
		File< Striped<_FileCount, TFile> > &file)
	{
		typedef typename Position<TFile>::Type pos_t;
		#ifdef SEQAN_VERBOSE
			printf("writePage: %x from page %d\n", (unsigned)(TValue*)pf.begin, pageNo);
		#endif
		pf.status = pf.WRITING;
		return awriteAt(
			file[pageNo % _FileCount], 
			(TValue*)pf.begin, 
			size(pf), 
			(pos_t)(pageNo / _FileCount) * (pos_t)pageSize(pf), 
			pf.request);
	}

	template < typename TValue, unsigned _FileCount, typename TFile, typename TSpec, typename TSize > inline
    bool readLastPage(
		int pageNo, 
		PageFrame<TValue, File< Striped<_FileCount, TFile> >, TSpec> &pf, 
		File< Striped<_FileCount, TFile> > &file,
		TSize size)
	{
		typedef typename Position<TFile>::Type pos_t;
		#ifdef SEQAN_VERBOSE
			printf("readPage:  %x from page %d size %d\n", (unsigned)(TValue*)pf.begin, pageNo, size);
		#endif
		pf.dirty = false;
		pf.status = pf.READY;
		return readAt(
			file[pageNo % _FileCount], 
			(TValue*)pf.begin, 
			size, 
			(pos_t)(pageNo / _FileCount) * (pos_t)pageSize(pf));
	}

	template < typename TValue, unsigned _FileCount, typename TFile, typename TSpec, typename TSize > inline
	bool writeLastPage(
		PageFrame<TValue, File< Striped<_FileCount, TFile> >, TSpec> &pf, 
		int pageNo, 
		File< Striped<_FileCount, TFile> > &file,
		TSize size)
	{
		typedef typename Position<TFile>::Type pos_t;
		#ifdef SEQAN_VERBOSE
			printf("writePage: %x from page %d size %d\n", (unsigned)(TValue*)pf.begin, pageNo, size);
		#endif
		pf.dirty = false;
		pf.status = pf.READY;
//        resize(pf, size);
		return writeAt(
			file[pageNo % _FileCount], 
			(TValue*)pf.begin, 
			size, 
			(pos_t)(pageNo / _FileCount) * (pos_t)pageSize(pf));
	}


	//////////////////////////////////////////////////////////////////////////////
	// bucket based read/write methods for striped files

	template < typename TValue, unsigned _FileCount, typename TFile > inline
	unsigned readBucket(
		PageBucket<TValue> &b, 
		int pageNo, 
		unsigned pageSize, 
		unsigned dataSize, 
		File< Striped<_FileCount, TFile> > &file) 
	{
		typedef typename Position<TFile>::Type pos_t;
        unsigned readSize = Min(dataSize - b.pageOfs, (unsigned)(b.end - b.begin));
		#ifdef SEQAN_VERBOSE
			printf("readBucket:  %x from page %d at %d size %d\n", (unsigned)b.begin, pageNo, pageNo * pageSize + b.pageOfs, readSize);
		#endif
        if (readSize && readAt(file[pageNo % _FileCount], b.begin, readSize, (pos_t)(pageNo / _FileCount) * (pos_t)pageSize + b.pageOfs)) {
            b.pageOfs += readSize;
            b.cur = b.begin;
            b.end = b.begin + readSize;
            return readSize;
        } else
            return 0;
	}

	template < typename TValue, unsigned _FileCount, typename TFile > inline
	bool writeBucket(
		PageBucket<TValue> &b,
		int pageNo, 
		unsigned pageSize, 
		File< Striped<_FileCount, TFile> > &file) 
	{
		typedef typename Position<TFile>::Type pos_t;
		#ifdef SEQAN_VERBOSE
			printf("writeBucket:  %x from page %d at %d size %d\n", (unsigned)b.begin, pageNo, pageNo * pageSize + b.pageOfs, b.cur - b.begin);
		#endif
        if ((b.cur == b.begin) || writeAt(file[pageNo % _FileCount], b.begin, b.cur - b.begin, (pos_t)(pageNo / _FileCount) * (pos_t)pageSize + b.pageOfs)) {
            b.pageOfs += b.cur - b.begin;
            b.cur = b.begin;
            return true;
        } else
            return false;
	}

	template < typename TValue, unsigned _FileCount, typename TFile, typename TSpec > inline
	bool writeBucket(
		PageFrame<TValue, File< Striped<_FileCount, TFile> >, Dynamic<TSpec> > &pf, 
		unsigned &pageOfs, 
		File< Striped<_FileCount, TFile> > &file) 
	{
		typedef typename Position<TFile>::Type pos_t;
		#ifdef SEQAN_VERBOSE
			printf("writeBucket:  %x from page %d at %d size %d\n", (unsigned)pf.begin, pf.pageNo, pf.pageNo * pageSize(pf) + pageOfs, size(pf));
		#endif
        if (pf.end == pf.begin) return true;
        if (awriteAt(file[pf.pageNo % _FileCount], pf.begin, size(pf), (pos_t)(pf.pageNo / _FileCount) * (pos_t)pageSize(pf) + pageOfs, pf.request)) {
            pf.status = pf.WRITING;
            pageOfs += size(pf);
            return true;
        } else
            return false;
	}

}

#endif
