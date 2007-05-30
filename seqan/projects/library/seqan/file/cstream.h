#ifndef SEQAN_HEADER_CSTREAM_H
#define SEQAN_HEADER_CSTREAM_H
 
#include <cstdio>

namespace SEQAN_NAMESPACE_MAIN
{
//////////////////////////////////////////////////////////////////////////////
/**
.Adaption."std::FILE *":
..summary:Standard library C style streams.
*/

//////////////////////////////////////////////////////////////////////////////
// Position is now defined in file/file_cstyle.h
/*
template <>
struct Position<FILE *>
{
	typedef long Type;
};
*/
//////////////////////////////////////////////////////////////////////////////

template <>
struct Value<FILE *>
{
	typedef char Type;
};

//////////////////////////////////////////////////////////////////////////////

template <>
struct Position<FILE *>
{
	typedef ::std::fpos_t Type;
};

//////////////////////////////////////////////////////////////////////////////

template <typename T>
struct _IsTellSeekStream;

template <>
struct _IsTellSeekStream<FILE *>
{
	typedef True Type;
};

//////////////////////////////////////////////////////////////////////////////

///.Internal._streamEOF.param.stream.type:Adaption."std::FILE *"

inline bool 
_streamEOF(::std::FILE * me)
{
SEQAN_CHECKPOINT
	return feof(me) || ferror(me);
}

//////////////////////////////////////////////////////////////////////////////

///.Internal._streamRead.param.stream.type:Adaption."std::FILE *"

template <typename TValue>
inline size_t 
_streamRead(TValue * target,
			::std::FILE * source,
			size_t limit)
{
SEQAN_CHECKPOINT
	return ::std::fread(target, sizeof(TValue), limit, source);
}

//////////////////////////////////////////////////////////////////////////////

///.Internal._streamGet.param.stream.type:Adaption."std::FILE *"

inline char 
_streamGet(::std::FILE * source)
{
SEQAN_CHECKPOINT
	return getc(source);
}

//////////////////////////////////////////////////////////////////////////////

///.Internal._streamPut.param.stream.type:Adaption."std::FILE *"

inline void
_streamPut(::std::FILE * target,
		   char character)
{
SEQAN_CHECKPOINT
	putc(character, target);
}


//////////////////////////////////////////////////////////////////////////////

///.Internal._streamPut.param.stream.type:Adaption."std::FILE *"

//////////////////////////////////////////////////////////////////////////////

///.Internal._streamTellG.param.stream.type:Adaption."std::FILE *"

inline Position<FILE *>::Type
_streamTellG(FILE * me)
{
SEQAN_CHECKPOINT
	Position<FILE *>::Type pos;
	::std::fgetpos(me, &pos);
	return pos;
}

//////////////////////////////////////////////////////////////////////////////

///.Internal._streamTellP.param.stream.type:Adaption."std::FILE *"

inline Position<FILE *>::Type
_streamTellP(FILE * me)
{
SEQAN_CHECKPOINT
	Position<FILE *>::Type pos;
	::std::fgetpos(me, &pos);
	return pos;
}

//////////////////////////////////////////////////////////////////////////////

///.Internal._streamSeekG.param.stream.type:Adaption."std::FILE *"

inline void
_streamSeekG(FILE * me,
			 Position<FILE *>::Type pos)
{
SEQAN_CHECKPOINT
	::std::fsetpos(me, &pos);
}

//////////////////////////////////////////////////////////////////////////////

///.Internal._streamSeekP.param.stream.type:Adaption."std::FILE *"

inline void
_streamSeekP(FILE * me,
			 Position<FILE *>::Type pos)
{
SEQAN_CHECKPOINT
	::std::fsetpos(me, &pos);
}

//////////////////////////////////////////////////////////////////////////////

///.Internal._streamSeek2G.param.stream.type:Adaption."std::FILE *"

inline void
_streamSeek2G(FILE * me,
	 int off)
{
SEQAN_CHECKPOINT
	::std::fseek(me, off, SEEK_CUR);
}

//////////////////////////////////////////////////////////////////////////////

///.Internal._streamUnget.param.stream.type:Adaption."std::FILE *"

inline void
_streamUnget(::std::FILE * stream)
{
SEQAN_CHECKPOINT
	_streamSeek2G(stream, -1);
}

//////////////////////////////////////////////////////////////////////////////
// Stream operators for FILE *
//////////////////////////////////////////////////////////////////////////////

// ISO C++ operators are only allowed for classes, not for pointers

/*
template <typename TSource>
inline FILE *
operator << (FILE * target, 
			 TSource & source)
{
SEQAN_CHECKPOINT
	write(target, source);
	return target;
}
template <typename TSource>
inline FILE *
operator << (FILE * target, 
			 TSource const & source)
{
SEQAN_CHECKPOINT
	write(target, source);
	return target;
}

//____________________________________________________________________________

template <typename TTarget>
inline FILE *
operator >> (FILE * source, 
			 TTarget & target)
{
SEQAN_CHECKPOINT
	read(source, target);
	return source;
}
template <typename TTarget>
inline FILE *
operator >> (FILE * source, 
			 TTarget const & target)
{
SEQAN_CHECKPOINT
	read(source, target);
	return source;
}
*/

//////////////////////////////////////////////////////////////////////////////

} //namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
