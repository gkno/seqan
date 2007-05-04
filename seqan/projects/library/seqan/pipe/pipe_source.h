/*
 *  pipe_source.h
 *  SeqAn
 *
 *  Created by David Weese on 17.07.05.
 *
 */

#ifndef SEQAN_HEADER_PIPE_SOURCE_H
#define SEQAN_HEADER_PIPE_SOURCE_H

namespace SEQAN_NAMESPACE_MAIN
{

	template < typename TSpec = void >
	struct Source;
/*
	template < typename TValue >
    struct Size< Pipe< TValue*, Source<> > > {
		typedef ptrdiff_t Type;
	};
*/

/**
.Spec.Source:
..cat:Pipelining
..general:Class.Pipe
..summary:Pipelining adaptor for arbitrary containers or iterators.
..signature:Pipe<TInput, Source<> >
..param.TInput:The type of container or iterator this module reads from.
*/


	//////////////////////////////////////////////////////////////////////////////
	// source class
	template < typename TInput, typename TSpec >
	struct Pipe< TInput, Source<TSpec> >
    {
		TInput const &in;
		typename Iterator<TInput const>::Type cur;

		Pipe(typename _RemoveConst<TInput>::Type &_cont):
			in(_cont) {}

		Pipe(TInput const &_cont):
			in(_cont) {}

        inline typename Value<TInput>::Type const & operator*() const {
            return *cur;
        }
    
        inline Pipe& operator++() {
            ++cur;
            return *this;
        }
    };

    template < typename TInput, typename TSpec >
	struct Iterator< Pipe< TInput, Source<TSpec> > > {
		typedef typename Iterator<TInput const>::Type Type;
	};

    template < typename TInput, typename TSpec >
	inline typename Iterator< Pipe< TInput, Source<TSpec> > >::Type
	begin(Pipe< TInput, Source<TSpec> > &pipe) {
		return begin(pipe.in);
	}

    template < typename TInput, typename TSpec >
	inline typename Iterator< Pipe< TInput, Source<TSpec> > >::Type
	end(Pipe< TInput, Source<TSpec> > &pipe) {
		return end(pipe.in);
	}



	//////////////////////////////////////////////////////////////////////////////
	// simple buffer adaption for fast random accessable containers
    template < typename TContainer >
	struct ContainerBuffer 
	{
		TContainer *cont;

        ContainerBuffer(): cont(NULL) {}
		ContainerBuffer(TContainer &_cont): cont(&_cont) {}

		template <typename TSize>
		inline typename Reference<TContainer>::Type 
		operator[] (TSize i) {
			return (*cont)[i];
		}
    };

    template < typename TContainer >
	struct Iterator< ContainerBuffer<TContainer> >:
		Iterator<TContainer> {};

    template < typename TContainer  >
	struct Value< ContainerBuffer<TContainer> >:
		Value<TContainer> {};

    template < typename TContainer >
	struct Size< ContainerBuffer<TContainer> >:
		Size<TContainer> {};


	template < typename TContainer >
    inline typename Size< ContainerBuffer<TContainer> >::Type
    length(ContainerBuffer<TContainer> const &me) {
		if (me.cont)
			return length(*me.cont);
		else
			return 0;
    }


	//////////////////////////////////////////////////////////////////////////////
	// simple buffer adaption for fast random accessable iterators
    template < typename TIterator >
	struct IteratorBuffer 
	{
		TIterator						begin;
		typename Size<TIterator>::Type	size;

        IteratorBuffer(): begin(NULL), size(0) {}
		template <typename TSize>
		IteratorBuffer(TIterator &_begin, TSize _size): begin(&_begin), size(_size) {}

		template <typename TSize>
		inline typename Reference<TIterator>::Type 
		operator[] (TSize i) {
			return *(begin + i);
		}
    };

    template < typename TIterator >
	struct Iterator< IteratorBuffer<TIterator> >:
		Iterator<TIterator> {};

    template < typename TIterator  >
	struct Value< IteratorBuffer<TIterator> >:
		Value<TIterator> {};

    template < typename TIterator >
	struct Size< IteratorBuffer<TIterator> >:
		Size<TIterator> {};


	template < typename TIterator >
    inline typename Size< IteratorBuffer<TIterator> >::Type
    length(IteratorBuffer<TIterator> const &me) {
		return me.size;
    }


    //////////////////////////////////////////////////////////////////////////////
	// handler that emulates a buffer by using a container buffer
    struct _SourceNonCachingSpec;
	typedef Tag<_SourceNonCachingSpec> SourceNonCachingSpec;

	template < typename TPipe >
	struct BufferHandler< TPipe, SourceNonCachingSpec >
    {
        typedef typename Source<TPipe>::Type	TSource;
        typedef ContainerBuffer<TSource const>	TBuffer;

		TPipe	&pipe;

		BufferHandler(TPipe &_pipe):
			pipe(_pipe) {}

		template <typename TSize>
		BufferHandler(TPipe &_pipe, TSize):
			pipe(_pipe) {}

        inline TBuffer first() {
            return TBuffer(pipe.in);
        }

        inline TBuffer next() {
            return TBuffer();
        }

        inline void process() {}
        inline void end() {}
        inline void cancel() {}
    };

	template < typename TPipe >
    struct Value< BufferHandler< TPipe, SourceNonCachingSpec > > {
        typedef ContainerBuffer< typename Source<TPipe>::Type const > Type;
    };


    //////////////////////////////////////////////////////////////////////////////
	// caching buffer handler
    struct _SourceCachingSpec;
	typedef Tag<_SourceCachingSpec> SourceCachingSpec;

	template < typename TPipe >
	struct BufferHandler< TPipe, SourceCachingSpec >
    {
        typedef typename Value<TPipe>::Type			TValue;
        typedef typename Size<TPipe>::Type			TSize;

        typedef SimpleBuffer<TValue>				TBuffer;
		typedef IPipeIterator<TPipe>				ISource;
		typedef typename Iterator<TBuffer>::Type	ITarget;

		TPipe		&pipe;
		size_t		bufferSize;
        TSize		rest;
        TBuffer		buffer;
		ISource		source;

		BufferHandler(TPipe &_pipe, size_t requestedSize):
			pipe(_pipe),
			bufferSize(requestedSize),
            rest(0) {}

        inline TBuffer& first() {
            rest = length(pipe);
			allocPage(buffer, Min(bufferSize, rest), *this);
			source = ISource(pipe);
			for(ITarget target = buffer.begin; target != buffer.end; ++target) {
				*target = *source;
				++source;
			}
            if (!(rest -= size(buffer))) source = ISource();
			return buffer;
        }

        inline TBuffer& next() {
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

	template < typename TPipe >
    struct Value< BufferHandler< TPipe, SourceCachingSpec > > {
        typedef SimpleBuffer< typename Value<TPipe>::Type > Type;
    };

    //////////////////////////////////////////////////////////////////////////////
	// buffer handler optimized for external string sources
    struct _ExtStringSourceCachingSpec;
	typedef Tag<_ExtStringSourceCachingSpec> ExtStringSourceCachingSpec;

	template < typename TSequence, typename TSpec >
	struct BufferHandler< Pipe<TSequence, TSpec>, ExtStringSourceCachingSpec >
    {
		typedef typename Value<TSequence>::Type		TValue;
        typedef typename Size<TSequence>::Type		TSize;
        typedef SimpleBuffer<TValue>				TBuffer;
        typedef Pipe<TSequence, TSpec>				TPipe;

		typedef typename Iterator<TSequence const>::Type	ISource;
		typedef typename Iterator<TBuffer>::Type			ITarget;

		TPipe		&pipe;
		unsigned	bufferSize;
        TSize		rest;
        TBuffer		buffer;
		ISource		source;

		BufferHandler(TPipe &_pipe, unsigned requestedSize):
			pipe(_pipe),
			bufferSize(requestedSize),
            rest(0) {}

        inline TBuffer& first() {
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

        inline TBuffer& next() {
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

	template < typename TSequence, typename TSpec >
    struct Value< BufferHandler< Pipe<TSequence, TSpec>, ExtStringSourceCachingSpec > > {
		typedef SimpleBuffer< typename Value<TSequence>::Type > Type;
    };

    //////////////////////////////////////////////////////////////////////////////
    // global functions

    // choose the most efficient buffer handler
    template < typename TInput, typename TSpec >
    struct BufReadHandler< Pipe<TInput, TSpec> > {
		typedef			
			typename IF< 
				AllowsFastRandomAccess<TInput>::VALUE,
				BufferHandler< Pipe<TInput, TSpec>, SourceNonCachingSpec>,
//				BufferHandler< Pipe<TInput, TSpec>, SourceCachingSpec>
				BufferHandler< Pipe<TInput, TSpec>, ExtStringSourceCachingSpec>
			>::Type Type;
    };


    template < typename TValue, typename TConfig, typename TSpec >
    struct BufReadHandler< Pipe< String<TValue, External<TConfig> >, TSpec > > {
		typedef BufferHandler< 
			Pipe< String<TValue, External<TConfig> > const, TSpec >, 
			ExtStringSourceCachingSpec 
		> Type;
    };

    template < typename TValue, typename TConfig, typename TSpec >
    struct BufReadHandler< Pipe< String<TValue, External<TConfig> > const, TSpec > > {
		typedef BufferHandler< 
			Pipe< String<TValue, External<TConfig> > const, TSpec >, 
			ExtStringSourceCachingSpec 
		> Type;
    };


	template < typename TInput, typename TSpec, typename TCommand >
	inline bool control(Pipe< TInput, Source<TSpec> > &me, TCommand const &command) {
        return true;
    }

	template < typename TInput, typename TSpec >
	inline bool control(Pipe< TInput, Source<TSpec> > &me, ControlBeginRead const &command) {
		me.cur = begin(me.in);
		return true;
	}
	
	template < typename TInput, typename TSpec >
	inline bool control(Pipe< TInput, Source<TSpec> > &me, ControlEndRead const &command) {
        me.cur = typename Iterator<TInput const>::Type();
		return true;
	}
	
	template < typename TInput, typename TSpec >
	inline bool control(Pipe< TInput, Source<TSpec> > &me, ControlEof const &command) {
		return me.cur == end(me.in);
	}

    template < typename TInput, typename TSpec >
    inline typename Size< Pipe< TInput, Source<TSpec> > >::Type
	length(Pipe< TInput, Source<TSpec> > &me) {
        return length(me.in);
    }

}

#endif
