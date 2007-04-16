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
	// simple class for fast random accessable containers buffers
    template < typename TIterator >
	struct IteratorBuffer 
	{
		TIterator           begin;
        TIterator           end;

        IteratorBuffer():
            begin(TIterator()),
            end(TIterator()) {}

		IteratorBuffer(TIterator _begin, TIterator _end):
            begin(_begin),
            end(_end) {}

		template <typename TSize>
        IteratorBuffer(TIterator _begin, TSize _size):
            begin(_begin),
            end(_begin + _size) {}

		template <typename TSize>
		inline typename Reference<TIterator>::Type 
		operator[] (TSize i) {
			return *(begin + i);
		}
    };

    template < typename TIterator >
    struct Iterator< IteratorBuffer<TIterator> > {
        typedef TIterator Type;
    };

    template < typename TIterator  >
	struct Value< IteratorBuffer<TIterator> >:
		Value<TIterator> {};

    template < typename TIterator >
	struct Size< IteratorBuffer<TIterator> >:
		Size<TIterator> {};


	template < typename TIterator >
    inline typename Size< IteratorBuffer<TIterator> >::Type
    length(IteratorBuffer<TIterator> const &me) {
        return me.end - me.begin;
    }


	//////////////////////////////////////////////////////////////////////////////

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

    //////////////////////////////////////////////////////////////////////////////
	// handler that manages a simple memory buffer
	template < typename TInput, typename TSpec >
	struct BufferHandler< Pipe< TInput, Source<TSpec> > >
    {
        typedef typename Iterator<TInput const>::Type	Iterator;
        typedef IteratorBuffer<Iterator>				Buffer;
        typedef Pipe< TInput, Source<TSpec> >			Pipe;

		Pipe	&pipe;

		BufferHandler(Pipe &_pipe):
			pipe(_pipe) {}

		template <typename TSize>
		BufferHandler(Pipe &_pipe, TSize):
			pipe(_pipe) {}

        inline Buffer first() {
            return Buffer(begin(pipe.in), length(pipe.in));
        }

        inline Buffer next() {
            return Buffer();
        }

        inline void process() {}
        inline void end() {}
        inline void cancel() {}
    };

	template < typename TInput, typename TSpec >
    struct Value< BufferHandler< Pipe< TInput, Source<TSpec> > > > {
        typedef IteratorBuffer<	typename Iterator< TInput const >::Type	> Type;
    };


    //////////////////////////////////////////////////////////////////////////////
    // global functions

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
