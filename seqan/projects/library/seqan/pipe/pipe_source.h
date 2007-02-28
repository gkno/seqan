/*
 *  pipe_source.h
 *  genindex
 *
 *  Created by David Weese on 17.07.05.
 *
 */

#ifndef SEQAN_HEADER_PIPE_SOURCE_H
#define SEQAN_HEADER_PIPE_SOURCE_H

namespace SEQAN_NAMESPACE_MAIN
{

	struct IteratorSpec;
	struct ContainerSpec;

	template < typename TSpec = ContainerSpec >
	struct Source;

/*
	template < typename TInput >
    struct Value< Pipe< TInput, Source<IteratorSpec> > > {
		typedef typename ::std::iterator_traits<TInput>::value_type Type;
	};

	template < typename TInput >
    struct Size< Pipe< TInput, Source<IteratorSpec> > > {
		typedef typename ::std::iterator_traits<TInput>::difference_type Type;
	};
*/
	template < typename TValue >
    struct Value< Pipe< TValue*, Source<IteratorSpec> > > {
		typedef TValue Type;
	};

	template < typename TValue >
    struct Size< Pipe< TValue*, Source<IteratorSpec> > > {
		typedef ptrdiff_t Type;
	};
/*
	template < typename TInput >
    struct Value< Pipe< TInput, Source<ContainerSpec> > > {
		typedef typename ::std::iterator_traits<typename TInput::const_iterator>::value_type Type;
	};

	template < typename TInput >
    struct Size< Pipe< TInput, Source<ContainerSpec> > > {
		typedef typename ::std::iterator_traits<typename TInput::const_iterator>::difference_type Type;
	};
*/

/**
.Spec.Source:
..cat:Pipelining
..general:Class.Pipe
..summary:Pipelining adaptor for arbitrary containers or iterators.
..signature:Pipe<TInput, Source<ContainerSpec> >
..signature:Pipe<TInput, Source<IteratorSpec> >
..param.TInput:The type of container or iterator this module reads from.
*/

    //////////////////////////////////////////////////////////////////////////////
    // source class
	template < typename TInput >
	struct Pipe< TInput, Source<IteratorSpec> >
    {
        TInput	begin, end;
        TInput	cur;
        
		Pipe(TInput _begin, TInput _end):
            begin(_begin),
            end(_end),
			cur(_begin) {}
        
        template < typename TSize >
		Pipe(TInput _begin, TSize _size):
//		Pipe(TInput _begin, typename Size<Pipe>::Type _size):
            begin(_begin),
            end(_begin + _size),
			cur(_begin) {}
        
        inline typename Value<Pipe>::Type const & operator*() const {
            return *cur;
        }
        
        inline Pipe& operator++() {
            ++cur;
            return *this;
        }
    };


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
	// handler that manages a simple memory buffer
	template < typename TInput >
	struct BufferHandler< Pipe< TInput, Source<IteratorSpec> >, void >
    {
        typedef typename Value< Pipe< TInput, Source<IteratorSpec> > >::Type Type;
        typedef typename Size< Pipe< TInput, Source<IteratorSpec> > >::Type  SizeType;

        typedef IteratorBuffer< TInput, Type, SizeType >    Buffer;
        typedef Pipe< TInput, Source<IteratorSpec> >        Pipe;

		Pipe	&pipe;
        Buffer	buffer;

		BufferHandler(Pipe &_pipe):
			pipe(_pipe) {}

		template <typename TSize>
		BufferHandler(Pipe &_pipe, TSize):
			pipe(_pipe) {}

        inline Buffer& begin() {            
            return buffer = Buffer(pipe.begin, pipe.end);
        }

        inline Buffer& next() {
            return buffer = Buffer();
        }

        inline void process() {}
        inline void end() {}
        inline void cancel() {}
    };

	template < typename TInput >
    struct Value< BufferHandler< Pipe< TInput, Source<IteratorSpec> > > > {
        IteratorBuffer< TInput,
                        typename Value< Pipe< TInput, Source<IteratorSpec> > >::Type,
                        typename Size< Pipe< TInput, Source<IteratorSpec> > >::Type > Type;
    };


	template < typename TInput >
	struct Pipe< TInput, Source<ContainerSpec> >
    {
		TInput const &in;
		typename Iterator<TInput const>::Type cur;

		Pipe(TInput &_cont):
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
	template < typename TInput >
	struct BufferHandler< Pipe< TInput, Source<ContainerSpec> > >
    {
        typedef typename Value< Pipe< TInput, Source<ContainerSpec> > >::Type   Type;
        typedef typename Size< Pipe< TInput, Source<ContainerSpec> > >::Type    SizeType;
        typedef typename Iterator<TInput const>::Type							Iterator;

        typedef IteratorBuffer< Iterator, Type, SizeType >	Buffer;
        typedef Pipe< TInput, Source<ContainerSpec> >		Pipe;

		Pipe	&pipe;
        Buffer	buffer;

		BufferHandler(Pipe &_pipe):
			pipe(_pipe) {}

		template <typename TSize>
		BufferHandler(Pipe &_pipe, TSize):
			pipe(_pipe) {}

        inline Buffer& first() {
            return buffer = Buffer(begin(pipe.in), length(pipe.in));
        }

        inline Buffer& next() {
            return buffer = Buffer();
        }

        inline void process() {}
        inline void end() {}
        inline void cancel() {}
    };

	template < typename TInput >
    struct Value< BufferHandler< Pipe< TInput, Source<ContainerSpec> > > > {
        IteratorBuffer< typename Iterator< TInput const >::Type,
                        typename Value< Pipe< TInput, Source<ContainerSpec> > >::Type,
                        typename Size< Pipe< TInput, Source<ContainerSpec> > >::Type > Type;
    };


    //////////////////////////////////////////////////////////////////////////////
    // global functions
    template < typename TInput, typename TCommand >
	inline bool control(Pipe< TInput, Source<IteratorSpec> > &me, TCommand const &command) {
        return true;
    }

	template < typename TInput >
	inline bool control(Pipe< TInput, Source<IteratorSpec> > &me, ControlBeginRead const &command) {
		me.cur = me.begin;
		return true;
	}
	
	template < typename TInput >
	inline bool control(Pipe< TInput, Source<IteratorSpec> > &me, ControlEndRead const &command) {
		me.cur = TInput();
		return true;
	}
	
	template < typename TInput >
	inline bool control(Pipe< TInput, Source<IteratorSpec> > &me, ControlEof const &command) {
		return me.cur == me.end;
	}

    template < typename TInput >
    inline typename Size< Pipe< TInput, Source<IteratorSpec> > >::Type
	length(Pipe< TInput, Source<IteratorSpec> > &me) {
        return me.end - me.begin;
    }


    template < typename TInput, typename TCommand >
	inline bool control(Pipe< TInput, Source<ContainerSpec> > &me, TCommand const &command) {
        return true;
    }

	template < typename TInput >
	inline bool control(Pipe< TInput, Source<ContainerSpec> > &me, ControlBeginRead const &command) {
		me.cur = begin(me.in);
		return true;
	}
	
	template < typename TInput >
	inline bool control(Pipe< TInput, Source<ContainerSpec> > &me, ControlEndRead const &command) {
        me.cur = typename Iterator<TInput const>::Type();
		return true;
	}
	
	template < typename TInput >
	inline bool control(Pipe< TInput, Source<ContainerSpec> > &me, ControlEof const &command) {
		return me.cur == end(me.in);
	}

    template < typename TInput >
    inline typename Size< Pipe< TInput, Source<ContainerSpec> > >::Type
	length(Pipe< TInput, Source<ContainerSpec> > &me) {
//        return size(me.in);
        return length(me.in);
    }

}

#endif
