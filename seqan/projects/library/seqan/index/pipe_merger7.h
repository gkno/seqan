/*
 *  merger7.h
 *  genindex
 *
 *  Created by David Weese on 17.07.05.
 *
 */

#ifndef SEQAN_HEADER_INDEX_MERGER7_H
#define SEQAN_HEADER_INDEX_MERGER7_H

namespace SEQAN_NAMESPACE_MAIN
{

//namespace SEQAN_NAMESPACE_PIPELINING
//{

    static char DC7_SHIFT[7][7] = {{0,6,5,6,3,3,5},
                                   {6,0,0,6,0,4,4},
                                   {5,0,0,1,0,1,5},
                                   {6,6,1,0,2,1,2},
                                   {3,0,0,2,0,3,2},
                                   {3,4,1,1,3,0,4},
                                   {5,4,5,2,2,4,0}};


    static char DC7_NINDX[7][7] = {{0,0,0,0,0,1,2},
                                   {0,0,0,0,1,1,2},
                                   {0,1,1,1,1,2,2},
                                   {0,0,1,1,1,1,2},
                                   {0,0,1,2,2,2,2},
                                   {0,0,0,1,2,2,2},
                                   {0,0,0,0,1,2,2}};


    template <typename TValue>
    struct SkewDCStream {
        TValue      i;
        unsigned    stream;
    };
    
    template <typename TValue>
	std::ostream& operator<<(std::ostream &out, const SkewDCStream<TValue> &s) {
		out << "< " << s.i.i1 << " , [ ";
		for(int i = 0; i < 3; ++i)
			out << s.i.i2[i] << " ";
		out << "] , [ ";
		for(int i = 0; i < 6; ++i)
			out << s.i.i3[i] << " ";
		out << "] , " << s.stream << " >";
		return out;
	}

    template <typename TValue>
    struct CompareSkewDCStream :
        public std::binary_function < SkewDCStream<TValue>,
                                      SkewDCStream<TValue>,
                                      bool >
    {
        inline bool operator()(const SkewDCStream<TValue> &a,
			                   const SkewDCStream<TValue> &b) const 
        {
            typedef typename TValue::T2::T SizeType;
            int shft = DC7_SHIFT[a.stream][b.stream];
//			printf("%d,%d___%d___%d:%d\n",a.stream,b.stream,shft,DC7_NINDX[a.stream][shft],DC7_NINDX[b.stream][shft]);
            for(int i = 0; i < shft; ++i) {
                if (a.i.i3[i] < b.i.i3[i]) return false;
                if (a.i.i3[i] > b.i.i3[i]) return true;
            }
            SizeType na = a.i.i2[DC7_NINDX[a.stream][shft]];
            SizeType nb = b.i.i2[DC7_NINDX[b.stream][shft]];
            if (na < nb) return false;
            if (na > nb) return true;
            SEQAN_ASSERT(false);  //names are unique -> we should never get here
            return (a.i.i1 < b.i.i1);
        }
    };

    template <typename T1, typename T2, typename T, const int _size>
    struct CompareSkewDCStream< Triple<T1,T2,Tuple<T,_size,Compressed>, Compressed> > :
        public std::binary_function < SkewDCStream<Triple<T1,T2,Tuple<T,_size,Compressed>,Compressed> >,
                                      SkewDCStream<Triple<T1,T2,Tuple<T,_size,Compressed>,Compressed> >,
                                      bool >
    {
        typedef Tuple<T,_size,Compressed> T3;
        inline bool operator()(const SkewDCStream<Triple<T1,T2,T3,Compressed> > &a,
			                   const SkewDCStream<Triple<T1,T2,T3,Compressed> > &b) const 
        {
            typedef typename T2::T SizeType;
            int shft = DC7_SHIFT[a.stream][b.stream];
            typename T3::CT mask = ~((1 << ((_size - shft) * T3::bitSize)) - 1);

            if ((a.i.i3.i & mask) < (b.i.i3.i & mask)) return false;
            if ((a.i.i3.i & mask) > (b.i.i3.i & mask)) return true;
            SizeType na = a.i.i2[DC7_NINDX[a.stream][shft]];
            SizeType nb = b.i.i2[DC7_NINDX[b.stream][shft]];
            if (na < nb) return false;
            if (na > nb) return true;
            SEQAN_ASSERT(false);  //names are unique -> we should never get here
            return (a.i.i1 < b.i.i1);
        }
    };

	struct Merger7;

    template < typename TInput0, typename TInput3, typename TInput5, typename TInput6, typename TInput124 >
    struct Value< Pipe< Bundle5< TInput0, TInput3, TInput5, TInput6, TInput124 >, Merger7 > > {
        typedef typename Size<TInput0>::Type Type;
    };


    //////////////////////////////////////////////////////////////////////////////
    // merger3 class
    template < typename TInput0, typename TInput3, typename TInput5, typename TInput6, typename TInput124 >
    struct Pipe< Bundle5< TInput0, TInput3, TInput5, TInput6, TInput124 >, Merger7 >
    {
        typedef typename Value<TInput0>::Type       InType0;
        typedef typename Value<TInput3>::Type       InType3;
        typedef typename Value<TInput5>::Type       InType5;
        typedef typename Value<TInput6>::Type       InType6;
        typedef typename Value<TInput124>::Type     InType124;
        typedef typename Size<Pipe>::Type           SizeType;

        typedef typename InType0::T3::T             Type;
        
        typedef SkewDCStream<InType0>               SkewDCStream;
        typedef CompareSkewDCStream<InType0>        CompareSkewDCStream;
        typedef ::std::priority_queue <
            SkewDCStream,
            ::std::vector<SkewDCStream>,
            CompareSkewDCStream >                   PQueue;

        Bundle5 <
            TInput0,
            TInput3,
            TInput5,
            TInput6,
            TInput124 > in;
        PQueue          queue;
        SizeType        N;
        
        Pipe(Bundle5< TInput0, TInput3, TInput5, TInput6, TInput124 > _in):
            in(_in) {}

        template <typename T1, typename T2, typename T3>
        inline static void _copy(SkewDCStream &dst, Triple<T1,T2,T3,Compressed> const &src) {
            memcpy(&dst.i.i2, &src.i2, sizeof(T2));
            memcpy(&dst.i.i3, &src.i3, sizeof(T3));
        }

        template <typename T1, typename T2, typename T, const int _size>
        inline static void _copy(SkewDCStream &dst, Triple<T1,T2,Tuple<T,_size,Compressed>,Compressed> const &src) {
            memcpy(&dst.i.i2, &src.i2, sizeof(T2));
            dst.i.i3.i = src.i3.i;
            dst.i.i3 <<= InType0::T3::size - _size;
        }

        template <typename TInput>
        inline void push(TInput &in, unsigned stream) {
			if (eof(in)) return;
            SkewDCStream s;
            s.i.i1 = N - (*in).i1;
            s.stream = stream;
            _copy(s, *in);
			queue.push(s);
            ++in;
        }

        inline void push124() {
			if (eof(in.in5)) return;
            SkewDCStream s;
			SizeType i1 = (*in.in5).i1;
            s.i.i1 = N - i1;
            s.stream = i1 % 7;
            _copy(s, *in.in5);
			queue.push(s);
            ++in.in5;
        }

        void fill() {
            push(in.in1, 0);
            push(in.in2, 3);
            push(in.in3, 5);
            push(in.in4, 6);
            push124();
        }

        inline typename Value<Pipe>::Type const operator*() {
            return queue.top().i.i1;
        }
        
        Pipe& operator++() {
            unsigned stream = queue.top().stream;
            queue.pop();
            switch (stream) {
                case 0:
                    push(in.in1, 0);
                    break;

				case 1:
				case 2:
				case 4:
					push124();
                    break;

                case 3:
                    push(in.in2, 3);
                    break;

                case 5:
                    push(in.in3, 5);
                    break;

                case 6:
                    push(in.in4, 6);
                    break;
            }
            return *this;
        }
    };
    

    //////////////////////////////////////////////////////////////////////////////
    // global pipe functions
    template < typename TInput >
	inline bool control(Pipe< TInput, Merger7 > &me, ControlBeginRead const &command) {
        if (!control(me.in, command)) return false;
        me.N = length(me);
        me.fill();
		return true;
	}
    
    template < typename TInput >
	inline bool control(Pipe< TInput, Merger7 > &me, ControlEof const &command) {
        return me.queue.size() == 0;
    }

    template < typename TInput >
    inline typename Size< Pipe< TInput, Merger7 > >::Type
    length(Pipe< TInput, Merger7 > const &me) {
        return length(me.in.in1) +
               length(me.in.in2) + 
               length(me.in.in3) +
               length(me.in.in4) +
               length(me.in.in5);
    }

//}

}

#endif
