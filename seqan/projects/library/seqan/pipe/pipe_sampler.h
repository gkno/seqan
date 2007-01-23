/*
 *  pipe_sampler.h
 *  genindex
 *
 *  Created by David Weese on 17.07.05.
 *
 */

#ifndef SEQAN_HEADER_PIPE_SAMPLER_H
#define SEQAN_HEADER_PIPE_SAMPLER_H

namespace SEQAN_NAMESPACE_MAIN
{

//namespace SEQAN_NAMESPACE_PIPELINING
//{

    template < const unsigned m, typename compress = void >
	struct Sampler;

    template < typename TInput, const unsigned m, typename compress >
    struct Value< Pipe< TInput, Sampler<m, compress> > > {
        typedef Tuple<typename Value<TInput>::Type, m, compress>		mTuple;
        typedef Pair<typename Size<TInput>::Type, mTuple, Compressed>	Type;
    };

	template <int I, typename T = void>
	struct _SkewDC;


/**
.Spec.Sampler:
..cat:Pipelining
..general:Class.Pipe
..summary:Outputs m-tuples beginning at a position of difference cover DC.
..signature:Pipe<TInput, Sampler<m, DC[, compress]> >
..param.TInput:The type of the pipeline module this module reads from.
..param.m:The tuple size.
..param.DC:A set of non-negative integers less than $m$.
..param.DC:$DC[0]$ contains the size of the set and $DC[1..DC[0]]$ contains the distinct and ordered elements.
..param.compress:Enable/Disable compression.
..param.compress:If $void$, no compression is used.
..param.compress:If $Compressed$, bit-compressed @Class.Tuple@s are used.
...default:void.
..example:The set ${1,2,4}$ is represented by $int DC[] = { 3, 1, 2, 4 }$.
..remarks:The output type is a @Class.Pair@ of size type and @Class.Tuple@ of input elements and length m (i.e. $Pair<Size<TInput>::Type, Tuple<Value<TInput>::Type, m, compress> >$).
..remarks:The first output field contains the beginning position of the m-tuple in the second field.
The m-tuples are substrings of the input stream beginning at positions $i$, with $i mod m$ is element of the set DC.
*/

    //////////////////////////////////////////////////////////////////////////////
    // sampler class
    template < typename TInput, const unsigned m, typename compress >
    struct Pipe< TInput, Sampler<m, compress> >
    {
        typedef typename Value<Pipe>::Type  OutType;
        typedef typename Size<Pipe>::Type   SizeType;

		TInput		&in;
        bool        filter[m];
        SizeType    idx, _size, _rest;
        unsigned    idxMod;
        OutType     tmp1, tmp2;
        OutType     *outRef, *tmpRef;
        bool        last;
        
        Pipe(TInput& _in):
            in(_in),
            outRef(&tmp1),
            tmpRef(&tmp2) {}
        
        inline void prepare() {
            memset<sizeof(filter), 0>(filter);
			for(unsigned i = 1; i <= _SkewDC<m>::VALUE[0]; i++)
                filter[_SkewDC<m>::VALUE[i]] = true;

            idx = length(in);
            idxMod = idx % m;

            while (!filter[idxMod] && !eof(in)) {
                ++in;
                if (idxMod == 0) idxMod = m;
                --idxMod; --idx;
            }
            _rest = length(*this);
            fill(m);
            swap();
        }
        
        inline void fill(int f) {
            int i;
            for(i = 0; i < f && !eof(in); ++i, ++in)
                tmpRef->i2.i[i] = *in;
            last = eof(in);
            for(; i < f; ++i)
                tmpRef->i2.i[i] = 0;
            tmpRef->i1 = idx;
        }
        
        inline void rotate(int r) {
            for(unsigned i = 0; i < m; ++i, ++r) {
                if (r == m) r = 0;
                tmpRef->i2.i[i] = outRef->i2.i[r];
            }
        }
        
        inline void swap() {
            OutType *newOutRef = tmpRef;
            tmpRef = outRef;
            outRef = newOutRef;
        }
        
        inline OutType const& operator*() {
            return *outRef;
        }
        
        Pipe& operator++() {
            unsigned skipped = 0;
            if (!last)
                do {
                    outRef->i2.i[skipped++] = *in;
                    ++in;
                    if (idxMod == 0) idxMod = m;
                    --idxMod; --idx;
                    if (eof(in)) {
                        last = true;
                        while (!filter[idxMod]) {
                            outRef->i2.i[skipped++] = 0;
                            if (idxMod == 0) idxMod = m;
                            --idxMod; --idx;
                        };
                        break;
                    }
                } while (!filter[idxMod]);
            else
                do {
                    outRef->i2.i[skipped++] = 0;
                    if (idxMod == 0) idxMod = m;
                    --idxMod; --idx;
                } while (!filter[idxMod]);
            rotate(skipped);
            --_rest;
            tmpRef->i1 = idx;
            swap();
            return *this;
        }        
    };


    //////////////////////////////////////////////////////////////////////////////
    // sampler class (uses bit compression)
    template < typename TInput, const unsigned m >
    struct Pipe< TInput, Sampler<m, Compressed> >
    {
        typedef typename Value<Pipe>::Type  OutType;
        typedef typename Size<Pipe>::Type   SizeType;
        typedef typename OutType::T2        TTuple;

		TInput		&in;
        bool        filter[m];
        SizeType    _size, _rest;
        unsigned    idxMod;
        OutType     tmp;
        bool        last;
        
        Pipe(TInput& _in):
            in(_in) {}
        
        inline void prepare() {
            memset<sizeof(filter), 0>(filter);
            for(unsigned i = 1; i <= _SkewDC<m>::VALUE[0]; i++)
                filter[_SkewDC<m>::VALUE[i]] = true;

            tmp.i1 = length(in);
            idxMod = tmp.i1 % m;

            while (!filter[idxMod] && !eof(in)) {
                ++in;
                if (idxMod == 0) idxMod = m;
                --idxMod; --tmp.i1;
            }
            _rest = length(*this);
            fill(m);
        }
        
        inline void fill(int f) {
            int i;
            tmp.i2.i = 0;
            for(i = 0; i < f && !eof(in); ++i, ++in) {
                tmp.i2 <<= 1;
                tmp.i2 |= *in;
            }
            last = eof(in);
            tmp.i2.i <<= (f - i);
        }
        
        inline OutType const& operator*() {
            return tmp;
        }
        
        Pipe& operator++() {
            if (!last)
                do {
                    tmp.i2 <<= 1;
                    tmp.i2 |= *in;
                    ++in;
                    if (idxMod == 0) idxMod = m;
                    --idxMod; --tmp.i1;
                    if (eof(in)) {
                        last = true;
                        while (!filter[idxMod]) {
                            tmp.i2 <<= 1;
                            if (idxMod == 0) idxMod = m;
                            --idxMod; --tmp.i1;
                        };
                        break;
                    }
                } while (!filter[idxMod]);
            else
                do {
                    tmp.i2 <<= 1;
                    if (idxMod == 0) idxMod = m;
                    --idxMod; --tmp.i1;
                } while (!filter[idxMod]);
            --_rest;
            return *this;
        }        
    };


    //////////////////////////////////////////////////////////////////////////////
    // global pipe functions
    template < typename TInput, const unsigned int m, typename compress >
	inline bool control(Pipe< TInput, Sampler<m, compress> > &me, ControlBeginRead const &command) {
        if (!control(me.in, command)) return false;
        me.prepare();
		return true;
    }

    template < typename TInput, const unsigned int m, typename compress >
	inline bool control(Pipe< TInput, Sampler<m, compress> > &me, ControlEof const &command) {
		return me._rest == 0;
    }

    template < typename TInput, const unsigned int m, typename compress >
    inline typename Size< Pipe< TInput, Sampler<m, compress> > >::Type
	length(Pipe< TInput, Sampler<m, compress> > &me) {
        typename Size< Pipe< TInput, Sampler<m> > >::Type _size = 0, n = length(me.in);
        for(unsigned i = 1; i <= _SkewDC<m>::VALUE[0]; i++)
            if (_SkewDC<m>::VALUE[i])
                _size += (n + m - _SkewDC<m>::VALUE[i]) / m;
            else
                _size += n / m;
        return _size;
    }
//}

}

#endif
