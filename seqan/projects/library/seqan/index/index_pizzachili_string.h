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
  $Id$
 ==========================================================================*/

#ifndef SEQAN_HEADER_INDEX_PIZZACHILI_STRING_H
#define SEQAN_HEADER_INDEX_PIZZACHILI_STRING_H

#include <cstdlib>

#include "index_pizzachili_include.h"

namespace SEQAN_NAMESPACE_MAIN {

/**
.Spec.Pizza & Chili String:
..summary:String used by the Pizza & Chili indices.
..remarks:The string is lazy in the sense that it holds a reference to the
compressed index structure it is associated with. Only when the text is
actually read, the index is queried for the text. If only a substring is
needed, this string tries to query only a substring.
..cat:Strings
..general:Class.String
..signature:String<TValue, PizzaChili<TSpec> >
..param.TValue:The value type, that is the type of them items/characters stored in the string.
...remarks:This type must be a simple type and it must hold that $sizeof(TValue) == 1$.
..param.TSpec:unused.
...default:$void$
*/

template <typename TSpec = void>
struct PizzaChili;

template <typename TValue, typename TSpec>
class String<TValue, PizzaChili<TSpec> > {
public:
    typedef TValue* TIter;
    typedef TValue const* TConstIter;

    impl::pizzachili::index_t index_handle;
    bool owned;
    mutable TIter data_begin;
    mutable TIter data_end;

    String() : index_handle(0), owned(true), data_begin(0), data_end(0) { }

    template <typename TText>
    String(TText& other)
        : index_handle(0), owned(true), data_begin(0), data_end(0)
    {
        assign(*this, other);
    }

    String(impl::pizzachili::index_t index_handle)
        : index_handle(index_handle), owned(true), data_begin(0), data_end(0) { }

    String(String const& other)
        : index_handle(0), owned(true), data_begin(0), data_end(0)
    {
        assign(*this, other);
    }

    String& operator =(String const& other) {
        if (this != &other)
            assign(*this, other);
        return *this;
    }

    ~String() {
        clear(*this);
    }

    // Convenience ...

    template <typename TPos>
    inline TValue& operator [](TPos index) {
        return value(*this, index);
    }

    template <typename TPos>
    inline TValue operator [](TPos index) const {
        return value(*this, index);
    }
};

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSpec, typename TTag>
struct Iterator<String<TValue, PizzaChili<TSpec> >, TTag> {
    typedef typename String<TValue, PizzaChili<TSpec> >::TIter Type;
};

template <typename TValue, typename TSpec, typename TTag>
struct Iterator<String<TValue, PizzaChili<TSpec> > const, TTag> {
    typedef typename String<TValue, PizzaChili<TSpec> >::TConstIter Type;
};

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSpec>
struct Infix<String<TValue, PizzaChili<TSpec> > > {
    typedef String<TValue, PizzaChili<TSpec> > Type;
};

template <typename TValue, typename TSpec>
struct Prefix<String<TValue, PizzaChili<TSpec> > >
    : Infix<String<TValue, PizzaChili<TSpec> > > { };

template <typename TValue, typename TSpec>
struct Suffix<String<TValue, PizzaChili<TSpec> > >
    : Infix<String<TValue, PizzaChili<TSpec> > > { };

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSpec>
struct DefaultOverflowImplicit<String<TValue, PizzaChili<TSpec> > > {
    typedef Exact Type;
};

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSpec>
inline void
clear(String<TValue, PizzaChili<TSpec> >& me) {
SEQAN_CHECKPOINT
    me.index_handle = 0;
    if (me.owned) {
        _deallocateStorage(me, me.data_begin, me.data_end - me.data_begin);
        me.data_begin = 0;
        me.data_end = 0;
    }
}

//////////////////////////////////////////////////////////////////////////////

template<typename TValue, typename TSpec, typename TSource, typename TExpand>
inline void
assign(
    String<TValue, PizzaChili<TSpec> >& target,
    TSource const& source,
    Tag<TExpand> const
) {
SEQAN_CHECKPOINT
    typedef String<TValue, PizzaChili<TSpec> > TTarget;

    target.owned = true;
    _Assign_String<Tag<TExpand> const>::assign_(target, source);
}

template<typename TValue, typename TSpec, typename TSource, typename TExpand>
inline void
assign(
    String<TValue, PizzaChili<TSpec> >& target,
    TSource const* source,
    Tag<TExpand> const
) {
SEQAN_CHECKPOINT
    typedef String<TValue, PizzaChili<TSpec> > TTarget;

    target.owned = true;
    _Assign_String<Tag<TExpand> const>::assign_(target, source);
}

template<typename TValue, typename TSpec, typename TExpand>
inline void
assign(
    String<TValue, PizzaChili<TSpec> >& target,
    String<TValue, PizzaChili<TSpec> > const& source,
    Tag<TExpand> const
) {
SEQAN_CHECKPOINT
    typedef String<TValue, PizzaChili<TSpec> > TTarget;

    target.owned = true;

    if (source.index_handle != 0) {
        clear(target);
        target.index_handle = source.index_handle;
    }
    else
        _Assign_String<Tag<TExpand> const>::assign_(target, source);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSpec>
inline typename Value<String<TValue, PizzaChili<TSpec> > >::Type*
_allocateStorage(
    String<TValue, PizzaChili<TSpec> >& me,
    typename Size<String<TValue, PizzaChili<TSpec> > >::Type new_capacity
) {
SEQAN_CHECKPOINT
    TValue* old = me.data_begin;
    me.data_begin = static_cast<TValue*>(::std::malloc(new_capacity));
    return old == me.data_begin ? 0 : old;
}

template <typename TValue, typename TSpec>
inline typename Value<String<TValue, PizzaChili<TSpec> > >::Type*
_reallocateStorage(
    String<TValue, PizzaChili<TSpec> >& me,
    typename Size<String<TValue, PizzaChili<TSpec> > >::Type new_capacity,
    Exact
) {
SEQAN_CHECKPOINT
    if (new_capacity <= capacity(me))
        return 0;

    if (me.data_begin == 0)
        return _allocateStorage(me, new_capacity);

    me.data_begin = static_cast<TValue*>(::std::realloc(me.data_begin, new_capacity));
    // ::std::realloc does the cleanup itself.
    return 0;
}

template <typename TValue, typename TSpec>
inline void
_deallocateStorage(
   String<TValue, PizzaChili<TSpec> >& /*me*/,
   TValue* begin,
   typename Size<String<TValue, PizzaChili<TSpec> > >::Type /*count*/
) {
    if (begin != 0)
        ::std::free(begin);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSpec>
inline typename Size<String<TValue, PizzaChili<TSpec> > >::Type
length(String<TValue, PizzaChili<TSpec> > const& me) {
SEQAN_CHECKPOINT
    if (me.data_begin != 0)
        return me.data_end - me.data_begin;
    else if (me.index_handle != 0) {
        impl::pizzachili::ulong len;
        impl::pizzachili::error_t e =
            impl::pizzachili::length(me.index_handle, &len);
        if (e != 0) {
            struct { } ex;
            throw ex;
        }

        return len;
    }
    else
        return 0;
}

template <typename TValue, typename TSpec>
inline void
_setLength(
    String<TValue, PizzaChili<TSpec> >& me,
    size_t new_length
) {
SEQAN_CHECKPOINT
    me.data_end = me.data_begin + new_length;
}

//////////////////////////////////////////////////////////////////////////////

namespace impl {
    template <typename TValue, typename TSpec>
    inline void
    queryText(String<TValue, PizzaChili<TSpec> > const& me) {
        if (me.data_begin != 0)
            return;
        if (me.index_handle == 0) {
            me.data_begin = static_cast<TValue*>(::std::malloc(1));
            me.data_begin[0] = '\0';
            me.data_end = me.data_begin;
        }
        else {
            pizzachili::uchar* snippet;
            pizzachili::ulong len;
            pizzachili::error_t e =
                pizzachili::extract(
                    me.index_handle,
                    0,
                    length(me) - 1,
                    &snippet,
                    &len
                );

            if (e != 0) {
    #if DEBUG || TRACE
                SEQAN_TREPORT(pizzachili::error_index(e));
    #endif
                struct { } ex;
                throw ex;
            }

            me.data_begin = reinterpret_cast<TValue*>(snippet);
            me.data_end = me.data_begin + len;
        }
    }
} // namespace impl

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSpec, typename TTag>
inline typename Iterator<String<TValue, PizzaChili<TSpec> >, Tag<TSpec> const>::Type 
begin(
    String<TValue, PizzaChili<TSpec> >& me,
    Tag<TTag> const
) {
SEQAN_CHECKPOINT
    impl::queryText(me);
    return me.data_begin;
}

template <typename TValue, typename TSpec, typename TTag>
inline typename Iterator<String<TValue, PizzaChili<TSpec> > const, Tag<TSpec> const>::Type 
begin(
    String<TValue, PizzaChili<TSpec> > const& me,
    Tag<TTag> const
) {
SEQAN_CHECKPOINT
    impl::queryText(me);
    return me.data_begin;
}

template <typename TValue, typename TSpec, typename TTag>
inline typename Iterator<String<TValue, PizzaChili<TSpec> >, Tag<TSpec> const>::Type 
end(
    String<TValue, PizzaChili<TSpec> >& me,
    Tag<TTag> const
) {
SEQAN_CHECKPOINT
    impl::queryText(me);
    return me.data_end;
}

template <typename TValue, typename TSpec, typename TTag>
inline typename Iterator<String<TValue, PizzaChili<TSpec> > const, Tag<TSpec> const>::Type 
end(
    String<TValue, PizzaChili<TSpec> > const& me,
    Tag<TTag> const
) {
SEQAN_CHECKPOINT
    impl::queryText(me);
    return me.data_end;
}

//////////////////////////////////////////////////////////////////////////////

namespace impl {
    template <typename TValue, typename TSpec, typename TPos>
    struct substringHelperPizzaChili {
        typedef typename Infix<String<TValue, PizzaChili<TSpec> > >::Type TResult;

        static inline TResult
        infix(String<TValue, PizzaChili<TSpec> > const& me, TPos begin, TPos end) {
SEQAN_CHECKPOINT
            TResult ret;

            if (me.data_begin != 0) {
                ret.owned = false;
                ret.data_begin = me.data_begin + begin;
                ret.data_end = me.data_begin + end;
            }
            else if (me.index_handle != 0) {
                pizzachili::uchar* snippet;
                pizzachili::ulong len;
                pizzachili::error_t e =
                    pizzachili::extract(
                        me.index_handle,
                        begin,
                        end - 1,
                        &snippet,
                        &len
                    );

                if (e != 0) {
        #if DEBUG || TRACE
                    SEQAN_TREPORT(pizzachili::error_index(e));
        #endif
                    struct { } ex;
                    throw ex;
                }

                ret.data_begin = reinterpret_cast<TValue*>(snippet);
                ret.data_end = ret.data_begin + len;
            }

            return ret;
        }

        static inline TResult
        prefix(String<TValue, PizzaChili<TSpec> > const& me, TPos end) {
            return infix(me, 0, end);
        }

        static inline TResult
        suffix(String<TValue, PizzaChili<TSpec> > const& me, TPos begin) {
            return infix(me, begin, length(me));
        }
    };

    template <typename TValue, typename TSpec>
    struct substringHelperPizzaChili<
        TValue,
        TSpec,
        typename Iterator<
            String<TValue, PizzaChili<TSpec> >,
            typename DefaultIteratorSpec<String<TValue, PizzaChili<TSpec> > >::Type
        >::Type
    > {
        typedef
            typename Iterator<
                String<TValue, PizzaChili<TSpec> >,
                typename DefaultIteratorSpec<String<TValue, PizzaChili<TSpec> > >::Type
            >::Type TPos;
        typedef typename Infix<String<TValue, PizzaChili<TSpec> > >::Type TResult;

        static inline TResult
        infix(String<TValue, PizzaChili<TSpec> > const& /*me*/, TPos begin, TPos end) {
SEQAN_CHECKPOINT
            // Iterators were used, therefore it's safe to assume that the
            // string is already in memory.
            TResult ret;
            ret.owned = false;
            ret.data_begin = begin;
            ret.data_end = end;
            return ret;
        }

        static inline TResult
        prefix(String<TValue, PizzaChili<TSpec> > const& me, TPos end) {
            return infix(me, begin(me), end);
        }

        static inline TResult
        prefix(String<TValue, PizzaChili<TSpec> >& me, TPos end) {
            return infix(me, begin(me), end);
        }

        static inline TResult
        suffix(String<TValue, PizzaChili<TSpec> > const& me, TPos begin) {
            return infix(me, begin, end(me));
        }

        static inline TResult
        suffix(String<TValue, PizzaChili<TSpec> >& me, TPos begin) {
            return infix(me, begin, end(me));
        }
    };

} // namespace impl

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSpec, typename TPosBegin, typename TPosEnd>
inline typename Infix<String<TValue, PizzaChili<TSpec> > >::Type
infix(
    String<TValue, PizzaChili<TSpec> > const& me,
    TPosBegin begin,
    TPosEnd end
) {
SEQAN_CHECKPOINT
    return impl::substringHelperPizzaChili<TValue, TSpec, TPosBegin>::infix(me, begin, end);
}

template <typename TValue, typename TSpec, typename TPosBegin, typename TPosEnd>
inline typename Infix<String<TValue, PizzaChili<TSpec> > >::Type
infix(
    String<TValue, PizzaChili<TSpec> >& me,
    TPosBegin begin,
    TPosEnd end
) {
SEQAN_CHECKPOINT
    return impl::substringHelperPizzaChili<TValue, TSpec, TPosBegin>::infix(me, begin, end);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSpec, typename TPos>
inline typename Prefix<String<TValue, PizzaChili<TSpec> > >::Type
prefix(
    String<TValue, PizzaChili<TSpec> > const& me,
    TPos end
) {
SEQAN_CHECKPOINT
    return impl::substringHelperPizzaChili<TValue, TSpec, TPos>::prefix(me, end);
}

template <typename TValue, typename TSpec, typename TPos>
inline typename Prefix<String<TValue, PizzaChili<TSpec> > >::Type
prefix(
    String<TValue, PizzaChili<TSpec> >& me,
    TPos end
) {
SEQAN_CHECKPOINT
    return impl::substringHelperPizzaChili<TValue, TSpec, TPos>::prefix(me, end);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSpec, typename TPos>
inline typename Suffix<String<TValue, PizzaChili<TSpec> > >::Type
suffix(
    String<TValue, PizzaChili<TSpec> > const& me,
    TPos begin
) {
SEQAN_CHECKPOINT
    return impl::substringHelperPizzaChili<TValue, TSpec, TPos>::suffix(me, begin);
}

template <typename TValue, typename TSpec, typename TPos>
inline typename Suffix<String<TValue, PizzaChili<TSpec> > >::Type
suffix(
    String<TValue, PizzaChili<TSpec> >& me,
    TPos begin
) {
SEQAN_CHECKPOINT
    return impl::substringHelperPizzaChili<TValue, TSpec, TPos>::suffix(me, begin);
}

//////////////////////////////////////////////////////////////////////////////

} // namespace SEQAN_NAMESPACE_MAIN

#endif // SEQAN_HEADER_INDEX_PIZZACHILI_STRING_H
