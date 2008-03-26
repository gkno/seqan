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

#ifndef SEQAN_HEADER_INDEX_PIZZACHILI_H
#define SEQAN_HEADER_INDEX_PIZZACHILI_H

#include "index_pizzachili_include.h"
#include "index_pizzachili_string.h"

namespace SEQAN_NAMESPACE_MAIN {

/**
.Tag.Pizza & Chili Index Fibres
..summary:Tag to select a specific fibre (e.g. table, object, ...) of a @Spec.Pizza & Chili Index@ index.
..cat:Index
..tag.PizzaChili_Text:The original text the index is based on.
..tag.PizzaChili_Compressed:The compressed suffix array.
...remarks:Pizza & Chili indices are compressed indices. Hence, this fibre is used for searching in the index.
..see:Metafunction.Fibre
..see:Function.getFibre
..see:Spec.Pizza & Chili Index
*/

struct _Fibre_PizzaChili_Compressed;

typedef Tag<_Fibre_Text> const Fibre_PizzaChili_Text;
typedef Tag<_Fibre_PizzaChili_Compressed> const Fibre_PizzaChili_Compressed;

typedef Fibre_PizzaChili_Text PizzaChili_Text;
typedef Fibre_PizzaChili_Compressed PizzaChili_Compressed;

/**
.Spec.Pizza & Chili Index:
..summary:An adapter for the Pizza & Chili index API.
..remarks:
..cat:Index
..general:Class.Index
..signature:Index<TText, PizzaChili<TSpec> >
..param.TText:The text type.
...type:Class.String
..param.TSpec:Algorithm-specific parameters.
...default:$void$
..demo:index_pizzachili.h
..see:Spec.Pizza & Chili String
..see:Tag.Pizza & Chili Index Fibres
..see:Tag.Index Find Algorithm
*/

// Already declared in included file "index_pizzachili_string.h".
//template <typename TSpec = void>
//struct PizzaChili;

template <typename TText, typename TSpec>
class Index<TText, PizzaChili<TSpec> > {
public:
    typedef typename Value<TText>::Type TValue;

    impl::pizzachili::index_t index_handle;
    Holder<String<TValue, PizzaChili<> > > text;

    Index() : index_handle(0), text() { }

    template <typename TOtherText>
    Index(TOtherText& text) : index_handle(0), text() {
        setIndexText(text);
    }

    ~Index() {
        clear(*this);
    }

private:
    Index(Index const& other);
    Index operator =(Index const& other);
};

//////////////////////////////////////////////////////////////////////////////

template <typename TText, typename TSpec>
struct Fibre<Index<TText, PizzaChili<TSpec> >, PizzaChili_Text> {
    typedef String<typename Value<TText>::Type, PizzaChili<> >& Type;
};

template <typename TText, typename TSpec>
struct Fibre<Index<TText, PizzaChili<TSpec> > const, PizzaChili_Text> {
    typedef String<typename Value<TText>::Type, PizzaChili<> > const& Type;
};

//////////////////////////////////////////////////////////////////////////////

namespace impl {
    template <typename TText, typename TSpec>
    inline void
    clearIndex(Index<TText, PizzaChili<TSpec> >& me) {
        if (me.index_handle != 0) {
            impl::pizzachili::error_t e =
                impl::pizzachili::free_index(me.index_handle);

            if (e != 0) {
                struct {} ex;
                throw ex;
            }

            me.index_handle = 0;
        }
    }
} // namespace impl

template <typename TText, typename TSpec>
inline void
clear(Index<TText, PizzaChili<TSpec> >& me) {
SEQAN_CHECKPOINT
    impl::clearIndex(me);
    clear(me.text);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TText, typename TSpec>
inline typename Fibre<Index<TText, PizzaChili<TSpec> >, PizzaChili_Text>::Type
indexText(Index<TText, PizzaChili<TSpec> >& me) {
SEQAN_CHECKPOINT
    return getFibre(me, PizzaChili_Text());
}

//////////////////////////////////////////////////////////////////////////////

template <typename TText, typename TSpec>
inline typename Fibre<Index<TText, PizzaChili<TSpec> >, PizzaChili_Text>::Type
getFibre(Index<TText, PizzaChili<TSpec> >& me, PizzaChili_Text const) {
SEQAN_CHECKPOINT
    return value(me.text);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TText, typename TSpec>
inline bool
indexSupplied(Index<TText, PizzaChili<TSpec> >& me, PizzaChili_Compressed const) {
SEQAN_CHECKPOINT
    return me.index_handle != 0;
}

template <typename TText, typename TSpec>
inline bool
indexSupplied(Index<TText, PizzaChili<TSpec> >& me, PizzaChili_Text const) {
SEQAN_CHECKPOINT
    return length(value(me.text)) > 0;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TText, typename TSpec>
inline bool
indexSolveDependencies(Index<TText, PizzaChili<TSpec> >& me, PizzaChili_Compressed const) {
SEQAN_CHECKPOINT
    return indexSupplied(me, PizzaChili_Text());
}

//////////////////////////////////////////////////////////////////////////////

namespace impl {
    template <typename TText, typename TSpec>
    inline bool
    createPizzaChiliIndex(
        Index<TText, PizzaChili<TSpec> >& me,
        pizzachili::uchar* textstart,
        pizzachili::ulong textlength
    ) {
        // Read-only access, therefore safe cast.
        char* options = const_cast<char*>("");
        pizzachili::error_t e =
            pizzachili::build_index(textstart, textlength, options, &me.index_handle);

        if (e != 0) {
    #if DEBUG || TRACE
            SEQAN_TREPORT(pizzachili::error_index(e));
    #endif
            me.index_handle = 0;
            return false;
        }

        value(me.text) = String<typename Value<TText>::Type, PizzaChili<> >(me.index_handle);

        return true;
    }
} // namespace impl

template <typename TText, typename TSpec>
inline bool
indexCreate(Index<TText, PizzaChili<TSpec> >& me, PizzaChili_Compressed const) {
SEQAN_CHECKPOINT
    typedef
        typename _RemoveConst<
            typename Index<TText, PizzaChili<TSpec> >::TValue
        >::Type alph_t;

    SEQAN_ASSERT(
        BitsPerValue<alph_t>::VALUE == 8 &&
        IsSimple<alph_t>::VALUE
    );

    impl::clearIndex(me);

    impl::pizzachili::uchar* textstart =
        reinterpret_cast<impl::pizzachili::uchar*>(
            const_cast<alph_t*>(indexText(me).data_begin)
        );
    impl::pizzachili::ulong textlength = length(indexText(me));
    return impl::createPizzaChiliIndex(me, textstart, textlength);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TText, typename TSpec, typename TOtherText>
inline void
setIndexText(Index<TText, PizzaChili<TSpec> >& me, TOtherText& text) {
SEQAN_CHECKPOINT
    clear(me);
    typedef
        typename _RemoveConst<
            typename Value<TOtherText>::Type
        >::Type alph_t;

    SEQAN_ASSERT(
        IsContiguous<TOtherText>::VALUE &&
        BitsPerValue<alph_t>::VALUE == 8 &&
        IsSimple<alph_t>::VALUE
    );

    String<alph_t, CStyle> cstr = text;
    impl::pizzachili::uchar* textstart =
        reinterpret_cast<impl::pizzachili::uchar*>(
            const_cast<alph_t*>(static_cast<alph_t const*>(cstr))
        );
    impl::pizzachili::ulong textlength = length(text);
    impl::createPizzaChiliIndex(me, textstart, textlength);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TText, typename TSpec>
inline bool open(
    Index<TText, PizzaChili<TSpec> >& me,
    char const* filename
) {
SEQAN_CHECKPOINT
    clear(me);
    impl::pizzachili::error_t e =
        impl::pizzachili::load_index(const_cast<char*>(filename), &me.index_handle);
    return e == 0;
}

template <typename TText, typename TSpec>
inline bool save(
    Index<TText, PizzaChili<TSpec> >& me,
    char const* filename
) {
SEQAN_CHECKPOINT
    impl::pizzachili::error_t e =
        impl::pizzachili::save_index(me.index_handle, const_cast<char*>(filename));
    return e == 0;
}

} // namespace SEQAN_NAMESPACE_MAIN

#endif // SEQAN_HEADER_INDEX_PIZZACHILI_H
