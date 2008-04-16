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

//SEQAN_xNO_GENERATED_FORWARDS: no forwards are generated for this file

#ifndef SEQAN_HEADER_PIZZACHILI_API_H
#define SEQAN_HEADER_PIZZACHILI_API_H

#include <seqan/basic.h>

namespace SEQAN_NAMESPACE_MAIN {

namespace impl {
    typedef unsigned char uchar_t;
    typedef unsigned long ulong_t;
    typedef void* index_t;
    typedef int error_t;
} // namespace impl

struct InvalidPizzaChiliSpec;

template <typename TSpec>
struct PizzaChiliCodeProvider {
    typedef InvalidPizzaChiliSpec Type;
};

/**
.Spec.PizzaChili_AF:
..summary:Pizza & Chili Alphabet-Friendly FM Index.
..cat:Index
..base:Spec.Pizza & Chili Index
..signature:Index<TText, PizzaChili<PizzaChili_AF> >
..param.TText:The text type.
...type:Class.String
..remarks:To be used in conjunction with Pizza & Chili's $af_index$.

.Spec.PizzaChili_CCSA:
..summary:Pizza & Chili Compressed Compact Suffix Array Index.
..cat:Index
..base:Spec.Pizza & Chili Index
..signature:Index<TText, PizzaChili<PizzaChili_CCSA> >
..param.TText:The text type.
...type:Class.String
..remarks:To be used in conjunction with Pizza & Chili's $CCSA$ index.

.Spec.PizzaChili_SA:
..summary:Pizza & Chili Suffix Array Index.
..cat:Index
..base:Spec.Pizza & Chili Index
..signature:Index<TText, PizzaChili<PizzaChili_SA> >
..param.TText:The text type.
...type:Class.String
..remarks:To be used in conjunction with Pizza & Chili's $SAu$ and $SAc$ index.

.Spec.PizzaChili_SSA:
..summary:Pizza & Chili Sufficient Suffix Array Index.
..cat:Index
..base:Spec.Pizza & Chili Index
..signature:Index<TText, PizzaChili<PizzaChili_SSA> >
..param.TText:The text type.
...type:Class.String
..remarks:To be used in conjunction with Pizza & Chili's $SSA*$ indices.
*/

// We need to declare these explicitly instead through macro expansion in order
// for them to be included in the forward generated declarations.

struct _PizzaChili_AF;
typedef Tag<_PizzaChili_AF> const PizzaChili_AF;

struct _PizzaChili_CCSA;
typedef Tag<_PizzaChili_CCSA> const PizzaChili_CCSA;

struct _PizzaChili_FM;
typedef Tag<_PizzaChili_FM> const PizzaChili_FM;

struct _PizzaChili_LZ;
typedef Tag<_PizzaChili_LZ> const PizzaChili_LZ;

struct _PizzaChili_RSA;
typedef Tag<_PizzaChili_RSA> const PizzaChili_RSA;

struct _PizzaChili_RLFM;
typedef Tag<_PizzaChili_RLFM> const PizzaChili_RLFM;

struct _PizzaChili_SA;
typedef Tag<_PizzaChili_SA> const PizzaChili_SA;

struct _PizzaChili_SADA;
typedef Tag<_PizzaChili_SADA> const PizzaChili_SADA;

struct _PizzaChili_SSA;
typedef Tag<_PizzaChili_SSA> const PizzaChili_SSA;

struct _PizzaChili_Test;
typedef Tag<_PizzaChili_Test> const PizzaChili_Test;

#define SEQAN_MAKE_PIZZACHILI_PROVIDER(name) \
    class PizzaChiliApi##name { \
    public: \
        static char* error_index(impl::error_t e); \
        static int build_index( \
            impl::uchar_t* text, \
            impl::ulong_t length, \
            char* build_options, \
            impl::index_t* index \
            ); \
        static int save_index(impl::index_t index, char* filename); \
        static int load_index(char* filename, impl::index_t* index); \
        static int free_index(impl::index_t index); \
        static int index_size(impl::index_t index, impl::ulong_t* size); \
        static int count( \
            impl::index_t index, \
            impl::uchar_t* pattern, \
            impl::ulong_t length, \
            impl::ulong_t* numocc \
        ); \
        static int locate( \
            impl::index_t index, \
            impl::uchar_t* pattern, \
            impl::ulong_t length, \
            impl::ulong_t** occ, \
            impl::ulong_t* numocc \
        ); \
        static int get_length(impl::index_t index, impl::ulong_t* length); \
        static int extract( \
            impl::index_t index, \
            impl::ulong_t from, \
            impl::ulong_t to, \
            impl::uchar_t** snippet, \
            impl::ulong_t* snippet_length \
        ); \
        static int display( \
            impl::index_t index, \
            impl::uchar_t* pattern, \
            impl::ulong_t length, \
            impl::ulong_t numc,  \
            impl::ulong_t* numocc, \
            impl::uchar_t** snippet_text, \
            impl::ulong_t** snippet_length \
        ); \
    }; \
    \
    /*struct _PizzaChili_##name; \
    typedef Tag<_PizzaChili_##name> const PizzaChili_##name;*/ \
    \
    template <> \
    struct PizzaChiliCodeProvider<PizzaChili_##name> { \
        typedef PizzaChiliApi##name Type; \
    };

SEQAN_MAKE_PIZZACHILI_PROVIDER(AF)
SEQAN_MAKE_PIZZACHILI_PROVIDER(CCSA)
SEQAN_MAKE_PIZZACHILI_PROVIDER(FM)
SEQAN_MAKE_PIZZACHILI_PROVIDER(LZ)
SEQAN_MAKE_PIZZACHILI_PROVIDER(RSA)
SEQAN_MAKE_PIZZACHILI_PROVIDER(RLFM)
SEQAN_MAKE_PIZZACHILI_PROVIDER(SA)
SEQAN_MAKE_PIZZACHILI_PROVIDER(SADA)
SEQAN_MAKE_PIZZACHILI_PROVIDER(SSA)
SEQAN_MAKE_PIZZACHILI_PROVIDER(Test)

#undef SEQAN_MAKE_PIZZACHILI_PROVIDER

} // namespace SEQAN_NAMESPACE_MAIN

#endif // SEQAN_HEADER_PIZZACHILI_API_H
