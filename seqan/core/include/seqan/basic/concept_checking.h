// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2011, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================
// A minimal subset of the Boost Concept Checking Library.  A lot of the code
// in the BCCL deals with support of non-conforming compilers and we cut this
// away.  The code here has been adjusted to work with the compilers supported
// by SeqAn and be as simple as possible while still creating useful compiler
// errors.
// ==========================================================================

// SEQAN_NO_GENERATED_FORWARDS

#ifndef CORE_INCLUDE_SEQAN_BASIC_CONCEPT_CHECKING_H_
#define CORE_INCLUDE_SEQAN_BASIC_CONCEPT_CHECKING_H_

namespace seqan {

// ---------------------------------------------------------------------------
// ==> boost/parameter/aux_/paranthesized_type.hpp <==
// ---------------------------------------------------------------------------

// Copyright David Abrahams 2006. Distributed under the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt).

template <class UnaryFunctionPointer>
struct unaryfunptr_arg_type;

template <class Arg>
struct unaryfunptr_arg_type<void(*)(Arg)>
{
    typedef Arg type; 
};

template <>
struct unaryfunptr_arg_type<void(*)(void)>
{
    typedef void type;
};

// ---------------------------------------------------------------------------
// ==> boost/concept_check/general.hpp <==
// ---------------------------------------------------------------------------

// Copyright David Abrahams 2006. Distributed under the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt).

namespace concept_checking
{
template <void(*)()> struct instantiate {};
}

template <class ModelFn> struct concept_check_;

template <class Model>
void concept_check_failed()
{
    ((Model*)0)->~Model();
}

template <class Model>
struct concept_check
{
    concept_checking::instantiate<concept_check_failed<Model> > x;
    enum { instantiate = 1 };
};

template <class Model>
struct concept_check_<void(*)(Model)>
        : concept_check<Model>
{};

#  define SEQAN_CONCEPT_ASSERT_FN( ModelFnPtr )             \
    typedef ::seqan::detail::instantiate<          \
    &::seqan::requirement_<ModelFnPtr>::failed>    \
      SEQAN_PP_CAT(seqan_concept_check,__LINE__)

// ---------------------------------------------------------------------------
// ==> boost/concept/assert.hpp <==
// ---------------------------------------------------------------------------

// Copyright David Abrahams 2006. Distributed under the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt).

// Usage, in class or function context:
//     SEQAN_CONCEPT_ASSERT((UnaryFunctionConcept<F,bool,int>));
# define SEQAN_CONCEPT_ASSERT(ModelInParens) \
    SEQAN_CONCEPT_ASSERT_FN(void(*)ModelInParens)

// usage.hpp

template <class Model>
struct usage_requirements
{
    ~usage_requirements() { ((Model*)0)->~Model(); }
};

#define SEQAN_CONCEPT_USAGE(model)                                      \
    SEQAN_CONCEPT_ASSERT((seqan::usage_requirements<model>));           \
    ~model()

// ---------------------------------------------------------------------------
// ==> boost/concept/detail/has_constraints.hpp <==
// ---------------------------------------------------------------------------

// Copyright David Abrahams 2006. Distributed under the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt).

namespace detail {
  typedef char yes;
  typedef char (&no)[2];

  template <class Model, void (Model::*)()>
  struct wrap_constraints {};
    
  template <class Model>
  inline yes has_constraints_(Model*, wrap_constraints<Model,&Model::constraints>* = 0);
  inline no has_constraints_(...);
}

// This would be called "detail::has_constraints," but it has a strong
// tendency to show up in error messages.
template <class Model>
struct not_satisfied
{
    enum {value = sizeof( detail::has_constraints_((Model*)0) ) == sizeof(detail::yes) };
    typedef typename Eval<value>::Type Type;
};

// ---------------------------------------------------------------------------
// ==> boost/concept_check/detail/general.hpp <==
// ---------------------------------------------------------------------------

// Copyright David Abrahams 2006. Distributed under the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt).

template <class ModelFn>
struct requirement_;

namespace detail
{
  template <void(*)()> struct instantiate {};
}

template <class Model>
struct requirement
{
    static void failed() { ((Model*)0)->~Model(); }
};

struct failed {};

template <class Model>
struct requirement<failed ************ Model::************>
{
    static void failed() { ((Model*)0)->~Model(); }
};

template <class Model>
struct constraint
{
    static void failed() { ((Model*)0)->constraints(); }
};
  
template <class Model>
struct requirement_<void(*)(Model)>
        : If<not_satisfied<Model>::Type::VALUE, /* should be called "has_constraints", see above */
             constraint<Model>,
             requirement<failed ************ Model::************>
             >::Type
{};

#  define SEQAN_CONCEPT_ASSERT_FN( ModelFnPtr )             \
    typedef ::seqan::detail::instantiate<          \
    &::seqan::requirement_<ModelFnPtr>::failed>    \
      SEQAN_PP_CAT(seqan_concept_check,__LINE__)

// ---------------------------------------------------------------------------
// ==> boost/concept_check/detail/requires.hpp <==
// ---------------------------------------------------------------------------

// Copyright David Abrahams 2006. Distributed under the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt).

// Template for use in handwritten assertions
template <class Model, class More>
struct requires_ : More
{
    SEQAN_CONCEPT_ASSERT((Model));
};

// Template for use by macros, where models must be wrapped in parens.
// This isn't in namespace detail to keep extra cruft out of resulting
// error messages.

template <class ModelFn>
struct _requires_
{
    enum { value = 0 };
    SEQAN_CONCEPT_ASSERT_FN(ModelFn);
};

template <int check, class Result>
struct Requires_ : unaryfunptr_arg_type<Result>
{};

#  define SEQAN_CONCEPT_REQUIRES_(r,data,t) + (::seqan::_requires_<void(*)t>::value)

#if defined(NDEBUG)

# define SEQAN_CONCEPT_REQUIRES(models, result)                      \
    typename unaryfunptr_arg_type<void(*)result>::type

#else  // #if defined(NDEBUG)

# define SEQAN_CONCEPT_REQUIRES(models, result)                                        \
    typename ::seqan::Requires_<                                                       \
      (0 SEQAN_PP_SEQ_FOR_EACH(SEQAN_CONCEPT_REQUIRES_, ~, models)),                   \
      void(*)result                                                                 \
    >::type

#endif  // #if defined(NDEBUG)

}  // namespace seqan

#endif  // #ifndef CORE_INCLUDE_SEQAN_BASIC_CONCEPT_CHECKING_H_
