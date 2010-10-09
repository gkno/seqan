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
 Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
 ============================================================================
 This module contains simple, generic support code for parallelism in SeqAn.
 It mainly defines the macros SEQAN_ENABLE_PARALLELISM and
 SEQAN_PRAGMA_IF_PARALLEL.
 ==========================================================================*/
#ifndef SEQAN_BASIC_BASIC_PARALLELISM_H_
#define SEQAN_BASIC_BASIC_PARALLELISM_H_

/* By default, parallelism is enabled if OpenMP is activated.
 */
#if !defined(SEQAN_ENABLE_PARALLELISM)
#if defined(_OPENMP)
#define SEQAN_ENABLE_PARALLELISM 1
#else  // defined(_OPENMP)
#define SEQAN_ENABLE_PARALLELISM 0
#endif  // defined(_OPENMP)
#endif  // !defined(SEQAN_ENABLE_PARALLELISM)


/* The SEQAN_PRAGMA_IF_PARALLEL macro expands to a pragma as in the argument
 * if SEQAN_ENABLE_PARALLELISM is 1.
 */
#if SEQAN_ENABLE_PARALLELISM
// SEQAN_MKSTRING_() comes from basic_testing.h, so include before this header.
#define SEQAN_PRAGMA_IF_PARALLEL(arg) \
  _Pragma(SEQAN_MKSTRING_(arg))
#else  // SEQAN_ENABLE_PARALLELISM
#define SEQAN_PRAGMA_IF_PARALLEL
#endif  // SEQAN_ENABLE_PARALLELISM

#endif  // SEQAN_BASIC_BASIC_PARALLELISM_H_
