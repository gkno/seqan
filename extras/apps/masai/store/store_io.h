// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2010, Knut Reinert, FU Berlin
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
// Author: Enrico Siragusa <enrico.siragusa@fu-berlin.de>
// ==========================================================================
// This file replaces store.h
// ==========================================================================

#ifndef SEQAN_EXTRAS_APPS_MASAI_STORE_IO_H_
#define SEQAN_EXTRAS_APPS_MASAI_STORE_IO_H_

// ===========================================================================
// Prerequisites.
// ===========================================================================

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/align.h>
#include <seqan/file.h>
#include <seqan/bam_io.h>

#include <sstream>
#include <iostream>
#include <algorithm>
#include <set>

#ifndef SEQAN_HAS_SAMTOOLS
#define SEQAN_HAS_SAMTOOLS 0
#endif  // #ifndef SEQAN_HAS_SAMTOOLS

#if SEQAN_HAS_SAMTOOLS
#include <sam.h>
#include <bam.h>
#endif  // #if SEQAN_HAS_SAMTOOLS


// ===========================================================================
// Fragment Store Sub-Containers.
// ===========================================================================

#include <seqan/store/store_base.h>
#include <seqan/store/store_read.h>
#include <seqan/store/store_matepair.h>
#include <seqan/store/store_library.h>
#include <seqan/store/store_contig.h>
#include <seqan/store/store_align.h>
#include <seqan/store/store_annotation.h>
#include <seqan/store/store_all.h>

#include <seqan/store/store_align_intervals.h>
#include <seqan/store/store_intervaltree.h>

#include <seqan/store/store_io.h>
//#include <seqan/store/store_io_sam.h>
#include <seqan/store/store_io_gff.h>
#include <seqan/store/store_io_ucsc.h>
#if SEQAN_HAS_SAMTOOLS
#include <seqan/store/store_io_bam.h>
#endif  // #if SEQAN_HAS_SAMTOOLS


// ===========================================================================
// Fragment Store New IO.
// ===========================================================================

#include "store_io_sam.h"

#endif //#ifndef SEQAN_EXTRAS_APPS_MASAI_STORE_IO_H_
