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
// ARE DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================

#ifndef SEQAN_HEADER_SEQUENCE_H
#define SEQAN_HEADER_SEQUENCE_H

//____________________________________________________________________________
// prerequisites

#include <seqan/basic.h>
#include <map>
#include <vector>

//____________________________________________________________________________

#include <seqan/sequence/sequence_forwards.h>

#ifdef SEQAN_SWITCH_USE_FORWARDS
#include <seqan/sequence/sequence_generated_forwards.h>
#endif

#include <seqan/sequence/sequence_interface.h>
#include <seqan/sequence/sequence_stream.h>
#include <seqan/sequence/lexical.h>

//____________________________________________________________________________
// segments (suffix, ...)

#include <seqan/sequence/segment_base.h>
#include <seqan/sequence/segment_infix.h>
#include <seqan/sequence/segment_suffix.h>
#include <seqan/sequence/segment_prefix.h>

//____________________________________________________________________________
// strings

#include <seqan/sequence/string_base.h>
#include <seqan/sequence/string_pointer.h>
#include <seqan/sequence/string_alloc.h>
#include <seqan/sequence/string_array.h>
#include <seqan/sequence/string_cstyle.h>
#include <seqan/sequence/string_stack.h>
#include <seqan/sequence/string_packed.h>
#include <seqan/sequence/string_value_expand.h>

#include <seqan/sequence/std_string.h>
#include <seqan/sequence/std_vector.h>

#include <seqan/sequence/sequence_multiple.h>
#include <seqan/sequence/sequence_shortcuts.h>

//____________________________________________________________________________
// Adaptions

#include <seqan/sequence/adapt_std_list.h>

#endif //#ifndef SEQAN_HEADER_...
