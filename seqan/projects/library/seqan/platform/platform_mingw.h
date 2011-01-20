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

//SEQAN_NO_GENERATED_FORWARDS: no forwards are generated for this file

#define PLATFORM_WINDOWS
#define PLATFORM_WINDOWS_MINGW

#include <inttypes.h>

#define finline __inline__

// default 64bit type
typedef int64_t __int64;
typedef uint64_t __uint64;

// default 32bit type
typedef int32_t __int32;
typedef uint32_t __uint32;

// default 16bit type
typedef int16_t __int16;
typedef uint16_t __uint16;

// default 8bit type
typedef int8_t __int8;
typedef uint8_t __uint8;

//define SEQAN_SWITCH_USE_FORWARDS to use generated forwards 
#define SEQAN_SWITCH_USE_FORWARDS

#include "platform_generated_forwards.h"

#ifndef SEQAN_HEADER_PLATFORM_GENERATED_FORWARDS_H
#error To use the SeqAn library you first have to execute 'make forwards' in the root directory
#endif
